# =============================================================================
# mod_occupancy -- Occupancy correction, reproducibility filter, tiering
#
# Receives norm_obj (from mod_normalization) and cfg.
# Returns a list with the SAME shape as the old combined norm_obj so that
# mod_de, mod_plots, and mod_export require no changes:
#   $y.sites, $y.tp          -- normalised ELists
#   $q.site, $q.prot, $q.occ -- dpcQuant + occupancy-corrected EList
#   $dpc.dg, $dpc.tp
#   $site_info, $link
#   $lines_present
# =============================================================================

mod_occupancy_ui <- function(id) {
  ns <- NS(id)
  tagList(

    # ── Run ─────────────────────────────────────────────────────────────────
    h4("Occupancy Correction & Reproducibility Filter"),
    card(
      card_body(
        tags$p(class = "text-muted small mb-2",
          "Computes occupancy = site log2 abundance \u2212 parent protein log2 abundance ",
          "(with SEs combined in quadrature). Then applies the per-line reproducibility ",
          "filter and assigns each site a tier."),
        tags$div(style = "display:flex; align-items:center; gap:12px; margin-bottom:12px;",
          actionButton(ns("run"), "Run occupancy correction",
                       class = "btn-primary"),
          spinner_ui(ns("sp_occ"))
        ),
        verbatimTextOutput(ns("occ_log"))
      )
    ),

    # ── Occupancy correction summary ─────────────────────────────────────────
    conditionalPanel(
      condition = sprintf("output['%s']", ns("done")),

      h4("Site\u2013Protein Matching"),
      card(
        card_body(
          tags$p(class = "text-muted small mb-1",
            tags$b("exact:"), " site protein group matched directly. ",
            tags$b("leading_acc:"), " matched via leading accession. ",
            tags$b("unmatched:"), " no protein group found \u2014 site tested on raw layer only."),
          verbatimTextOutput(ns("match_summary"))
        )
      ),

      # ── Tier table ──────────────────────────────────────────────────────────
      h4("Reproducibility Filter & Tier Summary"),
      card(
        card_body(
          tags$p(class = "text-muted small mb-2",
            tags$b("corrected:"), " passes repro filter + parent protein quantified \u2192 ",
            "occupancy-corrected DE. ",
            tags$b("no_protein_group / protein_not_quantified:"), " site-layer DE only. ",
            tags$b("low_reproducibility:"), " excluded from DE."),
          fluidRow(
            column(4, verbatimTextOutput(ns("tier_summary"))),
            column(8, DT::dataTableOutput(ns("tier_table")))
          )
        )
      ),

      # ── Occupancy matrix boxplot ─────────────────────────────────────────────
      h4("Occupancy Distribution (post-correction)"),
      card(
        card_body(
          plotOutput(ns("occ_boxplot"), height = "320px")
        )
      )
    )
  )
}


mod_occupancy_server <- function(id, norm_obj, cfg, status = reactiveVal("idle")) {
  moduleServer(id, function(input, output, session) {

    log_msg      <- reactiveVal("")
    occ_res      <- reactiveVal(NULL)
    sp_occ_state <- reactiveVal("idle")
    spinner_server("sp_occ", sp_occ_state)

    output$done <- reactive(!is.null(occ_res()))
    outputOptions(output, "done", suspendWhenHidden = FALSE)

    append_log <- function(...) {
      msg <- paste0(paste(...), "\n")
      log_msg(paste0(log_msg(), msg))
    }

    observeEvent(input$run, {
      nm <- norm_obj();  c <- cfg()
      req(nm, c)
      log_msg(""); occ_res(NULL)
      status("running"); sp_occ_state("running")

      withCallingHandlers(
        tryCatch({
          q.site  <- nm$q.site
          q.prot  <- nm$q.prot
          y.sites <- nm$y.sites

          prog <- make_progress_stepper(session, "Occupancy correction\u2026", c(
            "Site\u2013protein matching"  = 0.25,
            "Occupancy subtraction"       = 0.55,
            "Reproducibility filter"      = 0.80,
            "Finalising"                  = 1.00))
          on.exit(prog$close(), add = TRUE)

          # ── Site-protein link table ────────────────────────────────────────
          append_log("=== Site\u2013protein matching ===")
          prog$tick("Building site\u2013protein link table\u2026")

          site2pg <- dplyr::tibble(site_id = rownames(q.site)) |>
            dplyr::bind_cols(dplyr::as_tibble(q.site$genes)) |>
            dplyr::transmute(
              site_id, site_label, Genes, Protein.Group,
              lead = stringr::str_split_fixed(Protein.Group, ";", 2)[, 1])

          pg_index <- dplyr::tibble(
            pg   = rownames(q.prot),
            lead = stringr::str_split_fixed(rownames(q.prot), ";", 2)[, 1])

          link <- site2pg |>
            dplyr::left_join(
              pg_index |> dplyr::transmute(Protein.Group = pg, pg_exact = pg),
              by = "Protein.Group") |>
            dplyr::left_join(
              pg_index |> dplyr::transmute(lead, pg_lead = pg),
              by = "lead", relationship = "many-to-many") |>
            dplyr::group_by(site_id) |>
            dplyr::summarise(
              site_label    = dplyr::first(site_label),
              Genes         = dplyr::first(Genes),
              Protein.Group = dplyr::first(Protein.Group),
              pg_matched    = dplyr::coalesce(
                dplyr::first(na.omit(pg_exact)),
                dplyr::first(na.omit(pg_lead))),
              match_type    = dplyr::case_when(
                !is.na(dplyr::first(na.omit(pg_exact))) ~ "exact",
                !is.na(dplyr::first(na.omit(pg_lead)))  ~ "leading_acc",
                TRUE                                      ~ "unmatched"),
              .groups = "drop")

          append_log(sprintf("  Match types: %s",
            paste(names(table(link$match_type)),
                  table(link$match_type), sep = "=", collapse = " | ")))

          # ── Occupancy subtraction ──────────────────────────────────────────
          append_log("\n=== Occupancy correction ===")
          prog$tick("Computing site \u2212 protein log2 differences\u2026")

          mapped   <- link |>
            dplyr::filter(match_type != "unmatched",
                          pg_matched %in% rownames(q.prot))
          idx_site <- match(mapped$site_id,    rownames(q.site))
          idx_prot <- match(mapped$pg_matched, rownames(q.prot))

          E.occ  <- q.site$E[idx_site, , drop = FALSE] -
                    q.prot$E[idx_prot, , drop = FALSE]
          SE.occ <- sqrt(
            q.site$other$standard.error[idx_site, , drop = FALSE]^2 +
            q.prot$other$standard.error[idx_prot, , drop = FALSE]^2)
          N.occ  <- pmin(
            q.site$other$n.observations[idx_site, , drop = FALSE],
            q.prot$other$n.observations[idx_prot, , drop = FALSE])
          rownames(E.occ) <- rownames(SE.occ) <- rownames(N.occ) <- mapped$site_id

          q.occ <- methods::new("EList", list(
            E       = E.occ,
            genes   = as.data.frame(mapped),
            targets = y.sites$targets,
            other   = list(standard.error = SE.occ, n.observations = N.occ)))
          q.occ$dpc <- q.site$dpc

          append_log(sprintf("  Occupancy-corrected sites: %d of %d",
                             nrow(q.occ), nrow(q.site)))

          # ── Reproducibility filter + tiering ──────────────────────────────
          append_log("\n=== Reproducibility filter & tiering ===")
          prog$tick("Applying reproducibility filter\u2026")

          targets     <- y.sites$targets
          line_vec    <- as.character(targets$cell_line)
          lines_pres  <- levels(droplevels(targets$cell_line))
          target_line <- c$target_line
          rest_lines  <- c$reference_lines

          obs_per_line <- function(N)
            sapply(lines_pres, function(cl)
              rowSums(N[, line_vec == cl, drop = FALSE] > 0))

          obs_site <- obs_per_line(q.site$other$n.observations)

          site_info <- dplyr::tibble(site_id = rownames(q.site)) |>
            dplyr::left_join(link, by = "site_id") |>
            dplyr::mutate(
              n_obs_target   = obs_site[, target_line],
              n_obs_rest_max = apply(
                obs_site[, intersect(rest_lines, colnames(obs_site)), drop = FALSE],
                1, max),
              n_lines_repro  = rowSums(obs_site >= c$min_obs_per_line),
              passes_repro   = (n_obs_target   >= c$min_obs_per_line) |
                               (n_obs_rest_max >= c$min_obs_per_line),
              tier = dplyr::case_when(
                !passes_repro                 ~ "low_reproducibility",
                site_id %in% rownames(q.occ)  ~ "corrected",
                match_type == "unmatched"      ~ "no_protein_group",
                TRUE                           ~ "protein_not_quantified"))

          tier_tbl <- dplyr::count(site_info, tier, sort = TRUE)
          append_log(sprintf("  Tier counts: %s",
            paste(tier_tbl$tier, tier_tbl$n, sep = "=", collapse = " | ")))

          keep_sites <- site_info$site_id[site_info$passes_repro]
          q.site_f   <- q.site[rownames(q.site) %in% keep_sites, ]
          q.occ_f    <- q.occ[ rownames(q.occ)  %in% keep_sites, ]

          append_log(sprintf(
            "  After repro filter (>=%d obs/line): %d sites, %d corrected",
            c$min_obs_per_line, nrow(q.site_f), nrow(q.occ_f)))

          prog$tick("Finalising\u2026")
          append_log("\n[OK] Occupancy correction complete.")
          status("done"); sp_occ_state("done")

          # Return the full norm_obj-shaped list expected by DE/plots/export
          occ_res(list(
            y.sites       = y.sites,
            y.tp          = nm$y.tp,
            q.site        = q.site_f,
            q.prot        = q.prot,
            q.occ         = q.occ_f,
            dpc.dg        = nm$dpc.dg,
            dpc.tp        = nm$dpc.tp,
            site_info     = site_info,
            link          = link,
            lines_present = lines_pres
          ))

        }, error = function(e) {
          append_log(paste0("[ERR] ERROR: ", e$message))
          status("error"); sp_occ_state("error")
        }),
        message = function(m) {
          append_log(conditionMessage(m)); invokeRestart("muffleMessage")
        },
        warning = function(w) {
          append_log(paste("[!]", conditionMessage(w))); invokeRestart("muffleWarning")
        }
      )
    })

    # ── Outputs ───────────────────────────────────────────────────────────────
    output$occ_log <- renderPrint({ cat(log_msg()) })

    output$match_summary <- renderPrint({
      res <- occ_res(); req(res)
      print(dplyr::count(res$link, match_type, sort = TRUE))
    })

    output$tier_summary <- renderPrint({
      res <- occ_res(); req(res)
      cat("Tier counts:\n")
      print(dplyr::count(res$site_info, tier, sort = TRUE))
      cat(sprintf("\nPasses repro: %d / %d sites\n",
                  sum(res$site_info$passes_repro),
                  nrow(res$site_info)))
      cat(sprintf("Corrected (with occ): %d\n",
                  sum(res$site_info$tier == "corrected")))
    })

    output$tier_table <- DT::renderDataTable({
      res <- occ_res(); req(res)
      res$site_info |>
        dplyr::select(site_id, site_label, Genes, tier, match_type,
                      n_obs_target, n_obs_rest_max, n_lines_repro, passes_repro) |>
        dplyr::arrange(tier, dplyr::desc(n_obs_target))
    }, options  = list(scrollX = TRUE, pageLength = 15,
                       columnDefs = list(list(width = "130px", targets = 0))),
       rownames = FALSE,
       filter   = "top")

    output$occ_boxplot <- renderPlot({
      res <- occ_res(); req(res)
      targets <- res$y.sites$targets
      all_cl  <- levels(factor(as.character(targets$cell_line)))
      pal     <- grDevices::colorRampPalette(
        c("#e6ab02","#d95f02","#7570b3","#e7298a",
          "#1b9e77","#66a61e","#a6761d","#666666"))(length(all_cl))
      line_colors <- setNames(pal, all_cl)
      meta_df <- data.frame(name      = rownames(targets),
                            cell_line = as.character(targets$cell_line))

      as.data.frame(res$q.occ$E) |>
        tibble::rownames_to_column("site_id") |>
        tidyr::pivot_longer(-site_id, names_to = "name") |>
        dplyr::left_join(meta_df, by = "name") |>
        ggplot2::ggplot(ggplot2::aes(name, value, fill = cell_line)) +
        ggplot2::geom_boxplot(outlier.size = 0.2, na.rm = TRUE) +
        ggplot2::scale_fill_manual(values = line_colors) +
        ggplot2::labs(x = NULL, y = "log2 occupancy (site \u2212 protein)",
                      title = "Occupancy distribution per sample") +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    })

    return(occ_res)
  })
}
