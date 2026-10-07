# cRAP contaminant accessions (camprotR 2024-05-12 + cRAP 2019-04-01)
.CRAP_ACCESSIONS <- c(
  "cRAP-GluC", "cRAP-O00762", "cRAP-O43790", "cRAP-O76009", "cRAP-O76011", "cRAP-O76013",
  "cRAP-O76014", "cRAP-O76015", "cRAP-O77727", "cRAP-O82803", "cRAP-P00004", "cRAP-P00167",
  "cRAP-P00330", "cRAP-P00366", "cRAP-P00441", "cRAP-P00698", "cRAP-P00709", "cRAP-P00711",
  "cRAP-P00722", "cRAP-P00760", "cRAP-P00761", "cRAP-P00766", "cRAP-P00767", "cRAP-P00772",
  "cRAP-P00791", "cRAP-P00792", "cRAP-P00883", "cRAP-P00915", "cRAP-P00918", "cRAP-P00921",
  "cRAP-P01008", "cRAP-P01012", "cRAP-P01031", "cRAP-P01112", "cRAP-P01127", "cRAP-P01133",
  "cRAP-P01344", "cRAP-P01375", "cRAP-P02144", "cRAP-P02438", "cRAP-P02439", "cRAP-P02440",
  "cRAP-P02441", "cRAP-P02443", "cRAP-P02444", "cRAP-P02445", "cRAP-P02534", "cRAP-P02539",
  "cRAP-P02662", "cRAP-P02663", "cRAP-P02666", "cRAP-P02668", "cRAP-P02741", "cRAP-P02753",
  "cRAP-P02768", "cRAP-P02769", "cRAP-P02787", "cRAP-P02788", "cRAP-P04040", "cRAP-P04264",
  "cRAP-P04745", "cRAP-P05413", "cRAP-P06396", "cRAP-P06732", "cRAP-P07339", "cRAP-P08131",
  "cRAP-P08263", "cRAP-P08311", "cRAP-P08758", "cRAP-P09211", "cRAP-P0AEX9", "cRAP-P0C1U8",
  "cRAP-P10145", "cRAP-P10599", "cRAP-P10636", "cRAP-P12081", "cRAP-P12763", "cRAP-P13645",
  "cRAP-P15241", "cRAP-P15252", "cRAP-P15559", "cRAP-P16083", "cRAP-P21578", "cRAP-P22629",
  "cRAP-P25690", "cRAP-P25691", "cRAP-P26372", "cRAP-P30879", "cRAP-P32503", "cRAP-P35527",
  "cRAP-P35908", "cRAP-P41159", "cRAP-P42212", "cRAP-P51965", "cRAP-P55957", "cRAP-P61626",
  "cRAP-P61769", "cRAP-P62937", "cRAP-P62979", "cRAP-P63165", "cRAP-P63279", "cRAP-P68082",
  "cRAP-P68871", "cRAP-P69905", "cRAP-P78385", "cRAP-P78386", "cRAP-P81054", "cRAP-P99999",
  "cRAP-Q02958", "cRAP-Q06830", "cRAP-Q10735", "cRAP-Q14525", "cRAP-Q14532", "cRAP-Q14533",
  "cRAP-Q15323", "cRAP-Q15843", "cRAP-Q28153", "cRAP-Q29461", "cRAP-Q29463", "cRAP-Q58D62",
  "cRAP-Q7M135", "cRAP-Q92764", "cRAP-Q9NSB2", "cRAP-Q9NSB4", "cRAP-rLys-C"
)
.GFP_ACCESSION <- "cRAP-P42212"

# =============================================================================
# mod_filtering -- GG filtering, missingness filter, site mapping, DPC
# =============================================================================

mod_filtering_ui <- function(id) {
  ns <- NS(id)
  tagList(

    # ── Missingness filter options (set BEFORE running) ─────────────────────
    h4("Contaminant Filtering"),
    card(
      card_body(
        tags$p(class = "text-muted small mb-2",
          "Removes precursors whose Protein.Group matches contaminant accessions. ",
          "Applied before all other filters."),
        fluidRow(
          column(4,
            checkboxInput(ns("filter_crap"), "Remove cRAP contaminants", value = TRUE),
            helpText("Built-in list: 125 accessions (camprotR cRAP 2024-05-12). ",
                     "Matches any Protein.Group containing 'cRAP-'.")
          ),
          column(4,
            checkboxInput(ns("keep_gfp"), "Keep GFP (cRAP-P42212)", value = FALSE),
            helpText(icon("info-circle"),
              "When ", tags$b("Remove cRAP is ON"),
              ", GFP (cRAP-P42212) is excluded along with all other cRAP entries. ",
              "Tick this to ", tags$b("rescue GFP precursors"), " — all other cRAP ",
              "are still removed.")
          ),
          column(4,
            checkboxInput(ns("filter_custom_cont"), "Remove custom accessions",
                          value = FALSE),
            conditionalPanel(
              condition = sprintf("input['%s']", ns("filter_custom_cont")),
              textAreaInput(ns("custom_cont_accs"),
                            "Substrings to match (one per line)",
                            value = "", rows = 3,
                            placeholder = "e.g. REV__
CON__"),
              helpText("Any Protein.Group containing these strings is removed.")
            )
          )
        )
      )
    ),

    h4("Precursor Missingness Filter"),
    card(
      card_body(
        tags$p(class = "text-muted small mb-2",
          "Applied after identification filters and ", tags$b("before site mapping"),
          " — acts on raw precursor intensities. ",
          "Controls how many samples must have a valid (non-NA) intensity value ",
          "for a precursor to be kept. ",
          tags$b("Not the same as Min obs / line"),
          " (Configuration tab), which is a post-DPC-Quant reproducibility filter ",
          "applied per cell line on the quantified site matrix in the Occupancy tab."),
        fluidRow(
          column(4,
            radioButtons(ns("missing_mode"), "Filter mode",
              choices = c(
                "Minimum valid values (absolute)"  = "abs",
                "Minimum fraction of samples"      = "frac",
                "Allow all missing (no filter)"    = "none"),
              selected = "abs")
          ),
          column(4,
            conditionalPanel(
              condition = sprintf("input['%s'] == 'abs'", ns("missing_mode")),
              numericInput(ns("min_valid_abs"),
                           "Min valid values across all samples",
                           value = 1, min = 0, step = 1),
              helpText("0 = keep even if all samples are NA (same as 'allow all missing').")
            ),
            conditionalPanel(
              condition = sprintf("input['%s'] == 'frac'", ns("missing_mode")),
              sliderInput(ns("min_valid_frac"),
                          "Min fraction of samples with valid value",
                          min = 0, max = 1, value = 0.5, step = 0.05),
              helpText("E.g. 0.5 = present in at least 50% of samples.")
            ),
            conditionalPanel(
              condition = sprintf("input['%s'] == 'none'", ns("missing_mode")),
              helpText(icon("info-circle"),
                "All precursors are kept regardless of missingness. ",
                "DPC will model the missing values probabilistically.")
            )
          ),
          column(4,
            helpText(icon("info-circle"),
              tags$b("Note:"), " limpa uses a detection probability curve (DPC) ",
              "to represent missing values — you do not need to remove all NAs. ",
              "This filter removes precursors that are ", tags$i("entirely"), " absent ",
              "or below a reproducibility threshold before site mapping.")
          )
        )
      )
    ),

    # ── Run ─────────────────────────────────────────────────────────────────
    h4("Apply Filters & Map Sites"),
    card(
      card_body(
        tags$div(style = "display:flex; align-items:center; gap:12px; margin-bottom:12px;",
          actionButton(ns("run"), "Run filtering + site mapping",
                       class = "btn-primary"),
          spinner_ui(ns("sp_filt"))
        ),
        verbatimTextOutput(ns("filter_log"))
      )
    ),

    # ── QC outputs ──────────────────────────────────────────────────────────
    h4("Localisation Confidence Distribution"),
    card(
      card_body(
        plotOutput(ns("conf_hist"), height = "280px")
      )
    ),

    h4("Site Summary"),
    card(
      card_body(
        verbatimTextOutput(ns("site_summary"))
      )
    ),

    h4("Detection Probability Curves"),
    card(
      card_body(
        plotOutput(ns("dpc_plot"), height = "320px")
      )
    )
  )
}

mod_filtering_server <- function(id, data_obj, cfg, status = reactiveVal("idle")) {
  moduleServer(id, function(input, output, session) {

    log_msg       <- reactiveVal("")
    filt_res      <- reactiveVal(NULL)
    sp_filt_state <- reactiveVal("idle")
    spinner_server("sp_filt", sp_filt_state)

    append_log <- function(...) {
      msg <- paste0(paste(...), "\n")
      log_msg(paste0(log_msg(), msg))
    }

    add_filter <- function(keep, condition, label) {
      condition[is.na(condition)] <- FALSE
      new_keep <- keep & condition
      append_log(sprintf("  %-36s %7d -> %7d  (-%d)",
                         label, sum(keep), sum(new_keep), sum(keep) - sum(new_keep)))
      new_keep
    }

    gg_positions <- function(mod_seq, gg_tag = "UniMod:121") {
      m    <- gregexpr("\\((UniMod|Unimod):[0-9]+\\)", mod_seq, perl = TRUE)
      lens <- lapply(m, function(x) attr(x, "match.length"))
      full <- paste0("(", gg_tag, ")")
      mapply(function(s, st, ln) {
        if (st[1] == -1) return(integer(0))
        consumed <- cumsum(c(0, ln))[seq_along(ln)]
        res_idx  <- st - 1 - consumed
        tok      <- substring(s, st, st + ln - 1)
        idx <- res_idx[tok == full]
        idx[idx >= 1]
      }, mod_seq, m, lens, SIMPLIFY = FALSE, USE.NAMES = FALSE)
    }

    find_peptide_factory <- function(prot_seq) {
      function(pep, ids) {
        for (a in strsplit(ids, ";")[[1]]) {
          s <- prot_seq[[a]]
          if (is.null(s)) next
          h <- gregexpr(pep, s, fixed = TRUE)[[1]]
          if (h[1] != -1) return(c(a, h[1], length(h), nchar(s)))
        }
        c(NA, NA, 0L, NA)
      }
    }

    # ── Main computation ─────────────────────────────────────────────────────
    observeEvent(input$run, {
      d <- data_obj();  c <- cfg()
      req(d, c)
      log_msg(""); filt_res(NULL)
      status("running"); sp_filt_state("running")

      withCallingHandlers(
        tryCatch({
          y.dg <- d$y.dg;  y.tp <- d$y.tp
          ubi_genes <- c("UBB", "UBC", "RPS27A", "UBA52", "UBBP4")

          prog <- make_progress_stepper(session, "Filtering & mapping\u2026", c(
            "diGly precursor filters"   = 0.18,
            "Missingness filter (diGly)"= 0.28,
            "Proteome filters"          = 0.38,
            "Missingness filter (prot)" = 0.45,
            "Mapping GG sites to FASTA" = 0.72,
            "DPC estimation (diGly)"    = 0.86,
            "DPC estimation (proteome)" = 1.00))
          on.exit(prog$close(), add = TRUE)

          # ── Contaminant filtering (both layers) ────────────────────────────
          if (input$filter_crap || input$filter_custom_cont) {
            append_log("=== Contaminant filtering ===")
            cont_patterns <- character(0)
            if (input$filter_crap)
              cont_patterns <- c(cont_patterns, "cRAP-")
            if (input$filter_custom_cont) {
              user_accs <- Filter(nchar,
                                  trimws(strsplit(input$custom_cont_accs, "\n")[[1]]))
              cont_patterns <- c(cont_patterns, user_accs)
            }
            if (length(cont_patterns) > 0) {
              # Use fixed-string matching per pattern; rescue GFP if requested
              is_cont <- function(pg) {
                hits <- Reduce(`|`, lapply(cont_patterns, function(p)
                  stringr::str_detect(pg, stringr::fixed(p))))
                hits[is.na(hits)] <- FALSE
                # If Keep GFP is ticked, un-flag GFP precursors
                if (input$keep_gfp) {
                  is_gfp <- stringr::str_detect(pg, stringr::fixed(.GFP_ACCESSION))
                  is_gfp[is.na(is_gfp)] <- FALSE
                  hits <- hits & !is_gfp
                }
                hits
              }
              is_cont_dg <- is_cont(y.dg$genes$Protein.Group)
              n_rm_dg    <- sum(is_cont_dg)
              y.dg       <- y.dg[!is_cont_dg, ]
              append_log(sprintf("  diGly:    removed %d contaminant precursors (%d remaining)%s",
                                 n_rm_dg, nrow(y.dg),
                                 if (input$keep_gfp) " [GFP kept]" else ""))
              is_cont_tp <- is_cont(y.tp$genes$Protein.Group)
              n_rm_tp    <- sum(is_cont_tp)
              y.tp       <- y.tp[!is_cont_tp, ]
              append_log(sprintf("  Proteome: removed %d contaminant precursors (%d remaining)%s",
                                 n_rm_tp, nrow(y.tp),
                                 if (input$keep_gfp) " [GFP kept]" else ""))
            }
          }

          # Refresh genes slot after contaminant removal
          g <- y.dg$genes

          # ── diGly identification filters ────────────────────────────────
          append_log("=== diGly precursor filters ===")
          prog$tick("Applying diGly identification filters\u2026")
          keep <- rep(TRUE, nrow(y.dg))

          # Precursor q-value
          q_cols_dg <- intersect(
            c("Q.Value", "Lib.Q.Value", "Lib.PG.Q.Value", "PG.Q.Value"),
            names(g))
          if (length(q_cols_dg) > 0) {
            q_min <- apply(g[, q_cols_dg, drop = FALSE], 1,
                           function(x) min(as.numeric(x), na.rm = TRUE))
            keep <- add_filter(keep, q_min <= c$q_value,
                               sprintf("Precursor q-value <= %g", c$q_value))
          } else {
            append_log("  (no q-value columns found in diGly genes slot)")
          }

          keep <- add_filter(keep,
            stringr::str_detect(g$Modified.Sequence, stringr::fixed(c$gg_unimod)),
            "GG (UniMod:121)")

          if ("PTM.Q.Value" %in% names(g))
            keep <- add_filter(keep, g$PTM.Q.Value <= c$ptm_q_value,
                               sprintf("PTM q-value <= %g", c$ptm_q_value))

          if ("PTM.Site.Confidence" %in% names(g))
            keep <- add_filter(keep, g$PTM.Site.Confidence >= c$site_conf,
                               sprintf("Site confidence >= %g", c$site_conf))

          if (c$drop_cterm_gg)
            keep <- add_filter(keep,
              !stringr::str_detect(g$Modified.Sequence,
                                   paste0("\\(", c$gg_unimod, "\\)$")),
              "Not on peptide C-term")

          if (c$drop_ubiquitin)
            keep <- add_filter(keep, !(g$Genes %in% ubi_genes), "Not ubiquitin gene")

          if ("Proteotypic" %in% names(g))
            keep <- add_filter(keep, as.logical(g$Proteotypic) %in% TRUE, "Proteotypic")

          y.dg <- y.dg[keep, ]

          # ── diGly missingness filter ────────────────────────────────────
          append_log("\n=== diGly missingness filter ===")
          prog$tick("Applying diGly missingness filter\u2026")
          n_samples_dg <- ncol(y.dg$E)
          nobs_dg      <- rowSums(!is.na(y.dg$E))

          min_val_dg <- switch(input$missing_mode,
            abs  = input$min_valid_abs,
            frac = ceiling(input$min_valid_frac * n_samples_dg),
            none = 0L)

          if (input$missing_mode == "none") {
            append_log(sprintf(
              "  Mode: allow all missing — all %d precursors kept", nrow(y.dg)))
          } else {
            y.dg <- y.dg[nobs_dg >= max(1L, as.integer(min_val_dg)), ]
            append_log(sprintf(
              "  Mode: %s | threshold: %d / %d samples | retained: %d",
              input$missing_mode, max(1L, as.integer(min_val_dg)),
              n_samples_dg, nrow(y.dg)))
          }
          append_log(sprintf("  diGly precursors retained: %d", nrow(y.dg)))

          # ── Proteome identification filters ─────────────────────────────
          append_log("\n=== Proteome precursor filters ===")
          prog$tick("Applying proteome identification filters\u2026")
          gp      <- y.tp$genes
          keep_tp <- rep(TRUE, nrow(y.tp))

          q_cols_tp <- intersect(
            c("Q.Value", "Lib.Q.Value", "Lib.PG.Q.Value", "PG.Q.Value"),
            names(gp))
          if (length(q_cols_tp) > 0) {
            q_min_tp <- apply(gp[, q_cols_tp, drop = FALSE], 1,
                              function(x) min(as.numeric(x), na.rm = TRUE))
            keep_tp <- add_filter(keep_tp, q_min_tp <= c$q_value,
                                  sprintf("Precursor q-value <= %g", c$q_value))
          }

          if ("Proteotypic" %in% names(gp))
            keep_tp <- add_filter(keep_tp,
              as.logical(gp$Proteotypic) %in% TRUE, "Proteotypic")

          keep_tp <- add_filter(keep_tp,
            !stringr::str_detect(gp$Protein.Group, ";"), "Single-protein group")
          y.tp <- y.tp[keep_tp, ]

          # ── Proteome missingness filter ─────────────────────────────────
          append_log("\n=== Proteome missingness filter ===")
          prog$tick("Applying proteome missingness filter\u2026")
          n_samples_tp <- ncol(y.tp$E)
          nobs_tp      <- rowSums(!is.na(y.tp$E))

          min_val_tp <- switch(input$missing_mode,
            abs  = input$min_valid_abs,
            frac = ceiling(input$min_valid_frac * n_samples_tp),
            none = 0L)

          if (input$missing_mode == "none") {
            append_log(sprintf(
              "  Mode: allow all missing — all %d precursors kept", nrow(y.tp)))
          } else {
            y.tp <- y.tp[nobs_tp >= max(1L, as.integer(min_val_tp)), ]
            append_log(sprintf(
              "  Mode: %s | threshold: %d / %d samples | retained: %d",
              input$missing_mode, max(1L, as.integer(min_val_tp)),
              n_samples_tp, nrow(y.tp)))
          }
          append_log(sprintf("  Proteome precursors retained: %d", nrow(y.tp)))

          # ── Map GG sites onto protein sequences ─────────────────────────
          append_log("\n=== Mapping GG sites onto protein sequences ===")
          prog$tick("Mapping GG sites onto protein sequences\u2026")
          ann <- dplyr::as_tibble(y.dg$genes) |>
            dplyr::mutate(row_id     = seq_len(dplyr::n()),
                          gg_pos_pep = gg_positions(Modified.Sequence, c$gg_unimod))

          find_peptide <- find_peptide_factory(d$prot_seq)
          hits <- Map(find_peptide, ann$Stripped.Sequence, ann$Protein.Group)
          ann <- ann |>
            dplyr::mutate(
              map_acc   = vapply(hits, function(x) x[1], character(1)),
              pep_start = as.integer(vapply(hits, function(x) x[2], character(1))),
              n_hits    = as.integer(vapply(hits, function(x) x[3], character(1))),
              prot_len  = as.integer(vapply(hits, function(x) x[4], character(1))))

          append_log(sprintf("  unmapped: %d (%.2f%%)",
            sum(is.na(ann$map_acc)), 100 * mean(is.na(ann$map_acc))))

          sites <- ann |>
            dplyr::filter(!is.na(map_acc)) |>
            dplyr::mutate(site_pos = Map(function(st, p) st + p - 1L,
                                          pep_start, gg_pos_pep)) |>
            dplyr::select(row_id, Protein.Group, Genes, map_acc, prot_len,
                           n_hits, site_pos) |>
            tidyr::unnest(site_pos) |>
            dplyr::mutate(
              site_id       = paste0(map_acc, "_K", site_pos),
              site_label    = ifelse(is.na(Genes) | Genes == "",
                                      site_id, paste0(Genes, "_K", site_pos)),
              is_prot_cterm = site_pos == prot_len,
              ambiguous_pos = n_hits > 1)

          append_log(sprintf("  Unique diGly sites: %d",
                             dplyr::n_distinct(sites$site_id)))

          y.sites <- y.dg[sites$row_id, ]
          y.sites$genes <- cbind(y.sites$genes,
            sites[, c("site_id", "site_label", "map_acc",
                       "is_prot_cterm", "ambiguous_pos")])
          rownames(y.sites) <- make.unique(
            paste0(sites$site_id, "|", sites$row_id))
          y.sites$targets <- d$targets

          # ── DPC estimation ──────────────────────────────────────────────
          append_log("\n=== Detection probability curves ===")
          prog$tick("Estimating DPC (diGly layer) \u2014 may take 1\u20132 min\u2026")
          dpc.dg <- limpa::dpc(y.sites)
          prog$tick("Estimating DPC (proteome layer)\u2026")
          dpc.tp <- limpa::dpc(y.tp)
          append_log("  DPCs estimated.")
          append_log("\n[OK] Filtering + site mapping complete.")
          status("done"); sp_filt_state("done")

          filt_res(list(
            y.sites = y.sites, y.tp   = y.tp,
            sites   = sites,   ann    = ann,
            dpc.dg  = dpc.dg,  dpc.tp = dpc.tp
          ))

        }, error = function(e) {
          append_log(paste0("[ERR] ERROR: ", e$message))
          status("error"); sp_filt_state("error")
        }),
        message = function(m) {
          append_log(conditionMessage(m)); invokeRestart("muffleMessage")
        },
        warning = function(w) {
          append_log(paste("[!]", conditionMessage(w))); invokeRestart("muffleWarning")
        }
      )
    })

    # ── Outputs ──────────────────────────────────────────────────────────────
    output$filter_log <- renderPrint({ cat(log_msg()) })

    output$conf_hist <- renderPlot({
      res <- filt_res(); req(res)
      g   <- res$y.sites$genes
      if (!"PTM.Site.Confidence" %in% names(g)) {
        plot.new(); text(0.5, 0.5, "PTM.Site.Confidence not available"); return()
      }
      c <- cfg()
      hist(g$PTM.Site.Confidence, breaks = 50, col = "grey30", border = "white",
           main = "Localisation confidence (post-filter)",
           xlab = "PTM.Site.Confidence")
      abline(v = c$site_conf, col = "#e66101", lty = 2, lwd = 2)
    })

    output$site_summary <- renderPrint({
      res <- filt_res(); req(res)
      cat(sprintf("diGly precursors (post-filter)  : %d\n", nrow(res$y.sites)))
      cat(sprintf("Unique diGly sites               : %d\n",
                  dplyr::n_distinct(res$sites$site_id)))
      cat(sprintf("Proteome precursors (post-filter): %d\n", nrow(res$y.tp)))
    })

    output$dpc_plot <- renderPlot({
      res <- filt_res(); req(res)
      par(mfrow = c(1, 2))
      limpa::plotDPC(res$dpc.dg, main = "diGly layer")
      limpa::plotDPC(res$dpc.tp, main = "Total proteome")
      par(mfrow = c(1, 1))
    })

    return(filt_res)
  })
}
