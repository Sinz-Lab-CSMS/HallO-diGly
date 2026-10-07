# =============================================================================
# mod_de -- Differential analysis: one-vs-rest OR one-vs-control
# =============================================================================

mod_de_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h4(" Differential Analysis"),
    card(
      card_body(
        tags$div(style = "display:flex; align-items:center; gap:12px; margin-bottom:12px;",
          actionButton(ns("run"), "Run DE (design + contrasts + dpcDE)",
                       class = "btn-primary"),
          spinner_ui(ns("sp_de"))
        ),
        verbatimTextOutput(ns("de_log"))
      )
    ),

    h4(" Results Table"),
    card(
      card_body(
        fluidRow(
          column(3,
            selectInput(ns("view_layer"), "View",
                        choices  = c("Occupancy (full)"       = "occ",
                                     "Site (full)"            = "site",
                                     "Protein (full)"         = "prot",
                                     "Combined (all layers)"  = "combined"),
                        selected = "occ")),
          column(3,
            conditionalPanel(
              condition = sprintf("input['%s'] == 'combined'", ns("view_layer")),
              helpText("Shows the merged table with selected columns from all three layers.")
            )
          ),
          column(2, checkboxInput(ns("only_sig"),        "Significant only",  FALSE)),
          column(2, checkboxInput(ns("only_consistent"), "Consistent only",   FALSE)),
          column(2,
            helpText(tags$i(
              "'Consistent': direction holds across all reference lines."
            ))
          )
        ),
        DT::dataTableOutput(ns("res_table"))
      )
    ),

    h4(" Classification Counts"),
    card(
      card_body(verbatimTextOutput(ns("class_counts")))
    )
  )
}

mod_de_server <- function(id, norm_obj, cfg, status = reactiveVal('idle')) {
  moduleServer(id, function(input, output, session) {

    log_msg <- reactiveVal("")
    de_res  <- reactiveVal(NULL)
    sp_de_state <- reactiveVal("idle")
    spinner_server("sp_de", sp_de_state)

    append_log <- function(...) {
      msg <- paste0(paste(...), "\n")
      log_msg(paste0(log_msg(), msg))
    }

    observeEvent(input$run, {
      norm <- norm_obj();  c <- cfg()
      req(norm, c)
      log_msg(""); de_res(NULL)
      status("running")
      sp_de_state("running")

      withCallingHandlers(
        tryCatch({
          q.site <- norm$q.site
          q.occ  <- norm$q.occ
          q.prot <- norm$q.prot
          targets <- norm$y.sites$targets

          target_line  <- c$target_line
          ref_lines    <- c$reference_lines   # one line (vs_ctrl) or many (vs_rest)
          mode         <- c$comparison_mode

          # Lines that appear in the data (post-exclusion)
          cl           <- droplevels(targets$cell_line)
          lines_present <- levels(cl)

          # Validate: target and all reference lines must be in the data
          missing <- setdiff(c(target_line, ref_lines), lines_present)
          if (length(missing))
            stop(sprintf("Lines not found in data after exclusions: %s",
                         paste(missing, collapse = ", ")))

          prog <- make_progress_stepper(session, "Differential analysis…", c(
            "Design matrix"           = 0.10,
            "dpcDE (site layer)"      = 0.35,
            "dpcDE (occupancy layer)" = 0.60,
            "dpcDE (protein layer)"   = 0.75,
            "Results table"           = 0.88,
            "Consistency check"       = 1.00))
          on.exit(prog$close(), add = TRUE)

          # ---- Design matrix -----------------------------------------------
          append_log("=== Design matrix ===")
          prog$tick("Building design matrix & contrast…")
          design <- stats::model.matrix(~ 0 + cl)
          colnames(design) <- make.names(levels(cl))

          tgt <- make.names(target_line)
          rst <- make.names(ref_lines)

          # Contrast string differs by mode:
          #   vs_rest : target - (A + B + C) / 3
          #   vs_ctrl : target - control         (simple two-group)
          contr_str <- if (mode == "vs_ctrl") {
            sprintf("%s - %s", tgt, rst)
          } else {
            sprintf("%s - (%s)/%d", tgt, paste(rst, collapse = " + "), length(rst))
          }

          contr_name <- if (mode == "vs_ctrl") {
            sprintf("%s_vs_%s", make.names(target_line), make.names(ref_lines))
          } else {
            paste0(make.names(target_line), "_vs_rest")
          }

          cm <- limma::makeContrasts(contrasts = contr_str, levels = design)
          colnames(cm) <- contr_name

          append_log(sprintf("  Mode    : %s", mode))
          append_log(sprintf("  Contrast: %s", contr_str))
          append_log(sprintf("  Coef    : %s", contr_name))

          # ---- dpcDE -------------------------------------------------------
          append_log("\n=== dpcDE ===")
          run_de <- function(y, use_treat, lfc) {
            fit <- limpa::dpcDE(y, design, plot = FALSE)
            fit <- limma::contrasts.fit(fit, cm)
            if (use_treat) limma::treat(fit, lfc = lfc) else limma::eBayes(fit)
          }

          prog$tick("dpcDE: site layer…")
          fit.site <- run_de(q.site, c$use_treat, c$logfc_threshold)
          prog$tick("dpcDE: occupancy layer…")
          fit.occ  <- run_de(q.occ,  c$use_treat, c$logfc_threshold)
          prog$tick("dpcDE: protein layer…")
          fit.prot <- run_de(q.prot, c$use_treat, c$logfc_threshold)

          top_fn <- if (c$use_treat) limma::topTreat else limma::topTable

          tt.site <- top_fn(fit.site, coef = contr_name, number = Inf) |>
            tibble::rownames_to_column("site_id")
          tt.occ  <- top_fn(fit.occ,  coef = contr_name, number = Inf)
          if (!"site_id" %in% names(tt.occ))
            tt.occ <- tibble::rownames_to_column(tt.occ, "site_id")
          tt.prot <- top_fn(fit.prot, coef = contr_name, number = Inf) |>
            tibble::rownames_to_column("Protein.Group")

          append_log(sprintf("  Tested: %d sites | %d occ sites | %d proteins",
                             nrow(tt.site), nrow(tt.occ), nrow(tt.prot)))

          prog$tick("Assembling results table…")
          # ---- Combine results ---------------------------------------------
          fdr <- c$fdr_threshold

          res <- tt.site |>
            dplyr::select(site_id, site_label, Genes,
                          logFC_site = logFC, P_site  = P.Value,
                          FDR_site   = adj.P.Val,  t_site = t) |>
            dplyr::full_join(
              tt.occ |> dplyr::select(site_id,
                                      logFC_occ = logFC,  P_occ  = P.Value,
                                      FDR_occ   = adj.P.Val, t_occ = t),
              by = "site_id") |>
            dplyr::left_join(
              norm$site_info |>
                dplyr::select(site_id, Protein.Group, tier,
                              n_obs_target, n_obs_rest_max),
              by = "site_id") |>
            dplyr::left_join(
              tt.prot |>
                dplyr::select(Protein.Group,
                              logFC_prot = logFC, FDR_prot = adj.P.Val,
                              t_prot = t),
              by = "Protein.Group") |>
            dplyr::mutate(
              sig_site = !is.na(FDR_site) & FDR_site < fdr,
              sig_occ  = !is.na(FDR_occ)  & FDR_occ  < fdr,
              sig_prot = !is.na(FDR_prot) & FDR_prot < fdr,
              classification = dplyr::case_when(
                sig_occ &  sig_site & !sig_prot ~ "ubiquitination-specific",
                sig_occ &  sig_site &  sig_prot ~ "both layers change",
                sig_occ & !sig_site             ~ "occupancy only",
                !sig_occ & sig_site &  sig_prot ~ "abundance-driven",
                !sig_occ & sig_site             ~ "site change, not occupancy",
                TRUE                             ~ "ns"),
              direction = dplyr::case_when(
                is.na(logFC_occ) ~ NA_character_,
                logFC_occ > 0   ~ paste("up in", target_line),
                TRUE             ~ paste("down in", target_line))
            ) |>
            dplyr::arrange(FDR_occ, dplyr::desc(abs(logFC_occ)))

          # ---- Consistency check ------------------------------------------
          # vs_rest : direction consistent against every background line
          # vs_ctrl : direction consistent vs the single control
          #           (always TRUE by design, but min_abs_diff still informative)
          prog$tick("Consistency check across reference lines…")
          append_log("\n=== Consistency check ===")
          line_vec <- as.character(targets$cell_line)

          lm.occ <- sapply(lines_present, function(ln)
            rowMeans(q.occ$E[, line_vec == ln, drop = FALSE]))
          rownames(lm.occ) <- rownames(q.occ)

          ref_in_occ <- intersect(ref_lines, colnames(lm.occ))
          tgt_in_occ <- intersect(target_line, colnames(lm.occ))

          if (length(tgt_in_occ) == 0 || length(ref_in_occ) == 0) {
            append_log("  [!]  Cannot compute consistency -- target or reference absent from q.occ.")
            consistency <- dplyr::tibble(
              site_id     = rownames(q.occ),
              consistent  = NA,
              weakest_line= NA_character_,
              min_abs_diff= NA_real_,
              max_abs_diff= NA_real_)
          } else {
            d_ref <- sweep(
              -lm.occ[, ref_in_occ, drop = FALSE],
              1, lm.occ[, tgt_in_occ], "+")

            consistency <- dplyr::tibble(
              site_id          = rownames(q.occ),
              n_ref_same_sign  = pmax(rowSums(d_ref > 0), rowSums(d_ref < 0)),
              min_abs_diff     = apply(abs(d_ref), 1, min),
              max_abs_diff     = apply(abs(d_ref), 1, max)) |>
              dplyr::mutate(
                # For vs_ctrl with one reference the n_ref check is trivially 1/1
                consistent   = n_ref_same_sign == length(ref_in_occ) &
                               min_abs_diff    >= c$logfc_threshold / 2,
                weakest_line = ref_lines[apply(abs(d_ref), 1, which.min)])
          }

          res <- res |>
            dplyr::select(-dplyr::matches(
              "^(consistent|weakest_line|n_ref_same_sign|min_abs_diff|max_abs_diff)(\\.[xy])?$")) |>
            dplyr::left_join(
              consistency |> dplyr::select(site_id, consistent,
                                           weakest_line, min_abs_diff, max_abs_diff),
              by = "site_id")

          append_log(sprintf("  sig_occ + consistent: %d",
                             sum(res$sig_occ & res$consistent %in% TRUE, na.rm = TRUE)))
          append_log("\n[OK] DE analysis complete.")
          status("done")
          sp_de_state("done")

          de_res(list(
            res          = res,
            fit.site     = fit.site, fit.occ = fit.occ, fit.prot = fit.prot,
            tt.site      = tt.site,  tt.occ  = tt.occ,  tt.prot  = tt.prot,
            design       = design,   cm      = cm,
            contr_name   = contr_name,
            contr_str    = contr_str,
            mode         = mode,
            line_means_occ = lm.occ,
            lines_present  = lines_present
          ))
        }, error = function(e) {
          append_log(paste0("[ERR] ERROR: ", e$message)); status("error"); sp_de_state("error")
        }),
        message = function(m) { append_log(conditionMessage(m)); invokeRestart("muffleMessage") },
        warning = function(w) { append_log(paste("[!]", conditionMessage(w))); invokeRestart("muffleWarning") }
      )
    })

    # ---- outputs -----------------------------------------------------------
    output$de_log <- renderPrint({ cat(log_msg()) })

    res_display <- reactive({
      de <- de_res(); req(de)
      c  <- cfg()
      lyr <- input$view_layer

      # Pick the right table for the selected layer
      r <- switch(lyr,
        occ      = de$tt.occ,
        site     = de$tt.site,
        prot     = de$tt.prot,
        combined = de$res
      )
      req(r)

      # FDR column for filtering
      fdr_col <- switch(lyr,
        occ      = if ("adj.P.Val" %in% names(r)) "adj.P.Val" else "FDR_occ",
        site     = if ("adj.P.Val" %in% names(r)) "adj.P.Val" else "FDR_site",
        prot     = if ("adj.P.Val" %in% names(r)) "adj.P.Val" else "FDR_prot",
        combined = "FDR_occ"
      )

      if (input$only_sig && !is.null(fdr_col) && fdr_col %in% names(r))
        r <- r[!is.na(r[[fdr_col]]) & r[[fdr_col]] < c$fdr_threshold, ]

      if (input$only_consistent && "consistent" %in% names(r))
        r <- r[r$consistent %in% TRUE, ]

      r
    })

    output$res_table <- DT::renderDataTable({
      r   <- res_display()
      lyr <- input$view_layer

      # For combined view: show selected summary columns only
      if (lyr == "combined") {
        cols_show <- c("site_id", "site_label", "Genes",
                       "logFC_occ",  "P_occ",  "FDR_occ",
                       "logFC_site", "P_site", "FDR_site",
                       "logFC_prot", "P_prot", "FDR_prot",
                       "classification", "direction",
                       "consistent", "weakest_line", "tier")
        r <- r[, intersect(cols_show, names(r)), drop = FALSE]
      }
      # For layer views: show all topTable/topTreat columns
      # Round numeric columns for readability
      num_cols <- sapply(r, is.numeric)
      r[num_cols] <- lapply(r[num_cols], function(x) round(x, 5))
      r
    },
    options  = list(scrollX = TRUE, pageLength = 20),
    filter   = "top",
    rownames = FALSE)

    output$class_counts <- renderPrint({
      de <- de_res(); req(de)
      cat(sprintf("Contrast: %s\n\n", de$contr_str))
      print(dplyr::count(de$res, classification, sort = TRUE))
    })

    return(de_res)
  })
}
