# =============================================================================
# mod_export -- Save results, matrices, RDS objects, and figures (SVG)
# Figures: volcano, DPC curves, PCA/PLS-DA, voom plot
# =============================================================================

mod_export_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h4("Export Results"),
    card(
      card_body(
        tags$p(class = "text-muted small mb-3",
          "Outputs write to the Results / Figures / limpa directories set in the Data tab."),

        fluidRow(
          # --- TSV / RDS -------------------------------------------------------
          column(4,
            h5("TSV tables"),
            checkboxInput(ns("exp_res"),       "Results table (site + occ + protein combined)", TRUE),
            tags$hr(style="margin:4px 0;"),
            tags$small(class="text-muted", "Full DE tables per layer:"),
            checkboxInput(ns("exp_de_occ"),   "DE — Occupancy layer (full)",  TRUE),
            checkboxInput(ns("exp_de_site"),  "DE — Site layer (full)",       FALSE),
            checkboxInput(ns("exp_de_prot"),  "DE — Protein layer (full)",    FALSE),
            tags$hr(style="margin:4px 0;"),
            checkboxInput(ns("exp_occ_mat"), "Occupancy log2 matrix", TRUE),
            checkboxInput(ns("exp_site_mat"),"Site log2 matrix", TRUE),
            checkboxInput(ns("exp_meta"),    "Sample metadata", TRUE),
            hr(),
            h5("R objects"),
            checkboxInput(ns("exp_rds"),     "Full analysis RDS", TRUE)
          ),

          # --- Figures ---------------------------------------------------------
          column(5,
            h5("Figures (SVG)"),
            checkboxInput(ns("exp_volcano"),  "Volcano plot", TRUE),
            checkboxInput(ns("exp_dpc"),      "DPC curves (diGly + proteome)",   TRUE),
            checkboxInput(ns("exp_mv"),       "PCA / PLS-DA scores plot",        TRUE),
            tags$label("vooma mean-variance trend", class = "form-label mb-1 mt-1"),
            radioButtons(ns("exp_voom"), label = NULL,
                         choices = c(
                           "Site layer only"         = "site",
                           "Occupancy layer only"    = "occ",
                           "Site + Occupancy (both)" = "both",
                           "Do not export"           = "none"),
                         selected = "both"),
            checkboxInput(ns("exp_heatmap"),  "Heatmap PNG (top consistent)",    FALSE),
            hr(),
            h5("PCA / PLS-DA settings for export"),
            fluidRow(
              column(6,
                radioButtons(ns("fig_mv_layer"), "Layer",
                             choices  = c("diGly" = "digly", "Proteome" = "proteome"),
                             selected = "digly", inline = TRUE)
              ),
              column(6,
                radioButtons(ns("fig_mv_method"), "Method",
                             choices  = c("PCA" = "pca", "PLS-DA" = "plsda"),
                             selected = "pca", inline = TRUE)
              )
            ),
            numericInput(ns("fig_mv_ncomp"), "Components", 3, 2, 10, 1)
          ),
          tags$hr(),
          h5("Volcano settings for export"),
          fluidRow(
            column(4,
              selectInput(ns("fig_v_layer"), "Layer",
                          choices  = c("Occupancy" = "occ",
                                       "Site"      = "site",
                                       "Protein"   = "prot"),
                          selected = "occ")
            ),
            column(4,
              numericInput(ns("fig_v_label_n"), "Max labels per direction",
                           30, 5, 200, 5)
            ),
            column(4,
              checkboxInput(ns("fig_v_consistent"), "Consistent only", FALSE)
            )
          ),

          # --- Run -------------------------------------------------------------
          column(3,
            tags$div(style = "margin-top: 28px;",
              actionButton(ns("export"), "Export selected outputs",
                           class = "btn-success btn-lg w-100"),
              spinner_ui(ns("sp_export"))
            ),
            verbatimTextOutput(ns("export_log"))
          )
        )
      )
    ),

    # --- Browser downloads ---------------------------------------------------
    h4("Browser downloads"),
    card(
      card_body(
        fluidRow(
          column(3, downloadButton(ns("dl_res"),      "Results TSV")),
          column(3, downloadButton(ns("dl_occ"),      "Occupancy matrix TSV")),
          column(3, downloadButton(ns("dl_site"),     "Site matrix TSV"))
        ),
        fluidRow(
          column(3, downloadButton(ns("dl_de_occ"),  "DE occupancy TSV")),
          column(3, downloadButton(ns("dl_de_site"), "DE site TSV")),
          column(3, downloadButton(ns("dl_de_prot"), "DE protein TSV"))
        )
      )
    )
  )
}


mod_export_server <- function(id, de_obj, norm_obj, filt_obj,
                               cfg, status = reactiveVal("idle")) {
  moduleServer(id, function(input, output, session) {

    log_msg      <- reactiveVal("")
    sp_exp_state <- reactiveVal("idle")
    spinner_server("sp_export", sp_exp_state)

    append_log <- function(...) {
      msg <- paste0(paste(...), "\n")
      log_msg(paste0(log_msg(), msg))
    }

    safe_tsv <- function(df, path) {
      tryCatch({
        readr::write_tsv(df, path)
        append_log(sprintf("  [OK] %s", basename(path)))
      }, error = function(e)
        append_log(sprintf("  [ERR] %s: %s", basename(path), e$message)))
    }

    save_svg <- function(plot_expr, path, width = 9, height = 7) {
      tryCatch({
        svglite::svglite(path, width = width, height = height)
        force(plot_expr)
        grDevices::dev.off()
        append_log(sprintf("  [OK] %s", basename(path)))
      }, error = function(e) {
        if (grDevices::dev.cur() > 1) grDevices::dev.off()
        append_log(sprintf("  [ERR] %s: %s", basename(path), e$message))
      })
    }

    # ── Line-color palette helper ─────────────────────────────────────────────
    make_palette <- function(lines) {
      pal <- grDevices::colorRampPalette(
        c("#e6ab02","#d95f02","#7570b3","#e7298a",
          "#1b9e77","#66a61e","#a6761d","#666666"))(length(lines))
      setNames(pal, lines)
    }

    # ── Export to disk ────────────────────────────────────────────────────────
    observeEvent(input$export, {
      de   <- de_obj();   norm <- norm_obj()
      filt <- filt_obj(); c    <- cfg()
      req(de, norm, c)
      log_msg(""); status("running"); sp_exp_state("running")

      tgt <- make.names(c$target_line)

      for (d in c(c$results_dir, c$figures_dir, c$limpa_dir))
        if (nchar(d) > 0) dir.create(d, showWarnings = FALSE, recursive = TRUE)

      append_log("=== TSV tables ===")

      # Results table
      if (input$exp_res && nchar(c$results_dir) > 0)
        safe_tsv(de$res,
                 file.path(c$results_dir,
                           sprintf("diGly_%s_vs_rest_site_and_occupancy.tsv", tgt)))

      # Occupancy matrix
      # --- Full DE tables per layer -------------------------------------------
      if (input$exp_de_occ && nchar(c$results_dir) > 0 && !is.null(de$tt.occ))
        safe_tsv(de$tt.occ,
                 file.path(c$results_dir,
                           sprintf("diGly_%s_DE_occupancy_full.tsv", tgt)))

      if (input$exp_de_site && nchar(c$results_dir) > 0 && !is.null(de$tt.site))
        safe_tsv(de$tt.site,
                 file.path(c$results_dir,
                           sprintf("diGly_%s_DE_site_full.tsv", tgt)))

      if (input$exp_de_prot && nchar(c$results_dir) > 0 && !is.null(de$tt.prot))
        safe_tsv(de$tt.prot,
                 file.path(c$results_dir,
                           sprintf("diGly_%s_DE_protein_full.tsv", tgt)))

      if (input$exp_occ_mat && nchar(c$results_dir) > 0) {
        mat <- as.data.frame(norm$q.occ$E) |>
          tibble::rownames_to_column("site_id") |>
          dplyr::left_join(norm$site_info, by = "site_id")
        safe_tsv(mat,
                 file.path(c$results_dir,
                           sprintf("diGly_%s_occupancy_log2_matrix.tsv", tgt)))
      }

      # Site matrix
      if (input$exp_site_mat && nchar(c$results_dir) > 0)
        safe_tsv(
          as.data.frame(norm$q.site$E) |> tibble::rownames_to_column("site_id"),
          file.path(c$results_dir,
                    sprintf("diGly_%s_site_log2_matrix.tsv", tgt)))

      # Metadata
      if (input$exp_meta && nchar(c$results_dir) > 0)
        safe_tsv(
          as.data.frame(norm$y.sites$targets) |> tibble::rownames_to_column("sample_id"),
          file.path(c$results_dir, "samples_used.tsv"))

      # RDS
      append_log("\n=== R objects ===")
      if (input$exp_rds && nchar(c$limpa_dir) > 0) {
        path <- file.path(c$limpa_dir,
                          sprintf("limpa_occupancy_%s_vs_rest.rds", tgt))
        tryCatch({
          saveRDS(list(
            # Quantification
            q.site        = norm$q.site,
            q.prot        = norm$q.prot,
            q.occ         = norm$q.occ,
            dpc.dg        = norm$dpc.dg,
            dpc.tp        = norm$dpc.tp,
            # Normalised ELists
            y.sites       = norm$y.sites,
            y.tp          = norm$y.tp,
            # DE
            fit.site      = de$fit.site,
            fit.prot      = de$fit.prot,
            fit.occ       = de$fit.occ,
            tt.site       = de$tt.site,
            tt.prot       = de$tt.prot,
            tt.occ        = de$tt.occ,
            design        = de$design,
            contrasts     = de$cm,
            res           = de$res,
            # Occupancy metadata
            site_info     = norm$site_info,
            link          = norm$link,
            lines_present = norm$lines_present,
            # Config
            cfg           = c,
            # Version tag for compatibility checks
            hallo_version = "2.0"), path)
          append_log(sprintf("  [OK] %s", basename(path)))
        }, error = function(e) append_log(sprintf("  [ERR] RDS: %s", e$message)))
      }

      append_log("\n=== Figures (SVG) ===")

      # ── Volcano ─────────────────────────────────────────────────────────────
      if (input$exp_volcano && nchar(c$figures_dir) > 0) {
        lyr  <- input$fig_v_layer
        path <- file.path(c$figures_dir,
                          sprintf("DE_%s_%s_vs_rest.svg", lyr, tgt))
        p <- make_volcano(
          data             = de$res,
          layer            = lyr,
          fdr_thresh       = c$fdr_threshold,
          logfc_thresh     = c$logfc_threshold,
          selected_genes   = c$selected_genes,
          target_line      = c$target_line,
          require_consistent = input$fig_v_consistent,
          label_n          = input$fig_v_label_n)
        if (!is.null(p)) save_svg(print(p), path)
      }

      # ── DPC curves (diGly + proteome) ────────────────────────────────────────
      if (input$exp_dpc && !is.null(filt) && nchar(c$figures_dir) > 0) {
        path <- file.path(c$figures_dir,
                          sprintf("DPC_curves_%s.svg", tgt))
        save_svg({
          par(mfrow = c(1, 2))
          limpa::plotDPC(filt$dpc.dg, main = "DPC — diGly layer")
          limpa::plotDPC(filt$dpc.tp, main = "DPC — Proteome layer")
          par(mfrow = c(1, 1))
        }, path, width = 12, height = 6)
      } else if (input$exp_dpc && is.null(filt)) {
        append_log("  [SKIP] DPC curves: run Filtering first.")
      }

      # ── PCA / PLS-DA ─────────────────────────────────────────────────────────
      if (input$exp_mv && nchar(c$figures_dir) > 0) {
        path <- file.path(c$figures_dir,
                          sprintf("%s_%s_%s.svg",
                                  toupper(input$fig_mv_method),
                                  input$fig_mv_layer, tgt))
        tryCatch({
          q_obj       <- if (input$fig_mv_layer == "digly") norm$q.site else norm$q.prot
          layer_label <- if (input$fig_mv_layer == "digly") "diGly sites" else "Proteome"
          X  <- t(q_obj$E)
          Y  <- droplevels(q_obj$targets$cell_line)
          nc <- max(2L, as.integer(input$fig_mv_ncomp))

          # NA + zero-variance filter
          if (anyNA(X)) X <- X[, colSums(is.na(X)) == 0, drop = FALSE]
          feat_var <- apply(X, 2, var, na.rm = TRUE)
          X <- X[, feat_var > 0 & !is.na(feat_var), drop = FALSE]
          nc <- min(nc, nrow(X) - 1L, ncol(X))

          fit <- if (input$fig_mv_method == "pca")
            mixOmics::pca(X, ncomp = nc, center = TRUE, scale = TRUE)
          else
            mixOmics::plsda(X, Y, ncomp = nc, scale = TRUE)

          lines   <- levels(Y)
          col_map <- make_palette(lines)

          save_svg({
            mixOmics::plotIndiv(
              fit,
              comp          = c(1, 2),
              group         = Y,
              ind.names     = rownames(X),
              legend        = TRUE,
              ellipse       = FALSE,
              col.per.group = col_map[levels(Y)],
              title         = sprintf("%s — %s (%s vs rest)",
                                      toupper(input$fig_mv_method),
                                      layer_label, c$target_line))
          }, path, width = 8, height = 7)

        }, error = function(e)
          append_log(sprintf("  [ERR] PCA/PLS-DA: %s", e$message)))
      }

      # ── vooma mean-variance trend ─────────────────────────────────────────────
      if (input$exp_voom != "none" && nchar(c$figures_dir) > 0) {
        tryCatch({
          req(de$design)
          # Build the list of layers to export based on user selection
          layers <- switch(input$exp_voom,
            site = { l <- list(); if (!is.null(norm$q.site)) l[["site"]] <- norm$q.site; l },
            occ  = { l <- list(); if (!is.null(norm$q.occ))  l[["occupancy"]] <- norm$q.occ; l },
            both = {
              l <- list()
              if (!is.null(norm$q.site)) l[["site"]]      <- norm$q.site
              if (!is.null(norm$q.occ))  l[["occupancy"]] <- norm$q.occ
              l
            }
          )
          if (length(layers) == 0)
            stop("Selected layer(s) not available — run Occupancy tab first.")

          lyr_tag <- switch(input$exp_voom, site = "site", occ = "occ", both = "site_occ")
          path    <- file.path(c$figures_dir,
                               sprintf("vooma_meanvar_%s_%s.svg", lyr_tag, tgt))
          save_svg({
            par(mfrow = c(1, length(layers)))
            for (nm in names(layers)) {
              limpa::dpcDE(layers[[nm]], de$design, plot = TRUE)
              title(main = paste("vooma —", nm, "layer"), line = 0.5)
            }
            par(mfrow = c(1, 1))
          }, path, width = 5 * length(layers), height = 6)
        }, error = function(e)
          append_log(sprintf("  [ERR] vooma plot: %s", e$message)))
      }

      # ── Heatmap PNG ──────────────────────────────────────────────────────────
      if (input$exp_heatmap && nchar(c$figures_dir) > 0) {
        path <- file.path(c$figures_dir,
                          sprintf("heatmap_%s.png", tgt))
        tryCatch({
          top <- de$res |>
            dplyr::filter(!is.na(FDR_occ), consistent == TRUE) |>
            dplyr::slice_min(FDR_occ, n = 50)
          mat <- norm$q.occ$E[top$site_id[top$site_id %in% rownames(norm$q.occ$E)],
                               , drop = FALSE]
          if (nrow(mat) < 2) stop("Fewer than 2 consistent sites for heatmap.")
          grDevices::png(path, width = 2000, height = 2400, res = 300)
          pheatmap::pheatmap(mat, scale = "row",
                             main = sprintf("Consistent sites — %s vs rest", c$target_line))
          grDevices::dev.off()
          append_log(sprintf("  [OK] %s", basename(path)))
        }, error = function(e) {
          if (grDevices::dev.cur() > 1) grDevices::dev.off()
          append_log(sprintf("  [ERR] Heatmap: %s", e$message))
        })
      }

      append_log("\nExport done.")
      status("done"); sp_exp_state("done")
    })

    output$export_log <- renderPrint({ cat(log_msg()) })

    # ── Browser downloads ─────────────────────────────────────────────────────
    output$dl_res <- downloadHandler(
      filename = function()
        sprintf("diGly_%s_vs_rest_results.tsv", make.names(cfg()$target_line)),
      content = function(file) {
        de <- de_obj(); req(de); readr::write_tsv(de$res, file)
      }
    )

    output$dl_occ <- downloadHandler(
      filename = function()
        sprintf("diGly_%s_occupancy_matrix.tsv", make.names(cfg()$target_line)),
      content = function(file) {
        norm <- norm_obj(); req(norm)
        as.data.frame(norm$q.occ$E) |>
          tibble::rownames_to_column("site_id") |>
          dplyr::left_join(norm$site_info, by = "site_id") |>
          readr::write_tsv(file)
      }
    )

    output$dl_site <- downloadHandler(
      filename = function()
        sprintf("diGly_%s_site_matrix.tsv", make.names(cfg()$target_line)),
      content = function(file) {
        norm <- norm_obj(); req(norm)
        as.data.frame(norm$q.site$E) |>
          tibble::rownames_to_column("site_id") |>
          readr::write_tsv(file)
      }
    )

    output$dl_de_occ <- downloadHandler(
      filename = function()
        sprintf("diGly_%s_DE_occupancy_full.tsv", make.names(cfg()$target_line)),
      content = function(file) {
        de <- de_obj(); req(de, de$tt.occ)
        readr::write_tsv(de$tt.occ, file)
      }
    )

    output$dl_de_site <- downloadHandler(
      filename = function()
        sprintf("diGly_%s_DE_site_full.tsv", make.names(cfg()$target_line)),
      content = function(file) {
        de <- de_obj(); req(de, de$tt.site)
        readr::write_tsv(de$tt.site, file)
      }
    )

    output$dl_de_prot <- downloadHandler(
      filename = function()
        sprintf("diGly_%s_DE_protein_full.tsv", make.names(cfg()$target_line)),
      content = function(file) {
        de <- de_obj(); req(de, de$tt.prot)
        readr::write_tsv(de$tt.prot, file)
      }
    )
  })
}


