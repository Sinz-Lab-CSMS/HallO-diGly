# =============================================================================
# mod_normalization -- Median normalization + DPC-Quant + PCA/PLS-DA
# Stops BEFORE occupancy correction (that lives in mod_occupancy.R).
#
# Returns:
#   $y.sites, $y.tp          -- normalised ELists
#   $q.site, $q.prot         -- dpcQuant outputs
#   $dpc.dg, $dpc.tp         -- DPC objects
#   $lines_present           -- cell lines in data
# =============================================================================

mod_normalization_ui <- function(id) {
  ns <- NS(id)
  tagList(

    # ── Run ─────────────────────────────────────────────────────────────────
    h4("Normalization & DPC-Quant"),
    card(
      card_body(
        tags$div(style = "display:flex; align-items:center; gap:12px; margin-bottom:12px;",
          actionButton(ns("run"), "Run normalization + DPC-Quant",
                       class = "btn-primary"),
          spinner_ui(ns("sp_norm"))
        ),
        verbatimTextOutput(ns("norm_log"))
      )
    ),

    # ── Intensity distributions ──────────────────────────────────────────────
    h4("Intensity Distributions (post-normalization)"),
    card(
      card_body(
        plotOutput(ns("boxplot"), height = "380px")
      )
    ),

    # ── DPC-Quant summary tables ─────────────────────────────────────────────
    h4("DPC-Quant Summary"),
    card(
      card_body(
        tags$p(class = "text-muted small mb-2",
          tags$b("NPrec:"), " precursors contributing to the site/protein estimate. ",
          tags$b("PropObs:"), " fraction of non-imputed values (0\u20131). ",
          tags$b("mean_SE:"), " mean quantification SE \u2014 high = high uncertainty."),
        fluidRow(
          column(6,
            tags$b("Sites (diGly layer)"),
            DT::dataTableOutput(ns("dpcquant_site_table"))
          ),
          column(6,
            tags$b("Protein groups (proteome)"),
            DT::dataTableOutput(ns("dpcquant_prot_table"))
          )
        )
      )
    ),

    # ── Multivariate QC: PCA / PLS-DA ───────────────────────────────────────
    h4("Multivariate QC \u2014 PCA / PLS-DA"),
    card(
      card_body(
        fluidRow(
          column(3,
            radioButtons(ns("mv_layer"), "Layer",
                         choices  = c("diGly sites" = "digly",
                                      "Proteome"    = "proteome"),
                         selected = "digly", inline = TRUE)
          ),
          column(3,
            radioButtons(ns("mv_method"), "Method",
                         choices  = c("PCA" = "pca", "PLS-DA" = "plsda"),
                         selected = "pca", inline = TRUE)
          ),
          column(2,
            numericInput(ns("mv_ncomp"), "Components", value = 3, min = 2, max = 10, step = 1)
          ),
          column(2,
            selectInput(ns("mv_xcomp"), "X axis", choices = NULL)
          ),
          column(2,
            selectInput(ns("mv_ycomp"), "Y axis", choices = NULL)
          )
        ),
        tags$div(style = "margin-bottom:10px;",
          actionButton(ns("run_mv"), "Run", class = "btn-outline-primary btn-sm"),
          spinner_ui(ns("sp_mv"))
        ),
        conditionalPanel(
          condition = sprintf("input['%s'] == 'plsda'", ns("mv_method")),
          helpText(icon("info-circle"),
            "PLS-DA uses cell_line as the grouping factor. ",
            "It is a supervised method \u2014 separation is expected by design. ",
            "Use PCA first as the unsupervised reference.")
        ),
        plotOutput(ns("mv_plot"), height = "480px"),
        verbatimTextOutput(ns("mv_summary"))
      )
    )
  )
}


mod_normalization_server <- function(id, filt_obj, cfg, status = reactiveVal("idle")) {
  moduleServer(id, function(input, output, session) {

    log_msg       <- reactiveVal("")
    norm_res      <- reactiveVal(NULL)
    sp_norm_state <- reactiveVal("idle")
    sp_mv_state   <- reactiveVal("idle")
    spinner_server("sp_norm", sp_norm_state)
    spinner_server("sp_mv",   sp_mv_state)

    append_log <- function(...) {
      msg <- paste0(paste(...), "\n")
      log_msg(paste0(log_msg(), msg))
    }

    median_centre <- function(y) {
      y$E <- sweep(y$E, 2, apply(y$E, 2, median, na.rm = TRUE), "-")
      y
    }

    # ── Main run ─────────────────────────────────────────────────────────────
    observeEvent(input$run, {
      filt <- filt_obj();  c <- cfg()
      req(filt, c)
      log_msg(""); norm_res(NULL)
      status("running"); sp_norm_state("running")

      withCallingHandlers(
        tryCatch({
          y.sites <- filt$y.sites
          y.tp    <- filt$y.tp
          dpc.dg  <- filt$dpc.dg
          dpc.tp  <- filt$dpc.tp

          prog <- make_progress_stepper(session, "Normalization & DPC-Quant\u2026", c(
            "Median normalisation"  = 0.15,
            "DPC-Quant (sites)"     = 0.55,
            "DPC-Quant (proteins)"  = 0.90,
            "Finalising"            = 1.00))
          on.exit(prog$close(), add = TRUE)

          # Normalization
          append_log("=== Normalization ===")
          prog$tick("Median centring intensities\u2026")
          if (c$digly_norm == "median") {
            y.sites <- median_centre(y.sites)
            dpc.dg  <- limpa::dpc(y.sites)
            append_log("  diGly: median-centred; DPC re-estimated.")
          } else {
            append_log("  diGly: no normalization.")
          }
          if (c$proteome_norm == "median") {
            y.tp   <- median_centre(y.tp)
            dpc.tp <- limpa::dpc(y.tp)
            append_log("  Proteome: median-centred; DPC re-estimated.")
          } else {
            append_log("  Proteome: no normalization.")
          }

          # DPC-Quant
          append_log("\n=== DPC-Quant ===")
          prog$tick("DPC-Quant: site-level quantification\u2026")
          q.site <- limpa::dpcQuant(y.sites, protein.id = "site_id",
                                    dpc = dpc.dg, verbose = FALSE)
          prog$tick("DPC-Quant: protein-level quantification\u2026")
          q.prot <- limpa::dpcQuant(y.tp, protein.id = "Protein.Group",
                                    dpc = dpc.tp, verbose = FALSE)
          q.site$targets <- q.prot$targets <- y.sites$targets
          append_log(sprintf("  Sites: %d | Protein groups: %d",
                             nrow(q.site), nrow(q.prot)))

          prog$tick("Finalising\u2026")
          append_log("\n[OK] Normalization + DPC-Quant complete.")
          status("done"); sp_norm_state("done")

          lines_pres <- levels(droplevels(y.sites$targets$cell_line))

          norm_res(list(
            y.sites       = y.sites,
            y.tp          = y.tp,
            q.site        = q.site,
            q.prot        = q.prot,
            dpc.dg        = dpc.dg,
            dpc.tp        = dpc.tp,
            lines_present = lines_pres
          ))

          # Update component dropdowns
          nc <- input$mv_ncomp
          comp_choices <- setNames(seq_len(nc), paste0("Comp ", seq_len(nc)))
          updateSelectInput(session, "mv_xcomp",
                            choices = comp_choices, selected = 1)
          updateSelectInput(session, "mv_ycomp",
                            choices = comp_choices, selected = 2)

        }, error = function(e) {
          append_log(paste0("[ERR] ERROR: ", e$message))
          status("error"); sp_norm_state("error")
        }),
        message = function(m) {
          append_log(conditionMessage(m)); invokeRestart("muffleMessage")
        },
        warning = function(w) {
          append_log(paste("[!]", conditionMessage(w))); invokeRestart("muffleWarning")
        }
      )
    })

    # Update component dropdowns when ncomp changes
    observeEvent(input$mv_ncomp, {
      nc <- max(2L, as.integer(input$mv_ncomp))
      comp_choices <- setNames(seq_len(nc), paste0("Comp ", seq_len(nc)))
      updateSelectInput(session, "mv_xcomp",
                        choices = comp_choices, selected = 1)
      updateSelectInput(session, "mv_ycomp",
                        choices = comp_choices, selected = min(2, nc))
    })

    # ── Multivariate QC ──────────────────────────────────────────────────────
    mv_res <- reactiveVal(NULL)

    observeEvent(input$run_mv, {
      res <- norm_res(); req(res)
      sp_mv_state("running")
      mv_res(NULL)

      tryCatch({
        # Select layer
        q_obj       <- if (input$mv_layer == "digly") res$q.site else res$q.prot
        layer_label <- if (input$mv_layer == "digly") "diGly sites" else "Proteome"
        X  <- t(q_obj$E)
        Y  <- droplevels(q_obj$targets$cell_line)
        nc <- max(2L, as.integer(input$mv_ncomp))

        # Remove features with any NA
        n_na <- if (anyNA(X)) sum(colSums(is.na(X)) > 0) else 0L
        if (n_na > 0)
          X <- X[, colSums(is.na(X)) == 0, drop = FALSE]

        # Remove zero-variance features (constant across all samples after
        # DPC-Quant imputation — common for sites seen in only one condition)
        feat_var <- apply(X, 2, var, na.rm = TRUE)
        n_zv     <- sum(feat_var == 0 | is.na(feat_var))
        if (n_zv > 0)
          X <- X[, feat_var > 0 & !is.na(feat_var), drop = FALSE]

        if (ncol(X) == 0) {
          showNotification(
            "No usable features after NA and zero-variance filtering.",
            type = "error", duration = 8)
          sp_mv_state("error"); return()
        }

        if (n_na > 0 || n_zv > 0)
          showNotification(
            sprintf("Filtered: %d NA-containing, %d zero-variance features removed. Using %d.",
                    n_na, n_zv, ncol(X)),
            type = "warning", duration = 6)

        # Cap ncomp at feasible maximum
        nc <- min(nc, nrow(X) - 1L, ncol(X))

        if (input$mv_method == "pca") {
          fit <- mixOmics::pca(X, ncomp = nc, center = TRUE, scale = TRUE)
        } else {
          fit <- mixOmics::plsda(X, Y, ncomp = nc, scale = TRUE)
        }
        mv_res(list(fit         = fit,
                    Y           = Y,
                    method      = input$mv_method,
                    layer       = input$mv_layer,
                    layer_label = layer_label,
                    xcomp       = as.integer(input$mv_xcomp),
                    ycomp       = as.integer(input$mv_ycomp),
                    lines       = res$lines_present,
                    n_features  = ncol(X)))
        sp_mv_state("done")

      }, error = function(e) {
        showNotification(paste("MV error:", e$message), type = "error")
        sp_mv_state("error")
      })
    })

    output$mv_plot <- renderPlot({
      r <- mv_res(); req(r)
      res  <- norm_res(); req(res)
      xc   <- r$xcomp
      yc   <- r$ycomp

      all_cl <- r$lines
      pal    <- grDevices::colorRampPalette(
        c("#e6ab02","#d95f02","#7570b3","#e7298a",
          "#1b9e77","#66a61e","#a6761d","#666666"))(length(all_cl))
      col_map <- setNames(pal, all_cl)

      q_used <- if (r$layer == "digly") res$q.site else res$q.prot
      mixOmics::plotIndiv(
        r$fit,
        comp          = c(xc, yc),
        group         = r$Y,
        ind.names     = colnames(q_used),
        legend        = TRUE,
        ellipse       = FALSE,
        col.per.group = col_map[levels(r$Y)],
        title         = sprintf("%s \u2014 %s (Comp %d vs %d)",
                                toupper(r$method), r$layer_label, xc, yc)
      )
    })

    output$mv_summary <- renderPrint({
      r <- mv_res(); req(r)
      cat(sprintf("Layer    : %s\n", r$layer_label))
      cat(sprintf("Method   : %s\n", toupper(r$method)))
      cat(sprintf("Features : %d\n", r$n_features))
      cat(sprintf("Samples  : %d\n", length(r$Y)))
      cat("\n")
      if (r$method == "pca") {
        expl <- r$fit$prop_expl_var$X
        nc   <- min(length(expl), 5L)
        cat("Explained variance (first 5 components):\n")
        df <- data.frame(
          Component  = paste0("PC", seq_len(nc)),
          Individual = round(expl[seq_len(nc)], 3),
          Cumulative = round(cumsum(expl)[seq_len(nc)], 3))
        print(df, row.names = FALSE)
      } else {
        cat("PLS-DA: supervised method.\n")
        cat(sprintf("Components: %d\n", r$fit$ncomp))
        cat("Grouping factor: cell_line\n")
      }
    })

    # ── Static outputs ────────────────────────────────────────────────────────
    output$norm_log <- renderPrint({ cat(log_msg()) })

    output$boxplot <- renderPlot({
      res <- norm_res(); req(res)
      targets <- res$y.sites$targets
      all_cl  <- levels(factor(as.character(targets$cell_line)))
      pal     <- grDevices::colorRampPalette(
        c("#e6ab02","#d95f02","#7570b3","#e7298a",
          "#1b9e77","#66a61e","#a6761d","#666666"))(length(all_cl))
      line_colors <- setNames(pal, all_cl)
      meta_df <- data.frame(name      = rownames(targets),
                            cell_line = as.character(targets$cell_line))

      dplyr::bind_rows(
        as.data.frame(res$y.sites$E) |>
          tibble::rownames_to_column("prec") |>
          tidyr::pivot_longer(-prec, names_to = "name") |>
          dplyr::mutate(layer = "diGly"),
        as.data.frame(res$y.tp$E) |>
          tibble::rownames_to_column("prec") |>
          tidyr::pivot_longer(-prec, names_to = "name") |>
          dplyr::mutate(layer = "Proteome")
      ) |>
        dplyr::left_join(meta_df, by = "name") |>
        ggplot2::ggplot(ggplot2::aes(name, value, fill = cell_line)) +
        ggplot2::geom_boxplot(outlier.size = 0.2, na.rm = TRUE) +
        ggplot2::facet_wrap(~ layer, scales = "free_y", ncol = 1) +
        ggplot2::scale_fill_manual(values = line_colors) +
        ggplot2::labs(x = NULL, y = "log2 intensity (normalised)") +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    })

    # DPC-Quant tables helper
    dpcquant_summary <- function(q) {
      g       <- as.data.frame(q$genes)
      mean_se <- if (!is.null(q$other$standard.error))
                   rowMeans(q$other$standard.error, na.rm = TRUE) else NA_real_
      df <- dplyr::tibble(
        id      = rownames(q),
        NPrec   = if ("NPrec"   %in% names(g)) g$NPrec   else NA_integer_,
        PropObs = if ("PropObs" %in% names(g)) round(g$PropObs, 3) else NA_real_,
        mean_SE = round(mean_se, 4)
      )
      if ("site_label" %in% names(g))
        df <- dplyr::mutate(df, label = g$site_label, .after = id)
      else if ("Genes" %in% names(g))
        df <- dplyr::mutate(df, label = g$Genes, .after = id)
      df
    }

    output$dpcquant_site_table <- DT::renderDataTable({
      res <- norm_res(); req(res)
      dpcquant_summary(res$q.site)
    }, options = list(scrollX = TRUE, pageLength = 10), rownames = FALSE)

    output$dpcquant_prot_table <- DT::renderDataTable({
      res <- norm_res(); req(res)
      dpcquant_summary(res$q.prot)
    }, options = list(scrollX = TRUE, pageLength = 10), rownames = FALSE)

    return(norm_res)
  })
}
