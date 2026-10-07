# =============================================================================
# mod_plots -- Volcano, heatmap, acquisition QC
# make_volcano is top-level so mod_export.R can call it directly.
# =============================================================================

# ---- Top-level volcano helper (used by both plots tab and export) -----------
make_volcano <- function(data, layer, fdr_thresh, logfc_thresh,
                         selected_genes, target_line,
                         require_consistent, label_n) {
  fdr_col   <- paste0("FDR_", layer)
  logfc_col <- paste0("logFC_", layer)
  lab_col   <- if (layer == "prot") "Genes" else "site_label"
  gene_col  <- "Genes"

  if (!all(c(fdr_col, logfc_col) %in% names(data))) {
    return(ggplot2::ggplot() +
             ggplot2::annotate("text", x=0, y=0, label="Layer not available") +
             ggplot2::theme_void())
  }

  pd <- data |>
    dplyr::mutate(
      .fdr   = .data[[fdr_col]],
      .logfc = .data[[logfc_col]],
      .lab   = as.character(.data[[lab_col]]),
      .gene  = as.character(.data[[gene_col]])) |>
    dplyr::filter(!is.na(.fdr), !is.na(.logfc))

  if (nrow(pd) == 0) return(NULL)

  if ("site_id" %in% names(pd))
    pd <- pd |> dplyr::mutate(.lab = ifelse(is.na(.lab) | .lab == "", site_id, .lab))

  ok_consistent <- if (!require_consistent || !"consistent" %in% names(pd)) {
    TRUE
  } else { pd$consistent %in% TRUE }

  pd <- pd |>
    dplyr::mutate(
      neglog10FDR = -log10(.fdr + 1e-300),
      .passes     = .fdr < fdr_thresh & abs(.logfc) >= logfc_thresh & ok_consistent,
      dir = dplyr::case_when(.passes & .logfc > 0 ~ "UP",
                             .passes & .logfc < 0 ~ "DN",
                             TRUE                 ~ "NS"),
      is_sel  = .gene %in% selected_genes,
      category = dplyr::case_when(
        is_sel              ~ "Selected",
        dir == "UP"         ~ "Up",
        dir == "DN"         ~ "Down",
        TRUE                ~ "NS"),
      category = factor(category, levels = c("Up","Down","NS","Selected")))

  n_up <- sum(pd$dir == "UP")
  n_dn <- sum(pd$dir == "DN")

  top_labs <- pd |>
    dplyr::filter(dir != "NS") |>
    dplyr::group_by(dir) |>
    dplyr::arrange(.fdr, dplyr::desc(abs(.logfc)), .by_group = TRUE) |>
    dplyr::slice_head(n = label_n) |>
    dplyr::pull(.lab)

  pd <- pd |>
    dplyr::mutate(label = ifelse(.lab %in% top_labs | is_sel, .lab, ""))

  col_map <- c("Up" = "#e66101", "Down" = "#5e3c99",
               "NS" = "grey70",   "Selected" = "#1b9e77")

  ggplot2::ggplot(pd, ggplot2::aes(.logfc, neglog10FDR, color = category)) +
    ggplot2::geom_point(
      ggplot2::aes(alpha = ifelse(category == "NS", 0.3, 0.6),
                   size  = ifelse(category %in% c("NS","Up","Down"), 1.2, 2.5))) +
    ggplot2::geom_hline(yintercept = -log10(fdr_thresh),
                        linetype = "dashed", linewidth = 0.5) +
    ggplot2::geom_vline(xintercept = c(-logfc_thresh, logfc_thresh),
                        linetype = "dashed", linewidth = 0.5) +
    ggrepel::geom_text_repel(
      ggplot2::aes(label = label),
      size = 3, max.overlaps = 20,
      box.padding = 0.4, segment.size = 0.3, na.rm = TRUE) +
    ggplot2::scale_color_manual(values = col_map, drop = FALSE) +
    ggplot2::scale_alpha_identity() +
    ggplot2::scale_size_identity() +
    ggplot2::annotate("text", x = -Inf, y = Inf,
                      label = sprintf("UP: %d  DN: %d", n_up, n_dn),
                      hjust = -0.1, vjust = 1.5, size = 4.5, fontface = "bold") +
    ggplot2::labs(
      title    = sprintf("%s -- %s layer (FDR < %.2f, |logFC| > %.1f)",
                         target_line, toupper(layer), fdr_thresh, logfc_thresh),
      subtitle = sprintf("%s tested%s",
                         format(nrow(pd), big.mark = ","),
                         if (require_consistent) " | consistent only" else ""),
      x = expression(log[2]~"FC"),
      y = expression(-log[10]~"(FDR)")) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
}

# =============================================================================
# UI
# =============================================================================
mod_plots_ui <- function(id) {
  ns <- NS(id)
  tagList(

    # ---- Volcano ------------------------------------------------------------
    h4(" Volcano Plot"),
    card(
      card_body(
        fluidRow(
          column(3, selectInput(ns("v_layer"), "Layer",
                                choices  = c("Occupancy" = "occ",
                                             "Site (uncorrected)" = "site",
                                             "Protein" = "prot"),
                                selected = "occ")),
          column(2, checkboxInput(ns("v_consistent"), "Consistent only", FALSE)),
          column(2, numericInput(ns("v_label_n"), "Max labels per direction", 30, 5, 200, 5)),
          column(3, actionButton(ns("draw_volcano"), "Draw volcano",
                                 class = "btn-primary mt-4"))
        ),
        plotOutput(ns("volcano"), height = "550px")
      )
    ),

    # ---- Heatmap ------------------------------------------------------------
    h4(" Heatmap (top sig sites)"),
    card(
      card_body(
        fluidRow(
          column(3, numericInput(ns("h_n"), "Top N sites", 50, 5, 500, 5)),
          column(3, checkboxInput(ns("h_consistent"), "Consistent only", FALSE)),
          column(3, actionButton(ns("draw_heatmap"), "Draw heatmap",
                                 class = "btn-primary mt-4"))
        ),
        plotOutput(ns("heatmap"), height = "600px")
      )
    ),

    # ---- Acquisition / depth QC ------------------------------------------
    h4(" Acquisition & Depth QC"),
    card(
      card_body(
        fluidRow(
          column(6, plotOutput(ns("depth_qc"),  height = "280px")),
          column(6, plotOutput(ns("acq_order"), height = "280px"))
        )
      )
    )
  )
}

# =============================================================================
# Server
# =============================================================================
mod_plots_server <- function(id, de_obj, norm_obj, cfg) {
  moduleServer(id, function(input, output, session) {

    make_line_colors <- function(cell_lines) {
      cl  <- unique(as.character(cell_lines))
      pal <- grDevices::colorRampPalette(
               c("#e6ab02","#d95f02","#7570b3","#e7298a",
                 "#1b9e77","#66a61e","#a6761d","#666666"))(length(cl))
      setNames(pal, cl)
    }
    color_palette <- grDevices::colorRampPalette(
      c("#5e3c99", "#b2abd2", "#f7f7f7", "#fdb863", "#e66101"))(100)

    # ---- Volcano -----------------------------------------------------------
    observeEvent(input$draw_volcano, {
      de   <- de_obj();  c <- cfg()
      req(de, c)

      output$volcano <- renderPlot({
        make_volcano(
          data               = de$res,
          layer              = input$v_layer,
          fdr_thresh         = c$fdr_threshold,
          logfc_thresh       = c$logfc_threshold,
          selected_genes     = c$selected_genes,
          target_line        = c$target_line,
          require_consistent = input$v_consistent,
          label_n            = input$v_label_n)
      })
    })

    # ---- Heatmap -----------------------------------------------------------
    observeEvent(input$draw_heatmap, {
      de   <- de_obj();  norm <- norm_obj();  c <- cfg()
      req(de, norm, c)

      output$heatmap <- renderPlot({
        r <- de$res
        if (input$h_consistent) r <- r[r$consistent %in% TRUE, ]

        top_sites <- r |>
          dplyr::filter(sig_occ) |>
          dplyr::slice_head(n = input$h_n) |>
          dplyr::pull(site_id)

        if (length(top_sites) < 2) {
          plot.new(); text(0.5, 0.5, "Fewer than 2 sig sites for heatmap"); return()
        }

        M <- norm$q.occ$E[top_sites, , drop = FALSE]
        M <- t(scale(t(M)))

        rn_map <- setNames(de$res$site_label, de$res$site_id)
        rownames(M) <- ifelse(is.na(rn_map[top_sites]) | rn_map[top_sites] == "",
                              top_sites, rn_map[top_sites])

        targets <- norm$y.sites$targets
        ann_col <- data.frame(
          cell_line = as.character(targets$cell_line),
          row.names = rownames(targets))

        line_colors <- make_line_colors(targets$cell_line)
        ann_colors  <- list(
          cell_line = line_colors[intersect(names(line_colors),
                                             unique(ann_col$cell_line))])

        pheatmap::pheatmap(
          M, color = color_palette,
          annotation_col    = ann_col,
          annotation_colors = ann_colors,
          show_colnames     = TRUE,
          fontsize_row      = max(5, 9 - floor(nrow(M) / 20)),
          main = sprintf("Top %d occupancy sites -- %s vs rest (z-scored)",
                         length(top_sites), c$target_line))
      })
    })

    # ---- Acquisition QC plots (always rendered when data available) --------
    output$depth_qc <- renderPlot({
      norm <- norm_obj(); req(norm)
      targets <- norm$y.sites$targets
      meta_df <- data.frame(
        sample_id = rownames(targets),
        cell_line = as.character(targets$cell_line),
        inj_order = targets$inj_order,
        n_prec    = colSums(!is.na(norm$y.sites$E)))

      ggplot2::ggplot(meta_df, ggplot2::aes(inj_order, n_prec, color = cell_line)) +
        ggplot2::geom_line(ggplot2::aes(group = 1), color = "grey60", linewidth = 0.3) +
        ggplot2::geom_point(size = 3) +
        ggplot2::scale_color_manual(values = make_line_colors(targets$cell_line)) +
        ggplot2::labs(x = "Injection order", y = "diGly precursors ID'd",
                      title = "Identification depth") +
        ggplot2::theme_bw()
    })

    output$acq_order <- renderPlot({
      norm <- norm_obj(); req(norm)
      targets <- norm$y.sites$targets
      meta_df <- data.frame(
        sample_id = rownames(targets),
        cell_line = as.character(targets$cell_line),
        inj_order = targets$inj_order)

      ggplot2::ggplot(meta_df,
                      ggplot2::aes(inj_order, cell_line, color = cell_line)) +
        ggplot2::geom_point(size = 4) +
        ggplot2::scale_color_manual(values = make_line_colors(targets$cell_line)) +
        ggplot2::labs(x = "Injection order", y = NULL,
                      title = "Acquisition order by cell line") +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "none")
    })
  })
}
