# =============================================================================
# mod_rnk -- GSEA .rnk rank-file generator
#
# Builds gene-level rank files (GSEA Desktop / GSEA_4.4.0-style .rnk, two
# columns, no header) from the DE results already computed in the app
# (Occupancy / Site / Protein layers), with the same ranking methods used in
# the standalone 03a_prepare_rnk_files_proteomics.R script:
#
#   logFC             -- rank by logFC
#   logFC_pvalue      -- rank by logFC * -log10(pValue)
#   t_stat            -- rank by moderated t-statistic (one row per feature)
#   t_stat_stouffer   -- gene-level Stouffer-combined t: sum(t) / sqrt(n_sites)
#                        (recommended for the site/occupancy layers, where
#                        several diGly sites can map to the same gene)
#
# Receives de_obj (from mod_de) and cfg. Writes optionally to
# <results_dir>/GSEA_rnk/ and always offers browser downloads.
# =============================================================================

# ---- core ranking logic -----------------------------------------------------

# Pull a tidy gene/logFC/pvalue/tstat frame out of the DE object for one layer
rnk_layer_frame <- function(de, layer, gene_col) {
  if (layer %in% c("occ", "site")) {
    df <- de$res
    logfc_col <- paste0("logFC_", layer)
    p_col     <- paste0("P_",     layer)
    t_col     <- paste0("t_",     layer)
    req_cols  <- c(gene_col, logfc_col, p_col, t_col)
    missing   <- setdiff(req_cols, names(df))
    if (length(missing))
      stop(sprintf("Combined results table is missing column(s): %s",
                   paste(missing, collapse = ", ")))

    dplyr::transmute(df,
      gene   = .data[[gene_col]],
      logFC  = .data[[logfc_col]],
      pvalue = .data[[p_col]],
      tstat  = .data[[t_col]])

  } else if (layer == "prot") {
    tt <- de$tt.prot
    if (is.null(tt)) stop("Protein-layer DE table (tt.prot) not available.")

    if (gene_col == "Genes" && !"Genes" %in% names(tt)) {
      # Protein layer often only carries Protein.Group -- borrow the
      # protein -> gene mapping already resolved on the site layer.
      map <- de$res |>
        dplyr::filter(!is.na(Protein.Group)) |>
        dplyr::distinct(Protein.Group, Genes)
      tt <- dplyr::left_join(tt, map, by = "Protein.Group")
    }
    if (!gene_col %in% names(tt))
      stop(sprintf("Protein-layer table has no '%s' column.", gene_col))

    dplyr::transmute(tt,
      gene   = .data[[gene_col]],
      logFC  = logFC,
      pvalue = P.Value,
      tstat  = t)

  } else {
    stop(sprintf("Unknown layer '%s'.", layer))
  }
}

# Build a ranked (gene, score) table for one layer + method
build_rnk <- function(de, layer, gene_col = "Genes",
                       method = c("t_stat_stouffer", "t_stat",
                                  "logFC", "logFC_pvalue")) {
  method <- match.arg(method)
  df <- rnk_layer_frame(de, layer, gene_col)
  df <- df[!is.na(df$gene) & nzchar(trimws(as.character(df$gene))), , drop = FALSE]

  n_dup <- 0L

  if (method == "t_stat_stouffer") {
    out <- df |>
      dplyr::filter(!is.na(tstat)) |>
      dplyr::group_by(gene) |>
      dplyr::summarise(
        n_sites = dplyr::n(),
        score   = sum(tstat) / sqrt(n_sites),
        .groups = "drop") |>
      dplyr::arrange(dplyr::desc(score)) |>
      dplyr::select(gene, score, n_sites)

  } else if (method == "t_stat") {
    out <- df |>
      dplyr::filter(!is.na(tstat)) |>
      dplyr::transmute(gene, score = tstat) |>
      dplyr::arrange(dplyr::desc(score))

  } else if (method == "logFC") {
    out <- df |>
      dplyr::filter(!is.na(logFC)) |>
      dplyr::transmute(gene, score = logFC) |>
      dplyr::arrange(dplyr::desc(score))

  } else if (method == "logFC_pvalue") {
    out <- df |>
      dplyr::filter(!is.na(logFC), !is.na(pvalue)) |>
      dplyr::transmute(gene, score = logFC * -log10(pvalue + 1e-300)) |>
      dplyr::arrange(dplyr::desc(score))
  }

  # Stouffer already collapses to one row per gene; other methods can still
  # have several features (sites) per gene -- keep the best-ranked one, same
  # behaviour as the standalone script.
  if (method != "t_stat_stouffer") {
    n_dup <- sum(duplicated(out$gene))
    out   <- out[!duplicated(out$gene), , drop = FALSE]
  }

  list(data = out, method = method, layer = layer, n_dup = n_dup,
       n_genes = nrow(out))
}

# ---- UI ----------------------------------------------------------------------

mod_rnk_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h4(" GSEA Rank File (.rnk)"),
    card(
      card_body(
        tags$p(class = "text-muted small",
          "Generates GSEA Desktop-compatible .rnk files (gene", tags$code("\\t"), "score, ",
          "no header) from the DE results computed in the DE Analysis tab."),
        fluidRow(
          column(3,
            checkboxGroupInput(ns("layers"), "Layer(s)",
                                choices  = c("Occupancy" = "occ",
                                             "Site"       = "site",
                                             "Protein"    = "prot"),
                                selected = "occ")
          ),
          column(3,
            selectInput(ns("method"), "Ranking method",
                        choices = c(
                          "t-stat, Stouffer (gene-level, recommended)" = "t_stat_stouffer",
                          "t-statistic"                                = "t_stat",
                          "logFC"                                      = "logFC",
                          "logFC × -log10(p)"                     = "logFC_pvalue"),
                        selected = "t_stat_stouffer"),
            conditionalPanel(
              condition = sprintf("input['%s'] == 't_stat_stouffer'", ns("method")),
              helpText("z", tags$sub("combined"), " = Σ t", tags$sub("site"),
                       " / √ n", tags$sub("sites"),
                       " -- combines every diGly site mapped to a gene into one score.")
            )
          ),
          column(3,
            selectInput(ns("gene_col"), "Gene identifier column",
                        choices  = c("Gene symbol (Genes)" = "Genes",
                                     "Protein group"        = "Protein.Group"),
                        selected = "Genes")
          ),
          column(3,
            checkboxInput(ns("write_disk"), "Also write to Results/GSEA_rnk/", TRUE),
            tags$div(style = "margin-top:10px;",
              actionButton(ns("generate"), "Generate .rnk file(s)", class = "btn-primary"),
              spinner_ui(ns("sp_rnk"))
            )
          )
        ),
        verbatimTextOutput(ns("rnk_log"))
      )
    ),

    h4(" Preview & Download"),
    uiOutput(ns("rnk_tabs"))
  )
}

# ---- server --------------------------------------------------------------

mod_rnk_server <- function(id, de_obj, cfg, status = reactiveVal("idle")) {
  moduleServer(id, function(input, output, session) {

    log_msg      <- reactiveVal("")
    rnk_results  <- reactiveVal(NULL)   # named list, one entry per layer
    sp_rnk_state <- reactiveVal("idle")
    spinner_server("sp_rnk", sp_rnk_state)

    layer_label <- c(occ = "Occupancy", site = "Site", prot = "Protein")

    append_log <- function(...) {
      msg <- paste0(paste(...), "\n")
      log_msg(paste0(log_msg(), msg))
    }

    observeEvent(input$generate, {
      de <- de_obj(); c <- cfg()
      req(de, c, input$layers)
      log_msg(""); rnk_results(NULL)
      status("running"); sp_rnk_state("running")

      tgt <- make.names(c$target_line %||% "target")
      out_dir <- if (nchar(c$results_dir %||% "") > 0)
        file.path(c$results_dir, "GSEA_rnk") else NULL
      if (input$write_disk && !is.null(out_dir))
        dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

      results <- list()

      withCallingHandlers(
        tryCatch({
          append_log(sprintf("=== Ranking method: %s ===", input$method))
          for (lyr in input$layers) {
            append_log(sprintf("\n-- %s layer --", layer_label[[lyr]]))
            r <- build_rnk(de, lyr, gene_col = input$gene_col, method = input$method)
            results[[lyr]] <- r

            append_log(sprintf("  Ranked genes: %d", r$n_genes))
            if (r$n_dup > 0)
              append_log(sprintf("  [!] %d duplicate gene symbols collapsed (kept best-ranked).",
                                 r$n_dup))
            if (nrow(r$data) > 0)
              append_log(sprintf("  Range: %.3f to %.3f",
                                 min(r$data$score), max(r$data$score)))

            if (input$write_disk && !is.null(out_dir)) {
              fname <- sprintf("%s_%s_%s.rnk", tgt, lyr, input$method)
              fpath <- file.path(out_dir, fname)
              readr::write_tsv(r$data[, c("gene", "score")], fpath, col_names = FALSE)
              append_log(sprintf("  [OK] Saved: %s", fname))
            }
          }

          rnk_results(results)
          append_log("\n[OK] .rnk generation complete.")
          status("done"); sp_rnk_state("done")
        }, error = function(e) {
          append_log(paste0("[ERR] ERROR: ", e$message))
          status("error"); sp_rnk_state("error")
        }),
        warning = function(w) { append_log(paste("[!]", conditionMessage(w))); invokeRestart("muffleWarning") }
      )
    })

    output$rnk_log <- renderPrint({ cat(log_msg()) })

    # ---- dynamic preview / download tabs, one per generated layer ---------
    output$rnk_tabs <- renderUI({
      res <- rnk_results()
      if (is.null(res) || length(res) == 0)
        return(tags$p(class = "text-muted", "Generate a rank file to preview it here."))

      ns <- session$ns
      do.call(tabsetPanel, c(
        list(id = ns("rnk_tabset")),
        lapply(names(res), function(lyr) {
          tabPanel(
            title = layer_label[[lyr]],
            tags$div(style = "margin-top:12px;",
              downloadButton(ns(paste0("dl_", lyr)),
                             sprintf("Download %s .rnk", layer_label[[lyr]])),
              tags$span(class = "text-muted small", style = "margin-left:10px;",
                sprintf("%d genes", res[[lyr]]$n_genes))
            ),
            DT::dataTableOutput(ns(paste0("tbl_", lyr)))
          )
        })
      ))
    })

    # Render a preview table + download handler for each possible layer.
    # Guarded with req() so it's a no-op until that layer has been generated.
    lapply(c("occ", "site", "prot"), function(lyr) {
      local({
        layer <- lyr
        output[[paste0("tbl_", layer)]] <- DT::renderDataTable({
          res <- rnk_results(); req(res, res[[layer]])
          d <- res[[layer]]$data
          d$score <- round(d$score, 5)
          d
        }, options = list(pageLength = 15, scrollX = TRUE), rownames = FALSE)

        output[[paste0("dl_", layer)]] <- downloadHandler(
          filename = function() {
            c <- cfg()
            sprintf("%s_%s_%s.rnk",
                   make.names(c$target_line %||% "target"), layer, input$method)
          },
          content = function(file) {
            res <- rnk_results(); req(res, res[[layer]])
            readr::write_tsv(res[[layer]]$data[, c("gene", "score")], file, col_names = FALSE)
          }
        )
      })
    })

    `%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && nzchar(a[[1]])) a[[1]] else b

    invisible(rnk_results)
  })
}
