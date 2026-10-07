# =============================================================================
# mod_psp -- PhosphoSitePlus ubiquitination site annotation
#
# Cross-references DE results (de$res) against a PSP flat file
# (Ubiquitination_site_dataset_Human.tsv or equivalent).
# PSP files have a multi-line license preamble; the header row is detected
# automatically as the first line starting with "GENE\t".
#
# Adds to de$res:
#   res_pos         -- numeric lysine position extracted from site_id
#   key             -- ACC_ID_respos join key
#   reported_in_PSP -- logical
#   LT_LIT, MS_LIT, MS_CST, SITE_GRP_ID -- PSP evidence columns
# =============================================================================

mod_psp_ui <- function(id) {
  ns <- NS(id)
  tagList(

    h4("PhosphoSitePlus Ubiquitination Annotation"),
    card(
      card_body(
        tags$p(class = "text-muted small mb-2",
          "Cross-references DE sites against the PSP ubiquitination site dataset. ",
          "Download the human ubiquitination flat file from ",
          tags$a("PhosphoSitePlus",
                 href   = "https://www.phosphosite.org/staticDownloads",
                 target = "_blank"),
          " (free academic registration required). ",
          "The multi-line license preamble is skipped automatically."),

        # PSP file selector ---------------------------------------------------
        fluidRow(
          column(8,
            hybrid_path_ui(ns, "psp_file",
                           "PSP ubiquitination flat file (.tsv / .gz)",
                           type = "file",
                           placeholder = "/path/to/Ubiquitination_site_dataset.tsv")
          ),
          column(4,
            tags$div(style = "margin-top:28px;",
              actionButton(ns("run"), "Annotate with PSP",
                           class = "btn-primary w-100"),
              spinner_ui(ns("sp_psp"))
            )
          )
        ),
        verbatimTextOutput(ns("psp_log"))
      )
    ),

    # Results -----------------------------------------------------------------
    conditionalPanel(
      condition = sprintf("output['%s']", ns("done")),

      h4("Known / Novel Summary"),
      card(
        card_body(
          fluidRow(
            column(4, verbatimTextOutput(ns("summary_print"))),
            column(8, plotOutput(ns("summary_plot"), height = "200px"))
          )
        )
      ),

      h4("Annotated DE Results"),
      card(
        card_body(
          tags$p(class = "text-muted small mb-1",
            tags$b("reported_in_PSP:"), " site found in PSP database. ",
            tags$b("LT_LIT:"), " low-throughput literature references. ",
            tags$b("MS_LIT:"), " MS literature references. ",
            tags$b("MS_CST:"), " Cell Signaling Technology MS references."),
          DT::dataTableOutput(ns("annotated_table"))
        )
      ),

      h4("Novel Sites (not in PSP)"),
      card(
        card_body(
          tags$p(class = "text-muted small mb-1",
            "Significant sites (occupancy FDR < threshold) absent from PSP. ",
            "Sorted by FDR."),
          DT::dataTableOutput(ns("novel_table"))
        )
      ),

      h4("Export"),
      card(
        card_body(
          fluidRow(
            column(4,
              actionButton(ns("export_disk"), "Save TSV to Results directory",
                           class = "btn-success")
            ),
            column(4,
              downloadButton(ns("dl_annotated"), "Download annotated TSV")
            ),
            column(4,
              downloadButton(ns("dl_novel"),     "Download novel sites TSV")
            )
          ),
          verbatimTextOutput(ns("export_log"))
        )
      )
    )
  )
}


mod_psp_server <- function(id, de_obj, cfg) {
  moduleServer(id, function(input, output, session) {

    psp_res      <- reactiveVal(NULL)
    sp_psp_state <- reactiveVal("idle")
    log_msg      <- reactiveVal("")
    export_log   <- reactiveVal("")
    spinner_server("sp_psp", sp_psp_state)

    output$done <- reactive(!is.null(psp_res()))
    outputOptions(output, "done", suspendWhenHidden = FALSE)

    # shinyFiles root for PSP file
    roots <- c(Computer = if (.Platform$OS.type == "windows") "C:/" else "/")
    shinyFiles::shinyFileChoose(input, "btn_psp_file", roots = roots,
                                filetypes = c("tsv", "gz", "txt", "csv"))
    observeEvent(input$btn_psp_file, {
      p <- tryCatch(
        shinyFiles::parseFilePaths(roots, input$btn_psp_file)$datapath[[1]],
        error = function(e) NULL)
      if (!is.null(p) && nzchar(p))
        updateTextInput(session, "txt_psp_file", value = p)
    }, ignoreInit = TRUE, ignoreNULL = TRUE)

    append_log <- function(msg) log_msg(paste0(log_msg(), msg, "\n"))

    alog <- function(msg) log_msg(paste0(log_msg(), msg, "\n"))

    observeEvent(input$run, {
      de   <- de_obj();  c <- cfg()
      req(de, c)

      psp_path <- trimws(input$txt_psp_file %||% "")
      if (!nzchar(psp_path) || !file.exists(psp_path)) {
        showNotification("PSP file not found — select a valid path.",
                         type = "error"); return()
      }

      log_msg(""); psp_res(NULL); sp_psp_state("running")

      tryCatch({
        # ── Load PSP flat file (skip license preamble) ─────────────────────
        alog("=== Loading PSP flat file ===")
        psp_lines <- readr::read_lines(psp_path, n_max = 50)
        header_idx <- which(stringr::str_detect(psp_lines, "^GENE\\t"))[1]
        if (is.na(header_idx))
          stop("Could not find header row starting with 'GENE\\t'. ",
               "Check the file is a PSP ubiquitination dataset.")

        psp_raw <- readr::read_tsv(
          psp_path,
          skip          = header_idx - 1,
          show_col_types = FALSE,
          progress      = FALSE)

        # Filter to human and extract residue position
        psp <- psp_raw |>
          dplyr::filter(tolower(ORGANISM) == "human") |>
          dplyr::mutate(
            res_pos = stringr::str_extract(MOD_RSD, "(?<=K)\\d+(?=-ub)"),
            key     = paste(ACC_ID, res_pos, sep = "_"))

        alog(sprintf("  PSP entries (human): %d | unique sites: %d",
                     nrow(psp), dplyr::n_distinct(psp$key)))

        psp_keys     <- unique(psp$key)
        psp_evidence <- psp |>
          dplyr::select(key,
                        dplyr::any_of(c("LT_LIT","MS_LIT","MS_CST","SITE_GRP_ID"))) |>
          dplyr::distinct(key, .keep_all = TRUE)

        # ── Prepare DE results ─────────────────────────────────────────────
        alog("\n=== Matching DE sites against PSP ===")
        res <- de$res

        # site_id format examples: "P15121_K263|1", "Q9UHX1_K108|2"
        # site_label format:       "GENE_K263"
        # PSP key format:          "ACC_ID_res_pos" e.g. "P15121_263"
        #
        # Strategy:
        #   res_pos  -- extract digits after _K in site_id (before optional |)
        #   acc_id   -- everything before _K<digits> in site_id
        res <- res |>
          dplyr::mutate(
            # Position: digits after _K, optionally followed by | or end
            res_pos = stringr::str_extract(site_id, "(?<=_K)(\\d+)(?=[|]|$)"),
            # Accession: everything before the _K<digits> suffix
            acc_id  = stringr::str_remove(site_id, "_K\\d+([|].*)?$"),
            # Fall back to leading entry of Protein.Group if acc_id looks wrong
            acc_id  = dplyr::if_else(
              !is.na(Protein.Group) & nchar(Protein.Group) > 0,
              stringr::str_split_fixed(Protein.Group, ";", 2)[, 1],
              acc_id),
            key = dplyr::if_else(
              !is.na(acc_id) & !is.na(res_pos),
              paste(acc_id, res_pos, sep = "_"),
              NA_character_))

        res <- res |>
          dplyr::mutate(reported_in_PSP = key %in% psp_keys) |>
          dplyr::left_join(psp_evidence, by = "key")

        n_known <- sum(res$reported_in_PSP, na.rm = TRUE)
        n_novel <- sum(!res$reported_in_PSP, na.rm = TRUE)
        alog(sprintf("  Known (in PSP): %d  |  Novel: %d  |  Match rate: %.1f%%",
                     n_known, n_novel, 100 * n_known / nrow(res)))

        # ── Novel significant sites ────────────────────────────────────────
        fdr_col <- if ("FDR_occ" %in% names(res)) "FDR_occ" else
                   if ("FDR_site" %in% names(res)) "FDR_site" else NULL

        novel <- res |>
          dplyr::filter(!reported_in_PSP)
        if (!is.null(fdr_col))
          novel <- novel |>
            dplyr::filter(.data[[fdr_col]] < c$fdr_threshold) |>
            dplyr::arrange(.data[[fdr_col]])

        alog(sprintf("  Novel significant sites (FDR < %g): %d",
                     c$fdr_threshold, nrow(novel)))
        alog("\n[OK] PSP annotation complete.")
        sp_psp_state("done")

        psp_res(list(
          annotated = res,
          novel     = novel,
          n_known   = n_known,
          n_novel   = n_novel,
          fdr_col   = fdr_col
        ))

      }, error = function(e) {
        alog(paste0("[ERR] ", e$message))
        sp_psp_state("error")
      })
    })

    # ── Outputs ───────────────────────────────────────────────────────────────
    output$psp_log <- renderPrint({ cat(log_msg()) })

    output$summary_print <- renderPrint({
      r <- psp_res(); req(r)
      tot <- r$n_known + r$n_novel
      cat(sprintf("Total sites       : %d\n", tot))
      cat(sprintf("Known in PSP      : %d (%.1f%%)\n",
                  r$n_known, 100 * r$n_known / tot))
      cat(sprintf("Novel (not in PSP): %d (%.1f%%)\n",
                  r$n_novel, 100 * r$n_novel / tot))
    })

    output$summary_plot <- renderPlot({
      r <- psp_res(); req(r)
      df <- data.frame(
        Status = c("Known in PSP", "Novel"),
        n      = c(r$n_known, r$n_novel),
        fill   = c("#4e9a8c", "#c0392b"))
      ggplot2::ggplot(df, ggplot2::aes(x = "", y = n, fill = Status)) +
        ggplot2::geom_col(width = 0.5) +
        ggplot2::scale_fill_manual(
          values = c("Known in PSP" = "#4e9a8c", "Novel" = "#c0392b")) +
        ggplot2::geom_text(ggplot2::aes(label = paste0(n, "\n(",
                           round(100*n/sum(n),1), "%)")),
                           position = ggplot2::position_stack(vjust = 0.5),
                           size = 4, color = "white", fontface = "bold") +
        ggplot2::coord_flip() +
        ggplot2::labs(x = NULL, y = "Sites", title = "PSP annotation coverage") +
        ggplot2::theme_bw(base_size = 12) +
        ggplot2::theme(legend.position = "right")
    })

    output$annotated_table <- DT::renderDataTable({
      r <- psp_res(); req(r)
      disp_cols <- intersect(
        c("site_id","site_label","Genes","logFC_occ","FDR_occ",
          "logFC_site","FDR_site","classification","direction",
          "reported_in_PSP","LT_LIT","MS_LIT","MS_CST","SITE_GRP_ID"),
        names(r$annotated))
      r$annotated |>
        dplyr::select(dplyr::all_of(disp_cols)) |>
        dplyr::arrange(reported_in_PSP, dplyr::pick(dplyr::any_of("FDR_occ")))
    },
    options  = list(scrollX = TRUE, pageLength = 15),
    filter   = "top",
    rownames = FALSE)

    output$novel_table <- DT::renderDataTable({
      r <- psp_res(); req(r)
      disp_cols <- intersect(
        c("site_id","site_label","Genes","Protein.Group",
          "logFC_occ","P_occ","FDR_occ",
          "logFC_site","FDR_site","classification","direction"),
        names(r$novel))
      r$novel |> dplyr::select(dplyr::all_of(disp_cols))
    },
    options  = list(scrollX = TRUE, pageLength = 15),
    filter   = "top",
    rownames = FALSE)

    # ── Export ────────────────────────────────────────────────────────────────
    `%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && nzchar(a[[1]])) a[[1]] else b

    observeEvent(input$export_disk, {
      r <- psp_res(); c <- cfg(); req(r, c)
      export_log("")
      if (!nzchar(c$results_dir %||% "")) {
        export_log("Results directory not set in Configuration.")
        return()
      }
      dir.create(c$results_dir, showWarnings = FALSE, recursive = TRUE)
      tgt  <- make.names(c$target_line)
      path <- file.path(c$results_dir,
                        sprintf("diGly_%s_vs_rest_PSP_annotated.tsv", tgt))
      tryCatch({
        readr::write_tsv(r$annotated, path)
        export_log(sprintf("[OK] Saved: %s", basename(path)))
      }, error = function(e) export_log(paste("[ERR]", e$message)))
    })

    output$export_log <- renderPrint({ cat(export_log()) })

    output$dl_annotated <- downloadHandler(
      filename = function() {
        sprintf("diGly_%s_vs_rest_PSP_annotated.tsv",
                make.names(cfg()$target_line %||% "analysis"))
      },
      content = function(file) {
        r <- psp_res(); req(r)
        readr::write_tsv(r$annotated, file)
      }
    )

    output$dl_novel <- downloadHandler(
      filename = function() {
        sprintf("diGly_%s_vs_rest_PSP_novel_sites.tsv",
                make.names(cfg()$target_line %||% "analysis"))
      },
      content = function(file) {
        r <- psp_res(); req(r)
        readr::write_tsv(r$novel, file)
      }
    )
  })
}
