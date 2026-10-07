# =============================================================================
# mod_data -- File paths, run-name token map, scan + interactive run selection
#
# Tab order: paths are set here first, scan auto-fires once both parquets
# exist, then the user selects/annotates runs and loads. Configuration
# (design, filters, thresholds) comes in the next tab.
#
# Returns a reactiveVal list:
#   $y.dg, $y.tp           -- ELists (post-load)
#   $meta                  -- sample metadata data.frame
#   $targets               -- targets data.frame
#   $prot_seq              -- named character vector from FASTA
#   $paths                 -- file/dir paths forwarded to cfg and export
#   $cell_lines_present    -- character vector of unique cell lines loaded
# =============================================================================

# Hybrid path input helper ----------------------------------------------------
# Used in mod_data_ui; defined here so it is available when mod_data.R is sourced.
hybrid_path_ui <- function(ns, id, label, type = "file", placeholder = "") {
  btn <- if (type == "file") {
    shinyFiles::shinyFilesButton(
      ns(paste0("btn_", id)),
      label       = "\U0001F4C2",
      title       = paste("Select", label),
      multiple    = FALSE,
      buttonType  = "outline-secondary",
      class       = "btn-sm")
  } else {
    shinyFiles::shinyDirButton(
      ns(paste0("btn_", id)),
      label      = "\U0001F4C2",
      title      = paste("Select", label),
      buttonType = "outline-secondary",
      class      = "btn-sm")
  }
  tagList(
    tags$label(label, class = "form-label mb-1 mt-2"),
    tags$div(
      class = "input-group input-group-sm mb-1",
      tags$span(class = "input-group-text p-0 border-0 bg-transparent", btn),
      textInput(ns(paste0("txt_", id)),
                label       = NULL,
                placeholder = placeholder,
                width       = "100%") |>
        tagAppendAttributes(style = "flex:1;")
    ),
    uiOutput(ns(paste0("disp_", id)))
  )
}

# =============================================================================
# UI
# =============================================================================
mod_data_ui <- function(id) {
  ns <- NS(id)
  tagList(

    tags$style(HTML("
      .path-display {
        font-size: 0.82em; color: #555; background: #f8f9fa;
        border: 1px solid #dee2e6; border-radius: 4px;
        padding: 3px 8px; min-height: 24px; word-break: break-all;
        margin-bottom: 4px;
      }
      .path-display.ok   { border-color: #28a745; color: #155724; background: #d4edda; }
      .path-display.none { color: #aaa; font-style: italic; }
      .unparsed { color: #856404; font-weight: 600; }
    ")),

    # Step 1: File paths ------------------------------------------------------
    h4("Step 1 \u2014 Select files & directories"),
    card(
      card_body(
        tags$p(class = "text-muted small mb-3",
               "\U0001F4C2 Click Browse or paste a path. ",
               "Scanning fires automatically once both parquets are set."),
        fluidRow(
          column(6,
            hybrid_path_ui(ns, "digly",   "diGly parquet file",             "file",
                           placeholder = "/path/to/GG_report.parquet"),
            hybrid_path_ui(ns, "proteome","Proteome parquet file",           "file",
                           placeholder = "/path/to/Proteome_report.parquet"),
            hybrid_path_ui(ns, "fasta",   "FASTA file",                      "file",
                           placeholder = "/path/to/uniprotkb.fasta"),
            hybrid_path_ui(ns, "map_csv", "Proteome run map CSV (optional)", "file",
                           placeholder = "(optional)")
          ),
          column(6,
            hybrid_path_ui(ns, "results", "Results directory",   "dir",
                           placeholder = "/path/to/RESULTS"),
            hybrid_path_ui(ns, "figures", "Figures directory",   "dir",
                           placeholder = "/path/to/FIGURES"),
            hybrid_path_ui(ns, "limpa",   "limpa RDS directory", "dir",
                           placeholder = "/path/to/LIMPA")
          )
        )
      )
    ),

    # Step 2: Token map -------------------------------------------------------
    h4("Step 2 \u2014 Run-name token map"),
    card(
      card_body(
        fluidRow(
          column(6,
            textAreaInput(ns("line_map_raw"),
                          "Token = Cell line  (one per line, e.g.  ASPC=AsPC-1)",
                          value = "",
                          rows  = 5),
            helpText("Tokens are matched case-insensitively in run names. ",
                     "Leave blank to assign cell lines manually in the table below.")
          ),
          column(6,
            verbatimTextOutput(ns("scan_log"))
          )
        )
      )
    ),

    # Step 3: Run selection table (appears after scan) ------------------------
    conditionalPanel(
      condition = sprintf("output['%s']", ns("scanned")),

      h4("Step 3 \u2014 Select & annotate runs"),
      tags$p(class = "text-muted small mb-3",
        tags$b("Cell line"), " and ", tags$b("Replicate"),
        " cells are editable \u2014 double-click to change. ",
        tags$span(class = "unparsed", "Yellow rows"),
        " were not parsed automatically and must be assigned before loading."
      ),

      card(
        card_header(
          fluidRow(
            column(5, tags$b("diGly runs"),
                      uiOutput(ns("dg_parse_badge"), inline = TRUE)),
            column(2, actionButton(ns("dg_all"),      "All",
                                   class = "btn-sm btn-outline-secondary w-100")),
            column(2, actionButton(ns("dg_none"),     "None",
                                   class = "btn-sm btn-outline-secondary w-100")),
            column(3, actionButton(ns("dg_unparsed"), "Unparsed",
                                   class = "btn-sm btn-outline-warning w-100"))
          )
        ),
        card_body(DT::dataTableOutput(ns("dg_run_table")))
      ),

      card(
        card_header(
          fluidRow(
            column(5, tags$b("Proteome runs"),
                      uiOutput(ns("tp_parse_badge"), inline = TRUE)),
            column(2, actionButton(ns("tp_all"),      "All",
                                   class = "btn-sm btn-outline-secondary w-100")),
            column(2, actionButton(ns("tp_none"),     "None",
                                   class = "btn-sm btn-outline-secondary w-100")),
            column(3, actionButton(ns("tp_unparsed"), "Unparsed",
                                   class = "btn-sm btn-outline-warning w-100"))
          )
        ),
        card_body(DT::dataTableOutput(ns("tp_run_table")))
      ),

      # Step 4: Load ----------------------------------------------------------
      h4("Step 4 \u2014 Load selected runs"),
      card(
        card_body(
          fluidRow(
            column(9,
              tags$p(class = "text-muted small mb-2",
                "Calls readDIANN() on selected runs and reads the FASTA. ",
                "~1\u20132 min for large parquets.")
            ),
            column(3,
              actionButton(ns("load"), "Load selected runs",
                           class = "btn-primary w-100 mt-1"),
              spinner_ui(ns("sp_load"))
            )
          ),
          verbatimTextOutput(ns("load_log"))
        )
      ),

      conditionalPanel(
        condition = sprintf("output['%s']", ns("loaded")),
        h4("Final sample sheet"),
        card(card_body(DT::dataTableOutput(ns("meta_table"))))
      )
    )
  )
}

# =============================================================================
# Server
# =============================================================================
mod_data_server <- function(id, status = reactiveVal("idle")) {
  moduleServer(id, function(input, output, session) {

    # shinyFiles roots --------------------------------------------------------
    roots <- c(Computer = if (.Platform$OS.type == "windows") "C:/" else "/")

    shinyFiles::shinyFileChoose(input, "btn_digly",   roots = roots, filetypes = "parquet")
    shinyFiles::shinyFileChoose(input, "btn_proteome",roots = roots, filetypes = "parquet")
    shinyFiles::shinyFileChoose(input, "btn_fasta",   roots = roots,
                                filetypes = c("fasta", "fa", "faa", "gz"))
    shinyFiles::shinyFileChoose(input, "btn_map_csv", roots = roots, filetypes = "csv")
    shinyFiles::shinyDirChoose(input,  "btn_results", roots = roots)
    shinyFiles::shinyDirChoose(input,  "btn_figures", roots = roots)
    shinyFiles::shinyDirChoose(input,  "btn_limpa",   roots = roots)

    push_file <- function(btn_id, txt_id) {
      observeEvent(input[[btn_id]], {
        p <- tryCatch(
          shinyFiles::parseFilePaths(roots, input[[btn_id]])$datapath[[1]],
          error = function(e) NULL)
        if (!is.null(p) && nzchar(p)) updateTextInput(session, txt_id, value = p)
      }, ignoreInit = TRUE, ignoreNULL = TRUE)
    }
    push_dir <- function(btn_id, txt_id) {
      observeEvent(input[[btn_id]], {
        p <- tryCatch(
          shinyFiles::parseDirPath(roots, input[[btn_id]])[[1]],
          error = function(e) NULL)
        if (!is.null(p) && nzchar(p)) updateTextInput(session, txt_id, value = p)
      }, ignoreInit = TRUE, ignoreNULL = TRUE)
    }

    push_file("btn_digly",    "txt_digly")
    push_file("btn_proteome", "txt_proteome")
    push_file("btn_fasta",    "txt_fasta")
    push_file("btn_map_csv",  "txt_map_csv")
    push_dir("btn_results",   "txt_results")
    push_dir("btn_figures",   "txt_figures")
    push_dir("btn_limpa",     "txt_limpa")

    # Reactive paths ----------------------------------------------------------
    `%||%` <- function(a, b) {
      if (!is.null(a) && length(a) > 0 && nzchar(a[[1]])) a[[1]] else b
    }

    paths <- reactive({
      list(
        digly_parquet    = trimws(input$txt_digly    %||% ""),
        proteome_parquet = trimws(input$txt_proteome %||% ""),
        fasta            = trimws(input$txt_fasta    %||% ""),
        proteome_map_csv = trimws(input$txt_map_csv  %||% ""),
        results_dir      = trimws(input$txt_results  %||% ""),
        figures_dir      = trimws(input$txt_figures  %||% ""),
        limpa_dir        = trimws(input$txt_limpa    %||% "")
      )
    })

    # Path status displays ----------------------------------------------------
    fmt_path <- function(p, optional = FALSE) {
      if (is.null(p) || !nzchar(p)) {
        lbl <- if (optional) "optional" else ""
        return(tags$div(class = "path-display none", lbl))
      }
      ok  <- file.exists(p)
      cls <- paste("path-display", if (ok) "ok" else "")
      icn <- if (ok) icon("circle-check") else icon("triangle-exclamation")
      tags$div(class = cls, icn, p)
    }

    output$disp_digly    <- renderUI({ fmt_path(paths()$digly_parquet) })
    output$disp_proteome <- renderUI({ fmt_path(paths()$proteome_parquet) })
    output$disp_fasta    <- renderUI({ fmt_path(paths()$fasta) })
    output$disp_map_csv  <- renderUI({ fmt_path(paths()$proteome_map_csv, optional = TRUE) })
    output$disp_results  <- renderUI({ fmt_path(paths()$results_dir) })
    output$disp_figures  <- renderUI({ fmt_path(paths()$figures_dir) })
    output$disp_limpa    <- renderUI({ fmt_path(paths()$limpa_dir) })

    # Spinner states -----------------------------------------------------------
    sp_load_state <- reactiveVal("idle")
    spinner_server("sp_load", sp_load_state)

    # Shared reactive state ---------------------------------------------------
    scan_log  <- reactiveVal("")
    load_log  <- reactiveVal("")
    dg_df     <- reactiveVal(NULL)
    tp_df     <- reactiveVal(NULL)
    scan_info <- reactiveVal(NULL)
    data_obj  <- reactiveVal(NULL)

    output$scanned <- reactive(!is.null(dg_df()))
    output$loaded  <- reactive(!is.null(data_obj()))
    outputOptions(output, "scanned", suspendWhenHidden = FALSE)
    outputOptions(output, "loaded",  suspendWhenHidden = FALSE)

    alog <- function(rv, msg) rv(paste0(rv(), msg, "\n"))

    # Run-name parser helpers -------------------------------------------------
    parse_line_map <- function(raw) {
      pairs <- Filter(nchar, trimws(strsplit(raw, "\n")[[1]]))
      if (length(pairs) == 0) return(character(0))
      kv    <- strsplit(pairs, "=")
      valid <- sapply(kv, length) == 2
      if (!any(valid)) return(character(0))
      setNames(sapply(kv[valid], `[[`, 2), sapply(kv[valid], `[[`, 1))
    }

    parse_runs <- function(runs, line_map) {
      df <- dplyr::tibble(
        run       = runs,
        replicate = as.integer(
          stringr::str_match(runs, "_([0-9]+)_R[A-Z][0-9]+_1_")[, 2]),
        well      = stringr::str_match(runs, "_(R[A-Z][0-9]+)_1_")[, 2],
        inj_order = as.integer(stringr::str_match(runs, "_([0-9]+)$")[, 2]),
        cell_line = NA_character_,
        plate_row = NA_character_
      )
      if (length(line_map) > 0) {
        tokens    <- paste(names(line_map), collapse = "|")
        tok_match <- toupper(stringr::str_match(
          runs,
          stringr::regex(paste0("(", tokens, ")"), ignore_case = TRUE)
        )[, 2])
        df$cell_line <- as.character(unname(line_map[tok_match]))
        df$plate_row <- stringr::str_sub(df$well, 1, 2)
      }
      df
    }

    # Auto-scan when both parquet paths become valid --------------------------
    observeEvent(
      list(input$txt_digly, input$txt_proteome, input$line_map_raw),
      {
        p <- paths()
        if (!nzchar(p$digly_parquet) || !nzchar(p$proteome_parquet)) return()
        if (!file.exists(p$digly_parquet) || !file.exists(p$proteome_parquet)) return()

        scan_log(""); dg_df(NULL); tp_df(NULL); data_obj(NULL); scan_info(NULL)
        lm <- parse_line_map(input$line_map_raw)

        withCallingHandlers(
          tryCatch({
            alog(scan_log, "=== Scanning diGly parquet ===")
            ds_dg    <- arrow::open_dataset(p$digly_parquet)
            avail_dg <- names(ds_dg)
            dg_runs  <- ds_dg |>
              dplyr::select(dplyr::any_of("Run")) |>
              dplyr::collect() |>
              dplyr::distinct() |>
              dplyr::pull(1)
            dg_counts <- ds_dg |>
              dplyr::group_by(Run) |>
              dplyr::summarise(n_precursors = dplyr::n(), .groups = "drop") |>
              dplyr::collect()
            alog(scan_log, paste0("  ", length(dg_runs), " runs found."))

            alog(scan_log, "\n=== Scanning proteome parquet ===")
            ds_tp   <- arrow::open_dataset(p$proteome_parquet)
            tp_runs <- ds_tp |>
              dplyr::select(dplyr::any_of("Run")) |>
              dplyr::collect() |>
              dplyr::distinct() |>
              dplyr::pull(1)
            tp_counts <- ds_tp |>
              dplyr::group_by(Run) |>
              dplyr::summarise(n_precursors = dplyr::n(), .groups = "drop") |>
              dplyr::collect()
            alog(scan_log, paste0("  ", length(tp_runs), " runs found."))

            dg_parsed <- parse_runs(dg_runs, lm)
            tp_parsed <- parse_runs(tp_runs, lm)

            dg_df(dg_parsed |>
              dplyr::left_join(dg_counts, by = c("run" = "Run")) |>
              dplyr::arrange(inj_order))
            tp_df(tp_parsed |>
              dplyr::left_join(tp_counts, by = c("run" = "Run")) |>
              dplyr::arrange(inj_order))

            n_na_dg <- sum(is.na(dg_parsed$cell_line))
            n_na_tp <- sum(is.na(tp_parsed$cell_line))
            if (n_na_dg > 0)
              alog(scan_log, paste0("  WARNING: ", n_na_dg,
                                    " diGly run(s) not parsed - assign in table."))
            if (n_na_tp > 0)
              alog(scan_log, paste0("  WARNING: ", n_na_tp,
                                    " proteome run(s) not parsed - assign in table."))
            alog(scan_log, "\nScan complete. Edit cell lines if needed, then Load.")
            scan_info(list(avail_dg = avail_dg))

          }, error = function(e) alog(scan_log, paste0("ERROR: ", e$message))),
          message = function(m) {
            alog(scan_log, conditionMessage(m))
            invokeRestart("muffleMessage")
          },
          warning = function(w) {
            alog(scan_log, paste0("WARNING: ", conditionMessage(w)))
            invokeRestart("muffleWarning")
          }
        )
      },
      ignoreInit = TRUE
    )

    # Parse badges ------------------------------------------------------------
    make_badge <- function(df) {
      if (is.null(df)) return(NULL)
      n_na  <- sum(is.na(df$cell_line))
      n_tot <- nrow(df)
      if (n_na == 0)
        tags$span(class = "badge bg-success ms-2",
                  paste0(n_tot, " / ", n_tot, " parsed"))
      else
        tags$span(class = "badge bg-warning text-dark ms-2",
                  paste0(n_na, " unparsed"))
    }
    output$dg_parse_badge <- renderUI({ make_badge(dg_df()) })
    output$tp_parse_badge <- renderUI({ make_badge(tp_df()) })

    # Build DT ----------------------------------------------------------------
    to_display <- function(df) {
      df |> dplyr::transmute(
        Run            = run,
        `Cell line`    = ifelse(is.na(cell_line), "", as.character(cell_line)),
        Replicate      = replicate,
        Well           = well,
        `Inj order`    = inj_order,
        `N precursors` = n_precursors
      )
    }

    make_dt <- function(df) {
      disp <- to_display(df)
      na_0 <- which(is.na(df$cell_line)) - 1L
      cb   <- DT::JS(paste0(
        "function(row, data, index) {",
        "  if ([", paste(na_0, collapse = ","), "].indexOf(index) !== -1)",
        "    $('td', row).css('background-color', '#fff3cd');",
        "}"
      ))
      DT::datatable(
        disp,
        selection = list(mode = "multiple", selected = seq_len(nrow(disp))),
        editable  = list(target  = "cell",
                         disable = list(columns = c(0L, 3L, 4L, 5L))),
        rownames  = FALSE,
        options   = list(
          scrollX     = TRUE,
          pageLength  = 25,
          rowCallback = cb,
          columnDefs  = list(
            list(width = "38%", targets = 0),
            list(width = "12%", targets = 1)
          )
        )
      )
    }

    dg_proxy <- DT::dataTableProxy("dg_run_table")
    tp_proxy <- DT::dataTableProxy("tp_run_table")

    output$dg_run_table <- DT::renderDataTable({
      df <- dg_df(); req(df); make_dt(df)
    })
    output$tp_run_table <- DT::renderDataTable({
      df <- tp_df(); req(df); make_dt(df)
    })

    # Cell edits (cell line col = 1, replicate col = 2) -----------------------
    apply_edit <- function(rv, info) {
      df  <- rv()
      row <- info$row
      col <- info$col
      val <- info$value
      if (col == 1L) {
        df$cell_line[row] <- if (nzchar(trimws(val))) trimws(val) else NA_character_
      } else if (col == 2L) {
        n <- suppressWarnings(as.integer(val))
        if (is.na(n)) {
          showNotification("Replicate must be an integer.", type = "warning")
          return()
        }
        df$replicate[row] <- n
      }
      rv(df)
    }
    observeEvent(input$dg_run_table_cell_edit,
                 apply_edit(dg_df, input$dg_run_table_cell_edit))
    observeEvent(input$tp_run_table_cell_edit,
                 apply_edit(tp_df, input$tp_run_table_cell_edit))

    # Select all / none / unparsed buttons ------------------------------------
    observeEvent(input$dg_all, {
      df <- dg_df(); req(df); DT::selectRows(dg_proxy, seq_len(nrow(df)))
    })
    observeEvent(input$dg_none, { DT::selectRows(dg_proxy, NULL) })
    observeEvent(input$tp_all, {
      df <- tp_df(); req(df); DT::selectRows(tp_proxy, seq_len(nrow(df)))
    })
    observeEvent(input$tp_none, { DT::selectRows(tp_proxy, NULL) })
    observeEvent(input$dg_unparsed, {
      df <- dg_df(); req(df)
      DT::selectRows(dg_proxy, which(is.na(df$cell_line)))
    })
    observeEvent(input$tp_unparsed, {
      df <- tp_df(); req(df)
      DT::selectRows(tp_proxy, which(is.na(df$cell_line)))
    })

    # Load button -------------------------------------------------------------
    observeEvent(input$load, {
      si <- scan_info(); p <- paths()
      dg <- dg_df();     tp <- tp_df()
      req(si, p, dg, tp)

      sel_dg <- input$dg_run_table_rows_selected
      sel_tp <- input$tp_run_table_rows_selected

      if (length(sel_dg) == 0) {
        showNotification("Select at least one diGly run.", type = "error")
        return()
      }
      if (length(sel_tp) == 0) {
        showNotification("Select at least one proteome run.", type = "error")
        return()
      }
      if (any(is.na(dg$cell_line[sel_dg])) || any(is.na(tp$cell_line[sel_tp]))) {
        showNotification(
          "Some selected runs have no cell line assigned. Assign them first.",
          type = "error", duration = 8)
        return()
      }

      load_log(""); data_obj(NULL); status("running")
      sp_load_state("running")

      withCallingHandlers(
        tryCatch({
          chosen_dg <- dg$run[sel_dg]
          chosen_tp <- tp$run[sel_tp]
          alog(load_log, paste0("Selected: ", length(chosen_dg),
                                " diGly | ", length(chosen_tp), " proteome runs"))

          prog <- make_progress_stepper(session, "Loading data…", c(
            "Reading diGly parquet"   = 0.20,
            "Reading proteome parquet"= 0.40,
            "Building metadata"       = 0.60,
            "Reading FASTA"           = 0.85,
            "Finalising"              = 1.00))
          on.exit(prog$close(), add = TRUE)

          # read diGly --------------------------------------------------------
          alog(load_log, "\n=== Reading diGly parquet ===")
          prog$tick("Reading diGly parquet…")
          avail_dg   <- si$avail_dg
          # Read without q-value filtering so annotation columns
          # (Modified.Sequence, PTM.Q.Value, etc.) stay intact for all rows.
          # Precursor q-value filtering happens in the Filtering tab
          # using the threshold from Configuration.
          extra_cols <- c("Protein.Group", "Protein.Ids", "Genes", "Proteotypic",
                          "Stripped.Sequence", "Modified.Sequence",
                          "PTM.Site.Confidence", "PTM.Q.Value",
                          "Q.Value", "Lib.Q.Value", "Lib.PG.Q.Value", "PG.Q.Value")
          y.dg <- limpa::readDIANN(
            file          = p$digly_parquet,
            format        = "parquet",
            qty.column    = "Precursor.Normalised",
            q.columns     = character(0),
            q.cutoffs     = Inf,
            extra.columns = intersect(extra_cols, avail_dg))
          y.dg <- y.dg[, colnames(y.dg) %in% chosen_dg]
          alog(load_log, paste0("  ", nrow(y.dg), " precursors x ", ncol(y.dg),
                                " runs  (",
                                round(100 * mean(is.na(y.dg$E)), 1), "% missing)"))

          # read proteome -----------------------------------------------------
          alog(load_log, "\n=== Reading proteome parquet ===")
          prog$tick("Reading proteome parquet…")
          avail_tp  <- names(arrow::open_dataset(p$proteome_parquet))
          # Proteome: also read without q-filter; filtering handled downstream.
          y.tp <- limpa::readDIANN(
            file          = p$proteome_parquet,
            format        = "parquet",
            qty.column    = "Precursor.Normalised",
            q.columns     = character(0),
            q.cutoffs     = Inf,
            extra.columns = intersect(
              c("Protein.Group", "Genes", "Proteotypic",
                "Q.Value", "Lib.Q.Value", "Lib.PG.Q.Value", "PG.Q.Value"), avail_tp))
          y.tp <- y.tp[, colnames(y.tp) %in% chosen_tp]
          alog(load_log, paste0("  ", nrow(y.tp), " precursors x ", ncol(y.tp),
                                " runs  (",
                                round(100 * mean(is.na(y.tp$E)), 1), "% missing)"))

          # metadata ----------------------------------------------------------
          alog(load_log, "\n=== Building sample metadata ===")
          prog$tick("Building sample metadata…")
          meta_dg <- dg |>
            dplyr::filter(run %in% colnames(y.dg)) |>
            dplyr::mutate(cell_line = factor(cell_line), digly_run = run) |>
            dplyr::arrange(inj_order)

          tp_map <- tp |>
            dplyr::filter(run %in% colnames(y.tp)) |>
            dplyr::mutate(cell_line = factor(cell_line)) |>
            dplyr::select(cell_line, replicate, proteome_run = run)

          meta <- meta_dg |>
            dplyr::left_join(tp_map,
                             by      = c("cell_line", "replicate"),
                             relationship = "many-to-many") |>
            dplyr::filter(!is.na(proteome_run)) |>
            dplyr::group_by(digly_run) |>
            dplyr::slice(1) |>
            dplyr::ungroup()

          if (nrow(meta) == 0)
            stop(paste("No diGly-proteome sample pairs matched.",
                       "Check cell line and replicate assignments."))

          unmatched <- nrow(meta_dg) - nrow(meta)
          if (unmatched > 0)
            alog(load_log, paste0("  WARNING: ", unmatched,
                                  " diGly run(s) had no matching proteome run - dropped."))

          meta <- meta |>
            dplyr::mutate(
              cell_line = droplevels(cell_line),
              sample_id = paste0(
                stringr::str_replace_all(as.character(cell_line), "[-.\\s]", ""),
                "_R", replicate))

          y.dg <- y.dg[, meta$digly_run];   colnames(y.dg) <- meta$sample_id
          y.tp <- y.tp[, meta$proteome_run]; colnames(y.tp) <- meta$sample_id

          targets <- data.frame(
            cell_line = meta$cell_line,
            replicate = meta$replicate,
            inj_order = meta$inj_order,
            row.names = meta$sample_id)
          y.dg$targets <- y.tp$targets <- targets

          # FASTA -------------------------------------------------------------
          alog(load_log, "\n=== Reading FASTA ===")
          prog$tick("Reading FASTA…")
          fa  <- Biostrings::readAAStringSet(p$fasta)
          acc <- stringr::str_match(
            names(fa), "^(?:[a-z]{2}[|])?([A-Z0-9-]+)")[, 2]
          prot_seq <- setNames(as.character(fa), acc)
          prot_seq <- prot_seq[!duplicated(names(prot_seq))]
          alog(load_log, paste0("  ", length(prot_seq), " FASTA entries."))

          prog$tick("Finalising…")
          alog(load_log, paste0(
            "\nLoaded: ", nrow(meta), " samples | ",
            nrow(y.dg), " diGly prec | ",
            nrow(y.tp), " proteome prec"))

          data_obj(list(
            y.dg               = y.dg,
            y.tp               = y.tp,
            meta               = meta,
            targets            = targets,
            prot_seq           = prot_seq,
            paths              = p,
            cell_lines_present = levels(droplevels(meta$cell_line))
          ))
          status("done")
          sp_load_state("done")

        }, error = function(e) {
          alog(load_log, paste0("ERROR: ", e$message))
          status("error")
          sp_load_state("error")
        }),
        message = function(m) {
          alog(load_log, conditionMessage(m))
          invokeRestart("muffleMessage")
        },
        warning = function(w) {
          alog(load_log, paste0("WARNING: ", conditionMessage(w)))
          invokeRestart("muffleWarning")
        }
      )
    })

    # Outputs -----------------------------------------------------------------
    output$scan_log <- renderPrint({ cat(scan_log()) })
    output$load_log <- renderPrint({ cat(load_log()) })

    output$meta_table <- DT::renderDataTable({
      d <- data_obj(); req(d)
      d$meta |>
        dplyr::select(sample_id, cell_line, replicate,
                      well, inj_order, digly_run, proteome_run)
    }, options = list(scrollX = TRUE, pageLength = 20), rownames = FALSE)

    return(data_obj)
  })
}
