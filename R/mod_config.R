# =============================================================================
# mod_config  Analysis configuration panel
# Paths and run-name token map live in mod_data (Data tab).
# This tab covers: experimental design, filters, thresholds, highlighted genes.
# Cell-line dropdowns are populated from the loaded data object.
# =============================================================================

mod_config_ui <- function(id) {
  ns <- NS(id)
  tagList(

    h4("Experimental Design"),
    card(
      card_body(
        fluidRow(
          column(4,
            selectInput(ns("comparison_mode"), "Comparison mode",
                        choices  = c("One vs. rest (mean of all others)" = "vs_rest",
                                     "One vs. single control"             = "vs_ctrl"),
                        selected = "vs_rest"),
            helpText(tags$b("vs rest:"), "target vs average of all other lines.",
                     tags$br(),
                     tags$b("vs control:"), "target vs one designated control line.")
          ),
          column(4,
            selectInput(ns("target_line"), "Target cell line", choices = NULL),
            helpText("Numerator of the contrast (logFC > 0 = up in target).")
          ),
          column(4,
            conditionalPanel(
              condition = sprintf("input['%s'] == 'vs_ctrl'", ns("comparison_mode")),
              selectInput(ns("ctrl_line"), "Control cell line", choices = NULL),
              helpText("Single denominator line (vs_ctrl mode only).")
            )
          )
        ),
        helpText(icon("info-circle"),
                 "Dropdowns are populated once runs are loaded in the Data tab. ",
                 "Run exclusions: deselect unwanted runs there before loading.")
      )
    ),

    h4("Identification Filters"),
    card(
      card_body(
        fluidRow(
          column(3, numericInput(ns("q_value"),     "Precursor q-value",    0.01, 0, 1,   0.001)),
          column(3, numericInput(ns("ptm_q_value"), "PTM q-value",          0.05, 0, 1,   0.01)),
          column(3, numericInput(ns("site_conf"),   "Site confidence (≥)",  0.75, 0, 1,   0.05)),
          column(3,
            numericInput(ns("min_obs"), "Min obs / line (post-DPC-Quant)", 2, 1, 10, 1),
            helpText("Applied after DPC-Quant in the Occupancy tab. ",
                     "A site passes if the target OR any reference line has ",
                     "≥ N non-zero observation counts in the quantified site matrix. ",
                     tags$b("Different from the precursor missingness filter"),
                     " (Filtering tab), which acts on raw precursor intensities ",
                     "before site mapping."))
        ),
        fluidRow(
          column(4, checkboxInput(ns("drop_cterm"),     "Drop C-terminal GG",  TRUE)),
          column(4, checkboxInput(ns("drop_ubiquitin"), "Drop ubiquitin genes", TRUE))
        )
      )
    ),

    h4("Statistical Thresholds & Normalization"),
    card(
      card_body(
        fluidRow(
          column(3, numericInput(ns("fdr_threshold"),   "FDR threshold",   0.05, 0, 1,   0.01)),
          column(3, numericInput(ns("logfc_threshold"), "logFC threshold", 1.0,  0, Inf, 0.25)),
          column(3, selectInput(ns("digly_norm"),   "diGly normalization",
                                choices = c("median", "none"), selected = "median")),
          column(3, selectInput(ns("proteome_norm"), "Proteome normalization",
                                choices = c("median", "none"), selected = "median"))
        ),
        fluidRow(
          column(4,
            checkboxInput(ns("use_treat"), "Use treat() (vs eBayes)", TRUE),
            helpText(
              tags$b("treat()"), " tests whether |logFC| > threshold (recommended). ",
              tags$b("eBayes()"), " tests whether logFC ≠ 0 — any change, however small, ",
              "can be significant. ",
              tags$br(),
              "treat() builds the fold-change requirement into the null hypothesis: ",
              tags$i("H₀: |logFC| ≤ threshold"),
              ". This is more stringent than filtering eBayes results by logFC after the fact, ",
              "because it accounts for uncertainty in the logFC estimate. ",
              tags$br(),
              tags$b("For diGly occupancy data, treat() is strongly recommended"),
              " — occupancy values combine noise from two DPC-Quant estimates and ",
              "you want sites with robust, biologically meaningful changes."
            )
          )
        )
      )
    ),

    h4("Selected / Highlighted Genes"),
    card(
      card_body(
        textAreaInput(ns("selected_genes_raw"),
                      "Gene symbols to highlight in volcano (one per line)",
                      value = "", rows = 6),
        helpText("Leave blank for no highlights.")
      )
    ),

    div(class = "mt-3",
        actionButton(ns("apply"), "Apply Configuration",
                     class = "btn-primary btn-lg"),
        verbatimTextOutput(ns("cfg_summary"))
    )
  )
}


mod_config_server <- function(id, data_obj) {
  moduleServer(id, function(input, output, session) {

    split_lines <- function(x) Filter(nchar, trimws(strsplit(x, "\n")[[1]]))

    # ── Populate dropdowns from loaded data ───────────────────────────────
    observeEvent(data_obj(), {
      d <- data_obj(); req(d)
      lines <- d$cell_lines_present
      updateSelectInput(session, "target_line", choices = lines, selected = lines[1])
      sel_ctrl <- if (length(lines) > 1) lines[2] else lines[1]
      updateSelectInput(session, "ctrl_line",   choices = lines, selected = sel_ctrl)
    })

    observeEvent(input$target_line, {
      d <- data_obj(); req(d)
      other <- setdiff(d$cell_lines_present, input$target_line)
      cur   <- isolate(input$ctrl_line)
      sel   <- if (!is.null(cur) && cur %in% other) cur else other[1]
      updateSelectInput(session, "ctrl_line", choices = other, selected = sel)
    })

    # ── cfg reactive (fires on Apply) ─────────────────────────────────────
    cfg <- eventReactive(input$apply, {
      d <- data_obj()

      # Pull paths from the loaded data object (stored there at load time)
      p <- if (!is.null(d)) d$paths else list(
        digly_parquet = "", proteome_parquet = "", fasta = "",
        proteome_map_csv = "", results_dir = "", figures_dir = "", limpa_dir = "")

      all_lines <- if (!is.null(d)) d$cell_lines_present else character(0)
      target    <- input$target_line
      mode      <- input$comparison_mode
      ctrl      <- if (mode == "vs_ctrl") input$ctrl_line else NULL
      ref_lines <- if (mode == "vs_ctrl") ctrl else setdiff(all_lines, target)

      if (length(all_lines) == 0)
        showNotification("Load data first (Data tab) before applying configuration.",
                         type = "warning", duration = 6)

      list(
        # paths (forwarded from data)
        digly_parquet    = p$digly_parquet,
        proteome_parquet = p$proteome_parquet,
        fasta            = p$fasta,
        proteome_map_csv = p$proteome_map_csv,
        results_dir      = p$results_dir,
        figures_dir      = p$figures_dir,
        limpa_dir        = p$limpa_dir,
        # design
        target_line      = target,
        ctrl_line        = ctrl,
        comparison_mode  = mode,
        reference_lines  = ref_lines,
        all_lines        = all_lines,
        line_map         = character(0),   # token map lives in mod_data
        # filters
        q_value          = input$q_value,
        ptm_q_value      = input$ptm_q_value,
        site_conf        = input$site_conf,
        min_obs_per_line = input$min_obs,
        drop_cterm_gg    = input$drop_cterm,
        drop_ubiquitin   = input$drop_ubiquitin,
        gg_unimod        = "UniMod:121",
        # stats
        fdr_threshold    = input$fdr_threshold,
        logfc_threshold  = input$logfc_threshold,
        digly_norm       = input$digly_norm,
        proteome_norm    = input$proteome_norm,
        use_treat        = input$use_treat,
        selected_genes   = split_lines(input$selected_genes_raw)
      )
    }, ignoreNULL = FALSE)

    output$cfg_summary <- renderPrint({
      req(cfg())
      c <- cfg()
      cat(sprintf(paste(
        "Mode     : %s",
        "Target   : %s",
        "Ref      : %s",
        "FDR      : %.3f | logFC: %.2f | treat: %s",
        "diGly    : %s",
        "Proteome : %s",
        sep = "\n"),
        c$comparison_mode,
        c$target_line %||% "(not set)",
        paste(c$reference_lines, collapse = ", "),
        c$fdr_threshold, c$logfc_threshold, c$use_treat,
        basename(c$digly_parquet    %||% "(not loaded)"),
        basename(c$proteome_parquet %||% "(not loaded)")
      ), "\n")
    })

    `%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && nzchar(a[[1]])) a[[1]] else b

    return(cfg)
  })
}
