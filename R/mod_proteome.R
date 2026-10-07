# =============================================================================
# mod_proteome -- Proteome integration: degradative vs non-degradative events
#
# Combines gene-level diGly OCCUPANCY (from the DE Analysis tab, or from a
# site-level occupancy DE file) with a gene-level total-PROTEOME differential expression
# table measured WITHOUT proteasome inhibitor (MG-132), to separate:
#
#   non-degradative candidate  occupancy FDR < occ_cut  AND  proteome FDR > null_cut
#                              (occupancy shifts with no protein-level change)
#   same direction (strict)    occupancy FDR < occ_cut  AND  proteome FDR < strict_cut,
#                              logFC_prot and occupancy z of the same sign
#   same direction (relaxed)   as above with proteome FDR < loose_cut
#   opposite direction         occupancy FDR < occ_cut  AND  proteome FDR < loose_cut,
#                              opposite signs
#   intermediate               occupancy significant, proteome FDR between the
#                              relaxed and null cut-offs (not clearly either)
#
# Gene-level occupancy = Stouffer-combined z over all diGly sites of a gene:
#   z_site = sign(t_occ) * qnorm(P_occ / 2, lower.tail = FALSE)
#   z_gene = sum(z_site) / sqrt(n_sites);  BH-adjusted across genes.
#
# One run = one comparison (e.g. one target line vs rest). Runs accumulate, so
# several comparisons can be shown side by side (faceted by comparison label).
# Thresholds are reactive: changing them does not require re-running.
#
# Receives de_obj (from mod_de / restore) and cfg. Optional E3-substrate
# literature table (e.g. literature.E3_H.sapiens.txt) annotates candidates.
# =============================================================================

# ---- core logic (pure functions, no Shiny) ----------------------------------

PI_CLASSES <- c("non-degradative candidate",
                "same direction (strict)",
                "same direction (relaxed)",
                "opposite direction",
                "intermediate",
                "occupancy n.s.")

# Read a delimited table; .csv -> comma, everything else -> tab
pi_read_table <- function(path, n_max = Inf) {
  if (grepl("\\.csv(\\.gz)?$", path, ignore.case = TRUE))
    readr::read_csv(path, n_max = n_max, show_col_types = FALSE, progress = FALSE)
  else
    readr::read_tsv(path, n_max = n_max, show_col_types = FALSE, progress = FALSE)
}

# First candidate name found (case-insensitive) in `cols`, else ""
pi_guess_col <- function(cols, candidates) {
  hit <- cols[match(tolower(candidates), tolower(cols), nomatch = 0L)]
  if (length(hit)) hit[[1]] else ""
}

# Site-level occupancy DE -> gene-level Stouffer occupancy
pi_gene_level_occupancy <- function(df, gene_col = "Genes",
                                    p_col = "P_occ", t_col = "t_occ") {
  miss <- setdiff(c(gene_col, p_col, t_col), names(df))
  if (length(miss))
    stop(sprintf("Occupancy table is missing column(s): %s", paste(miss, collapse = ", ")))

  d <- dplyr::tibble(
    Genes = trimws(as.character(df[[gene_col]])),
    P     = suppressWarnings(as.numeric(df[[p_col]])),
    t     = suppressWarnings(as.numeric(df[[t_col]]))) |>
    dplyr::filter(!is.na(P), !is.na(t), !is.na(Genes), Genes != "") |>
    dplyr::mutate(z_site = stats::qnorm(pmax(P, 1e-300) / 2, lower.tail = FALSE) * sign(t))

  d |>
    dplyr::group_by(Genes) |>
    dplyr::summarise(n_sites = dplyr::n(),
                     occ_z   = sum(z_site) / sqrt(dplyr::n()),
                     .groups = "drop") |>
    dplyr::mutate(occ_p   = 2 * stats::pnorm(-abs(occ_z)),
                  occ_FDR = stats::p.adjust(occ_p, method = "BH"))
}

# Proteome DE -> one row per gene (multi-mapped entries: keep largest |t|)
pi_prepare_proteome <- function(df, gene_col, logfc_col, t_col, fdr_col) {
  cols <- c(gene_col, logfc_col, t_col, fdr_col)
  if (any(!nzchar(cols)))
    stop("Select the gene, logFC, t-statistic and FDR columns of the proteome table.")
  miss <- setdiff(cols, names(df))
  if (length(miss))
    stop(sprintf("Proteome table is missing column(s): %s", paste(miss, collapse = ", ")))

  dplyr::tibble(
    Genes      = trimws(as.character(df[[gene_col]])),
    logFC_prot = suppressWarnings(as.numeric(df[[logfc_col]])),
    t_prot     = suppressWarnings(as.numeric(df[[t_col]])),
    FDR_prot   = suppressWarnings(as.numeric(df[[fdr_col]]))) |>
    dplyr::filter(!is.na(t_prot), !is.na(logFC_prot), !is.na(FDR_prot),
                  !is.na(Genes), Genes != "") |>
    dplyr::group_by(Genes) |>
    dplyr::slice_max(abs(t_prot), n = 1, with_ties = FALSE) |>
    dplyr::ungroup()
}

# Threshold-based classification of every gene with both layers
pi_classify <- function(m, occ_fdr = 0.05, prot_strict = 0.05,
                        prot_loose = 0.25, prot_null = 0.5) {
  prot_loose <- max(prot_loose, prot_strict)
  m |>
    dplyr::mutate(
      same_sign = sign(occ_z) == sign(logFC_prot),
      classification = dplyr::case_when(
        !(occ_FDR < occ_fdr)                        ~ "occupancy n.s.",
        FDR_prot > prot_null                        ~ "non-degradative candidate",
        FDR_prot < prot_strict &  same_sign         ~ "same direction (strict)",
        FDR_prot < prot_loose  &  same_sign         ~ "same direction (relaxed)",
        FDR_prot < prot_loose  & !same_sign         ~ "opposite direction",
        TRUE                                        ~ "intermediate"),
      classification = factor(classification, levels = PI_CLASSES)) |>
    dplyr::select(-same_sign)
}

pi_out_cols <- c("Experiment", "Genes", "n_sites", "occ_z", "occ_FDR",
                 "logFC_prot", "t_prot", "FDR_prot", "classification")

pi_load_e3 <- function(path, species = "H.sapiens") {
  e3 <- pi_read_table(path)
  need <- c("Gene Symbol (E3)", "Gene Symbol (Substrate)")
  miss <- setdiff(need, names(e3))
  if (length(miss))
    stop(sprintf("E3 table is missing column(s): %s", paste(miss, collapse = ", ")))
  names(e3)[names(e3) == "Gene Symbol (E3)"]       <- "E3"
  names(e3)[names(e3) == "Gene Symbol (Substrate)"] <- "Substrate"
  if (nzchar(species) && "species" %in% names(e3))
    e3 <- e3[e3[["species"]] %in% species, , drop = FALSE]
  dplyr::distinct(e3, E3, Substrate, .keep_all = TRUE)
}

pi_e3_hits <- function(e3, cand) {
  if (is.null(e3) || nrow(cand) == 0) return(NULL)
  e3 |>
    dplyr::filter(Substrate %in% cand$Genes) |>
    dplyr::inner_join(cand, by = c("Substrate" = "Genes"),
                      relationship = "many-to-many")
}

# ---- plots (ggplot2) --------------------------------------------------------

pi_empty_plot <- function(msg) {
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = 0, label = msg, colour = "grey40", size = 5) +
    ggplot2::theme_void()
}

pi_theme <- function() {
  ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "white", color = "white"),
      strip.text.y     = ggplot2::element_text(face = "bold", angle = 0),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position  = "right",
      plot.margin      = ggplot2::margin(10, 15, 10, 10))
}

# Keep the top-n genes (smallest occupancy FDR) per comparison
pi_top_n <- function(d, n) {
  if (nrow(d) == 0 || is.na(n) || n < 1) return(d)
  d |>
    dplyr::group_by(Experiment) |>
    dplyr::slice_min(occ_FDR, n = n, with_ties = FALSE) |>
    dplyr::ungroup()
}

# Lollipop of gene-level occupancy z (non-degradative candidates)
pi_plot_nondeg <- function(d, e3_genes = character(0)) {
  if (nrow(d) == 0) return(pi_empty_plot("No non-degradative candidates at the current thresholds"))
  d <- d |>
    dplyr::mutate(
      direction   = ifelse(occ_z > 0, "Occupancy up", "Occupancy down"),
      neglog10FDR = -log10(pmax(occ_FDR, 1e-14)),
      lab         = ifelse(Genes %in% e3_genes, paste0(Genes, " *"), Genes)) |>
    dplyr::arrange(Experiment, occ_z) |>
    dplyr::mutate(row_id = factor(paste(Experiment, lab, sep = " || "),
                                  levels = unique(paste(Experiment, lab, sep = " || "))))

  cap <- if (length(e3_genes))
    "* documented substrate of a literature E3 ligase" else NULL

  ggplot2::ggplot(d, ggplot2::aes(x = occ_z, y = row_id)) +
    ggplot2::geom_vline(xintercept = 0, color = "grey60", linewidth = 0.4) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = occ_z, yend = row_id, color = direction),
                          linewidth = 0.7) +
    ggplot2::geom_point(ggplot2::aes(size = neglog10FDR, color = direction)) +
    ggplot2::scale_color_manual(values = c("Occupancy up" = "#e66101",
                                           "Occupancy down" = "#5e3c99"), name = NULL) +
    ggplot2::scale_size_continuous(name = "-log10(occupancy FDR)", range = c(2, 7)) +
    ggplot2::scale_y_discrete(labels = function(x) sub(".* \\|\\| ", "", x)) +
    ggplot2::facet_grid(Experiment ~ ., scales = "free_y", space = "free_y") +
    ggplot2::labs(x = "Gene-level occupancy Stouffer z-score", y = NULL, caption = cap) +
    pi_theme() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 9))
}

# Two-column dot plot: occupancy z vs proteome t (same / opposite direction)
pi_plot_direction <- function(d) {
  if (nrow(d) == 0) return(pi_empty_plot("No genes in this category at the current thresholds"))
  d <- d |>
    dplyr::arrange(Experiment, occ_z) |>
    dplyr::mutate(row_id = paste(Experiment, Genes, sep = " || "),
                  row_id = factor(row_id, levels = rev(unique(row_id))))

  lyr <- c("Occupancy\n(z-score)", "Proteome\n(t-statistic)")
  long <- dplyr::bind_rows(
    d |> dplyr::transmute(Experiment, row_id, layer = lyr[1], value = occ_z,  FDR = occ_FDR),
    d |> dplyr::transmute(Experiment, row_id, layer = lyr[2], value = t_prot, FDR = FDR_prot)) |>
    dplyr::mutate(layer = factor(layer, levels = lyr),
                  neglog10FDR = -log10(pmax(FDR, 1e-6)))

  ggplot2::ggplot(long, ggplot2::aes(x = layer, y = row_id)) +
    ggplot2::geom_point(ggplot2::aes(size = neglog10FDR, fill = value),
                        shape = 21, color = "grey30", stroke = 0.3) +
    ggplot2::scale_fill_gradient2(low = "#5e3c99", mid = "#f7f7f7", high = "#e66101",
                                  midpoint = 0, name = "z-score /\nt-statistic") +
    ggplot2::scale_size_continuous(name = "FDR", range = c(3, 10),
                                   breaks = c(0.7, 1.3, 2, 3),
                                   labels = c("0.2", "0.05", "0.01", "0.001")) +
    ggplot2::scale_y_discrete(labels = function(x) sub(".* \\|\\| ", "", x)) +
    ggplot2::facet_grid(Experiment ~ ., scales = "free_y", space = "free_y") +
    ggplot2::labs(x = NULL, y = NULL) +
    pi_theme() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10),
                   axis.text.x = ggplot2::element_text(size = 9.5, face = "bold"))
}

pi_plot_height <- function(n_rows, n_facets, per_row = 19) {
  as.integer(max(320, 110 + per_row * n_rows + 16 * n_facets))
}

# ---- UI ---------------------------------------------------------------------

mod_proteome_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h4("Proteome Integration \u2014 degradative vs non-degradative events"),

    # ── Inputs ────────────────────────────────────────────────────────────────
    card(
      card_body(
        tags$p(class = "text-muted small mb-2",
          "Combines gene-level diGly occupancy (Stouffer z over all sites of a gene) with a ",
          "gene-level total-proteome DE table measured ", tags$b("without MG-132"), ". ",
          "Occupancy shifts with no protein change are candidate non-degradative events; ",
          "shifts that track protein abundance are candidate degradative/concordant events. ",
          tags$b("The proteome logFC must be oriented like the occupancy contrast "),
          "(target vs reference) or the same/opposite-direction classes are inverted."),

        fluidRow(
          # 1. Occupancy -------------------------------------------------------
          column(4,
            tags$h6("1. Occupancy (diGly)"),
            radioButtons(ns("occ_source"), NULL,
                         choices  = c("Current session (DE Analysis tab)" = "session",
                                      "Site-level Occupancy DE file"                = "file"),
                         selected = "session"),
            conditionalPanel(
              condition = sprintf("input['%s'] == 'file'", ns("occ_source")),
              hybrid_path_ui(ns, "occ_file", "Site-level Occupancy DE table (.tsv / .csv)",
                             type = "file",
                             placeholder = "diGly_..._PSP_annotated.tsv"),
              selectInput(ns("occ_gene"), "Gene column",     choices = c("(select)" = "")),
              selectInput(ns("occ_p"),    "P-value column",  choices = c("(select)" = "")),
              selectInput(ns("occ_t"),    "t-statistic column", choices = c("(select)" = ""))
            ),
            conditionalPanel(
              condition = sprintf("input['%s'] == 'session'", ns("occ_source")),
              helpText("Uses the occupancy layer (Genes, P_occ, t_occ) of the current DE results.")
            )
          ),

          # 2. Proteome --------------------------------------------------------
          column(4,
            tags$h6("2. Proteome DE (no MG-132)"),
            hybrid_path_ui(ns, "prot_file", "Proteome DE table (.tsv / .csv)",
                           type = "file",
                           placeholder = "DE_<line>_vs_Background_proteins.tsv"),
            selectInput(ns("prot_gene"),  "Gene column",         choices = c("(select)" = "")),
            selectInput(ns("prot_logfc"), "logFC column",        choices = c("(select)" = "")),
            selectInput(ns("prot_t"),     "t-statistic column",  choices = c("(select)" = "")),
            selectInput(ns("prot_fdr"),   "FDR column",          choices = c("(select)" = ""))
          ),

          # 3. Labels / E3 -----------------------------------------------------
          column(4,
            tags$h6("3. Comparison & annotation"),
            textInput(ns("label"), "Comparison label",
                      placeholder = "e.g. AsPC-1 (target vs rest)"),
            helpText("Shown as the facet/Experiment name. Re-using a label replaces that comparison."),
            hybrid_path_ui(ns, "e3_file", "E3-substrate literature table (optional)",
                           type = "file",
                           placeholder = "literature.E3_H.sapiens.txt"),
            textInput(ns("e3_species"), "E3 table species filter (blank = none)",
                      value = "H.sapiens")
          )
        ),

        tags$div(style = "display:flex; align-items:center; gap:12px; margin-top:8px;",
          actionButton(ns("run"),   "Add / update comparison", class = "btn-primary"),
          actionButton(ns("clear"), "Clear all comparisons",   class = "btn-outline-secondary"),
          spinner_ui(ns("sp_prot"))
        ),
        verbatimTextOutput(ns("log"))
      )
    ),

    # ── Results ───────────────────────────────────────────────────────────────
    conditionalPanel(
      condition = sprintf("output['%s']", ns("done")),

      h4("Thresholds & comparisons"),
      card(
        card_body(
          fluidRow(
            column(2, numericInput(ns("occ_fdr"),     "Occupancy FDR <",          0.05, min = 0, max = 1, step = 0.01)),
            column(2, numericInput(ns("prot_strict"), "Proteome FDR < (strict)",  0.05, min = 0, max = 1, step = 0.01)),
            column(2, numericInput(ns("prot_loose"),  "Proteome FDR < (relaxed)", 0.25, min = 0, max = 1, step = 0.05)),
            column(2, numericInput(ns("prot_null"),   "Proteome FDR > (null)",    0.5,  min = 0, max = 1, step = 0.05)),
            column(2, numericInput(ns("top_n"),       "Max genes / comparison in plots", 50, min = 1, step = 5)),
            column(2, checkboxGroupInput(ns("exp_show"), "Comparisons shown", choices = NULL))
          ),
          helpText(
            tags$b("non-degradative candidate:"), " occupancy FDR < cut-off and proteome FDR > null. ",
            tags$b("same / opposite direction:"), " occupancy and proteome both significant ",
            "(strict / relaxed proteome FDR), with equal / opposite sign of occupancy z and proteome logFC. ",
            tags$b("intermediate:"), " occupancy significant, proteome neither clearly changed nor unchanged."),
          verbatimTextOutput(ns("summary"))
        )
      ),

      tabsetPanel(
        id = ns("tabs"),

        tabPanel("Non-degradative candidates",
          tags$div(style = "margin-top:12px;",
            downloadButton(ns("dl_nd_tbl"),  "Table (.csv)"),
            downloadButton(ns("dl_nd_e3"),   "E3 pairs (.csv)"),
            downloadButton(ns("dl_nd_svg"),  "Figure (.svg)")),
          plotOutput(ns("plot_nd"), height = "auto"),
          h5("Candidates"),
          DT::dataTableOutput(ns("tbl_nd")),
          h5("Documented E3 ligases (literature)"),
          DT::dataTableOutput(ns("tbl_e3"))
        ),

        tabPanel("Same / opposite direction",
          tags$div(style = "margin-top:12px;",
            radioButtons(ns("dir_view"), NULL, inline = TRUE,
                         choices = c("Same direction (strict + relaxed)" = "same",
                                     "Opposite direction"                = "opp",
                                     "Both"                              = "both"),
                         selected = "same"),
            downloadButton(ns("dl_dir_tbl"), "Table (.csv)"),
            downloadButton(ns("dl_dir_svg"), "Figure (.svg)")),
          plotOutput(ns("plot_dir"), height = "auto"),
          h5("Genes"),
          DT::dataTableOutput(ns("tbl_dir"))
        ),

        tabPanel("All genes",
          tags$div(style = "margin-top:12px;",
            downloadButton(ns("dl_all_tbl"), "Table (.csv)")),
          DT::dataTableOutput(ns("tbl_all"))
        )
      ),

      h4("Save to disk"),
      card(
        card_body(
          tags$p(class = "text-muted small mb-2",
            "Writes every table (.csv) and figure (.svg) to ",
            tags$code("<Results directory>/Proteome_integration/"), "."),
          actionButton(ns("save_all"), "Save all tables + figures", class = "btn-success"),
          verbatimTextOutput(ns("save_log"))
        )
      )
    )
  )
}

# ---- server -----------------------------------------------------------------

mod_proteome_server <- function(id, de_obj, cfg, status = reactiveVal("idle")) {
  moduleServer(id, function(input, output, session) {

    `%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && nzchar(a[[1]])) a[[1]] else b

    log_msg  <- reactiveVal("")
    save_log <- reactiveVal("")
    store    <- reactiveVal(list())     # label -> list(merged, n_occ, n_prot, ...)
    e3_db    <- reactiveVal(NULL)
    sp_state <- reactiveVal("idle")
    spinner_server("sp_prot", sp_state)

    append_log <- function(...) log_msg(paste0(log_msg(), paste(...), "\n"))

    output$done <- reactive(length(store()) > 0)
    outputOptions(output, "done", suspendWhenHidden = FALSE)

    # ── File pickers + header-driven column auto-detection ────────────────────
    roots <- c(Computer = if (.Platform$OS.type == "windows") "C:/" else "/")
    wire_file <- function(key, filetypes, on_path = NULL) {
      btn <- paste0("btn_", key); txt <- paste0("txt_", key)
      shinyFiles::shinyFileChoose(input, btn, roots = roots, filetypes = filetypes)
      observeEvent(input[[btn]], {
        p <- tryCatch(shinyFiles::parseFilePaths(roots, input[[btn]])$datapath[[1]],
                      error = function(e) NULL)
        if (!is.null(p) && nzchar(p)) updateTextInput(session, txt, value = p)
      }, ignoreInit = TRUE, ignoreNULL = TRUE)

      output[[paste0("disp_", key)]] <- renderUI({
        p <- trimws(input[[txt]] %||% "")
        if (!nzchar(p)) return(NULL)
        if (file.exists(p))
          tags$small(style = "color:#28a745;", "\u2714 found")
        else
          tags$small(style = "color:#dc3545;", "\u2716 file not found")
      })

      if (!is.null(on_path))
        observeEvent(input[[txt]], {
          p <- trimws(input[[txt]] %||% "")
          if (nzchar(p) && file.exists(p))
            tryCatch(on_path(names(pi_read_table(p, n_max = 0))),
                     error = function(e) showNotification(
                       paste("Could not read header:", e$message), type = "warning"))
        }, ignoreInit = TRUE)
    }

    set_choices <- function(inputId, cols, candidates) {
      updateSelectInput(session, inputId,
        choices  = c("(select)" = "", cols),
        selected = pi_guess_col(cols, candidates))
    }
    gene_cands <- c("Genes.x", "Genes", "Gene", "Gene.Symbol", "Gene_Symbol",
                    "Gene.Names", "GeneName", "gene_name", "PG.Genes", "SYMBOL")

    wire_file("occ_file",  c("tsv", "txt", "csv", "gz"), function(cols) {
      set_choices("occ_gene", cols, gene_cands)
      set_choices("occ_p",    cols, c("P_occ", "P.Value", "P_value", "pvalue", "PValue"))
      set_choices("occ_t",    cols, c("t_occ", "t", "t_stat", "statistic"))
    })
    wire_file("prot_file", c("tsv", "txt", "csv", "gz"), function(cols) {
      set_choices("prot_gene",  cols, gene_cands)
      set_choices("prot_logfc", cols, c("logFC", "log2FC", "log2FoldChange", "logFC_prot"))
      set_choices("prot_t",     cols, c("t", "t_prot", "t.value", "statistic", "stat"))
      set_choices("prot_fdr",   cols, c("FDR", "adj.P.Val", "adj.P.Value", "padj",
                                         "qvalue", "q.value", "FDR_prot", "adj_p"))
    })
    wire_file("e3_file",   c("tsv", "txt", "csv", "gz"))

    # Default comparison label = configured target line
    observe({
      c <- cfg()
      tl <- if (is.null(c)) "" else (c$target_line %||% "")
      if (nzchar(tl) && !nzchar(isolate(input$label %||% "")))
        updateTextInput(session, "label", value = tl)
    })

    # ── Run: build gene-level table for one comparison ────────────────────────
    observeEvent(input$run, {
      log_msg(""); status("running"); sp_state("running")

      withCallingHandlers(
        tryCatch({
          label <- trimws(input$label %||% "")
          if (!nzchar(label)) stop("Enter a comparison label.")

          # -- occupancy --------------------------------------------------------
          append_log("=== Occupancy (gene level) ===")
          if (input$occ_source == "session") {
            de <- de_obj()
            if (is.null(de))
              stop("No DE results in this session. Run the DE Analysis tab (or restore a session), ",
                   "or choose 'Site-level Occupancy DE file'.")
            occ_df <- de$res
            g <- "Genes"; p <- "P_occ"; tt <- "t_occ"
            append_log(sprintf("  Source: session DE results (%d sites)", nrow(occ_df)))
          } else {
            path <- trimws(input$txt_occ_file %||% "")
            if (!nzchar(path) || !file.exists(path)) stop("Occupancy file not found.")
            occ_df <- pi_read_table(path)
            g <- input$occ_gene; p <- input$occ_p; tt <- input$occ_t
            if (!all(nzchar(c(g, p, tt)))) stop("Select the gene, P-value and t-statistic columns of the occupancy file.")
            append_log(sprintf("  Source: %s (%d rows)", basename(path), nrow(occ_df)))
          }
          occ_gene <- pi_gene_level_occupancy(occ_df, g, p, tt)
          append_log(sprintf("  Genes with occupancy: %d | significant (FDR<%.2f): %d",
                             nrow(occ_gene), input$occ_fdr %||% 0.05,
                             sum(occ_gene$occ_FDR < (input$occ_fdr %||% 0.05))))

          # -- proteome ---------------------------------------------------------
          append_log("\n=== Proteome DE ===")
          ppath <- trimws(input$txt_prot_file %||% "")
          if (!nzchar(ppath) || !file.exists(ppath)) stop("Proteome DE file not found.")
          prot_raw <- pi_read_table(ppath)
          prot <- pi_prepare_proteome(prot_raw, input$prot_gene, input$prot_logfc,
                                      input$prot_t, input$prot_fdr)
          append_log(sprintf("  %s: %d rows -> %d genes", basename(ppath),
                             nrow(prot_raw), nrow(prot)))

          # -- merge ------------------------------------------------------------
          merged <- dplyr::inner_join(occ_gene, prot, by = "Genes")
          append_log(sprintf("\n  Genes in both layers: %d", nrow(merged)))
          if (nrow(merged) == 0)
            stop("No overlapping gene symbols between occupancy and proteome tables. ",
                 "Check the gene columns / identifier type.")
          if (nrow(merged) < 0.2 * min(nrow(occ_gene), nrow(prot)))
            append_log("  [!] Low overlap (<20%) -- check that both tables use the same gene symbols.")

          merged$Experiment <- label
          s <- store()
          if (label %in% names(s)) append_log(sprintf("  Replaced existing comparison '%s'.", label))
          s[[label]] <- list(merged = merged, n_occ = nrow(occ_gene), n_prot = nrow(prot))
          store(s)
          updateCheckboxGroupInput(session, "exp_show",
                                   choices = names(s), selected = names(s))

          # -- optional E3 table -----------------------------------------------
          e3path <- trimws(input$txt_e3_file %||% "")
          if (nzchar(e3path)) {
            if (!file.exists(e3path)) stop("E3 literature file not found.")
            e3_db(pi_load_e3(e3path, trimws(input$e3_species %||% "")))
            append_log(sprintf("\n  E3 table loaded: %d E3-substrate pairs", nrow(e3_db())))
          }

          append_log("\n[OK] Comparison added.")
          status("done"); sp_state("done")
        }, error = function(e) {
          append_log(paste0("[ERR] ERROR: ", conditionMessage(e)))
          status("error"); sp_state("error")
        }),
        message = function(m) { append_log(conditionMessage(m)); invokeRestart("muffleMessage") },
        warning = function(w) { append_log(paste("[!]", conditionMessage(w))); invokeRestart("muffleWarning") }
      )
    })

    observeEvent(input$clear, {
      store(list()); e3_db(NULL); log_msg(""); status("idle"); sp_state("idle")
      updateCheckboxGroupInput(session, "exp_show", choices = character(0))
    })

    # ── Reactive tables ────────────────────────────────────────────────────────
    num <- function(x, default) if (is.null(x) || is.na(x)) default else x

    classified <- reactive({
      s <- store(); req(length(s) > 0, input$exp_show)
      m <- dplyr::bind_rows(lapply(s[intersect(names(s), input$exp_show)], `[[`, "merged"))
      req(nrow(m) > 0)
      m$Experiment <- factor(m$Experiment, levels = names(s))
      pi_classify(m,
                  occ_fdr     = num(input$occ_fdr,     0.05),
                  prot_strict = num(input$prot_strict, 0.05),
                  prot_loose  = num(input$prot_loose,  0.25),
                  prot_null   = num(input$prot_null,   0.5)) |>
        dplyr::select(dplyr::all_of(pi_out_cols)) |>
        dplyr::arrange(Experiment, occ_FDR)
    })

    nondeg <- reactive(dplyr::filter(classified(), classification == "non-degradative candidate"))
    same_d <- reactive(dplyr::filter(classified(),
                         classification %in% c("same direction (strict)", "same direction (relaxed)")))
    same_s <- reactive(dplyr::filter(classified(), classification == "same direction (strict)"))
    opp_d  <- reactive(dplyr::filter(classified(), classification == "opposite direction"))
    dir_d  <- reactive(switch(input$dir_view %||% "same",
                              same = same_d(), opp = opp_d(), both = dplyr::bind_rows(same_d(), opp_d())))

    e3_hits <- reactive(pi_e3_hits(e3_db(), nondeg()))

    nd_plot_df  <- reactive(pi_top_n(nondeg(), num(input$top_n, 50)))
    dir_plot_df <- reactive(pi_top_n(dir_d(),  num(input$top_n, 50)))

    n_facets <- function(d) dplyr::n_distinct(d$Experiment)

    # ── Outputs ────────────────────────────────────────────────────────────────
    output$log <- renderPrint({ cat(log_msg()) })

    output$summary <- renderPrint({
      s <- store(); req(length(s) > 0)
      cat("Comparisons (genes with both layers):\n")
      print(as.data.frame(dplyr::bind_rows(lapply(names(s), function(k)
        data.frame(Experiment = k, genes_occ = s[[k]]$n_occ, genes_prot = s[[k]]$n_prot,
                   genes_both = nrow(s[[k]]$merged))))), row.names = FALSE)
      cat("\nGenes per class:\n")
      tab <- table(Experiment = classified()$Experiment, class = classified()$classification)
      print(tab)
    })

    # -- non-degradative
    output$plot_nd <- renderPlot({
      d <- nd_plot_df()
      pi_plot_nondeg(d, e3_genes = unique(e3_hits()$Substrate))
    }, height = function() tryCatch(
         pi_plot_height(nrow(nd_plot_df()), n_facets(nd_plot_df())),
         error = function(e) 320L))

    round_tbl <- function(d) dplyr::mutate(d, dplyr::across(where(is.numeric), ~ signif(.x, 4)))
    dt_opts <- list(scrollX = TRUE, pageLength = 15)

    output$tbl_nd  <- DT::renderDataTable(round_tbl(nondeg()), options = dt_opts,
                                          filter = "top", rownames = FALSE)
    output$tbl_e3  <- DT::renderDataTable({
      h <- e3_hits()
      if (is.null(h)) {
        msg <- if (is.null(e3_db())) "Provide an E3 literature table and re-run to see documented E3 ligases."
               else "No documented E3 ligases for the current candidates."
        return(data.frame(Note = msg))
      }
      round_tbl(h)
    }, options = list(scrollX = TRUE, pageLength = 10), rownames = FALSE)

    # -- same / opposite direction
    output$plot_dir <- renderPlot(pi_plot_direction(dir_plot_df()),
      height = function() tryCatch(
        pi_plot_height(nrow(dir_plot_df()), n_facets(dir_plot_df()), per_row = 24),
        error = function(e) 320L))
    output$tbl_dir <- DT::renderDataTable(round_tbl(dir_d()), options = dt_opts,
                                          filter = "top", rownames = FALSE)

    # -- all
    output$tbl_all <- DT::renderDataTable(round_tbl(classified()), options = dt_opts,
                                          filter = "top", rownames = FALSE)

    # ── Figure writers (shared by downloads and save-to-disk) ─────────────────
    save_svg <- function(file, plot, n_rows, n_fac, per_row = 19) {
      h <- pi_plot_height(n_rows, n_fac, per_row) / 96
      ggplot2::ggsave(file, plot, width = 9.5, height = h,
                      device = svglite::svglite, limitsize = FALSE)
    }
    fig_nd  <- function() { d <- nd_plot_df()
      list(p = pi_plot_nondeg(d, unique(e3_hits()$Substrate)), n = nrow(d), f = n_facets(d)) }
    fig_dir <- function() { d <- dir_plot_df()
      list(p = pi_plot_direction(d), n = nrow(d), f = n_facets(d)) }

    output$dl_nd_tbl <- downloadHandler(
      filename = function() "nondegradative_candidates.csv",
      content  = function(file) readr::write_csv(nondeg(), file))
    output$dl_nd_e3 <- downloadHandler(
      filename = function() "nondegradative_candidates_E3_literature.csv",
      content  = function(file) {
        h <- e3_hits()
        readr::write_csv(if (is.null(h)) data.frame() else h, file)
      })
    output$dl_nd_svg <- downloadHandler(
      filename = function() "nondegradative_candidates.svg",
      content  = function(file) { f <- fig_nd(); save_svg(file, f$p, f$n, f$f) })

    output$dl_dir_tbl <- downloadHandler(
      filename = function() sprintf("%s_direction_candidates.csv", input$dir_view %||% "same"),
      content  = function(file) readr::write_csv(dir_d(), file))
    output$dl_dir_svg <- downloadHandler(
      filename = function() sprintf("%s_direction_candidates.svg", input$dir_view %||% "same"),
      content  = function(file) { f <- fig_dir(); save_svg(file, f$p, f$n, f$f, 24) })

    output$dl_all_tbl <- downloadHandler(
      filename = function() "proteome_integration_all_genes.csv",
      content  = function(file) readr::write_csv(classified(), file))

    # ── Save everything to <results_dir>/Proteome_integration ─────────────────
    observeEvent(input$save_all, {
      c <- cfg(); req(classified())
      rd <- if (is.null(c)) "" else (c$results_dir %||% "")
      if (!nzchar(rd)) {
        save_log("Results directory not set (Data tab). Use the download buttons instead.")
        return()
      }
      out <- file.path(rd, "Proteome_integration")
      dir.create(out, showWarnings = FALSE, recursive = TRUE)
      tryCatch({
        w <- function(d, name) {
          readr::write_csv(d, file.path(out, name))
          paste0("  [OK] ", name, sprintf(" (%d rows)", nrow(d)))
        }
        msgs <- c(
          w(classified(), "proteome_integration_all_genes.csv"),
          w(nondeg(),     "nondegradative_candidates.csv"),
          w(same_d(),     "same_direction_candidates.csv"),
          w(same_s(),     "same_direction_strict.csv"),
          w(opp_d(),      "opposite_direction_candidates.csv"))
        h <- e3_hits()
        if (!is.null(h)) msgs <- c(msgs, w(h, "nondegradative_candidates_E3_literature.csv"))

        f <- fig_nd();  save_svg(file.path(out, "nondegradative_candidates.svg"), f$p, f$n, f$f)
        f <- fig_dir(); save_svg(file.path(out, "direction_candidates.svg"),      f$p, f$n, f$f, 24)
        msgs <- c(msgs, "  [OK] nondegradative_candidates.svg", "  [OK] direction_candidates.svg")
        save_log(paste(c(paste("Saved to", out), msgs), collapse = "\n"))
      }, error = function(e) save_log(paste("[ERR]", conditionMessage(e))))
    })
    output$save_log <- renderPrint({ cat(save_log()) })

    invisible(classified)
  })
}
