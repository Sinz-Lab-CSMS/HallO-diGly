# =============================================================================
# diGly limpa Occupancy Analysis -- Shiny App
# Tab order: Data → Configuration → Filtering → Normalization → DE → Plots → Export
#            → PSP → GSEA .rnk → Proteome integration
# =============================================================================

library(shiny)
library(bslib)
library(shinyFiles)

# Increase max upload to 2 GB for large RDS session files
options(shiny.maxRequestSize = 2 * 1024^3)

source("R/progress_helpers.R")  # spinner + progress bar helpers
source("R/mod_data.R")          # includes hybrid_path_ui helper
source("R/mod_config.R")
source("R/mod_filtering.R")
source("R/mod_normalization.R")
source("R/mod_occupancy.R")
source("R/mod_de.R")
source("R/mod_plots.R")
source("R/mod_export.R")
source("R/mod_psp.R")
source("R/mod_rnk.R")
source("R/mod_proteome.R")
source("R/mod_restore.R")

# Status indicator UI/server helper ------------------------------------------
status_ui <- function(id) {
  ns <- NS(id)
  tags$span(
    id    = ns("wrap"),
    style = "display:inline-flex; align-items:center; gap:6px;",
    uiOutput(ns("icon"),  inline = TRUE),
    uiOutput(ns("label"), inline = TRUE)
  )
}

status_server <- function(id, state) {
  moduleServer(id, function(input, output, session) {
    output$icon <- renderUI({
      switch(state(),
        idle    = tags$span(style = "color:#6c757d;", ""),
        running = tags$span(style = "color:#fd7e14;", ""),
        done    = tags$span(style = "color:#28a745;", "[OK]"),
        error   = tags$span(style = "color:#dc3545;", "[ERR]"),
        tags$span("")
      )
    })
    output$label <- renderUI({
      tags$small(style = "color:#6c757d;",
        switch(state(),
          idle    = "idle",
          running = tags$span(style = "color:#fd7e14; font-weight:600;", "running"),
          done    = tags$span(style = "color:#28a745;", "done"),
          error   = tags$span(style = "color:#dc3545;", "error"),
          ""
        )
      )
    })
  })
}

# UI -------------------------------------------------------------------------
ui <- page_navbar(
  title    = tags$span(
               tags$img(src = "logo.png",
                        height = "42px",
                        style  = "margin-right:8px; vertical-align:middle;"),
               tags$span("HallO-diGly",
                         style = "vertical-align:middle; font-weight:700; font-size:1.1em;")),
  theme    = bs_theme(bootswatch = "flatly", base_font = font_google("Inter")),
  fillable = FALSE,

  # Global progress + spinner styles
  header = tags$head(tags$style(HTML("
    /* ---- Shiny built-in progress bar: thicker + branded ---- */
    #shiny-progress-container .progress {
      height: 6px !important;
    }
    #shiny-progress-container .progress-bar {
      background-color: #2c7be5 !important;
    }
    #shiny-progress-container .shiny-progress-message {
      font-size: 0.92em;
      font-weight: 600;
      color: #2c7be5;
    }
    #shiny-progress-container .shiny-progress-detail {
      font-size: 0.82em;
      color: #555;
    }

    /* ---- Inline spinner shown next to action buttons ---- */
    .step-spinner {
      display: inline-flex;
      align-items: center;
      gap: 7px;
      font-size: 0.88em;
      color: #fd7e14;
      font-weight: 600;
      margin-left: 12px;
    }
    .step-spinner .spinner-border {
      width: 1rem;
      height: 1rem;
      border-width: 0.15em;
    }
    /* ---- Log box: fixed height + auto-scroll ---- */
    .log-box {
      background: #1e1e2e;
      color: #cdd6f4;
      font-family: 'Courier New', monospace;
      font-size: 0.80em;
      padding: 10px 14px;
      border-radius: 6px;
      max-height: 340px;
      overflow-y: auto;
      white-space: pre-wrap;
      word-break: break-all;
    }
    .log-box .log-ok    { color: #a6e3a1; }
    .log-box .log-warn  { color: #f9e2af; }
    .log-box .log-err   { color: #f38ba8; }
    .log-box .log-head  { color: #89b4fa; font-weight: bold; }
  "))),

  nav_panel(
    title = tagList("🔄 Restore session"),
    mod_restore_ui("restore")
  ),
  nav_panel(
    title = tagList("📂 Data ", status_ui("st_data")),
    mod_data_ui("data")
  ),
  nav_panel(
    title = tagList("⚙️ Configuration"),
    mod_config_ui("cfg")
  ),
  nav_panel(
    title = tagList("🔬 Filtering ", status_ui("st_filt")),
    mod_filtering_ui("filt")
  ),
  nav_panel(
    title = tagList("📐 Normalization ", status_ui("st_norm")),
    mod_normalization_ui("norm")
  ),
  nav_panel(
    title = tagList("🔗 Occupancy ", status_ui("st_occ")),
    mod_occupancy_ui("occ")
  ),
  nav_panel(
    title = tagList("📊 DE Analysis ", status_ui("st_de")),
    mod_de_ui("de")
  ),
  nav_panel(
    title = tagList("🌋 Plots"),
    mod_plots_ui("plots")
  ),
  nav_panel(
    title = tagList("💾 Export ", status_ui("st_export")),
    mod_export_ui("export")
  ),

  nav_panel(
    title = tagList("📖 PSP Annotation"),
    mod_psp_ui("psp")
  ),
  nav_panel(
    title = tagList("🧬 GSEA .rnk ", status_ui("st_rnk")),
    mod_rnk_ui("rnk")
  ),
  nav_panel(
    title = tagList("\U0001F9EA Proteome integration ", status_ui("st_prot")),
    mod_proteome_ui("prot")
  ),

  nav_spacer(),
  nav_item(tags$small(class = "text-muted me-2", "limpa | DIA-NN parquet"))
)

# Server ---------------------------------------------------------------------
server <- function(input, output, session) {

  st_data   <- reactiveVal("idle")
  st_filt   <- reactiveVal("idle")
  st_norm   <- reactiveVal("idle")
  st_occ    <- reactiveVal("idle")
  st_de     <- reactiveVal("idle")
  st_export <- reactiveVal("idle")
  st_rnk    <- reactiveVal("idle")
  st_prot   <- reactiveVal("idle")

  # Pipeline modules — each writes to its own reactiveVal
  data_obj <- mod_data_server("data",  st_data)
  cfg_pipe <- mod_config_server("cfg", data_obj)
  filt_obj <- mod_filtering_server("filt",  data_obj, cfg_pipe, st_filt)

  # Shared writable reactiveVals — written by pipeline OR restore
  norm_obj_rv <- reactiveVal(NULL)
  occ_obj_rv  <- reactiveVal(NULL)
  de_obj_rv   <- reactiveVal(NULL)
  cfg_rv      <- reactiveVal(NULL)

  # Forward pipeline outputs into the shared reactiveVals
  norm_pipe <- mod_normalization_server("norm", filt_obj, cfg_pipe, st_norm)
  occ_pipe  <- mod_occupancy_server("occ",  norm_pipe, cfg_pipe, st_occ)
  de_pipe   <- mod_de_server("de",  occ_pipe, cfg_pipe, st_de)

  observe({ if (!is.null(norm_pipe())) norm_obj_rv(norm_pipe()) })
  observe({ if (!is.null(occ_pipe()))  occ_obj_rv(occ_pipe())  })
  observe({ if (!is.null(de_pipe()))   de_obj_rv(de_pipe())    })
  observe({ if (!is.null(cfg_pipe()))  cfg_rv(cfg_pipe())      })

  # Restore module can overwrite the same reactiveVals
  mod_restore_server("restore", norm_obj_rv, occ_obj_rv, de_obj_rv, cfg_rv)

  # Downstream consumers use the shared reactiveVals
  mod_plots_server("plots",  de_obj_rv,  occ_obj_rv, cfg_rv)
  mod_export_server("export", de_obj_rv, occ_obj_rv, filt_obj, cfg_rv, st_export)
  mod_psp_server("psp", de_obj_rv, cfg_rv)
  mod_rnk_server("rnk", de_obj_rv, cfg_rv, st_rnk)
  mod_proteome_server("prot", de_obj_rv, cfg_rv, st_prot)

  status_server("st_data",   st_data)
  status_server("st_filt",   st_filt)
  status_server("st_norm",   st_norm)
  status_server("st_occ",    st_occ)
  status_server("st_de",     st_de)
  status_server("st_export", st_export)
  status_server("st_rnk",    st_rnk)
  status_server("st_prot",   st_prot)
}

shinyApp(ui, server)
