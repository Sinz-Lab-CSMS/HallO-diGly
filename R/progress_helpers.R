# =============================================================================
# progress_helpers.R
# Shared UI + server helpers for inline spinners and progress reporting.
#
# spinner_ui(ns, "run_btn_id")  -- place after an actionButton
# run_with_progress(session, steps, expr)  -- wraps withProgress
# =============================================================================

# Inline spinner shown next to a button while a step runs -------------------
# state_rv : reactiveVal returning "idle" | "running" | "done" | "error"
spinner_ui <- function(id) {
  ns <- NS(id)
  uiOutput(ns("spinner_out"), inline = TRUE)
}

spinner_server <- function(id, state_rv) {
  moduleServer(id, function(input, output, session) {
    output$spinner_out <- renderUI({
      s <- state_rv()
      switch(s,
        running = tags$span(
          class = "step-spinner",
          tags$span(class = "spinner-border spinner-border-sm",
                    role  = "status",
                    `aria-hidden` = "true"),
          "Running\u2026"
        ),
        done  = tags$span(
          style = "color:#28a745; font-weight:600; font-size:0.88em; margin-left:12px;",
          "\u2714 Done"),
        error = tags$span(
          style = "color:#dc3545; font-weight:600; font-size:0.88em; margin-left:12px;",
          "\u2716 Error"),
        NULL   # idle: show nothing
      )
    })
  })
}

# Named-step progress bar wrapper -------------------------------------------
# steps : named character vector, e.g.
#   c("Reading parquets" = 0.15, "DPC estimation" = 0.35, ...)
#   names = message shown; values = cumulative fraction (0-1)
# Each step emits setProgress(value, detail = name).
# Call incr_progress(label) inside the computation to advance.

make_progress_stepper <- function(session, title, steps) {
  # Returns a function: call it with the step name to advance the bar.
  # Also starts the progress display.
  progress <- shiny::Progress$new(session, min = 0, max = 1)
  progress$set(message = title, value = 0)

  step_names   <- names(steps)
  step_values  <- unname(steps)
  current_step <- 0L

  list(
    tick = function(label = NULL) {
      current_step <<- current_step + 1L
      v <- if (current_step <= length(step_values)) step_values[[current_step]] else 1
      lbl <- if (!is.null(label)) label else
             if (current_step <= length(step_names)) step_names[[current_step]] else ""
      progress$set(value = v, detail = lbl)
    },
    close = function() progress$close()
  )
}
