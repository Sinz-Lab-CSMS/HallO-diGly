# =============================================================================
# mod_restore -- Load a previously saved HallO-diGly RDS environment
#
# Repopulates norm_obj, occ_obj, and de_obj reactives so the user can jump
# straight to Plots, PSP Annotation, or Export without re-running the pipeline.
#
# The RDS must have been saved by mod_export.R (hallo_version >= "2.0").
# Older RDS files may partially restore — the module reports what was found.
# =============================================================================

mod_restore_ui <- function(id) {
  ns <- NS(id)
  tagList(

    h4("Restore a Previous Analysis"),
    card(
      card_body(
        tags$p(class = "text-muted small mb-2",
          "Upload an RDS file saved by HallO-diGly (Export tab ",
          tags$b("\u2192 Full analysis RDS"), "). ",
          "This restores quantification, occupancy correction, and DE results ",
          "so you can go directly to Plots, PSP Annotation, or Export without ",
          "re-running the full pipeline."),

        fluidRow(
          column(8,
            fileInput(ns("rds_upload"),
                      label     = "Select RDS file",
                      accept    = ".rds",
                      multiple  = FALSE,
                      width     = "100%")
          ),
          column(4,
            tags$div(style = "margin-top:24px;",
              actionButton(ns("restore"), "Restore session",
                           class = "btn-primary btn-lg w-100"),
              spinner_ui(ns("sp_restore"))
            )
          )
        ),
        verbatimTextOutput(ns("restore_log"))
      )
    ),

    conditionalPanel(
      condition = sprintf("output['%s']", ns("done")),
      card(
        card_header(tags$b("Restored session summary")),
        card_body(
          verbatimTextOutput(ns("session_summary"))
        )
      )
    )
  )
}


mod_restore_server <- function(id,
                                norm_obj_rv,
                                occ_obj_rv,
                                de_obj_rv,
                                cfg_rv) {
  moduleServer(id, function(input, output, session) {

    log_msg         <- reactiveVal("")
    sp_restore_state <- reactiveVal("idle")
    spinner_server("sp_restore", sp_restore_state)

    output$done <- reactive(nchar(log_msg()) > 0 &&
                              !grepl("\\[ERR\\]", log_msg(), fixed = TRUE) &&
                              grepl("\\[OK\\]", log_msg(), fixed = TRUE))
    outputOptions(output, "done", suspendWhenHidden = FALSE)

    alog <- function(msg) log_msg(paste0(log_msg(), msg, "\n"))

    observeEvent(input$restore, {
      f <- input$rds_upload
      if (is.null(f) || !file.exists(f$datapath)) {
        showNotification("No file uploaded.", type = "error"); return()
      }
      log_msg(""); sp_restore_state("running")

      tryCatch({
        alog(sprintf("Loading: %s (%.1f MB)", f$name, f$size / 1e6))
        env <- readRDS(f$datapath)

        if (!is.list(env))
          stop("RDS is not a list — does not appear to be a HallO-diGly session.")

        ver <- env$hallo_version %||% "(unknown)"
        alog(sprintf("HallO-diGly version tag: %s", ver))

        # ── Restore cfg ─────────────────────────────────────────────────────
        if (!is.null(env$cfg)) {
          cfg_rv(env$cfg)
          alog(sprintf("  [OK] Configuration restored (target: %s)",
                       env$cfg$target_line %||% "?"))
        } else {
          alog("  [WARN] cfg not found in RDS — Configuration tab values unchanged.")
        }

        # ── Restore norm_obj ────────────────────────────────────────────────
        required_norm <- c("q.site", "q.prot", "dpc.dg", "dpc.tp")
        missing_norm  <- setdiff(required_norm, names(env))
        if (length(missing_norm) == 0) {
          norm_obj_rv(list(
            q.site        = env$q.site,
            q.prot        = env$q.prot,
            dpc.dg        = env$dpc.dg,
            dpc.tp        = env$dpc.tp,
            y.sites       = env$y.sites,
            y.tp          = env$y.tp,
            lines_present = env$lines_present %||%
                            levels(droplevels(env$y.sites$targets$cell_line))
          ))
          alog(sprintf("  [OK] Normalisation restored (%d sites, %d proteins)",
                       nrow(env$q.site), nrow(env$q.prot)))
        } else {
          alog(sprintf("  [WARN] Normalisation incomplete — missing: %s",
                       paste(missing_norm, collapse = ", ")))
        }

        # ── Restore occ_obj ─────────────────────────────────────────────────
        required_occ <- c("q.occ", "site_info")
        missing_occ  <- setdiff(required_occ, names(env))
        if (length(missing_occ) == 0) {
          occ_obj_rv(list(
            q.site        = env$q.site,
            q.prot        = env$q.prot,
            q.occ         = env$q.occ,
            dpc.dg        = env$dpc.dg,
            dpc.tp        = env$dpc.tp,
            y.sites       = env$y.sites,
            y.tp          = env$y.tp,
            site_info     = env$site_info,
            link          = env$link,
            lines_present = env$lines_present %||%
                            levels(droplevels(env$y.sites$targets$cell_line))
          ))
          n_corrected <- if (!is.null(env$site_info))
            sum(env$site_info$tier == "corrected", na.rm = TRUE) else NA
          alog(sprintf("  [OK] Occupancy restored (%d corrected sites)",
                       n_corrected))
        } else {
          alog(sprintf("  [WARN] Occupancy incomplete — missing: %s",
                       paste(missing_occ, collapse = ", ")))
        }

        # ── Restore de_obj ──────────────────────────────────────────────────
        required_de <- c("fit.occ", "res")
        missing_de  <- setdiff(required_de, names(env))
        if (length(missing_de) == 0) {
          de_obj_rv(list(
            fit.site = env$fit.site,
            fit.occ  = env$fit.occ,
            fit.prot = env$fit.prot,
            tt.site  = env$tt.site,
            tt.occ   = env$tt.occ,
            tt.prot  = env$tt.prot,
            design   = env$design,
            cm       = env$contrasts,
            res      = env$res
          ))
          n_sig <- if (!is.null(env$res) && "FDR_occ" %in% names(env$res))
            sum(env$res$FDR_occ < (env$cfg$fdr_threshold %||% 0.05),
                na.rm = TRUE) else NA
          alog(sprintf("  [OK] DE results restored (%s significant at occupancy layer)",
                       if (is.na(n_sig)) "?" else format(n_sig, big.mark = ",")))
        } else {
          alog(sprintf("  [WARN] DE incomplete — missing: %s",
                       paste(missing_de, collapse = ", ")))
        }

        alog("\n[OK] Session restored. Proceed to Plots, PSP, or Export.")
        sp_restore_state("done")

      }, error = function(e) {
        alog(paste0("[ERR] ", e$message))
        sp_restore_state("error")
      })
    })

    output$restore_log <- renderPrint({ cat(log_msg()) })

    output$session_summary <- renderPrint({
      req(nchar(log_msg()) > 0)
      occ  <- occ_obj_rv()
      de   <- de_obj_rv()
      cfg  <- cfg_rv()
      if (is.null(occ) || is.null(de)) {
        cat("Restoration incomplete.\n"); return()
      }
      cat(sprintf("Target line   : %s\n", cfg$target_line %||% "?"))
      cat(sprintf("Ref lines     : %s\n",
                  paste(cfg$reference_lines %||% "?", collapse = ", ")))
      cat(sprintf("Cell lines    : %s\n",
                  paste(occ$lines_present %||% "?", collapse = ", ")))
      cat(sprintf("Sites (occ)   : %d\n", nrow(occ$q.occ)))
      cat(sprintf("Sites (total) : %d\n", nrow(occ$q.site)))
      cat(sprintf("Proteins      : %d\n", nrow(occ$q.prot)))
      if (!is.null(de$res)) {
        fdr_t <- cfg$fdr_threshold %||% 0.05
        if ("FDR_occ" %in% names(de$res))
          cat(sprintf("Sig occ sites : %d (FDR < %g)\n",
                      sum(de$res$FDR_occ < fdr_t, na.rm = TRUE), fdr_t))
      }
    })

    `%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b
  })
}
