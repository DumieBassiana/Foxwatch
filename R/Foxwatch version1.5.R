# =============================================================================
#
#    ███████╗ ██████╗ ██╗  ██╗██╗    ██╗ █████╗ ████████╗ ██████╗██╗  ██╗
#    ██╔════╝██╔═══██╗╚██╗██╔╝██║    ██║██╔══██╗╚══██╔══╝██╔════╝██║  ██║
#    █████╗  ██║   ██║ ╚███╔╝ ██║ █╗ ██║███████║   ██║   ██║     ███████║
#    ██╔══╝  ██║   ██║ ██╔██╗ ██║███╗██║██╔══██║   ██║   ██║     ██╔══██║
#    ██║     ╚██████╔╝██╔╝ ██╗╚███╔███╔╝██║  ██║   ██║   ╚██████╗██║  ██║
#    ╚═╝      ╚═════╝ ╚═╝  ╚═╝ ╚══╝╚══╝ ╚═╝  ╚═╝   ╚═╝    ╚═════╝╚═╝  ╚═╝
#
#                  🦊  Smart Data Management & Abundance Estimation  🦊
#
# =============================================================================
#  Title   : FoxWatch — A User-Friendly Tool for Estimating Fox Abundance
#            from Camera Traps Using Simultaneous-Count Models
#  Version : 1.5
#  Authors : Duminda S.B. Dissanayake & Graeme Armstrong
#  Year    : 2026
#  Contact : dumie.dissanayake@dcceew.nsw.gov.au
# -----------------------------------------------------------------------------
#  Programme    : NSW Saving Our Species
#  Paper        : <DOI link when available>
#  Repository   : https://github.com/DumieBassiana/Foxwatch
# -----------------------------------------------------------------------------
#  Description:
#    FoxWatch is an R Shiny application that enables field ecologists and
#    conservation practitioners to estimate red fox (Vulpes vulpes) abundance
#    from camera-trap data using a Bayesian Simultaneous-Count Model (SCM),
#    without requiring statistical programming expertise.
#
#    Key features:
#      • Automated EXIF metadata parsing from camera-trap text files
#      • Interactive spatial visualisation of fox detections per site per year
#      • Bayesian effort-scaled Simultaneous-Count Model (SCM) fitted via
#        JAGS. Detection probability per replicate scales with survey effort
#        via p_eff = 1 - (1 - p)^effort, and all N[i] share a common Gamma
#        prior on lambda, so blocks with sparse data borrow strength from the
#        others (partial pooling). Default priors Beta(1, 19) on p and
#        Gamma(2, 0.25) on lambda encode a rare-detection, low-abundance
#        regime; all priors and MCMC settings are user-adjustable in the UI.
#      • Gelman-Rubin R-hat, effective sample size (ESS) and Monte-Carlo
#        standard error (MCSE) reported for every parameter
#      • Colour-coded results tables and uncertainty ribbon plots
#      • One-click CSV export of all abundance estimates
# -----------------------------------------------------------------------------
#  Dependencies:
#    shiny, dplyr, readr, stringr, DT, tools, ggplot2, plotly,
#    shinycssloaders, lubridate, rjags, modeest, coda,
#    sf, ggspatial, viridis, ggrepel, hms, colorspace,
#    RColorBrewer, forcats, grid
#    External: JAGS >= 4.3.1  (https://mcmc-jags.sourceforge.io)
# -----------------------------------------------------------------------------
#  Citation:
#    Dissanayake D.S.B. & Armstrong G. (2026) FoxWatch: A User-Friendly Tool
#    for Estimating Fox Abundance from Camera Traps Using Simultaneous-Count
#    Models. <Journal>. <DOI>
#
#    Armstrong G. & McSorley A. (2024) Estimating fox abundance using the
#    simultaneous count method. Wildlife Letters, 2, 102-109.
#    https://doi.org/10.1002/wll2.12038
# -----------------------------------------------------------------------------
#  Licence : MIT — see LICENSE file in repository root
#  Notes   : Requires R >= 4.2.0 and JAGS >= 4.3.1
# =============================================================================

library(shiny)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(DT)
library(tools)
library(ggplot2)
library(plotly)
library(shinycssloaders)
library(lubridate)
library(rjags)
library(modeest)
library(coda)
library(sf)
library(ggspatial)
library(viridis)
library(ggrepel)
library(hms)
library(colorspace)
library(RColorBrewer)
library(forcats)
library(grid)

# ---------------------------------------------------------------------------
# Parse EXIF lines — fox detections only
# ---------------------------------------------------------------------------
parse_exif_lines <- function(lines) {
  data <- lapply(lines, function(line) {
    match <- str_match(
      line,
      "^(IMG_\\d+)-(\\d{1,2}/\\d{1,2}/\\d{4})\\s+(\\d{1,2}:\\d{2}:\\d{2})\\s+[APM]{2}.*-(Fox)"
    )
    if (!is.na(match[1, 1])) {
      return(data.frame(
        Date      = match[1, 3],
        Time      = match[1, 4],
        Detection = match[1, 5],
        Image_ID  = match[1, 2],
        stringsAsFactors = FALSE
      ))
    }
    return(NULL)
  })
  bind_rows(data)
}

# ---------------------------------------------------------------------------
# Vector layer reader
# ---------------------------------------------------------------------------
read_vector_layer <- function(upload) {
  req(upload)
  vdir <- file.path(tempdir(), paste0("v_", as.integer(runif(1, 1e9))))
  dir.create(vdir, showWarnings = FALSE, recursive = TRUE)

  zip_idx <- grepl("\\.zip$", upload$name, ignore.case = TRUE)
  if (any(zip_idx)) {
    zpath    <- upload$datapath[which(zip_idx)[1]]
    utils::unzip(zpath, exdir = vdir)
    shp_path <- list.files(vdir, pattern = "\\.shp$", full.names = TRUE,
                           ignore.case = TRUE)
    if (length(shp_path) > 0) {
      Sys.setenv(GDAL_SHAPE_RESTORE_SHX = "YES")
      return(sf::st_read(shp_path[1], quiet = TRUE))
    }
  }
  out_paths <- file.path(vdir, upload$name)
  file.copy(upload$datapath, out_paths, overwrite = TRUE)
  shp_path  <- out_paths[grepl("\\.shp$", upload$name, ignore.case = TRUE)]
  if (length(shp_path) > 0) {
    Sys.setenv(GDAL_SHAPE_RESTORE_SHX = "YES")
    return(sf::st_read(shp_path[1], quiet = TRUE))
  }
  return(NULL)
}

# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------
ui <- fluidPage(

  tags$head(
    # Google Fonts
    tags$link(
      rel  = "stylesheet",
      href = "https://fonts.googleapis.com/css2?family=Inter:wght@400;600;700&family=Montserrat:wght@700;800&display=swap"
    ),

    tags$style(HTML("

      /* ── Base ─────────────────────────────────────────────────── */
      body {
        background: linear-gradient(135deg, #f0f4f8 0%, #e8edf2 100%);
        font-family: 'Inter', 'Segoe UI', sans-serif;
        font-size: 16px;
        line-height: 1.65;
        color: #1a2332;
        min-height: 100vh;
      }

      /* ── Banner ───────────────────────────────────────────────── */
      .fw-banner {
        background: linear-gradient(135deg, #0d1b2a 0%, #1b3a5c 45%, #0d1b2a 100%);
        padding: 32px 40px 28px 40px;
        border-radius: 16px;
        text-align: center;
        margin-bottom: 28px;
        position: relative;
        overflow: hidden;
        box-shadow: 0 8px 32px rgba(13,27,42,0.45);
      }
      .fw-banner::before {
        content: '';
        position: absolute;
        top: -60px; left: -60px;
        width: 220px; height: 220px;
        background: radial-gradient(circle, rgba(230,126,34,0.18) 0%, transparent 70%);
        border-radius: 50%;
      }
      .fw-banner::after {
        content: '';
        position: absolute;
        bottom: -40px; right: -40px;
        width: 180px; height: 180px;
        background: radial-gradient(circle, rgba(52,152,219,0.15) 0%, transparent 70%);
        border-radius: 50%;
      }
      .fw-title {
        font-family: 'Montserrat', 'Inter', sans-serif;
        font-size: 36px;
        font-weight: 800;
        color: #ffffff;
        letter-spacing: -0.5px;
        margin: 0 0 6px 0;
        line-height: 1.2;
        position: relative; z-index: 1;
      }
      .fw-title .fox-accent {
        color: #e67e22;
      }
      .fw-subtitle {
        font-family: 'Inter', sans-serif;
        font-size: 15px;
        font-weight: 400;
        color: #a8c0d6;
        margin: 0 0 16px 0;
        letter-spacing: 0.3px;
        position: relative; z-index: 1;
      }
      .fw-badge-row {
        display: flex;
        justify-content: center;
        align-items: center;
        gap: 10px;
        flex-wrap: wrap;
        position: relative; z-index: 1;
      }
      .fw-badge {
        display: inline-flex;
        align-items: center;
        gap: 6px;
        padding: 5px 14px;
        border-radius: 20px;
        font-size: 12px;
        font-weight: 600;
        letter-spacing: 0.4px;
      }
      .fw-badge-version {
        background: rgba(230,126,34,0.20);
        border: 1px solid rgba(230,126,34,0.50);
        color: #f0a257;
      }
      .fw-badge-method {
        background: rgba(52,152,219,0.18);
        border: 1px solid rgba(52,152,219,0.40);
        color: #7ec8f5;
      }
      .fw-badge-org {
        background: rgba(255,255,255,0.08);
        border: 1px solid rgba(255,255,255,0.18);
        color: #c8d8e8;
      }

      /* ── Sidebar ──────────────────────────────────────────────── */
      .well {
        background: #ffffff;
        border: none;
        border-radius: 14px;
        box-shadow: 0 4px 18px rgba(0,0,0,0.08);
        padding: 20px 18px;
      }
      .sidebar-title {
        font-family: 'Montserrat', sans-serif;
        font-size: 15px;
        font-weight: 700;
        color: #0d1b2a;
        text-transform: uppercase;
        letter-spacing: 1px;
        margin-bottom: 16px;
        padding-bottom: 10px;
        border-bottom: 2px solid #e67e22;
      }

      /* ── Upload cards ─────────────────────────────────────────── */
      .upload-card {
        background: #f8fafd;
        border-radius: 10px;
        padding: 14px 16px;
        margin-bottom: 14px;
        border-left: 4px solid #1b3a5c;
        transition: border-color 0.2s;
      }
      .upload-card:hover { border-left-color: #e67e22; }
      .upload-card-title {
        font-size: 13px;
        font-weight: 700;
        color: #1b3a5c;
        text-transform: uppercase;
        letter-spacing: 0.6px;
        margin-bottom: 10px;
      }
      .upload-card .form-group { margin-bottom: 8px; }
      .upload-card label {
        font-size: 13px;
        color: #4a5568;
        font-weight: 500;
      }
      .help-block { font-size: 12px; color: #718096; }

      /* ── Tabs ─────────────────────────────────────────────────── */
      .nav-tabs {
        border-bottom: 2px solid #e2e8f0;
        margin-bottom: 0;
      }
      .nav-tabs > li > a {
        font-family: 'Inter', sans-serif;
        font-size: 14px;
        font-weight: 600;
        color: #4a5568;
        background: #f1f5f9;
        border: 1px solid #e2e8f0;
        border-bottom: none;
        border-radius: 10px 10px 0 0;
        margin-right: 4px;
        padding: 10px 20px;
        transition: all 0.2s;
      }
      .nav-tabs > li > a:hover {
        background: #e8edf5;
        color: #1b3a5c;
      }
      .nav-tabs > li.active > a,
      .nav-tabs > li.active > a:hover {
        background: linear-gradient(135deg, #1b3a5c, #0d5c8a);
        color: #ffffff;
        font-weight: 700;
        border-color: #1b3a5c;
        box-shadow: 0 -2px 0 #e67e22 inset;
      }
      .tab-content {
        background: #ffffff;
        border: 1px solid #e2e8f0;
        border-top: none;
        border-radius: 0 0 14px 14px;
        padding: 24px;
        box-shadow: 0 4px 18px rgba(0,0,0,0.06);
      }

      /* ── Section cards ────────────────────────────────────────── */
      .section-card {
        background: #ffffff;
        border-radius: 12px;
        padding: 22px 24px;
        margin: 18px 0;
        border: 1px solid #e8edf5;
        box-shadow: 0 2px 10px rgba(0,0,0,0.05);
        transition: box-shadow 0.2s;
      }
      .section-card:hover { box-shadow: 0 4px 20px rgba(0,0,0,0.09); }
      .section-card h3 {
        font-family: 'Montserrat', sans-serif;
        font-size: 17px;
        font-weight: 700;
        color: #0d1b2a;
        margin-top: 0;
        margin-bottom: 16px;
        padding-bottom: 10px;
        border-bottom: 2px solid #f0f4f8;
      }

      /* ── Headings ─────────────────────────────────────────────── */
      h3 {
        font-family: 'Montserrat', sans-serif;
        font-size: 17px;
        font-weight: 700;
        color: #0d1b2a;
      }
      h4 {
        font-family: 'Inter', sans-serif;
        font-size: 15px;
        font-weight: 700;
        color: #1b3a5c;
      }

      /* ── Buttons ──────────────────────────────────────────────── */
      .btn {
        font-family: 'Inter', sans-serif;
        font-size: 14px;
        font-weight: 600;
        border-radius: 8px;
        padding: 9px 20px;
        transition: all 0.2s;
        border: none;
        letter-spacing: 0.2px;
      }
      .btn-success {
        background: linear-gradient(135deg, #27ae60, #1e8449);
        color: #ffffff;
        box-shadow: 0 3px 10px rgba(39,174,96,0.30);
      }
      .btn-success:hover {
        background: linear-gradient(135deg, #1e8449, #196f3d);
        box-shadow: 0 4px 14px rgba(39,174,96,0.40);
        transform: translateY(-1px);
        color: #ffffff;
      }
      .btn-primary {
        background: linear-gradient(135deg, #1b3a5c, #0d5c8a);
        color: #ffffff;
        box-shadow: 0 3px 10px rgba(27,58,92,0.30);
      }
      .btn-primary:hover {
        background: linear-gradient(135deg, #0d5c8a, #0a4f7a);
        box-shadow: 0 4px 14px rgba(27,58,92,0.40);
        transform: translateY(-1px);
        color: #ffffff;
      }
      .btn-lg { padding: 12px 28px; font-size: 15px; }

      /* ── Download links ───────────────────────────────────────── */
      .shiny-download-link {
        display: inline-block;
        margin-top: 8px;
        background: linear-gradient(135deg, #2471a3, #1a5276);
        color: #ffffff !important;
        padding: 9px 20px;
        border-radius: 8px;
        font-size: 14px;
        font-weight: 600;
        text-decoration: none;
        transition: all 0.2s;
        box-shadow: 0 3px 10px rgba(36,113,163,0.30);
      }
      .shiny-download-link:hover {
        background: linear-gradient(135deg, #1a5276, #154360);
        transform: translateY(-1px);
        box-shadow: 0 4px 14px rgba(36,113,163,0.40);
      }

      /* ── Convergence legend ───────────────────────────────────── */
      .convergence-legend {
        font-size: 13px;
        margin-top: 12px;
        padding: 12px 16px;
        background: #f8fafd;
        border-radius: 8px;
        border-left: 4px solid #1b3a5c;
      }
      .convergence-legend ul { margin: 6px 0 0 0; padding-left: 20px; }
      .convergence-legend li { margin-bottom: 4px; }

      /* ── DT tables ────────────────────────────────────────────── */
      .dataTables_wrapper .dataTables_filter input,
      .dataTables_wrapper .dataTables_length select {
        font-size: 13px;
        padding: 5px 10px;
        border: 1px solid #cbd5e0;
        border-radius: 6px;
      }
      .dataTables_wrapper .dataTables_info,
      .dataTables_wrapper .dataTables_paginate { font-size: 13px; }
      table.dataTable thead th {
        background: #1b3a5c;
        color: #ffffff;
        font-family: 'Inter', sans-serif;
        font-size: 13px;
        font-weight: 600;
      }

      /* ── Spinner ──────────────────────────────────────────────── */
      .shiny-spinner-output-container .shiny-spinner-placeholder {
        color: #e67e22 !important;
      }

      /* ── Footer ───────────────────────────────────────────────── */
      .fw-footer {
        text-align: center;
        padding: 24px 0 16px 0;
        color: #718096;
        font-size: 12px;
        border-top: 1px solid #e2e8f0;
        margin-top: 20px;
      }

      /* ── HR ───────────────────────────────────────────────────── */
      hr { border-color: #e2e8f0; }

    "))
  ),

  # ── BANNER ──────────────────────────────────────────────────────────────
  div(class = "fw-banner",
      div(class = "fw-title",
          HTML("&#x1F98A; Fox<span class='fox-accent'>Watch</span>")
      ),
      div(class = "fw-subtitle",
          "Estimating Fox Abundance from Camera Traps Using Simultaneous-Count Models"
      ),
      div(class = "fw-badge-row",
          span(class = "fw-badge fw-badge-version",
               HTML("&#x25C6; Version 1.1"))
      )
  ),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      div(class = "sidebar-title", "Data Upload"),

      # Spatial Data
      div(class = "upload-card",
          div(class = "upload-card-title",
              icon("map"), " Spatial Data"),
          fileInput("fox_csv",   "Fox Detection Records (CSV)",
                    accept = ".csv"),
          fileInput("sites_csv", "Site Locations (CSV)",
                    accept = ".csv"),
          fileInput("park_shp",  "Park Boundary (shapefile / zip)",
                    multiple = TRUE,
                    accept   = c(".shp", ".shx", ".dbf", ".prj", ".zip"))
      ),

      # Camera Data
      div(class = "upload-card",
          div(class = "upload-card-title",
              icon("camera"), " Camera Data"),
          fileInput("txtfiles", "Camera TXT Files (EXIF)",
                    multiple = TRUE,
                    accept   = ".txt"),
          helpText("Multiple files can be selected")
      ),

      # Abundance Data
      div(class = "upload-card",
          div(class = "upload-card-title",
              icon("chart-line"), " Abundance Data"),
          fileInput("abundance_csv",
                    "Daily Counts (TXT) or Fox Detections (CSV)",
                    accept = c(".txt", ".csv")),
          helpText("TXT: one integer count per line. CSV: SiteID, Date, Time, Detection")
      )
    ),

    mainPanel(
      width = 9,

      tabsetPanel(
        id   = "main_tabs",
        type = "tabs",

        # ── TAB 1: Spatial Analysis ────────────────────────────────────
        tabPanel(
          title = tagList(icon("map-marked-alt"), " Spatial Analysis"),
          br(),
          div(class = "section-card",
              h3(icon("map"), " Fox Detections Per Site Per Year"),
              fluidRow(
                column(4,
                       selectInput("min_year", "Minimum year:",
                                   choices  = 2019:2030,
                                   selected = 2023)),
                column(4,
                       textInput("drop_year", "Exclude year:",
                                 value = "2022")),
                column(4,
                       br(),
                       actionButton("make_map", " Generate Map",
                                    class = "btn btn-success",
                                    icon  = icon("map")))
              ),
              hr(),
              withSpinner(plotOutput("detections_map", height = "650px"),
                          color = "#e67e22"),
              br(),
              downloadButton("download_map_png", " Download Map (PNG)",
                             class = "btn btn-primary")
          )
        ),

        # ── TAB 2: Camera Data & Detection ────────────────────────────
        tabPanel(
          title = tagList(icon("camera-retro"), " Camera Data & Detection"),
          br(),
          div(class = "section-card",
              h3(icon("cogs"), " Process Camera Files"),
              uiOutput("camera_tabs")
          ),
          br(),
          div(class = "section-card",
              h3(icon("table"), " Summary Statistics"),
              withSpinner(DTOutput("summary_table"), color = "#e67e22")
          ),
          br(),
          div(class = "section-card",
              h3(icon("chart-bar"), " Detection Visualisation"),
              withSpinner(plotlyOutput("detection_plot", height = "480px"),
                          color = "#e67e22"),
              br(),
              div(style = "text-align:right;",
                  tags$small(style = "color:#718096; margin-right:10px;",
                             "Publication-quality export:"),
                  downloadButton("download_detection_png",
                                 label = tagList(icon("image"), " PNG (300 dpi)"),
                                 class = "btn btn-primary"),
                  span(" "),
                  downloadButton("download_detection_pdf",
                                 label = tagList(icon("file-pdf"), " PDF"),
                                 class = "btn btn-primary")
              )
          ),
          br(),
          div(class = "section-card",
              h3(icon("download"), " Export Detection Data"),
              div(style = "text-align: center;",
                  downloadButton("download_fox_only",
                                 label = tagList(icon("paw"),
                                                 " Download Fox Detections"))
              )
          ),
          br(),
          div(class = "section-card",
              style = "border-left: 4px solid #e67e22;",
              h3(icon("arrow-right"), " Generate Abundance Model Input"),
              p(style = "font-size:14px; color:#4a5568; margin-bottom:14px;",
                "Click below to aggregate processed fox detections into daily
                camera-presence counts ready for the Abundance Model tab.
                For each day, FoxWatch counts the number of camera stations
                that detected at least one fox, fills days with no detections
                as zero, and outputs one integer per line — the exact
                format required by the SCM."),
              div(style = "text-align: center;",
                  actionButton("generate_daily_counts",
                               label = tagList(icon("calculator"),
                                               " Generate Daily Counts"),
                               class = "btn btn-success"),
                  span("  "),
                  downloadButton("download_daily_counts",
                                 label = tagList(icon("file-alt"),
                                                 " Download for Model"))
              ),
              br(),
              uiOutput("daily_counts_preview")
          )
        ),

        # -- TAB 3: Fox Abundance Model ----------------------------------
        tabPanel(
          title = tagList(icon("chart-line"),
                          " Fox Abundance Model"),
          br(),
          div(class = "section-card",
              style = "border-left: 4px solid #0d5c8a;",
              h3(icon("sliders-h"),
                 " Fox Abundance Model - Effort-Scaled SCM with Partial Pooling"),
              br(),
              p(style = "font-size:14px; color:#4a5568;",
                "FoxWatch fits a Bayesian Simultaneous-Count Model (SCM; Armstrong &
                McSorley 2024) to your daily fox detection counts. Per-visit
                detection probability scales with survey effort via the
                geometric-complement link ",
                tags$em("p_eff = 1 - (1 - p)^effort"),
                " -- i.e. if the daily detection probability is p, then
                over k trap-nights the probability of at least one detection
                is 1 - (1 - p)^k. All occasions are fitted jointly and
                share a common Poisson mean lambda, so blocks with sparse
                data borrow strength from the others (partial pooling).
                Default priors Beta(1, 19) on p (mean ~ 0.05) and
                Gamma(2, 0.25) on lambda (mean 8) encode a rare-detection,
                low-abundance regime. All priors and MCMC settings below
                are user-adjustable."),
              hr(),
              h4(icon("bezier-curve"), " Prior hyperparameters"),
              fluidRow(
                column(3, numericInput("mb_pa",
                                       "Beta prior alpha on p",
                                       value = 1, min = 0.01, step = 0.1)),
                column(3, numericInput("mb_pb",
                                       "Beta prior beta on p",
                                       value = 19, min = 0.01, step = 0.5)),
                column(3, numericInput("mb_shape",
                                       "Gamma shape on lambda",
                                       value = 2, min = 0.01, step = 0.1)),
                column(3, numericInput("mb_rate",
                                       "Gamma rate on lambda",
                                       value = 0.25, min = 0.001,
                                       step = 0.05))
              ),
              hr(),
              h4(icon("random"), " MCMC settings"),
              fluidRow(
                column(2, numericInput("mb_chains", "Chains",
                                       value = 3, min = 1, max = 8,
                                       step = 1)),
                column(2, numericInput("mb_adapt", "Adapt",
                                       value = 1000, min = 100,
                                       step = 100)),
                column(3, numericInput("mb_burnin", "Burn-in",
                                       value = 5000, min = 500,
                                       step = 500)),
                column(3, numericInput("mb_iter", "Samples",
                                       value = 10000, min = 1000,
                                       step = 1000)),
                column(2, numericInput("mb_thin", "Thin",
                                       value = 1, min = 1, step = 1))
              ),
              hr(),
              h4(icon("video"), " Camera array size"),
              fluidRow(
                column(4, numericInput("mb_n_cameras",
                                       "Number of active camera stations",
                                       value = 15, min = 1, step = 1)),
                column(8,
                       helpText("Total number of cameras deployed across the array. Used for a data-quality warning (if any daily count exceeds this number, the app flags it as a probable data-entry issue) and for a derived foxes-per-camera metric shown alongside the abundance table. Does not enter the JAGS likelihood."))
              ),
              hr(),
              h4(icon("clock"), " Effort matrix (optional)"),
              fileInput("mb_effort_csv",
                        "Effort CSV (rows = occasions, cols = replicates)",
                        accept = ".csv"),
              helpText("If omitted, effort = 1 for every cell (equivalent to
                        one trap-night per replicate). Provide a CSV with no
                        header, one row per 10-day block, one column per day
                        (10 columns), where each value is the number of
                        active trap-nights that day."),
              hr(),
              div(style = "text-align: center;",
                  actionButton("run_model_b",
                               label = tagList(icon("play-circle"),
                                               " Run Abundance Model"),
                               class = "btn btn-success btn-lg"))
          ),
          br(),
          div(class = "section-card",
              style = "border-left: 4px solid #27ae60;",
              h3(icon("balance-scale"),
                 " Study-wide Summary"),
              withSpinner(uiOutput("mb_summary_text"),
                          color = "#27ae60")
          ),
          br(),
          div(class = "section-card",
              h3(icon("table"),
                 " Per-Occasion Abundance Estimates (N[i])"),
              withSpinner(DTOutput("mb_N_table"), color = "#0d5c8a")
          ),
          br(),
          div(class = "section-card",
              h3(icon("cogs"), " Global Parameters (lambda, p)"),
              withSpinner(DTOutput("mb_param_table"), color = "#0d5c8a")
          ),
          br(),
          div(class = "section-card",
              h3(icon("stethoscope"),
                 " Convergence Diagnostics (R-hat, ESS, MCSE)"),
              p(style = "font-size:13px; color:#4a5568;",
                "R-hat colour key: green <= 1.05 (good),
                amber 1.05-1.10 (acceptable), red > 1.10 (poor
                -- treat that parameter with caution)."),
              withSpinner(DTOutput("mb_diag_table"), color = "#0d5c8a")
          ),
          br(),
          div(class = "section-card",
              h3(icon("chart-area"),
                 " Posterior Mode Abundance with 95% HDI"),
              withSpinner(plotlyOutput("mb_mode_plot",
                                       height = "500px"),
                          color = "#0d5c8a"),
              br(),
              div(style = "text-align:right;",
                  tags$small(style = "color:#718096; margin-right:8px;",
                             "Publication-quality export:"),
                  downloadButton("mb_download_mode_png",
                                 label = tagList(icon("image"),
                                                 " PNG (300 dpi)"),
                                 class = "btn btn-primary"),
                  span(" "),
                  downloadButton("mb_download_mode_pdf",
                                 label = tagList(icon("file-pdf"),
                                                 " PDF"),
                                 class = "btn btn-primary"))
          ),
          br(),
          div(class = "section-card",
              h3(icon("chart-area"),
                 " Posterior Mean Abundance with 95% HDI"),
              withSpinner(plotlyOutput("mb_mean_plot",
                                       height = "500px"),
                          color = "#0d5c8a"),
              br(),
              div(style = "text-align:right;",
                  tags$small(style = "color:#718096; margin-right:8px;",
                             "Publication-quality export:"),
                  downloadButton("mb_download_mean_png",
                                 label = tagList(icon("image"),
                                                 " PNG (300 dpi)"),
                                 class = "btn btn-primary"),
                  span(" "),
                  downloadButton("mb_download_mean_pdf",
                                 label = tagList(icon("file-pdf"),
                                                 " PDF"),
                                 class = "btn btn-primary"))
          ),
          br(),
          div(class = "section-card",
              h3(icon("download"), " Export Results"),
              div(style = "text-align: center;",
                  downloadButton("mb_download_csv",
                                 label = tagList(icon("file-csv"),
                                                 " Download Estimates (CSV)"),
                                 class = "btn btn-primary"),
                  span("  "),
                  downloadButton("mb_download_params_csv",
                                 label = tagList(icon("file-csv"),
                                                 " Download lambda, p Table (CSV)"),
                                 class = "btn btn-primary"))
          )
        )
      ),

      # ── FOOTER ──────────────────────────────────────────────────────
      div(class = "fw-footer",
          actionButton(
            inputId = "contact_button",
            label   = tagList(icon("envelope"), " Contact for Support"),
            class   = "btn btn-primary",
            onclick = "window.location.href='mailto:dumie.dissanayake@dcceew.nsw.gov.au?subject=FoxWatch%20v1.1%20Support';"
          ),
          br(), br(),
          HTML("FoxWatch v1.1 &nbsp;&bull;&nbsp;
                Dissanayake &amp; Armstrong 2026 &nbsp;&bull;&nbsp;
                <a href='https://github.com/DumieBassiana/Foxwatch'
                   target='_blank'
                   style='color:#718096;'>GitHub</a>")
      )
    )
  )
)

# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------
server <- function(input, output, session) {

  camera_files      <- reactiveValues(files = list(), cleaned = list())
  master_detections <- reactiveVal(data.frame())
  summary_table_rv  <- reactiveVal(data.frame())
  map_plot_obj      <- reactiveVal(NULL)

  # ── Camera file upload ────────────────────────────────────────────────
  observeEvent(input$txtfiles, {
    req(input$txtfiles)
    for (i in 1:nrow(input$txtfiles)) {
      name      <- input$txtfiles$name[i]
      path      <- input$txtfiles$datapath[i]
      camera_id <- str_extract(name, "Ool_\\d+")
      if (is.na(camera_id) || is.null(camera_id))
        camera_id <- tools::file_path_sans_ext(name)
      if (!camera_id %in% names(camera_files$files))
        camera_files$files[[camera_id]] <- list(name = name, path = path)
    }
    showNotification(paste("Loaded", nrow(input$txtfiles), "camera file(s)"),
                     type = "message", duration = 3)
  })

  # ── Camera tabs UI ────────────────────────────────────────────────────
  output$camera_tabs <- renderUI({
    req(camera_files$files)
    tabs <- lapply(names(camera_files$files), function(camera_id) {
      tabPanel(
        title = camera_id,
        br(),
        fluidRow(column(12,
                        h4(paste("Camera File:", camera_files$files[[camera_id]]$name)),
                        DTOutput(paste0("table_", camera_id)),
                        br(),
                        actionButton(paste0("process_", camera_id),
                                     paste("Process", camera_id),
                                     class = "btn btn-success",
                                     icon  = icon("cogs")),
                        span("  "),
                        downloadButton(paste0("download_", camera_id),
                                       paste("Download", camera_id, ".txt"))
        ))
      )
    })
    do.call(tabsetPanel, c(tabs, id = "camera_tabset"))
  })

  # ── Camera processing ─────────────────────────────────────────────────
  observe({
    lapply(names(camera_files$files), function(camera_id) {

      output[[paste0("table_", camera_id)]] <- renderDT({
        req(camera_files$files[[camera_id]])
        lines <- read_lines(camera_files$files[[camera_id]]$path)
        df    <- parse_exif_lines(lines)
        datatable(df, options = list(pageLength = 10, scrollX = TRUE))
      })

      observeEvent(input[[paste0("process_", camera_id)]], {
        lines <- read_lines(camera_files$files[[camera_id]]$path)
        df    <- parse_exif_lines(lines)
        df$SiteID <- camera_id
        df <- df %>% select(SiteID, Date, Time, Detection, Image_ID)
        colnames(df)[colnames(df) == "Image_ID"] <- "Image.ID"

        cleaned_file <- file.path(tempdir(), paste0(camera_id, ".txt"))
        write_tsv(df, cleaned_file, quote = "all")
        camera_files$cleaned[[camera_id]] <- list(data = df,
                                                  file = cleaned_file)
        current <- master_detections()
        master_detections(bind_rows(current, df))

        showNotification(paste("Processed", camera_id, "successfully"),
                         type = "message", duration = 3)

        latest_summary <- master_detections() %>%
          filter(Detection == "Fox") %>%
          group_by(SiteID, Detection) %>%
          summarise(Total = n(), .groups = "drop")
        summary_table_rv(latest_summary)
      })

      output[[paste0("download_", camera_id)]] <- downloadHandler(
        filename = function() paste0(camera_id, ".txt"),
        content  = function(file) {
          req(camera_files$cleaned[[camera_id]]$file)
          file.copy(camera_files$cleaned[[camera_id]]$file, file)
        }
      )
    })
  })

  # ── Summary table ─────────────────────────────────────────────────────
  output$summary_table <- renderDT({
    df <- summary_table_rv()
    if (!is.null(df) && nrow(df) > 0) {
      datatable(df, options = list(pageLength = 15))
    } else {
      datatable(data.frame(
        Message = "Process camera files to see summary"))
    }
  })

  # ── Detection plot ─────────────────────────────────────────────────────
  output$detection_plot <- renderPlotly({
    df <- master_detections()
    if (nrow(df) == 0) return(plotly_empty())
    df      <- df %>% filter(Detection == "Fox")
    df$SiteID <- factor(df$SiteID)
    p <- ggplot(df, aes(x = SiteID)) +
      geom_bar(fill = "#e67e22", width = 0.65) +
      theme_minimal(base_family = "sans") +
      labs(title = "Total Fox Detections per Camera",
           x = "Camera Station", y = "Detection Count") +
      theme(
        plot.title   = element_text(face = "bold", size = 14,
                                    color = "#0d1b2a"),
        axis.text.x  = element_text(angle = 45, hjust = 1, size = 11),
        axis.title   = element_text(face = "bold", size = 12),
        panel.grid.major.x = element_blank()
      )
    ggplotly(p)
  })

  # ── Fox download ──────────────────────────────────────────────────────
  output$download_fox_only <- downloadHandler(
    filename = function() paste0("Fox_Detections_", Sys.Date(), ".txt"),
    content  = function(file) {
      df <- master_detections()
      write.table(df %>% filter(Detection == "Fox"),
                  file, sep = "\t", row.names = FALSE, quote = FALSE)
    }
  )

  # ── Daily counts reactive ─────────────────────────────────────────────
  daily_counts_rv <- reactiveVal(NULL)

  # ── Generate daily camera-presence counts ─────────────────────────────
  observeEvent(input$generate_daily_counts, {
    df <- master_detections()

    if (is.null(df) || nrow(df) == 0) {
      showNotification(
        "No processed detections found. Please process camera files first.",
        type = "error", duration = 5)
      return()
    }

    df <- df %>%
      filter(Detection == "Fox") %>%
      mutate(Date_parsed = suppressWarnings(lubridate::dmy(Date))) %>%
      filter(!is.na(Date_parsed))

    if (nrow(df) == 0) {
      showNotification("No valid fox detections found after parsing dates.",
                       type = "error", duration = 5)
      return()
    }

    # Full date sequence from first to last detection
    date_range <- seq(min(df$Date_parsed), max(df$Date_parsed), by = "day")

    # For each day: count unique cameras with >= 1 fox detection (binary per camera)
    daily <- df %>%
      group_by(Date_parsed) %>%
      summarise(Count = n_distinct(SiteID), .groups = "drop")

    # Fill missing days with zero
    daily_full <- data.frame(Date_parsed = date_range) %>%
      left_join(daily, by = "Date_parsed") %>%
      mutate(Count = tidyr::replace_na(Count, 0L)) %>%
      arrange(Date_parsed)

    daily_counts_rv(daily_full)

    showNotification(
      paste0("Generated ", nrow(daily_full), " daily counts across ",
             length(unique(df$SiteID)), " camera stations. ",
             "Download the file and upload it to the Abundance Model tab."),
      type = "message", duration = 8)
  })

  # ── Daily counts preview UI ───────────────────────────────────────────
  output$daily_counts_preview <- renderUI({
    df <- daily_counts_rv()
    if (is.null(df)) return(NULL)
    tagList(
      hr(),
      p(style = "font-size:13px; color:#4a5568;",
        strong(paste0(nrow(df), " daily values generated")),
        paste0(" — date range: ",
               format(min(df$Date_parsed), "%d/%m/%Y"),
               " to ",
               format(max(df$Date_parsed), "%d/%m/%Y"),
               ". Non-detection days filled with 0.
               Approx. ",
               floor(nrow(df) / 10),
               " ten-day blocks available for modelling.")),
      DTOutput("daily_counts_table")
    )
  })

  output$daily_counts_table <- renderDT({
    df <- daily_counts_rv()
    if (is.null(df)) return(NULL)
    datatable(
      head(df %>% rename(Date = Date_parsed, `Cameras Detecting Fox` = Count),
           20),
      options  = list(pageLength = 10, dom = "tp", scrollX = TRUE),
      rownames = FALSE,
      caption  = "Preview (first 20 rows)"
    ) %>%
      formatStyle("Cameras Detecting Fox",
                  background = styleColorBar(c(0, max(df$Count)), "#aed6f1"),
                  backgroundSize = "100% 80%",
                  backgroundRepeat = "no-repeat",
                  backgroundPosition = "center")
  })

  # ── Download daily counts for model ──────────────────────────────────
  output$download_daily_counts <- downloadHandler(
    filename = function()
      paste0("DailyCounts_ForModel_", Sys.Date(), ".txt"),
    content = function(file) {
      df <- daily_counts_rv()
      if (is.null(df)) {
        showNotification("Generate daily counts first.", type = "error")
        return()
      }
      writeLines(as.character(df$Count), file)
    }
  )

  # ── High-resolution detection plot exports ───────────────────────────
  # Shared ggplot2 builder for static export (mirrors renderPlotly logic)
  detection_ggplot <- reactive({
    df <- master_detections()
    req(nrow(df) > 0)
    df <- df %>% filter(Detection == "Fox")
    df$SiteID <- factor(df$SiteID)
    ggplot(df, aes(x = SiteID)) +
      geom_bar(fill = "#e67e22", width = 0.65) +
      theme_minimal(base_size = 14) +
      labs(title = "Total Fox Detections per Camera Station",
           x = "Camera Station", y = "Detection Count") +
      theme(
        plot.title   = element_text(face = "bold", size = 16,
                                    color = "#0d1b2a", hjust = 0.5),
        axis.text.x  = element_text(angle = 45, hjust = 1, size = 12),
        axis.title   = element_text(face = "bold", size = 13),
        panel.grid.major.x = element_blank()
      )
  })

  output$download_detection_png <- downloadHandler(
    filename = function()
      paste0("FoxWatch_Detections_", Sys.Date(), ".png"),
    content = function(file) {
      p <- detection_ggplot()
      req(!is.null(p))
      ggplot2::ggsave(file, plot = p,
                      width = 12, height = 6,
                      dpi = 300, units = "in",
                      bg = "white")
    }
  )

  output$download_detection_pdf <- downloadHandler(
    filename = function()
      paste0("FoxWatch_Detections_", Sys.Date(), ".pdf"),
    content = function(file) {
      p <- detection_ggplot()
      req(!is.null(p))
      ggplot2::ggsave(file, plot = p,
                      width = 12, height = 6,
                      device = "pdf", units = "in")
    }
  )

  # =========================================================================
  # -- Fox Abundance Model: effort-scaled SCM with partial pooling ------
  # =========================================================================
  #
  # JAGS specification (built dynamically from the UI hyperparameter inputs):
  #
  #   p      ~ dbeta(mb_pa, mb_pb)          # default Beta(1, 19)
  #   lambda ~ dgamma(mb_shape, mb_rate)    # default Gamma(2, 0.25)
  #   for (i in 1:nOcc) {
  #     N[i] ~ dpois(lambda)
  #     for (j in 1:nRep) {
  #       p_eff[i,j] <- 1 - (1 - p)^effort[i,j]
  #       y[i,j]     ~ dbin(p_eff[i,j], N[i])
  #     }
  #   }
  #
  # Data reshaping: the same daily-counts vector uploaded on the sidebar is
  # cut into 10-day blocks (matrix(counts, nrow = 10)), then transposed so
  # each row is one occasion (10-day block) and each column is one day
  # within that occasion. Effort matrix is user-uploaded (same dimensions)
  # or defaulted to 1s.
  # =========================================================================

  model_b_rv <- reactiveValues(
    N_summary     = NULL,
    param_summary = NULL,
    diag          = NULL,
    n_occ         = 0L,
    n_rep         = 0L,
    start_date    = as.Date(NA),
    n_cameras     = NA_integer_,
    summary_text  = NULL
  )

  observeEvent(input$run_model_b, {
    req(input$abundance_csv)

    withProgress(message = "Running the abundance model...",
                 value = 0.05, {

      # Detect input file type: TXT = one integer per line;
      # CSV = SiteID, Date, Time, Detection (as on the Spatial
      # Analysis tab). When CSV, aggregate to daily unique-camera
      # counts and record the first-day date so the plot can use
      # a real date axis.
      file_name  <- input$abundance_csv$name
      is_csv     <- grepl("\\.csv$", file_name, ignore.case = TRUE)
      start_date <- as.Date(NA)

      if (is_csv) {
        fox_df <- tryCatch(
          readr::read_csv(input$abundance_csv$datapath,
                          show_col_types = FALSE),
          error = function(e) NULL
        )
        if (is.null(fox_df) || nrow(fox_df) == 0 ||
            !all(c("SiteID", "Date", "Detection") %in% names(fox_df))) {
          showNotification(
            "CSV must have columns SiteID, Date, Detection.",
            type = "error", duration = 6)
          return()
        }
        fox_df <- fox_df %>%
          dplyr::filter(Detection == "Fox") %>%
          dplyr::mutate(Date = suppressWarnings(lubridate::dmy(Date))) %>%
          dplyr::filter(!is.na(Date))
        if (nrow(fox_df) == 0) {
          showNotification(
            "No valid Fox detections found in the uploaded CSV.",
            type = "error", duration = 6)
          return()
        }
        date_range <- seq(min(fox_df$Date), max(fox_df$Date), by = "day")
        daily <- fox_df %>%
          dplyr::group_by(Date) %>%
          dplyr::summarise(Count = dplyr::n_distinct(SiteID),
                           .groups = "drop")
        daily_full <- data.frame(Date = date_range) %>%
          dplyr::left_join(daily, by = "Date") %>%
          dplyr::mutate(Count = tidyr::replace_na(Count, 0L)) %>%
          dplyr::arrange(Date)
        counts     <- daily_full$Count
        start_date <- min(daily_full$Date)
        showNotification(
          paste0("CSV parsed: ", nrow(daily_full), " days from ",
                 format(start_date, "%d %b %Y"), " to ",
                 format(max(daily_full$Date), "%d %b %Y"),
                 " across ", dplyr::n_distinct(fox_df$SiteID),
                 " camera stations."),
          type = "message", duration = 6)
      } else {
        counts <- read_lines(input$abundance_csv$datapath)
        counts <- suppressWarnings(as.numeric(counts))
        counts <- counts[!is.na(counts)]
      }

      if (length(counts) < 10) {
        showNotification("Need at least 10 daily counts to fit the model.",
                         type = "error", duration = 5)
        return()
      }

      # -- Camera-array-size data-quality check (new in v1.5) --
      n_cameras <- suppressWarnings(as.integer(input$mb_n_cameras))
      if (is.na(n_cameras) || n_cameras < 1) n_cameras <- 1L
      max_daily <- max(counts)
      if (max_daily > n_cameras) {
        showNotification(
          paste0("Data-quality warning: the maximum daily count is ",
                 max_daily, ", but you entered only ", n_cameras,
                 " active camera(s). A daily count above the number of cameras usually indicates a data-entry issue (e.g. duplicated SiteID rows or wrong aggregation)."),
          type = "warning", duration = 12)
      }

      usable_length <- floor(length(counts) / 10) * 10
      detection_vec <- counts[1:usable_length]
      days_mat      <- matrix(detection_vec, nrow = 10)
      y_mat         <- t(days_mat)
      n_occ         <- nrow(y_mat)
      n_rep         <- ncol(y_mat)

      incProgress(0.10, detail = "Loading effort matrix...")

      effort_mat <- NULL
      if (!is.null(input$mb_effort_csv)) {
        effort_df <- tryCatch(
          utils::read.csv(input$mb_effort_csv$datapath,
                          header = FALSE, stringsAsFactors = FALSE),
          error = function(e) NULL
        )
        if (!is.null(effort_df) &&
            nrow(effort_df) == n_occ &&
            ncol(effort_df) == n_rep) {
          effort_mat <- as.matrix(effort_df)
          storage.mode(effort_mat) <- "double"
          if (any(!is.finite(effort_mat) | effort_mat < 0)) {
            showNotification(
              "Effort CSV contains non-finite or negative values; falling
               back to effort = 1.",
              type = "warning", duration = 6)
            effort_mat <- NULL
          }
        } else {
          showNotification(
            paste0("Effort CSV dimensions (",
                   ifelse(is.null(effort_df), "?", nrow(effort_df)), " x ",
                   ifelse(is.null(effort_df), "?", ncol(effort_df)),
                   ") do not match daily-count blocks (",
                   n_occ, " x ", n_rep,
                   "). Falling back to constant effort = 1."),
            type = "warning", duration = 8)
        }
      }
      if (is.null(effort_mat))
        effort_mat <- matrix(1, nrow = n_occ, ncol = n_rep)

      incProgress(0.20, detail = "Building JAGS model...")

      modelString <- paste0(
        "model {\n",
        "  p      ~ dbeta(",  input$mb_pa,    ", ", input$mb_pb,   ")\n",
        "  lambda ~ dgamma(", input$mb_shape, ", ", input$mb_rate, ")\n",
        "  for (i in 1:nOcc) {\n",
        "    N[i] ~ dpois(lambda)\n",
        "    for (j in 1:nRep) {\n",
        "      p_eff[i, j] <- 1 - pow(1 - p, effort[i, j])\n",
        "      y[i, j]     ~ dbin(p_eff[i, j], N[i])\n",
        "    }\n",
        "  }\n",
        "}\n"
      )
      writeLines(modelString, con = "model_b.txt")

      dataList <- list(
        y      = y_mat,
        effort = effort_mat,
        nOcc   = n_occ,
        nRep   = n_rep
      )

      incProgress(0.35, detail = "Adapting MCMC...")

      make_inits <- function(chain_idx) {
        set.seed(42L + chain_idx)
        list(
          N      = apply(y_mat, 1, max) + sample(0:3, n_occ, replace = TRUE),
          lambda = runif(1, min = 1,   max = 30),
          p      = runif(1, min = 0.02, max = 0.30)
        )
      }
      n_chains  <- max(1L, as.integer(input$mb_chains))
      inits_lst <- lapply(seq_len(n_chains), make_inits)

      jm <- tryCatch(
        rjags::jags.model(
          file     = "model_b.txt",
          data     = dataList,
          inits    = inits_lst,
          n.chains = n_chains,
          n.adapt  = as.integer(input$mb_adapt),
          quiet    = TRUE
        ),
        error = function(e) {
          showNotification(paste("Model failed to compile:", e$message),
                           type = "error", duration = 8)
          NULL
        }
      )
      req(!is.null(jm))

      incProgress(0.55, detail = "Burning in...")
      update(jm, as.integer(input$mb_burnin), progress.bar = "none")

      incProgress(0.75, detail = "Sampling from posterior...")
      samples <- coda.samples(
        jm,
        variable.names = c("lambda", "p", "N"),
        n.iter         = as.integer(input$mb_iter),
        thin           = max(1L, as.integer(input$mb_thin)),
        progress.bar   = "none"
      )

      incProgress(0.90, detail = "Summarising posterior...")

      post   <- as.matrix(samples)
      N_cols <- grep("^N\\[", colnames(post))

      N_summary <- do.call(rbind, lapply(seq_along(N_cols), function(i) {
        x   <- post[, N_cols[i]]
        hpd <- coda::HPDinterval(coda::mcmc(x), prob = 0.95)
        mode_est <- tryCatch(
          as.numeric(modeest::mlv(x, method = "mfv")[[1]]),
          error = function(e)
            as.numeric(names(sort(table(x), decreasing = TRUE)[1]))
        )
        data.frame(
          Occasion  = i,
          DayNumber = i * 10L,
          Mean      = mean(x),
          Median    = as.numeric(median(x)),
          Mode      = mode_est,
          SD        = sd(x),
          HDI_Lower = hpd[1, "lower"],
          HDI_Upper = hpd[1, "upper"],
          stringsAsFactors = FALSE
        )
      }))

      hpd_lambda <- coda::HPDinterval(coda::mcmc(post[, "lambda"]),
                                      prob = 0.95)
      hpd_p      <- coda::HPDinterval(coda::mcmc(post[, "p"]),
                                      prob = 0.95)
      param_summary <- data.frame(
        Parameter = c("lambda", "p"),
        Mean      = c(mean(post[, "lambda"]),   mean(post[, "p"])),
        Median    = c(median(post[, "lambda"]), median(post[, "p"])),
        SD        = c(sd(post[, "lambda"]),     sd(post[, "p"])),
        HDI_Lower = c(hpd_lambda[1, "lower"],   hpd_p[1, "lower"]),
        HDI_Upper = c(hpd_lambda[1, "upper"],   hpd_p[1, "upper"]),
        stringsAsFactors = FALSE
      )

      gelman_res <- tryCatch(
        coda::gelman.diag(samples, multivariate = FALSE, autoburnin = FALSE),
        error = function(e) NULL
      )
      rhat <- if (!is.null(gelman_res)) {
        setNames(round(gelman_res$psrf[, "Point est."], 3),
                 rownames(gelman_res$psrf))
      } else {
        setNames(rep(NA_real_, ncol(post)), colnames(post))
      }
      ess  <- round(as.numeric(coda::effectiveSize(samples)))
      names(ess) <- colnames(post)
      mcse <- apply(post, 2, function(x) sd(x) / sqrt(length(x)))

      diag_df <- data.frame(
        Parameter = colnames(post),
        Rhat      = rhat[colnames(post)],
        ESS       = ess[colnames(post)],
        MCSE      = round(mcse, 4),
        stringsAsFactors = FALSE
      )

      incProgress(1.00, detail = "Complete!")

      # -- Attach per-occasion R-hat & convergence category --
      # For each N[i], pull the R-hat computed above and classify.
      n_rhat_names        <- paste0("N[", seq_along(N_cols), "]")
      N_summary$Rhat_N    <- as.numeric(rhat[n_rhat_names])
      N_summary$Convergence <- factor(
        ifelse(is.na(N_summary$Rhat_N),        "Unknown",
          ifelse(N_summary$Rhat_N <= 1.05,     "Good",
            ifelse(N_summary$Rhat_N <= 1.10,   "Acceptable",
                                               "Poor"))),
        levels = c("Good", "Acceptable", "Poor", "Unknown")
      )

      # Attach real-date columns when we have a start_date
      if (!is.na(start_date)) {
        N_summary$BlockStart <- start_date +
          (N_summary$Occasion - 1L) * 10L
        N_summary$BlockEnd   <- start_date +
          N_summary$Occasion * 10L - 1L
        N_summary$BlockMid   <- start_date +
          (N_summary$Occasion - 1L) * 10L + 4L
      }

      # -- Foxes-per-camera derived columns (new in v1.5) --
      N_summary$FoxPerCamera_Mode <- round(N_summary$Mode / n_cameras, 3)
      N_summary$FoxPerCamera_Mean <- round(N_summary$Mean / n_cameras, 3)

      # -- Build the study-wide summary text (v1.5) --
      mean_abundance <- round(mean(N_summary$Mode), 1)
      fox_per_cam    <- round(mean_abundance / n_cameras, 2)
      date_bit <- if (!is.na(start_date)) {
        paste0(" over the period ",
               format(start_date, "%d %b %Y"), " to ",
               format(max(N_summary$BlockEnd), "%d %b %Y"))
      } else { "" }
      model_b_rv$summary_text <- HTML(paste0(
        "<div style=\"padding:10px 14px; font-size:14px; color:#0d1b2a;\">",
        "<strong>Camera array:</strong> ", n_cameras,
        " active camera station(s).<br>",
        "<strong>Fitted blocks:</strong> ", nrow(N_summary),
        " ten-day occasion(s)", date_bit, ".<br>",
        "<strong>Study-wide mean posterior mode:</strong> ",
        mean_abundance, " foxes across the array &rarr; approximately ",
        "<strong>", fox_per_cam, " foxes per camera</strong>.",
        "</div>"
      ))

      model_b_rv$n_cameras     <- n_cameras
      model_b_rv$N_summary     <- N_summary
      model_b_rv$param_summary <- param_summary
      model_b_rv$diag          <- diag_df
      model_b_rv$n_occ         <- n_occ
      model_b_rv$n_rep         <- n_rep
      model_b_rv$start_date    <- start_date
    })

    showNotification(
      paste0("Model completed: fitted ", model_b_rv$n_occ,
             " occasions x ", model_b_rv$n_rep, " replicates."),
      type = "message", duration = 6)

    if (!is.null(model_b_rv$diag)) {
      n_poor <- sum(model_b_rv$diag$Rhat > 1.10, na.rm = TRUE)
      if (n_poor > 0) {
        showNotification(
          paste0(n_poor,
                 " parameter(s) show poor MCMC convergence (R-hat > 1.10). ",
                 "Consider increasing burn-in, samples, or number of chains."),
          type = "warning", duration = 10)
      }
    }
  })

  # -- Model B: Study-wide summary text (new in v1.5) --
  output$mb_summary_text <- renderUI({
    txt <- model_b_rv$summary_text
    if (is.null(txt)) {
      return(HTML(paste0(
        "<div style=\"padding:10px 14px; color:#4a5568; font-style:italic;\">",
        "Run the model to see a study-wide summary.",
        "</div>"
      )))
    }
    txt
  })

  output$mb_N_table <- renderDT({
    df <- model_b_rv$N_summary
    if (is.null(df) || nrow(df) == 0) {
      return(datatable(
        data.frame(Message = "Run the model to see per-occasion estimates."),
        options = list(dom = "t")
      ))
    }
    datatable(
      df,
      options  = list(pageLength = 12, dom = "tp", scrollX = TRUE),
      rownames = FALSE,
      caption  = "N[i] posterior summaries -- one row per 10-day block."
    ) %>%
      formatRound(c("Mean", "Median", "Mode", "SD",
                    "HDI_Lower", "HDI_Upper"), 2) %>%
      formatRound(c("FoxPerCamera_Mode", "FoxPerCamera_Mean"), 3)
  })

  output$mb_param_table <- renderDT({
    df <- model_b_rv$param_summary
    if (is.null(df)) {
      return(datatable(
        data.frame(Message = "Run the model to see lambda, p estimates."),
        options = list(dom = "t")
      ))
    }
    datatable(
      df,
      options  = list(dom = "t"),
      rownames = FALSE,
      caption  = "Global parameters shared across all occasions."
    ) %>%
      formatRound(c("Mean", "Median", "SD",
                    "HDI_Lower", "HDI_Upper"), 3)
  })

  output$mb_diag_table <- renderDT({
    df <- model_b_rv$diag
    if (is.null(df) || nrow(df) == 0) {
      return(datatable(
        data.frame(Message = "Run the model to see convergence diagnostics."),
        options = list(dom = "t")
      ))
    }
    datatable(
      df,
      options  = list(pageLength = 20, dom = "tp", scrollX = TRUE),
      rownames = FALSE
    ) %>%
      formatRound("Rhat", 3) %>%
      formatRound("MCSE", 4) %>%
      formatStyle(
        "Rhat",
        backgroundColor = styleInterval(
          c(1.05, 1.10),
          c("#d4edda", "#fff3cd", "#f8d7da")
        )
      )
  })

  mb_mode_ggplot <- reactive({
    df <- model_b_rv$N_summary
    req(!is.null(df) && nrow(df) > 0)
    use_dates   <- "BlockMid" %in% names(df)

    # Palette for convergence traffic-light and plot backdrops
    conv_cols <- c(
      "Good"       = "#27AE60",
      "Acceptable" = "#F39C12",
      "Poor"       = "#E74C3C",
      "Unknown"    = "#95A5A6"
    )
    ribbon_fill <- "#7EB6D9"
    hdi_edge    <- "#2E5C82"
    trend_col   <- "#0F2A44"   # LOESS trend line - dark navy

    # LOESS span - tighter for long series, wider for short ones
    span_use <- max(0.15, min(0.5, 30 / nrow(df)))

    if (use_dates) {
      p <- ggplot(df, aes(x = BlockMid, y = Mode)) +
        geom_ribbon(aes(ymin = HDI_Lower, ymax = HDI_Upper),
                    fill = ribbon_fill, alpha = 0.6) +
        geom_line(aes(y = HDI_Lower), colour = hdi_edge,
                  linetype = "22", linewidth = 0.55) +
        geom_line(aes(y = HDI_Upper), colour = hdi_edge,
                  linetype = "22", linewidth = 0.55) +
        geom_smooth(method = "loess", span = span_use, se = FALSE,
                    colour = trend_col, linewidth = 0.9) +
        geom_point(aes(fill = Convergence, shape = Convergence),
                   colour = trend_col, size = 2.2, stroke = 0.35) +
        scale_x_date(date_labels = "%d %b %Y",
                     date_breaks = if (nrow(df) > 12) "2 months" else "1 month") +
        labs(x = "10-day block (midpoint date)")
    } else {
      n_occ    <- nrow(df)
      mid_occ  <- ceiling(n_occ / 2)
      brk_x    <- unique(c(1L, as.integer(mid_occ), as.integer(n_occ)))
      lab_x    <- c("Start of\ndeployment", "Mid",
                    "End of\ndeployment")[seq_along(brk_x)]
      p <- ggplot(df, aes(x = Occasion, y = Mode)) +
        geom_ribbon(aes(ymin = HDI_Lower, ymax = HDI_Upper),
                    fill = ribbon_fill, alpha = 0.6) +
        geom_line(aes(y = HDI_Lower), colour = hdi_edge,
                  linetype = "22", linewidth = 0.55) +
        geom_line(aes(y = HDI_Upper), colour = hdi_edge,
                  linetype = "22", linewidth = 0.55) +
        geom_smooth(method = "loess", span = span_use, se = FALSE,
                    colour = trend_col, linewidth = 0.9) +
        geom_point(aes(fill = Convergence, shape = Convergence),
                   colour = trend_col, size = 2.2, stroke = 0.35) +
        scale_x_continuous(breaks = brk_x, labels = lab_x,
                           expand = expansion(mult = c(0.02, 0.02))) +
        labs(x = NULL)
    }
    p +
      scale_fill_manual(name   = "Convergence\n(R-hat)",
                        values = conv_cols, drop = FALSE) +
      scale_shape_manual(name  = "Convergence\n(R-hat)",
                         values = c("Good"       = 21,
                                    "Acceptable" = 22,
                                    "Poor"       = 24,
                                    "Unknown"    = 23),
                         drop   = FALSE) +
      # Force every convergence category to appear in the legend with a
      # coloured symbol, even if the current data has no points at that level.
      guides(
        fill = guide_legend(
          override.aes = list(
            shape  = c(21, 22, 24, 23),
            fill   = c("#27AE60", "#F39C12", "#E74C3C", "#95A5A6"),
            colour = "#0F2A44",
            size   = 3.2,
            stroke = 0.4
          )
        ),
        shape = "none"
      ) +
      scale_y_continuous(expand = expansion(mult = c(0.01, 0.05))) +
      labs(
        title    = "Posterior Mode Abundance",
        subtitle = "Shaded band = 95% HDI (dashed edges); points coloured & shaped by R-hat convergence.",
        y        = "Estimated Abundance (N)",
        caption  = "FoxWatch — Bayesian SCM with effort-scaled p and shared λ"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title         = element_text(face = "bold", hjust = 0.5,
                                          colour = "#0d1b2a", size = 16),
        plot.subtitle      = element_text(hjust = 0.5, colour = "#4a5568",
                                          size = 11),
        plot.caption       = element_text(hjust = 0.5, colour = "#7a8a9a",
                                          size = 9, face = "italic"),
        axis.title         = element_text(face = "bold", size = 13,
                                          colour = "#0d1b2a"),
        axis.text          = element_text(colour = "#0d1b2a"),
        axis.text.x        = element_text(angle = if (use_dates) 45 else 0,
                                          hjust = if (use_dates) 1 else 0.5,
                                          lineheight = 0.9),
        panel.grid.major.y = element_line(colour = "grey88", linewidth = 0.35),
        panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank(),
        plot.margin        = ggplot2::margin(15, 20, 12, 15),
        legend.position    = "right",
        legend.title       = element_text(face = "bold", size = 11),
        legend.text        = element_text(size = 10)
      )
  })

  mb_mean_ggplot <- reactive({
    df <- model_b_rv$N_summary
    req(!is.null(df) && nrow(df) > 0)
    use_dates   <- "BlockMid" %in% names(df)

    conv_cols <- c(
      "Good"       = "#27AE60",
      "Acceptable" = "#F39C12",
      "Poor"       = "#E74C3C",
      "Unknown"    = "#95A5A6"
    )
    ribbon_fill <- "#F3BE86"   # soft peach for mean plot
    hdi_edge    <- "#C46A2E"
    trend_col   <- "#0F2A44"   # dark navy (consistent with mode plot)

    span_use <- max(0.15, min(0.5, 30 / nrow(df)))

    if (use_dates) {
      p <- ggplot(df, aes(x = BlockMid, y = Mean)) +
        geom_ribbon(aes(ymin = HDI_Lower, ymax = HDI_Upper),
                    fill = ribbon_fill, alpha = 0.6) +
        geom_line(aes(y = HDI_Lower), colour = hdi_edge,
                  linetype = "22", linewidth = 0.55) +
        geom_line(aes(y = HDI_Upper), colour = hdi_edge,
                  linetype = "22", linewidth = 0.55) +
        geom_smooth(method = "loess", span = span_use, se = FALSE,
                    colour = trend_col, linewidth = 0.9) +
        geom_point(aes(fill = Convergence, shape = Convergence),
                   colour = trend_col, size = 2.2, stroke = 0.35) +
        scale_x_date(date_labels = "%d %b %Y",
                     date_breaks = if (nrow(df) > 12) "2 months" else "1 month") +
        labs(x = "10-day block (midpoint date)")
    } else {
      n_occ    <- nrow(df)
      mid_occ  <- ceiling(n_occ / 2)
      brk_x    <- unique(c(1L, as.integer(mid_occ), as.integer(n_occ)))
      lab_x    <- c("Start of\ndeployment", "Mid",
                    "End of\ndeployment")[seq_along(brk_x)]
      p <- ggplot(df, aes(x = Occasion, y = Mean)) +
        geom_ribbon(aes(ymin = HDI_Lower, ymax = HDI_Upper),
                    fill = ribbon_fill, alpha = 0.6) +
        geom_line(aes(y = HDI_Lower), colour = hdi_edge,
                  linetype = "22", linewidth = 0.55) +
        geom_line(aes(y = HDI_Upper), colour = hdi_edge,
                  linetype = "22", linewidth = 0.55) +
        geom_smooth(method = "loess", span = span_use, se = FALSE,
                    colour = trend_col, linewidth = 0.9) +
        geom_point(aes(fill = Convergence, shape = Convergence),
                   colour = trend_col, size = 2.2, stroke = 0.35) +
        scale_x_continuous(breaks = brk_x, labels = lab_x,
                           expand = expansion(mult = c(0.02, 0.02))) +
        labs(x = NULL)
    }
    p +
      scale_fill_manual(name   = "Convergence\n(R-hat)",
                        values = conv_cols, drop = FALSE) +
      scale_shape_manual(name  = "Convergence\n(R-hat)",
                         values = c("Good"       = 21,
                                    "Acceptable" = 22,
                                    "Poor"       = 24,
                                    "Unknown"    = 23),
                         drop   = FALSE) +
      # Force every convergence category to appear in the legend with a
      # coloured symbol, even if the current data has no points at that level.
      guides(
        fill = guide_legend(
          override.aes = list(
            shape  = c(21, 22, 24, 23),
            fill   = c("#27AE60", "#F39C12", "#E74C3C", "#95A5A6"),
            colour = "#0F2A44",
            size   = 3.2,
            stroke = 0.4
          )
        ),
        shape = "none"
      ) +
      scale_y_continuous(expand = expansion(mult = c(0.01, 0.05))) +
      labs(
        title    = "Posterior Mean Abundance",
        subtitle = "Shaded band = 95% HDI (dashed edges); points coloured & shaped by R-hat convergence.",
        y        = "Estimated Abundance (N)",
        caption  = "FoxWatch — Bayesian SCM with effort-scaled p and shared λ"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title         = element_text(face = "bold", hjust = 0.5,
                                          colour = "#0d1b2a", size = 16),
        plot.subtitle      = element_text(hjust = 0.5, colour = "#4a5568",
                                          size = 11),
        plot.caption       = element_text(hjust = 0.5, colour = "#7a8a9a",
                                          size = 9, face = "italic"),
        axis.title         = element_text(face = "bold", size = 13,
                                          colour = "#0d1b2a"),
        axis.text          = element_text(colour = "#0d1b2a"),
        axis.text.x        = element_text(angle = if (use_dates) 45 else 0,
                                          hjust = if (use_dates) 1 else 0.5,
                                          lineheight = 0.9),
        panel.grid.major.y = element_line(colour = "grey88", linewidth = 0.35),
        panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank(),
        plot.margin        = ggplot2::margin(15, 20, 12, 15),
        legend.position    = "right",
        legend.title       = element_text(face = "bold", size = 11),
        legend.text        = element_text(size = 10)
      )
  })

  output$mb_mode_plot <- renderPlotly({
    if (is.null(model_b_rv$N_summary)) return(plotly_empty())
    ggplotly(mb_mode_ggplot())
  })

  output$mb_mean_plot <- renderPlotly({
    if (is.null(model_b_rv$N_summary)) return(plotly_empty())
    ggplotly(mb_mean_ggplot())
  })

  output$mb_download_mode_png <- downloadHandler(
    filename = function()
      paste0("FoxWatch_ModelB_Mode_", Sys.Date(), ".png"),
    content = function(file) {
      p <- mb_mode_ggplot()
      req(!is.null(p))
      ggplot2::ggsave(file, plot = p, width = 12, height = 6,
                      dpi = 300, units = "in", bg = "white")
    }
  )
  output$mb_download_mode_pdf <- downloadHandler(
    filename = function()
      paste0("FoxWatch_ModelB_Mode_", Sys.Date(), ".pdf"),
    content = function(file) {
      p <- mb_mode_ggplot()
      req(!is.null(p))
      ggplot2::ggsave(file, plot = p, width = 12, height = 6,
                      device = "pdf", units = "in")
    }
  )
  output$mb_download_mean_png <- downloadHandler(
    filename = function()
      paste0("FoxWatch_ModelB_Mean_", Sys.Date(), ".png"),
    content = function(file) {
      p <- mb_mean_ggplot()
      req(!is.null(p))
      ggplot2::ggsave(file, plot = p, width = 12, height = 6,
                      dpi = 300, units = "in", bg = "white")
    }
  )
  output$mb_download_mean_pdf <- downloadHandler(
    filename = function()
      paste0("FoxWatch_ModelB_Mean_", Sys.Date(), ".pdf"),
    content = function(file) {
      p <- mb_mean_ggplot()
      req(!is.null(p))
      ggplot2::ggsave(file, plot = p, width = 12, height = 6,
                      device = "pdf", units = "in")
    }
  )

  output$mb_download_csv <- downloadHandler(
    filename = function()
      paste0("FoxWatch_ModelB_Abundance_", Sys.Date(), ".csv"),
    content = function(file) {
      df <- model_b_rv$N_summary
      if (is.null(df)) {
        showNotification("Run the model first.", type = "error")
        return()
      }
      utils::write.csv(df, file, row.names = FALSE)
    }
  )
  output$mb_download_params_csv <- downloadHandler(
    filename = function()
      paste0("FoxWatch_ModelB_Parameters_", Sys.Date(), ".csv"),
    content = function(file) {
      df <- model_b_rv$param_summary
      if (is.null(df)) {
        showNotification("Run the model first.", type = "error")
        return()
      }
      utils::write.csv(df, file, row.names = FALSE)
    }
  )

  # ── Map section ────────────────────────────────────────────────────────
  observeEvent(input$make_map, {
    req(input$fox_csv, input$sites_csv, input$park_shp)

    withProgress(message = "Creating map visualisation...", value = 0.1, {
      tryCatch({
        fox   <- readr::read_csv(input$fox_csv$datapath,
                                 show_col_types = FALSE)
        sites <- readr::read_csv(input$sites_csv$datapath,
                                 show_col_types = FALSE)
        incProgress(0.3, detail = "Processing detection data...")

        fox <- fox %>% mutate(
          Date = suppressWarnings(lubridate::dmy(Date)),
          Time = suppressWarnings(hms::as_hms(Time)),
          Year = lubridate::year(Date)
        )
        fox_sites  <- fox %>% left_join(sites, by = "SiteID")
        min_y      <- as.integer(input$min_year)
        drop_y     <- suppressWarnings(as.integer(input$drop_year))
        years_keep <- sort(unique(
          fox_sites$Year[fox_sites$Year >= min_y]))
        if (!is.na(drop_y)) years_keep <- setdiff(years_keep, drop_y)
        validate(need(length(years_keep) > 0,
                      "No years remain after filtering."))

        site_levels <- fox_sites %>%
          filter(Year %in% years_keep) %>%
          distinct(SiteID) %>% arrange(SiteID) %>% pull(SiteID)
        fox_sites <- fox_sites %>%
          mutate(SiteID = factor(SiteID, levels = site_levels))
        N         <- length(site_levels)
        site_cols <- setNames(
          qualitative_hcl(N, palette = "Dynamic", c = 90, l = 55),
          site_levels
        )

        incProgress(0.5, detail = "Loading spatial layers...")
        park <- tryCatch(read_vector_layer(input$park_shp),
                         error = function(e) NULL)
        if (!is.null(park)) {
          if (is.na(sf::st_crs(park))) {
            sf::st_crs(park) <- 4326
          } else if (!isTRUE(sf::st_crs(park)$epsg == 4326)) {
            park <- sf::st_transform(park, 4326)
          }
        }

        incProgress(0.7, detail = "Creating visualisation...")
        detections_year <- fox_sites %>%
          filter(Year %in% years_keep, Detection == "Fox") %>%
          group_by(SiteID, Year, latitude, longitude) %>%
          summarise(TotalDetections = n(), .groups = "drop") %>%
          mutate(Year   = factor(Year,   levels = years_keep),
                 SiteID = factor(SiteID, levels = site_levels))
        validate(need(nrow(detections_year) > 0,
                      "No fox detections found for selected years."))

        detections_sf <- st_as_sf(
          detections_year,
          coords = c("longitude", "latitude"), crs = 4326)

        ref_label   <- "Site Locations"
        year_levels <- c(ref_label, as.character(years_keep))
        detections_sf <- detections_sf %>%
          dplyr::mutate(Year = factor(as.character(Year),
                                      levels = year_levels))
        sites_sf <- st_as_sf(sites,
                             coords = c("longitude", "latitude"),
                             crs    = 4326) %>%
          dplyr::mutate(Year = factor(ref_label, levels = year_levels))

        if (!is.null(park)) {
          bb <- sf::st_bbox(park)
        } else {
          validate(need(
            all(c("latitude", "longitude") %in% names(sites)),
            "Sites CSV must contain 'latitude' and 'longitude' columns."))
          pts <- sf::st_as_sf(sites,
                              coords = c("longitude", "latitude"),
                              crs = 4326, remove = FALSE)
          bb  <- sf::st_bbox(pts)
          showNotification(
            "Park layer could not be opened; using site extent instead.",
            type = "warning")
        }
        xlim_bb <- c(bb["xmin"], bb["xmax"])
        ylim_bb <- c(bb["ymin"], bb["ymax"])

        sites_lab       <- sf::st_transform(sites_sf, 3577)
        lab_xy          <- sf::st_coordinates(sites_lab)
        sites_lab       <- cbind(sites_lab,
                                 lab_x = lab_xy[, 1], lab_y = lab_xy[, 2])
        sites_lab       <- sf::st_transform(sites_lab, 4326)
        lab_xy_ll       <- sf::st_coordinates(sites_lab)
        sites_lab$lab_lon <- lab_xy_ll[, 1]
        sites_lab$lab_lat <- lab_xy_ll[, 2]

        max_det <- max(detections_sf$TotalDetections, na.rm = TRUE)
        brks    <- c(1, 5, 10, 20, 50, 100)
        brks    <- brks[brks <= max_det]
        if (length(brks) == 0) brks <- max_det

        incProgress(0.9, detail = "Rendering map...")
        p <- ggplot()
        if (!is.null(park))
          p <- p + geom_sf(data = park, fill = "darkgrey", color = "black")

        p <- p +
          geom_sf(data  = detections_sf,
                  aes(size = TotalDetections, colour = TotalDetections),
                  alpha = 0.98) +
          geom_sf(data   = sites_sf,
                  shape  = 24, fill = "darkblue", color = "black",
                  size   = 3.2, stroke = 0.6, inherit.aes = FALSE) +
          ggrepel::geom_text_repel(
            data               = sites_lab,
            aes(x = lab_lon, y = lab_lat, label = SiteID),
            seed               = 42, size = 3,
            box.padding        = grid::unit(0.18, "lines"),
            point.padding      = grid::unit(0.10, "lines"),
            min.segment.length = 0,
            inherit.aes        = FALSE
          ) +
          scale_size_continuous(range = c(2, 10)) +
          scale_colour_viridis_c(option = "inferno", direction = 1,
                                 begin = 0.05, end = 0.95,
                                 trans = "sqrt", breaks = brks) +
          facet_wrap(~ Year, nrow = 1, scales = "fixed") +
          coord_sf(xlim = xlim_bb, ylim = ylim_bb, expand = FALSE) +
          annotation_scale(location = "bl",
                           pad_x = grid::unit(0.1, "cm"),
                           pad_y = grid::unit(0.1, "cm")) +
          annotation_north_arrow(
            location = "tl", which_north = "true",
            style    = north_arrow_fancy_orienteering,
            height   = grid::unit(9, "mm"), width = grid::unit(9, "mm"),
            pad_x    = grid::unit(-0.3, "cm"),
            pad_y    = grid::unit(0.2, "cm")
          ) +
          theme_minimal(base_size = 14) +
          theme(plot.margin = ggplot2::margin(5.5, 5.5, 5.5, 7, "pt")) +
          labs(title = "Fox Detections Per Site Per Year",
               size = "Detections", colour = "Detections",
               x = NULL, y = NULL)

        map_plot_obj(p)
        incProgress(1, detail = "Map complete!")
        showNotification("Map generated successfully.",
                         type = "message", duration = 3)

      }, error = function(e) {
        showNotification(paste("Error creating map:", e$message),
                         type = "error", duration = 5)
      })
    })
  })

  output$detections_map <- renderPlot({
    p <- map_plot_obj()
    if (is.null(p)) {
      plot(1, type = "n", xlab = "", ylab = "",
           main = "Upload spatial data files and click 'Generate Map'")
    } else { p }
  })

  output$download_map_png <- downloadHandler(
    filename = function()
      sprintf("FoxWatch_Map_%s.png", Sys.Date()),
    content = function(file) {
      p <- map_plot_obj()
      if (!is.null(p))
        ggplot2::ggsave(file, plot = p, width = 14, height = 6, dpi = 300)
    }
  )
}

shinyApp(ui = ui, server = server)
