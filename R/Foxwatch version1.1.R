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
#  Version : 1.1
#  Authors : Duminda S.B. Dissanayake & Graeme Armstrong
#  Year    : 2026
#  Contact : dumie.dissanayake@dcceew.nsw.gov.au
# -----------------------------------------------------------------------------
#  Programme    : NSW Saving Our Species
#  Paper        : <DOI link when available>
#  Repository   : <GitHub URL>
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
#      • Bayesian SCM fitted via JAGS with two-chain MCMC sampling
#      • Gelman-Rubin convergence diagnostics (R-hat) per 10-day block
#      • Colour-coded results table and uncertainty ribbon plots
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
          fileInput("abundance_csv", "Daily Detection Counts (TXT)",
                    accept = ".txt"),
          helpText("One integer count per line")
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
                as zero, and outputs one integer per line \u2014 the exact
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
        
        # ── TAB 3: Fox Abundance Model ─────────────────────────────────
        tabPanel(
          title = tagList(icon("chart-line"), " Fox Abundance Model"),
          br(),
          div(class = "section-card",
              h3(icon("calculator"), " Fox Abundance Estimation"),
              br(),
              p(style = "font-size:14px; color:#4a5568;",
                "FoxWatch fits a Bayesian Simultaneous-Count Model (SCM;
                Armstrong & McSorley 2024) to daily fox detection counts.
                Two independent MCMC chains are run per 10-day block;
                Gelman-Rubin R-hat statistics are computed automatically
                to assess convergence."),
              hr(),
              actionButton("run_abundance_model",
                           label = tagList(icon("play-circle"),
                                           " Run Abundance Model"),
                           class = "btn btn-success btn-lg"),
              br(), br(),
              div(class = "convergence-legend",
                  tags$strong("Convergence key \u2014 Gelman-Rubin R-hat:"),
                  tags$ul(
                    tags$li(HTML(
                      "<span style='background:#d4edda; padding:2px 10px;
                       border-radius:4px; font-size:13px;'>
                       <strong>Good</strong> (&le; 1.05) \u2014
                       chains agree; estimate reliable</span>"
                    )),
                    tags$li(HTML(
                      "<span style='background:#fff3cd; padding:2px 10px;
                       border-radius:4px; font-size:13px;'>
                       <strong>Acceptable</strong> (1.05 &ndash; 1.10) \u2014
                       minor discrepancy; use with caution</span>"
                    )),
                    tags$li(HTML(
                      "<span style='background:#f8d7da; padding:2px 10px;
                       border-radius:4px; font-size:13px;'>
                       <strong>Poor</strong> (&gt; 1.10) \u2014
                       treat estimate with caution; check for sparse
                       detection counts in that period</span>"
                    ))
                  )
              ),
              br(),
              withSpinner(DTOutput("abundance_preview"), color = "#e67e22")
          ),
          br(),
          div(class = "section-card",
              h3(icon("chart-area"), " Abundance Visualisation"),
              withSpinner(plotlyOutput("abundance_plot", height = "580px"),
                          color = "#e67e22"),
              br(),
              div(style = "display:flex; justify-content:space-between; align-items:center; flex-wrap:wrap; gap:8px;",
                  downloadButton("download_abundance_data",
                                 label = tagList(icon("file-csv"),
                                                 " Download Abundance Estimates (CSV)"),
                                 class = "btn btn-primary"),
                  div(
                    tags$small(style = "color:#718096; margin-right:8px;",
                               "Publication-quality export:"),
                    downloadButton("download_abundance_png",
                                   label = tagList(icon("image"), " PNG (300 dpi)"),
                                   class = "btn btn-primary"),
                    span(" "),
                    downloadButton("download_abundance_pdf",
                                   label = tagList(icon("file-pdf"), " PDF"),
                                   class = "btn btn-primary")
                  )
              )
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
# Server  (unchanged from v1.0 — all logic identical)
# ---------------------------------------------------------------------------
server <- function(input, output, session) {
  
  camera_files      <- reactiveValues(files = list(), cleaned = list())
  master_detections <- reactiveVal(data.frame())
  summary_table_rv  <- reactiveVal(data.frame())
  abundance_results <- reactiveVal()
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
        paste0(" \u2014 date range: ",
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
  
  # ── Abundance model ───────────────────────────────────────────────────
  observeEvent(input$run_abundance_model, {
    req(input$abundance_csv)
    
    withProgress(message = "Running fox abundance model...", value = 0.1, {
      
      counts <- read_lines(input$abundance_csv$datapath)
      counts <- as.numeric(counts)
      counts <- counts[!is.na(counts)]
      
      if (length(counts) < 10) {
        showNotification("Need at least 10 daily counts", type = "error")
        return()
      }
      
      usable_length <- floor(length(counts) / 10) * 10
      detection_vec <- counts[1:usable_length]
      days          <- matrix(detection_vec, nrow = 10)
      
      incProgress(0.2, detail = "Building Bayesian SCM...")
      
      modelString <- "
      model {
        for (i in 1:length(y)) {
          y[i] ~ dbin(p, N)
        }
        p      ~ dbeta(2, 2)
        N      ~ dpois(lambda)
        lambda ~ dgamma(2, 0.05)
      }"
      writeLines(modelString, con = "model.txt")
      
      incProgress(0.4, detail = "Running MCMC chains...")
      
      jagsModel <- lapply(1:ncol(days), function(a) {
        y <- days[, a]
        jags.model(
          file     = "model.txt",
          data     = list(y = y),
          inits    = list(
            list(lambda = 100, p = 0.5),
            list(lambda = 50,  p = 0.3)
          ),
          n.chains = 2,
          n.adapt  = 1000,
          quiet    = TRUE
        )
      })
      
      incProgress(0.6, detail = "Sampling from posterior...")
      
      codaSamples <- lapply(jagsModel, function(mod) {
        update(mod, 2000, progress.bar = "none")
        coda.samples(mod,
                     variable.names = c("N", "p"),
                     n.iter         = 10000,
                     progress.bar   = "none")
      })
      
      incProgress(0.8, detail = "Summarising results and checking convergence...")
      
      modN <- sapply(codaSamples, function(samp) {
        modeest::mlv(as.vector(as.matrix(samp)[, "N"]), method = "mfv")[[1]]
      })
      upperhdi <- sapply(codaSamples, function(samp) {
        HPDinterval(mcmc(as.vector(as.matrix(samp)[, "N"])))[1, 2]
      })
      lowerhdi <- sapply(codaSamples, function(samp) {
        HPDinterval(mcmc(as.vector(as.matrix(samp)[, "N"])))[1, 1]
      })
      rhat_N <- sapply(codaSamples, function(samp) {
        tryCatch({
          coda::gelman.diag(samp)$psrf["N", "Point est."]
        }, error = function(e) NA_real_)
      })
      converged <- dplyr::case_when(
        is.na(rhat_N)  ~ "Unknown",
        rhat_N <= 1.05 ~ "Good",
        rhat_N <= 1.10 ~ "Acceptable",
        TRUE           ~ "Poor"
      )
      
      incProgress(0.9, detail = "Finalising...")
      
      abundance_results(data.frame(
        Period      = 1:length(modN),
        Abundance   = modN,
        LowerHDI    = lowerhdi,
        UpperHDI    = upperhdi,
        Rhat        = round(rhat_N, 3),
        Convergence = converged
      ))
      
      incProgress(1, detail = "Complete!")
    })
    
    showNotification("Abundance model completed successfully.",
                     type    = "message",
                     duration = 5)
    
    res    <- abundance_results()
    n_poor <- sum(res$Convergence == "Poor", na.rm = TRUE)
    if (n_poor > 0) {
      showNotification(
        paste0(n_poor, " block(s) showed poor MCMC convergence (R-hat > 1.10). ",
               "Abundance estimates for those periods should be treated with ",
               "caution. Review detection data for sparse or zero counts."),
        type     = "warning",
        duration = 12
      )
    }
  })
  
  # ── Abundance preview table ───────────────────────────────────────────
  output$abundance_preview <- renderDT({
    df <- abundance_results()
    if (is.null(df) || nrow(df) == 0) {
      return(datatable(
        data.frame(Message = "Run the model to see abundance estimates"),
        options = list(dom = "t")
      ))
    }
    datatable(
      head(df, 10),
      options  = list(pageLength = 10, dom = "t", scrollX = TRUE),
      rownames = FALSE
    ) %>%
      formatStyle(
        "Convergence",
        backgroundColor = styleEqual(
          c("Good",    "Acceptable", "Poor",    "Unknown"),
          c("#d4edda", "#fff3cd",    "#f8d7da", "#e2e3e5")
        )
      ) %>%
      formatRound(columns = c("Abundance", "LowerHDI", "UpperHDI", "Rhat"),
                  digits  = 2)
  })
  
  # ── Abundance plot ────────────────────────────────────────────────────
  output$abundance_plot <- renderPlotly({
    df <- abundance_results()
    if (is.null(df) || nrow(df) == 0) return(plotly_empty())
    
    df$DayNumber <- df$Period * 10
    conv_colours <- c(
      "Good"       = "#27ae60",
      "Acceptable" = "#f39c12",
      "Poor"       = "#e74c3c",
      "Unknown"    = "#95a5a6"
    )
    
    p <- ggplot(df, aes(x = DayNumber)) +
      geom_ribbon(aes(ymin = LowerHDI, ymax = UpperHDI),
                  fill = "#aed6f1", alpha = 0.35) +
      geom_line(aes(y = UpperHDI),
                color = "#7fb3d3", linetype = "dashed", linewidth = 0.5) +
      geom_line(aes(y = LowerHDI),
                color = "#7fb3d3", linetype = "dashed", linewidth = 0.5) +
      geom_smooth(aes(y = Abundance),
                  method = "loess", se = FALSE,
                  color = "#1b3a5c", linewidth = 1.2, span = 0.2) +
      geom_point(aes(y = Abundance, colour = Convergence), size = 3) +
      scale_colour_manual(
        name   = "Convergence\n(R-hat)",
        values = conv_colours
      ) +
      labs(
        title = "Fox Abundance — 10-Day Sampling Periods",
        x     = "Day",
        y     = "Estimated Abundance (N)"
      ) +
      scale_x_continuous(breaks = seq(0, max(df$DayNumber), by = 50)) +
      theme_minimal(base_family = "sans") +
      theme(
        plot.title   = element_text(face = "bold", size = 15,
                                    hjust = 0.5, color = "#0d1b2a"),
        axis.title   = element_text(face = "bold", size = 12),
        axis.text.x  = element_text(size = 10, angle = 45, hjust = 1),
        axis.text.y  = element_text(size = 11),
        legend.title = element_text(size = 11, face = "bold"),
        panel.grid.minor = element_blank()
      )
    ggplotly(p)
  })
  
  # ── Abundance download ─────────────────────────────────────────────────
  output$download_abundance_data <- downloadHandler(
    filename = function()
      paste0("FoxWatch_Abundance_", Sys.Date(), ".csv"),
    content = function(file)
      write.csv(abundance_results(), file, row.names = FALSE)
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
  
  # ── High-resolution abundance plot exports ────────────────────────────
  abundance_ggplot <- reactive({
    df <- abundance_results()
    req(!is.null(df) && nrow(df) > 0)
    df$DayNumber <- df$Period * 10
    conv_colours <- c(
      "Good"       = "#27ae60",
      "Acceptable" = "#f39c12",
      "Poor"       = "#e74c3c",
      "Unknown"    = "#95a5a6"
    )
    ggplot(df, aes(x = DayNumber)) +
      geom_ribbon(aes(ymin = LowerHDI, ymax = UpperHDI),
                  fill = "#aed6f1", alpha = 0.35) +
      geom_line(aes(y = UpperHDI),
                color = "#7fb3d3", linetype = "dashed", linewidth = 0.5) +
      geom_line(aes(y = LowerHDI),
                color = "#7fb3d3", linetype = "dashed", linewidth = 0.5) +
      geom_smooth(aes(y = Abundance),
                  method = "loess", se = FALSE,
                  color = "#1b3a5c", linewidth = 1.2, span = 0.2) +
      geom_point(aes(y = Abundance, colour = Convergence), size = 2.5) +
      scale_colour_manual(
        name   = "Convergence
(R-hat)",
        values = conv_colours
      ) +
      labs(
        title = "Fox Abundance — 10-Day Sampling Periods",
        x     = "Day",
        y     = "Estimated Abundance (N)"
      ) +
      scale_x_continuous(breaks = seq(0, max(df$DayNumber), by = 50)) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title   = element_text(face = "bold", size = 16,
                                    hjust = 0.5, color = "#0d1b2a"),
        axis.title   = element_text(face = "bold", size = 13),
        axis.text.x  = element_text(size = 11, angle = 45, hjust = 1),
        axis.text.y  = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold"),
        panel.grid.minor = element_blank()
      )
  })
  
  output$download_abundance_png <- downloadHandler(
    filename = function()
      paste0("FoxWatch_Abundance_", Sys.Date(), ".png"),
    content = function(file) {
      p <- abundance_ggplot()
      req(!is.null(p))
      ggplot2::ggsave(file, plot = p,
                      width = 14, height = 6,
                      dpi = 300, units = "in",
                      bg = "white")
    }
  )
  
  output$download_abundance_pdf <- downloadHandler(
    filename = function()
      paste0("FoxWatch_Abundance_", Sys.Date(), ".pdf"),
    content = function(file) {
      p <- abundance_ggplot()
      req(!is.null(p))
      ggplot2::ggsave(file, plot = p,
                      width = 14, height = 6,
                      device = "pdf", units = "in")
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