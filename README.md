# FoxWatch

**A User-Friendly Tool for Estimating Fox Abundance from Camera Traps Using Simultaneous-Count Models**

Version 1.1 • Authors: Duminda S.B. Dissanayake & Graeme Armstrong • 2026
Programme: NSW Saving Our Species

---

## Overview

FoxWatch is an R Shiny application that enables field ecologists and conservation
practitioners to estimate red fox (*Vulpes vulpes*) abundance from camera-trap
data using a Bayesian Simultaneous-Count Model (SCM), without requiring
statistical programming expertise.

Key features:

- Automated EXIF metadata parsing from camera-trap text files
- Interactive spatial visualisation of fox detections per site per year
- Bayesian SCM fitted via JAGS with two-chain MCMC sampling
- Gelman-Rubin convergence diagnostics (R-hat) per 10-day block
- Colour-coded results table and uncertainty ribbon plots
- One-click CSV export of all abundance estimates

## Repository structure

```
Foxwatch/
├── R/                          # Shiny application source
│   └── Foxwatch version1.1.R
├── data/
│   └── example/                # Example data shipped with the tool
│       ├── camera_txt/         # Per-camera EXIF text files (CamOol_*)
│       ├── detections/         # Fox detection records (CSV / TXT)
│       ├── sites/              # Site location coordinates (CSV)
│       └── park_boundary/      # Oolambeyan National Park shapefile
├── docs/                       # Supplementary documentation
│   └── Guide walks you through installing FoxWatch_SupplementaryFile_*.pdf
├── README.md
├── LICENSE                     # MIT
├── CITATION.cff                # Citation metadata
└── .gitignore
```

## Requirements

- **R** ≥ 4.2.0
- **JAGS** ≥ 4.3.1 — install separately from <https://mcmc-jags.sourceforge.io>
- The following R packages:

```r
install.packages(c(
  "shiny", "dplyr", "tidyr", "readr", "stringr", "DT", "tools",
  "ggplot2", "plotly", "shinycssloaders", "lubridate",
  "rjags", "modeest", "coda",
  "sf", "ggspatial", "viridis", "ggrepel", "hms",
  "colorspace", "RColorBrewer", "forcats"
))
```

Note: `grid` ships with base R and does not need to be installed separately.

## Running the app

From an R session in the repository root:

```r
shiny::runApp("R/Foxwatch version1.1.R")
```

Or open the file in RStudio and click **Run App**.

## Using the example data

Once the app is running, upload the files from `data/example/` via the
sidebar inputs:

| App input                          | File(s) to upload                                        |
| ---------------------------------- | -------------------------------------------------------- |
| Fox Detection Records (CSV)        | `data/example/detections/Foxdetection_Example_*.csv`     |
| Site Locations (CSV)               | `data/example/sites/SiteLocation_Example_*.csv`          |
| Park Boundary (shapefile)          | All files in `data/example/park_boundary/` (select all)  |
| Camera TXT Files (EXIF)            | All `CamOol_*.txt` files in `data/example/camera_txt/`   |
| Daily Detection Counts (TXT)       | `data/example/detections/Foxonlydata.txt`                |

A step-by-step walkthrough with screenshots is in
`docs/Guide walks you through installing FoxWatch_SupplementaryFile_*.pdf`.

## How to cite

If you use FoxWatch, please cite:

> Dissanayake D.S.B. & Armstrong G. (2026) *FoxWatch: A User-Friendly Tool for
> Estimating Fox Abundance from Camera Traps Using Simultaneous-Count Models.*
> [Journal name to be inserted on acceptance]. DOI: *pending*.

The simultaneous-count method itself is described in:

> Armstrong G. & McSorley A. (2024) Estimating fox abundance using the
> simultaneous count method. *Wildlife Letters*, **2**, 102–109.
> <https://doi.org/10.1002/wll2.12038>

A `CITATION.cff` file is included so GitHub renders a "Cite this repository"
button automatically.

## Data availability

The data in `data/example/` are example records from cameras at Oolambeyan
National Park, NSW. They are provided to allow reviewers and end-users to
reproduce the worked example in the manuscript and supplementary material.

## Licence

Released under the MIT Licence — see [`LICENSE`](LICENSE).

## Contact

Duminda S.B. Dissanayake — <dumie.dissanayake@dcceew.nsw.gov.au>
NSW Saving Our Species, Department of Climate Change, Energy, the Environment
and Water (DCCEEW).

---

*Repository status: private during peer review. It will be made public on
publication of the associated manuscript.*
