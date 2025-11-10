# ==============================================================================
# Script: Install Required R Packages
# Purpose: Install all packages needed for GBIF Kenya bias assessment
# Author: Automated Research Pipeline
# Date: 2025-11-10
# ==============================================================================

cat("\n==================================================================\n")
cat("  Installing Required Packages for GBIF Kenya Analysis\n")
cat("==================================================================\n\n")

# List of required packages ---------------------------------------------------
required_packages <- c(
  # Data manipulation and utilities
  "tidyverse",        # Data wrangling and visualization
  "here",             # Path management
  "lubridate",        # Date handling
  "scales",           # Formatting numbers and dates

  # Spatial analysis
  "sf",               # Spatial data handling
  "terra",            # Raster data
  "rnaturalearth",    # Country boundaries
  "rnaturalearthdata",# Natural Earth data
  "spdep",            # Spatial dependence
  "ape",              # Phylogenetics and spatial stats

  # GBIF and biodiversity data
  "rgbif",            # GBIF API interface
  "CoordinateCleaner",# Coordinate cleaning
  "geodata",          # Geographic data download

  # Statistical analysis
  "MASS",             # Negative binomial models
  "lme4",             # Mixed models
  "mgcv",             # Generalized additive models
  "DHARMa",           # Model diagnostics
  "performance",      # Model performance metrics
  "effects",          # Model effect plots
  "MuMIn",            # Model selection

  # Diversity analysis
  "vegan",            # Community ecology
  "iNEXT",            # Interpolation/extrapolation
  "ineq",             # Inequality metrics (Gini)

  # Time series
  "Kendall",          # Mann-Kendall tests
  "forecast",         # Time series forecasting
  "tseries",          # Time series analysis

  # Visualization
  "ggplot2",          # Included in tidyverse but listed explicitly
  "viridis",          # Color palettes
  "patchwork",        # Combine plots
  "treemap",          # Treemap visualizations
  "ggalluvial",       # Alluvial plots

  # Reproducibility
  "renv",             # Package management
  "rmarkdown",        # Dynamic documents
  "knitr",            # Literate programming
  "kableExtra",       # Table formatting

  # Optional but useful
  "usethis",          # Project setup utilities
  "devtools"          # Package development tools
)

# GBIF and biodiversity data
if (!"occAssess" %in% installed.packages()) devtools::install_github("https://github.com/robboyd/occAssess") # Occurrence bias assessment

# Function to install packages ------------------------------------------------
install_if_missing <- function(packages) {
  installed <- installed.packages()[, "Package"]
  to_install <- packages[!packages %in% installed]

  if (length(to_install) == 0) {
    cat("All required packages are already installed!\n")
    return(invisible(NULL))
  }

  cat("The following packages will be installed:\n")
  cat(paste(" -", to_install, collapse = "\n"), "\n\n")

  cat("Installing", length(to_install), "packages...\n\n")

  # Try installing from CRAN
  for (pkg in to_install) {
    cat("Installing", pkg, "... ")

    tryCatch({
      install.packages(pkg, dependencies = TRUE, quiet = TRUE)
      cat("SUCCESS\n")
    }, error = function(e) {
      cat("FAILED\n")
      cat("  Error:", as.character(e), "\n")
    })
  }
}

# Install packages ------------------------------------------------------------
install_if_missing(required_packages)

# Verify installation ---------------------------------------------------------
cat("\n==================================================================\n")
cat("  Verifying Installation\n")
cat("==================================================================\n\n")

installed <- installed.packages()[, "Package"]
status <- required_packages %in% installed

summary_table <- data.frame(
  Package = required_packages,
  Status = ifelse(status, "Installed", "Missing")
)

missing_count <- sum(!status)
installed_count <- sum(status)

cat("Installed:", installed_count, "/", length(required_packages), "\n")
if (missing_count > 0) {
  cat("Missing:", missing_count, "\n\n")
  cat("Missing packages:\n")
  print(summary_table[!status, ], row.names = FALSE)
  cat("\nPlease install missing packages manually or try running this script again.\n")
} else {
  cat("\nAll packages installed successfully!\n")
}

# Initialize renv (optional) --------------------------------------------------
cat("\n==================================================================\n")
cat("  Initialize renv? (optional)\n")
cat("==================================================================\n\n")

cat("Would you like to initialize renv for reproducibility?\n")
cat("This will create a project-specific package library.\n\n")

if (interactive()) {
  response <- readline(prompt = "Initialize renv? (y/n): ")

  if (tolower(response) == "y") {
    if ("renv" %in% installed) {
      cat("\nInitializing renv...\n")
      renv::init()
      cat("\nrenv initialized successfully!\n")
      cat("Run renv::snapshot() to save the current package state.\n")
    } else {
      cat("\nrenv is not installed. Install it first:\n")
      cat("  install.packages('renv')\n")
    }
  } else {
    cat("\nSkipping renv initialization.\n")
  }
} else {
  cat("Run renv::init() manually if you want to use renv.\n")
}

# Session info ----------------------------------------------------------------
cat("\n==================================================================\n")
cat("  Session Information\n")
cat("==================================================================\n\n")

cat("R version:", R.version.string, "\n")
cat("Platform:", R.version$platform, "\n")
cat("Running under:", R.version$os, "\n\n")

cat("Installation complete!\n")
cat("You can now run the analysis scripts.\n\n")

cat("Quick start:\n")
cat("  source('run_all_analyses.R')\n\n")

cat("Or run scripts individually:\n")
cat("  source('scripts/01_data_download.R')\n")
cat("  source('scripts/02_spatial_bias.R')\n")
cat("  # ... and so on\n\n")

# Return summary invisibly
invisible(summary_table)
