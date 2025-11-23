# ==============================================================================
# Master Script: Run Complete GBIF Kenya Bias Assessment Pipeline
# Author: Kwiz Computing Technologies
# Date: 2025-11-10
# ==============================================================================

# This script runs the complete analysis pipeline in the correct order
# It includes error handling and progress reporting

cat("\n")
cat("==============================================================================\n")
cat("  GBIF Kenya Bias Assessment Pipeline\n")
cat("  Starting analysis at:", as.character(Sys.time()), "\n")
cat("==============================================================================\n\n")

# Set up environment ----------------------------------------------------------
cat("Setting up environment...\n")

# Load required packages
required_packages <- c("here", "tidyverse")
missing_packages <- required_packages[!required_packages %in% installed.packages()[,"Package"]]

if(length(missing_packages) > 0) {
  cat("Installing missing packages:", paste(missing_packages, collapse = ", "), "\n")
  install.packages(missing_packages)
}

library(here)
library(tidyverse)

# Source utility functions
source(here("R", "utils.R"))
source(here("scripts", "results_config.R"))

# Create log file
log_file <- here("analysis_log.txt")
sink(log_file, split = TRUE)

cat("Analysis started at:", as.character(Sys.time()), "\n")
cat("Working directory:", here(), "\n")
cat("R version:", R.version.string, "\n\n")

# Define analysis scripts -----------------------------------------------------
scripts <- c(
  "01_data_download.R",
  "01b_data_quality_assessment.R",
  "02_spatial_bias.R",
  "03_temporal_bias.R",
  "04_taxonomic_bias.R",
  "05_statistical_models.R"
)

script_names <- c(
  "Data Download and Cleaning",
  "Data Quality Assessment",
  "Spatial Bias Assessment",
  "Temporal Bias Assessment",
  "Taxonomic Bias Assessment",
  "Statistical Modeling"
)

# Define results directories and required files
results_configs <- list(
  list(dir = NULL, files = NULL),  # 01_data_download has no results folder yet
  list(dir = "data_quality", files = DATA_QUALITY_FILES),
  list(dir = "spatial_bias", files = SPATIAL_BIAS_FILES),
  list(dir = "temporal_bias", files = TEMPORAL_BIAS_FILES),
  list(dir = "taxonomic_bias", files = TAXONOMIC_BIAS_FILES),
  list(dir = "statistical_models", files = STATISTICAL_MODELS_FILES)
)

# Run analysis pipeline -------------------------------------------------------
n_scripts <- length(scripts)
success <- logical(n_scripts)
timings <- numeric(n_scripts)
errors <- vector("list", n_scripts)

for (i in seq_along(scripts)) {
  cat("\n")
  cat("==============================================================================\n")
  cat("Step", i, "of", n_scripts, ":", script_names[i], "\n")
  cat("==============================================================================\n")
  cat("Script:", scripts[i], "\n")
  cat("Started at:", as.character(Sys.time()), "\n\n")

  script_path <- here("scripts", scripts[i])

  # Check if script exists
  if (!file.exists(script_path)) {
    cat("ERROR: Script not found:", script_path, "\n")
    success[i] <- FALSE
    errors[[i]] <- paste("Script not found:", script_path)
    next
  }

  # Check if results already exist for this analysis
  config <- results_configs[[i]]
  skip_analysis <- FALSE

  if (!is.null(config$dir) && !is.null(config$files)) {
    results_path <- here("results", config$dir)
    if (check_results_exist(results_path, config$files)) {
      cat("✓ All required results already exist for", script_names[i], "\n")
      cat("  Skipping analysis. To re-run, delete files in:", results_path, "\n")
      skip_analysis <- TRUE
      success[i] <- TRUE
      errors[[i]] <- NA
      timings[i] <- 0
    }
  }

  if (!skip_analysis) {
    # Run script with error handling
    start_time <- Sys.time()

    tryCatch({
      source(script_path, echo = FALSE, verbose = FALSE)
      success[i] <- TRUE
      errors[[i]] <- NA
    }, error = function(e) {
      cat("\nERROR in", scripts[i], ":\n")
      cat(as.character(e), "\n")
      success[i] <<- FALSE
      errors[[i]] <<- as.character(e)
    }, warning = function(w) {
      cat("\nWARNING in", scripts[i], ":\n")
      cat(as.character(w), "\n")
    })

    end_time <- Sys.time()
    timings[i] <- as.numeric(difftime(end_time, start_time, units = "mins"))

    cat("\nCompleted at:", as.character(end_time), "\n")
    cat("Time elapsed:", round(timings[i], 2), "minutes\n")

    if (success[i]) {
      cat("Status: SUCCESS\n")
    } else {
      cat("Status: FAILED\n")
    }
  }
}

# Generate manuscript ---------------------------------------------------------
cat("\n")
cat("==============================================================================\n")
cat("Compiling R Markdown Manuscript\n")
cat("==============================================================================\n\n")

manuscript_path <- here("docs", "kenya_gbif_bias_assessment.Rmd")

if (file.exists(manuscript_path) && all(success)) {
  cat("Rendering manuscript to HTML...\n")

  tryCatch({
    rmarkdown::render(
      manuscript_path,
      output_format = "html_document",
      output_dir = here("docs"),
      quiet = FALSE
    )
    cat("Manuscript compiled successfully!\n")
    manuscript_success <- TRUE
  }, error = function(e) {
    cat("\nERROR rendering manuscript:\n")
    cat(as.character(e), "\n")
    manuscript_success <- FALSE
  })
} else if (!all(success)) {
  cat("Skipping manuscript compilation due to errors in analysis scripts\n")
  manuscript_success <- FALSE
} else {
  cat("Manuscript file not found:", manuscript_path, "\n")
  manuscript_success <- FALSE
}

# Summary Report --------------------------------------------------------------
cat("\n")
cat("==============================================================================\n")
cat("  PIPELINE SUMMARY\n")
cat("==============================================================================\n\n")

summary_table <- data.frame(
  Step = seq_along(scripts),
  Script = script_names,
  Status = ifelse(success, "SUCCESS", "FAILED"),
  Time_mins = round(timings, 2)
)

print(summary_table, row.names = FALSE)

cat("\nTotal scripts run:", n_scripts, "\n")
cat("Successful:", sum(success), "\n")
cat("Failed:", sum(!success), "\n")
cat("Total time:", round(sum(timings), 2), "minutes\n")

if (exists("manuscript_success")) {
  cat("Manuscript compilation:", ifelse(manuscript_success, "SUCCESS", "FAILED"), "\n")
}

# Report errors
if (any(!success)) {
  cat("\n")
  cat("==============================================================================\n")
  cat("  ERRORS\n")
  cat("==============================================================================\n\n")

  for (i in which(!success)) {
    cat("Script:", script_names[i], "\n")
    cat("Error:", errors[[i]], "\n\n")
  }
}

# Final message ---------------------------------------------------------------
cat("\n")
cat("==============================================================================\n")
cat("  PIPELINE COMPLETE\n")
cat("==============================================================================\n")
cat("Finished at:", as.character(Sys.time()), "\n")
cat("Log saved to:", log_file, "\n")

if (all(success)) {
  cat("\nAll analyses completed successfully!\n")
  cat("\nOutputs:\n")
  cat("  - Processed data: data/processed/\n")
  cat("  - Analysis results: data/outputs/\n")
  cat("  - Figures: figures/\n")
  cat("  - Manuscript: docs/\n")
} else {
  cat("\nSome analyses failed. Please check the log file for details.\n")
}

cat("\n")

# Close log file
sink()

# Return summary invisibly
invisible(list(
  success = success,
  timings = timings,
  errors = errors,
  summary_table = summary_table
))
