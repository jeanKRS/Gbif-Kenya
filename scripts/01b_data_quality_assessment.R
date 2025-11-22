# ==============================================================================
# Script: 01b_data_quality_assessment.R
# Purpose: Comprehensive data quality assessment and quantification
# Author: Automated Research Pipeline
# Date: 2025-11-11
# ==============================================================================
# This script performs detailed tracking and quantification of all data quality
# issues identified during the GBIF data cleaning process, documenting the
# number and percentage of records affected by each issue relative to the
# original dataset.
# ==============================================================================

# Load required packages -------------------------------------------------------
suppressPackageStartupMessages({
  library(rgbif)
  library(tidyverse)
  library(sf)
  library(here)
  library(CoordinateCleaner)
  library(lubridate)
  library(knitr)
})

# Source utility functions
source(here("R", "utils.R"))

# Setup paths ------------------------------------------------------------------
data_raw <- here("data", "raw")
data_processed <- here("data", "processed")
data_outputs <- here("data", "outputs")  # Keep for backward compatibility
figures_dir <- here("figures")
results_dir <- here("results", "data_quality")

# Create directories if needed
dir.create(data_outputs, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# Define required output files ------------------------------------------------
required_files <- c(
  "data_quality_assessment.rds",
  "data_quality_tracking.csv",
  "coordinate_issues_summary.csv"
)

# Check if analysis should be run ----------------------------------------------
if (!skip_if_complete("Data Quality Assessment", results_dir, required_files)) {

  # Load data ------------------------------------------------------------------
  message("Loading data...")
  kenya_raw <- readRDS(file.path(data_raw, "kenya_gbif_raw.rds"))
  kenya_final <- readRDS(file.path(data_processed, "kenya_gbif_clean.rds"))
  metadata <- readRDS(file.path(data_processed, "metadata.rds"))

  # Initialize tracking data frame ---------------------------------------------
  n_original <- nrow(kenya_raw)
  message("Original dataset: ", format(n_original, big.mark = ","), " records\n")

  # Create quality tracking data frame
  quality_tracking <- tibble(
    step = character(),
    description = character(),
    records_removed = numeric(),
    records_remaining = numeric(),
    percent_removed = numeric(),
    cumulative_percent_removed = numeric()
  )

  # Function to add tracking step
  add_step <- function(step_name, description, data_before, data_after) {
    n_before <- nrow(data_before)
    n_after <- nrow(data_after)
    n_removed <- n_before - n_after
    pct_removed <- (n_removed / n_original) * 100
    cumulative_removed <- ((n_original - n_after) / n_original) * 100

    quality_tracking <<- quality_tracking %>%
      add_row(
        step = step_name,
        description = description,
        records_removed = n_removed,
        records_remaining = n_after,
        percent_removed = pct_removed,
        cumulative_percent_removed = cumulative_removed
      )

    if (n_removed > 0) {
      message(sprintf("  %s: %s records removed (%.2f%%)",
                      step_name,
                      format(n_removed, big.mark = ","),
                      pct_removed))
    }

    return(data_after)
  }

  # Step 1: Missing coordinates ------------------------------------------------
  message("\n1. Checking for missing coordinates...")
  kenya_step1 <- kenya_raw %>%
    filter(!is.na(decimalLongitude), !is.na(decimalLatitude))

  kenya_raw <- add_step(
    "Missing coordinates",
    "Records without decimal latitude or longitude",
    kenya_raw,
    kenya_step1
  )

  # Step 2: High coordinate uncertainty ----------------------------------------
  message("\n2. Checking coordinate uncertainty...")

  # Analyze uncertainty distribution
  uncertainty_summary <- kenya_step1 %>%
    filter(!is.na(coordinateUncertaintyInMeters)) %>%
    summarise(
      n_with_uncertainty = n(),
      median_uncertainty = median(coordinateUncertaintyInMeters),
      mean_uncertainty = mean(coordinateUncertaintyInMeters),
      q25_uncertainty = quantile(coordinateUncertaintyInMeters, 0.25),
      q75_uncertainty = quantile(coordinateUncertaintyInMeters, 0.75),
      max_uncertainty = max(coordinateUncertaintyInMeters),
      n_gt_1km = sum(coordinateUncertaintyInMeters > 1000),
      n_gt_5km = sum(coordinateUncertaintyInMeters > 5000),
      n_gt_10km = sum(coordinateUncertaintyInMeters > 10000),
      n_gt_50km = sum(coordinateUncertaintyInMeters > 50000)
    )

  kenya_step2 <- kenya_step1 %>%
    filter(coordinateUncertaintyInMeters < 10000 | is.na(coordinateUncertaintyInMeters))

  kenya_step1 <- add_step(
    "High uncertainty",
    "Coordinate uncertainty > 10 km",
    kenya_step1,
    kenya_step2
  )

  # Step 3: Basis of record ----------------------------------------------------
  message("\n3. Checking basis of record...")

  # Summarize basis of record
  basis_summary <- kenya_step2 %>%
    count(basisOfRecord, sort = TRUE) %>%
    mutate(percent = (n / nrow(kenya_step2)) * 100)

  kenya_step3 <- kenya_step2 %>%
    filter(!basisOfRecord %in% c("FOSSIL_SPECIMEN", "LIVING_SPECIMEN"))

  kenya_step2 <- add_step(
    "Inappropriate basis",
    "Fossil or living specimens (not natural occurrences)",
    kenya_step2,
    kenya_step3
  )

  # Step 4: Species-level identification --------------------------------------
  message("\n4. Checking taxonomic identification...")

  # Analyze taxonomic completeness
  tax_completeness <- kenya_step3 %>%
    summarise(
      n_with_species = sum(!is.na(species) & species != ""),
      n_with_genus = sum(!is.na(genus) & genus != ""),
      n_with_family = sum(!is.na(family) & family != ""),
      n_with_order = sum(!is.na(order) & order != ""),
      n_with_class = sum(!is.na(class) & class != ""),
      n_with_phylum = sum(!is.na(phylum) & phylum != ""),
      n_with_kingdom = sum(!is.na(kingdom) & kingdom != "")
    ) %>%
    mutate(across(everything(), ~ (.x / nrow(kenya_step3)) * 100, .names = "pct_{.col}"))

  kenya_step4 <- kenya_step3 %>%
    filter(!is.na(species), species != "")

  kenya_step3 <- add_step(
    "Missing species ID",
    "Records without species-level identification",
    kenya_step3,
    kenya_step4
  )

  # Step 5: Date parsing and validation ----------------------------------------
  message("\n5. Parsing and validating dates...")

  kenya_step5 <- kenya_step4 %>%
    mutate(
      eventDate = ymd(eventDate),
      year = year(eventDate),
      month = month(eventDate),
      day = day(eventDate)
    )

  # Check for invalid dates
  date_issues <- kenya_step5 %>%
    summarise(
      n_missing_date = sum(is.na(eventDate)),
      n_missing_year = sum(is.na(year)),
      n_future_dates = sum(year > year(Sys.Date()), na.rm = TRUE),
      n_pre_1950 = sum(year < 1950, na.rm = TRUE),
      n_valid_year = sum(year >= 1950 & year <= year(Sys.Date()), na.rm = TRUE)
    )

  kenya_step6 <- kenya_step5 %>%
    filter(year >= 1950, year <= year(Sys.Date()))

  kenya_step5 <- add_step(
    "Invalid dates",
    "Records with dates before 1950 or in the future",
    kenya_step5,
    kenya_step6
  )

  # Step 6: Duplicate records --------------------------------------------------
  message("\n6. Identifying duplicate records...")

  # Identify duplicates
  duplicates <- kenya_step6 %>%
    group_by(species, decimalLongitude, decimalLatitude, eventDate) %>%
    filter(n() > 1) %>%
    ungroup()

  n_duplicate_sets <- duplicates %>%
    group_by(species, decimalLongitude, decimalLatitude, eventDate) %>%
    summarise(n = n(), .groups = "drop") %>%
    nrow()

  kenya_step7 <- kenya_step6 %>%
    distinct(species, decimalLongitude, decimalLatitude, eventDate, .keep_all = TRUE)

  kenya_step6 <- add_step(
    "Duplicates",
    "Exact duplicate records (same species, location, and date)",
    kenya_step6,
    kenya_step7
  )

  message("\n7. Applying CoordinateCleaner tests...")

  # Individual CoordinateCleaner tests -----------------------------------------
  # Run each test separately to track individual issues

  coord_issues <- list()

  # Test 1: Capitals
  message("  - Testing capitals...")
  flags_capitals <- cc_cap(
    kenya_step7,
    lon = "decimalLongitude",
    lat = "decimalLatitude",
    buffer = 10000,
    value = "flagged",
    verbose = FALSE
  )
  coord_issues$capitals <- sum(flags_capitals)

  # Test 2: Centroids
  message("  - Testing centroids...")
  flags_centroids <- cc_cen(
    kenya_step7,
    lon = "decimalLongitude",
    lat = "decimalLatitude",
    buffer = 5000,
    value = "flagged",
    verbose = FALSE
  )
  coord_issues$centroids <- sum(flags_centroids)

  # Test 3: Equal coordinates
  message("  - Testing equal coordinates...")
  flags_equal <- cc_equ(
    kenya_step7,
    lon = "decimalLongitude",
    lat = "decimalLatitude",
    value = "flagged",
    verbose = FALSE
  )
  coord_issues$equal <- sum(flags_equal)

  # Test 4: GBIF headquarters
  message("  - Testing GBIF headquarters...")
  flags_gbif <- cc_gbif(
    kenya_step7,
    lon = "decimalLongitude",
    lat = "decimalLatitude",
    value = "flagged",
    verbose = FALSE
  )
  coord_issues$gbif <- sum(flags_gbif)

  # Test 5: Zeros
  message("  - Testing zeros...")
  flags_zeros <- cc_zero(
    kenya_step7,
    lon = "decimalLongitude",
    lat = "decimalLatitude",
    buffer = 0.5,
    value = "flagged",
    verbose = FALSE
  )
  coord_issues$zeros <- sum(flags_zeros)

  # Test 6: Urban areas
  message("  - Testing urban areas...")
  flags_urban <- cc_urb(
    kenya_step7,
    lon = "decimalLongitude",
    lat = "decimalLatitude",
    value = "flagged",
    verbose = FALSE
  )
  coord_issues$urban <- sum(flags_urban)

  # Test 7: Outliers
  message("  - Testing outliers...")
  flags_outliers <- cc_outl(
    kenya_step7,
    lon = "decimalLongitude",
    lat = "decimalLatitude",
    species = "species",
    method = "quantile",
    mltpl = 5,
    tdi = 1000,
    min_occs = 10,
    value = "flagged",
    verbose = FALSE
  )
  coord_issues$outliers <- sum(flags_outliers)

  # Combined coordinate cleaning
  flags_all <- clean_coordinates(
    x = kenya_step7,
    lon = "decimalLongitude",
    lat = "decimalLatitude",
    species = "species",
    countries = "countryCode",
    tests = c("capitals", "centroids", "equal", "gbif", "outliers", "zeros", "urban"),
    capitals_rad = 10000,
    centroids_rad = 5000,
    zeros_rad = 0.5,
    outliers_method = "quantile",
    outliers_mtp = 5,
    outliers_td = 1000,
    outliers_size = 10,
    verbose = FALSE,
    value = "flagged"
  )

  kenya_step8 <- kenya_step7 %>%
    filter(!flags_all)

  kenya_step7 <- add_step(
    "Coordinate flags",
    "Records flagged by CoordinateCleaner (any test)",
    kenya_step7,
    kenya_step8
  )

  # Create detailed coordinate issues table
  # NOTE: Individual test counts may overlap - one record can be flagged by multiple tests
  # Therefore, sum(records_flagged) may exceed the total number of unique records flagged
  coord_issues_df <- tibble(
    test = c("Capitals", "Centroids", "Equal coordinates", "GBIF HQ",
             "Zeros", "Urban areas", "Outliers"),
    description = c(
      "Within 10 km of country/province capitals",
      "Within 5 km of country/province centroids",
      "Identical latitude and longitude values",
      "Coordinates matching GBIF headquarters",
      "Coordinates at (0,0) or near equator/prime meridian",
      "Coordinates in major urban centers",
      "Statistical outliers based on species distribution"
    ),
    records_flagged = c(
      coord_issues$capitals,
      coord_issues$centroids,
      coord_issues$equal,
      coord_issues$gbif,
      coord_issues$zeros,
      coord_issues$urban,
      coord_issues$outliers
    ),
    percent_of_original = (records_flagged / n_original) * 100
  ) %>%
    arrange(desc(records_flagged))

  # Verify final count matches
  if (nrow(kenya_step8) != nrow(kenya_final)) {
    warning("Final count mismatch! Expected: ", nrow(kenya_final),
            ", Got: ", nrow(kenya_step8))
  }

  # Summary statistics ---------------------------------------------------------
  message("\n=== Data Quality Summary ===")
  message("Original records: ", format(n_original, big.mark = ","))
  message("Final clean records: ", format(nrow(kenya_final), big.mark = ","))
  message("Total removed: ", format(n_original - nrow(kenya_final), big.mark = ","),
          " (", round(((n_original - nrow(kenya_final)) / n_original) * 100, 2), "%)")

  # Calculate additional quality metrics ---------------------------------------
  quality_metrics <- list(
    # Overall retention
    retention_rate = (nrow(kenya_final) / n_original) * 100,
    removal_rate = ((n_original - nrow(kenya_final)) / n_original) * 100,

    # Data completeness
    completeness = list(
      coordinates = (nrow(kenya_step1) / n_original) * 100,
      species_id = tax_completeness$pct_n_with_species,
      genus_id = tax_completeness$pct_n_with_genus,
      family_id = tax_completeness$pct_n_with_family,
      dates = (date_issues$n_valid_year / n_original) * 100
    ),

    # Coordinate quality
    coordinate_quality = list(
      uncertainty_summary = uncertainty_summary,
      coord_issues_summary = coord_issues_df,
      # Count of unique records flagged by ANY CoordinateCleaner test
      total_coord_flags = sum(flags_all),
      # Percentage of records at step 7 that were flagged (should be 0-100%)
      percent_coord_flags = (sum(flags_all) / nrow(kenya_step7)) * 100
    ),

    # Taxonomic completeness
    taxonomic_completeness = tax_completeness,

    # Temporal quality
    temporal_quality = date_issues,

    # Duplicates
    duplicate_info = list(
      n_duplicate_records = nrow(duplicates),
      n_duplicate_sets = n_duplicate_sets,
      percent_duplicates = (nrow(duplicates) / nrow(kenya_step6)) * 100
    ),

    # Basis of record
    basis_of_record = basis_summary
  )

  # Compile final results ------------------------------------------------------
  quality_results <- list(
    tracking = quality_tracking,
    coord_issues = coord_issues_df,
    metrics = quality_metrics,
    n_original = n_original,
    n_final = nrow(kenya_final),
    assessment_date = Sys.Date()
  )

  # Save results ---------------------------------------------------------------
  saveRDS(quality_results, file.path(results_dir, "data_quality_assessment.rds"))
  write_csv(quality_tracking, file.path(results_dir, "data_quality_tracking.csv"))
  write_csv(coord_issues_df, file.path(results_dir, "coordinate_issues_summary.csv"))

  # Also save to old location for backward compatibility
  saveRDS(quality_results, file.path(data_outputs, "data_quality_assessment.rds"))
  write_csv(quality_tracking, file.path(data_outputs, "data_quality_tracking.csv"))
  write_csv(coord_issues_df, file.path(data_outputs, "coordinate_issues_summary.csv"))

  message("\n=== Data quality assessment complete ===")
  message("Results saved to: ", results_dir)
}

# Load results for visualization -----------------------------------------------
if (!exists("quality_results")) {
  message("Loading saved data quality results...")
  quality_results <- readRDS(file.path(results_dir, "data_quality_assessment.rds"))
}

# Create visualizations --------------------------------------------------------
message("\nGenerating data quality visualizations...")

# 1. Waterfall/Cascade plot showing filtering steps
cascade_plot <- quality_results$tracking %>%
  mutate(
    step = fct_reorder(factor(step), row_number()),
    step_label = paste0(step, "\n(", format(records_removed, big.mark = ","), ")")
  ) %>%
  ggplot(aes(x = step, y = records_removed, fill = percent_removed)) +
  geom_col() +
  geom_text(aes(label = format(records_removed, big.mark = ",")),
            vjust = -0.5, size = 3) +
  scale_fill_viridis_c(option = "magma", direction = -1,
                       name = "% of original") +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = "Data Quality Filtering Cascade",
    subtitle = paste0("From ", format(quality_results$n_original, big.mark = ","),
                      " to ", format(quality_results$n_final, big.mark = ","),
                      " records (",
                      round(quality_results$metrics$retention_rate, 1), "% retained)"),
    x = "Quality Filter",
    y = "Records Removed"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "right"
  )

ggsave(file.path(figures_dir, "data_quality_cascade.png"),
       cascade_plot, width = 10, height = 6, dpi = 300)

# 2. Cumulative retention plot
retention_plot <- quality_results$tracking %>%
  mutate(
    step = fct_reorder(factor(step), row_number()),
    retention_pct = (records_remaining / quality_results$n_original) * 100
  ) %>%
  ggplot(aes(x = step, y = retention_pct, group = 1)) +
  geom_line(color = "#2C3E50", size = 1) +
  geom_point(color = "#E74C3C", size = 3) +
  geom_hline(yintercept = 100, linetype = "dashed", alpha = 0.5) +
  geom_text(aes(label = paste0(round(retention_pct, 1), "%")),
            vjust = -1, size = 3) +
  scale_y_continuous(limits = c(0, 105), breaks = seq(0, 100, 20)) +
  labs(
    title = "Cumulative Data Retention Through Quality Filters",
    subtitle = "Percentage of original records remaining after each filter",
    x = "Quality Filter",
    y = "Records Retained (%)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", size = 14)
  )

ggsave(file.path(figures_dir, "data_quality_retention.png"),
       retention_plot, width = 10, height = 6, dpi = 300)

# 3. CoordinateCleaner issues breakdown
coord_plot <- quality_results$coord_issues %>%
  filter(records_flagged > 0) %>%
  mutate(test = fct_reorder(test, records_flagged)) %>%
  ggplot(aes(x = test, y = records_flagged, fill = percent_of_original)) +
  geom_col() +
  geom_text(aes(label = paste0(format(records_flagged, big.mark = ","),
                                "\n(", round(percent_of_original, 2), "%)")),
            hjust = -0.1, size = 3) +
  scale_fill_viridis_c(option = "plasma", direction = -1,
                       name = "% of original") +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.15))) +
  coord_flip() +
  labs(
    title = "Coordinate Quality Issues Detected",
    subtitle = "Number of records flagged by each CoordinateCleaner test",
    x = "Quality Test",
    y = "Records Flagged"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "right"
  )

ggsave(file.path(figures_dir, "coordinate_issues_breakdown.png"),
       coord_plot, width = 10, height = 6, dpi = 300)

# 4. Summary pie chart
summary_data <- tibble(
  category = c("Retained (Clean)", "Removed (Quality Issues)"),
  count = c(quality_results$n_final,
            quality_results$n_original - quality_results$n_final),
  percentage = c(quality_results$metrics$retention_rate,
                 quality_results$metrics$removal_rate)
)

summary_pie <- summary_data %>%
  ggplot(aes(x = "", y = count, fill = category)) +
  geom_col(width = 1, color = "white") +
  geom_text(aes(label = paste0(format(count, big.mark = ","), "\n",
                                "(", round(percentage, 1), "%)")),
            position = position_stack(vjust = 0.5), size = 5) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = c("Retained (Clean)" = "#27AE60",
                                "Removed (Quality Issues)" = "#E74C3C")) +
  labs(
    title = "Overall Data Quality Assessment",
    subtitle = paste0("GBIF Kenya: ", format(quality_results$n_original, big.mark = ","),
                      " original records"),
    fill = ""
  ) +
  theme_void(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "bottom"
  )

ggsave(file.path(figures_dir, "data_quality_summary_pie.png"),
       summary_pie, width = 8, height = 8, dpi = 300)

message("Visualizations saved to: ", figures_dir)
message("\n=== All data quality assessments complete ===")
