# ==============================================================================
# Script: 01b_data_quality_assessment.R
# Purpose: Comprehensive data quality assessment with flagging (not filtering)
# Author: Kwiz Computing Technologies
# Date: 2025-11-22
# ==============================================================================
# This script identifies and FLAGS all data quality issues without removing
# records. Each downstream analysis can then filter based on relevant flags.
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
  "coordinate_issues_summary.csv",
  "kenya_gbif_flagged.rds"  # New: data with quality flags
)

# Check if analysis should be run ----------------------------------------------
if (!skip_if_complete("Data Quality Assessment", results_dir, required_files)) {

  message("=== Comprehensive Data Quality Assessment (Flagging Mode) ===\n")

  # Load data ------------------------------------------------------------------
  message("Loading data...")
  kenya_raw <- readRDS(file.path(data_raw, "kenya_gbif_raw.rds"))

  # Initialize tracking
  n_original <- nrow(kenya_raw)
  message("Original dataset: ", format(n_original, big.mark = ","), " records\n")

  # Create quality tracking data frame
  quality_tracking <- tibble(
    issue = character(),
    description = character(),
    records_flagged = numeric(),
    percent_flagged = numeric()
  )

  # Function to add tracking step
  add_flag_tracking <- function(issue_name, description, flag_vector) {
    n_flagged <- sum(flag_vector, na.rm = TRUE)
    pct_flagged <- (n_flagged / n_original) * 100

    quality_tracking <<- quality_tracking %>%
      add_row(
        issue = issue_name,
        description = description,
        records_flagged = n_flagged,
        percent_flagged = pct_flagged
      )

    if (n_flagged > 0) {
      message(sprintf("  %s: %s records flagged (%.2f%%)",
                      issue_name,
                      format(n_flagged, big.mark = ","),
                      pct_flagged))
    }

    return(flag_vector)
  }

  # Start with all data, add flag columns ---------------------------------------
  kenya_flagged <- kenya_raw

  # Flag 1: Missing coordinates ------------------------------------------------
  message("\n1. Flagging missing coordinates...")
  kenya_flagged$flag_missing_coords <- add_flag_tracking(
    "Missing coordinates",
    "Records without decimal latitude or longitude",
    is.na(kenya_raw$decimalLongitude) | is.na(kenya_raw$decimalLatitude)
  )

  # Flag 2: High coordinate uncertainty ----------------------------------------
  message("\n2. Flagging coordinate uncertainty...")

  # Analyze uncertainty distribution
  uncertainty_summary <- kenya_raw %>%
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

  kenya_flagged$flag_high_uncertainty <- add_flag_tracking(
    "High uncertainty",
    "Coordinate uncertainty > 10 km",
    !is.na(kenya_raw$coordinateUncertaintyInMeters) &
      kenya_raw$coordinateUncertaintyInMeters >= 10000
  )

  # Flag 3: Basis of record ----------------------------------------------------
  message("\n3. Flagging basis of record...")

  # Summarize basis of record
  basis_summary <- kenya_raw %>%
    count(basisOfRecord, sort = TRUE) %>%
    mutate(percent = (n / n_original) * 100)

  kenya_flagged$flag_inappropriate_basis <- add_flag_tracking(
    "Inappropriate basis",
    "Fossil or living specimens (not natural occurrences)",
    kenya_raw$basisOfRecord %in% c("FOSSIL_SPECIMEN", "LIVING_SPECIMEN")
  )

  # Flag 4: Species-level identification --------------------------------------
  message("\n4. Flagging taxonomic identification...")

  # Analyze taxonomic completeness
  tax_completeness <- kenya_raw %>%
    summarise(
      n_with_species = sum(!is.na(species) & species != ""),
      n_with_genus = sum(!is.na(genus) & genus != ""),
      n_with_family = sum(!is.na(family) & family != ""),
      n_with_order = sum(!is.na(order) & order != ""),
      n_with_class = sum(!is.na(class) & class != ""),
      n_with_phylum = sum(!is.na(phylum) & phylum != ""),
      n_with_kingdom = sum(!is.na(kingdom) & kingdom != "")
    ) %>%
    mutate(across(everything(), ~ (.x / n_original) * 100, .names = "pct_{.col}"))

  kenya_flagged$flag_missing_species <- add_flag_tracking(
    "Missing species ID",
    "Records without species-level identification",
    is.na(kenya_raw$species) | kenya_raw$species == ""
  )

  # Flag 5: Date parsing and validation ----------------------------------------
  message("\n5. Flagging date issues...")

  # Parse dates
  kenya_flagged <- kenya_flagged %>%
    mutate(
      eventDate = ymd(eventDate),
      year = year(eventDate),
      month = month(eventDate),
      day = day(eventDate)
    )

  # Check for invalid dates
  date_issues <- kenya_flagged %>%
    summarise(
      n_missing_date = sum(is.na(eventDate)),
      n_missing_year = sum(is.na(year)),
      n_future_dates = sum(year > year(Sys.Date()), na.rm = TRUE),
      n_pre_1950 = sum(year < 1950, na.rm = TRUE),
      n_valid_year = sum(year >= 1950 & year <= year(Sys.Date()), na.rm = TRUE)
    )

  kenya_flagged$flag_invalid_date <- add_flag_tracking(
    "Invalid dates",
    "Records with dates before 1950 or in the future",
    is.na(kenya_flagged$year) |
      kenya_flagged$year < 1950 |
      kenya_flagged$year > year(Sys.Date())
  )

  # Flag 6: Duplicate records --------------------------------------------------
  message("\n6. Flagging duplicate records...")

  # Identify duplicates
  kenya_flagged <- kenya_flagged %>%
    group_by(species, decimalLongitude, decimalLatitude, eventDate) %>%
    mutate(
      duplicate_group_size = n(),
      flag_duplicate = duplicate_group_size > 1
    ) %>%
    ungroup()

  duplicates_summary <- kenya_flagged %>%
    filter(flag_duplicate) %>%
    summarise(
      n_duplicate_records = n(),
      n_duplicate_sets = n_distinct(species, decimalLongitude, decimalLatitude, eventDate)
    )

  add_flag_tracking(
    "Duplicates",
    "Exact duplicate records (same species, location, and date)",
    kenya_flagged$flag_duplicate
  )

  # Flag 7: CoordinateCleaner tests --------------------------------------------
  message("\n7. Applying CoordinateCleaner tests...")

  # Need to work with records that have coordinates
  kenya_with_coords <- kenya_flagged %>%
    filter(!flag_missing_coords)

  coord_issues <- list()

  # Test 1: Capitals
  message("  - Testing capitals...")
  flags_capitals <- !cc_cap(
    kenya_with_coords,
    lon = "decimalLongitude",
    lat = "decimalLatitude",
    buffer = 10000,
    value = "flagged",
    verbose = FALSE
  )
  coord_issues$capitals <- sum(flags_capitals)

  # Test 2: Centroids
  message("  - Testing centroids...")
  flags_centroids <- !cc_cen(
    kenya_with_coords,
    lon = "decimalLongitude",
    lat = "decimalLatitude",
    buffer = 5000,
    value = "flagged",
    verbose = FALSE
  )
  coord_issues$centroids <- sum(flags_centroids)

  # Test 3: Equal coordinates
  message("  - Testing equal coordinates...")
  flags_equal <- !cc_equ(
    kenya_with_coords,
    lon = "decimalLongitude",
    lat = "decimalLatitude",
    value = "flagged",
    verbose = FALSE
  )
  coord_issues$equal <- sum(flags_equal)

  # Test 4: GBIF headquarters
  message("  - Testing GBIF headquarters...")
  flags_gbif <- !cc_gbif(
    kenya_with_coords,
    lon = "decimalLongitude",
    lat = "decimalLatitude",
    value = "flagged",
    verbose = FALSE
  )
  coord_issues$gbif <- sum(flags_gbif)

  # Test 5: Zeros
  message("  - Testing zeros...")
  flags_zeros <- !cc_zero(
    kenya_with_coords,
    lon = "decimalLongitude",
    lat = "decimalLatitude",
    buffer = 0.5,
    value = "flagged",
    verbose = FALSE
  )
  coord_issues$zeros <- sum(flags_zeros)

  # Test 6: Urban areas
  message("  - Testing urban areas...")
  flags_urban <- !cc_urb(
    kenya_with_coords,
    lon = "decimalLongitude",
    lat = "decimalLatitude",
    value = "flagged",
    verbose = FALSE
  )
  coord_issues$urban <- sum(flags_urban)

  # Test 7: Outliers
  message("  - Testing outliers...")
  flags_outliers <- !cc_outl(
    kenya_with_coords,
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

  # Add coordinate flags to main dataset
  # Initialize all coordinate flags as FALSE
  kenya_flagged$flag_coord_capitals <- FALSE
  kenya_flagged$flag_coord_centroids <- FALSE
  kenya_flagged$flag_coord_equal <- FALSE
  kenya_flagged$flag_coord_gbif <- FALSE
  kenya_flagged$flag_coord_zeros <- FALSE
  kenya_flagged$flag_coord_urban <- FALSE
  kenya_flagged$flag_coord_outliers <- FALSE

  # Set flags for records with coordinates
  valid_indices <- which(!kenya_flagged$flag_missing_coords)
  kenya_flagged$flag_coord_capitals[valid_indices] <- flags_capitals
  kenya_flagged$flag_coord_centroids[valid_indices] <- flags_centroids
  kenya_flagged$flag_coord_equal[valid_indices] <- flags_equal
  kenya_flagged$flag_coord_gbif[valid_indices] <- flags_gbif
  kenya_flagged$flag_coord_zeros[valid_indices] <- flags_zeros
  kenya_flagged$flag_coord_urban[valid_indices] <- flags_urban
  kenya_flagged$flag_coord_outliers[valid_indices] <- flags_outliers

  # Add tracking for coordinate flags
  add_flag_tracking("Coord: Capitals", "Within 10km of capitals", kenya_flagged$flag_coord_capitals)
  add_flag_tracking("Coord: Centroids", "Within 5km of centroids", kenya_flagged$flag_coord_centroids)
  add_flag_tracking("Coord: Equal", "Identical lat/lon", kenya_flagged$flag_coord_equal)
  add_flag_tracking("Coord: GBIF HQ", "At GBIF headquarters", kenya_flagged$flag_coord_gbif)
  add_flag_tracking("Coord: Zeros", "At (0,0) or near", kenya_flagged$flag_coord_zeros)
  add_flag_tracking("Coord: Urban", "In urban areas", kenya_flagged$flag_coord_urban)
  add_flag_tracking("Coord: Outliers", "Statistical outliers", kenya_flagged$flag_coord_outliers)

  # Create combined coordinate quality flag
  kenya_flagged$flag_any_coord_issue <- (
    kenya_flagged$flag_coord_capitals |
    kenya_flagged$flag_coord_centroids |
    kenya_flagged$flag_coord_equal |
    kenya_flagged$flag_coord_gbif |
    kenya_flagged$flag_coord_zeros |
    kenya_flagged$flag_coord_urban |
    kenya_flagged$flag_coord_outliers
  )

  add_flag_tracking("Any coord issue", "Any CoordinateCleaner flag", kenya_flagged$flag_any_coord_issue)

  # Create detailed coordinate issues table
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

  # Summary statistics ---------------------------------------------------------
  message("\n=== Data Quality Summary ===")
  message("Total records: ", format(n_original, big.mark = ","))

  # Calculate how many records would be retained with different filter strategies
  clean_counts <- kenya_flagged %>%
    summarise(
      strict_clean = sum(
        !flag_missing_coords & !flag_high_uncertainty &
        !flag_inappropriate_basis & !flag_missing_species &
        !flag_invalid_date & !flag_duplicate & !flag_any_coord_issue
      ),
      moderate_clean = sum(
        !flag_missing_coords & !flag_missing_species &
        !flag_invalid_date & !flag_any_coord_issue
      ),
      minimal_clean = sum(
        !flag_missing_coords & !flag_missing_species
      )
    )

  message("\nRecords passing different filter levels:")
  message("  Strict (all flags pass): ", format(clean_counts$strict_clean, big.mark = ","),
          " (", round(100 * clean_counts$strict_clean / n_original, 1), "%)")
  message("  Moderate (core flags): ", format(clean_counts$moderate_clean, big.mark = ","),
          " (", round(100 * clean_counts$moderate_clean / n_original, 1), "%)")
  message("  Minimal (coords + species): ", format(clean_counts$minimal_clean, big.mark = ","),
          " (", round(100 * clean_counts$minimal_clean / n_original, 1), "%)")

  # Calculate additional quality metrics ---------------------------------------
  quality_metrics <- list(
    # Overall stats
    n_records = n_original,
    n_strict_clean = clean_counts$strict_clean,
    n_moderate_clean = clean_counts$moderate_clean,
    n_minimal_clean = clean_counts$minimal_clean,

    # Data completeness
    completeness = list(
      coordinates = sum(!kenya_flagged$flag_missing_coords),
      species_id = tax_completeness$pct_n_with_species,
      genus_id = tax_completeness$pct_n_with_genus,
      family_id = tax_completeness$pct_n_with_family,
      dates = 100 * sum(!kenya_flagged$flag_invalid_date, na.rm = TRUE) / n_original
    ),

    # Coordinate quality
    coordinate_quality = list(
      uncertainty_summary = uncertainty_summary,
      coord_issues_summary = coord_issues_df,
      total_coord_flags = sum(kenya_flagged$flag_any_coord_issue),
      percent_coord_flags = 100 * sum(kenya_flagged$flag_any_coord_issue) / n_original
    ),

    # Taxonomic completeness
    taxonomic_completeness = tax_completeness,

    # Temporal quality
    temporal_quality = date_issues,

    # Duplicates
    duplicate_info = duplicates_summary,

    # Basis of record
    basis_of_record = basis_summary
  )

  # Compile final results ------------------------------------------------------
  quality_results <- list(
    tracking = quality_tracking,
    coord_issues = coord_issues_df,
    metrics = quality_metrics,
    n_original = n_original,
    clean_counts = clean_counts,
    assessment_date = Sys.Date(),
    flag_columns = names(kenya_flagged)[grepl("^flag_", names(kenya_flagged))]
  )

  # Save results ---------------------------------------------------------------
  saveRDS(quality_results, file.path(results_dir, "data_quality_assessment.rds"))
  write_csv(quality_tracking, file.path(results_dir, "data_quality_tracking.csv"))
  write_csv(coord_issues_df, file.path(results_dir, "coordinate_issues_summary.csv"))
  saveRDS(kenya_flagged, file.path(results_dir, "kenya_gbif_flagged.rds"))

  # Also save to old location for backward compatibility
  saveRDS(quality_results, file.path(data_outputs, "data_quality_assessment.rds"))
  write_csv(quality_tracking, file.path(data_outputs, "data_quality_tracking.csv"))
  write_csv(coord_issues_df, file.path(data_outputs, "coordinate_issues_summary.csv"))
  saveRDS(kenya_flagged, file.path(data_outputs, "kenya_gbif_flagged.rds"))

  # Also save to processed for easy access by other scripts
  saveRDS(kenya_flagged, file.path(data_processed, "kenya_gbif_flagged.rds"))

  message("\n=== Data quality assessment complete ===")
  message("Results saved to: ", results_dir)
  message("Flagged dataset saved to: ", file.path(data_processed, "kenya_gbif_flagged.rds"))
}

# Load results for visualization -----------------------------------------------
if (!exists("quality_results")) {
  message("Loading saved data quality results...")
  quality_results <- readRDS(file.path(results_dir, "data_quality_assessment.rds"))
}

# Create visualizations --------------------------------------------------------
message("\nGenerating data quality visualizations...")

# 1. Flag frequency plot
flag_plot <- quality_results$tracking %>%
  filter(records_flagged > 0) %>%
  mutate(issue = fct_reorder(factor(issue), records_flagged)) %>%
  ggplot(aes(x = issue, y = records_flagged, fill = percent_flagged)) +
  geom_col() +
  geom_text(aes(label = paste0(format(records_flagged, big.mark = ","), "\n",
                                "(", round(percent_flagged, 1), "%)")),
            hjust = -0.1, size = 3) +
  scale_fill_viridis_c(option = "magma", direction = -1,
                       name = "% of records") +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.15))) +
  coord_flip() +
  labs(
    title = "Data Quality Flags by Issue Type",
    subtitle = paste0("Total records: ", format(quality_results$n_original, big.mark = ",")),
    x = "Quality Issue",
    y = "Records Flagged"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "right"
  )

ggsave(file.path(figures_dir, "data_quality_flags.png"),
       flag_plot, width = 12, height = 10, dpi = 300)

# 2. Filter level comparison
filter_comparison <- tibble(
  filter_level = c("Original\n(no filtering)",
                   "Minimal\n(coords + species)",
                   "Moderate\n(+ dates + coord quality)",
                   "Strict\n(all flags pass)"),
  n_records = c(quality_results$n_original,
                quality_results$clean_counts$minimal_clean,
                quality_results$clean_counts$moderate_clean,
                quality_results$clean_counts$strict_clean),
  percent = 100 * n_records / quality_results$n_original
) %>%
  mutate(filter_level = factor(filter_level, levels = filter_level))

filter_plot <- ggplot(filter_comparison, aes(x = filter_level, y = n_records, fill = percent)) +
  geom_col(alpha = 0.8) +
  geom_text(aes(label = paste0(format(n_records, big.mark = ","), "\n",
                                "(", round(percent, 1), "%)")),
            vjust = -0.5, size = 4) +
  scale_fill_viridis_c(option = "plasma", direction = -1,
                       name = "% of original") +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "Records Retained Under Different Filtering Strategies",
    subtitle = "Downstream analyses can choose appropriate filter level",
    x = "Filter Level",
    y = "Records Retained"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "right",
    axis.text.x = element_text(hjust = 0.5)
  )

ggsave(file.path(figures_dir, "data_quality_filter_levels.png"),
       filter_plot, width = 12, height = 8, dpi = 300)

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

message("Visualizations saved to: ", figures_dir)
message("\n=== All data quality assessments complete ===")
message("\nNOTE: This script FLAGS quality issues but does NOT filter data.")
message("Downstream analyses should load 'kenya_gbif_flagged.rds' and apply")
message("filters appropriate to their specific requirements.")
