# ==============================================================================
# Script: 01_data_download.R
# Purpose: Download and clean GBIF occurrence data for Kenya
# Author: Automated Research Pipeline
# Date: 2025-11-10
# ==============================================================================

# Load required packages -------------------------------------------------------
suppressPackageStartupMessages({
  library(rgbif)
  library(tidyverse)
  library(data.table)
  library(sf)
  library(here)
  library(CoordinateCleaner)
  library(lubridate)
})

# Setup paths ------------------------------------------------------------------
data_raw <- here("data", "raw")
data_processed <- here("data", "processed")
dir.create(data_raw, showWarnings = FALSE, recursive = TRUE)
dir.create(data_processed, showWarnings = FALSE, recursive = TRUE)

# Check for existing data ------------------------------------------------------
raw_data_file <- file.path(data_raw, "kenya_gbif_raw.rds")
processed_data_file <- file.path(data_processed, "kenya_gbif_clean.rds")

# Set force_download to TRUE to force re-download even if data exists
force_download <- FALSE

# Flag to track whether we need to clean data
skip_cleaning <- FALSE

if (file.exists(processed_data_file) && !force_download) {
  message("=== PROCESSED DATA ALREADY EXISTS ===")
  message("Found existing cleaned data at: ", processed_data_file)
  message("Loading existing data and skipping download/cleaning steps.")
  message("To force re-download, set force_download <- TRUE")

  # Load processed data and skip to summary generation
  kenya_final <- readRDS(processed_data_file)
  message("Loaded ", nrow(kenya_final), " cleaned records from existing file")
  skip_cleaning <- TRUE

} else if (file.exists(raw_data_file) && !force_download) {
  message("=== RAW DATA ALREADY EXISTS ===")
  message("Found existing raw data at: ", raw_data_file)
  message("Loading existing raw data instead of downloading...")
  kenya_raw <- readRDS(raw_data_file)
  message("Loaded ", nrow(kenya_raw), " records from existing file")
  message("To force re-download, set force_download <- TRUE")

} else {
  # Set GBIF credentials (set as environment variables) -----------------------
  # Use: usethis::edit_r_environ() to set GBIF_USER, GBIF_PWD, GBIF_EMAIL

  # Define download parameters -------------------------------------------------
  message("=== INITIATING NEW GBIF DOWNLOAD ===")
  message("No existing data found. Downloading from GBIF for Kenya...")

  # Option 1: Use existing download key if available
  # Uncomment and replace with actual download key if data already downloaded
  # download_key <- "YOUR_DOWNLOAD_KEY_HERE"
  # kenya_zip <- occ_download_get(download_key, path = data_raw, overwrite = TRUE)
  # kenya_raw <- occ_download_import(x = kenya_zip)

  # Option 2: Request new download (requires GBIF login)
  # Note: This creates a download request. Check status at www.gbif.org
  kenya_download <- occ_download(
    pred("country", "KE"),
    pred("hasCoordinate", TRUE),
    pred("hasGeospatialIssue", FALSE),
    pred_gte("year", 1950),
    pred("occurrenceStatus", "PRESENT"),
    format = "SIMPLE_CSV",
    user = Sys.getenv("GBIF_USER"),
    pwd = Sys.getenv("GBIF_PWD"),
    email = Sys.getenv("GBIF_EMAIL")
  )

  message("Download key: ", kenya_download)
  message("Check download status at: https://www.gbif.org/occurrence/download/", kenya_download)
  options(timeout = 1000000)

  # Wait for download to complete
  message("Waiting for download to complete...")
  occ_download_wait(kenya_download, status_ping = 10, curlopts = list(), quiet = FALSE)

  # Download the data
  message("Downloading data...")
  kenya_zip <- occ_download_get(kenya_download, path = data_raw, overwrite = TRUE)

  # Import data
  message("Importing data...")
  kenya_raw <- occ_download_import(x = kenya_zip)

  # Save raw data
  saveRDS(kenya_raw, file.path(data_raw, "kenya_gbif_raw.rds"))
  message("Raw data saved: ", nrow(kenya_raw), " records")
}

# Data cleaning ----------------------------------------------------------------
if (!skip_cleaning) {
  message("\n=== Data Cleaning ===")
  message("Initial dataset: ", nrow(kenya_raw), " records")

# Convert to data.table for memory-efficient operations
library(data.table)
kenya_dt <- as.data.table(kenya_raw)

# Clear original to free memory
rm(kenya_raw)
gc()

# Filter 1: Remove records with missing coordinates (most restrictive first)
message("Filtering missing coordinates...")
kenya_dt <- kenya_dt[!is.na(decimalLongitude) & !is.na(decimalLatitude)]
message("  After coordinate filter: ", nrow(kenya_dt), " records")
gc()

# Filter 2: Remove records without species-level identification
message("Filtering species identification...")
kenya_dt <- kenya_dt[!is.na(species) & species != ""]
message("  After species filter: ", nrow(kenya_dt), " records")
gc()

# Filter 3: Remove fossil and living specimens
message("Filtering basis of record...")
kenya_dt <- kenya_dt[!basisOfRecord %in% c("FOSSIL_SPECIMEN", "LIVING_SPECIMEN")]
message("  After basis filter: ", nrow(kenya_dt), " records")
gc()

# Filter 4: Remove records with high coordinate uncertainty (>10km)
message("Filtering coordinate uncertainty...")
kenya_dt <- kenya_dt[coordinateUncertaintyInMeters < 10000 | is.na(coordinateUncertaintyInMeters)]
message("  After uncertainty filter: ", nrow(kenya_dt), " records")
gc()

# Handle date columns efficiently
message("Processing date information...")
current_year <- year(Sys.Date())

# Check if year/month/day columns already exist in GBIF data
if ("year" %in% names(kenya_dt) && !"eventDate_parsed" %in% names(kenya_dt)) {
  # Use existing year column and validate
  message("  Using existing year column from GBIF data")
  kenya_dt[, year := as.integer(year)]

  # Only parse eventDate for records with valid years or missing years
  kenya_dt[is.na(year) & !is.na(eventDate),
           eventDate_parsed := suppressWarnings(ymd(eventDate, quiet = TRUE))]
  kenya_dt[is.na(year) & !is.na(eventDate_parsed),
           year := year(eventDate_parsed)]

  # Parse month and day if they don't exist
  if (!"month" %in% names(kenya_dt)) {
    kenya_dt[, month := as.integer(month)]
  }
  if (!"day" %in% names(kenya_dt)) {
    kenya_dt[, day := as.integer(day)]
  }

} else {
  # Parse dates from eventDate column
  message("  Parsing eventDate column (this may take a moment)...")

  # Only attempt to parse non-empty, non-NA dates
  kenya_dt[!is.na(eventDate) & eventDate != "",
           eventDate_parsed := suppressWarnings(ymd(eventDate, quiet = TRUE))]

  # Extract year, month, day from parsed dates
  kenya_dt[!is.na(eventDate_parsed), `:=`(
    year = year(eventDate_parsed),
    month = month(eventDate_parsed),
    day = day(eventDate_parsed)
  )]
}

# Count parsing failures
n_failed <- kenya_dt[is.na(year) | year < 1950 | year > current_year, .N]
message("  Records with invalid/missing years: ", n_failed)

# Filter 5: Remove records with invalid years
message("Filtering invalid years...")
kenya_dt <- kenya_dt[!is.na(year) & year >= 1950 & year <= current_year]
message("  After year filter: ", nrow(kenya_dt), " records")
gc()

# Filter 6: Remove duplicate records
message("Removing duplicates...")
# Create a date column for duplicate detection (handle NA dates)
kenya_dt[, event_date := as.Date(ifelse(!is.na(eventDate_parsed),
                                         as.character(eventDate_parsed),
                                         NA_character_))]

# Remove duplicates based on species, location, and date
initial_n <- nrow(kenya_dt)
kenya_dt <- unique(kenya_dt, by = c("species", "decimalLongitude", "decimalLatitude", "event_date"))
message("  Removed ", initial_n - nrow(kenya_dt), " duplicate records")
message("  After deduplication: ", nrow(kenya_dt), " records")
gc()

# Convert back to tibble for compatibility with downstream code
kenya_clean <- as_tibble(kenya_dt)
rm(kenya_dt)
gc()

message("After initial cleaning: ", nrow(kenya_clean), " records")

# Coordinate cleaning using CoordinateCleaner ----------------------------------
message("Applying coordinate cleaning filters...")

# Convert to spatial object
kenya_sf <- kenya_clean %>%
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)

# Flag potentially problematic coordinates
flags <- clean_coordinates(
  x = kenya_clean,
  lon = "decimalLongitude",
  lat = "decimalLatitude",
  species = "species",
  countries = "countryCode",
  tests = c("capitals", "centroids", "equal", "gbif",
            "outliers", "zeros", "urban"),
  capitals_ref = NULL,
  centroids_ref = NULL,
  outliers_method = "quantile",
  outliers_mtp = 5,
  outliers_td = 1000,
  outliers_size = 10,
  range_ref = NULL,
  zeros_rad = 0.5,
  capitals_rad = 10000,
  centroids_rad = 5000,
  verbose = TRUE,
  value = "flagged"
)

# Add flags to dataset
# Note: value = "flagged" returns TRUE for problematic records, so negate for clean records
kenya_clean <- kenya_clean %>%
  mutate(coord_clean = !flags)

# Summary of flagged records
message("Coordinate cleaning results:")
message("  Records flagged: ", sum(flags))
message("  Clean records: ", sum(!flags))

# Keep only clean records
kenya_final <- kenya_clean %>%
  filter(coord_clean == TRUE) %>%
  dplyr::select(-coord_clean)

message("Final dataset: ", nrow(kenya_final), " records")

  # Save cleaned data
  saveRDS(kenya_final, file.path(data_processed, "kenya_gbif_clean.rds"))
  write_csv(kenya_final, file.path(data_processed, "kenya_gbif_clean.csv"))
  message("Clean data saved to: ", file.path(data_processed, "kenya_gbif_clean.rds"))

} # End of cleaning section

# Data summary -----------------------------------------------------------------
message("\n=== Generating Summary Statistics ===")
summary_stats <- kenya_final %>%
  summarise(
    n_records = n(),
    n_species = n_distinct(species),
    n_genera = n_distinct(genus),
    n_families = n_distinct(family),
    n_orders = n_distinct(order),
    n_classes = n_distinct(class),
    n_phyla = n_distinct(phylum),
    year_min = min(year, na.rm = TRUE),
    year_max = max(year, na.rm = TRUE),
    coord_uncertainty_median = median(coordinateUncertaintyInMeters, na.rm = TRUE)
  )

print(summary_stats)

# Taxonomic summary
tax_summary <- kenya_final %>%
  group_by(kingdom, phylum, class) %>%
  summarise(
    n_records = n(),
    n_species = n_distinct(species),
    .groups = "drop"
  ) %>%
  arrange(desc(n_records))

print(head(tax_summary, 20))

# Save summary statistics
saveRDS(summary_stats, file.path(data_processed, "summary_stats.rds"))
write_csv(tax_summary, file.path(data_processed, "taxonomic_summary.csv"))

# Create metadata file ---------------------------------------------------------
# Check if metadata already exists to preserve original download info
metadata_file <- file.path(data_processed, "metadata.rds")

if (file.exists(metadata_file)) {
  message("Loading existing metadata...")
  metadata <- readRDS(metadata_file)
  # Update clean records count in case reprocessing happened
  metadata$clean_records <- nrow(kenya_final)
  metadata$last_processed <- Sys.Date()
} else {
  # Create new metadata
  # kenya_download will exist if we did a new download
  if (exists("kenya_download")) {
    download_key <- kenya_download
    download_doi <- paste0("https://doi.org/10.15468/dl.", kenya_download)
  } else {
    # If loaded from existing raw data, try to get info from file
    download_key <- "Unknown (loaded from existing file)"
    download_doi <- "Unknown (loaded from existing file)"
  }

  # Get original record count from raw data file
  if (file.exists(raw_data_file)) {
    kenya_raw_meta <- readRDS(raw_data_file)
    raw_count <- nrow(kenya_raw_meta)
    rm(kenya_raw_meta)
  } else if (exists("kenya_raw")) {
    raw_count <- nrow(kenya_raw)
  } else {
    raw_count <- NA
  }

  metadata <- list(
    download_key = download_key,
    download_date = Sys.Date(),
    gbif_doi = download_doi,
    raw_records = raw_count,
    clean_records = nrow(kenya_final),
    cleaning_steps = c(
      "Removed records with missing coordinates",
      "Removed records with uncertainty > 10km",
      "Removed fossil and living specimens",
      "Removed records without species-level ID",
      "Removed records outside 1950-present",
      "Removed duplicate records",
      "Applied CoordinateCleaner filters"
    )
  )
}

saveRDS(metadata, metadata_file)

message("\n=== Data download and cleaning complete ===")
message("Clean data saved to: ", file.path(data_processed, "kenya_gbif_clean.rds"))
message("Citation DOI: ", metadata$gbif_doi)
