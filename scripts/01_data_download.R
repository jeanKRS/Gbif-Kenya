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

# Set GBIF credentials (set as environment variables) -------------------------
# Use: usethis::edit_r_environ() to set GBIF_USER, GBIF_PWD, GBIF_EMAIL

# Define download parameters ---------------------------------------------------
message("Initiating GBIF download for Kenya...")

# Option 1: Use existing download key if available
# Uncomment and replace with actual download key if data already downloaded
# download_key <- 0024317-251025141854904 # "YOUR_DOWNLOAD_KEY_HERE"

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

# Wait for download to complete
message("Waiting for download to complete...")
occ_download_wait(kenya_download, status_ping = 10, curlopts = list(), quiet = FALSE)

# Download the data
message("Downloading data...")
options(timeout = 1000000)
kenya_zip <- occ_download_get(kenya_download, path = data_raw, overwrite = TRUE)

# Import data
message("Importing data...")
kenya_raw <- occ_download_import(kenya_zip)

# Save raw data
saveRDS(kenya_raw, file.path(data_raw, "kenya_gbif_raw.rds"))
message("Raw data saved: ", nrow(kenya_raw), " records")

# Data cleaning ----------------------------------------------------------------
message("Cleaning data...")

kenya_clean <- kenya_raw %>%
  # Remove records with missing coordinates
  filter(!is.na(decimalLongitude), !is.na(decimalLatitude)) %>%
  # Remove records with high coordinate uncertainty (>10km)
  filter(coordinateUncertaintyInMeters < 10000 | is.na(coordinateUncertaintyInMeters)) %>%
  # Remove fossil and living specimens
  filter(!basisOfRecord %in% c("FOSSIL_SPECIMEN", "LIVING_SPECIMEN")) %>%
  # Remove records without species-level identification
  filter(!is.na(species), species != "") %>%
  # Parse dates
  mutate(
    eventDate = ymd(eventDate),
    year = year(eventDate),
    month = month(eventDate),
    day = day(eventDate)
  ) %>%
  # Remove records with invalid years
  filter(year >= 1950, year <= year(Sys.Date())) %>%
  # Remove duplicate records (same species, location, date)
  distinct(species, decimalLongitude, decimalLatitude, eventDate, .keep_all = TRUE)

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
            "outliers", "zeros", "seas", "urban"),
  seas_ref = buffland,
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
kenya_clean <- kenya_clean %>%
  mutate(coord_clean = flags$.summary)

# Summary of flagged records
message("Coordinate cleaning results:")
message("  Records flagged: ", sum(!flags$.summary))
message("  Clean records: ", sum(flags$.summary))

# Keep only clean records
kenya_final <- kenya_clean %>%
  filter(coord_clean == TRUE) %>%
  select(-coord_clean)

message("Final dataset: ", nrow(kenya_final), " records")

# Data summary -----------------------------------------------------------------
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

# Save cleaned data ------------------------------------------------------------
saveRDS(kenya_final, file.path(data_processed, "kenya_gbif_clean.rds"))
write_csv(kenya_final, file.path(data_processed, "kenya_gbif_clean.csv"))

# Save summary statistics
saveRDS(summary_stats, file.path(data_processed, "summary_stats.rds"))
write_csv(tax_summary, file.path(data_processed, "taxonomic_summary.csv"))

# Create metadata file
metadata <- list(
  download_key = kenya_download,
  download_date = Sys.Date(),
  gbif_doi = paste0("https://doi.org/10.15468/dl.", kenya_download),
  raw_records = nrow(kenya_raw),
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

saveRDS(metadata, file.path(data_processed, "metadata.rds"))

message("\n=== Data download and cleaning complete ===")
message("Clean data saved to: ", file.path(data_processed, "kenya_gbif_clean.rds"))
message("Citation DOI: ", metadata$gbif_doi)
