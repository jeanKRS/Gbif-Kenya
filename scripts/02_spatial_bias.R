# ==============================================================================
# Script: 02_spatial_bias.R
# Purpose: Assess spatial biases in GBIF data for Kenya
# Author: Automated Research Pipeline
# Date: 2025-11-10
# ==============================================================================

# Load required packages -------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
  library(here)
  library(occAssess)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(terra)
  library(geodata)
  library(osmdata)
  library(ape)
  library(spdep)
  library(viridis)
  library(patchwork)
})

# Setup paths ------------------------------------------------------------------
data_processed <- here("data", "processed")
data_outputs <- here("data", "outputs")
figures_dir <- here("figures")
dir.create(data_outputs, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# Load cleaned data ------------------------------------------------------------
message("Loading cleaned GBIF data...")
kenya_data <- readRDS(file.path(data_processed, "kenya_gbif_clean.rds"))

# Get Kenya boundary -----------------------------------------------------------
message("Loading Kenya boundary...")
kenya_sf <- ne_countries(country = "Kenya", returnclass = "sf", scale = "large")
kenya_boundary <- st_union(kenya_sf)

# Create spatial object from occurrence data
kenya_occurrences <- kenya_data %>%
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)

# 1. SAMPLING EFFORT ANALYSIS --------------------------------------------------
message("\n=== Analyzing sampling effort ===")

# Create hexagonal grid for analysis (10km resolution)
kenya_grid <- st_make_grid(
  kenya_boundary,
  cellsize = 0.1,  # ~10km at equator
  square = FALSE,   # hexagonal grid
  what = "polygons"
) %>%
  st_sf() %>%
  mutate(grid_id = row_number()) %>%
  st_intersection(kenya_boundary)

# Calculate records per grid cell
grid_effort <- kenya_grid %>%
  st_join(kenya_occurrences) %>%
  group_by(grid_id) %>%
  summarise(
    n_records = sum(!is.na(species)),
    n_species = n_distinct(species, na.rm = TRUE),
    .groups = "drop"
  )

# Add grid cells with zero records
kenya_grid_complete <- kenya_grid %>%
  left_join(st_drop_geometry(grid_effort), by = "grid_id") %>%
  mutate(
    n_records = replace_na(n_records, 0),
    n_species = replace_na(n_species, 0),
    log_records = log10(n_records + 1)
  )

# Save grid data
saveRDS(kenya_grid_complete, file.path(data_outputs, "spatial_grid_effort.rds"))

# Summary statistics
effort_summary <- kenya_grid_complete %>%
  st_drop_geometry() %>%
  summarise(
    total_cells = n(),
    cells_with_records = sum(n_records > 0),
    percent_covered = 100 * cells_with_records / total_cells,
    mean_records = mean(n_records),
    median_records = median(n_records),
    sd_records = sd(n_records),
    max_records = max(n_records),
    mean_species = mean(n_species),
    median_species = median(n_species)
  )

print(effort_summary)

# 2. SPATIAL AUTOCORRELATION ---------------------------------------------------
message("\n=== Testing spatial autocorrelation ===")

# Create neighbors list (queen contiguity)
coords_centroids <- st_centroid(kenya_grid_complete) %>%
  st_coordinates()

# Create distance-based neighbors (10km threshold)
neighbors <- dnearneigh(coords_centroids, 0, 15000)  # 15km
weights <- nb2listw(neighbors, style = "W", zero.policy = TRUE)

# Calculate Moran's I for sampling effort
moran_test <- moran.test(kenya_grid_complete$n_records, weights, zero.policy = TRUE)

moran_results <- data.frame(
  variable = "sampling_effort",
  morans_i = moran_test$estimate["Moran I statistic"],
  expectation = moran_test$estimate["Expectation"],
  variance = moran_test$estimate["Variance"],
  p_value = moran_test$p.value,
  interpretation = ifelse(moran_test$statistic > 0, "Clustered", "Dispersed")
)

print(moran_results)
saveRDS(moran_results, file.path(data_outputs, "spatial_autocorrelation.rds"))

# 3. ACCESSIBILITY BIAS ASSESSMENT ---------------------------------------------
message("\n=== Assessing accessibility bias ===")

# Download environmental data for Kenya
message("Downloading elevation data...")
elevation_data <- elevation_30s(country = "KEN", path = tempdir())
elevation <- if (is.character(elevation_data)) terra::rast(elevation_data) else elevation_data
elevation = terra::extract(elevation, kenya_coords)[, 1]

# Get climate data
message("Downloading climate data...")
climate_data <- worldclim_country(country = "KEN", var = "bio", path = tempdir())
climate <- if (is.character(climate_data)) terra::rast(climate_data) else climate_data

# Extract environmental values for occurrences
kenya_coords <- kenya_data %>%
  select(decimalLongitude, decimalLatitude) %>%
  as.matrix()

# Extract elevation
kenya_data_env <- kenya_data %>%
  mutate(
    elevation = elevation,
    bio1_temp = terra::extract(climate[[1]], kenya_coords)[, 1],  # Annual mean temp
    bio12_precip = terra::extract(climate[[12]], kenya_coords)[, 1]  # Annual precip
  )

# Calculate distance to nearest record for each grid cell
grid_centroids <- st_centroid(kenya_grid_complete)
occurrence_coords <- st_coordinates(kenya_occurrences)

dist_to_nearest <- apply(
  st_coordinates(grid_centroids),
  1,
  function(x) {
    dists <- sqrt((occurrence_coords[, 1] - x[1])^2 +
                  (occurrence_coords[, 2] - x[2])^2)
    min(dists)
  }
)

kenya_grid_complete <- kenya_grid_complete %>%
  mutate(dist_to_nearest_km = dist_to_nearest * 111)  # Convert degrees to km

# 4. OCCASSESS ANALYSES --------------------------------------------------------
message("\n=== Running occAssess assessments ===")

# Prepare data for occAssess with taxonomic levels (data.frame format required)
kenya_df <- kenya_data %>%
  select(species, genus, family, order, class, phylum, kingdom,
         decimalLongitude, decimalLatitude, coordinateUncertaintyInMeters,
         eventDate, year) %>%
  as.data.frame()

# Create 10-year period breaks
message("Creating 10-year period breaks...")
year_range <- range(kenya_df$year, na.rm = TRUE)
period_breaks <- seq(
  floor(year_range[1] / 10) * 10,  # Round down to nearest decade
  ceiling(year_range[2] / 10) * 10,  # Round up to nearest decade
  by = 10
)

message(sprintf("Period breaks: %s to %s",
                min(period_breaks), max(period_breaks)))
message(sprintf("Number of periods: %d", length(period_breaks) - 1))

# Add period labels to data
kenya_df <- kenya_df %>%
  mutate(
    period = cut(year,
                breaks = period_breaks,
                include.lowest = TRUE,
                right = FALSE,
                labels = paste0(head(period_breaks, -1), "-",
                              tail(period_breaks, -1) - 1))
  )

# Summary of records by period
period_summary <- kenya_df %>%
  group_by(period) %>%
  summarise(
    n_records = n(),
    n_species = n_distinct(species),
    .groups = "drop"
  )

print(period_summary)
saveRDS(period_summary, file.path(data_outputs, "records_by_period.rds"))

# Define taxonomic levels and other grouping features for analysis
identifiers <- list(
  taxonomic = c("species", "genus", "family", "order", "class", "phylum", "kingdom"),
  grouping = c("period")  # Can add more grouping features here
)

# Record number assessment - Run for each taxonomic level
message("\n--- Assessing record numbers by taxonomic level ---")
record_assessments <- list()

for (tax_level in identifiers$taxonomic) {
  message(sprintf("  Running assessRecordNumber for identifier: %s", tax_level))

  # Skip if column has all NA values
  if (all(is.na(kenya_df[[tax_level]]))) {
    message(sprintf("  Skipping %s - all values are NA", tax_level))
    next
  }

  tryCatch({
    record_assessments[[tax_level]] <- assessRecordNumber(
      dat = kenya_df,
      periods = period_breaks,
      species = "species",
      x = "decimalLongitude",
      y = "decimalLatitude",
      year = "year",
      spatialUncertainty = "coordinateUncertaintyInMeters",
      identifier = tax_level
    )
    message(sprintf("  ✓ Completed assessment for %s", tax_level))
  }, error = function(e) {
    message(sprintf("  ✗ Error in assessment for %s: %s", tax_level, e$message))
  })
}

# Also run assessment grouped by period
message("  Running assessRecordNumber for identifier: period")
tryCatch({
  record_assessments[["period"]] <- assessRecordNumber(
    dat = kenya_df %>% filter(!is.na(period)),
    periods = period_breaks,
    species = "species",
    x = "decimalLongitude",
    y = "decimalLatitude",
    year = "year",
    spatialUncertainty = "coordinateUncertaintyInMeters",
    identifier = "period"
  )
  message("  ✓ Completed assessment for period")
}, error = function(e) {
  message(sprintf("  ✗ Error in assessment for period: %s", e$message))
})

# Species number assessment - Run for different taxonomic levels
message("\n--- Assessing species numbers by taxonomic level ---")
species_assessments <- list()

for (tax_level in identifiers$taxonomic) {
  message(sprintf("  Running assessSpeciesNumber for %s", tax_level))

  if (all(is.na(kenya_df[[tax_level]]))) {
    message(sprintf("  Skipping %s - all values are NA", tax_level))
    next
  }

  tryCatch({
    species_assessments[[tax_level]] <- assessSpeciesNumber(
      dat = kenya_df,
      xCol = "decimalLongitude",
      yCol = "decimalLatitude",
      gridRes = 10000,
      logCount = TRUE,
      speciesCol = tax_level
    )
    message(sprintf("  ✓ Completed assessment for %s", tax_level))
  }, error = function(e) {
    message(sprintf("  ✗ Error in assessment for %s: %s", tax_level, e$message))
  })
}

# Species coverage assessment - Run for different taxonomic levels
message("\n--- Assessing species coverage by taxonomic level ---")
coverage_assessments <- list()

for (tax_level in identifiers$taxonomic) {
  message(sprintf("  Running assessSpeciesCoverage for %s", tax_level))

  if (all(is.na(kenya_df[[tax_level]]))) {
    message(sprintf("  Skipping %s - all values are NA", tax_level))
    next
  }

  tryCatch({
    coverage_assessments[[tax_level]] <- assessSpeciesCoverage(
      dat = kenya_df,
      xCol = "decimalLongitude",
      yCol = "decimalLatitude",
      speciesCol = tax_level,
      minCoords = 5  # Minimum 5 occurrences per taxonomic unit
    )
    message(sprintf("  ✓ Completed assessment for %s", tax_level))
  }, error = function(e) {
    message(sprintf("  ✗ Error in assessment for %s: %s", tax_level, e$message))
  })
}

# Save occAssess results
message("\n--- Saving occAssess results ---")
saveRDS(record_assessments, file.path(data_outputs, "occassess_records_by_taxlevel.rds"))
saveRDS(species_assessments, file.path(data_outputs, "occassess_species_by_taxlevel.rds"))
saveRDS(coverage_assessments, file.path(data_outputs, "occassess_coverage_by_taxlevel.rds"))

# Also save the main species-level assessment for backward compatibility
if (!is.null(record_assessments[["species"]])) {
  saveRDS(record_assessments[["species"]], file.path(data_outputs, "occassess_records.rds"))
}
if (!is.null(species_assessments[["species"]])) {
  saveRDS(species_assessments[["species"]], file.path(data_outputs, "occassess_species.rds"))
}
if (!is.null(coverage_assessments[["species"]])) {
  saveRDS(coverage_assessments[["species"]], file.path(data_outputs, "occassess_coverage.rds"))
}

# Create summary of assessments completed
assessment_summary <- data.frame(
  taxonomic_level = identifiers$taxonomic,
  record_assessment = sapply(identifiers$taxonomic, function(x) !is.null(record_assessments[[x]])),
  species_assessment = sapply(identifiers$taxonomic, function(x) !is.null(species_assessments[[x]])),
  coverage_assessment = sapply(identifiers$taxonomic, function(x) !is.null(coverage_assessments[[x]]))
)

print(assessment_summary)
saveRDS(assessment_summary, file.path(data_outputs, "assessment_summary.rds"))

# 5. ENVIRONMENTAL BIAS --------------------------------------------------------
message("\n=== Assessing environmental bias ===")

# Sample available environmental space
set.seed(123)
background_points <- st_sample(kenya_boundary, size = 10000) %>%
  st_coordinates()

# Extract environmental values for background
bg_env <- data.frame(
  type = "available",
  elevation = terra::extract(elevation, background_points)[, 2],
  temperature = terra::extract(climate[[1]], background_points)[, 2],
  precipitation = terra::extract(climate[[12]], background_points)[, 2]
)

# Environmental values for occurrences
occ_env <- kenya_data_env %>%
  select(elevation, temperature = bio1_temp, precipitation = bio12_precip) %>%
  filter(!is.na(elevation)) %>%
  mutate(type = "sampled")

# Combine for comparison
env_comparison <- bind_rows(
  occ_env %>% select(type, elevation, temperature, precipitation),
  bg_env %>% filter(!is.na(elevation))
)

# Statistical tests for environmental bias
env_tests <- list(
  elevation = ks.test(
    occ_env$elevation,
    bg_env$elevation[!is.na(bg_env$elevation)]
  ),
  temperature = ks.test(
    occ_env$temperature[!is.na(occ_env$temperature)],
    bg_env$temperature[!is.na(bg_env$temperature)]
  ),
  precipitation = ks.test(
    occ_env$precipitation[!is.na(occ_env$precipitation)],
    bg_env$precipitation[!is.na(bg_env$precipitation)]
  )
)

env_bias_summary <- map_df(names(env_tests), ~{
  tibble(
    variable = .x,
    D_statistic = env_tests[[.x]]$statistic,
    p_value = env_tests[[.x]]$p.value,
    significant = p_value < 0.05
  )
})

print(env_bias_summary)
saveRDS(env_comparison, file.path(data_outputs, "environmental_comparison.rds"))
saveRDS(env_bias_summary, file.path(data_outputs, "environmental_bias_tests.rds"))

# 6. VISUALIZATION -------------------------------------------------------------
message("\n=== Creating visualizations ===")

# Plot 1: Sampling effort map
p1 <- ggplot() +
  geom_sf(data = kenya_grid_complete, aes(fill = log_records), color = NA) +
  geom_sf(data = kenya_boundary, fill = NA, color = "black", linewidth = 0.5) +
  scale_fill_viridis_c(
    name = "log10(Records + 1)",
    option = "plasma",
    na.value = "grey90"
  ) +
  labs(title = "Sampling Effort Across Kenya",
       subtitle = paste0("Total records: ", nrow(kenya_data))) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text = element_text(size = 8)
  )

ggsave(file.path(figures_dir, "01_sampling_effort_map.png"),
       p1, width = 10, height = 8, dpi = 300)

# Plot 2: Species richness map
p2 <- ggplot() +
  geom_sf(data = kenya_grid_complete, aes(fill = n_species), color = NA) +
  geom_sf(data = kenya_boundary, fill = NA, color = "black", linewidth = 0.5) +
  scale_fill_viridis_c(
    name = "Number of Species",
    option = "viridis",
    na.value = "grey90",
    trans = "log10"
  ) +
  labs(title = "Species Richness Across Kenya",
       subtitle = paste0("Total species: ", n_distinct(kenya_data$species))) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text = element_text(size = 8)
  )

ggsave(file.path(figures_dir, "02_species_richness_map.png"),
       p2, width = 10, height = 8, dpi = 300)

# Plot 3: Environmental space comparison
env_long <- env_comparison %>%
  pivot_longer(cols = c(elevation, temperature, precipitation),
               names_to = "variable", values_to = "value") %>%
  filter(!is.na(value))

p3 <- ggplot(env_long, aes(x = value, fill = type)) +
  geom_density(alpha = 0.6) +
  facet_wrap(~variable, scales = "free", ncol = 1) +
  scale_fill_manual(values = c("available" = "grey50", "sampled" = "red"),
                    name = "Environmental Space") +
  labs(title = "Environmental Bias Assessment",
       subtitle = "Comparison of sampled vs. available environmental space",
       x = "Environmental Value", y = "Density") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(file.path(figures_dir, "03_environmental_bias.png"),
       p3, width = 10, height = 10, dpi = 300)

# Plot 4: Combined occAssess visualizations for species level
if (!is.null(record_assessments[["species"]]) && !is.null(species_assessments[["species"]])) {
  png(file.path(figures_dir, "04_occassess_assessments_species.png"),
      width = 3000, height = 2000, res = 300)
  par(mfrow = c(2, 2))
  plot(record_assessments[["species"]], main = "Record Number Assessment (Species)")
  plot(species_assessments[["species"]], main = "Species Number Assessment")
  dev.off()
}

# Plot 5-11: occAssess visualizations for each taxonomic level
message("Creating visualizations for each taxonomic level...")
for (tax_level in identifiers$taxonomic) {
  if (!is.null(record_assessments[[tax_level]])) {
    tryCatch({
      png(file.path(figures_dir, sprintf("04_occassess_%s.png", tax_level)),
          width = 2400, height = 1600, res = 300)
      par(mfrow = c(1, 1))
      plot(record_assessments[[tax_level]],
           main = sprintf("Record Number Assessment (%s)", tax_level))
      dev.off()
      message(sprintf("  ✓ Created visualization for %s", tax_level))
    }, error = function(e) {
      message(sprintf("  ✗ Error creating visualization for %s: %s", tax_level, e$message))
    })
  }
}

# Create comparative summary plot of record counts by taxonomic level and period
message("Creating comparative summary plots...")
record_summary_data <- data.frame()

for (tax_level in identifiers$taxonomic) {
  if (!is.null(record_assessments[[tax_level]])) {
    tryCatch({
      # Extract summary data from assessment object
      summary_df <- data.frame(
        taxonomic_level = tax_level,
        n_taxa = length(unique(kenya_df[[tax_level]])),
        stringsAsFactors = FALSE
      )
      record_summary_data <- bind_rows(record_summary_data, summary_df)
    }, error = function(e) {
      message(sprintf("  Warning: Could not extract summary for %s", tax_level))
    })
  }
}

if (nrow(record_summary_data) > 0) {
  p_tax_summary <- ggplot(record_summary_data,
                          aes(x = reorder(taxonomic_level, n_taxa), y = n_taxa)) +
    geom_col(fill = "steelblue") +
    geom_text(aes(label = n_taxa), hjust = -0.2, size = 3.5) +
    coord_flip() +
    scale_y_log10(labels = scales::comma) +
    labs(title = "Number of Taxa by Taxonomic Level",
         x = "Taxonomic Level",
         y = "Number of Taxa (log scale)") +
    theme_minimal() +
    theme(panel.grid.minor = element_blank())

  ggsave(file.path(figures_dir, "05_taxonomic_summary.png"),
         p_tax_summary, width = 10, height = 6, dpi = 300)
}

# 7. SUMMARY REPORT ------------------------------------------------------------
spatial_summary <- list(
  effort_summary = effort_summary,
  moran_results = moran_results,
  env_bias_summary = env_bias_summary,
  period_summary = period_summary,
  period_breaks = period_breaks,
  assessment_summary = assessment_summary,
  taxonomic_levels_analyzed = identifiers$taxonomic,
  coverage_stats = if (!is.null(coverage_assessments[["species"]])) {
    summary(coverage_assessments[["species"]])
  } else {
    "No species-level coverage assessment available"
  }
)

saveRDS(spatial_summary, file.path(data_outputs, "spatial_bias_summary.rds"))

message("\n=== Spatial bias assessment complete ===")
message("Results saved to: ", data_outputs)
message("Figures saved to: ", figures_dir)

