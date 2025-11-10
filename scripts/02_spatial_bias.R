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

# Prepare data for occAssess (data.frame format required)
kenya_df <- kenya_data %>%
  select(species, decimalLongitude, decimalLatitude, eventDate, year) %>%
  as.data.frame()

# Record number assessment
message("Assessing record numbers...")
record_assessment <- assessRecordNumber(
  dat = kenya_df,
  xCol = "decimalLongitude",
  yCol = "decimalLatitude",
  gridRes = 10000,  # 10km
  logCount = TRUE
)

# Species number assessment
message("Assessing species numbers...")
species_assessment <- assessSpeciesNumber(
  dat = kenya_df,
  xCol = "decimalLongitude",
  yCol = "decimalLatitude",
  gridRes = 10000,
  logCount = TRUE,
  speciesCol = "species"
)

# Species coverage assessment
message("Assessing species coverage...")
coverage_assessment <- assessSpeciesCoverage(
  dat = kenya_df,
  xCol = "decimalLongitude",
  yCol = "decimalLatitude",
  speciesCol = "species",
  minCoords = 5  # Minimum 5 occurrences per species
)

# Save occAssess results
saveRDS(record_assessment, file.path(data_outputs, "occassess_records.rds"))
saveRDS(species_assessment, file.path(data_outputs, "occassess_species.rds"))
saveRDS(coverage_assessment, file.path(data_outputs, "occassess_coverage.rds"))

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

# Plot 4: Combined occAssess visualizations
png(file.path(figures_dir, "04_occassess_assessments.png"),
    width = 3000, height = 2000, res = 300)
par(mfrow = c(2, 2))
plot(record_assessment, main = "Record Number Assessment")
plot(species_assessment, main = "Species Number Assessment")
dev.off()

# 7. SUMMARY REPORT ------------------------------------------------------------
spatial_summary <- list(
  effort_summary = effort_summary,
  moran_results = moran_results,
  env_bias_summary = env_bias_summary,
  coverage_stats = summary(coverage_assessment)
)

saveRDS(spatial_summary, file.path(data_outputs, "spatial_bias_summary.rds"))

message("\n=== Spatial bias assessment complete ===")
message("Results saved to: ", data_outputs)
message("Figures saved to: ", figures_dir)

