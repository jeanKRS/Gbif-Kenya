# ==============================================================================
# Script: 05_statistical_models.R
# Purpose: Statistical modeling of sampling bias predictors
# Author: Automated Research Pipeline
# Date: 2025-11-10
# ==============================================================================

# Load required packages -------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
  library(here)
  library(terra)
  library(geodata)
  library(MASS)
  library(lme4)
  library(DHARMa)
  library(performance)
  library(effects)
  library(MuMIn)
  library(mgcv)
})

# Setup paths ------------------------------------------------------------------
data_processed <- here("data", "processed")
data_outputs <- here("data", "outputs")
figures_dir <- here("figures")

# Load data --------------------------------------------------------------------
message("Loading data...")
kenya_data <- readRDS(file.path(data_processed, "kenya_gbif_clean.rds"))
spatial_grid <- readRDS(file.path(data_outputs, "spatial_grid_effort.rds"))

# 1. PREPARE ENVIRONMENTAL COVARIATES ------------------------------------------
message("\n=== Preparing environmental covariates ===")

# Download environmental data
message("Downloading elevation data...")
elevation_data <- elevation_30s(country = "KEN", path = tempdir())
elevation <- if (is.character(elevation_data)) terra::rast(elevation_data) else elevation_data

message("Downloading climate data...")
climate_data <- worldclim_country(country = "KEN", var = "bio", path = tempdir())
climate <- if (is.character(climate_data)) terra::rast(climate_data) else climate_data

# Extract bio1 (temperature) and bio12 (precipitation)
temperature <- climate[[1]]
precipitation <- climate[[12]]

# Get grid centroids for extraction
grid_centroids <- st_centroid(spatial_grid)
grid_coords <- st_coordinates(grid_centroids)

# Extract environmental values
message("Extracting environmental values for grid cells...")
grid_env <- spatial_grid %>%
  mutate(
    elevation = terra::extract(elevation, grid_coords)[, 2],
    temperature = terra::extract(temperature, grid_coords)[, 2],
    precipitation = terra::extract(precipitation, grid_coords)[, 2],
    lon = grid_coords[, 1],
    lat = grid_coords[, 2]
  )

# 2. CALCULATE GEOGRAPHIC COVARIATES -------------------------------------------
message("\n=== Calculating geographic covariates ===")

# Distance to major cities (Nairobi, Mombasa, Kisumu, Nakuru, Eldoret)
major_cities <- data.frame(
  city = c("Nairobi", "Mombasa", "Kisumu", "Nakuru", "Eldoret"),
  lon = c(36.8219, 39.6682, 34.7617, 36.0667, 35.2698),
  lat = c(-1.2921, -4.0435, -0.1022, -0.3031, 0.5143)
) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

# Calculate minimum distance to cities
dist_to_city <- st_distance(grid_centroids, major_cities)
min_dist_to_city <- apply(dist_to_city, 1, min) / 1000  # Convert to km

# Distance from equator (absolute latitude)
dist_from_equator <- abs(grid_coords[, 2])

# Distance from coast (approximate - eastern boundary)
dist_from_coast <- abs(grid_coords[, 1] - 40.5) * 111  # Rough km conversion

# Add to grid data
grid_model <- st_drop_geometry(grid_env) %>%
  mutate(
    dist_to_city_km = as.numeric(min_dist_to_city),
    dist_from_equator = dist_from_equator,
    dist_from_coast_km = dist_from_coast,
    # Scale environmental variables
    elevation_scaled = scale(elevation)[, 1],
    temperature_scaled = scale(temperature)[, 1],
    precipitation_scaled = scale(precipitation)[, 1],
    dist_to_city_scaled = scale(dist_to_city_km)[, 1],
    # Create presence/absence for logistic model
    presence = as.numeric(n_records > 0)
  ) %>%
  filter(!is.na(elevation), !is.na(temperature), !is.na(precipitation))

# 3. GENERALIZED LINEAR MODELS -------------------------------------------------
message("\n=== Fitting GLMs for sampling effort ===")

# Model 1: Poisson GLM for record count
message("Fitting Poisson GLM...")
glm_poisson <- glm(
  n_records ~ elevation_scaled + temperature_scaled + precipitation_scaled +
    dist_to_city_scaled + dist_from_equator + dist_from_coast_km,
  data = grid_model,
  family = poisson(link = "log")
)

# Check for overdispersion
overdispersion <- sum(residuals(glm_poisson, type = "pearson")^2) / glm_poisson$df.residual
message("Overdispersion parameter: ", round(overdispersion, 2))

# Model 2: Negative binomial GLM (if overdispersed)
if (overdispersion > 2) {
  message("Data is overdispersed. Fitting Negative Binomial GLM...")
  glm_nb <- glm.nb(
    n_records ~ elevation_scaled + temperature_scaled + precipitation_scaled +
      dist_to_city_scaled + dist_from_equator + dist_from_coast_km,
    data = grid_model
  )
  primary_model <- glm_nb
  model_type <- "Negative Binomial"
} else {
  primary_model <- glm_poisson
  model_type <- "Poisson"
}

# Model 3: Binomial GLM for presence/absence
message("Fitting Binomial GLM for sampling presence...")
glm_binomial <- glm(
  presence ~ elevation_scaled + temperature_scaled + precipitation_scaled +
    dist_to_city_scaled + dist_from_equator + dist_from_coast_km,
  data = grid_model,
  family = binomial(link = "logit")
)

# 4. MODEL DIAGNOSTICS ---------------------------------------------------------
message("\n=== Running model diagnostics ===")

# DHARMa residuals
sim_resid_count <- simulateResiduals(primary_model, n = 250)
sim_resid_presence <- simulateResiduals(glm_binomial, n = 250)

# Save diagnostic plots
png(file.path(figures_dir, "16_model_diagnostics_count.png"),
    width = 3000, height = 1500, res = 300)
plot(sim_resid_count)
dev.off()

png(file.path(figures_dir, "17_model_diagnostics_presence.png"),
    width = 3000, height = 1500, res = 300)
plot(sim_resid_presence)
dev.off()

# Model performance metrics
perf_count <- performance::model_performance(primary_model)
perf_presence <- performance::model_performance(glm_binomial)

print(perf_count)
print(perf_presence)

# 5. MODEL SUMMARIES AND INFERENCE ---------------------------------------------
message("\n=== Extracting model summaries ===")

# Summary for count model
summary_count <- broom::tidy(primary_model, conf.int = TRUE) %>%
  mutate(
    model = "Count Model",
    significant = p.value < 0.05,
    effect_direction = ifelse(estimate > 0, "Positive", "Negative")
  )

# Summary for presence model
summary_presence <- broom::tidy(glm_binomial, conf.int = TRUE) %>%
  mutate(
    model = "Presence Model",
    significant = p.value < 0.05,
    effect_direction = ifelse(estimate > 0, "Positive", "Negative")
  )

# Combine summaries
model_summaries <- bind_rows(summary_count, summary_presence)
print(model_summaries)

saveRDS(model_summaries, file.path(data_outputs, "model_summaries.rds"))

# 6. MODEL SELECTION -----------------------------------------------------------
message("\n=== Performing model selection ===")

# Fit models with different covariate combinations
models_count <- list(
  full = primary_model,
  geographic = update(primary_model, . ~ dist_to_city_scaled + dist_from_equator + dist_from_coast_km),
  environmental = update(primary_model, . ~ elevation_scaled + temperature_scaled + precipitation_scaled),
  elevation_only = update(primary_model, . ~ elevation_scaled),
  accessibility = update(primary_model, . ~ dist_to_city_scaled)
)

# AIC comparison
aic_comparison <- data.frame(
  model = names(models_count),
  AIC = sapply(models_count, AIC),
  BIC = sapply(models_count, BIC)
) %>%
  arrange(AIC) %>%
  mutate(
    delta_AIC = AIC - min(AIC),
    weight = exp(-0.5 * delta_AIC) / sum(exp(-0.5 * delta_AIC))
  )

print(aic_comparison)
saveRDS(aic_comparison, file.path(data_outputs, "model_selection_aic.rds"))

# 7. GENERALIZED ADDITIVE MODELS -----------------------------------------------
message("\n=== Fitting GAMs for non-linear relationships ===")

# GAM for record count
gam_count <- gam(
  n_records ~ s(elevation, k = 5) + s(temperature, k = 5) +
    s(precipitation, k = 5) + s(dist_to_city_km, k = 5) +
    s(lon, lat, k = 20),
  data = grid_model %>% filter(n_records > 0),
  family = nb(),
  method = "REML"
)

# GAM summary
summary_gam <- summary(gam_count)
print(summary_gam)

# Save GAM
saveRDS(gam_count, file.path(data_outputs, "gam_model.rds"))

# GAM diagnostic plots
png(file.path(figures_dir, "18_gam_diagnostics.png"),
    width = 3000, height = 2000, res = 300)
par(mfrow = c(2, 3))
gam.check(gam_count)
dev.off()

# 8. EFFECT PLOTS --------------------------------------------------------------
message("\n=== Creating effect plots ===")

# Effect plots for GLM
effects_count <- allEffects(primary_model)

png(file.path(figures_dir, "19_glm_effects.png"),
    width = 3600, height = 2400, res = 300)
plot(effects_count, main = paste("Effects on Sampling Effort -", model_type, "GLM"))
dev.off()

# GAM smooth plots
png(file.path(figures_dir, "20_gam_smooths.png"),
    width = 3600, height = 2400, res = 300)
par(mfrow = c(2, 3))
plot(gam_count, shade = TRUE, pages = 1, scale = 0,
     main = "GAM Smooth Effects on Sampling Effort")
dev.off()

# 9. PREDICTIONS AND SPATIAL PATTERNS ------------------------------------------
message("\n=== Generating predictions ===")

# Predict sampling effort
grid_model <- grid_model %>%
  mutate(
    predicted_count = predict(primary_model, newdata = ., type = "response"),
    predicted_presence = predict(glm_binomial, newdata = ., type = "response")
  )

# Calculate residuals
grid_model <- grid_model %>%
  mutate(
    residual_count = n_records - predicted_count,
    residual_std = residual_count / sd(residual_count),
    under_sampled = residual_std < -1,
    over_sampled = residual_std > 1
  )

# Identify under-sampled areas
under_sampled_summary <- grid_model %>%
  filter(under_sampled) %>%
  summarise(
    n_cells = n(),
    mean_elevation = mean(elevation, na.rm = TRUE),
    mean_temp = mean(temperature, na.rm = TRUE),
    mean_precip = mean(precipitation, na.rm = TRUE),
    mean_dist_city = mean(dist_to_city_km, na.rm = TRUE)
  )

print("Under-sampled areas characteristics:")
print(under_sampled_summary)

saveRDS(grid_model, file.path(data_outputs, "grid_with_predictions.rds"))

# 10. VISUALIZATION ------------------------------------------------------------
message("\n=== Creating modeling visualizations ===")

# Rejoin with spatial data for mapping
grid_spatial_pred <- spatial_grid %>%
  left_join(grid_model %>% select(grid_id, predicted_count, residual_std,
                                  under_sampled, over_sampled),
            by = "grid_id")

# Plot: Predicted sampling effort
p1 <- ggplot() +
  geom_sf(data = grid_spatial_pred, aes(fill = log10(predicted_count + 1)), color = NA) +
  scale_fill_viridis_c(option = "plasma", name = "log10(Predicted Records)") +
  labs(title = "Predicted Sampling Effort Based on Environmental & Geographic Covariates",
       subtitle = paste("Model:", model_type, "GLM")) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(file.path(figures_dir, "21_predicted_effort.png"),
       p1, width = 12, height = 8, dpi = 300)

# Plot: Standardized residuals (sampling bias map)
p2 <- ggplot() +
  geom_sf(data = grid_spatial_pred, aes(fill = residual_std), color = NA) +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = 0,
    name = "Standardized Residual",
    limits = c(-3, 3),
    oob = scales::squish
  ) +
  labs(title = "Sampling Bias Map: Observed vs Expected Effort",
       subtitle = "Blue = Under-sampled, Red = Over-sampled") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(file.path(figures_dir, "22_sampling_bias_map.png"),
       p2, width = 12, height = 8, dpi = 300)

# Plot: Coefficient plot
coef_plot_data <- summary_count %>%
  filter(term != "(Intercept)") %>%
  mutate(term = str_remove(term, "_scaled"))

p3 <- ggplot(coef_plot_data, aes(x = reorder(term, estimate), y = estimate)) +
  geom_point(size = 4, color = "steelblue") +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  coord_flip() +
  labs(
    title = "Predictors of Sampling Effort",
    subtitle = paste("Coefficients from", model_type, "GLM with 95% CI"),
    x = "Predictor",
    y = "Coefficient Estimate"
  ) +
  theme_minimal()

ggsave(file.path(figures_dir, "23_coefficient_plot.png"),
       p3, width = 10, height = 6, dpi = 300)

# 11. TAXONOMIC-SPECIFIC MODELS -----------------------------------------------
message("\n=== Fitting taxonomic-specific models ===")

# Create grid data with taxonomic composition
grid_taxonomy <- kenya_data %>%
  filter(!is.na(class)) %>%
  mutate(
    cell_id = paste0(
      floor(decimalLongitude * 10), "_",
      floor(decimalLatitude * 10)
    )
  ) %>%
  group_by(grid_id = cell_id, class) %>%
  summarise(
    n_records_class = n(),
    .groups = "drop"
  )

# Get top 5 classes for modeling
top_5_classes <- kenya_data %>%
  count(class, sort = TRUE) %>%
  head(5) %>%
  pull(class)

message(sprintf("Fitting models for top 5 classes: %s", paste(top_5_classes, collapse = ", ")))

# Join taxonomic data with environmental grid
grid_tax_env <- grid_model %>%
  select(grid_id, elevation_scaled, temperature_scaled, precipitation_scaled,
         dist_to_city_scaled, dist_from_equator, dist_from_coast_km,
         elevation, temperature, precipitation, dist_to_city_km) %>%
  inner_join(
    grid_taxonomy %>% filter(class %in% top_5_classes),
    by = "grid_id"
  )

# Fit separate models for each major taxonomic group
taxonomic_models <- list()
taxonomic_summaries <- list()

for (tax_class in top_5_classes) {
  message(sprintf("  Fitting model for %s...", tax_class))

  class_data <- grid_tax_env %>%
    filter(class == tax_class)

  if (nrow(class_data) < 30) {
    message(sprintf("  Skipping %s - insufficient data (n = %d)", tax_class, nrow(class_data)))
    next
  }

  tryCatch({
    # Fit negative binomial model for this taxonomic group
    model_class <- glm.nb(
      n_records_class ~ elevation_scaled + temperature_scaled + precipitation_scaled +
        dist_to_city_scaled + dist_from_equator + dist_from_coast_km,
      data = class_data
    )

    taxonomic_models[[tax_class]] <- model_class

    # Extract summary
    taxonomic_summaries[[tax_class]] <- broom::tidy(model_class, conf.int = TRUE) %>%
      mutate(
        taxonomic_class = tax_class,
        significant = p.value < 0.05,
        effect_direction = ifelse(estimate > 0, "Positive", "Negative")
      )

    message(sprintf("  ✓ Completed model for %s", tax_class))
  }, error = function(e) {
    message(sprintf("  ✗ Error in model for %s: %s", tax_class, e$message))
  })
}

# Combine taxonomic summaries
if (length(taxonomic_summaries) > 0) {
  taxonomic_model_summaries <- bind_rows(taxonomic_summaries)
  print(taxonomic_model_summaries)

  saveRDS(taxonomic_models, file.path(data_outputs, "taxonomic_specific_models.rds"))
  saveRDS(taxonomic_model_summaries, file.path(data_outputs, "taxonomic_model_summaries.rds"))

  # Create coefficient comparison plot
  coef_comparison <- taxonomic_model_summaries %>%
    filter(term != "(Intercept)") %>%
    mutate(term = str_remove(term, "_scaled"))

  p_tax_coef <- ggplot(coef_comparison,
                       aes(x = term, y = estimate, color = taxonomic_class)) +
    geom_point(size = 3, position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                  width = 0.2, position = position_dodge(width = 0.5)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    scale_color_viridis_d(option = "turbo") +
    coord_flip() +
    labs(
      title = "Sampling Effort Predictors by Taxonomic Class",
      subtitle = "Coefficient estimates with 95% CI from Negative Binomial GLMs",
      x = "Predictor",
      y = "Coefficient Estimate",
      color = "Taxonomic Class"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")

  ggsave(file.path(figures_dir, "24_taxonomic_coefficients_comparison.png"),
         p_tax_coef, width = 12, height = 8, dpi = 300)
}

# 12. SUMMARY REPORT -----------------------------------------------------------
modeling_summary <- list(
  model_type = model_type,
  overdispersion = overdispersion,
  performance_count = perf_count,
  performance_presence = perf_presence,
  aic_comparison = aic_comparison,
  model_summaries = model_summaries,
  under_sampled_summary = under_sampled_summary,
  gam_summary = summary_gam,
  taxonomic_models_fitted = names(taxonomic_models),
  taxonomic_model_summaries = if (exists("taxonomic_model_summaries")) taxonomic_model_summaries else NULL
)

saveRDS(modeling_summary, file.path(data_outputs, "modeling_summary.rds"))

# Save final models
saveRDS(primary_model, file.path(data_outputs, "final_count_model.rds"))
saveRDS(glm_binomial, file.path(data_outputs, "final_presence_model.rds"))

message("\n=== Statistical modeling complete ===")
message("Results saved to: ", data_outputs)
message("Figures saved to: ", figures_dir)
