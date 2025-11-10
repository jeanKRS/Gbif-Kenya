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
    elevation = terra::extract(elevation, grid_coords)[, 1],
    temperature = terra::extract(temperature, grid_coords)[, 1],
    precipitation = terra::extract(precipitation, grid_coords)[, 1],
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
  # Filter out rows with NA values in any predictor
  filter(!is.na(elevation), !is.na(temperature), !is.na(precipitation),
         !is.na(dist_to_city_km), !is.na(dist_from_equator), !is.na(dist_from_coast_km)) %>%
  # Filter out rows with NA/NaN/Inf in scaled variables
  filter(!is.na(elevation_scaled), !is.na(temperature_scaled), !is.na(precipitation_scaled),
         !is.na(dist_to_city_scaled),
         is.finite(elevation_scaled), is.finite(temperature_scaled),
         is.finite(precipitation_scaled), is.finite(dist_to_city_scaled),
         is.finite(dist_from_equator), is.finite(dist_from_coast_km))

# Data validation report
message(sprintf("  Data after cleaning: %d grid cells retained", nrow(grid_model)))
message(sprintf("  Summary of predictor ranges:"))
message(sprintf("    Elevation: %.1f to %.1f m", min(grid_model$elevation), max(grid_model$elevation)))
message(sprintf("    Temperature: %.1f to %.1f °C", min(grid_model$temperature), max(grid_model$temperature)))
message(sprintf("    Precipitation: %.1f to %.1f mm", min(grid_model$precipitation), max(grid_model$precipitation)))
message(sprintf("    Distance to city: %.1f to %.1f km", min(grid_model$dist_to_city_km), max(grid_model$dist_to_city_km)))

# 3. GENERALIZED LINEAR MODELS -------------------------------------------------
message("\n=== Fitting GLMs for sampling effort ===")

# Check if models already exist
models_exist <- file.exists(file.path(data_outputs, "final_count_model.rds")) &&
                file.exists(file.path(data_outputs, "final_presence_model.rds"))

if (models_exist) {
  message("  ℹ GLM models already fitted. Loading from files...")
  primary_model <- readRDS(file.path(data_outputs, "final_count_model.rds"))
  glm_binomial <- readRDS(file.path(data_outputs, "final_presence_model.rds"))

  # Determine model type
  if (inherits(primary_model, "negbin")) {
    model_type <- "Negative Binomial"
    overdispersion <- primary_model$theta
  } else {
    model_type <- "Poisson"
    overdispersion <- sum(residuals(primary_model, type = "pearson")^2) / primary_model$df.residual
  }
  message(sprintf("  Model type: %s", model_type))
} else {
  message("  → Fitting GLM models...")

  # Model 1: Poisson GLM for record count
  message("  Fitting Poisson GLM...")
  glm_poisson <- glm(
    n_records ~ elevation_scaled + temperature_scaled + precipitation_scaled +
      dist_to_city_scaled + dist_from_equator + dist_from_coast_km,
    data = grid_model,
    family = poisson(link = "log"),
    control = glm.control(epsilon = 1e-10, maxit = 200, trace = FALSE)
  )

  # Check for overdispersion
  overdispersion <- sum(residuals(glm_poisson, type = "pearson")^2) / glm_poisson$df.residual
  message(sprintf("  Overdispersion parameter: %.2f", overdispersion))

  # Model 2: Negative binomial GLM (if overdispersed)
  if (overdispersion > 2) {
    message("  Data is overdispersed. Fitting Negative Binomial GLM...")

    # Add data diagnostics
    message(sprintf("    Data summary: mean = %.2f, variance = %.2f, zeros = %d (%.1f%%)",
                    mean(grid_model$n_records),
                    var(grid_model$n_records),
                    sum(grid_model$n_records == 0),
                    100 * sum(grid_model$n_records == 0) / nrow(grid_model)))

    # Try fitting with error handling and fallback options
    glm_nb <- NULL
    fit_successful <- FALSE

    # Attempt 1: Full model
    tryCatch({
      glm_nb <- glm.nb(
        n_records ~ elevation_scaled + temperature_scaled + precipitation_scaled +
          dist_to_city_scaled + dist_from_equator + dist_from_coast_km,
        data = grid_model,
        control = glm.control(epsilon = 1e-10, maxit = 200, trace = FALSE)
      )
      fit_successful <- TRUE
      message("    ✓ Full model converged")
    }, error = function(e) {
      message(sprintf("    First attempt failed: %s", e$message))
    })

    # Attempt 2: Simpler model (main environmental variables only)
    if (!fit_successful) {
      message("    Trying with simpler model (environmental variables only)...")
      tryCatch({
        glm_nb <- glm.nb(
          n_records ~ elevation_scaled + temperature_scaled + precipitation_scaled,
          data = grid_model,
          control = glm.control(epsilon = 1e-10, maxit = 200, trace = FALSE)
        )
        fit_successful <- TRUE
        message("    ✓ Simpler model converged")
      }, error = function(e) {
        message(sprintf("    Simpler model also failed: %s", e$message))
      })
    }

    # Attempt 3: Use starting values from Poisson model
    if (!fit_successful) {
      message("    Trying with starting values from Poisson model...")
      tryCatch({
        # Get starting values from Poisson model
        poisson_coefs <- coef(glm_poisson)

        glm_nb <- glm.nb(
          n_records ~ elevation_scaled + temperature_scaled + precipitation_scaled +
            dist_to_city_scaled + dist_from_equator + dist_from_coast_km,
          data = grid_model,
          init.theta = 1,
          control = glm.control(epsilon = 1e-10, maxit = 200, trace = FALSE)
        )
        fit_successful <- TRUE
        message("    ✓ Model with starting values converged")
      }, error = function(e) {
        message(sprintf("    Model with starting values failed: %s", e$message))
      })
    }

    # Attempt 4: GAM with negative binomial family as final fallback
    if (!fit_successful) {
      message("    Trying GAM with negative binomial family as fallback...")
      tryCatch({
        require(mgcv)
        glm_nb <- gam(
          n_records ~ s(elevation_scaled, k = 4) + s(temperature_scaled, k = 4) +
            s(precipitation_scaled, k = 4) + s(dist_to_city_scaled, k = 4),
          data = grid_model,
          family = nb(),
          method = "REML"
        )
        fit_successful <- TRUE
        message("    ✓ GAM model converged")
      }, error = function(e) {
        message(sprintf("    GAM also failed: %s", e$message))
      })
    }

    # If all attempts failed, fall back to Poisson
    if (!fit_successful || is.null(glm_nb)) {
      message("    ⚠ All Negative Binomial attempts failed. Using Poisson model.")
      message("    Note: Results may be unreliable due to extreme overdispersion.")
      primary_model <- glm_poisson
      model_type <- "Poisson (fallback)"
    } else {
      primary_model <- glm_nb
      if (inherits(glm_nb, "gam")) {
        model_type <- "Negative Binomial GAM"
      } else {
        model_type <- "Negative Binomial"
      }
    }
  } else {
    primary_model <- glm_poisson
    model_type <- "Poisson"
  }

  # Model 3: Binomial GLM for presence/absence
  message("  Fitting Binomial GLM for sampling presence...")
  glm_binomial <- glm(
    presence ~ elevation_scaled + temperature_scaled + precipitation_scaled +
      dist_to_city_scaled + dist_from_equator + dist_from_coast_km,
    data = grid_model,
    family = binomial(link = "logit"),
    control = glm.control(epsilon = 1e-10, maxit = 200, trace = FALSE)
  )

  message("  ✓ GLM models fitted")
}

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

# Helper function to safely extract model summaries with Wald CIs
safe_tidy_with_ci <- function(model, model_name) {
  # Get basic tidy output without CI
  tidy_result <- broom::tidy(model)

  # Calculate Wald confidence intervals (more robust than profile CIs)
  tryCatch({
    ci <- confint.default(model)  # Uses Wald method instead of profile
    tidy_result$conf.low <- ci[, 1]
    tidy_result$conf.high <- ci[, 2]
    message(sprintf("  ✓ Calculated Wald confidence intervals for %s", model_name))
  }, error = function(e) {
    message(sprintf("  ⚠ Could not calculate CIs for %s: %s", model_name, e$message))
    # Calculate approximate CIs using standard errors
    tidy_result$conf.low <- tidy_result$estimate - 1.96 * tidy_result$std.error
    tidy_result$conf.high <- tidy_result$estimate + 1.96 * tidy_result$std.error
  })

  return(tidy_result)
}

# Summary for count model
summary_count <- safe_tidy_with_ci(primary_model, "Count Model") %>%
  mutate(
    model = "Count Model",
    significant = p.value < 0.05,
    effect_direction = ifelse(estimate > 0, "Positive", "Negative")
  )

# Summary for presence model
summary_presence <- safe_tidy_with_ci(glm_binomial, "Presence Model") %>%
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

# Only perform model selection for GLM models (not GAM)
if (!inherits(primary_model, "gam")) {
  # Pre-flight check: Verify data integrity before model selection
  message("  Checking data integrity before model selection...")

  predictor_cols <- c("elevation_scaled", "temperature_scaled", "precipitation_scaled",
                     "dist_to_city_scaled", "dist_from_equator", "dist_from_coast_km")

  # Check for any NA/NaN/Inf values
  data_check <- sapply(predictor_cols, function(col) {
    values <- grid_model[[col]]
    list(
      has_na = any(is.na(values)),
      has_inf = any(is.infinite(values)),
      n_valid = sum(is.finite(values) & !is.na(values))
    )
  })

  # Report any issues
  if (any(sapply(data_check["has_na", ], identity))) {
    warning("NA values detected in predictors")
    print(data_check)
  }

  if (any(sapply(data_check["has_inf", ], identity))) {
    warning("Inf values detected in predictors")
    print(data_check)
  }

  # Create clean dataset for model selection
  grid_model_clean <- grid_model %>%
    filter(if_all(all_of(predictor_cols), ~ is.finite(.) & !is.na(.)))

  message(sprintf("  Using %d of %d observations for model selection (%.1f%%)",
                  nrow(grid_model_clean), nrow(grid_model),
                  100 * nrow(grid_model_clean) / nrow(grid_model)))

  # Refit primary model on clean data if needed
  if (nrow(grid_model_clean) < nrow(grid_model)) {
    message("  Refitting primary model on cleaned dataset...")

    if (inherits(primary_model, "negbin")) {
      primary_model <- glm.nb(
        n_records ~ elevation_scaled + temperature_scaled + precipitation_scaled +
          dist_to_city_scaled + dist_from_equator + dist_from_coast_km,
        data = grid_model_clean,
        control = glm.control(epsilon = 1e-10, maxit = 200, trace = FALSE)
      )
    } else {
      primary_model <- glm(
        n_records ~ elevation_scaled + temperature_scaled + precipitation_scaled +
          dist_to_city_scaled + dist_from_equator + dist_from_coast_km,
        data = grid_model_clean,
        family = poisson(link = "log"),
        control = glm.control(epsilon = 1e-10, maxit = 200, trace = FALSE)
      )
    }
  }

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
} else {
  message("  ℹ Skipping model selection (primary model is GAM)")
  message("  Using AIC from primary model: ", AIC(primary_model))
}

# 7. GENERALIZED ADDITIVE MODELS -----------------------------------------------
message("\n=== Fitting GAMs for non-linear relationships ===")

# Check if GAM already exists
if (file.exists(file.path(data_outputs, "gam_model.rds"))) {
  message("  ℹ GAM model already fitted. Loading from file...")
  gam_count <- readRDS(file.path(data_outputs, "gam_model.rds"))
  summary_gam <- summary(gam_count)
  print(summary_gam)
} else {
  message("  → Fitting GAM model...")

  # Try fitting GAM with error handling
  gam_fit_successful <- FALSE

  # Attempt 1: Full GAM with negative binomial
  tryCatch({
    gam_count <- gam(
      n_records ~ s(elevation, k = 5) + s(temperature, k = 5) +
        s(precipitation, k = 5) + s(dist_to_city_km, k = 5) +
        s(lon, lat, k = 20),
      data = grid_model %>% filter(n_records > 0),
      family = nb(),
      method = "REML"
    )
    gam_fit_successful <- TRUE
    message("    ✓ Negative binomial GAM converged")
  }, error = function(e) {
    message(sprintf("    Negative binomial GAM failed: %s", e$message))
  })

  # Attempt 2: Try with Poisson family
  if (!gam_fit_successful) {
    message("    Trying GAM with Poisson family...")
    tryCatch({
      gam_count <- gam(
        n_records ~ s(elevation, k = 5) + s(temperature, k = 5) +
          s(precipitation, k = 5) + s(dist_to_city_km, k = 5) +
          s(lon, lat, k = 20),
        data = grid_model %>% filter(n_records > 0),
        family = poisson(link = "log"),
        method = "REML"
      )
      gam_fit_successful <- TRUE
      message("    ✓ Poisson GAM converged (note: may be overdispersed)")
    }, error = function(e) {
      message(sprintf("    Poisson GAM failed: %s", e$message))
    })
  }

  # Attempt 3: Simpler GAM with fewer knots
  if (!gam_fit_successful) {
    message("    Trying simpler GAM with fewer knots...")
    tryCatch({
      gam_count <- gam(
        n_records ~ s(elevation, k = 3) + s(temperature, k = 3) +
          s(precipitation, k = 3) + s(dist_to_city_km, k = 3),
        data = grid_model %>% filter(n_records > 0),
        family = nb(),
        method = "REML"
      )
      gam_fit_successful <- TRUE
      message("    ✓ Simpler GAM converged")
    }, error = function(e) {
      message(sprintf("    Simpler GAM failed: %s", e$message))
    })
  }

  if (gam_fit_successful) {
    # GAM summary
    summary_gam <- summary(gam_count)
    print(summary_gam)

    # Save GAM
    saveRDS(gam_count, file.path(data_outputs, "gam_model.rds"))
    message("  ✓ GAM model fitted")
  } else {
    message("  ⚠ All GAM fitting attempts failed. Skipping GAM analysis.")
    gam_count <- NULL
  }
}

# GAM diagnostic plots (only if GAM was successfully fitted)
if (!is.null(gam_count)) {
  png(file.path(figures_dir, "18_gam_diagnostics.png"),
      width = 3000, height = 2000, res = 300)
  par(mfrow = c(2, 3))
  gam.check(gam_count)
  dev.off()
} else {
  message("  ℹ Skipping GAM diagnostics (no GAM model available)")
}

# 8. EFFECT PLOTS --------------------------------------------------------------
message("\n=== Creating effect plots ===")

# Effect plots for GLM (only if primary model is not a GAM)
if (!inherits(primary_model, "gam")) {
  effects_count <- allEffects(primary_model)

  png(file.path(figures_dir, "19_glm_effects.png"),
      width = 3600, height = 2400, res = 300)
  plot(effects_count, main = paste("Effects on Sampling Effort -", model_type, "GLM"))
  dev.off()
} else {
  message("  ℹ Skipping allEffects (primary model is GAM - see GAM smooth plots instead)")
}

# GAM smooth plots (only if GAM was successfully fitted)
if (!is.null(gam_count)) {
  png(file.path(figures_dir, "20_gam_smooths.png"),
      width = 3600, height = 2400, res = 300)
  par(mfrow = c(2, 3))
  plot(gam_count, shade = TRUE, pages = 1, scale = 0,
       main = "GAM Smooth Effects on Sampling Effort")
  dev.off()
} else {
  message("  ℹ Skipping GAM smooth plots (no GAM model available)")
}

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

# Check if taxonomic models already exist
taxonomic_models_exist <- file.exists(file.path(data_outputs, "taxonomic_specific_models.rds")) &&
                          file.exists(file.path(data_outputs, "taxonomic_model_summaries.rds"))

if (taxonomic_models_exist) {
  message("  ℹ Taxonomic-specific models already fitted. Loading from files...")
  taxonomic_models <- readRDS(file.path(data_outputs, "taxonomic_specific_models.rds"))
  taxonomic_model_summaries <- readRDS(file.path(data_outputs, "taxonomic_model_summaries.rds"))
  message(sprintf("  Loaded models for %d taxonomic classes", length(taxonomic_models)))
  print(names(taxonomic_models))
} else {
  message("  → Fitting taxonomic-specific models...")

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

  message(sprintf("  Fitting models for top 5 classes: %s", paste(top_5_classes, collapse = ", ")))

  # Join taxonomic data with environmental grid
  grid_tax_env <- grid_model %>%
    select(grid_id, elevation_scaled, temperature_scaled, precipitation_scaled,
           dist_to_city_scaled, dist_from_equator, dist_from_coast_km,
           elevation, temperature, precipitation, dist_to_city_km) %>%
    # Ensure clean data for taxonomic models
    filter(is.finite(elevation_scaled), is.finite(temperature_scaled),
           is.finite(precipitation_scaled), is.finite(dist_to_city_scaled),
           is.finite(dist_from_equator), is.finite(dist_from_coast_km),
           !is.na(elevation_scaled), !is.na(temperature_scaled),
           !is.na(precipitation_scaled), !is.na(dist_to_city_scaled)) %>%
    inner_join(
      grid_taxonomy %>% filter(class %in% top_5_classes),
      by = "grid_id"
    )

  message(sprintf("  Grid cells with taxonomic data: %d", nrow(grid_tax_env)))

  # Fit separate models for each major taxonomic group
  taxonomic_models <- list()
  taxonomic_summaries <- list()

  for (tax_class in top_5_classes) {
    message(sprintf("    Fitting model for %s...", tax_class))

    class_data <- grid_tax_env %>%
      filter(class == tax_class)

    if (nrow(class_data) < 30) {
      message(sprintf("    Skipping %s - insufficient data (n = %d)", tax_class, nrow(class_data)))
      next
    }

    # Try fitting with multiple fallback strategies
    model_class <- NULL
    fit_successful <- FALSE

    # Attempt 1: Full negative binomial model
    tryCatch({
      model_class <- glm.nb(
        n_records_class ~ elevation_scaled + temperature_scaled + precipitation_scaled +
          dist_to_city_scaled + dist_from_equator + dist_from_coast_km,
        data = class_data,
        control = glm.control(epsilon = 1e-10, maxit = 200, trace = FALSE)
      )
      fit_successful <- TRUE
    }, error = function(e) {
      message(sprintf("      Full model failed: %s", e$message))
    })

    # Attempt 2: Simpler model
    if (!fit_successful) {
      message(sprintf("      Trying simpler model for %s...", tax_class))
      tryCatch({
        model_class <- glm.nb(
          n_records_class ~ elevation_scaled + temperature_scaled + precipitation_scaled,
          data = class_data,
          control = glm.control(epsilon = 1e-10, maxit = 200, trace = FALSE)
        )
        fit_successful <- TRUE
      }, error = function(e) {
        message(sprintf("      Simpler model failed: %s", e$message))
      })
    }

    # Attempt 3: Poisson model as fallback
    if (!fit_successful) {
      message(sprintf("      Trying Poisson model for %s...", tax_class))
      tryCatch({
        model_class <- glm(
          n_records_class ~ elevation_scaled + temperature_scaled + precipitation_scaled +
            dist_to_city_scaled + dist_from_equator + dist_from_coast_km,
          data = class_data,
          family = poisson(link = "log"),
          control = glm.control(epsilon = 1e-10, maxit = 200, trace = FALSE)
        )
        fit_successful <- TRUE
        message(sprintf("      Note: Using Poisson model for %s (NB failed to converge)", tax_class))
      }, error = function(e) {
        message(sprintf("      Poisson model failed: %s", e$message))
      })
    }

    # If successful, store the model and extract summary
    if (fit_successful && !is.null(model_class)) {
      taxonomic_models[[tax_class]] <- model_class

      # Extract summary using safe method with Wald CIs
      taxonomic_summaries[[tax_class]] <- safe_tidy_with_ci(model_class, paste("Taxonomic model:", tax_class)) %>%
        mutate(
          taxonomic_class = tax_class,
          significant = p.value < 0.05,
          effect_direction = ifelse(estimate > 0, "Positive", "Negative")
        )

      message(sprintf("    ✓ Completed model for %s", tax_class))
    } else {
      message(sprintf("    ✗ All modeling attempts failed for %s", tax_class))
    }
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

  message("  ✓ Taxonomic-specific models fitted")
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
