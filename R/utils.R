# ==============================================================================
# Utility Functions for GBIF Kenya Bias Assessment
# Author: Automated Research Pipeline
# Date: 2025-11-10
# ==============================================================================

#' Load all required packages
#'
#' @param quietly Logical, suppress package startup messages
#' @return NULL (invisible)
load_packages <- function(quietly = TRUE) {
  packages <- c(
    "tidyverse", "sf", "here", "rgbif", "occAssess",
    "rnaturalearth", "rnaturalearthdata", "terra", "geodata",
    "CoordinateCleaner", "lubridate", "vegan", "iNEXT",
    "MASS", "lme4", "DHARMa", "performance", "mgcv",
    "viridis", "patchwork", "scales", "Kendall"
  )

  if (quietly) {
    suppressPackageStartupMessages({
      invisible(lapply(packages, library, character.only = TRUE))
    })
  } else {
    invisible(lapply(packages, library, character.only = TRUE))
  }

  message("All packages loaded successfully")
  invisible(NULL)
}


#' Calculate grid-based sampling metrics
#'
#' @param data sf object with occurrence data
#' @param grid_size numeric, grid cell size in degrees
#' @param boundary sf object with study area boundary
#' @return sf object with grid and metrics
calculate_grid_metrics <- function(data, grid_size = 0.1, boundary = NULL) {

  if (is.null(boundary)) {
    boundary <- st_union(data) %>% st_convex_hull()
  }

  # Create grid
  grid <- st_make_grid(boundary, cellsize = grid_size, square = FALSE) %>%
    st_sf() %>%
    mutate(grid_id = row_number()) %>%
    st_intersection(boundary)

  # Calculate metrics
  grid_metrics <- grid %>%
    st_join(data) %>%
    group_by(grid_id) %>%
    summarise(
      n_records = sum(!is.na(species)),
      n_species = n_distinct(species, na.rm = TRUE),
      n_families = n_distinct(family, na.rm = TRUE),
      .groups = "drop"
    )

  # Join back to complete grid
  grid_complete <- grid %>%
    left_join(st_drop_geometry(grid_metrics), by = "grid_id") %>%
    mutate(
      n_records = replace_na(n_records, 0),
      n_species = replace_na(n_species, 0),
      n_families = replace_na(n_families, 0)
    )

  return(grid_complete)
}


#' Extract environmental data for points
#'
#' @param coords matrix of longitude, latitude coordinates
#' @param country character, ISO3 country code
#' @return data.frame with environmental values
extract_env_data <- function(coords, country = "KEN") {

  # Download data (geodata functions return SpatRaster objects directly)
  elev_path <- geodata::elevation_30s(country = country, path = tempdir())
  clim_path <- geodata::worldclim_country(country = country, var = "bio", path = tempdir())

  # Load rasters properly - check if result is path or SpatRaster
  if (is.character(elev_path)) {
    elev <- terra::rast(elev_path)
  } else {
    elev <- elev_path
  }

  if (is.character(clim_path)) {
    clim <- terra::rast(clim_path)
  } else {
    clim <- clim_path
  }

  # Extract values
  env_data <- data.frame(
    elevation = terra::extract(elev, coords)[, 2],
    bio1_temp = terra::extract(clim[[1]], coords)[, 2],
    bio12_precip = terra::extract(clim[[12]], coords)[, 2],
    bio15_precip_seasonality = terra::extract(clim[[15]], coords)[, 2]
  )

  return(env_data)
}


#' Calculate Gini coefficient
#'
#' @param x numeric vector
#' @return numeric, Gini coefficient
calculate_gini <- function(x) {
  x <- x[!is.na(x)]
  x <- sort(x)
  n <- length(x)
  G <- (2 * sum((1:n) * x)) / (n * sum(x)) - (n + 1) / n
  return(G)
}


#' Calculate completeness index
#'
#' @param observed numeric, observed richness
#' @param expected numeric, expected richness
#' @return numeric, completeness index (0-1)
calculate_completeness <- function(observed, expected) {
  completeness <- observed / expected
  completeness[completeness > 1] <- 1
  return(completeness)
}


#' Format numbers for reporting
#'
#' @param x numeric vector
#' @param digits integer, number of decimal places
#' @return character vector
format_number <- function(x, digits = 2) {
  format(round(x, digits), big.mark = ",", scientific = FALSE)
}


#' Create summary table
#'
#' @param data data.frame
#' @param group_vars character vector of grouping variables
#' @param summary_vars character vector of variables to summarize
#' @return data.frame
create_summary_table <- function(data, group_vars = NULL, summary_vars = NULL) {

  if (is.null(summary_vars)) {
    summary_vars <- names(data)[sapply(data, is.numeric)]
  }

  if (is.null(group_vars)) {
    summary <- data %>%
      summarise(across(all_of(summary_vars),
                      list(mean = ~mean(., na.rm = TRUE),
                           median = ~median(., na.rm = TRUE),
                           sd = ~sd(., na.rm = TRUE),
                           min = ~min(., na.rm = TRUE),
                           max = ~max(., na.rm = TRUE)),
                      .names = "{.col}_{.fn}"))
  } else {
    summary <- data %>%
      group_by(across(all_of(group_vars))) %>%
      summarise(across(all_of(summary_vars),
                      list(mean = ~mean(., na.rm = TRUE),
                           median = ~median(., na.rm = TRUE),
                           sd = ~sd(., na.rm = TRUE)),
                      .names = "{.col}_{.fn}"),
               .groups = "drop")
  }

  return(summary)
}


#' Plot spatial distribution
#'
#' @param data sf object
#' @param fill_var character, variable to map
#' @param title character, plot title
#' @param palette character, viridis palette option
#' @return ggplot object
plot_spatial <- function(data, fill_var, title = "", palette = "viridis") {

  p <- ggplot(data) +
    geom_sf(aes(fill = .data[[fill_var]]), color = NA) +
    scale_fill_viridis_c(option = palette, name = fill_var) +
    labs(title = title) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      axis.text = element_text(size = 8)
    )

  return(p)
}


#' Run data quality checks
#'
#' @param data data.frame with GBIF data
#' @return list with quality metrics
quality_check <- function(data) {

  checks <- list(
    total_records = nrow(data),
    missing_coords = sum(is.na(data$decimalLongitude) | is.na(data$decimalLatitude)),
    missing_species = sum(is.na(data$species) | data$species == ""),
    missing_dates = sum(is.na(data$eventDate)),
    duplicate_records = nrow(data) - nrow(distinct(data, species, decimalLongitude,
                                                   decimalLatitude, eventDate)),
    invalid_years = sum(data$year < 1700 | data$year > year(Sys.Date()), na.rm = TRUE),
    high_uncertainty = sum(data$coordinateUncertaintyInMeters > 10000, na.rm = TRUE)
  )

  checks$quality_score <- 100 * (1 - (checks$missing_coords + checks$missing_species +
                                     checks$missing_dates + checks$duplicate_records) /
                                   (checks$total_records * 4))

  return(checks)
}


#' Calculate spatial autocorrelation
#'
#' @param values numeric vector of values
#' @param coords matrix of coordinates
#' @param k integer, number of nearest neighbors
#' @return list with Moran's I results
calculate_morans_i <- function(values, coords, k = 8) {

  # Create neighbors
  nb <- spdep::knn2nb(spdep::knearneigh(coords, k = k))
  listw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)

  # Calculate Moran's I
  moran <- spdep::moran.test(values, listw, zero.policy = TRUE)

  result <- list(
    morans_i = moran$estimate["Moran I statistic"],
    expected = moran$estimate["Expectation"],
    variance = moran$estimate["Variance"],
    p_value = moran$p.value,
    statistic = moran$statistic
  )

  return(result)
}


#' Create publication-ready table
#'
#' @param data data.frame
#' @param caption character, table caption
#' @param digits integer, decimal places
#' @return kable object
create_table <- function(data, caption = "", digits = 2) {

  if (requireNamespace("knitr", quietly = TRUE) &&
      requireNamespace("kableExtra", quietly = TRUE)) {

    table <- data %>%
      knitr::kable(caption = caption, digits = digits, format = "html") %>%
      kableExtra::kable_styling(bootstrap_options = c("striped", "hover", "condensed"))

    return(table)
  } else {
    return(data)
  }
}


#' Save plot with standard settings
#'
#' @param plot ggplot object
#' @param filename character, output filename
#' @param width numeric, plot width in inches
#' @param height numeric, plot height in inches
#' @param dpi numeric, resolution
#' @return NULL (invisible)
save_plot <- function(plot, filename, width = 10, height = 8, dpi = 300) {

  ggsave(filename, plot, width = width, height = height, dpi = dpi)
  message("Plot saved to: ", filename)
  invisible(NULL)
}


#' Generate citation for GBIF data
#'
#' @param download_key character, GBIF download key
#' @param access_date Date, date of access
#' @return character, formatted citation
gbif_citation <- function(download_key, access_date = Sys.Date()) {

  citation <- paste0(
    "GBIF.org (", format(access_date, "%d %B %Y"), ") ",
    "GBIF Occurrence Download https://doi.org/10.15468/dl.", download_key
  )

  return(citation)
}


#' Validate coordinate data
#'
#' @param data data.frame with coordinates
#' @param lon_col character, longitude column name
#' @param lat_col character, latitude column name
#' @return logical vector of valid coordinates
validate_coordinates <- function(data, lon_col = "decimalLongitude",
                                lat_col = "decimalLatitude") {

  valid <- !is.na(data[[lon_col]]) &
           !is.na(data[[lat_col]]) &
           data[[lon_col]] >= -180 &
           data[[lon_col]] <= 180 &
           data[[lat_col]] >= -90 &
           data[[lat_col]] <= 90 &
           !(data[[lon_col]] == 0 & data[[lat_col]] == 0)

  return(valid)
}


#' Calculate temporal completeness
#'
#' @param years numeric vector of years
#' @return list with temporal metrics
temporal_metrics <- function(years) {

  years <- years[!is.na(years)]
  year_range <- seq(min(years), max(years))

  metrics <- list(
    min_year = min(years),
    max_year = max(years),
    span = max(years) - min(years) + 1,
    n_years_sampled = length(unique(years)),
    completeness = length(unique(years)) / (max(years) - min(years) + 1),
    n_gaps = length(year_range) - length(unique(years)),
    cv = sd(table(years)) / mean(table(years))
  )

  return(metrics)
}

message("Utility functions loaded successfully")
