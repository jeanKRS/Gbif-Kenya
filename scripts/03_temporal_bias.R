# ==============================================================================
# Script: 03_temporal_bias.R
# Purpose: Assess temporal biases in GBIF data for Kenya
# Author: Kwiz Computing Technologies
# Date: 2025-11-10
# ==============================================================================

# Load required packages -------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(occAssess)
  library(lubridate)
  library(Kendall)
  library(forecast)
  library(tseries)
  library(viridis)
  library(patchwork)
})

# Setup paths ------------------------------------------------------------------
data_processed <- here("data", "processed")
data_outputs <- here("data", "outputs")
figures_dir <- here("figures")

# Load cleaned data ------------------------------------------------------------
message("Loading cleaned GBIF data...")
kenya_data <- readRDS(file.path(data_processed, "kenya_gbif_clean.rds"))

# 1. TEMPORAL COVERAGE ANALYSIS ------------------------------------------------
message("\n=== Analyzing temporal coverage ===")

# Aggregate data by year
temporal_summary <- kenya_data %>%
  filter(!is.na(year)) %>%
  group_by(year) %>%
  summarise(
    n_records = n(),
    n_species = n_distinct(species),
    n_genera = n_distinct(genus),
    n_families = n_distinct(family),
    n_collectors = n_distinct(recordedBy, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(year)

# Overall temporal statistics
temporal_stats <- temporal_summary %>%
  summarise(
    year_min = min(year),
    year_max = max(year),
    year_range = year_max - year_min,
    total_years = n_distinct(year),
    mean_records_per_year = mean(n_records),
    median_records_per_year = median(n_records),
    sd_records_per_year = sd(n_records),
    cv_records = sd_records_per_year / mean_records_per_year,
    mean_species_per_year = mean(n_species),
    median_species_per_year = median(n_species)
  )

print(temporal_stats)
saveRDS(temporal_stats, file.path(data_outputs, "temporal_statistics.rds"))

# 2. TEMPORAL TRENDS -----------------------------------------------------------
message("\n=== Testing temporal trends ===")

# Mann-Kendall trend test for records
mk_records <- MannKendall(temporal_summary$n_records)
mk_species <- MannKendall(temporal_summary$n_species)

trend_tests <- data.frame(
  variable = c("n_records", "n_species"),
  tau = c(mk_records$tau, mk_species$tau),
  p_value = c(mk_records$sl, mk_species$sl),
  trend = c(
    ifelse(mk_records$sl < 0.05, ifelse(mk_records$tau > 0, "Increasing", "Decreasing"), "No trend"),
    ifelse(mk_species$sl < 0.05, ifelse(mk_species$tau > 0, "Increasing", "Decreasing"), "No trend")
  )
)

print(trend_tests)
saveRDS(trend_tests, file.path(data_outputs, "temporal_trends.rds"))

# 3. TEMPORAL GAPS AND COMPLETENESS -------------------------------------------
message("\n=== Identifying temporal gaps ===")

# Identify years with no records
year_range <- seq(min(temporal_summary$year), max(temporal_summary$year))
missing_years <- setdiff(year_range, temporal_summary$year)

# Create complete temporal sequence
temporal_complete <- tibble(year = year_range) %>%
  left_join(temporal_summary, by = "year") %>%
  mutate(
    n_records = replace_na(n_records, 0),
    n_species = replace_na(n_species, 0),
    has_records = n_records > 0
  )

# Gap analysis
gap_summary <- temporal_complete %>%
  summarise(
    total_years = n(),
    years_with_records = sum(has_records),
    years_without_records = sum(!has_records),
    temporal_completeness = 100 * years_with_records / total_years,
    max_gap_length = max(rle(!has_records)$lengths[rle(!has_records)$values]),
    n_gaps = sum(rle(!has_records)$values)
  )

print(gap_summary)
saveRDS(gap_summary, file.path(data_outputs, "temporal_gaps.rds"))

# 4. OCCASSESS TEMPORAL ASSESSMENTS --------------------------------------------
message("\n=== Running occAssess temporal assessments ===")

# Prepare data for occAssess with taxonomic levels
kenya_df <- kenya_data %>%
  filter(!is.na(eventDate)) %>%
  select(species, genus, family, order, class, phylum, kingdom,
         decimalLongitude, decimalLatitude, eventDate, year) %>%
  as.data.frame()

# Create 10-year period breaks (consistent with spatial analysis)
year_range <- range(kenya_df$year, na.rm = TRUE)
period_breaks <- seq(
  floor(year_range[1] / 10) * 10,
  ceiling(year_range[2] / 10) * 10,
  by = 10
)

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

# Define taxonomic levels for analysis
taxonomic_levels <- c("species", "genus", "family", "order", "class", "phylum", "kingdom")

# Check if all temporal assessments already exist
all_temporal_assessments_exist <- file.exists(file.path(data_outputs, "occassess_record_time_by_taxlevel.rds")) &&
                                  file.exists(file.path(data_outputs, "occassess_species_time_by_taxlevel.rds")) &&
                                  file.exists(file.path(data_outputs, "temporal_assessment_summary.rds"))

if (all_temporal_assessments_exist) {
  message("  ℹ All temporal occAssess assessments already completed. Loading from files...")
  record_time_assessments <- readRDS(file.path(data_outputs, "occassess_record_time_by_taxlevel.rds"))
  species_time_assessments <- readRDS(file.path(data_outputs, "occassess_species_time_by_taxlevel.rds"))
  temporal_assessment_summary <- readRDS(file.path(data_outputs, "temporal_assessment_summary.rds"))
  print(temporal_assessment_summary)
} else {
  message("  → Running temporal occAssess assessments...")

  # Record time assessment - Run for each taxonomic level
  message("\n--- Assessing record temporal distribution by taxonomic level ---")
  record_time_assessments <- list()

  for (tax_level in taxonomic_levels) {
    message(sprintf("  Running assessRecordTime for %s", tax_level))

    if (all(is.na(kenya_df[[tax_level]]))) {
      message(sprintf("  Skipping %s - all values are NA", tax_level))
      next
    }

    tryCatch({
      record_time_assessments[[tax_level]] <- assessRecordTime(
        dat = kenya_df,
        dateCol = "eventDate",
        intervals = "year"
      )
      message(sprintf("  ✓ Completed time assessment for %s", tax_level))
    }, error = function(e) {
      message(sprintf("  ✗ Error in time assessment for %s: %s", tax_level, e$message))
    })
  }

  # Species temporal coverage - Run for different taxonomic levels
  message("\n--- Assessing temporal coverage by taxonomic level ---")
  species_time_assessments <- list()

  for (tax_level in taxonomic_levels) {
    message(sprintf("  Running assessSpeciesTime for %s", tax_level))

    if (all(is.na(kenya_df[[tax_level]]))) {
      message(sprintf("  Skipping %s - all values are NA", tax_level))
      next
    }

    tryCatch({
      species_time_assessments[[tax_level]] <- assessSpeciesTime(
        dat = kenya_df,
        speciesCol = tax_level,
        dateCol = "eventDate",
        threshold = 5  # Minimum 5 years between first and last record
      )
      message(sprintf("  ✓ Completed temporal coverage for %s", tax_level))
    }, error = function(e) {
      message(sprintf("  ✗ Error in temporal coverage for %s: %s", tax_level, e$message))
    })
  }

  # Save occAssess results
  message("\n--- Saving occAssess temporal results ---")
  saveRDS(record_time_assessments, file.path(data_outputs, "occassess_record_time_by_taxlevel.rds"))
  saveRDS(species_time_assessments, file.path(data_outputs, "occassess_species_time_by_taxlevel.rds"))

  # Also save species-level assessment for backward compatibility
  if (!is.null(record_time_assessments[["species"]])) {
    saveRDS(record_time_assessments[["species"]], file.path(data_outputs, "occassess_record_time.rds"))
  }
  if (!is.null(species_time_assessments[["species"]])) {
    saveRDS(species_time_assessments[["species"]], file.path(data_outputs, "occassess_species_time.rds"))
  }

  # Create summary of temporal assessments
  temporal_assessment_summary <- data.frame(
    taxonomic_level = taxonomic_levels,
    record_time_assessment = sapply(taxonomic_levels, function(x) !is.null(record_time_assessments[[x]])),
    species_time_assessment = sapply(taxonomic_levels, function(x) !is.null(species_time_assessments[[x]]))
  )

  print(temporal_assessment_summary)
  saveRDS(temporal_assessment_summary, file.path(data_outputs, "temporal_assessment_summary.rds"))
  message("  ✓ Temporal occAssess assessments completed")
}

# 5. TEMPORAL BIAS BY TAXONOMIC GROUP ------------------------------------------
message("\n=== Analyzing temporal bias by taxonomic group ===")

# Temporal trends by class
temporal_by_class <- kenya_data %>%
  filter(!is.na(year), !is.na(class)) %>%
  group_by(year, class) %>%
  summarise(
    n_records = n(),
    n_species = n_distinct(species),
    .groups = "drop"
  )

# Get top 10 classes by total records
top_classes <- kenya_data %>%
  count(class, sort = TRUE) %>%
  head(10) %>%
  pull(class)

temporal_top_classes <- temporal_by_class %>%
  filter(class %in% top_classes)

# Calculate sampling effort change over decades
decade_analysis <- kenya_data %>%
  filter(!is.na(year)) %>%
  mutate(decade = floor(year / 10) * 10) %>%
  group_by(decade) %>%
  summarise(
    n_records = n(),
    n_species = n_distinct(species),
    n_families = n_distinct(family),
    .groups = "drop"
  )

print(decade_analysis)
saveRDS(decade_analysis, file.path(data_outputs, "decade_analysis.rds"))

# 6. SEASONAL PATTERNS ---------------------------------------------------------
message("\n=== Analyzing seasonal patterns ===")

# Monthly patterns
monthly_patterns <- kenya_data %>%
  filter(!is.na(month)) %>%
  group_by(month) %>%
  summarise(
    n_records = n(),
    n_species = n_distinct(species),
    .groups = "drop"
  ) %>%
  mutate(
    month_name = month.abb[month],
    season = case_when(
      month %in% c(12, 1, 2) ~ "DJF",  # Long dry season
      month %in% c(3, 4, 5) ~ "MAM",   # Long rains
      month %in% c(6, 7, 8) ~ "JJA",   # Short dry season
      month %in% c(9, 10, 11) ~ "SON"  # Short rains
    )
  )

# Test for seasonal bias (Chi-square test)
seasonal_test <- chisq.test(monthly_patterns$n_records)

seasonal_summary <- data.frame(
  chi_square = seasonal_test$statistic,
  p_value = seasonal_test$p.value,
  seasonal_bias = ifelse(seasonal_test$p.value < 0.05, "Significant", "Not significant")
)

print(seasonal_summary)
saveRDS(monthly_patterns, file.path(data_outputs, "monthly_patterns.rds"))
saveRDS(seasonal_summary, file.path(data_outputs, "seasonal_bias_test.rds"))

# 7. SPECIES-LEVEL TEMPORAL COVERAGE -------------------------------------------
message("\n=== Analyzing species-level temporal coverage ===")

# Calculate temporal span for each species
species_temporal <- kenya_data %>%
  filter(!is.na(year)) %>%
  group_by(species) %>%
  summarise(
    n_records = n(),
    year_first = min(year),
    year_last = max(year),
    temporal_span = year_last - year_first,
    n_years = n_distinct(year),
    temporal_completeness = n_years / (temporal_span + 1),
    .groups = "drop"
  ) %>%
  filter(n_records >= 10)  # Only species with at least 10 records

# Summary statistics
species_temporal_summary <- species_temporal %>%
  summarise(
    n_species = n(),
    mean_temporal_span = mean(temporal_span),
    median_temporal_span = median(temporal_span),
    mean_completeness = mean(temporal_completeness),
    median_completeness = median(temporal_completeness)
  )

print(species_temporal_summary)
saveRDS(species_temporal, file.path(data_outputs, "species_temporal_coverage.rds"))

# 8. VISUALIZATION -------------------------------------------------------------
message("\n=== Creating temporal visualizations ===")

# Plot 1: Temporal trends in records and species
p1 <- ggplot(temporal_summary, aes(x = year)) +
  geom_line(aes(y = n_records, color = "Records"), linewidth = 1) +
  geom_line(aes(y = n_species * 10, color = "Species (×10)"), linewidth = 1) +
  scale_color_manual(values = c("Records" = "blue", "Species (×10)" = "red")) +
  scale_y_continuous(
    name = "Number of Records",
    sec.axis = sec_axis(~./10, name = "Number of Species")
  ) +
  labs(
    title = "Temporal Trends in GBIF Data for Kenya",
    subtitle = paste0(min(temporal_summary$year), " - ", max(temporal_summary$year)),
    x = "Year",
    color = ""
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.title.y.right = element_text(color = "red"),
    axis.text.y.right = element_text(color = "red"),
    axis.title.y.left = element_text(color = "blue"),
    axis.text.y.left = element_text(color = "blue")
  )

ggsave(file.path(figures_dir, "05_temporal_trends.png"),
       p1, width = 12, height = 6, dpi = 300)

# Plot 2: Decadal comparison
p2 <- ggplot(decade_analysis, aes(x = factor(decade), y = n_records)) +
  geom_col(fill = "steelblue", alpha = 0.8) +
  geom_text(aes(label = scales::comma(n_records)), vjust = -0.5, size = 3) +
  labs(
    title = "Decadal Distribution of Records",
    x = "Decade",
    y = "Number of Records"
  ) +
  scale_y_continuous(labels = scales::comma) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(figures_dir, "06_decadal_records.png"),
       p2, width = 10, height = 6, dpi = 300)

# Plot 3: Seasonal patterns
p3 <- ggplot(monthly_patterns, aes(x = factor(month), y = n_records)) +
  geom_col(aes(fill = season), alpha = 0.8) +
  scale_fill_viridis_d(option = "plasma") +
  scale_x_discrete(labels = month.abb) +
  labs(
    title = "Seasonal Distribution of Records",
    subtitle = "Grouped by East African Seasonal Patterns",
    x = "Month",
    y = "Number of Records",
    fill = "Season"
  ) +
  scale_y_continuous(labels = scales::comma) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(file.path(figures_dir, "07_seasonal_patterns.png"),
       p3, width = 10, height = 6, dpi = 300)

# Plot 4: Temporal trends by taxonomic class
p4 <- ggplot(temporal_top_classes, aes(x = year, y = n_records, color = class)) +
  geom_line(linewidth = 1, alpha = 0.7) +
  scale_color_viridis_d(option = "turbo") +
  scale_y_log10(labels = scales::comma) +
  labs(
    title = "Temporal Trends by Taxonomic Class",
    subtitle = "Top 10 classes by total records (log scale)",
    x = "Year",
    y = "Number of Records (log scale)",
    color = "Class"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

ggsave(file.path(figures_dir, "08_temporal_by_class.png"),
       p4, width = 12, height = 6, dpi = 300)

# Plot 5: Species temporal coverage histogram
p5 <- ggplot(species_temporal, aes(x = temporal_span)) +
  geom_histogram(bins = 50, fill = "darkgreen", alpha = 0.7) +
  labs(
    title = "Distribution of Species Temporal Coverage",
    subtitle = paste0("Species with ≥10 records (n = ", nrow(species_temporal), ")"),
    x = "Temporal Span (years)",
    y = "Number of Species"
  ) +
  theme_minimal()

ggsave(file.path(figures_dir, "09_species_temporal_span.png"),
       p5, width = 10, height = 6, dpi = 300)

# 9. SUMMARY REPORT ------------------------------------------------------------

# Create period summary
period_summary <- kenya_df %>%
  filter(!is.na(period)) %>%
  group_by(period) %>%
  summarise(
    n_records = n(),
    n_species = n_distinct(species),
    n_genera = n_distinct(genus),
    n_families = n_distinct(family),
    .groups = "drop"
  )

print(period_summary)
saveRDS(period_summary, file.path(data_outputs, "temporal_period_summary.rds"))

temporal_bias_summary <- list(
  temporal_summary = temporal_summary,  # Year-by-year data (year, n_records, n_species, etc.)
  temporal_stats = temporal_stats,      # Overall statistics (year_min, year_max, means, etc.)
  trend_tests = trend_tests,
  gap_summary = gap_summary,
  seasonal_summary = seasonal_summary,
  species_temporal_summary = species_temporal_summary,
  period_breaks = period_breaks,
  period_summary = period_summary,
  temporal_assessment_summary = temporal_assessment_summary,
  taxonomic_levels_analyzed = taxonomic_levels
)

saveRDS(temporal_bias_summary, file.path(data_outputs, "temporal_bias_summary.rds"))

message("\n=== Temporal bias assessment complete ===")
message("Results saved to: ", data_outputs)
message("Figures saved to: ", figures_dir)

# Copy results to results folder -----------------------------------------------
source(here("scripts", "results_config.R"))
copy_to_results(TEMPORAL_BIAS_FILES, from_dir = data_outputs, to_subdir = "temporal_bias")
