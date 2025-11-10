# ==============================================================================
# Script: 04_taxonomic_bias.R
# Purpose: Assess taxonomic biases in GBIF data for Kenya
# Author: Automated Research Pipeline
# Date: 2025-11-10
# ==============================================================================

# Load required packages -------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(occAssess)
  library(vegan)
  library(iNEXT)
  library(treemapify)
  library(viridis)
  library(ggalluvial)
  library(scales)
})

# Setup paths ------------------------------------------------------------------
data_processed <- here("data", "processed")
data_outputs <- here("data", "outputs")
figures_dir <- here("figures")

# Load cleaned data ------------------------------------------------------------
message("Loading cleaned GBIF data...")
kenya_data <- readRDS(file.path(data_processed, "kenya_gbif_clean.rds"))

# 1. TAXONOMIC COVERAGE OVERVIEW -----------------------------------------------
message("\n=== Analyzing taxonomic coverage ===")

# Overall taxonomic summary
tax_overview <- kenya_data %>%
  summarise(
    n_records = n(),
    n_kingdoms = n_distinct(kingdom, na.rm = TRUE),
    n_phyla = n_distinct(phylum, na.rm = TRUE),
    n_classes = n_distinct(class, na.rm = TRUE),
    n_orders = n_distinct(order, na.rm = TRUE),
    n_families = n_distinct(family, na.rm = TRUE),
    n_genera = n_distinct(genus, na.rm = TRUE),
    n_species = n_distinct(species, na.rm = TRUE)
  )

print(tax_overview)
saveRDS(tax_overview, file.path(data_outputs, "taxonomic_overview.rds"))

# 2. HIERARCHICAL TAXONOMIC SUMMARY --------------------------------------------
message("\n=== Creating hierarchical taxonomic summaries ===")

# Kingdom level
kingdom_summary <- kenya_data %>%
  filter(!is.na(kingdom)) %>%
  group_by(kingdom) %>%
  summarise(
    n_records = n(),
    n_phyla = n_distinct(phylum, na.rm = TRUE),
    n_classes = n_distinct(class, na.rm = TRUE),
    n_orders = n_distinct(order, na.rm = TRUE),
    n_families = n_distinct(family, na.rm = TRUE),
    n_species = n_distinct(species, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_records))

# Phylum level
phylum_summary <- kenya_data %>%
  filter(!is.na(phylum)) %>%
  group_by(kingdom, phylum) %>%
  summarise(
    n_records = n(),
    n_classes = n_distinct(class, na.rm = TRUE),
    n_orders = n_distinct(order, na.rm = TRUE),
    n_families = n_distinct(family, na.rm = TRUE),
    n_species = n_distinct(species, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_records))

# Class level
class_summary <- kenya_data %>%
  filter(!is.na(class)) %>%
  group_by(kingdom, phylum, class) %>%
  summarise(
    n_records = n(),
    n_orders = n_distinct(order, na.rm = TRUE),
    n_families = n_distinct(family, na.rm = TRUE),
    n_species = n_distinct(species, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_records))

# Order level
order_summary <- kenya_data %>%
  filter(!is.na(order)) %>%
  group_by(kingdom, phylum, class, order) %>%
  summarise(
    n_records = n(),
    n_families = n_distinct(family, na.rm = TRUE),
    n_species = n_distinct(species, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_records))

# Save summaries
saveRDS(kingdom_summary, file.path(data_outputs, "kingdom_summary.rds"))
saveRDS(phylum_summary, file.path(data_outputs, "phylum_summary.rds"))
saveRDS(class_summary, file.path(data_outputs, "class_summary.rds"))
saveRDS(order_summary, file.path(data_outputs, "order_summary.rds"))

# Print top results
message("\nTop 10 Classes:")
print(head(class_summary, 10))

# 3. TAXONOMIC BIAS METRICS ----------------------------------------------------
message("\n=== Calculating taxonomic bias metrics ===")

# Calculate proportion of records per major group
tax_proportions <- class_summary %>%
  mutate(
    prop_records = n_records / sum(n_records),
    prop_species = n_species / sum(n_species),
    cumsum_records = cumsum(prop_records),
    cumsum_species = cumsum(prop_species)
  )

# Calculate Lorenz curve and Gini coefficient for records
gini_records <- ineq::Gini(class_summary$n_records)
gini_species <- ineq::Gini(class_summary$n_species)

# Simpson's diversity index
simpson_class <- diversity(class_summary$n_records, index = "simpson")
shannon_class <- diversity(class_summary$n_records, index = "shannon")

# Pielou's evenness
evenness_class <- shannon_class / log(nrow(class_summary))

diversity_metrics <- data.frame(
  metric = c("Gini_records", "Gini_species", "Simpson_diversity",
             "Shannon_diversity", "Pielou_evenness"),
  value = c(gini_records, gini_species, simpson_class, shannon_class, evenness_class),
  interpretation = c(
    paste0("Records concentration (0=equal, 1=concentrated): ", round(gini_records, 3)),
    paste0("Species concentration: ", round(gini_species, 3)),
    paste0("Class diversity: ", round(simpson_class, 3)),
    paste0("Class diversity: ", round(shannon_class, 3)),
    paste0("Evenness (0=uneven, 1=even): ", round(evenness_class, 3))
  )
)

print(diversity_metrics)
saveRDS(diversity_metrics, file.path(data_outputs, "taxonomic_diversity_metrics.rds"))

# 4. SPECIES ACCUMULATION CURVES -----------------------------------------------
message("\n=== Computing species accumulation curves ===")

# Overall species accumulation
set.seed(123)

# Create species-by-site matrix (sites = 1km grid cells)
species_matrix <- kenya_data %>%
  mutate(
    cell_id = paste0(
      floor(decimalLongitude * 10), "_",
      floor(decimalLatitude * 10)
    )
  ) %>%
  group_by(cell_id, species) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = species, values_from = n, values_fill = 0) %>%
  column_to_rownames("cell_id")

# Rarefaction curve
spec_accum <- specaccum(species_matrix, method = "random", permutations = 100)

# Save accumulation data
accum_data <- data.frame(
  sites = spec_accum$sites,
  richness = spec_accum$richness,
  sd = spec_accum$sd
)
saveRDS(accum_data, file.path(data_outputs, "species_accumulation.rds"))

# iNEXT analysis for coverage-based rarefaction
# Prepare frequency data (only first 100 species for computational efficiency)
freq_data <- kenya_data %>%
  count(species, sort = TRUE) %>%
  head(100) %>%
  pull(n)

inext_result <- iNEXT(freq_data, q = 0, datatype = "abundance", endpoint = max(freq_data) * 2)
saveRDS(inext_result, file.path(data_outputs, "inext_rarefaction.rds"))

# 5. TAXONOMIC REPRESENTATION BY SPATIAL COVERAGE ------------------------------
message("\n=== Analyzing taxonomic spatial representation ===")

# Calculate spatial extent and sampling intensity for each class
class_spatial <- kenya_data %>%
  filter(!is.na(class)) %>%
  group_by(class) %>%
  summarise(
    n_records = n(),
    n_species = n_distinct(species),
    n_locations = n_distinct(paste(decimalLongitude, decimalLatitude)),
    lon_range = max(decimalLongitude) - min(decimalLongitude),
    lat_range = max(decimalLatitude) - min(decimalLatitude),
    spatial_extent = lon_range * lat_range,
    .groups = "drop"
  ) %>%
  mutate(
    records_per_location = n_records / n_locations,
    species_per_location = n_species / n_locations
  ) %>%
  arrange(desc(n_records))

saveRDS(class_spatial, file.path(data_outputs, "class_spatial_coverage.rds"))

# 6. RARE VS COMMON SPECIES ANALYSIS -------------------------------------------
message("\n=== Analyzing rarity patterns ===")

# Species frequency distribution
species_freq <- kenya_data %>%
  count(species, name = "n_records") %>%
  mutate(
    rarity_class = case_when(
      n_records == 1 ~ "Singleton (1 record)",
      n_records == 2 ~ "Doubleton (2 records)",
      n_records <= 5 ~ "Very rare (3-5 records)",
      n_records <= 10 ~ "Rare (6-10 records)",
      n_records <= 50 ~ "Uncommon (11-50 records)",
      n_records <= 100 ~ "Common (51-100 records)",
      TRUE ~ "Very common (>100 records)"
    )
  )

rarity_summary <- species_freq %>%
  group_by(rarity_class) %>%
  summarise(
    n_species = n(),
    total_records = sum(n_records),
    .groups = "drop"
  ) %>%
  mutate(
    prop_species = n_species / sum(n_species),
    prop_records = total_records / sum(total_records)
  )

print(rarity_summary)
saveRDS(rarity_summary, file.path(data_outputs, "rarity_summary.rds"))

# Identify most recorded species
top_species <- species_freq %>%
  arrange(desc(n_records)) %>%
  head(50)

saveRDS(top_species, file.path(data_outputs, "top_50_species.rds"))

# 7. TEMPORAL CHANGES IN TAXONOMIC COMPOSITION ---------------------------------
message("\n=== Analyzing temporal changes in taxonomic composition ===")

# Create 10-year period breaks
year_range <- range(kenya_data$year, na.rm = TRUE)
period_breaks <- seq(
  floor(year_range[1] / 10) * 10,
  ceiling(year_range[2] / 10) * 10,
  by = 10
)

# Add period labels to data
kenya_data_periods <- kenya_data %>%
  mutate(
    period = cut(year,
                breaks = period_breaks,
                include.lowest = TRUE,
                right = FALSE,
                labels = paste0(head(period_breaks, -1), "-",
                              tail(period_breaks, -1) - 1))
  )

# Taxonomic composition by period
period_tax_summary <- kenya_data_periods %>%
  filter(!is.na(period), !is.na(class)) %>%
  group_by(period, class) %>%
  summarise(
    n_records = n(),
    n_species = n_distinct(species),
    .groups = "drop"
  )

# Get top 10 classes overall
top_10_classes <- kenya_data %>%
  count(class, sort = TRUE) %>%
  head(10) %>%
  pull(class)

# Filter period summary for top classes
period_tax_top <- period_tax_summary %>%
  filter(class %in% top_10_classes)

saveRDS(period_tax_summary, file.path(data_outputs, "taxonomic_by_period.rds"))
saveRDS(period_tax_top, file.path(data_outputs, "taxonomic_top10_by_period.rds"))

# Calculate diversity metrics by period
period_diversity <- kenya_data_periods %>%
  filter(!is.na(period)) %>%
  group_by(period) %>%
  summarise(
    n_records = n(),
    n_species = n_distinct(species),
    n_genera = n_distinct(genus),
    n_families = n_distinct(family),
    n_orders = n_distinct(order),
    n_classes = n_distinct(class),
    .groups = "drop"
  )

print(period_diversity)
saveRDS(period_diversity, file.path(data_outputs, "diversity_by_period.rds"))

# 8. OCCASSESS TAXONOMIC ASSESSMENTS -------------------------------------------
message("\n=== Running occAssess taxonomic assessments ===")

# Prepare data for occAssess with taxonomic levels
kenya_df <- kenya_data %>%
  select(species, genus, family, order, class, phylum, kingdom,
         decimalLongitude, decimalLatitude, eventDate, year) %>%
  filter(!is.na(species)) %>%
  as.data.frame()

# Add period information
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

# Check if taxonomic assessments already exist
taxonomic_assessments_exist <- file.exists(file.path(data_outputs, "occassess_species_identity_by_taxlevel.rds")) &&
                               file.exists(file.path(data_outputs, "taxonomic_assessment_summary.rds"))

if (taxonomic_assessments_exist) {
  message("  ℹ Taxonomic occAssess assessments already completed. Loading from files...")
  species_identity_assessments <- readRDS(file.path(data_outputs, "occassess_species_identity_by_taxlevel.rds"))
  taxonomic_assessment_summary <- readRDS(file.path(data_outputs, "taxonomic_assessment_summary.rds"))
  print(taxonomic_assessment_summary)
} else {
  message("  → Running taxonomic occAssess assessments...")

  # Species identity assessment - Run for each taxonomic level
  message("\n--- Assessing species identity by taxonomic level ---")
  species_identity_assessments <- list()

  for (tax_level in taxonomic_levels) {
    message(sprintf("  Running assessSpeciesID for %s", tax_level))

    if (all(is.na(kenya_df[[tax_level]]))) {
      message(sprintf("  Skipping %s - all values are NA", tax_level))
      next
    }

    tryCatch({
      species_identity_assessments[[tax_level]] <- assessSpeciesID(
        dat = kenya_df,
        speciesCol = tax_level,
        assumptions = c("singletonSpecies", "highlyCommonSpecies"),
        singletonThreshold = 1,
        commonThreshold = 100
      )
      message(sprintf("  ✓ Completed identity assessment for %s", tax_level))
    }, error = function(e) {
      message(sprintf("  ✗ Error in identity assessment for %s: %s", tax_level, e$message))
    })
  }

  # Save occAssess results
  message("\n--- Saving occAssess taxonomic results ---")
  saveRDS(species_identity_assessments, file.path(data_outputs, "occassess_species_identity_by_taxlevel.rds"))

  # Also save species-level assessment for backward compatibility
  if (!is.null(species_identity_assessments[["species"]])) {
    saveRDS(species_identity_assessments[["species"]], file.path(data_outputs, "occassess_species_identity.rds"))
  }

  # Create summary of taxonomic assessments
  taxonomic_assessment_summary <- data.frame(
    taxonomic_level = taxonomic_levels,
    identity_assessment = sapply(taxonomic_levels, function(x) !is.null(species_identity_assessments[[x]]))
  )

  print(taxonomic_assessment_summary)
  saveRDS(taxonomic_assessment_summary, file.path(data_outputs, "taxonomic_assessment_summary.rds"))
  message("  ✓ Taxonomic occAssess assessments completed")
}

# 9. COMPARISON WITH EXPECTED DIVERSITY ----------------------------------------
message("\n=== Comparing with expected diversity (if reference data available) ===")

# Note: This section would compare with Kenya's known species lists
# For now, we document what we have vs potential sources

comparison_notes <- data.frame(
  taxonomic_group = c("Birds", "Mammals", "Reptiles", "Amphibians", "Plants"),
  gbif_species = c(
    sum(class_summary$n_species[class_summary$class == "Aves"]),
    sum(class_summary$n_species[class_summary$class == "Mammalia"]),
    sum(class_summary$n_species[class_summary$class == "Reptilia"]),
    sum(class_summary$n_species[class_summary$class == "Amphibia"]),
    sum(class_summary$n_species[class_summary$kingdom == "Plantae"], na.rm = TRUE)
  ),
  estimated_total_kenya = c(1136, 390, 200, 120, 7000),  # Approximate known totals
  reference = c(
    "Birdlife Kenya",
    "Kenya Wildlife Service",
    "National Museums of Kenya",
    "National Museums of Kenya",
    "Kenya Forest Service"
  )
) %>%
  mutate(
    coverage_percent = 100 * gbif_species / estimated_total_kenya,
    gap = estimated_total_kenya - gbif_species
  )

print(comparison_notes)
saveRDS(comparison_notes, file.path(data_outputs, "taxonomic_coverage_comparison.rds"))

# 10. VISUALIZATION -------------------------------------------------------------
message("\n=== Creating taxonomic visualizations ===")

# Plot 1: Treemap of taxonomic groups
top_classes <- head(class_summary, 15)

p1 <- ggplot(top_classes, aes(area = n_records, fill = n_species,
                              label = paste0(class, "\n",
                                           scales::comma(n_records), " records\n",
                                           n_species, " species"))) +
  geom_treemap() +
  geom_treemap_text(color = "white", place = "centre", size = 10) +
  scale_fill_viridis_c(option = "plasma", name = "Number of Species") +
  labs(title = "Taxonomic Representation in GBIF Kenya Data",
       subtitle = "Top 15 Classes by Number of Records") +
  theme_minimal()

ggsave(file.path(figures_dir, "10_taxonomic_treemap.png"),
       p1, width = 14, height = 10, dpi = 300)

# Plot 2: Lorenz curve (taxonomic inequality)
p2 <- ggplot(tax_proportions, aes(x = seq_along(cumsum_records) / n(),
                                  y = cumsum_records)) +
  geom_line(linewidth = 1.5, color = "blue") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  annotate("text", x = 0.7, y = 0.3,
           label = paste0("Gini coefficient = ", round(gini_records, 3)),
           size = 5) +
  labs(
    title = "Taxonomic Inequality: Lorenz Curve",
    subtitle = "Distribution of records across taxonomic classes",
    x = "Cumulative Proportion of Classes",
    y = "Cumulative Proportion of Records"
  ) +
  theme_minimal() +
  coord_fixed()

ggsave(file.path(figures_dir, "11_lorenz_curve.png"),
       p2, width = 8, height = 8, dpi = 300)

# Plot 3: Species accumulation curve
p3 <- ggplot(accum_data, aes(x = sites, y = richness)) +
  geom_ribbon(aes(ymin = richness - 2*sd, ymax = richness + 2*sd),
              alpha = 0.3, fill = "blue") +
  geom_line(linewidth = 1.5, color = "darkblue") +
  labs(
    title = "Species Accumulation Curve",
    subtitle = "Random accumulation with 95% confidence interval",
    x = "Number of Sites (1km grid cells)",
    y = "Cumulative Number of Species"
  ) +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = comma) +
  theme_minimal()

ggsave(file.path(figures_dir, "12_species_accumulation.png"),
       p3, width = 10, height = 6, dpi = 300)

# Plot 4: Rarity distribution
rarity_summary$rarity_class <- factor(
  rarity_summary$rarity_class,
  levels = c("Singleton (1 record)", "Doubleton (2 records)",
             "Very rare (3-5 records)", "Rare (6-10 records)",
             "Uncommon (11-50 records)", "Common (51-100 records)",
             "Very common (>100 records)")
)

p4 <- ggplot(rarity_summary, aes(x = rarity_class, y = n_species)) +
  geom_col(aes(fill = rarity_class), alpha = 0.8) +
  geom_text(aes(label = comma(n_species)), vjust = -0.5, size = 3) +
  scale_fill_viridis_d(option = "magma") +
  labs(
    title = "Distribution of Species by Rarity",
    subtitle = paste0("Total species: ", sum(rarity_summary$n_species)),
    x = "Rarity Class",
    y = "Number of Species"
  ) +
  scale_y_continuous(labels = comma) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

ggsave(file.path(figures_dir, "13_rarity_distribution.png"),
       p4, width = 12, height = 6, dpi = 300)

# Plot 5: Top taxa barplot
p5 <- ggplot(head(order_summary, 20), aes(x = reorder(order, n_records), y = n_records)) +
  geom_col(aes(fill = class), alpha = 0.8) +
  coord_flip() +
  scale_fill_viridis_d(option = "turbo") +
  scale_y_continuous(labels = comma) +
  labs(
    title = "Top 20 Orders by Number of Records",
    x = "Order",
    y = "Number of Records",
    fill = "Class"
  ) +
  theme_minimal()

ggsave(file.path(figures_dir, "14_top_orders.png"),
       p5, width = 12, height = 8, dpi = 300)

# Plot 6: Records vs Species scatter
p6 <- ggplot(head(class_summary, 30),
             aes(x = n_records, y = n_species, label = class)) +
  geom_point(aes(size = n_families, color = kingdom), alpha = 0.7) +
  geom_text(vjust = -1, size = 2.5, check_overlap = TRUE) +
  scale_x_log10(labels = comma) +
  scale_y_log10(labels = comma) +
  scale_size_continuous(name = "Number of Families") +
  scale_color_viridis_d() +
  labs(
    title = "Records vs Species Richness by Class",
    subtitle = "Top 30 classes (log scale)",
    x = "Number of Records (log scale)",
    y = "Number of Species (log scale)",
    color = "Kingdom"
  ) +
  theme_minimal()

ggsave(file.path(figures_dir, "15_records_vs_species.png"),
       p6, width = 12, height = 8, dpi = 300)

# Plot 7: Taxonomic composition changes over time periods
p7 <- ggplot(period_tax_top, aes(x = period, y = n_records, fill = class)) +
  geom_col(position = "fill", alpha = 0.8) +
  scale_fill_viridis_d(option = "turbo") +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "Taxonomic Composition Changes Over Time",
    subtitle = "Proportion of records by class across 10-year periods (Top 10 classes)",
    x = "Time Period",
    y = "Proportion of Records",
    fill = "Class"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

ggsave(file.path(figures_dir, "16_taxonomic_composition_over_time.png"),
       p7, width = 14, height = 8, dpi = 300)

# Plot 8: Diversity metrics by period
p8 <- ggplot(period_diversity, aes(x = period)) +
  geom_line(aes(y = n_species, group = 1, color = "Species"), linewidth = 1.2) +
  geom_line(aes(y = n_families * 10, group = 1, color = "Families (×10)"), linewidth = 1.2) +
  geom_point(aes(y = n_species, color = "Species"), size = 3) +
  geom_point(aes(y = n_families * 10, color = "Families (×10)"), size = 3) +
  scale_color_manual(values = c("Species" = "blue", "Families (×10)" = "red")) +
  scale_y_continuous(
    name = "Number of Species",
    sec.axis = sec_axis(~./10, name = "Number of Families")
  ) +
  labs(
    title = "Taxonomic Diversity by Time Period",
    subtitle = "Species and family richness across 10-year periods",
    x = "Time Period",
    color = ""
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    axis.title.y.right = element_text(color = "red"),
    axis.text.y.right = element_text(color = "red"),
    axis.title.y.left = element_text(color = "blue"),
    axis.text.y.left = element_text(color = "blue")
  )

ggsave(file.path(figures_dir, "17_diversity_by_period.png"),
       p8, width = 12, height = 6, dpi = 300)

# 11. SUMMARY REPORT -----------------------------------------------------------
taxonomic_summary <- list(
  overview = tax_overview,
  diversity_metrics = diversity_metrics,
  rarity_summary = rarity_summary,
  top_classes = head(class_summary, 20),
  coverage_comparison = comparison_notes,
  period_breaks = period_breaks,
  period_diversity = period_diversity,
  taxonomic_assessment_summary = taxonomic_assessment_summary,
  taxonomic_levels_analyzed = taxonomic_levels
)

saveRDS(taxonomic_summary, file.path(data_outputs, "taxonomic_bias_summary.rds"))

message("\n=== Taxonomic bias assessment complete ===")
message("Results saved to: ", data_outputs)
message("Figures saved to: ", figures_dir)
