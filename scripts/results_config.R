# ==============================================================================
# Results Configuration File
# Defines output files for each analysis to enable smart caching
# ==============================================================================

# Data Quality Assessment ------------------------------------------------------
DATA_QUALITY_FILES <- c(
  "data_quality_assessment.rds",
  "data_quality_tracking.csv",
  "coordinate_issues_summary.csv",
  "kenya_gbif_flagged.rds"  # Flagged data for downstream analyses
)

# Spatial Bias Analysis --------------------------------------------------------
SPATIAL_BIAS_FILES <- c(
  "spatial_grid_effort.rds",
  "spatial_autocorrelation.rds",
  "records_by_period.rds",
  "occassess_records_by_taxlevel.rds",
  "occassess_species_by_taxlevel.rds",
  "occassess_coverage_by_taxlevel.rds",
  "assessment_summary.rds",
  "environmental_comparison.rds",
  "environmental_bias_tests.rds",
  "spatial_bias_summary.rds"
)

# Temporal Bias Analysis -------------------------------------------------------
TEMPORAL_BIAS_FILES <- c(
  "temporal_statistics.rds",
  "temporal_trends.rds",
  "temporal_gaps.rds",
  "occassess_record_time_by_taxlevel.rds",
  "occassess_species_time_by_taxlevel.rds",
  "temporal_assessment_summary.rds",
  "decade_analysis.rds",
  "monthly_patterns.rds",
  "seasonal_bias_test.rds",
  "species_temporal_coverage.rds",
  "temporal_period_summary.rds",
  "temporal_bias_summary.rds"
)

# Taxonomic Bias Analysis ------------------------------------------------------
TAXONOMIC_BIAS_FILES <- c(
  "taxonomic_overview.rds",
  "kingdom_summary.rds",
  "phylum_summary.rds",
  "class_summary.rds",
  "order_summary.rds",
  "taxonomic_diversity_metrics.rds",
  "species_accumulation.rds",
  "inext_rarefaction.rds",
  "class_spatial_coverage.rds",
  "rarity_summary.rds",
  "top_50_species.rds",
  "taxonomic_by_period.rds",
  "taxonomic_top10_by_period.rds",
  "diversity_by_period.rds",
  "occassess_species_identity_by_taxlevel.rds",
  "taxonomic_assessment_summary.rds",
  "taxonomic_coverage_comparison.rds",
  "taxonomic_bias_summary.rds"
)

# Statistical Models -----------------------------------------------------------
STATISTICAL_MODELS_FILES <- c(
  "final_count_model.rds",
  "final_presence_model.rds",
  "model_summaries.rds",
  "model_selection_aic.rds",
  "grid_with_predictions.rds",
  "modeling_summary.rds"
)

# Optional files (GAM may not always be fitted)
STATISTICAL_MODELS_OPTIONAL <- c(
  "gam_model.rds",
  "taxonomic_specific_models.rds",
  "taxonomic_model_summaries.rds"
)
