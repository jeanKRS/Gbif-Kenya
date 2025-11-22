# Results Folder Structure

This folder contains the output files from each analysis in the GBIF Kenya Bias Assessment pipeline. The analyses are organized into subfolders by type, enabling efficient caching and selective re-running of analyses.

## How It Works

When you run `run_all_analyses.R`, the pipeline will:

1. Check if all required output files exist for each analysis
2. Skip analyses that have already been completed (all outputs present)
3. Run only analyses with missing outputs
4. Save new results to the appropriate subfolder

This approach significantly reduces computation time when re-running the pipeline, especially useful during iterative development or when only specific analyses need updating.

## Folder Structure

```
results/
├── data_quality/           # Data quality assessment outputs
├── spatial_bias/           # Spatial bias analysis outputs
├── temporal_bias/          # Temporal bias analysis outputs
├── taxonomic_bias/         # Taxonomic bias analysis outputs
└── statistical_models/     # Statistical modeling outputs
```

## Contents by Folder

### data_quality/

**Analysis Script:** `scripts/01b_data_quality_assessment.R`

**Key Outputs:**
- `data_quality_assessment.rds` - Complete quality assessment results
- `data_quality_tracking.csv` - Step-by-step filtering cascade
- `coordinate_issues_summary.csv` - CoordinateCleaner test results

**Purpose:** Comprehensive tracking of data quality issues including missing data, coordinate quality, duplicates, and spatial flags.

---

### spatial_bias/

**Analysis Script:** `scripts/02_spatial_bias.R`

**Key Outputs:**
- `spatial_grid_effort.rds` - Hexagonal grid with sampling effort
- `spatial_autocorrelation.rds` - Moran's I test results
- `environmental_comparison.rds` - Sampled vs. available environmental space
- `environmental_bias_tests.rds` - Kolmogorov-Smirnov test results
- `occassess_*.rds` - occAssess package assessments by taxonomic level
- `spatial_bias_summary.rds` - Overall spatial bias summary

**Purpose:** Quantify spatial sampling biases including clustering, accessibility effects, and environmental representation.

---

### temporal_bias/

**Analysis Script:** `scripts/03_temporal_bias.R`

**Key Outputs:**
- `temporal_statistics.rds` - Temporal coverage metrics
- `temporal_trends.rds` - Mann-Kendall trend test results
- `temporal_gaps.rds` - Years with missing data
- `monthly_patterns.rds` - Seasonal sampling patterns
- `seasonal_bias_test.rds` - Chi-square test for seasonality
- `decade_analysis.rds` - Records aggregated by decade
- `species_temporal_coverage.rds` - Species-level temporal spans
- `temporal_bias_summary.rds` - Overall temporal bias summary

**Purpose:** Assess temporal biases including trends, gaps, seasonality, and species-level temporal coverage.

---

### taxonomic_bias/

**Analysis Script:** `scripts/04_taxonomic_bias.R`

**Key Outputs:**
- `taxonomic_overview.rds` - Overall taxonomic summary
- `{kingdom,phylum,class,order}_summary.rds` - Hierarchical summaries
- `taxonomic_diversity_metrics.rds` - Gini, Simpson, Shannon indices
- `species_accumulation.rds` - Rarefaction curve data
- `rarity_summary.rds` - Singleton/doubleton analysis
- `diversity_by_period.rds` - Temporal changes in diversity
- `taxonomic_coverage_comparison.rds` - GBIF vs. expected diversity
- `taxonomic_bias_summary.rds` - Overall taxonomic bias summary

**Purpose:** Quantify taxonomic representation, diversity, rarity patterns, and coverage gaps.

---

### statistical_models/

**Analysis Script:** `scripts/05_statistical_models.R`

**Key Outputs:**
- `final_count_model.rds` - GLM/GAM for record counts
- `final_presence_model.rds` - Binomial GLM for sampling presence
- `model_summaries.rds` - Coefficient estimates and significance
- `model_selection_aic.rds` - AIC comparison of candidate models
- `grid_with_predictions.rds` - Predicted vs. observed sampling effort
- `modeling_summary.rds` - Overall modeling summary

**Optional Outputs (may not exist):**
- `gam_model.rds` - Generalized Additive Model (if fitted)
- `taxonomic_specific_models.rds` - Class-specific models
- `taxonomic_model_summaries.rds` - Taxonomic model coefficients

**Purpose:** Statistical models predicting sampling effort from environmental and geographic covariates, identifying over/under-sampled areas.

## Re-running Specific Analyses

To force re-run a specific analysis:

```r
# Delete the results folder for that analysis
unlink("results/spatial_bias", recursive = TRUE)

# Then run the full pipeline
source("run_all_analyses.R")
```

Or run individual scripts:

```r
source("scripts/02_spatial_bias.R")
```

## Backward Compatibility

All results are also saved to `data/outputs/` for backward compatibility with existing code and manuscripts. This dual-save approach ensures:

- New pipeline can cache and skip completed analyses
- Existing analysis code continues to work without modification
- RMarkdown documents can still find their expected data files

## Configuration

The list of required files for each analysis is defined in `scripts/results_config.R`. This file is used by `run_all_analyses.R` to determine if an analysis can be skipped.

To modify what files are required for an analysis, edit the corresponding list in `results_config.R`:

```r
SPATIAL_BIAS_FILES <- c(
  "spatial_grid_effort.rds",
  "spatial_autocorrelation.rds",
  # ... add more required files
)
```

## Utility Functions

The following utility functions in `R/utils.R` support the results caching system:

- `check_results_exist(results_dir, required_files)` - Check if all required files exist
- `skip_if_complete(analysis_name, results_dir, required_files)` - Skip analysis if complete
- `copy_to_results(files, from_dir, to_subdir)` - Copy files to results folder

## Tips

1. **Large files:** `.gitignore` is configured to exclude `*.rds` files, so results won't be committed to git
2. **Clean rebuild:** Delete the entire `results/` folder to force a complete re-analysis
3. **Partial rebuild:** Delete individual subfolders to re-run specific analyses
4. **Development mode:** When developing a new analysis, delete just that subfolder frequently

## Related Files

- `scripts/results_config.R` - Defines required output files
- `run_all_analyses.R` - Main pipeline with caching logic
- `data/outputs/` - Backward-compatible output location
- `figures/` - Visualization outputs (separate from cached data)
