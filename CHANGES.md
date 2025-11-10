# Analysis Updates - Taxonomic Levels and Temporal Periods

## Summary
Updated all bias assessment scripts to incorporate 10-year temporal periods and taxonomic hierarchy analysis across multiple taxonomic levels.

## Changes Made

### 1. Spatial Bias Analysis (`scripts/02_spatial_bias.R`)

**Added:**
- 10-year period breaks calculation based on year range
- Period labels for temporal grouping
- Taxonomic backbone extraction (kingdom → species)
- Multiple taxonomic levels as identifiers in `assessRecordNumber()`
- Taxonomic-level specific runs of:
  - `assessRecordNumber()` for 7 taxonomic levels + period grouping
  - `assessSpeciesNumber()` for all taxonomic levels
  - `assessSpeciesCoverage()` for all taxonomic levels
- Error handling with try-catch blocks for robustness
- Assessment summary tracking which taxonomic levels completed successfully
- Backward compatibility by saving species-level results separately
- Enhanced visualizations for each taxonomic level
- Taxonomic summary plot showing number of taxa at each level

**New Output Files:**
- `occassess_records_by_taxlevel.rds` - Record assessments for all taxonomic levels
- `occassess_species_by_taxlevel.rds` - Species number assessments by taxonomic level
- `occassess_coverage_by_taxlevel.rds` - Coverage assessments by taxonomic level
- `records_by_period.rds` - Summary of records across 10-year periods
- `assessment_summary.rds` - Tracking table of completed assessments
- Individual visualization PNGs for each taxonomic level

### 2. Temporal Bias Analysis (`scripts/03_temporal_bias.R`)

**Added:**
- 10-year period breaks (consistent with spatial analysis)
- Period labels integrated with data
- Taxonomic levels (kingdom → species) in data preparation
- Taxonomic-level specific runs of:
  - `assessRecordTime()` for all 7 taxonomic levels
  - `assessSpeciesTime()` for all 7 taxonomic levels
- Period-based summary statistics
- Temporal assessment summary tracking
- Enhanced summary report with period information

**New Output Files:**
- `occassess_record_time_by_taxlevel.rds` - Temporal distribution by taxonomic level
- `occassess_species_time_by_taxlevel.rds` - Temporal coverage by taxonomic level
- `temporal_period_summary.rds` - Records, species, genera, families by period
- `temporal_assessment_summary.rds` - Tracking table of temporal assessments

### 3. Taxonomic Bias Analysis (`scripts/04_taxonomic_bias.R`)

**Added:**
- 10-year period breaks calculation
- Period-based taxonomic composition analysis
- Taxonomic composition changes over time
- Diversity metrics by period (species, genera, families, orders, classes)
- Taxonomic-level specific runs of `assessSpeciesID()` for all 7 levels
- Visualizations:
  - Taxonomic composition over time (stacked proportions)
  - Diversity metrics by period (species and families)

**New Output Files:**
- `taxonomic_by_period.rds` - Full taxonomic composition by period
- `taxonomic_top10_by_period.rds` - Top 10 classes by period
- `diversity_by_period.rds` - Diversity metrics across periods
- `occassess_species_identity_by_taxlevel.rds` - Identity assessments by taxonomic level
- `taxonomic_assessment_summary.rds` - Tracking table

**New Visualizations:**
- `16_taxonomic_composition_over_time.png` - Stacked area chart
- `17_diversity_by_period.png` - Time series of diversity metrics

### 4. Statistical Models (`scripts/05_statistical_models.R`)

**Added:**
- Taxonomic-specific GLM models for top 5 taxonomic classes
- Environmental predictor effects by taxonomic group
- Comparative coefficient plot across taxonomic classes
- Separate negative binomial models for each major class

**New Output Files:**
- `taxonomic_specific_models.rds` - Models for each major taxonomic class
- `taxonomic_model_summaries.rds` - Model coefficients and statistics

**New Visualizations:**
- `24_taxonomic_coefficients_comparison.png` - Coefficient comparison across classes

## Key Features

### Temporal Periods
- **Period Breaks**: Automatically calculated based on data year range
- **Granularity**: 10-year intervals (e.g., 1950-1959, 1960-1969, etc.)
- **Consistency**: Same period breaks used across all scripts
- **Integration**: Period information included in all relevant assessments

### Taxonomic Levels
- **Hierarchy**: 7 levels analyzed (kingdom, phylum, class, order, family, genus, species)
- **Flexibility**: Each level used as identifier in occAssess functions
- **Robustness**: Automatic skipping of levels with all NA values
- **Error Handling**: Try-catch blocks prevent script failure if one level fails

### Assessment Functions Enhanced
All occAssess functions now run with multiple identifiers:
1. `assessRecordNumber()` - with periods parameter and 7 taxonomic levels + period
2. `assessSpeciesNumber()` - for all 7 taxonomic levels
3. `assessSpeciesCoverage()` - for all 7 taxonomic levels
4. `assessRecordTime()` - for all 7 taxonomic levels
5. `assessSpeciesTime()` - for all 7 taxonomic levels
6. `assessSpeciesID()` - for all 7 taxonomic levels

### Backward Compatibility
- Species-level results saved separately in original file names
- Existing downstream code will continue to work
- New taxonomic-level results available for enhanced analysis

## Data Structure
All analysis scripts now include:
```r
kenya_df <- kenya_data %>%
  select(species, genus, family, order, class, phylum, kingdom,
         decimalLongitude, decimalLatitude, coordinateUncertaintyInMeters,
         eventDate, year) %>%
  mutate(
    period = cut(year,
                breaks = period_breaks,
                include.lowest = TRUE,
                right = FALSE,
                labels = paste0(head(period_breaks, -1), "-",
                              tail(period_breaks, -1) - 1))
  )
```

## Benefits
1. **Temporal Insights**: Understand how sampling patterns change across decades
2. **Taxonomic Breadth**: Assess biases at multiple taxonomic resolutions
3. **Comparative Analysis**: Compare assessment metrics across taxonomic groups
4. **Robust Analysis**: Error handling ensures partial failures don't stop analysis
5. **Comprehensive Output**: All results tracked and summarized systematically

## Testing Recommendations
1. Run scripts sequentially: 01 → 02 → 03 → 04 → 05
2. Check output files for completeness
3. Review assessment summary files to see which taxonomic levels completed
4. Inspect new visualizations for quality and interpretability
5. Verify period breaks match expected year ranges

## Notes
- All scripts maintain backward compatibility
- Original file outputs preserved for existing workflows
- New outputs use `_by_taxlevel` suffix for clarity
- Period breaks automatically adjust to data year range
- Missing taxonomic data handled gracefully (automatic skip with message)
