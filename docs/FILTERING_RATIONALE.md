# Data Quality Filtering Rationale by Analysis

This document explains which quality flags each analysis applies and why, based on the specific requirements of each bias assessment.

## Available Quality Flags

The data quality assessment (`01b_data_quality_assessment.R`) adds the following flags to each record:

### Core Flags
- `flag_missing_coords` - Missing decimal latitude or longitude
- `flag_high_uncertainty` - Coordinate uncertainty > 10 km
- `flag_inappropriate_basis` - Fossil or living specimens
- `flag_missing_species` - No species-level identification
- `flag_invalid_date` - Date before 1950 or in the future
- `flag_duplicate` - Exact duplicate (same species, location, date)

### CoordinateCleaner Flags
- `flag_coord_capitals` - Within 10 km of capitals
- `flag_coord_centroids` - Within 5 km of country/province centroids
- `flag_coord_equal` - Identical latitude and longitude values
- `flag_coord_gbif` - At GBIF headquarters location
- `flag_coord_zeros` - At (0,0) or near equator/prime meridian
- `flag_coord_urban` - In major urban areas
- `flag_coord_outliers` - Statistical outliers for the species
- `flag_any_coord_issue` - Any of the above coordinate issues

## Filtering by Analysis

### 02_spatial_bias.R - Spatial Bias Assessment

**Filters Applied:**
```r
kenya_clean <- kenya_flagged %>%
  filter(
    !flag_missing_coords,      # REQUIRED: Need coordinates for spatial analysis
    !flag_missing_species,      # REQUIRED: Need species for diversity metrics
    !flag_invalid_date,         # RECOMMENDED: For temporal comparisons
    !flag_coord_equal,          # STRICT: Clearly erroneous coordinates
    !flag_coord_gbif,           # STRICT: Clearly erroneous coordinates
    !flag_coord_zeros           # STRICT: Clearly erroneous coordinates
  )
```

**Rationale:**
- **Missing coords**: Impossible to analyze spatial patterns without coordinates
- **Missing species**: Needed for species richness and diversity metrics
- **Invalid dates**: Temporal trends in spatial bias require valid dates
- **Equal/GBIF/Zeros coords**: These are clear data errors that would create artificial hotspots

**Flags NOT Applied:**
- `flag_high_uncertainty`: High uncertainty doesn't invalidate the location, just reduces precision
- `flag_coord_capitals/centroids/urban`: These ARE the biases we're trying to quantify! Filtering them out would remove the signal we want to detect
- `flag_coord_outliers`: Outliers might represent real range expansions or rare habitats
- `flag_duplicate`: Duplicates indicate sampling intensity, relevant for bias assessment
- `flag_inappropriate_basis`: Living/fossil specimens still have spatial information

**Impact Documentation:**
Each analysis should report:
- Number of records before filtering
- Number removed by each flag
- Percentage of records retained
- Comparison to original dataset

---

### 03_temporal_bias.R - Temporal Bias Assessment

**Filters Applied:**
```r
kenya_clean <- kenya_flagged %>%
  filter(
    !flag_invalid_date,         # REQUIRED: Need valid dates for temporal analysis
    !flag_missing_species       # REQUIRED: Need species for diversity trends
  )
```

**Rationale:**
- **Invalid dates**: Core requirement for temporal analysis
- **Missing species**: Needed for taxonomic trends over time

**Flags NOT Applied:**
- `flag_missing_coords`: Temporal patterns don't require coordinates
- `flag_coord_*`: Coordinate quality doesn't affect temporal patterns
- `flag_high_uncertainty`: Doesn't affect date validity
- `flag_duplicate`: Multiple records on same date indicate sampling effort
- `flag_inappropriate_basis`: Fossils have temporal value; living specimens show collection trends

**Impact Documentation:**
Report temporal coverage changes after filtering and how it affects trend detection.

---

### 04_taxonomic_bias.R - Taxonomic Bias Assessment

**Filters Applied:**
```r
kenya_clean <- kenya_flagged %>%
  filter(
    !flag_missing_species       # REQUIRED: Need species for taxonomic analysis
  )
```

**Rationale:**
- **Missing species**: Absolutely required for taxonomic assessments

**Flags NOT Applied:**
- `flag_missing_coords`: Taxonomic coverage doesn't require spatial data
- `flag_invalid_date`: Taxonomic representation is relevant across all time periods
- `flag_coord_*`: Coordinate quality doesn't affect taxonomic identification
- `flag_high_uncertainty`: Doesn't affect taxonomy
- `flag_duplicate`: Duplicates show which species are well-sampled
- `flag_inappropriate_basis`: Fossils and living specimens contribute to taxonomic coverage

**Impact Documentation:**
Report how filtering affects representation across taxonomic groups.

---

### 05_statistical_models.R - Statistical Modeling

**Filters Applied:**
```r
kenya_clean <- kenya_flagged %>%
  filter(
    !flag_missing_coords,       # REQUIRED: Need coords for environmental extraction
    !flag_missing_species,      # REQUIRED: For species-specific models
    !flag_invalid_date,         # RECOMMENDED: For temporal covariates
    !flag_coord_equal,          # STRICT: Clear errors
    !flag_coord_gbif,           # STRICT: Clear errors
    !flag_coord_zeros,          # STRICT: Clear errors
    !flag_coord_outliers        # STRICT: Would bias environmental relationships
  )
```

**Rationale:**
- **Missing coords**: Need coordinates to extract environmental predictors
- **Missing species**: Required for species-level modeling
- **Invalid dates**: Temporal covariates require valid dates
- **Clear coord errors**: Would introduce noise in environmental relationships
- **Outliers**: Geographic outliers could drive spurious environmental correlations

**Flags NOT Applied:**
- `flag_high_uncertainty`: Can be included as a covariate in models
- `flag_coord_capitals/centroids/urban`: These biases are what we're modeling!
- `flag_duplicate`: Could be weighted or included as a covariate

**Additional Considerations:**
- Consider using `flag_high_uncertainty` as a model weight or covariate
- `flag_coord_urban` and `flag_coord_capitals` should be **predictors**, not filters
- Document how filtering affects model fit and coefficient estimates

---

## Sensitivity Analysis

Each analysis should include sensitivity checks:

1. **Filter Sensitivity**: Re-run with more/less strict filters and compare results
2. **Flag as Covariate**: Include flags as model predictors instead of filters
3. **Stratified Analysis**: Analyze flagged vs. unflagged records separately

## Documentation Template

Each analysis should document filtering decisions:

```r
# Applied Filters
cat("\\n=== Data Quality Filtering ===\\n")
cat("Original records:", nrow(kenya_flagged), "\\n")
cat("Filters applied:\\n")
cat("  - Missing coordinates:", sum(kenya_flagged$flag_missing_coords), "removed\\n")
cat("  - Missing species:", sum(kenya_flagged$flag_missing_species), "removed\\n")
# ... etc

cat("\\nRecords retained:", nrow(kenya_clean),
    "(", round(100 * nrow(kenya_clean) / nrow(kenya_flagged), 1), "%)\\n")

cat("\\nRationale: [Explain why these specific filters]\\n")
```

## Flag Combinations

Some analyses might want specific combinations:

```r
# Example: Moderate filtering for robust analysis
moderate_clean <- kenya_flagged %>%
  filter(
    !flag_missing_coords,
    !flag_missing_species,
    !flag_invalid_date,
    !flag_any_coord_issue  # Any coordinate quality issue
  )

# Example: Minimal filtering to preserve maximum data
minimal_clean <- kenya_flagged %>%
  filter(
    !flag_missing_coords,
    !flag_missing_species
  )

# Example: Strict filtering for high-quality subset
strict_clean <- kenya_flagged %>%
  filter(
    !flag_missing_coords,
    !flag_high_uncertainty,
    !flag_inappropriate_basis,
    !flag_missing_species,
    !flag_invalid_date,
    !flag_duplicate,
    !flag_any_coord_issue
  )
```

## Key Principle

**The golden rule**: Only filter out flags that represent:
1. **Data errors** (impossible values, clear mistakes)
2. **Missing required information** for that specific analysis
3. **Noise** that would obscure the pattern you're studying

**Do NOT filter**:
- Biases you're trying to measure
- Issues that can be modeled as covariates
- Problems that don't affect your specific analysis

Different analyses have different requirements - there is no "one size fits all" cleaning strategy!
