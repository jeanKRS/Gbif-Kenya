# Comprehensive Data Quality Assessment and Bias Analysis of GBIF Data in Kenya

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

This repository contains a comprehensive, reproducible analysis of **data quality issues** and **systematic biases** in Global Biodiversity Information Facility (GBIF) occurrence data for Kenya. The analysis provides transparent documentation of all quality filters applied during data cleaning, quantifying the number and percentage of records affected by each issue. It also uses the `occAssess` R package and modern statistical methods to quantify spatial, temporal, and taxonomic biases and identify gaps in biodiversity data coverage.

## Key Features

- **Comprehensive data quality tracking**: Systematic documentation and quantification of all quality issues (missing coordinates, coordinate uncertainty, duplicates, taxonomic completeness, etc.)
- **Transparent quality reporting**: Detailed breakdown of records removed at each filtering step with percentages relative to original dataset
- **CoordinateCleaner integration**: Seven independent coordinate quality tests with individual issue quantification
- **Interactive visualizations**: Dynamic tables, charts, and maps for hands-on exploration (educational)
- **Function flow documentation**: Complete traceability from raw data to final results
- **Reproducible workflow**: All analyses are fully reproducible using R scripts and R Markdown
- **Multi-dimensional bias assessment**: Spatial, temporal, and taxonomic dimensions analyzed
- **Statistical modeling**: GLMs and GAMs to identify predictors of sampling effort
- **Publication-ready outputs**: Figures, tables, and compiled manuscript
- **Educational design**: Perfect for teaching reproducible research and biodiversity informatics
- **Open science**: All code, data, and documentation openly available

## Project Structure

```
Gbif-Kenya/
├── scripts/               # Analysis scripts (run in order)
│   ├── 01_data_download.R            # Download and clean GBIF data
│   ├── 01b_data_quality_assessment.R # Comprehensive quality tracking & metrics
│   ├── 02_spatial_bias.R             # Spatial bias assessment
│   ├── 03_temporal_bias.R            # Temporal bias assessment
│   ├── 04_taxonomic_bias.R           # Taxonomic bias assessment
│   └── 05_statistical_models.R       # GLM/GAM modeling
├── R/                     # Utility functions
│   └── utils.R                       # Helper functions
├── data/                  # Data directory
│   ├── raw/                          # Raw GBIF downloads
│   ├── processed/                    # Cleaned data
│   └── outputs/                      # Analysis outputs & quality metrics
├── figures/               # Generated figures
│   ├── data_quality_cascade.png      # Quality filtering cascade
│   ├── data_quality_retention.png    # Cumulative retention plot
│   ├── coordinate_issues_breakdown.png # Coordinate quality issues
│   └── ...                           # Bias assessment figures
├── docs/                  # Documentation and manuscripts
│   ├── kenya_gbif_bias_assessment.Rmd  # Main manuscript
│   └── references.bib                   # Bibliography
└── README.md             # This file
```

## Requirements

### Software

- R (≥ 4.2.0)
- RStudio (recommended)
- Git

### R Packages

Main packages required:

- **Data acquisition**: `rgbif`, `geodata`, `rnaturalearth`
- **Data cleaning**: `CoordinateCleaner`, `tidyverse`, `sf`
- **Bias assessment**: `occAssess`, `vegan`, `iNEXT`
- **Statistical analysis**: `MASS`, `lme4`, `mgcv`, `DHARMa`, `performance`
- **Visualization**: `ggplot2`, `viridis`, `patchwork`, `treemap`
- **Interactive visualization**: `DT`, `plotly`, `leaflet`, `htmltools` (NEW!)
- **Reproducibility**: `here`, `renv`, `rmarkdown`, `knitr`

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/username/Gbif-Kenya.git
cd Gbif-Kenya
```

### 2. Set up R environment

Open R/RStudio in the project directory and install required packages:

```r
# Install renv for package management
install.packages("renv")

# Restore package environment (if renv.lock exists)
renv::restore()

# Or install packages manually
install.packages(c("tidyverse", "sf", "here", "rgbif", "occAssess",
                   "CoordinateCleaner", "rnaturalearth", "rnaturalearthdata",
                   "terra", "geodata", "vegan", "iNEXT", "MASS", "lme4",
                   "mgcv", "DHARMa", "performance", "viridis", "patchwork",
                   "Kendall", "lubridate", "scales", "knitr", "rmarkdown",
                   "treemap", "ggalluvial", "ineq", "spdep", "ape", "osmdata",
                   "forecast", "tseries", "effects", "MuMIn", "kableExtra",
                   "DT", "plotly", "leaflet", "htmltools"))  # Interactive viz packages
```

### 3. Set GBIF credentials

Create a GBIF account at [https://www.gbif.org/](https://www.gbif.org/) if you don't have one.

Set environment variables for GBIF authentication:

```r
# Edit R environment file
usethis::edit_r_environ()

# Add these lines (replace with your credentials):
GBIF_USER="your_username"
GBIF_PWD="your_password"
GBIF_EMAIL="your_email@example.com"

# Restart R session
```

## Usage

### Quick Start

Run the complete analysis pipeline:

```r
# Set working directory to project root
setwd("path/to/Gbif-Kenya")

# Run all analysis scripts in sequence
source("scripts/01_data_download.R")
source("scripts/01b_data_quality_assessment.R")  # NEW: Quality tracking
source("scripts/02_spatial_bias.R")
source("scripts/03_temporal_bias.R")
source("scripts/04_taxonomic_bias.R")
source("scripts/05_statistical_models.R")

# Compile the manuscript
rmarkdown::render("docs/kenya_gbif_bias_assessment.Rmd")
```

### Step-by-Step Workflow

#### Step 1: Data Download and Cleaning

```r
source("scripts/01_data_download.R")
```

This script:
- Downloads Kenya occurrence data from GBIF
- Applies quality filters (coordinates, uncertainty, dates)
- Performs coordinate cleaning using CoordinateCleaner
- Saves cleaned data to `data/processed/`

**Note**: GBIF downloads can take time (minutes to hours). The script will wait for download completion.

#### Step 1b: Data Quality Assessment (NEW)

```r
source("scripts/01b_data_quality_assessment.R")
```

This script:
- Systematically tracks each quality filter applied during cleaning
- Quantifies records removed at each step (counts and percentages)
- Runs individual CoordinateCleaner tests to quantify specific coordinate issues
- Calculates comprehensive quality metrics (taxonomic completeness, uncertainty distributions, duplicate rates)
- Generates quality visualizations:
  - Filtering cascade plot
  - Cumulative retention plot
  - Coordinate issues breakdown
  - Summary pie chart

Outputs:
- `data/outputs/data_quality_assessment.rds` - Comprehensive quality metrics
- `data/outputs/data_quality_tracking.csv` - Step-by-step filtering results
- `data/outputs/coordinate_issues_summary.csv` - CoordinateCleaner test results
- `figures/data_quality_*.png` - Quality visualization figures

#### Step 2: Spatial Bias Assessment

```r
source("scripts/02_spatial_bias.R")
```

Analyzes:
- Grid-based sampling effort and species richness
- Spatial autocorrelation (Moran's I)
- Environmental bias (sampled vs. available space)
- occAssess spatial metrics

Outputs:
- Spatial grid with metrics
- Maps of sampling effort and richness
- Environmental bias plots

#### Step 3: Temporal Bias Assessment

```r
source("scripts/03_temporal_bias.R")
```

Analyzes:
- Temporal trends (Mann-Kendall tests)
- Temporal completeness and gaps
- Seasonal patterns
- Species temporal coverage

Outputs:
- Temporal trend plots
- Seasonal distribution plots
- Gap analysis results

#### Step 4: Taxonomic Bias Assessment

```r
source("scripts/04_taxonomic_bias.R")
```

Analyzes:
- Taxonomic coverage across hierarchical levels
- Diversity indices (Gini, Simpson, Shannon)
- Species accumulation curves
- Rarity patterns

Outputs:
- Taxonomic treemaps
- Lorenz curves
- Rarity distribution plots

#### Step 5: Statistical Modeling

```r
source("scripts/05_statistical_models.R")
```

Fits models to identify predictors of sampling effort:
- Negative Binomial GLM for record counts
- Binomial GLM for sampling presence
- GAMs for non-linear relationships

Outputs:
- Model summaries and diagnostics
- Effect plots
- Sampling bias maps

#### Step 6: Compile Manuscript

```r
rmarkdown::render("docs/kenya_gbif_bias_assessment.Rmd",
                  output_format = "html_document") # "all"
```

Generates publication-ready documents in HTML, PDF, and Word formats.

## Key Results

Our analysis of GBIF data for Kenya reveals:

### Spatial Bias
- **Strong clustering**: Moran's I indicates significant spatial autocorrelation
- **Accessibility bias**: Sampling concentrated near cities and roads
- **Environmental bias**: Sampled areas differ from available environmental space

### Temporal Bias
- **Increasing trend**: Data collection has increased over time
- **High variability**: Large coefficient of variation in annual records
- **Temporal gaps**: Many years with no or few records

### Taxonomic Bias
- **Uneven coverage**: High Gini coefficient indicates concentration in few groups
- **Charismatic bias**: Birds and mammals over-represented
- **Rarity issues**: High proportion of singleton species

### Predictors of Sampling
- **Distance to cities**: Strongest predictor (negative effect)
- **Elevation**: Non-linear relationship
- **Environmental factors**: Moderate effects

## Interactive Features & Educational Design

### NEW: Interactive Visualizations

The manuscript now includes **interactive tables, charts, and maps** designed for educational purposes and hands-on exploration:

#### Interactive Tables (DT)
- **Click to sort** any column (ascending/descending)
- **Search and filter** data in real-time
- **Export** selections to CSV or Excel
- **Pagination** for browsing large datasets

**Example**: Sensitivity Analysis table allows sorting by retention percentage to compare filtering strategies

#### Interactive Maps (Leaflet)
- **Zoom and pan** to explore geographic patterns at any scale
- **Click cells** to see detailed statistics (records, species, coordinates)
- **Hover** to highlight specific regions
- **Basemap** provides geographic context

**Example**: Sampling effort map lets you zoom into Nairobi (high sampling) vs. northern Kenya (sampling gaps)

#### Interactive Charts (Plotly)
- **Zoom** into specific time periods or value ranges
- **Hover** for exact values at any data point
- **Toggle** data series on/off via legend
- **Download** charts as PNG for presentations

**Example**: Temporal trends chart allows zooming into specific decades to see annual variation

### Function Flow Documentation

Every output is annotated with **complete traceability**:

- **Source file** and line numbers showing where code lives
- **Function names** indicating which transformations were applied
- **Collapsible code blocks** revealing the exact R code
- **Educational notes** explaining the analysis pipeline

See **`docs/FUNCTION_FLOW.md`** for comprehensive documentation mapping every table and figure to its source code.

### Educational Use Cases

This project is designed for teaching:

1. **Reproducible Research**: Students can trace results back to raw data
2. **Biodiversity Informatics**: Demonstrates data quality assessment workflows
3. **Spatial Statistics**: Shows spatial autocorrelation and bias analysis
4. **Data Visualization**: Examples of static and interactive visualizations
5. **R Programming**: Well-documented code with reusable functions

**For Educators**: The interactive features allow students to:
- Explore patterns hands-on rather than passively viewing
- Test "what-if" scenarios by filtering and sorting
- Understand the impact of quality control decisions
- See exactly which code generated which results

**Function Flow Reference**: `docs/FUNCTION_FLOW.md` provides a complete learning pathway from beginner to advanced, with code explanations and educational notes.

## Outputs

All analysis outputs are saved to:

- **Processed data**: `data/processed/`
- **Analysis results**: `data/outputs/`
- **Figures**: `figures/` (PNG, 300 dpi)
- **Manuscript**: `docs/` (HTML, PDF, Word)

## Customization

### Adapting for Other Countries

To adapt this workflow for another country:

1. Modify `01_data_download.R`:
   ```r
   # Change country code (ISO 2-letter)
   kenya_download <- occ_download(
     pred("country", "TZ"),  # Tanzania
     # ... rest of predicates
   )
   ```

2. Update country boundaries:
   ```r
   country_sf <- ne_countries(country = "Tanzania", returnclass = "sf")
   ```

3. Adjust environmental data extraction:
   ```r
   elevation <- elevation_30s(country = "TZA", path = tempdir())
   ```

### Modifying Grid Resolution

Change grid cell size in spatial analyses:

```r
# In 02_spatial_bias.R
kenya_grid <- st_make_grid(
  kenya_boundary,
  cellsize = 0.05,  # ~5km instead of 10km
  square = FALSE
)
```

## Citation

If you use this code or analysis framework, please cite:

```
[Authors] (2025). Spatial, Temporal, and Taxonomic Bias Assessment of GBIF
Biodiversity Data in Kenya. GitHub repository:
https://github.com/username/Gbif-Kenya

GBIF data citation:
GBIF.org (DD Month YYYY) GBIF Occurrence Download https://doi.org/10.15468/dl.XXXXXX
```

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/improvement`)
3. Commit changes (`git commit -am 'Add improvement'`)
4. Push to branch (`git push origin feature/improvement`)
5. Create Pull Request

## License

This project is licensed under the MIT License - see LICENSE file for details.

## Acknowledgments

- **GBIF** and all data contributors for open biodiversity data
- **occAssess** developers for the bias assessment framework
- **R community** for developing open-source statistical tools
- **Funding agencies**: [Add your funders]

## Contact

For questions or issues:
- Open an issue on GitHub
- Email: [your-email@example.com]

## References

Key references:

- Marsh et al. (2023). occAssess: An R package for assessing potential biases in species occurrence data. *Ecography*, 2023(1), e06299.
- Beck et al. (2014). Spatial bias in the GBIF database and its effect on modeling species' geographic distributions. *Ecological Informatics*, 19, 10-15.
- Hortal et al. (2015). Seven shortfalls that beset large-scale knowledge of biodiversity. *Annual Review of Ecology, Evolution, and Systematics*, 46, 523-549.

See `docs/references.bib` for complete bibliography.

## Version History

- **v1.0.0** (2025-11-10): Initial release
  - Complete spatial, temporal, and taxonomic bias assessment
  - Statistical modeling of sampling predictors
  - Reproducible R Markdown manuscript

---

**Last Updated**: November 10, 2025
