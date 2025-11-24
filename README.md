# Comprehensive Data Quality Assessment and Cleaning: A Case Study Using GBIF Biodiversity Data from Kenya

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

This repository presents a **comprehensive, reproducible framework for data quality assessment and cleaning**, demonstrated through application to biodiversity occurrence data. Using Global Biodiversity Information Facility (GBIF) data from Kenya as a case study, this framework provides transparent documentation of quality control procedures, systematic bias assessment, and gap analysis methods that can be adapted to any large-scale observational dataset.

The framework emphasizes:
- **Transparent data cleaning workflows** with systematic documentation of all quality filters and their impacts
- **Quantitative bias assessment** across spatial, temporal, and taxonomic dimensions
- **Reproducible methods** that can be applied to other datasets, regions, or research domains
- **Educational design** suitable for teaching data quality principles and reproducible research practices

While demonstrated using GBIF biodiversity data from Kenya, the methods, visualizations, and quality control procedures implemented here serve as a general template for assessing and improving data quality in any large observational dataset with spatial, temporal, and categorical components.

## Key Features

### Data Quality Framework
- **Systematic quality tracking**: Comprehensive documentation and quantification of data quality issues at each cleaning step
- **Transparent quality reporting**: Detailed breakdown of records affected by each filter with counts and percentages
- **Modular quality checks**: Independent tests for coordinate validity, taxonomic completeness, temporal consistency, and duplicates
- **Sensitivity analysis**: Comparison of different filtering strategies to assess robustness

### Bias Assessment Methods
- **Multi-dimensional bias analysis**: Spatial, temporal, and taxonomic dimensions assessed independently
- **Statistical modeling**: GLMs and GAMs to identify systematic predictors of sampling patterns
- **Gap identification**: Quantitative methods to detect undersampled regions, time periods, and taxa
- **Accessibility bias quantification**: Analysis of infrastructure and environmental influences on sampling

### Reproducibility & Transparency
- **Complete function flow documentation**: Every output traced back to source code with line numbers
- **Reproducible workflow**: All analyses fully reproducible using documented R scripts
- **Interactive visualizations**: Dynamic tables, charts, and maps for hands-on exploration
- **Educational design**: Suitable for teaching data quality principles, bias assessment, and reproducible research

### Adaptability & Reusability
- **Framework approach**: Methods designed for adaptation to other datasets, regions, or domains
- **Well-documented code**: Modular functions with clear parameters and comments
- **Customization guides**: Instructions for adapting to different countries, grid resolutions, or taxonomic groups
- **Open science**: All code, data, and documentation openly available under MIT license

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

## Key Results from the Kenya Case Study

Applying the framework to GBIF data from Kenya reveals typical patterns found in large observational datasets:

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

This framework is designed for teaching fundamental data science concepts:

1. **Data Quality Principles**: Demonstrates systematic approaches to identifying and addressing data quality issues
2. **Reproducible Research**: Students can trace every result back to raw data and source code
3. **Bias Assessment Methods**: Shows how to quantify and visualize systematic biases in observational data
4. **Spatial Statistics**: Examples of spatial autocorrelation, clustering analysis, and environmental niche modeling
5. **Data Visualization**: Both static (publication-ready) and interactive (exploratory) visualization techniques
6. **Statistical Modeling**: GLMs and GAMs for identifying predictors and quantifying relationships
7. **R Programming**: Well-documented, modular code with reusable functions

**Beyond Biodiversity**: While using biodiversity data, the methods apply broadly to any observational dataset:
- Public health surveillance data (spatial, temporal, and demographic patterns)
- Environmental monitoring networks (sensor placement bias, temporal gaps)
- Social science surveys (geographic sampling bias, response patterns)
- Citizen science projects (participation patterns, quality control)

**For Educators**: The interactive features allow students to:
- Explore patterns hands-on rather than passively viewing
- Test "what-if" scenarios by filtering and sorting data
- Understand the impact of different quality control decisions
- See exactly which code generated which results
- Learn transferable skills applicable to any observational dataset

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

If you use this framework, code, or methods, please cite:

```
Kwiz Computing Technologies (2025). Comprehensive Data Quality Assessment and Cleaning:
A Case Study Using GBIF Biodiversity Data from Kenya. GitHub repository:
https://github.com/username/Gbif-Kenya

GBIF data citation (for Kenya case study):
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
