# Spatial, Temporal, and Taxonomic Bias Assessment of GBIF Data in Kenya

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

This repository contains a comprehensive, reproducible analysis of spatial, temporal, and taxonomic biases in Global Biodiversity Information Facility (GBIF) occurrence data for Kenya. The analysis uses the `occAssess` R package and modern statistical methods to quantify biases and identify gaps in biodiversity data coverage.

## Key Features

- **Reproducible workflow**: All analyses are fully reproducible using R scripts and R Markdown
- **Comprehensive bias assessment**: Spatial, temporal, and taxonomic dimensions analyzed
- **Statistical modeling**: GLMs and GAMs to identify predictors of sampling effort
- **Publication-ready outputs**: Figures, tables, and compiled manuscript
- **Open science**: All code, data, and documentation openly available

## Project Structure

```
Gbif-Kenya/
├── scripts/               # Analysis scripts (run in order)
│   ├── 01_data_download.R           # Download and clean GBIF data
│   ├── 02_spatial_bias.R            # Spatial bias assessment
│   ├── 03_spatial_bias.R           # Temporal bias assessment
│   ├── 04_taxonomic_bias.R          # Taxonomic bias assessment
│   └── 05_statistical_models.R      # GLM/GAM modeling
├── R/                     # Utility functions
│   └── utils.R                      # Helper functions
├── data/                  # Data directory
│   ├── raw/                         # Raw GBIF downloads
│   ├── processed/                   # Cleaned data
│   └── outputs/                     # Analysis outputs
├── figures/               # Generated figures
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
                   "forecast", "tseries", "effects", "MuMIn", "kableExtra"))
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
