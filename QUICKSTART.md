# Quick Start Guide: GBIF Kenya Bias Assessment

This guide will help you get started with the analysis in under 15 minutes.

## Prerequisites

1. **R** (version ≥ 4.2.0) - Download from [CRAN](https://cran.r-project.org/)
2. **RStudio** (recommended) - Download from [RStudio](https://posit.co/download/rstudio-desktop/)
3. **GBIF account** - Register at [GBIF.org](https://www.gbif.org/)

## Step 1: Clone Repository

```bash
git clone https://github.com/username/Gbif-Kenya.git
cd Gbif-Kenya
```

## Step 2: Install R Packages

Open R/RStudio and run:

```r
# Set working directory to project folder
setwd("path/to/Gbif-Kenya")

# Install all required packages
source("install_packages.R")
```

This will install approximately 40 packages. Takes ~5-10 minutes depending on your system.

## Step 3: Configure GBIF Credentials

**IMPORTANT**: You need GBIF credentials to download data.

1. Create a free account at [GBIF.org](https://www.gbif.org/)
2. Set environment variables:

```r
# Edit your R environment file
usethis::edit_r_environ()

# Add these lines (replace with YOUR credentials):
GBIF_USER="your_username"
GBIF_PWD="your_password"
GBIF_EMAIL="your_email@example.com"

# Save the file and restart R
```

**Alternative method** (temporary session):

```r
Sys.setenv(GBIF_USER = "your_username")
Sys.setenv(GBIF_PWD = "your_password")
Sys.setenv(GBIF_EMAIL = "your_email@example.com")
```

## Step 4: Run the Analysis

### Option A: Run Everything at Once

```r
source("run_all_analyses.R")
```

This will:
- Download GBIF data for Kenya (~30-60 min for download)
- Run spatial bias analysis (~10 min)
- Run temporal bias analysis (~5 min)
- Run taxonomic bias analysis (~10 min)
- Run statistical modeling (~15 min)
- Compile the manuscript (~5 min)

**Total time**: 1-2 hours (mostly waiting for GBIF download)

### Option B: Run Step-by-Step

```r
# Step 1: Download data (takes longest)
source("scripts/01_data_download.R")

# Step 2: Spatial bias
source("scripts/02_spatial_bias.R")

# Step 3: Temporal bias
source("scripts/03_temporal_bias.R")

# Step 4: Taxonomic bias
source("scripts/04_taxonomic_bias.R")

# Step 5: Statistical models
source("scripts/05_statistical_models.R")

# Step 6: Compile manuscript
rmarkdown::render("docs/kenya_gbif_bias_assessment.Rmd")
```

## Step 5: View Results

After analysis completes:

### Figures

```r
# List all generated figures
list.files("figures", pattern = "\\.png$")

# View a figure
system("open figures/01_sampling_effort_map.png")  # macOS
# OR
shell.exec("figures/01_sampling_effort_map.png")   # Windows
```

### Data Outputs

```r
# Load processed data
kenya_data <- readRDS("data/processed/kenya_gbif_clean.rds")

# View summary
summary_stats <- readRDS("data/processed/summary_stats.rds")
print(summary_stats)

# Load analysis results
spatial_summary <- readRDS("data/outputs/spatial_bias_summary.rds")
print(spatial_summary)
```

### Manuscript

Open the compiled manuscript:

```r
# Open HTML version
browseURL("docs/kenya_gbif_bias_assessment.html")
```

Or find it in the `docs/` folder.

## Troubleshooting

### Issue: GBIF download fails

**Solution**:
- Check your GBIF credentials are correct
- Verify internet connection
- Try downloading with an existing download key:

```r
# In scripts/01_data_download.R, uncomment and use existing key:
download_key <- "0123456-250101120000000"  # Example key
kenya_data <- occ_download_get(download_key)
```

### Issue: Package installation fails

**Solution**:
- Update R to latest version
- Install system dependencies (Linux):

```bash
# Ubuntu/Debian
sudo apt-get install libgdal-dev libproj-dev libgeos-dev libudunits2-dev

# macOS
brew install gdal proj geos udunits
```

### Issue: Memory errors

**Solution**:
- Increase R memory limit:

```r
# Increase memory (Windows)
memory.limit(size = 16000)  # 16 GB
```

- Use a subset of data for testing:

```r
# In scripts, add after loading data:
kenya_data <- kenya_data %>% slice_sample(n = 10000)
```

### Issue: Spatial packages fail to load

**Solution**:
- Install missing system libraries
- Try installing packages individually:

```r
install.packages("sf", type = "binary")  # Use pre-compiled version
```

## Working with the Data

### Explore the data

```r
library(tidyverse)

# Load cleaned data
kenya_data <- readRDS("data/processed/kenya_gbif_clean.rds")

# Basic exploration
glimpse(kenya_data)
summary(kenya_data)

# Top species
kenya_data %>%
  count(species, sort = TRUE) %>%
  head(20)

# Top classes
kenya_data %>%
  count(class, sort = TRUE)

# Temporal distribution
kenya_data %>%
  count(year) %>%
  ggplot(aes(x = year, y = n)) +
  geom_line() +
  labs(title = "Records per Year")
```

### Customize analyses

Modify parameters in scripts:

```r
# Change grid resolution (02_spatial_bias.R)
kenya_grid <- st_make_grid(
  kenya_boundary,
  cellsize = 0.05,  # 5km instead of 10km
  square = FALSE
)

# Change time period (01_data_download.R)
kenya_download <- occ_download(
  pred("country", "KE"),
  pred("hasCoordinate", TRUE),
  pred_gte("year", 2000),  # Only data from 2000 onwards
  # ...
)
```

## Next Steps

### Adapt for Another Country

1. Change country code in `01_data_download.R`:
   ```r
   pred("country", "TZ")  # Tanzania
   ```

2. Update boundary and environmental data downloads

3. Run the pipeline

### Add New Analyses

1. Create new script in `scripts/` directory
2. Add to `run_all_analyses.R`
3. Update manuscript to include results

### Share Your Work

1. Push to GitHub
2. Create DOI on Zenodo
3. Submit manuscript
4. Share data outputs

## Resources

- **GBIF API Documentation**: https://www.gbif.org/developer/summary
- **occAssess Package**: https://docs.ropensci.org/occAssess/
- **R Spatial**: https://r-spatial.org/
- **R for Data Science**: https://r4ds.had.co.nz/

## Getting Help

1. Check the main [README.md](README.md)
2. Review script comments
3. Open an issue on GitHub
4. Email: [your-email@example.com]

## Citation

If you use this workflow, please cite:

```
Kwiz Computing Technologies (2025). Spatial, Temporal, and Taxonomic Bias Assessment of
GBIF Biodiversity Data in Kenya. GitHub: https://github.com/username/Gbif-Kenya

GBIF data citation:
GBIF.org (Date) GBIF Occurrence Download https://doi.org/10.15468/dl.XXXXXX
```

---

**Happy analyzing!** 🦁🦒🐘

Last updated: November 10, 2025
