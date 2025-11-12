# Kenya GBIF Data Cleaning Optimization Guide

## Problem Summary

The original data cleaning pipeline encountered:
- **Memory Error**: Vector memory limit of 16.0 GB reached
- **Date Parsing Failures**: 204,534 records failed to parse with `ymd()`

## Optimizations Implemented

### 1. Use data.table for Memory Efficiency

**Before:**
```r
kenya_clean <- kenya_raw %>%
  filter(...) %>%
  mutate(...)
```

**After:**
```r
library(data.table)
kenya_dt <- as.data.table(kenya_raw)
kenya_dt <- kenya_dt[condition]  # In-place operations
```

**Benefits:**
- data.table uses 50-70% less memory than dplyr for large datasets
- In-place operations avoid creating intermediate copies
- Faster execution (5-10x speedup on large datasets)

### 2. Filter Early and Often

**Strategy**: Apply most restrictive filters first to reduce dataset size before expensive operations.

**Order of operations:**
1. Filter missing coordinates (removes many records)
2. Filter missing species (removes many records)
3. Filter basis of record
4. Filter coordinate uncertainty
5. Filter invalid years
6. Remove duplicates

Each filter is followed by `gc()` to free memory immediately.

### 3. Smart Date Parsing

**Problem**: Parsing 200K+ invalid dates wastes memory and time.

**Solution**:
- Use existing `year`, `month`, `day` columns from GBIF data when available
- Only parse `eventDate` for records with missing years
- Use `suppressWarnings()` and `quiet = TRUE` to avoid console spam
- Filter invalid years after parsing, not during

**Benefits:**
- Avoids parsing ~200K invalid dates
- Uses pre-computed GBIF year field when available
- Reduces memory allocation from failed parsing attempts

### 4. Memory Management

Added explicit memory management:
```r
rm(kenya_raw)  # Remove original after conversion
gc()           # Garbage collection after each major operation
```

### 5. Efficient Duplicate Removal

**Before:**
```r
distinct(species, decimalLongitude, decimalLatitude, eventDate, .keep_all = TRUE)
```

**After:**
```r
unique(kenya_dt, by = c("species", "decimalLongitude", "decimalLatitude", "event_date"))
```

data.table's `unique()` is faster and more memory-efficient for large datasets.

## Expected Performance Improvements

| Metric | Before | After | Improvement |
|--------|--------|-------|-------------|
| Memory Usage | ~16 GB (crashed) | ~8-10 GB (est.) | ~40-50% reduction |
| Processing Time | Unknown (crashed) | Faster | 5-10x speedup |
| Date Parsing | All records | Only records with missing year | ~80-90% fewer operations |

## Additional Optimization Strategies

If you still encounter memory issues, try these:

### Option 1: Increase Memory Limit (if available)
```r
# Increase memory limit to 32GB (if system has RAM available)
memory.limit(size = 32000)  # Windows only
# On Linux/Mac, set via: ulimit -v
```

### Option 2: Process in Chunks
```r
# Process data in chunks of 500,000 records
chunk_size <- 500000
n_chunks <- ceiling(nrow(kenya_raw) / chunk_size)

kenya_clean_list <- list()
for (i in 1:n_chunks) {
  start_idx <- (i - 1) * chunk_size + 1
  end_idx <- min(i * chunk_size, nrow(kenya_raw))

  chunk <- kenya_raw[start_idx:end_idx, ]
  chunk_clean <- clean_data_chunk(chunk)  # Apply cleaning pipeline
  kenya_clean_list[[i]] <- chunk_clean

  rm(chunk)
  gc()
}

kenya_clean <- rbindlist(kenya_clean_list)
```

### Option 3: Select Only Required Columns Early
```r
# Keep only essential columns before processing
essential_cols <- c("species", "decimalLongitude", "decimalLatitude",
                    "eventDate", "year", "month", "day",
                    "coordinateUncertaintyInMeters", "basisOfRecord",
                    "countryCode", "genus", "family", "order",
                    "class", "phylum", "kingdom")

kenya_dt <- kenya_dt[, ..essential_cols]
```

### Option 4: Use Disk-Based Processing
```r
# Use disk.frame for datasets too large for memory
library(disk.frame)

kenya_df <- as.disk.frame(kenya_raw, outdir = "data/temp")
kenya_clean <- kenya_df %>%
  filter(...) %>%
  collect()
```

### Option 5: Optimize CoordinateCleaner Step
The CoordinateCleaner operations can also be memory-intensive:

```r
# Process coordinate cleaning in batches
batch_size <- 100000
flags <- logical(nrow(kenya_clean))

for (i in seq(1, nrow(kenya_clean), by = batch_size)) {
  end_i <- min(i + batch_size - 1, nrow(kenya_clean))
  batch_flags <- clean_coordinates(
    x = kenya_clean[i:end_i, ],
    # ... parameters ...
  )
  flags[i:end_i] <- batch_flags
  gc()
}
```

## Monitoring Memory Usage

Add memory monitoring to track usage:

```r
# Function to report memory usage
report_memory <- function(label) {
  mem_used <- pryr::mem_used()
  message(sprintf("%s: %.2f GB used", label, mem_used / 1e9))
}

# Use throughout pipeline
report_memory("After loading data")
# ... operations ...
report_memory("After cleaning")
```

## Best Practices for Large GBIF Datasets

1. **Download strategically**: Use GBIF download predicates to pre-filter data
2. **Avoid repeated downloads**: Cache cleaned data and work with subsets during development
3. **Test with subsets**: Develop pipeline on 10% sample, then run on full dataset
4. **Monitor resources**: Use `htop` (Linux/Mac) or Task Manager (Windows) to watch memory
5. **Use appropriate hardware**: Large GBIF downloads may require >32 GB RAM

## Testing the Optimizations

To test the optimized pipeline:

```r
# Run on a subset first
source("scripts/01_data_download.R")

# Monitor for memory issues and parse failures
# Check output messages for:
# - "Records with invalid/missing years" count
# - Memory usage stays below limit
# - No vector memory errors
```

## Citation

If these optimizations help your research, please cite:
- GBIF.org (for the data)
- data.table package: Dowle M, Srinivasan A (2024). data.table: Extension of `data.frame`
