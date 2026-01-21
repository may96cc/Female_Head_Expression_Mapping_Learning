#!/usr/bin/env Rscript

library(DSPRqtlDataB)
library(data.table)

# Get all datasets in the B package
all_datasets <- data(package = "DSPRqtlDataB")$results[, "Item"]

# Initialize list for data.tables
dt_list <- vector("list", length(all_datasets))
names(dt_list) <- all_datasets

for(i in seq_along(all_datasets)) {
  ds_name <- all_datasets[i]
  
  # Load dataset
  data(list = ds_name, package = "DSPRqtlDataB")
  df <- get(ds_name)
  
  # Keep only ril + AA1–AA8
  df_sub <- df[, c("ril", paste0("BB", 1:8))]
  
  # Rename columns: chrom_position_AAi
  base_name <- sub("^B_", "", ds_name)
  setnames(df_sub, old = paste0("BB", 1:8), new = paste0(base_name, "_BB", 1:8))
  
  # Convert to data.table
  dt_list[[i]] <- as.data.table(df_sub)
}

# Merge all datasets by ril efficiently
merged_dt <- Reduce(function(x, y) merge(x, y, by = "ril", all = TRUE), dt_list)

# Save as CSV with requested name
fwrite(merged_dt, "B_RILS_Hap_Probs.csv")

