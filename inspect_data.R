# inspect_data.R

# Load necessary libraries
library(tidyverse)

cat("=== Verification Script Output ===\n")

# Find the latest .rds file in Results
# This assumes the results are in a directory structure like Results/full-<timestamp>/
# and the RDS file is study_data_list.rds within that.

results_base_dir <- file.path("Results")
# List all subdirectories in Results
all_result_dirs <- list.dirs(path = results_base_dir, full.names = TRUE, recursive = FALSE)

# Filter for directories that likely match our settings type (e.g., "full")
# This is a heuristic; a more robust way might involve parsing settings.type from RunSimulation.R
# or having RunSimulation.R output the exact path.
# For now, let's assume "full" is part of the directory name.
# And pick the most recent one.
relevant_dirs <- grep("full-", all_result_dirs, value = TRUE)

if (length(relevant_dirs) == 0) {
  stop("No relevant result directories found in Results/. Looked for 'full-*'. Current dirs: ",
       paste(all_result_dirs, collapse=", "))
}

# Get the most recently created directory
dir_info <- file.info(relevant_dirs)
latest_dir <- relevant_dirs[which.max(dir_info$mtime)]

rds_file_path <- file.path(latest_dir, "study_data_list.rds")

if (!file.exists(rds_file_path)) {
  stop(paste("RDS file not found at:", rds_file_path))
}

cat(paste("\nLoading RDS file from:", rds_file_path, "\n"))
study_list <- readRDS(rds_file_path)

# Extract the first (and only, due to sim.reps=1) study data
study_df <- study_list[[1]]

cat("\n--- Structure of a few rows: ---\n")
print(head(study_df))
print(str(study_df))


# Verify Cluster Assignment
cat("\n--- Cluster Assignment Check: ---\n")
print(table(study_df$ClusterID, useNA = "ifany"))
cat("Proportions:\n")
print(prop.table(table(study_df$ClusterID, useNA = "ifany")))

# For by-cluster summaries, it's better to use the wide format (one row per patient)
# The current study_df is long (multiple rows per patient if multiple visits)
# We need to get unique patient data for covariates.
# Most covariates (ClusterID, B1, U1, N1) are patient-level.
patient_df <- study_df %>%
  distinct(ID, .keep_all = TRUE) # Keep first row for each patient for patient-level covars

cat("\n--- Patient-level data (first row per ID) structure: ---\n")
print(head(patient_df))
print(str(patient_df))

cat("\nNumber of unique patients for covariate checks:", nrow(patient_df), "\n")
cat("Cluster assignment for unique patients:\n")
print(table(patient_df$ClusterID, useNA = "ifany"))


# Verify Binary Covariate B1 Distribution
cat("\n--- Binary Covariate B1 Proportions by Cluster (Expected C1: ~0.35, C2: ~0.6): ---\n")
# Use patient_df for this as B1 is a patient-level covariate
if ("B1" %in% names(patient_df) && "ClusterID" %in% names(patient_df)) {
  proportions_b1 <- patient_df %>%
    group_by(ClusterID) %>%
    summarise(
      count_0 = sum(B1 == 0, na.rm = TRUE),
      count_1 = sum(B1 == 1, na.rm = TRUE),
      prop_1 = mean(B1, na.rm = TRUE),
      .groups = 'drop'
    )
  print(proportions_b1)
} else {
  cat("B1 or ClusterID column not found in patient_df.\n")
}

# Verify Continuous Covariate U1 Mean
cat("\n--- Continuous Covariate U1 Mean by Cluster (Expected C1: ~0, C2: ~0.2): ---\n")
if ("U1" %in% names(patient_df) && "ClusterID" %in% names(patient_df)) {
  means_u1 <- patient_df %>%
    group_by(ClusterID) %>%
    summarise(
      mean_U1 = mean(U1, na.rm = TRUE),
      n = n(),
      .groups = 'drop'
    )
  print(means_u1)
} else {
  cat("U1 or ClusterID column not found in patient_df.\n")
}

# Verify Continuous Covariate N1 SD
cat("\n--- Continuous Covariate N1 SD by Cluster (Expected C1: ~1, C2: ~1.2): ---\n")
if ("N1" %in% names(patient_df) && "ClusterID" %in% names(patient_df)) {
  sds_n1 <- patient_df %>%
    group_by(ClusterID) %>%
    summarise(
      sd_N1 = sd(N1, na.rm = TRUE),
      n = n(),
      .groups = 'drop'
    )
  print(sds_n1)
} else {
  cat("N1 or ClusterID column not found in patient_df.\n")
}

# Verify Outcome Linear Predictor Muij Mean (this is visit-level, so use study_df)
cat("\n--- Mean of Muij by Cluster (should differ if true.params differ): ---\n")
if ("Muij" %in% names(study_df) && "ClusterID" %in% names(study_df)) {
  means_muij <- study_df %>%
    group_by(ClusterID) %>%
    summarise(
      mean_Muij = mean(Muij, na.rm = TRUE),
      n_visits = n(), # Number of visits in this cluster
      .groups = 'drop'
    )
  print(means_muij)
} else {
  cat("Muij or ClusterID column not found in study_df (long format).\n")
}

cat("\n=== End of Verification Script Output ===\n")
