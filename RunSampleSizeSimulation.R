# Load necessary libraries
library(DeclareDesign)
library(dplyr)
library(tidyr) # For pivot_wider if complex reshaping is needed, though not in current plan
library(ggplot2)

# Source the script that defines the CreateSampleDesign function
# This also ensures LoadAndProcessData.R is sourced, making site_props etc. available
source("Scripts/DefineDesign.R")

# --- Simulation Parameters ---
# Define a vector of sample sizes to test
sample_sizes_to_evaluate <- c(250, 500, 1000, 2000, 4000, 8000)
# Define the number of simulations to run for each sample size in diagnose_design
num_simulations_per_n <- 500 # Increase for more stable estimates, 500 is moderate

cat(paste("Starting sample size simulation for N values:", paste(sample_sizes_to_evaluate, collapse=", "), "\n"))
cat(paste("Number of simulations per N:", num_simulations_per_n, "\n\n"))

# Initialize an empty list to store diagnostics summaries from each N
all_diagnostics_summary <- list()

# --- Loop through Sample Sizes ---
for (N_val in sample_sizes_to_evaluate) {
  cat(paste("Processing N =", N_val, "...\n"))

  set.seed(123 + N_val)

  # Create the design for the current sample size, setting outcome_scenario to 1
  current_design <- CreateSampleDesign(current_N = N_val, outcome_scenario = 1)

  diagnostics_output <- DeclareDesign::diagnose_design(
    design = current_design,
    sims = num_simulations_per_n
  )

  # Extract relevant information from the diagnostics object
  # The column for estimator label in diagnose_design output is 'estimator_label'.
  # The column for estimand label is 'estimand_label'.
  summary_for_N <- diagnostics_output |>
    dplyr::filter(estimand_label %in% c("T1", "T2", "T1:T2") &
                  estimator_label == "OLS_clustered" & # UPDATED to OLS_clustered
                  diagnosand %in% c("power", "bias", "rmse", "coverage", "mean_estimate")) |>
    dplyr::select(estimand_label, diagnosand, estimate) |>
    dplyr::mutate(N = N_val)

  all_diagnostics_summary[[length(all_diagnostics_summary) + 1]] <- summary_for_N

  cat(paste("Finished processing N =", N_val, "\n\n"))
}

# Combine the list of data frames into a single data frame
final_summary_df <- dplyr::bind_rows(all_diagnostics_summary)

# --- Print and Plot Results ---
cat("Summary of Diagnostics Across Sample Sizes:\n")
# Ensure the full tibble is printed if it's long
print(final_summary_df, n = nrow(final_summary_df))


# Plot Power vs. N
power_plot_data <- final_summary_df |> dplyr::filter(diagnosand == "power")
power_plot <- ggplot2::ggplot(power_plot_data,
                              ggplot2::aes(x = N, y = estimate, color = estimand_label, group = estimand_label)) +
  ggplot2::geom_line(linewidth = 1) +
  ggplot2::geom_point(size = 2) +
  ggplot2::labs(title = "Power vs. Sample Size (Scenario 1, Clustered SEs)",
                x = "Sample Size (N persons)",
                y = "Power",
                color = "Estimand") +
  ggplot2::theme_minimal(base_size = 14) +
  ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  ggplot2::scale_x_continuous(breaks = sample_sizes_to_evaluate, labels = sample_sizes_to_evaluate)

print(power_plot)
ggsave("power_vs_N_scenario1_clustered.png", plot = power_plot, width = 8, height = 6)


# Plot Bias vs. N
bias_plot_data <- final_summary_df |> dplyr::filter(diagnosand == "bias")
bias_plot <- ggplot2::ggplot(bias_plot_data,
                             ggplot2::aes(x = N, y = estimate, color = estimand_label, group = estimand_label)) +
  ggplot2::geom_line(linewidth = 1) +
  ggplot2::geom_point(size = 2) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  ggplot2::labs(title = "Bias vs. Sample Size (Scenario 1, Clustered SEs)",
                x = "Sample Size (N persons)",
                y = "Bias",
                color = "Estimand") +
  ggplot2::theme_minimal(base_size = 14) +
  ggplot2::scale_x_continuous(breaks = sample_sizes_to_evaluate, labels = sample_sizes_to_evaluate)

print(bias_plot)
ggsave("bias_vs_N_scenario1_clustered.png", plot = bias_plot, width = 8, height = 6)


# Plot RMSE vs. N
rmse_plot_data <- final_summary_df |> dplyr::filter(diagnosand == "rmse")
rmse_plot <- ggplot2::ggplot(rmse_plot_data,
                             ggplot2::aes(x = N, y = estimate, color = estimand_label, group = estimand_label)) +
  ggplot2::geom_line(linewidth = 1) +
  ggplot2::geom_point(size = 2) +
  ggplot2::labs(title = "RMSE vs. Sample Size (Scenario 1, Clustered SEs)",
                x = "Sample Size (N persons)",
                y = "RMSE",
                color = "Estimand") +
  ggplot2::theme_minimal(base_size = 14) +
  ggplot2::scale_x_continuous(breaks = sample_sizes_to_evaluate, labels = sample_sizes_to_evaluate) +
  ggplot2::scale_y_continuous(limits = c(0, NA))

print(rmse_plot)
ggsave("rmse_vs_N_scenario1_clustered.png", plot = rmse_plot, width = 8, height = 6)

cat("\nSimulations complete. Plots saved to power_vs_N_scenario1_clustered.png, bias_vs_N_scenario1_clustered.png, and rmse_vs_N_scenario1_clustered.png.\n")
