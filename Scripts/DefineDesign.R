# Load necessary libraries
library(DeclareDesign)
library(fabricatr)
library(dplyr)
library(tidyr)
library(estimatr)

# Source data loading script for proportions
# This makes site_props, race_props_by_site, and sex_indication_props_by_site globally available
# for use within CreateSampleDesign
source("Scripts/LoadAndProcessData.R")

# Define a function to create the design, parameterized by sample size N
CreateSampleDesign <- function(current_N) {

  design <- DeclareDesign::declare_design(

    # == Population Step ==
    population = DeclareDesign::declare_population(
      N = current_N, # Use the function argument for sample size

      site_draw = fabricatr::draw_categorical(
        prob = site_props$Proportion,
        category_labels = site_props$Site
      ),

      race_draw = fabricatr::draw_conditional(
        condition_on = site_draw,
        FUN = function(data_so_far) {
          conditioned_draws <- character(nrow(data_so_far))
          current_sites <- as.character(data_so_far$site_draw)
          unique_sites_in_data <- unique(current_sites)
          for (s_unique in unique_sites_in_data) {
            indices <- which(current_sites == s_unique)
            site_specific_race_props <- race_props_by_site[as.character(race_props_by_site$Site) == s_unique, ]
            if(nrow(site_specific_race_props) > 0 && sum(site_specific_race_props$Proportion, na.rm = TRUE) > 0) {
              probs <- site_specific_race_props$Proportion / sum(site_specific_race_props$Proportion, na.rm = TRUE)
              if (nrow(site_specific_race_props) == 1) {
                 conditioned_draws[indices] <- rep(site_specific_race_props$Variable, length(indices))
              } else {
                 conditioned_draws[indices] <- sample(
                    x = site_specific_race_props$Variable,
                    size = length(indices),
                    replace = TRUE,
                    prob = probs
                 )
              }
            } else {
              conditioned_draws[indices] <- NA_character_
            }
          }
          return(conditioned_draws)
        }
      ),

      sex_indication_joint_temp_draw = fabricatr::draw_conditional(
        condition_on = site_draw,
        FUN = function(data_so_far) {
          conditioned_draws <- character(nrow(data_so_far))
          current_sites <- as.character(data_so_far$site_draw)
          unique_sites_in_data <- unique(current_sites)
          for (s_unique in unique_sites_in_data) {
            indices <- which(current_sites == s_unique)
            site_specific_sex_indication_props <- sex_indication_props_by_site[as.character(sex_indication_props_by_site$Site) == s_unique, ]
            if(nrow(site_specific_sex_indication_props) > 0 && sum(site_specific_sex_indication_props$Proportion, na.rm = TRUE) > 0) {
              combined_labels <- paste(site_specific_sex_indication_props$Sex, site_specific_sex_indication_props$Indication, sep = "_")
              probs <- site_specific_sex_indication_props$Proportion / sum(site_specific_sex_indication_props$Proportion, na.rm = TRUE)
              if (length(unique(combined_labels)) == 1) {
                conditioned_draws[indices] <- rep(unique(combined_labels)[1], length(indices))
              } else {
                conditioned_draws[indices] <- sample(
                  x = combined_labels,
                  size = length(indices),
                  replace = TRUE,
                  prob = probs
                )
              }
            } else {
              conditioned_draws[indices] <- NA_character_
            }
          }
          return(conditioned_draws)
        }
      ),

      post_processing = function(data) {
        data |>
          tidyr::separate(sex_indication_joint_temp_draw, into = c("Sex", "Indication"), sep = "_", remove = TRUE, fill = "right") |>
          dplyr::rename(Site = site_draw, Race = race_draw) |>
          dplyr::mutate(
            Site = factor(Site),
            Race = factor(Race),
            Sex = factor(Sex),
            Indication = factor(Indication)
          )
      }
    ),

    # Assignment Step infers N from the population data
    assignment = DeclareDesign::declare_assignment(
      T1 = simple_ra(prob = 0.5),
      T2 = simple_ra(prob = 0.5)
    ),

    # Potential Outcomes Step also infers N (number of units in the step's data)
    potential_outcomes = DeclareDesign::declare_potential_outcomes(
      Y_latent ~ 5 + 1.0 * T1 + 1.5 * T2 + 0.5 * T1*T2 + rnorm(N, mean = 0, sd = 3),
      Y_ord = as.integer(as.character(base::cut(
        Y_latent,
        breaks = c(-Inf, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, Inf),
        labels = 0:10,
        include.lowest = TRUE,
        right = FALSE
      )))
    ),

    estimands = DeclareDesign::declare_estimand(
      `T1` = mean( (Y_ord_T1_1_T2_0 + Y_ord_T1_1_T2_1)/2 - (Y_ord_T1_0_T2_0 + Y_ord_T1_0_T2_1)/2 ),
      `T2` = mean( (Y_ord_T1_0_T2_1 + Y_ord_T1_1_T2_1)/2 - (Y_ord_T1_0_T2_0 + Y_ord_T1_1_T2_0)/2 ),
      `T1:T2` = mean( (Y_ord_T1_1_T2_1 - Y_ord_T1_1_T2_0) - (Y_ord_T1_0_T2_1 - Y_ord_T1_0_T2_0) )
    ),

    estimators = DeclareDesign::declare_estimator(
      Y_ord ~ T1 * T2,
      model = estimatr::lm_robust,
      term = c("T1", "T2", "T1:T2"),
      label = "OLS"
    )
  ) # End of declare_design

  return(design)
} # End of CreateSampleDesign function

# To make the script runnable and testable when sourced or run directly:
if (sys.nframe() == 0) {

  N_for_test <- 1000 # Define a sample size for this test run

  cat(paste("\nCreating a test design instance with N =", N_for_test, "\n"))
  test_design_instance <- CreateSampleDesign(current_N = N_for_test)

  cat("Summary of the test design instance:\n")
  print(test_design_instance)

  cat("\nDiagnosing the test design instance (sims = 100, this may take a moment)...\n")
  set.seed(456) # For reproducibility of diagnose_design simulations
  diagnostics_for_test_design <- DeclareDesign::diagnose_design(test_design_instance, sims = 100)

  cat("Design Diagnostics for N =", N_for_test, ":\n")
  print(diagnostics_for_test_design)

  # The detailed simulation and printing of one data draw can be added here if desired
  # For example:
  # cat("\nSimulating one draw from the test design (N =", N_for_test, ")...\n")
  # set.seed(123) # For reproducibility of this specific simulation draw
  # simulated_data_one_draw <- DeclareDesign::simulate_design(test_design_instance)
  # cat("Head of simulated data for N =", N_for_test, ":\n")
  # print(head(simulated_data_one_draw))
  # cat("\nOutcome (Y_ord) summary for N =", N_for_test, ":\n")
  # print(table(simulated_data_one_draw$Y_ord, useNA = "ifany"))
}
