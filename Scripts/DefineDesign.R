# Load necessary libraries
library(DeclareDesign)
library(fabricatr)
library(dplyr)
library(tidyr)
library(estimatr)
library(stats) # For qlogis, plogis (though often available without explicit load)

# Source data loading script for proportions
source("Scripts/LoadAndProcessData.R") # Also loads outcome_dist_scenario1 and outcome_dist_scenario2

# Define a function to create the design, parameterized by sample size N and outcome scenario
CreateSampleDesign <- function(current_N, outcome_scenario) {

  # Define Beta Coefficients for treatment effects
  beta_T1 <- 1.0
  beta_T2 <- 1.5
  beta_interaction <- 0.5

  population_step <- DeclareDesign::declare_population(
      population = fabricatr::fabricate(
          persons = fabricatr::level(
              N = current_N,
              site_draw_pers = fabricatr::draw_categorical(
                  prob = site_props$Proportion,
                  category_labels = site_props$Site
              ),
              race_draw_pers = fabricatr::draw_conditional(
                  condition_on = site_draw_pers,
                  FUN = function(data) {
                      conditioned_draws <- character(nrow(data))
                      current_sites <- as.character(data$site_draw_pers)
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
                        } else { conditioned_draws[indices] <- NA_character_ }
                      }
                      return(conditioned_draws)
                  }
              ),
              sex_indication_joint_temp_pers = fabricatr::draw_conditional(
                  condition_on = site_draw_pers,
                  FUN = function(data) {
                      conditioned_draws <- character(nrow(data))
                      current_sites <- as.character(data$site_draw_pers)
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
                        } else { conditioned_draws[indices] <- NA_character_ }
                      }
                      return(conditioned_draws)
                  }
              ),
              Indication_derived_pers = sapply(strsplit(sex_indication_joint_temp_pers, "_"), `[`, 2),
              num_visits_total_pers = fabricatr::draw_conditional(
                  condition_on = Indication_derived_pers,
                  FUN = function(data) {
                      ifelse(data$Indication_derived_pers == "Bladder Cancer",
                             sample(1:4, size = nrow(data), replace = TRUE, prob = c(0.25, 0.35, 0.25, 0.15)),
                             1)
                  }
              )
          ),
          visits = fabricatr::add_level(
              N_by = "num_visits_total_pers",
              visit_number_by_person = 1:.N,
              is_first_visit_vis = (visit_number_by_person == 1),
              naivety_if_first_vis = fabricatr::draw_binary(prob = 0.4),
              Naivety_vis = fabricatr::modify_if(is_first_visit_vis, naivety_if_first_vis, FALSE)
          )
      ),
      data_transform = function(fabricated_data) {
          transformed_data <- fabricated_data |>
            tidyr::separate(persons.sex_indication_joint_temp_pers,
                            into = c("Sex_pers_char", "Indication_dummy"),
                            sep = "_", remove = FALSE, fill = "right") |>
            dplyr::rename(
                person_id = persons.persons_id,
                Site = persons.site_draw_pers,
                Race = persons.race_draw_pers,
                Sex = Sex_pers_char,
                Indication = persons.Indication_derived_pers,
                num_visits_total = persons.num_visits_total_pers,
                visit_id = visits.visits_id,
                visit_number = visits.visit_number_by_person,
                Naivety = visits.Naivety_vis
            ) |>
            dplyr::select(
                person_id, Site, Race, Sex, Indication, num_visits_total,
                visit_id, visit_number, Naivety
            ) |>
            dplyr::mutate(
                Site = factor(Site),
                Race = factor(Race),
                Sex = factor(Sex),
                Indication = factor(Indication),
                Naivety = factor(Naivety, levels = c(TRUE, FALSE))
            )
          return(transformed_data)
      }
  )

  assignment_step <- DeclareDesign::declare_assignment(
    handler = function(data) {
      all_treatment_combinations_df <- data.frame(T1 = c(0,0,1,1), T2 = c(0,1,0,1))
      processed_data_list <- split(data, data$person_id)
      assigned_data_list <- lapply(processed_data_list, function(person_df) {
        n_visits <- person_df$num_visits_total[1]
        if (n_visits == 1) {
          person_df$T1 <- sample(0:1, size = 1, replace = TRUE)
          person_df$T2 <- sample(0:1, size = 1, replace = TRUE)
        } else {
          sampled_indices <- sample(1:nrow(all_treatment_combinations_df), size = n_visits, replace = FALSE)
          distinct_treatments_for_person <- all_treatment_combinations_df[sampled_indices, ]
          person_df <- person_df[order(person_df$visit_number), ]
          person_df$T1 <- distinct_treatments_for_person$T1
          person_df$T2 <- distinct_treatments_for_person$T2
        }
        return(person_df)
      })
      final_assigned_data <- dplyr::bind_rows(assigned_data_list)
      return(final_assigned_data)
    }
  )

  potential_outcomes_step <- DeclareDesign::declare_potential_outcomes(
    Y_ord = {
      .DrawYordInternal <- function(t1_value, t2_value, N_units, scenario_idx_internal,
                                    b_T1, b_T2, b_int, dist1, dist2) {
        baseline_probs <- if (scenario_idx_internal == 1) dist1 else dist2
        baseline_cum_probs <- cumsum(baseline_probs)
        if (abs(baseline_cum_probs[length(baseline_cum_probs)] - 1.0) > 1e-9) {
        }
        baseline_cum_probs[length(baseline_cum_probs)] <- 1.0
        alpha_k_baseline <- stats::qlogis(baseline_cum_probs[1:(length(baseline_cum_probs)-1)])
        effect_score <- b_T1 * t1_value + b_T2 * t2_value + b_int * t1_value * t2_value
        alpha_k_treated <- alpha_k_baseline - effect_score
        treated_cum_probs_intermediate <- stats::plogis(alpha_k_treated)
        full_treated_cum_probs <- c(treated_cum_probs_intermediate, 1.0)
        treated_individual_probs <- diff(c(0, full_treated_cum_probs))
        treated_individual_probs[treated_individual_probs < 0] <- 0
        treated_individual_probs <- treated_individual_probs / sum(treated_individual_probs)
        outcomes <- sample(0:10, size = N_units, replace = TRUE, prob = treated_individual_probs)
        return(outcomes)
      }
      .DrawYordInternal(T1, T2, N, outcome_scenario,
                        beta_T1, beta_T2, beta_interaction,
                        outcome_dist_scenario1, outcome_dist_scenario2)
    }
  )

  estimators_step <- DeclareDesign::declare_estimator(
    formula = Y_ord ~ T1 * T2,
    model = estimatr::lm_robust,
    clusters = person_id,
    term = c("T1", "T2", "T1:T2"),
    label = "OLS_clustered"
  )

  design <- DeclareDesign::declare_design(
    population = population_step,
    assignment = assignment_step,
    potential_outcomes = potential_outcomes_step,
    estimands = DeclareDesign::declare_estimand(
      `T1` = mean( (Y_ord_T1_1_T2_0 + Y_ord_T1_1_T2_1)/2 - (Y_ord_T1_0_T2_0 + Y_ord_T1_0_T2_1)/2 ),
      `T2` = mean( (Y_ord_T1_0_T2_1 + Y_ord_T1_1_T2_1)/2 - (Y_ord_T1_0_T2_0 + Y_ord_T1_1_T2_0)/2 ),
      `T1:T2` = mean( (Y_ord_T1_1_T2_1 - Y_ord_T1_1_T2_0) - (Y_ord_T1_0_T2_1 - Y_ord_T1_0_T2_0) )
    ),
    estimators = estimators_step
  )

  return(design)
}

# To make the script runnable and testable when sourced or run directly:
if (sys.nframe() == 0) {

  N_persons_for_test <- 250
  outcome_scenario_for_test <- 1

  cat(paste("\nCreating a test design instance with N_persons =", N_persons_for_test,
              ", outcome_scenario =", outcome_scenario_for_test, "\n"))
  test_design_instance <- CreateSampleDesign(current_N = N_persons_for_test,
                                           outcome_scenario = outcome_scenario_for_test)

  # cat("Summary of the test design instance:\n") # Optional: print design summary
  # print(test_design_instance)

  cat("\n--- Simulating one data draw for verification (N_persons =", N_persons_for_test, ")...\n")
  set.seed(12345) # Consistent seed for verification data
  simulated_data_one_draw <- DeclareDesign::simulate_design(test_design_instance) # sims=1 is default

  cat("--- Data Sample (first 10 rows) ---\n")
  print(head(simulated_data_one_draw, 10))

  if (!is.null(simulated_data_one_draw)) {
    cat("\n--- Verification: Number of Visits ---\n")
    person_summary_visits <- simulated_data_one_draw |>
      dplyr::distinct(person_id, Indication, num_visits_total)

    cat("Table of num_visits_total by Indication (for persons):\n")
    print(table(person_summary_visits$Indication, person_summary_visits$num_visits_total, useNA = "ifany", dnn = c("Indication", "NumVisits")))

    bc_visits <- person_summary_visits |>
                   dplyr::filter(Indication == "Bladder Cancer")
    if(nrow(bc_visits) > 0) {
      cat("Summary of num_visits_total for Bladder Cancer indication:\n")
      print(summary(bc_visits$num_visits_total))
      print(prop.table(table(bc_visits$num_visits_total, dnn = "NumVisits_BC")))
    } else {
      cat("No Bladder Cancer patients in this sample for detailed visit check.\n")
    }

    other_visits <- person_summary_visits |>
                      dplyr::filter(Indication != "Bladder Cancer" & !is.na(Indication))
    if(nrow(other_visits) > 0) {
      cat("Num_visits_total for Other Indications (should all be 1):\n")
      all_one <- all(other_visits$num_visits_total == 1)
      print(paste("All non-BC Indications have 1 visit:", all_one))
      if (!all_one) print(table(other_visits$Indication, other_visits$num_visits_total, dnn = c("Indication", "NumVisits_Other")))
    } else {
      cat("No 'Other' indication patients in this sample for visit check.\n")
    }

    cat("\n--- Verification: Treatment Assignment Uniqueness (Multi-Visit) ---\n")
    multi_visit_person_ids <- simulated_data_one_draw |>
      dplyr::filter(num_visits_total > 1) |>
      dplyr::distinct(person_id) |>
      dplyr::pull(person_id)

    # Check first few (e.g., up to 3) multi-visit persons
    persons_to_check_treatment <- head(multi_visit_person_ids, 3)

    if (length(persons_to_check_treatment) > 0) {
      cat("Treatment assignments for first few multi-visit persons:\n")
      for (pid_val in persons_to_check_treatment) {
        cat("Person ID:", pid_val, "\n")
        treatments_for_pid <- simulated_data_one_draw |>
          dplyr::filter(person_id == pid_val) |>
          dplyr::select(visit_number, T1, T2, num_visits_total) |>
          dplyr::arrange(visit_number)
        print(treatments_for_pid)
        distinct_rows <- treatments_for_pid |> dplyr::distinct(T1, T2)
        cat("  Number of visits:", nrow(treatments_for_pid),
            "| Number of unique (T1,T2) pairs:", nrow(distinct_rows), "\n")
        if(nrow(treatments_for_pid) != nrow(distinct_rows)) {
          cat("  WARNING: Treatments may not be unique for person_id", pid_val, "\n")
        }
      }
    } else {
      cat("No multi-visit persons in this sample to check treatment uniqueness.\n")
    }

    cat("\n--- Verification: Naivety Distribution ---\n")
    first_visits_data <- simulated_data_one_draw |> dplyr::filter(visit_number == 1)
    if (nrow(first_visits_data) > 0) {
      naivety_prop_first_visit <- mean(first_visits_data$Naivety == TRUE, na.rm = TRUE)
      cat(sprintf("Proportion of Naivety=TRUE for first visits: %.3f (target approx 0.40)\n",
                  naivety_prop_first_visit))
    } else {
      cat("No first visits in data? (Should not happen if N_persons_for_test > 0)\n")
    }

    later_visits_data <- simulated_data_one_draw |> dplyr::filter(visit_number > 1)
    if (nrow(later_visits_data) > 0) {
      all_naivety_false_later_visits <- all(later_visits_data$Naivety == FALSE, na.rm = TRUE)
      cat(sprintf("All Naivety=FALSE for visits > 1: %s\n", all_naivety_false_later_visits))
      if(!all_naivety_false_later_visits){
          cat("  WARNING: Some later visits have Naivety = TRUE.\n")
          print(table(later_visits_data$visit_number, later_visits_data$Naivety, dnn=c("VisitNum", "Naivety")))
      }
    } else {
      cat("No visits > 1 in this sample (e.g., no Bladder Cancer patients or all had 1 visit).\n")
    }

    cat("\n--- Verification: Y_ord Distribution for Control Group ---\n")
    # cat("Marginal Y_ord distribution (all data in this draw):\n") # This was printed before, can be less verbose here
    # print(prop.table(table(simulated_data_one_draw$Y_ord, useNA = "ifany")))

    control_group_data <- simulated_data_one_draw |>
                            dplyr::filter(T1 == 0 & T2 == 0) # Control group defined by T1=0, T2=0

    if (nrow(control_group_data) > 0) {
      cat(sprintf("Marginal Y_ord for Control Group (T1=0, T2=0) (Scenario %d):\n",
                  outcome_scenario_for_test))

      observed_dist_control <- prop.table(table(factor(control_group_data$Y_ord, levels=0:10))) # Ensure all levels 0-10 shown
      print(round(observed_dist_control,3))

      target_dist_vector <- if (outcome_scenario_for_test == 1) {
        outcome_dist_scenario1
      } else {
        outcome_dist_scenario2
      }
      target_dist_named <- setNames(target_dist_vector, 0:10)

      cat(sprintf("Compare with input Scenario %d distribution:\n", outcome_scenario_for_test))
      print(round(target_dist_named,3))

      comparison_df <- data.frame(
          Outcome = 0:10,
          Observed_Control = as.numeric(observed_dist_control),
          Target_Scenario = as.numeric(target_dist_vector)
      )
      comparison_df$Difference <- comparison_df$Observed_Control - comparison_df$Target_Scenario
      cat("Differences (Observed Control - Target Scenario):\n")
      print(round(comparison_df,3))
    } else {
      cat("No data for Control Group (T1=0, T2=0) in this sample to verify Y_ord distribution.\n")
    }
  } # End of if(!is.null(simulated_data_one_draw))

  # The diagnose_design call can be kept for a quick check of diagnostics,
  # but it's very slow with complex data generation for many sims.
  # Using low sims here just to ensure it runs without error.
  # cat("\n--- Quick Diagnostics Check (sims = 5) ---\n")
  # set.seed(456)
  # diagnostics_for_test_design <- DeclareDesign::diagnose_design(test_design_instance, sims = 5)
  # cat("Design Diagnostics for N_persons =", N_persons_for_test, ":\n")
  # print(diagnostics_for_test_design)
  cat("\n--- End of Design Script Test Block ---\n")
}
