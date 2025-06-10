# Load necessary libraries
library(readr)
library(dplyr)

# Read the data
demo_by_site_raw <- read_csv("data/anon_demo.csv")
sex_indication_by_site_raw <- read_csv("data/anon_sex.csv")

# Calculate race proportions by site
race_props_by_site <- demo_by_site_raw |>
  dplyr::group_by(Site, Variable) |>
  dplyr::summarise(Count = sum(Count), .groups = 'drop') |>
  dplyr::group_by(Site) |>
  dplyr::mutate(Proportion = Count / sum(Count)) |>
  dplyr::ungroup() |>
  dplyr::select(Site, Variable, Proportion)

# Calculate joint sex/indication proportions by site
sex_indication_props_by_site <- sex_indication_by_site_raw |>
  dplyr::group_by(Site, Sex, Indication) |>
  dplyr::summarise(Count = sum(Count), .groups = 'drop') |>
  dplyr::group_by(Site) |>
  dplyr::mutate(Proportion = Count / sum(Count)) |>
  dplyr::ungroup() |>
  dplyr::select(Site, Sex, Indication, Proportion)

# Calculate site proportions based on total procedures
site_counts <- sex_indication_by_site_raw |>
  dplyr::group_by(Site) |>
  dplyr::summarise(total_procedures_site = sum(Count), .groups = 'drop')

site_props <- site_counts |>
  dplyr::mutate(Proportion = total_procedures_site / sum(total_procedures_site)) |>
  dplyr::select(Site, Proportion)

# Print the head of the resulting data frames
print("Race Proportions by Site:")
print(head(race_props_by_site))

print("Sex/Indication Proportions by Site:")
print(head(sex_indication_props_by_site))

print("Site Proportions:")
print(head(site_props))

# Define and normalize outcome distributions
# These represent probabilities for an ordinal outcome with 11 categories (e.g., 0 to 10)
outcome_dist_scenario1_raw <- c(0.30, 0.09, 0.33, 0.06, 0.11, 0.03, 0.04, 0.02, 0.03, 0.01, 0.0)
outcome_dist_scenario2_raw <- c(0.115, 0.180, 0.220, 0.230, 0.150, 0.065, 0.018, 0.007, 0.005, 0.005, 0.005)

# Normalize them to ensure they sum to 1
outcome_dist_scenario1 <- outcome_dist_scenario1_raw / sum(outcome_dist_scenario1_raw)
outcome_dist_scenario2 <- outcome_dist_scenario2_raw / sum(outcome_dist_scenario2_raw)

# Optional: Print to verify (can be commented out in final script for cleaner sourcing)
# print("Normalized Outcome Distribution - Scenario 1 (scores 0-10):")
# print(setNames(outcome_dist_scenario1, 0:10))
# print(paste("Sum Scenario 1:", sum(outcome_dist_scenario1)))
#
# print("Normalized Outcome Distribution - Scenario 2 (scores 0-10):")
# print(setNames(outcome_dist_scenario2, 0:10))
# print(paste("Sum Scenario 2:", sum(outcome_dist_scenario2)))
