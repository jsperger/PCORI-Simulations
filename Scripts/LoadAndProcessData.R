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
