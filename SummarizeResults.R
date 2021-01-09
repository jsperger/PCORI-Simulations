library(tidyverse)
source("./Scripts/AnalyzeSimulationRuns.R")

sc1_pvals <- read_csv("Results/nsg_final-2021-01-05 13_48_21/gee_pvals.csv") %>% 
  AdjustPvals(., "holm")

sc2_pvals <- read_csv("Results/dwell_sex_sg-2021-01-05 14_19_09/gee_pvals.csv") %>% 
  AdjustPvals(., "holm")

sc3_pvals <- read_csv("Results/full-2021-01-05 15_52_40/gee_pvals.csv") %>% 
  AdjustPvals(., "holm")

null_sc1_pvals <- read_csv("Results/nsg_null-2021-01-05 14_29_11/gee_pvals.csv") %>% 
  AdjustPvals(., "holm")

null_sc2_pvals <- read_csv("Results/null_dwell_sex_sg-2021-01-05 14_24_02/gee_pvals.csv") %>% 
  AdjustPvals(., "holm")

null_sc3_pvals <- read_csv("Results/null_full-2021-01-05 14_08_15/gee_pvals.csv") %>% 
  AdjustPvals(., "holm")

results_summary <- tibble(Scenario = c("1", "2", "3"),
                          PowerDwellMainEffect = c(mean(sc1_pvals$Dwell < .05),
                                                   mean(sc2_pvals$Dwell < .05),
                                                   mean(sc3_pvals$Dwell < .05)),
                          PowerMusicMainEffect = c(mean(sc1_pvals$Music < .05),
                                                   mean(sc2_pvals$Music < .05),
                                                   mean(sc3_pvals$Music < .05)),
                          LowestPower = c(colMeans(sc1_pvals < .05) %>% min,
                                          colMeans(sc2_pvals < .05) %>% min,
                                          colMeans(sc3_pvals < .05) %>% min),
                          Type1Error = c(mean(rowSums(null_sc1_pvals < .05) >= 1),
                                         mean(rowSums(null_sc2_pvals < .05) >= 1),
                                         mean(rowSums(null_sc3_pvals < .05) >= 1)))


            