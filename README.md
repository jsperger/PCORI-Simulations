# PCORI-Simulations

This repository contains the code used to conduct the sample size simulations for the PCORI APDTO research grant application. 

## Data Generation

Based on a review of the published literature, we assumed treatment effect sizes of -0.2 for long relative to short lidocaine dwell time, and -0.4 for music relative to visualization. No estimates were available for the effect sizes of interaction terms, and so we powered to detect interaction effects of the same magnitude as the main effects  because these have clinical significance because the best treatment for the subgroup may change (depending on other covariates). The exact settings used can be seen in the Settings sub-folder. The data were generated using a random intercept model to induce the correlation between observations for the same participant, and a latent normal RV with cut-points based on quantiles of the standard normal distribution was used to generate the observed ordinal responses. The cut-points defining the bins are controlled in the settings file. The latent continuous outcome is only used for data generation; the model fitting for the power simulations only uses the observed ordinal response. Note that the coefficients in the CSV files are on the scale of this latent variable, not on the scale of the 0-10 VAS. 

For returning patients, an exchangeable correlation structure with an intra-class correlation of 0.7 was assumed. No studies measuring the ICC specific to cystoscopy exist to our knowledge, but the literature on ICC for pain scores reports that they are typically very high, above 0.84, with a higher ICC for less painful stimuli. In the context of our pain score distribution, an ICC of .70 corresponds to an average session-to-session change in the absolute VAS pain score of 1.03 and a 32% probability of having the same score.

### Analysis

 GEE with an exchangeable correlation structure was used to fit the model for all scenarios. For all scenarios, the hypotheses tests were tests of the null hypothesis that the relevant parameters were equal to zero, and a Bonferroni-Holm correction applied. 

## Results

The results directory contains a sub-directory for every simulation run. The simulation sub-directories 
contain the parameters for the data generating model as a CSV,  text files with a snapshot of RunSimulation.R 
at the time of the simulation and of the settings R script used for the simulation, and the resulting
unadjusted p-values in a CSV. 

The resulting power and family-wise Type I error rate (adjusting for multiple comparisons using the Holm procedure) for the three scenarios are summarized in the table below

| Scenario | Power -  Dwell Time > 10 ME | Power -  Music | Power -  Smallest Interaction | Type I Error |
|:--------:|:---------------------------:|:--------------:|:-----------------------------:|:------------:|
|     1    |             97%             |      >99%      |              84%              |      4%      |
|     2    |             90%             |      >99%      |              81%              |      5%      |
|     3    |             91%             |      >99%      |              80%              |      4%      |

## Running a Simulation Scenario

1. Change 'settings.path <- "./Settings/full_scenario.R"' in the RunSimulation.R file to the path of the scenario you wish to run
2. Change the number of mc.cores in RunSimulation.R to a system-appropriate number
3. Source the RunSimulation.R file and the simulation will run and the results will be saved. 

Additional options are available in the settings files, and the scenarios can be changed by modifying the CSV files in the Settings directory. 