# PCORI-Simulations

This repository contains the code used to conduct the sample size simulations for the PCORI APDTO research grant application. 

## Data Generation

Based on a review of the published literature we estimated treatment effect sizes of -0.4, -1.2, -0.8, and -.9 for dwell time, music, visualization, and bag squeeze, respectively. For the combination interventions no effect size estimates were available, and we assume that the effect of adding a second intervention is sub-additive with the combination intervention’s effect approximately equal to the effect of the stronger intervention plus one half the effect of the weaker intervention. The exact settings used can be seen in the Settings sub-folder. For returning patients, an exchangeable correlation structure with ρ=.5 was assumed. 

The data were generated using a random intercept model to induce the correlation between observations for the same participant. To generate the VAS pain scale outcome an underlying normal random error is added and then the 0-10 pain score is created by binning the normal response. The cut-points definining the bins are controlled in the settings file. The normal response is only used for data generation; the intermediate working model for the adaptive randomization and all the final models used to determine power use the ordinal response and not the normal response. 

### Simplifications and Approximations

An approximate method is used for Thompson Sampling (the adaptive randomization). The working model is fit using a mixed model and REML, and then a parametric bootstrap draw is used to approximate sampling from the posterior distribution. 

## Results

|Scenario     |Method   | Power| PercentBest|
|:------------|:--------|-----:|-----------:|
|No Subgroups |Fixed    | 1.000|       0.091|
|No Subgroups |Adaptive | 0.980|       0.206|
|1            |Fixed    | 0.861|       0.087|
|1            |Adaptive | 0.922|       0.196|
|2            |Fixed    | 0.736|       0.094|
|2            |Adaptive | 0.795|       0.179|