# PCORI-Simulations

This repository contains the code used to conduct the sample size simulations for the PCORI APDTO research grant application. 

## Data Generation

Based on a review of the published literature we estimated treatment effect sizes of -0.4, -1.2, -0.8, and -.9 for dwell time, music, visualization, and bag squeeze, respectively. For the combination interventions no effect size estimates were available, and we assume that the effect of adding a second intervention is sub-additive with the combination intervention’s effect approximately equal to the effect of the stronger intervention plus one half the effect of the weaker intervention. The exact settings used can be seen in the Settings sub-folder. For returning patients, an exchangeable correlation structure with ρ=.5 was assumed. 

The data were generated using a random intercept model to induce the correlation between observations for the same participant. To generate the VAS pain scale outcome an underlying normal random error is added and then the 0-10 pain score is created by binning the normal response. The cut-points definining the bins are controlled in the settings file. The normal response is only used for data generation; the intermediate working model for the adaptive randomization and all the final models used to determine power use the ordinal response and not the normal response. 

### Adaptive Randomization

As a simplication we assumed that we've observed all observations for the first 1,500 patients so for example if a patient has 4 visits throughout the course of the study we assumed that all 4 responses were available. 

An approximate method is used for Thompson Sampling (the adaptive randomization). The working model is fit using a mixed model and REML, and then a parametric bootstrap draw is used to approximate sampling from the posterior distribution. 


### Analysis

For all scenarios GEE with an exchangeable correlation structure was used to fit the model. For the no subgroup scenario a Wald test was used for each hypothesis test and a Bonferroni correction applied. The null hypotheses in the no subgroup scenario were of the form total intervention effect = 0 for each intervention. Thus for a combination intervention the null hypothesis was that the sum of the individual intervention parameters and the interaction parameter was equal to zero, not that the interaction parameter was equal to zero. For the scenarios with subgroups, all two-way interaction terms are included in the model even if the true parameter is zero (which would not be known in reality). No  regularization was done, and for the actual study we would use residual weighted learning which is a more efficient method for discovering treatment rules than fitting a model and then finding the best treatment using the model predicted outcomes. 

## Results

|Scenario     |Method   | Power| PercentBest|
|:------------|:--------|-----:|-----------:|
|No Subgroups |Fixed    | 1.000|       0.091|
|No Subgroups |Adaptive | 0.980|       0.206|
|1            |Fixed    | 0.861|       0.087|
|1            |Adaptive | 0.922|       0.196|
|2            |Fixed    | 0.736|       0.094|
|2            |Adaptive | 0.795|       0.179|