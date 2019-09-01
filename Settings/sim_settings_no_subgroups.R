#################################################
#### Simulation Settings
#################################################

require(tidyverse)

####### Simulation Settings
# Save intermediate simulation run data files
# If true, saves the data for each simulation into a CSV
save.intermediate.data.files <- FALSE
# Save final simulation result summaries - includes pvals, treatment allocation summary, and settings files
save.results <- TRUE

seed <- 1011

# Number of simulation runs
sim.reps <- 1000

# Used as part of the name for creating the results directory
# Directory is Run-type Date Time
settings.type <- "nsg-ts"

# Whether hypothesis tests should be conducted
# Should only be set to TRUE if there are no subgroups
test.hypotheses.flag <- TRUE

# Whether percentage of oracle value and other out of sample comparison metrics should be calculated
# Should be set to TRUE when there are subgroups
calc.oos.metrics.flag <- FALSE

# Whether the percentage of patients in-sample who received the best treatment 
# for them should be calculated
# 
calc.percentage.best.treat.flag <- TRUE
#################################################
#### Study Design Settings
#################################################

# Total sample size
N <- 3000

# Visit proportions
# 1 Visit, 2 visits, 3 visits, 4 visits
visit.proportions <- c(.6, .1, .1, .2)

category.numbers <- N * visit.proportions
nobs.to.sample <- c(rep(1, times = category.numbers[1]),
                    rep(2, times = category.numbers[2]),
                    rep(3, times = category.numbers[3]),
                    rep(4, times = category.numbers[4]))

# First batch sample size
first.batch <- 1500

# Subsequent batch sizes
batch.increment <- 500

# Number of uniform (-1,1) covariates
k.unif <- 1

# Number of standard normal covariates
k.norm <- 1

# Number of binary covariates
# Not currently used except to check against the length of bin.props
k.bin <- 2

# Must have length equal to the number of binary covariates
# First is gender - I think about 75% men 25% women
bin.props <- c(.25, .5)

#################################################
#### Randomization Settings
#################################################
# Random Assignment method
# Options: "simple", "block", "trajectory", "ttts"
 randomization.method <- "ttts"
 # Set the block size when use block randomization
 block.size <-4
 # Flag that determines whether subjects who will have four visits will be randomized separately
 # It will be known at enrollment whether a subject will be expected to return four times during the study
 # If TRUE, subjects who return four times will be randomized to a trajectory where a trajectory is defined as
 # Soc, Treatment A, Treatment B, Treatment A&B in a random order
 # To facilitate within subject comparisons when possible
 
 trt.probs.for.non.traj <- c(.08, .05, .05, .05, .05, .12, .12, .12, .12, .12, .12)
 # randomize.to.trajectory <- FALSE

 # Tuning parameter for Top-two Thompson Sampling
 # Probability that the best result from TS is used
 # Second best is used with probability 1-pi.param
 pi.param <- .6
 
 #################################################
 #### Treatment Effect Model Settings
 #################################################
# Treatment effects
#param.df <- read.csv("./Settings/null_parameters.csv")
param.df <- read.csv("./Settings/parameter_ests.csv")
true.params <- param.df$Coefficient
names(true.params) <- param.df$Parameter

#One sided formula because it's used to generate the outcome data
true.model.formula <- formula(~ 1 + (Dwell + Music + Viz + Squeeze)^2)
# Formula used for the working model for the TS fits
#working.model.formula <- formula(Obsij ~ 1 + (Dwell + Music + Viz + Squeeze)^2)
working.model.formula <- "Obsij ~ (1|ID) + (Dwell + Music + Viz + Squeeze)^2"
gee.model.formula <- formula(Obsij ~ 1 + (Dwell + Music + Viz + Squeeze)^2)

# Noise parameters
# Standard deviation for the random intercepts
sigma.intercept <- 1

# Standard deviation of the within-subject (residual) variance component
sigma.noise <- 1

# Ordinal cutpoints
ordinal.breaks <- c(-Inf, 0, 1, 2, 3, 4, 5, 6, 7, 7.5, 8, Inf)

#################################################
#### Hypothesis Test Setup
#################################################
if(test.hypotheses.flag == TRUE){
GenTreatmentIndicators <- function(study.data){
  # Generate treatment indicators - 1 if a treatment was received (either by itself or in combination) 0 else
  
  treat.indicators <- tibble(Dwell = ifelse(study.data$Treatment %in% c(2,6, 7, 8), 1, 0),
                             Music = ifelse(study.data$Treatment %in% c(3,6, 9, 10), 1, 0),
                             Viz = ifelse(study.data$Treatment %in% c(4,7, 9, 11), 1, 0),
                             Squeeze = ifelse(study.data$Treatment %in% c(5,8, 10, 11), 1, 0))
  return(treat.indicators)
  
}
# Create the contrast matrix for the hypothesis tests
# Note: not actually a contrast matrix
# Each row is tested separately, not a joint hypothesis
temp.study <- tibble(Treatment = 1:11)
temp.study <- bind_cols(temp.study, GenTreatmentIndicators(temp.study))
treat.mat <- model.frame(true.model.formula, data = temp.study)
contrast.mat <- model.matrix(true.model.formula, treat.mat)[-1,]
contrast.mat[,1] <- 0

# alternative contrasts - test every coefficient except intercept
#contrast.mat <- diag(11)
#contrast.mat <- contrast.mat[2:11,]

rm(temp.study, treat.mat)
}

#################################################
#### Treatment Assignment Table Setup
#################################################

trt.table.cols <- 11
table.call <- expr(table(current.study.data$Treatment))


#################################################
#### Out of Sample Comparisons
#################################################

Noos <- 5000