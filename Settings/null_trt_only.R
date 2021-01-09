#################################################
#### Simulation Settings
#################################################

require(tidyverse)

####### Simulation Settings


# https://nclottery.com/Pick4 Thursday 12-31-20 Daytime Draw
sim.seed <- 2373

# Number of simulation runs
sim.reps <- 2500

# Used as part of the name for creating the results directory
# Directory is Run-type Date Time
settings.type <- "nsg_null"

# Is the scenario a null scenario?
# If true, then the FWER is computed
null.scenario <- TRUE

# Save final simulation result summaries - includes pvals, treatment allocation summary, and settings files
save.results <- TRUE
# Save intermediate simulation run data files
# If true, saves the data for each simulation into a CSV
save.sim.study.data <- FALSE

#################################################
#### Study Design Settings
#################################################

# Total sample size
N <- 2900

# Visit proportions
# 1 Visit, 2 visits, 3 visits, 4 visits
visit.proportions <-  c(.7, .12, .1, .08)

category.numbers <- N * visit.proportions
nobs.to.sample <- c(rep(1, times = category.numbers[1]),
                    rep(2, times = category.numbers[2]),
                    rep(3, times = category.numbers[3]),
                    rep(4, times = category.numbers[4]))

# First batch sample size
# legacy setting from the adaptive design
first.batch <- N

# Subsequent batch sizes
batch.increment <- 0

# Number of uniform (-1,1) covariates
k.unif <- 1

# Number of standard normal covariates
k.norm <- 1

# Number of binary covariates
# Not currently used except to check against the length of bin.props
k.bin <- 2

# Must have length equal to the number of binary covariates
# First is gender - I think about 75% men 25% women
# Second is cancer
bin.props <- c(.35, .4)

# What proportion of people are having their first cystoscopy
first.cysto.prop <- .6

arm.var.name <- "Arm"

#################################################
#### Treatment Effect Model Settings
#################################################
# Treatment effects
param.df <- read.csv("./Settings/trt_only_null.csv")
true.params <- param.df$Coefficient
names(true.params) <- param.df$Parameter

treatment.arm.map <- expand_grid(Dwell = c(0, 1),
                                 Music = c(0, 1)) %>% 
  mutate(Arm = 1:n())

#One sided formula because it's used to generate the outcome data
true.model.formula <- formula(~ 1 + (Dwell + Music )^2)


mixed.model.string <- "Yij ~ (1|ID) + (Dwell + Music)^2"
gee.model.formula <- formula(Yij ~ 1 + (Dwell + Music)^2)

# Noise parameters
# Because of the discretization of the response, to achieve an ICC of .70 on the observed
# response scale the ICC of the latent normal RV is slightly different
# Standard deviation for the random intercepts
sigma.intercept <- sqrt(.74)

# Standard deviation of the within-subject (residual) variance component
sigma.noise <- sqrt(.26)

# Ordinal cutpoints
ref_prob <- c(.115, 0.18, 0.22, 0.23, 0.150, 0.065, 0.0175, 0.0075, 0.005, 0.005, 0.005)
ref_cum <- cumsum(ref_prob)
ordinal.breaks <- c(-Inf, qnorm(ref_cum))

################################################################################
#### Hypothesis Testing
##
################################################################################

contrast.mat <- diag(4)[-1,]
rownames(contrast.mat) <- c("Dwell", "Music", "Dwell:Music")
