#################################################
#### Simulation Settings
#################################################

####### Simulation Settings
# Save intermediate simulation run data files
# If true, saves the data for each simulation into a CSV
save.data.files <- FALSE

# Number of simulation runs
sim.reps <- 1500

# Number of uniform (-1,1) covariates
p.unif <- 1

# Number of standard normal covariates
p.norm <- 1

# Number of binary covariates
# Not currently used except to check against the length of bin.props
p.bin <- 2

# Must have length equal to the number of binary covariates
# First is gender - I think about 75% men 25% women
bin.props <- c(.25, .5)
#################################################
#### Study Design Settings
#################################################

# Total sample size
N <- 3000

# Visit proportions
# 1 Visit, 2 visits, 3 visits, 4 visits
visit.proportions <- c(.6, .05, .1, .25)

category.numbers <- N * visit.proportions
nobs.to.sample <- c(rep(1, times = category.numbers[1]),
                    rep(2, times = category.numbers[2]),
                    rep(3, times = category.numbers[3]),
                    rep(4, times = category.numbers[4]))

# First batch sample size
first.batch <- 1500

# Subsequent batch sizes
batch.increment <- 500

#################################################
#### Randomization Settings
#################################################
# Random Assignment method
# Options: "simple", "block", "trajectory"
 randomization.method <- "trajectory"
 # Set the block size when use block randomization
 block.size <-4
 # Flag that determines whether subjects who will have four visits will be randomized separately
 # It will be known at enrollment whether a subject will be expected to return four times during the study
 # If TRUE, subjects who return four times will be randomized to a trajectory where a trajectory is defined as
 # Soc, Treatment A, Treatment B, Treatment A&B in a random order
 # To facilitate within subject comparisons when possible
 
 trt.probs.for.non.traj <- c(.08, .05, .05, .05, .05, .12, .12, .12, .12, .12, .12)
 # randomize.to.trajectory <- FALSE

 #################################################
 #### Treatment Effect Model Settings
 #################################################
# Treatment effects
param.df <- read.csv("./Settings/parameter_ests.csv")
trt.design <- read.csv("./Settings/treatment_design_matrix.csv", header=TRUE)
trt.names <- trt.design$X
trt.design <- data.matrix(trt.design[,-1])
rownames(trt.design) <- trt.names

trt.effects <- trt.design %*% param.df$Coefficient

# Noise parameters
# Standard deviation for the random intercepts
sigma.intercept <- 1

# Standard deviation of the within-subject (residual) variance component
sigma.noise <- 1

# Ordinal cutpoints
ordinal.breaks <- c(-Inf, 0, 1, 2, 3, 4, 5, 6, 7, 7.5, 8, Inf)
