#################################################
#### Data Generation
#################################################

# Note: Assumes that a settings file has already been sourced so that 
# variables like `N` are defined

require(tidyverse)

# Create tibble for study data
# We will generate data for four visits for every subject but for many subjects the visits after the first 
# are considered missing
study.data <- tibble(ID = rep(1:N, times = 4), Visit = c(rep(1, times = N), rep(2, times = N),
                                                         rep(3, times = N), rep(4, times = N)))

#### Generate random effects
random.intercepts <- rnorm(N, mean = 0, sd = sigma.intercept)

# Add random intercepts to the dataset
study.data$Ui <- rep(random.intercepts, times = 4)


# Generate errors
study.data$Errij <- rnorm(n = nrow(study.data), mean = 0, sd = sigma.noise)

# N obs per subject
study.data$Nobs <- rep(nobs.to.sample, times = 4)



#################################################
#### Generate Covariates
#################################################

#################################################
#### Assign Treatments
#################################################
# Generate Treatment Assignments


BlockRandomize <- function(n.to.randomize, treatments = 1:11, bs = block.size){
  trt.assignments <- vector(length = n.to.randomize)
  
  GenBlock <- function(trts, block){
    # Generate a single block randomization
    temp.trts <- rep(trts, times = block)
    assignments <- sample(temp.trts, replace = FALSE)
  }
  n.blocks <- ceiling(n.to.randomize/bs)
  # Subsetting is to deal with when the number of assignments to make is not perfectly divisble by the block size
  trt.assignments <- c(replicate(n = n.blocks, GenBlock(treatments, bs)))[1:n.to.randomize]
  
  return(trt.assignments)
}


if(randomization.method == "simple"){
  study.data$Treatment <- sample(seq(1,11, by = 1), size = nrow(study.data), replace=TRUE)
}

if(randomization.method == "block"){

  study.data$Treatment <- BlockRandomize(nrow(study.data))
}

if(randomization.method == "trajectory"){
  #Split the data up
  non.trajectory <- study.data %>% filter(Nobs <4)
  to.randomize <- study.data %>% filter(Nobs == 4) %>% arrange(ID)
  
  # Randomize subjects with 4 observations to trajectories
  # A trajectory is Soc, Treatment 1, Treatment 2, Treatment 1 & 2 in a random order
  n.trajectories <- nrow(to.randomize)/4
  interaction.assignment <- BlockRandomize(n.trajectories, treatments = 6:11)
  
  GenerateTrajectory <- function(interaction.trt){
    if(interaction.trt == 6){
      #Dwell x Music
      temp.trajectory <- c(1, 2, 3, 6)
    }
    if(interaction.trt == 7){
      #Dwell x Viz
      temp.trajectory <- c(1, 2, 4, 7)
    }
    if(interaction.trt == 8){
      #Dwell x Squeeze
      temp.trajectory <- c(1, 2, 5, 8)
    }
    if(interaction.trt == 9){
      #Music x Viz
      temp.trajectory <- c(1, 3, 4, 9)
    }
    if(interaction.trt == 10){
      #Music x Squeeze
      temp.trajectory <- c(1, 3, 5, 10)
    }
    if(interaction.trt == 11){
      #Viz x Squeeze
      temp.trajectory <- c(1, 4, 5, 11)
    }
    
    trajectory <- sample(temp.trajectory)
    return(trajectory)
  }
  
  assignments <- unlist(lapply(interaction.assignment, GenerateTrajectory))
  to.randomize$Treatment <- assignments
  
  non.trajectory$Treatment <- sample(seq(1,11, by = 1), size = nrow(non.trajectory),prob = trt.probs.for.non.traj, replace=TRUE)
  study.data <- bind_rows(non.trajectory, to.randomize)
  remove(non.trajectory, to.randomize)
}

#TODO: Add ability to tweak  randomization probabilities for other groups when the
#     4 visit subjects are randomized to trajectories

# Generate treatment indicators - 1 if a treatment was received (either by itself or in combination) 0 else
study.data$Dwell <- ifelse(study.data$Treatment %in% c(2,6, 7, 8), 1, 0)
study.data$Music <- ifelse(study.data$Treatment %in% c(3,6, 9, 10), 1, 0)
study.data$Viz <- ifelse(study.data$Treatment %in% c(4,7, 9, 11), 1, 0)
study.data$Squeeze <- ifelse(study.data$Treatment %in% c(5,8, 10, 11), 1, 0)
#################################################
#### Generate Potential Outcomes
#################################################

# Calculate the treatment effects based on the treatment assignment
study.data$TEij <- NA
for(i in 1:nrow(study.data)){
  study.data$TEij[i] <- trt.effects[study.data$Treatment[i]]
}

# Generate the outcomes
study.data$Yij <- study.data$TEij + study.data$Ui + study.data$Errij

#################################################
#### Generate Observed Outcomes
#################################################
CalcObsij <- function(in.obs){
  in.obs <- as.data.frame(t(in.obs))
  if(in.obs$Nobs == 4){
    obsij <- in.obs$Yij
  }
  
  if(in.obs$Nobs == 3){
    if(in.obs$Visit <= 3){
      obsij <- in.obs$Yij
    } else{
      obsij <- NA
    }
  }
  if(in.obs$Nobs == 2){
    if(in.obs$Visit <= 2){
      obsij <- in.obs$Yij
    } else{
      obsij <- NA
    }
  }
  
  if(in.obs$Nobs == 1){
    if(in.obs$Visit == 1){
      obsij <- in.obs$Yij
    } else{
      obsij <- NA
    }
  }
  
  return(obsij)
}

# Generate normal outcome
study.data$ObsNormij <- apply(study.data, 1, CalcObsij)

# Generate ordinal outcome
# breaks are defined in the settings file
# Subtract 1 because the R default begins numbering at 1
study.data$Obsij <- cut(study.data$ObsNormij, breaks = ordinal.breaks, labels=FALSE) -1

# Subset to include only observed data
# The arranging is necessary for geepack to work correctly downstream
observed.data <- study.data %>% filter(!is.na(Obsij)) %>% arrange(ID)
#################################################
#### Export data
#################################################

#TODO: Add data exporting
