source("AnalyzeStudyData.R")

#################################################
#### Data Generation
#################################################

# Note: Assumes that a settings file has already been sourced so that 
# variables like `N` are defined

require(tidyverse)

GenContinuousCovariates <- function(n.subj = N, d.unif = k.unif, d.norm = k.norm){
  # Generate Continuous Covariates
  # Continuous covariates are either standard normal or uniform (-1,1)
  # d.unif is the number of uniform RVs to generate for each subject
  # d.norm is the numebr of standard normal RVs to generate for each subject
  # Output is an N x (d.unif + d.norm) tibble
  U <- as_tibble(replicate(n = d.unif, runif(n=n.subj, min = -1, max = 1)))
  names(U) <- paste0("U", 1:d.unif)
  Nd <- as_tibble(replicate(n = d.norm, rnorm(n=n.subj, mean = 0, sd = 1)))
  names(Nd) <- paste0("N", 1:d.unif)
  X <- bind_cols(U, Nd)
  X$ID <- 1:n.subj
  return(X)
}

GenBinaryCovariates <- function(n.subj = N, d.bin = k.bin, props = .5){
  # Generate Binary Covariates
  # d.bin is the number of binary covariates to generate
  # props is a vector of length equal to d.bin or length 1 which 
  #     specifies the probability of a binary covariate being equal to 1
  # Output is an N x d.bin  tibble
  if(length(props) != d.bin & length(props) != 1){
    warning("The number of proportions specified does not match the stated number of binary covariates")
  }
  # props determines what proportion of the subjects have a value of 1
  X <- props %>% map_dfc(rbinom, n = n.subj, size = 1)
  names(X) <- paste0("B", 1:d.bin)
  X$ID <- 1:n.subj
  return(X)
}

GenTreatmentIndicators <- function(study.data){
  # Generate treatment indicators - 1 if a treatment was received (either by itself or in combination) 0 else
  
  treat.indicators <- tibble(Dwell = ifelse(study.data$Treatment %in% c(2,5, 6), 1, 0),
                             Music = ifelse(study.data$Treatment %in% c(3,5, 7), 1, 0),
                             Squeeze = ifelse(study.data$Treatment %in% c(4,6,7), 1, 0))
  return(treat.indicators)
  
}

#################################################
#### Randomization Methods
#################################################
BlockRandomize <- function(n.to.randomize, treatments = 1:7, bs = block.size){
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

TTTSRandomize <- function(study.data, first.batch, batch.increment, working.model.formula, pi.param){
  
  # Determine the waves
  # Waves are by ID
  # Simplifying assumption - we observe all observations for an individual in their wave
  ids <- unique(study.data$ID)
  n.waves <- 1 + ceiling((length(ids)-first.batch)/batch.increment)
  wave <- c(rep(1, times = first.batch), rep(seq(2, n.waves, by = 1), times = batch.increment))
  wave.assignment <- tibble(ID = ids, Wave = sample(wave, size = length(ids), replace = FALSE))
  study.data <- inner_join(study.data, wave.assignment, by = "ID")
  
  ####### Generate Data for the first wave and fit the working model based on the first wave
  wave.one <- study.data %>% filter(Wave == 1)
  #wave.one <- TrajectoryRandomize(study.data = wave.one, trt.probs.for.non.traj = trt.probs.for.non.traj)
  wave.one$Treatment <- BlockRandomize(n.to.randomize = nrow(wave.one),bs = 4)
  wave.one <- bind_cols(wave.one, GenTreatmentIndicators(wave.one))
  wave.one <- GenOutcomes(study.data = wave.one, true.model.formula = true.model.formula, true.params = param.df$Coefficient)
  working.model <- FitMixedModel(observed.data = wave.one, model.form = working.model.formula)
  
  observed.data <- wave.one
  
  TTTSAssignmentHelper <- function(cur.people, trajectory.flag, pi.param){
    # Helper function for assigning treatments to batches after the first
    # For those with 4 observations the trajectory flag should be set to TRUE
    # Trajectory flag indicates whether the patients should be randomize to a trajectory or to individual treatments
    if(trajectory.flag == TRUE){
      all.visits <- cur.people %>% group_by(ID) %>% arrange(ID)
      cur.people <- cur.people %>% filter(Visit == 1)
      # Only look at the interaction treatments
      data.dupe <- cur.people[rep(1:nrow(cur.people),each=6),]
      data.dupe$Treatment <- rep(6:11, times = nrow(cur.people))
    } else{
      data.dupe <- cur.people[rep(1:nrow(cur.people),each=7),]
      data.dupe$Treatment <- rep(1:7, times = nrow(cur.people))
    }
    
    data.dupe <- bind_cols(data.dupe, GenTreatmentIndicators(data.dupe))
    data.dupe$Pred <- PredictMM(working.model, newdat = data.dupe)
    
    if(trajectory.flag == TRUE){
      trt.visit.assign.helper <- data.dupe %>% group_by(ID, Visit) %>% summarise(BestVal = min(Pred), 
                                                                                 BestTrt = which.min(Pred)+5,
                                                                                 SecondTrt = which.min(Pred[-(BestTrt-5)])+5 +ifelse(which.min(Pred[-(BestTrt-5)]) < BestTrt-5, 0, 1),
                                                                                 CoinFlip = rbinom(n = 1, size = 1, p = pi.param),
                                                                                 Treatment = ifelse(CoinFlip == 1, BestTrt, SecondTrt))

      if(any(trt.visit.assign.helper$Treatment > 11) | any(trt.visit.assign.helper$Treatment < 6)){
        tp <- trt.visit.assign.helper %>% filter(Treatment > 11 | Treatment < 6)
        print(trt.visit.assign.helper)
      }
    } else{
      #The if-else is because by removing the best treatment we change the array indices for everything above the best treatment 
      trt.visit.assign.helper <- data.dupe %>% group_by(ID, Visit) %>% summarise(BestVal = min(Pred), 
                                                                                 BestTrt = which.min(Pred),
                                                                                 SecondTrt = which.min(Pred[-BestTrt]) + ifelse(which.min(Pred[-BestTrt]) < BestTrt, 0, 1),
                                                                                 CoinFlip = rbinom(n = 1, size = 1, p = pi.param),
                                                                                 Treatment = ifelse(CoinFlip == 1, BestTrt, SecondTrt))
    }
    
    cur.people <- inner_join(cur.people, trt.visit.assign.helper, by =c("ID", "Visit"))
    
    if(trajectory.flag == TRUE){
      assignments <- unlist(lapply(cur.people$Treatment, GenIndividualTrajectory))
      all.visits$Treatment <- assignments
      to.join <- cur.people %>% select(ID, BestTrt, SecondTrt, CoinFlip)
      
      to.return <- left_join(all.visits, to.join, by = "ID")
      return(to.return)
    } else{
      # Check if fewer than 5% are assigned to SoC
      soc.flag <- ifelse(sum(cur.people$Treatment == 1) < .05*nrow(cur.people), 1, 0)
      if(soc.flag == TRUE){
        assign.to.soc <- rbinom(n = nrow(cur.people), size = 1, p = .05)
        cur.people$Treatment <- ifelse(assign.to.soc == 1, 1, cur.people$Treatment)
      }
      return(cur.people)
    }
  }
  for(i in 2:n.waves){
    wave.to.randomize <- study.data %>% filter(Wave == i)
    #TODO Don't hard code trajectory flag
    trajectory.people <- TTTSAssignmentHelper(cur.people= wave.to.randomize %>% filter(Nobs == 4), 
                                              trajectory.flag = FALSE, pi.param=pi.param)
    non.trajectory.people <- TTTSAssignmentHelper(cur.people= wave.to.randomize %>% filter(Nobs != 4), 
                                              trajectory.flag = FALSE, pi.param=pi.param)
    
    current.wave <- bind_rows(trajectory.people, non.trajectory.people)
    current.wave <- bind_cols(current.wave, GenTreatmentIndicators(current.wave))
    current.wave <- GenOutcomes(current.wave, true.model.formula = true.model.formula, true.params = param.df$Coefficient)
    observed.data <- bind_rows(observed.data, current.wave)
  }
  return(observed.data)
}
#################################################
#### Generate Potential and Realized Outcomes
#################################################
GenOutcomes <- function(study.data, true.model.formula, true.params, oos.flag = FALSE){
  model.X <- model.matrix(true.model.formula, data=study.data)
  
  # Check that the design matrix and the parameter names match
  if(all(colnames(model.X) == names(true.params)) == FALSE){
    warning("Model formula variable names do not match parameter names")
  }
  
  # Calculate the outcome with no noise
  study.data$TEij <- c(model.X %*% true.params)
  
  if(oos.flag == TRUE){
    # No noise for the out of sample data
    random.ints <- tibble(ID = 1:max(study.data$ID), 
                          Ui = rnorm(n = length(unique(study.data$ID)), mean = 0, sd = 1.4))
    study.data <- inner_join(study.data, random.ints, by = "ID")
    study.data$Normij <- study.data$TEij + study.data$Ui
    study.data$Obsij <- cut(study.data$Normij, breaks = ordinal.breaks, labels=FALSE) -1
    return(study.data)
  }
  # Generate the observed outcome
  study.data$Yij <- study.data$TEij + study.data$Ui + study.data$Errij
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
  return(study.data)
  
}

#################################################
#### Generate Study Data
#################################################

GenerateStudyData <- function(N, nobs.to.sample, sigma.intercept, sigma.noise, randomization.method,
                              block.size, trt.probs.for.non.traj,ordinal.breaks,
                              k.unif, k.norm,
                              k.bin, bin.props,
                              true.model.formula, true.params, working.model.formula, ...){
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
  
  study.data <- inner_join(study.data, 
                           GenContinuousCovariates(n.subj = N, d.unif = k.unif, d.norm = k.norm), 
                           by = "ID")
  
  study.data <- inner_join(study.data, 
                           GenBinaryCovariates(n.subj = N, d.bin = k.bin, props = bin.props), 
                           by = "ID")
  
  
  #################################################
  #### Assign Treatments
  #################################################
  # Generate Treatment Assignments
  
  
  
  
  if(randomization.method == "simple"){
    study.data$Treatment <- sample(seq(1,7, by = 1), size = nrow(study.data), replace=TRUE)
  }
  
  if(randomization.method == "block"){
    
    study.data$Treatment <- BlockRandomize(nrow(study.data))
  }
  
  if(randomization.method == "ttts"){
    # Filter out the observed values ahead of time
    # TODO: Add the ability to not remove the unobserved values ahead of time for missing data analysis
    study.data <- TTTSRandomize(study.data = study.data %>% filter(Visit <= Nobs), 
                                first.batch = first.batch, 
                                batch.increment = batch.increment, 
                                working.model.formula = working.model.formula, pi.param = pi.param)
  } else{
    study.data <- bind_cols(study.data, GenTreatmentIndicators(study.data))
    study.data <- GenOutcomes(study.data = study.data, 
                              true.model.formula = true.model.formula, 
                              true.params = true.params)
  }
  #################################################
  #### Export data
  #################################################
  
  #TODO: Add the ability to write study data to a file 
  
  
  # Subset to include only observed data
  # TODO: Someday dump all the data for comparing missing data methods
  # The arranging is necessary for geepack to work correctly downstream
  observed.data <- study.data %>% filter(!is.na(Obsij)) %>% arrange(ID)
  
  return(observed.data)
}

