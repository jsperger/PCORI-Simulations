#################################################
#### Run a simulation
#################################################
settings.path <- "./Settings/sim_settings_sg1.R"
# load settings
source(settings.path)
source("./GenerateStudyData.R")
source("./AnalyzeStudyData.R")

set.seed(seed)
cur.time <- gsub(":", "_", x=Sys.time() ,fixed=TRUE)

if(test.hypotheses.flag == TRUE){
  # Create the containers for where p-values and treatment assignments 
  # from the simulation runs will be stored
  mixed.model.pvals.unadj <- matrix(nrow = sim.reps, ncol = 10)
  #colnames(mixed.model.pvals.unadj) <- names(true.params)
  # TODO: Decide what to do with the naming schemes esp. for the subgroup cases
  gee.pvals.unadj <- matrix(nrow = sim.reps, ncol = 10)
  #colnames(gee.pvals.unadj) <- trt.names[-1]
  
}

if(calc.percentage.best.treat.flag == TRUE){
  # Create the empty vectors for holding treatment assignments from each simulation rep
  treatment.assignments <- matrix(nrow = sim.reps, ncol = trt.table.cols)
  percentage.best.treatment <- vector(length = sim.reps)
}

if(calc.oos.metrics.flag == TRUE){
  # Create the empty matrix for the results summaries for each simulation rep
  percentage.of.oracle.value <- matrix(nrow = sim.reps, ncol = 5)
  colnames(percentage.of.oracle.value) <- c("PercOracleORD", "MSE-Norm", "MSE-SoC", "PerImpOverTO", "MissclassPerc")
}

start.time <- Sys.time()

#################################################
#### Run the Simulation
#################################################
for(i in 1:sim.reps){
  current.study.data <- GenerateStudyData(N=N, nobs.to.sample=nobs.to.sample, 
                                          randomization.method=randomization.method, 
                                          sigma.intercept=sigma.intercept, sigma.noise = sigma.noise,
                                          block.size=block.size, trt.probs.for.non.traj=trt.probs.for.non.traj,
                                          ordinal.breaks = ordinal.breaks,
                                          k.unif=k.unif, k.norm=k.norm,
                                          k.bin=k.bin, bin.props=bin.props,
                                          true.model.formula = true.model.formula, true.params, working.model.formula)
  gee.mod <- FitGEEModel(current.study.data, model.form = gee.model.formula)
  
  if(calc.oos.metrics.flag == TRUE){
    # Generate out of sample data for evaluating model performance
    oos.data <- GenerateOOSData(Noos = Noos, k.unif, k.norm, k.bin, bin.props,
                                true.model.formula, true.params)
    # Fit a treatment indicators only model for comparison
    trt.only.mod <- FitGEEModel(current.study.data, 
                                model.form = formula(Obsij ~ (Dwell + Music + Viz + Squeeze)^2))
    percentage.of.oracle.value[i,] <- CalcOOSValue(oos.data, working.model= gee.mod, 
                                                   trt.only.model= trt.only.mod,
                                                   model.formula = gee.model.formula)
  }
  

  
  if(calc.percentage.best.treat.flag == TRUE){
    # Calculate the proportion of patients in sample who receieved the best treatment for them 
    percentage.best.treatment[i] <- CalcISBestTreatment(current.study.data,
                                                        true.model.formula,
                                                        true.params)
    treatment.assignments[i,] <- eval(table.call)
  }

  # Different setups have different assignments of interest
  # With subtypes we're interested in the allocation by important subtype
  if(test.hypotheses.flag == TRUE){
    mixed.model.pvals.unadj[i,] <- CalcMMPval(FitMixedModel(current.study.data, working.model.formula), contrast.mat)
    gee.pvals.unadj[i,] <- CalcGEEPval(gee.mod, contrast.mat)
    
  }
}

end.time <- Sys.time()
#################################################
#### Write out simulation results
#################################################
if(save.results == TRUE){
  path.name <- paste0(getwd(), "/Results/",settings.type,"-", cur.time)
  # Create the results path
  dir.create(path.name)
  
  
  # Write out the parameters file used
  write.csv(param.df, file = paste0(path.name, "/parameter_values.csv"), 
            row.names = FALSE)
  
  if(test.hypotheses.flag == TRUE){
  # only write out the p-values if we did hypothesis tests
  # Write out the p-values
  write.csv(gee.pvals.unadj, file = paste0(path.name, "/unadjusted-gee-pvals.csv"), 
            row.names = FALSE)
  write.csv(mixed.model.pvals.unadj, file = paste0(path.name, "/unadjusted-mm-pvals.csv"), 
            row.names = FALSE)
  }
  # Write out the treatment assignment summary
  write.csv(treatment.assignments, file = paste0(path.name, "/treatment-assignment-summary.csv"), 
            row.names = FALSE)
  write.csv(percentage.best.treatment, file = paste0(path.name, "/percentage-best-treatment.csv"))
  # Write out the percentage of oracle value matrix
  write.csv(percentage.of.oracle.value, file = paste0(path.name, "/perc-oracle.csv"), 
            row.names = FALSE)
  
  
  # Write out the settings file used
  fileConn <- file(paste0(path.name,"/settings.txt"))
  settings.file <- c(readLines(settings.path),
                     paste("#Run time:", end.time-start.time))
  
  writeLines(settings.file, fileConn)
  close(fileConn)
  
}
