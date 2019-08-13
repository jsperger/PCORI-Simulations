#################################################
#### Run a simulation
#################################################

# load settings
source("./Settings/sim_settings_no_subgroups.R")
source("./GenerateStudyData.R")
source("./AnalyzeStudyData.R")

set.seed(seed)
cur.time <- Sys.time()

# Create the containers for where p-values and treatment assignments 
# from the simulation runs will be stored
mixed.model.pvals.unadj <- matrix(nrow = sim.reps, ncol = 10)
colnames(mixed.model.pvals.unadj) <- trt.names[-1]
gee.pvals.unadj <- matrix(nrow = sim.reps, ncol = 10)
colnames(gee.pvals.unadj) <- trt.names[-1]

treatment.assignments <- matrix(nrow = sim.reps, ncol = 11)
# Create the contrast matrix for the hypothesis tests
# Note: not actually a contrast matrix
# Each row is tested separately, not a joint hypothesis
temp.study <- tibble(Treatment = 1:11)
temp.study <- bind_cols(temp.study, GenTreatmentIndicators(temp.study))
treat.mat <- model.frame(true.model.formula, data = temp.study)
contrast.mat <- model.matrix(true.model.formula, treat.mat)[-1,]
contrast.mat[,1] <- 0
rm(temp.study, treat.mat)

#################################################
#### Run the Simulation
#################################################
for(i in 1:sim.reps){
  current.study.data <- GenerateStudyData(N=N, nobs.to.sample=nobs.to.sample, 
                                          randomization.method=randomization.method, 
                                          sigma.intercept=sigma.intercept, sigma.noise = sigma.noise,
                                          block.size=block.size, trt.probs.for.non.traj=trt.probs.for.non.traj,
                                          trt.effects=trt.effects, ordinal.breaks = ordinal.breaks)
  mixed.model.pvals.unadj[i,] <- CalcMMPval(FitMixedModel(current.study.data), contrast.mat)
  treatment.assignments[i,] <- table(current.study.data$Treatment)
  gee.pvals.unadj[i,] <- CalcGEEPval(FitGEEModel(current.study.data), contrast.mat)
}

#Bonferroni correction
adjusted.mm.pvals <- mixed.model.pvals.unadj*10
adjusted.gee.pvals <- gee.pvals.unadj*10

#################################################
#### Write out simulation results
#################################################
path.name <- paste0(getwd(), "/Results/",settings.type,"-", cur.time)
# Create the results path
dir.create(path.name)

# Write out the settings file used
fileConn <- file(paste0(path.name,"/settings.txt"))
settings.file <- readLines("./Settings/sim_settings_no_subgroups.R")
writeLines(settings.file, fileConn)
close(fileConn)

# Write out the p-values
write.csv(gee.pvals.unadj, file = paste0(path.name, "/unadjusted-gee-pvals.csv"), 
          row.names = FALSE)
write.csv(mixed.model.pvals.unadj, file = paste0(path.name, "/unadjusted-mm-pvals.csv"), 
          row.names = FALSE)
# Write out the treatment assignment summary
write.csv(treatment.assignments, file = paste0(path.name, "/treatment-assignment-summary.csv"), 
          row.names = FALSE)
