#################################################
#### Run a simulation
#################################################

# load settings
source("./Settings/sim_settings_no_subgroups.R")
source("./GenerateStudyData.R")
source("./AnalyzeStudyData.R")

mixed.model.pvals <- matrix(nrow = sim.reps, ncol = 10)
colnames(mixed.model.pvals) <- trt.names[-1]
gee.pvals <- matrix(nrow = sim.reps, ncol = 10)
colnames(gee.pvals) <- trt.names[-1]

for(i in 1:sim.reps){
  current.study.data <- GenerateStudyData(N=N, nobs.to.sample=nobs.to.sample, 
                                          randomization.method=randomization.method, 
                                          sigma.intercept=sigma.intercept, sigma.noise = sigma.noise,
                                          block.size=block.size, trt.probs.for.non.traj=trt.probs.for.non.traj,
                                          trt.effects=trt.effects, ordinal.breaks = ordinal.breaks)
  mixed.model.pvals[i,] <- FitMixedModel(current.study.data)
  gee.pvals[i,] <- FitGEEModel(current.study.data)
}


#Bonferroni correction
adjusted.mm.pvals <- mixed.model.pvals*10
adjusted.gee.pvals <- gee.pvals*10
