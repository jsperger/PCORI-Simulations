#################################################
#### Run a simulation
#################################################

library(tidyverse)
library(geepack)
library(parallel)
library(tictoc)
library(flexiblas)


# My knowledge of parallelization is a  little spotty, but we don't want the
# forked R processes spinning up multiple threads

serial_backend <- flexiblas_load_backend("OPENBLAS-SERIAL")
flexiblas_switch(serial_backend)

# I think this is redundant with setting the backend to serial, but just in case
flexiblas_set_num_threads(1)


settings.path <- "./Settings/full_scenario.R"
# load settings
source(settings.path)
source("./Scripts/GenerateStudyData.R")
source("./Scripts/AnalyzeStudyData.R")
source("./Scripts/AnalyzeSimulationRuns.R")
source("./Scripts/utilities.R")

set.seed(sim.seed)
cur.time <- gsub(":", "_", x=Sys.time() ,fixed=TRUE)

start.time <- Sys.time()
tic()

#################################################
#### Run the Simulation
#################################################
simWrapper <- function(i,  N, 
                       visit.proportions,
                       sigma.intercept, sigma.noise,
                       treatment.arm.map ,
                       arm.var.name = "Arm",
                       ordinal.breaks,
                       k.unif, 
                       k.norm,
                       k.bin, 
                       bin.props,
                       true.model.formula = true.model.formula, true.params = true.params,
                       first.cysto.prop = first.cysto.prop) {
  
  current.study.data <- CreateStudyData(N=N, 
                                        visit.proportions = visit.proportions,
                                        sigma.intercept=sigma.intercept, sigma.noise = sigma.noise,
                                        treatment.arm.map = treatment.arm.map,
                                        arm.var.name = "Arm",
                                        ordinal.breaks = ordinal.breaks,
                                        k.unif=k.unif, k.norm=k.norm,
                                        k.bin=k.bin, bin.props=bin.props,
                                        true.model.formula = true.model.formula, 
                                        true.params = true.params,
                                        first.cysto.prop = first.cysto.prop)
  #g1 <- FitGEEModel(current.study.data, gee.model.formula)
 
    return(current.study.data)
  
}


tic()
  
study.list <- parallel::mclapply(X = 1:sim.reps, FUN = simWrapper,
                                    N=N, 
                       visit.proportions = visit.proportions,
                       sigma.intercept=sigma.intercept, sigma.noise = sigma.noise,
                       treatment.arm.map = treatment.arm.map,
                       arm.var.name = "Arm",
                       ordinal.breaks = ordinal.breaks,
                       k.unif=k.unif, k.norm=k.norm,
                       k.bin=k.bin, bin.props=bin.props,
                       true.model.formula = true.model.formula, true.params = true.params,
                       first.cysto.prop = first.cysto.prop,
                       mc.cores = 8)

gee.mods <- parallel::mclapply(X = study.list, 
                               FUN = FitGEEModel,
                               model.form = gee.model.formula,
                               mc.cores = 8)

gee.pvals.unadj <- parallel::mclapply(X = gee.mods, FUN = CalcGEEPval,
                                      contrast.mat = contrast.mat,
                                      mc.cores = 8) %>% 
  bind_rows(.)



toc()

gee.pvals.adj <- AdjustPvals(gee.pvals.unadj)

gee_coefs <- map_dfr(gee.mods, coef)
colMeans(gee_coefs)


if (null.scenario == TRUE) print(paste0("FWER: ", round(mean(rowSums(gee.pvals.adj < .05) >= 1), 3)))


if (null.scenario == FALSE) print(colMeans(gee.pvals.adj < .05))


end.time <- Sys.time()

#################################################
#### Write out simulation results
#################################################
if(save.results == TRUE){
  path.name <- paste0(getwd(), "/Results/",settings.type,"-", cur.time)
  # Create the results path
  dir.create(path.name,
             recursive = TRUE,
             showWarnings = TRUE)
  
  
  # Write out the parameters file used
  write.csv(param.df, file = paste0(path.name, "/parameter_values.csv"), 
            row.names = FALSE)
  
  write_csv(gee.pvals.unadj,
          file = paste0(path.name, "/gee_pvals.csv"))
  
  if(save.sim.study.data == TRUE){
    saveRDS(study.list,
            file = paste0(path.name, "/study_data_list.rds"))

  }
  

  # Write out the settings file used
  fileConn <- file(paste0(path.name,"/settings.txt"))
  settings.file <- c(readLines(settings.path),
                     paste("#Run time:", end.time-start.time))
  
  writeLines(settings.file, fileConn)
  close(fileConn)
  
  # Write out the script used to run the simulations
  fileConn <- file(paste0(path.name,"/run_script.txt"))
  settings.file <- c(readLines("RunSimulation.R"),
                     paste("#Run time:", end.time-start.time))
  
  writeLines(settings.file, fileConn)
  close(fileConn)
  
}
