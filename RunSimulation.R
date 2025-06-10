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

# Attempt to load OpenBLAS, hopefully FlexiBLAS can manage thread settings.
# If "OPENBLAS" itself isn't found, this will error.
# If "OPENBLAS-SERIAL" was specific and this doesn't work,
# we might need to comment out flexiblas usage for testing in this environment.
serial_backend <- tryCatch({
  flexiblas_load_backend("OPENBLAS")
}, error = function(e) {
  warning("OPENBLAS backend not found, will try default. Error: ", e$message)
  return(NULL)
})

if (!is.null(serial_backend)) {
  flexiblas_switch(serial_backend)
} else {
  warning("FlexiBLAS backend 'OPENBLAS' could not be loaded. Using R's default BLAS.")
}

# I think this is redundant with setting the backend to serial, but just in case
flexiblas_set_num_threads(1)


settings.path <- "./Settings/full_scenario.R"
# load settings
source(settings.path)

# Check for n_clusters and load cluster-specific params
if (!exists("n_clusters")) {
  n_clusters <- 1
}

if (n_clusters > 1) {
  cluster_params_list <- list()
  for (i in 1:n_clusters) {
    cluster_file_path <- paste0("./Settings/cluster", i, "_params.csv")
    if (file.exists(cluster_file_path)) {
      cluster_params_list[[i]] <- read.csv(cluster_file_path)
    } else {
      # Fallback to default if a specific cluster file is missing
      warning(paste("Cluster parameter file not found:", cluster_file_path, ". Using default_param.df."))
      temp_df <- default_param.df
      if ("Coefficient" %in% names(temp_df)) {
        names(temp_df)[names(temp_df) == "Coefficient"] <- "Value"
      }
      cluster_params_list[[i]] <- temp_df
    }
  }
} else {
  # If n_clusters is 1 or not defined, use default_param.df
  # Ensure default_param.df is available (it should be loaded from the settings file)
  if (!exists("default_param.df") && exists("param.df")) {
    # If the old param.df exists (e.g. from older settings files not yet updated)
    default_param.df <- param.df
  }

  # Make sure default_param.df uses "Value" column
  if (exists("default_param.df")) {
    if ("Coefficient" %in% names(default_param.df)) {
      names(default_param.df)[names(default_param.df) == "Coefficient"] <- "Value"
    }
    cluster_params_list <- list(default_param.df)
  } else {
    # This case should ideally not be reached if settings are correct
    warning("default_param.df not found. Initializing cluster_params_list with an empty list or default structure.")
    # Create a minimal structure to avoid errors downstream, though this indicates a setup problem.
    cluster_params_list <- list(data.frame(Parameter = character(), Value = numeric()))
  }
}

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
                       true.model.formula = true.model.formula,
                       true.params = true.params, # This will likely be overridden by cluster_params_list_sim
                       first.cysto.prop = first.cysto.prop,
                       n_clusters_sim = 1,
                       cluster_proportions_sim = c(1),
                       cluster_params_list_sim = list()) {
  
  current.study.data <- CreateStudyData(N=N, 
                                        visit.proportions = visit.proportions,
                                        sigma.intercept=sigma.intercept, sigma.noise = sigma.noise,
                                        treatment.arm.map = treatment.arm.map,
                                        arm.var.name = arm.var.name, # Ensure correct variable is passed
                                        ordinal.breaks = ordinal.breaks,
                                        k.unif=k.unif, k.norm=k.norm,
                                        k.bin=k.bin, bin.props=bin.props,
                                        true.model.formula = true.model.formula, 
                                        true.params = true.params, # Fallback, might be overridden by cluster logic in CreateStudyData
                                        first.cysto.prop = first.cysto.prop,
                                        n_clusters = n_clusters_sim,
                                        cluster_proportions = cluster_proportions_sim,
                                        cluster_params_list = cluster_params_list_sim)
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
                       true.model.formula = true.model.formula,
                       true.params = true.params, # Pass the global true.params as a default/fallback
                       first.cysto.prop = first.cysto.prop,
                       n_clusters = n_clusters, # Pass the loaded n_clusters
                       cluster_proportions = cluster_proportions, # Pass the loaded cluster_proportions
                       cluster_params_list = cluster_params_list, # Pass the loaded list of cluster parameters
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
  write.csv(default_param.df, file = paste0(path.name, "/default_parameter_values.csv"),
            row.names = FALSE)
  
  # TODO: Consider saving cluster_params_list as well, perhaps as multiple files or an RDS object

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
