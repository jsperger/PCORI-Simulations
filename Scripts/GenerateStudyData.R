##


################################################################################
## Data Generation Wrapper
##
################################################################################

CreateStudyData <- function(N, visit.proportions,
                            sigma.intercept, sigma.noise,
                            treatment.arm.map, arm.var.name,
                            true.model.formula, # true.params removed, will come from cluster_params_list
                            ordinal.breaks, first.cysto.prop,
                            n_clusters = 1,
                            cluster_proportions = c(1),
                            cluster_params_list = list(),
                            ...){
  
  # study_data <- GeneratePatientsAndVisitTotals(N = N, visit.proportions = visit.proportions) %>%
  #   GenContinuousCovariates(indf = ., d.unif = k.unif, d.norm = k.norm) %>% # k.unif, k.norm removed
  #   GenBinaryCovariates(indf = ., k.bin, props = bin.props) %>% # k.bin, bin.props removed
  #   GenLatentVar(., sigma.intercept) %>%

  # Initial patient and visit generation
  study_data <- GeneratePatientsAndVisitTotals(N = N, visit.proportions = visit.proportions)

  # Assign to clusters
  study_data <- AssignToClusters(indf = study_data, N_total = N,
                                 n_clusters = n_clusters,
                                 cluster_proportions = cluster_proportions)

  # Generate Covariates based on clusters
  study_data <- study_data %>%
    GenContinuousCovariates(indf = .,
                            n_clusters = n_clusters,
                            cluster_params_list = cluster_params_list) %>%
    GenBinaryCovariates(indf = .,
                        n_clusters = n_clusters,
                        cluster_params_list = cluster_params_list) %>%
    GenLatentVar(., sigma.intercept) %>% 
    GenAge(.) %>% 
    GenFirstCystoIndicator(., first.cysto.prop = first.cysto.prop) %>% 
    ElongateAndRandomizeStudyDataStratified(., possible.arms = treatment.arm.map$Arm) %>% 
    GenTreatmentIndicators(., treatment.arm.map, arm.var.name) %>% 
    GenOutcomes(study.data = ., 
                true.model.formula = true.model.formula,
                n_clusters = n_clusters,
                cluster_params_list = cluster_params_list,
                sigma.noise = sigma.noise,
                ordinal.breaks = ordinal.breaks)
  
  return(study_data)
}


################################################################################
## Cluster Assignment Function
################################################################################
AssignToClusters <- function(indf, N_total, n_clusters, cluster_proportions) {
  # N_total is used here instead of nrow(indf) in case indf is already grouped or modified
  if (n_clusters == 1) {
    indf$ClusterID <- 1
  } else {
    if (length(cluster_proportions) != n_clusters) {
      stop("Length of cluster_proportions must be equal to n_clusters.")
    }
    if (sum(cluster_proportions) != 1) {
      stop("cluster_proportions must sum to 1.")
    }
    # Ensure that IDs are unique for sampling if indf contains only unique IDs
    # If indf is just the initial ID generation, N_total should match nrow(indf)
    indf$ClusterID <- sample(1:n_clusters, size = N_total, replace = TRUE, prob = cluster_proportions)
  }
  return(indf)
}

################################################################################
## Covariate Generation Functions
##
################################################################################

#' Generate Continuous Covariates
#' Continuous covariates are either standard normal or uniform (-1,1)
#' @param d.unif is the number of uniform RVs to generate for each subject (now per cluster)
#' @param d.norm is the number of standard normal RVs to generate for each subject (now per cluster)
#' @return Output is an N x (total d.unif + total d.norm) tibble, with covariates generated cluster-specifically
GenContinuousCovariates <- function(indf, n_clusters, cluster_params_list, ...){
  
  # Helper to get parameter value or return a default
  get_param_value <- function(params_df, param_name, default_value) {
    val <- params_df$Value[params_df$Parameter == param_name]
    if (length(val) == 0 || is.na(val)) {
      return(default_value)
    }
    return(as.numeric(val))
  }

  all_covars_list <- list() # To store covariates for each subject

  # Determine max k.unif and k.norm across clusters for consistent column structure if needed,
  # though current approach generates them per cluster and names them U1, U2.. N1, N2..
  # This means if cluster 1 has U1 and cluster 2 has U1, they are different variables.
  # The problem description implies U1, N1 are specific (e.g. mean.U1), so this should be fine.

  for (cl_id in 1:n_clusters) {
    cluster_data <- indf %>% filter(ClusterID == cl_id)
    if (nrow(cluster_data) == 0) next # Skip if no subjects in this cluster

    params_for_cluster <- cluster_params_list[[cl_id]]

    # Get k.unif and k.norm for this cluster (assuming they might not exist, default to 0)
    # These specific parameter names "k.unif" and "k.norm" might not be in cluster_params.csv yet.
    # The prompt implies mean.U1, sd.U1, so we assume k.unif=1 if mean.U1 is present, etc.
    # For now, let's look for specific U covariates like U1, U2 and N covariates N1, N2.

    # Max number of uniform/normal covariates mentioned in any cluster file (e.g. up to U<max_k_unif>)
    # This needs to be more robust by checking which specific Ux or Nx params are defined.
    # Example: check for "mean.U1", "sd.U1", then "mean.U2", "sd.U2" etc.

    # Let's assume for now that the cluster CSVs will define parameters for U1, N1 if they exist.
    # And that k.unif, k.norm in the *original* settings dictate the *number* of such variables globally.
    # This part of the prompt is a bit ambiguous. Let's re-read:
    # "For cluster 1, use values similar to the existing global settings (e.g., k.unif <- 1 and k.norm <- 1 imply U1 and N1)"
    # This suggests k.unif/k.norm are still somewhat global in defining *how many* U and N variables there are.
    # Let's assume the settings file still has global k.unif and k.norm for # of covariates,
    # but their distributions (mean/sd/props) are cluster specific.
    # This means CreateStudyData needs to pass k.unif and k.norm to this function.
    # The prompt said to remove them from CreateStudyData's signature. This is contradictory.
    #
    # RESOLUTION: I will assume that the *number* of covariates (e.g., one U, one N) is fixed globally for now
    # as suggested by "k.unif <- 1 and k.norm <- 1 imply U1 and N1".
    # The cluster CSVs will then provide parameters for these fixed covariates (U1, N1).
    # If a cluster CSV does *not* provide params for U1, it defaults.

    # Let's assume global_k_unif and global_k_norm are passed or accessible.
    # For now, I'll hardcode to 1 U and 1 N as per example, assuming U1 and N1.
    # This needs clarification if more Us or Ns are expected.

    current_k_unif = 1 # Assuming one uniform covariate U1 based on example
    current_k_norm = 1 # Assuming one normal covariate N1 based on example

    cluster_covars <- tibble(ID = cluster_data$ID)

    if (current_k_unif > 0) {
      for (i in 1:current_k_unif) {
        mean_u_val <- get_param_value(params_for_cluster, paste0("mean.U", i), 0) # U(-1,1) has mean 0
        sd_u_val <- get_param_value(params_for_cluster, paste0("sd.U", i), sqrt(1/3)) # U(-1,1) has var 1/3 -> sd sqrt(1/3) ~ 0.577
        # To generate from a specific mean/sd for uniform, it's easier to scale a U(0,1) or U(-1,1)
        # runif(n, min, max) -> mean (min+max)/2, var (max-min)^2/12
        # Given mean M and sd S for U(a,b): M=(a+b)/2, S^2=(b-a)^2/12.
        # b-a = sqrt(12)*S. a+b = 2M.  2b = 2M+sqrt(12)S -> b = M + sqrt(3)S. 2a = 2M-sqrt(12)S -> a = M - sqrt(3)S
        min_val = mean_u_val - sqrt(3)*sd_u_val
        max_val = mean_u_val + sqrt(3)*sd_u_val
        U_vals <- runif(n=nrow(cluster_data), min = min_val, max = max_val)
        cluster_covars[[paste0("U", i)]] <- U_vals
      }
    }

    if (current_k_norm > 0) {
      for (i in 1:current_k_norm) {
        mean_n_val <- get_param_value(params_for_cluster, paste0("mean.N", i), 0)
        sd_n_val <- get_param_value(params_for_cluster, paste0("sd.N", i), 1)
        Nd_vals <- rnorm(n=nrow(cluster_data), mean = mean_n_val, sd = sd_n_val)
        cluster_covars[[paste0("N", i)]] <- Nd_vals
      }
    }
    all_covars_list[[cl_id]] <- cluster_covars
  }
  
  # Combine all generated covariates
  if (length(all_covars_list) > 0) {
    all_generated_covars <- bind_rows(all_covars_list)
    # Need to handle cases where some clusters might not generate all covariates if k.unif/k.norm differs by cluster.
    # Current assumption: k.unif/k.norm are fixed globally for now.
    indf <- indf %>% left_join(all_generated_covars, by = "ID")
  }
  
  return(indf)
}

#' Generate Binary Covariates
#' @param n.subj number of subjects to generate covariates for
#' @param d.bin the number of binary covariates to generate
#' @param props is a vector of length equal to d.bin or length 1 which (now per-cluster, e.g. bin.prop.B1)
#'     specifies the probability of a binary covariate being equal to 1
#' @return an n.subj x d.bin tibble
GenBinaryCovariates <- function(indf, n_clusters, cluster_params_list, ...){
  
  # Helper to get parameter value or return a default (can be defined globally or passed)
  get_param_value <- function(params_df, param_name, default_value) {
    val <- params_df$Value[params_df$Parameter == param_name]
    if (length(val) == 0 || is.na(val)) {
      return(default_value)
    }
    # For props, it could be a string "c(0.5, 0.3)" or single numeric.
    # For now, assume individual props like bin.prop.B1, bin.prop.B2 are separate rows.
    return(as.numeric(val))
  }

  all_covars_list <- list()

  # Assume global k.bin dictates the number of binary covariates, e.g., B1, B2.
  # This should ideally come from a global setting or be inferred.
  # From problem: "bin.props <- c(.35, .4) means bin.prop.B1 could be 0.35 and bin.prop.B2 could be 0.4"
  # This suggests k.bin=2 is known.
  global_k_bin <- 2 # Assuming two binary covariates: B1, B2 based on prompt.

  for (cl_id in 1:n_clusters) {
    cluster_data <- indf %>% filter(ClusterID == cl_id)
    if (nrow(cluster_data) == 0) next

    params_for_cluster <- cluster_params_list[[cl_id]]
    n_subj_cluster <- nrow(cluster_data)

    cluster_bin_covars <- tibble(ID = cluster_data$ID)

    if (global_k_bin > 0) {
      for (i in 1:global_k_bin) {
        prop_val <- get_param_value(params_for_cluster, paste0("bin.prop.B", i), 0.5) # Default to 0.5

        B_vals <- rbinom(n = n_subj_cluster, size = 1, prob = prop_val)
        cluster_bin_covars[[paste0("B", i)]] <- B_vals
      }
    }
    all_covars_list[[cl_id]] <- cluster_bin_covars
  }

  if (length(all_covars_list) > 0) {
    all_generated_covars <- bind_rows(all_covars_list)
    # Ensure columns are correctly aligned if some clusters define fewer binary vars (not current assumption)
    indf <- indf %>% left_join(all_generated_covars, by = "ID")
  }
  
  return(indf)
}

#' Generate First Cystoscopy time varying covariate
#' @param indf.long Long 
#' @param d.bin the number of binary covariates to generate
#' @param props is a vector of length equal to d.bin or length 1 which 
#'     specifies the probability of a binary covariate being equal to 1
#' @return an n.subj x d.bin  tibble
GenFirstCystoscopyVar <- function(indf.long, first.cysto.prop = .5, ...){
  
  
  
  return(df_with_covars)
}

#' Generate the random intercepts which are used to induce correlation across 
#' observations. This is equivalent to an excahngeable correlation structure
#' @param indf (wide) data frame. The data frame should be wide to ensure that the 
#' random intercept for a patient is constant across visits. 
#' @param sigma.intercept standard deviation of the random intercepts
GenLatentVar <- function(indf, sigma.intercept){
  return(indf %>% mutate(Ui = rnorm(n = n(), sd = sigma.intercept)))
}


#' Generate treatment indicators - TRUE if a treatment was received 
#' FALSE otherwise
#' @param study.data (long) data frame containing a column with the arm designation
GenTreatmentIndicators <- function(study.data,
                                   treatment.arm.map,
                                   arm.var.name = "Arm"){
  
  study.data <- left_join(study.data, treatment.arm.map, by = arm.var.name)
  return(study.data)
}

#' Generate an indicator for first cystoscopy
#' 
GenFirstCystoIndicator <- function(study.data, first.cysto.prop){
 study.data <-  study.data %>% 
    group_by(ID) %>% 
    nest %>% 
    mutate(FirstCysto = rbinom(n = n(), size = 1, prob = parent.frame()$first.cysto.prop)) %>% 
   unnest(cols = c(data)) #%>% 
   #mutate(FirstCysto = if_else(Visit >= 2, 0L, FirstCysto))
 
 return(study.data)
}

#' Generate the age covariate. Age is scaled so that the ages lies in (0,1)
#' The actual generation is done using the beta distribution
#' 
GenAge <- function(study.data){
  study.data <- study.data %>% mutate(Age = rbeta(n = n(), 5, 3))
  
  return(study.data)
}


################################################################################
## Outcome Generation Functions
##
################################################################################

GenOutcomes <- function(study.data, 
                        true.model.formula, 
                        n_clusters, # Added
                        cluster_params_list, # Added
                        sigma.noise,
                        ordinal.breaks,
                        ...){ # true.params removed

  study.data$Muij <- NA # Initialize Muij column
  
  # Ensure study.data has an original row order identifier if not already present, to reassemble Muij correctly
  # However, if we are modifying study.data by cluster subsets, then combining results, order should be preserved by ID.
  # Let's assume study.data has unique IDs and ClusterID is already merged.

  all_mu_ij_list <- list()

  for (cl_id in 1:n_clusters) {
    cluster_data_subset <- study.data %>% filter(ClusterID == cl_id)
    if (nrow(cluster_data_subset) == 0) next

    params_for_cluster_df <- cluster_params_list[[cl_id]]

    # Extract true.params for this cluster:
    # These are parameters like "Dwell", "Music", "(Intercept)" etc.
    # They are NOT prefixed like bin.prop.B1 or mean.U1
    # So, we need to filter params_for_cluster_df for rows that are coefficient names.
    # This assumes coefficient names in the CSV do not clash with other param names like "mean.U1".
    # A safer way would be to have a specific section or prefix for true.params in CSV, or a pre-defined list of coeff names.
    # For now, assume true.params are all rows that are NOT "bin.prop.B1", "bin.prop.B2", "mean.U1", etc.
    # This is fragile. A better approach: list expected coefficient names from true.model.formula.

    # Generate model matrix for this subset
    model.X_cluster <- model.matrix(true.model.formula, data = cluster_data_subset)

    # Extract true.params as a named vector for this cluster
    # The cluster_params.csv has "Parameter" and "Value" columns.
    # Need to select rows that are actual coefficients.
    # A simple heuristic: they don't contain "." like "bin.prop.B1" or "mean.U1".
    # Or, more robustly, they are present as column names in model.X_cluster.

    expected_coeffs <- colnames(model.X_cluster)
    cluster_true_params_df <- params_for_cluster_df %>% filter(Parameter %in% expected_coeffs)

    if(nrow(cluster_true_params_df) == 0 && cl_id == 1 && exists("default_param.df")) {
        # Fallback for single cluster scenario or if cluster1 specifically uses default
        # This logic might be better placed in RunSimulation.R when populating cluster_params_list
        warning(paste("No specific true.params found for Cluster", cl_id, "in its CSV. Trying to use default_param.df if available."))
        # default_param.df should have Parameter and Coefficient columns
        # The main settings file originally loaded into 'default_param.df' with columns 'Parameter' and 'Coefficient'
        # The new cluster CSVs use 'Parameter' and 'Value'. This is an inconsistency.
        # For now, I'll assume cluster_params_list[[cl_id]] always has 'Parameter' and 'Value'.
        # And that RunSimulation.R ensures default_param.df is presented in this list with these col names if used as fallback.

        # Let's assume cluster_params_list[[cl_id]] for default IS default_param.df, but need to handle column name diff.
        # This part is tricky due to differing column names (Coefficient vs Value).
        # For now, I'll assume cluster_params_list elements *always* have Parameter, Value.
        # The CSVs for clusters were created with Parameter,Value.
        # default_param.df (from full.csv) had Parameter,Coefficient.
        # This needs to be harmonized, probably when cluster_params_list is created in RunSimulation.R.
        # For this function, I will strictly expect Parameter, Value in cluster_params_list[[cl_id]].
    }

    cluster_true_params <- setNames(cluster_true_params_df$Value, cluster_true_params_df$Parameter)

    # Ensure order of params matches model.X_cluster columns
    cluster_true_params <- cluster_true_params[expected_coeffs]
    # Handle NA for any missing params - should ideally default or error
    cluster_true_params[is.na(cluster_true_params)] <- 0 # Default to 0 if a param is missing, or error.

    if(length(cluster_true_params) != ncol(model.X_cluster) || !all(names(cluster_true_params) == colnames(model.X_cluster))) {
      warning(paste("Mismatch or missing true.params for cluster", cl_id, ". Check CSV and formula. Using 0 for missing effects."))
      # Ensure it's a correctly named vector of the right length
      temp_params <- setNames(rep(0, ncol(model.X_cluster)), colnames(model.X_cluster))
      # Fill with what we have
      common_names <- intersect(names(temp_params), names(cluster_true_params))
      temp_params[common_names] <- cluster_true_params[common_names]
      cluster_true_params <- temp_params
    }

    Muij_cluster <- c(model.X_cluster %*% cluster_true_params)

    # Store Muij with original IDs to ensure correct merging later
    all_mu_ij_list[[cl_id]] <- tibble(ID = cluster_data_subset$ID, Muij_cluster = Muij_cluster, Visit = cluster_data_subset$Visit)
  }

  # Combine all Muij values. Need to handle cases of multiple visits per ID.
  # The study.data is long, so ID+Visit should be unique key for merging.
  if (length(all_mu_ij_list) > 0) {
    all_mu_ij_df <- bind_rows(all_mu_ij_list)
    # This merge might be tricky if study.data doesn't have unique ID-Visit combination before this step.
    # Assuming study.data is already elongated and contains ID and Visit columns.
    study.data <- study.data %>%
      left_join(all_mu_ij_df, by = c("ID", "Visit")) %>%
      mutate(Muij = ifelse(!is.na(Muij_cluster), Muij_cluster, Muij)) %>% # Replace NA Muij with calculated
      select(-Muij_cluster) # Clean up temporary column
  } else {
    # This case should not happen if there's data
    study.data$Muij <- 0 # Fallback if no Muij calculated
  }
  
  study.data <- study.data %>% 
    mutate(Errij = rnorm(n = n(), mean = 0, sd = sigma.noise),
           Zij = Muij + Ui + Errij,
           Yij = .GenOutcomesOrdinalize(Zij, cutoffs = ordinal.breaks))
}

#' Calculate the observed ordinal outcome on a 0-10 based on the latent 
#' normal response and the cutoffs (quantiles of the standard normal distribution)
#' @param invec vector of 
#' @param cutoffs length 11 vector of quantiles. The first element shoudl be -Inf and the
#' last element should be Inf
.GenOutcomesOrdinalize <- function(invec,
                       cutoffs){
  # as.integer will return 1-11 instead of 0-10
  ordinal_vec <- as.integer(cut(invec, breaks = cutoffs, include.lowest = TRUE,
                                labels = 0:10)) - 1
  return(ordinal_vec)
}

################################################################################
## Patient Generation Functions
##
################################################################################
#' Generate a N x 2 shell with an ID column and a column indicating the total
#' number of visits the participant has over the study
#' @param N
#' @param visit.proportions
GeneratePatientsAndVisitTotals <- function(N, visit.proportions){
  stopifnot(sum(visit.proportions) == 1)
  
  
  visit_number_set <- map2(.x = 1:length(visit.proportions),
                               .y = ceiling(N*visit.proportions), 
                               ~rep(.x, length.out = .y)) %>% 
    unlist
  
  study_data <- tibble(ID = 1:N) %>% 
    mutate(TotalVisits = sample(visit_number_set,
                                size = n(),
                                replace = FALSE))
  
  return(study_data)
}

#' Create a long data frame and randomize patients
ElongateAndRandomizeStudyData <- function(study.data, possible.arms){
  study_long <- study.data %>% 
    group_by(ID, TotalVisits) %>% 
    nest %>%  
    mutate(Visit = list(1:TotalVisits)) %>% 
    unnest(., cols = Visit) %>% 
    group_by(ID) %>% 
    mutate(Arm = sample(possible.arms, size = TotalVisits)) %>% 
    unnest(., cols = data) 
  
  return(study_long)
}


#' Create a long data frame and randomize patients
ElongateAndRandomizeStudyDataStratified <- function(study.data, possible.arms,
                                                    strata.vars.syms = c(sym("B1"), sym("FirstCysto"))){
  
  study_long <- study.data %>%
    group_by(!!!strata.vars.syms) %>%
    mutate(Arm1 = sample(rep(possible.arms, each = ceiling(n()/length(possible.arms))), size = n())) %>%
    ungroup() %>%
    group_by(ID, TotalVisits, Arm1) %>% # Arm1 is the first assigned arm for this ID. TotalVisits is also fixed per ID.
    nest() %>%  # Nests all other columns (B1, FirstCysto, U1, N1, ClusterID etc.) into 'data'
    mutate(Visit = map(TotalVisits, ~1:.x)) %>% # Create list of visits 1:TotalVisits for each ID
    unnest(cols = c(Visit)) %>% # Expand rows for each visit
    group_by(ID) %>% # Group by ID to assign full arm sequence per patient
    mutate(Arm = {
      # Within each ID group, Arm1 and TotalVisits are effectively scalar
      # (they are repeated for each visit, but Arm1[[1]] and TotalVisits[[1]] give the unique value for that ID)
      current_id_arm1 <- Arm1[[1]]
      current_id_total_visits <- TotalVisits[[1]]

      # Calculate number of follow-up visits. Default to 0 if TotalVisits is NA or less than 1.
      num_follow_up_visits <- ifelse(is.na(current_id_total_visits) || current_id_total_visits < 1,
                                     0L,
                                     as.integer(current_id_total_visits) - 1L)

      follow_up_arms <- if (num_follow_up_visits > 0) {
        # Arms available for follow-up visits (all arms except the first one assigned)
        available_follow_up_arms <- setdiff(possible.arms, current_id_arm1)

        # Ensure sample size doesn't exceed available arms for sampling without replacement
        # (should be 3 available arms, num_follow_up_visits should be <=3)
        sample_size <- min(num_follow_up_visits, length(available_follow_up_arms))

        if (sample_size > 0) {
          sample(available_follow_up_arms, size = sample_size, replace = FALSE)
        } else {
          integer(0) # No follow-up arms to sample
        }
      } else {
        integer(0) # No follow-up visits
      }

      # Combine the first arm with the sampled follow-up arms
      # This vector will be recycled by dplyr to fill all visit rows for this ID.
      c(current_id_arm1, follow_up_arms)
    }) %>%
    unnest(cols = data) %>% # Unnest the original characteristics (B1, FirstCysto, covariates, etc.)
    select(-Arm1) %>% # Remove the temporary Arm1 column
    mutate(FirstCysto = if_else(Visit >= 2, 0L, FirstCysto)) # Ensure FirstCysto is 0 for Visit >= 2
  
  return(study_long)
}