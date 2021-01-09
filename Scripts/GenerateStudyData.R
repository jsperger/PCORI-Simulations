##


################################################################################
## Data Generation Wrapper
##
################################################################################

CreateStudyData <- function(N, visit.proportions,
                            k.unif, k.norm, k.bin, bin.props,
                            sigma.intercept, sigma.noise,
                            treatment.arm.map, arm.var.name,
                            true.model.formula, true.params,
                            ordinal.breaks, first.cysto.prop,
                            ...){
  
  study_data <- GeneratePatientsAndVisitTotals(N = N, visit.proportions = visit.proportions) %>% 
    GenContinuousCovariates(indf = ., d.unif = k.unif, d.norm = k.norm) %>% 
    GenBinaryCovariates(indf = ., k.bin, props = bin.props) %>% 
    GenLatentVar(., sigma.intercept) %>% 
    GenAge(.) %>% 
    GenFirstCystoIndicator(., first.cysto.prop = first.cysto.prop) %>% 
    ElongateAndRandomizeStudyDataStratified(., possible.arms = treatment.arm.map$Arm) %>% 
    GenTreatmentIndicators(., treatment.arm.map, arm.var.name) %>% 
    GenOutcomes(study.data = ., 
                true.model.formula, 
                true.params,
                sigma.noise,
                ordinal.breaks)
  
  return(study_data)
}


################################################################################
## Covariate Generation Functions
##
################################################################################

#' Generate Continuous Covariates
#' Continuous covariates are either standard normal or uniform (-1,1)
#' @param d.unif is the number of uniform RVs to generate for each subject
#' @param d.norm is the numebr of standard normal RVs to generate for each subject
#' @reutrn Output is an N x (d.unif + d.norm) tibble
GenContinuousCovariates <- function(indf, d.unif, d.norm, ...){
  if (d.unif == 0 & d.norm == 0) return(indf)
  
  n.subj <- nrow(indf)
  
  U <- replicate(n = d.unif, runif(n=n.subj, min = -1, max = 1))
  colnames(U) <- paste0("U", 1:d.unif)
  
  U <- as_tibble(U)
  
  Nd <- replicate(n = d.norm, rnorm(n=n.subj, mean = 0, sd = 1))
  colnames(Nd) <- paste0("N", 1:d.norm)
  
  Nd <- as_tibble(Nd)
  
  X <- bind_cols(U, Nd)
  X$ID <- 1:n.subj
  
  df_with_covars <- indf %>% left_join(., y = X, by = "ID")
  
  return(df_with_covars)
}

#' Generate Binary Covariates
#' @param n.subj number of subjects to generate covariates for
#' @param d.bin the number of binary covariates to generate
#' @param props is a vector of length equal to d.bin or length 1 which 
#'     specifies the probability of a binary covariate being equal to 1
#' @return an n.subj x d.bin  tibble
GenBinaryCovariates <- function(indf, d.bin = k.bin, props = .5, ...){
  
  n.subj <- nrow(indf)
  
  if(length(props) != d.bin & length(props) != 1){
    warning("The number of proportions specified does not match the stated number of binary covariates")
  }
  
  # props determines what proportion of the subjects have a value of 1
  X <- props %>% map(rbinom, n = n.subj, size = 1)
  names(X) <- paste0("B", 1:d.bin)
  X <- bind_cols(X)
  X$ID <- 1:n.subj
  
  df_with_covars <- indf %>% left_join(., y = X, by = "ID")
  
  return(df_with_covars)
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
                        true.params,
                        sigma.noise,
                        ordinal.breaks){
  
  model.X <- model.matrix(true.model.formula, data=study.data)
  model.X <- model.X[, names(true.params)]
  # Check that the design matrix and the parameter names match
  if(all(colnames(model.X) == names(true.params)) == FALSE){
    warning("Model formula variable names do not match parameter names")
  }
  
  # Calculate the outcome with no noise
  study.data$Muij <- c(model.X %*% true.params)
  
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
  
  study_long <- study.data %>% group_by(!!!strata.vars.syms) %>% 
    mutate(Arm1 = sample(rep(possible.arms, 
                            each = ceiling(n()/length(possible.arms))),
                        size = n())) %>% 
    ungroup %>% 
    group_by(ID, TotalVisits, Arm1) %>% 
    nest %>%  
    mutate(Visit = list(1:TotalVisits)) %>% 
    unnest(., cols = Visit) %>% 
    group_by(ID) %>% 
    mutate(Arm = c(Arm1[[1]], sample(setdiff(possible.arms, Arm1[[1]]), 
                                size = TotalVisits-1))) %>% 
    unnest(., cols = data) %>% 
    select(-Arm1) %>% 
    mutate(FirstCysto = if_else(Visit >= 2, 0L, FirstCysto))
  
  
  return(study_long)
}