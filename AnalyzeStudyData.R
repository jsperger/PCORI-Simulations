#################################################
#### Data Analysis
#################################################
require(geepack)
require(rstan)
require(lme4)
rstan_options(auto_write = TRUE)
options(mc.cores = 4)
require(MASS)
select <- dplyr::select
#library(DynTxRegime)
#library(multcomp)

###########################################################
### Mixed Model
###########################################################


FitMixedModel <- function(observed.data, model.form = "Obsij ~ (Dwell + Music + Viz + Squeeze)^2 + (1 | ID)"){
  mixed.fit <- lmer(model.form, 
                    data = observed.data, 
                    na.action = na.omit)
  return(mixed.fit)

}

FitGLMM <- function(observed.data, model.form = formula(Obsij ~ (Dwell + Music + Viz + Squeeze)^2)){
  # Not working yet
  # Uses ordinal package
  # 
  current.study.data$Resp <- factor(current.study.data$Obsij, ordered = TRUE)
  
  fmm1 <- clmm(Resp ~ (Dwell + Music + Viz + Squeeze)^2 + (1 | ID), data = current.study.data, link = "logit")

  return(mixed.fit)
  
}

CalcMMPval <- function(mixed.fit, contrast.mat){
  
  sigma.err <- attr(VarCorr(mixed.fit), "sc")
  Sigma <- sigma.err*summary(mixed.fit)$vcov
  
  wald.stats <- vector(length = nrow(contrast.mat))
  wald.pvals <- vector(length = nrow(contrast.mat))
  beta <- fixef(mixed.fit)
  for(i in 1:length(wald.stats)){
    current.contrast <- matrix(contrast.mat[i,], nrow=1)
    wald.outer <- (current.contrast %*% beta)
    wald.inner <- solve(current.contrast %*% Sigma %*% t(current.contrast))
    wald.stats[i] <- t(wald.outer) %*% wald.inner %*% wald.outer
    wald.pvals[i] <- pchisq(q = wald.stats[i], df = 1, lower.tail = FALSE)
  }
  
  
  
  #Exclude intercept p-vals
  
  #mixed.wald.chisq <- summary(mixed.fit)$coefficients[-1,3]^2
  #mixed.pvals <- pchisq(mixed.wald.chisq, df = 1, lower.tail = FALSE)
  return(wald.pvals)
}

PredictMM <- function(mixed.fit, newdat){
  # Not a general prediction function
  # Specifically for an approximation to Thompson Sampling
  # Ignores random intercepts because all comparisons will be within-subject
  
  model.form <- formula(mixed.fit, fixed.only = TRUE)
  # Delete the response from the terms otherwise it will get angry that the new data doesn't have a response
  model.form.nr <- delete.response(terms(model.form))
  model.frame.nr <- model.frame(model.form.nr, as.data.frame(newdat))
  
  design.mat <- model.matrix(model.form.nr, model.frame.nr)

  beta <- fixef(mixed.fit)
  sigma.err <- attr(VarCorr(mixed.fit), "sc")
  Sigma <- sigma.err*summary(mixed.fit)$vcov
  sigma.b <- sqrt(c(summary(mixed.fit)$var[[1]]))

  parameter.draws <- mvrnorm(n = nrow(newdat), mu = beta, Sigma = Sigma)
  # Result of the matrix multiplication is the predicted value under every parameter draw
  # The diagonals are the prediction for the j-th observation from the j-th parameter draw
  pred.vals <- diag(design.mat %*% t(parameter.draws))
  
  #pred.vals <- vector(length = nrow(newdat))
  #for(i in 1:nrow(newdat)){
  #  parameter.draw <- mvrnorm(n = 1, mu = beta, Sigma = Sigma)
  #  pred.vals[i] <- design.mat[i,] %*% parameter.draw
  #}
  
  return(pred.vals)
}

###########################################################
### GEE
###########################################################
FitGEEModel <- function(observed.data, model.form = formula(Obsij ~ (Dwell + Music + Viz + Squeeze)^2)){
  gee.fit <- geeglm(formula = model.form, id=ID, data = observed.data, 
                    corstr = "exchangeable", na.action="na.omit")
  return(gee.fit)
}

CalcGEEPval <- function(gee.fit, contrast.mat){
  Sigma <- summary(gee.fit)$cov.scaled
  
  wald.stats <- vector(length = nrow(contrast.mat))
  wald.pvals <- vector(length = nrow(contrast.mat))
  beta <- coef(gee.fit)
  for(i in 1:length(wald.stats)){
    current.contrast <- matrix(contrast.mat[i,], nrow=1)
    wald.outer <- (current.contrast %*% beta)
    wald.inner <- solve(current.contrast %*% Sigma %*% t(current.contrast))
    wald.stats[i] <- t(wald.outer) %*% wald.inner %*% wald.outer
    wald.pvals[i] <- pchisq(q = wald.stats[i], df = 1, lower.tail = FALSE)
  }
  #Exclude intercept pvals
  #gee.pvals <- summary(gee.fit)$coefficients[-1,4]
  return(wald.pvals)
}

###########################################################
### Bayesian Model
###########################################################

FitStanModel <- function(observed.data, model.form, prior.int, prior.int.var, prior.beta, prior.beta.var){
  #No intercept because that's done separately not because there is no intercept in the model
  model.form.stan <- formula(Obsij ~ (Dwell + Music + Viz + Squeeze)^2 - 1)
  mf <- model.frame(model.form.stan, observed.data)
  design.mat <- model.matrix(model.form.stan, mf)
  stanDat <- list(M = nrow(observed.data),
                  N = length(unique(observed.data$ID)),
                  K = ncol(design.mat),
                  x = design.mat,
                  subj = observed.data$ID,
                  y = observed.data$Obsij, 
                  a = 4,
                  a_var = 1.5,
                  b = rep(0, length = ncol(design.mat)),
                  b_var = rep(2, length = ncol(design.mat)))
  riFit <- stan(file = "random_ints.stan",
                data = stanDat,
                iter = 1000,
                chains = 4)
  return(riFit)
}

###########################################################
### Test Out of Sample
###########################################################
GenerateOOSData <- function(Noos, k.unif, k.norm, k.bin, bin.props,
                            true.model.formula, true.params){
  oos.data <- inner_join(GenContinuousCovariates(n.subj = Noos,
                                                 d.unif=k.unif,
                                                 d.norm = k.norm),
                         GenBinaryCovariates(n.subj=Noos, d.bin = k.bin, props=bin.props), by = "ID")
  data.dupe <- oos.data[rep(1:nrow(oos.data),each=11),]
  data.dupe$Treatment <- rep(1:11, times = nrow(oos.data))
  data.dupe <- bind_cols(data.dupe, GenTreatmentIndicators(data.dupe))
  # No random intercepts and no random noise for out of sample data
  data.dupe <- GenOutcomes(data.dupe, true.model.formula=true.model.formula, 
                           true.params=true.params, oos.flag=TRUE)
  data.dupe$Visit <- 1
  return(data.dupe)
}


CalcOOSValue <- function(oos.data, working.model, trt.only.model, 
                         model.formula, 
                         trt.only.formula = formula(Obsij ~ (Dwell + Music + Viz + Squeeze)^2)){
  #Debugging
  #working.model <- gee.mod
  #trt.only.model <- trt.only.mod
  #model.formula <- gee.model.formula
  #trt.only.formula <- formula(Obsij ~ (Dwell + Music + Viz + Squeeze)^2)
  #
  model.X <- model.matrix(delete.response(terms(model.formula)), data=oos.data)
  trt.X <- model.matrix(delete.response(terms(trt.only.formula)), data=oos.data)
  if(all(colnames(model.X) == names(coef(gee.mod))) == FALSE){
    warning("Prediction model parameter mismatch")
  }
  oos.data$GPred <- c(model.X %*% coef(working.model))
  oos.data$TOPred <- c(trt.X %*% coef(trt.only.model))
  
  value.summary <- oos.data %>% group_by(ID) %>% 
    summarise(
      BestValNorm = min(Normij),
      BestValOrd = min(Obsij),
      BestValPred = min(GPred),
      BestTrt = which.min(Normij),
      BestTrtOrd = which.min(Obsij),
      BestTrtPred = which.min(GPred),
      TOPred = which.min(TOPred),
      OrdValfromPredTreat = Obsij[BestTrtPred],
      OrdValfromTOTreat = Obsij[TOPred],
      NormValfromPredTreat = Normij[BestTrtPred],
      ValueDelta = BestValOrd - Obsij[BestTrtPred],
      AbsDifNorm = abs(BestValNorm - Normij[BestTrtPred]),
      PerfectInd = ifelse(BestTrt ==BestTrtPred,1,0),
      SocNorm = Normij[1])
  percentage.of.oracle.value <- sum(value.summary$BestValOrd)/sum(value.summary$OrdValfromPredTreat)
  mse <- sum((value.summary$BestValNorm-value.summary$NormValfromPredTreat)^2)
  mse.soc <- sum((value.summary$BestValNorm-value.summary$SocNorm)^2)
  perc.improvement.to <- sum(value.summary$OrdValfromTOTreat)/sum(value.summary$OrdValfromPredTreat)
  misclass.rate <- sum(value.summary$PerfectInd)/nrow(value.summary)
  to.return <- c(percentage.of.oracle.value, mse, mse.soc,perc.improvement.to, misclass.rate)
  #names(to.return) <- c("PercOracle", "PercOverTrtOnly")
  return(to.return)
}

CalcOracleValue <- function(oos.data){
  # Expects OOS data generated by the GenerateOOSData function
  value.summary <- oos.data %>% group_by(ID) %>% 
    summarise(
      BestValNorm = min(TEij),
      BestValOrd = min(Obsij),
      BestTrt = which.min(TEij),
      BestTrtOrd = which.min(Obsij))
}

CalcISBestTreatment <- function(current.data, true.model.formula, true.params){
  current.data$ActualTreatment <- current.data$Treatment
  data.dupe <- current.data[rep(1:nrow(current.data),each=11),]
  data.dupe$Treatment <- rep(1:11, times = nrow(current.data))
  data.dupe <- data.dupe %>% select(-Dwell, -Music, -Viz, -Squeeze, -TEij, -Yij, -ObsNormij, -Obsij)
  data.dupe <- bind_cols(data.dupe, GenTreatmentIndicators(data.dupe))
  data.dupe <- GenOutcomes(data.dupe, true.model.formula, true.params)
  value.summary <- data.dupe %>% group_by(ID, Visit) %>% 
    summarise(
      BestTrt = which.min(TEij),
      GotBestTrt = ifelse(BestTrt ==ActualTreatment[1],1,0))
  return(mean(value.summary$GotBestTrt))
}

CalcISBestTreatmentNoSG <- function(){
  # Hard coded to match the parameters in the subgroup-1 setting
}