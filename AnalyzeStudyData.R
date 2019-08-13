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
### Extract p-values of interest
###########################################################

#K <- diag(length(coef(m1)))[-1,]
#rownames(K) <- names(coef(m1))[-1]
#test1 <-glht(m1, linfct=K)




#gee.adjusted.pvals <- gee.pvals*10
