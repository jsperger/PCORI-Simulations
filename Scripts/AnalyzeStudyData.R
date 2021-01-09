# Functions for Analyzing a Single Study
# Author: John Sperger
# Last modified on 04 December 2020


################################################################################
## Mixed Model Functions
##
################################################################################
#' Wrapper around the lmer function from the \code{lme4} package
FitLMM <- function(study.data, 
                   model.string = "Yij ~ (Dwell + Music)^2 + (1 | ID)"){
  
  mixed.fit <- lme4::lmer(model.string, 
                          REML = FALSE,
                          data = study.data, 
                          na.action = na.omit)
  
  return(mixed.fit)
}


FitGLMM <- function(observed.data, 
                    model.form = formula(Obsij ~ (Dwell + Music + Viz + Squeeze)^2)){
  # Not working yet
  # Uses ordinal package
  # 
  current.study.data$Ordij <- factor(current.study.data$Yij, ordered = TRUE)
  
  fmm1 <- ordinal::clmm(Ordij ~ (Dwell + Music)^2 + (1 | ID), 
                        data = study.data, link = "logit")
  
  return(mixed.fit)
}

#' Calculate Wald p-values for 
#' @param mixed.fit fitted model from \code{lmer}
#' @param 
CalcMMPval <- function(mixed.fit, contrast.mat){
  
  sigma.err <- attr(lme4::VarCorr(mixed.fit), "sc")
  Sigma <- sigma.err*summary(mixed.fit)$vcov
  
  wald.stats <- vector(length = nrow(contrast.mat))
  wald.pvals <- vector(length = nrow(contrast.mat))
  
  beta_hat <- lme4::fixef(mixed.fit)
  
  for(i in 1:length(wald.stats)){
    current.contrast <- matrix(contrast.mat[i,], nrow=1)
    wald.outer <- (current.contrast %*% beta_hat)
    wald.inner <- solve(current.contrast %*% Sigma %*% t(current.contrast))
    wald.stats[i] <- t(wald.outer) %*% wald.inner %*% wald.outer
    wald.pvals[i] <- pchisq(q = wald.stats[i], df = 1, lower.tail = FALSE)
  }
  
  names(wald.pvals) <- rownames(contrast.mat)
  

  return(wald.pvals)
}



################################################################################
## GEE Functions
##
################################################################################

#' Fit a linear model using GEE with an exchangeable correlation structure
#' Essentially a wrapper around \code{geeglm} from the \code{geepack} package
#' @param study.data
#' @param model.form
#' @return fitted model
FitGEEModel <- function(study.data, 
                        model.form = formula(Yij ~ (Dwell + Music)^2)){
  
  gee.fit <- geepack::geeglm(formula = model.form, id=ID, data = study.data, 
                    corstr = "exchangeable", na.action="na.omit")
  return(gee.fit)
}

FitOrdGLMM <- function(observed.data, 
                    model.form = formula(Yij ~ (Dwell + Music + Viz + Squeeze)^2)){

  current.study.data$Yij <- factor(current.study.data$Yij, ordered = TRUE)

  
  gee.ord.fit <- geepack::ordgee(formula = model.form, data = study.data,
                           id = ID, corstr = "exchangeable")
  
  return(gee.ord.fit)
  
}

#' Calculate Wald statistics and associated p-values
CalcGEEPval <- function(gee.fit, contrast.mat){
  Sigma <- summary(gee.fit)$cov.scaled
  
  wald.stats <- vector(length = nrow(contrast.mat))
  wald.pvals <- vector(length = nrow(contrast.mat))
  beta_hat <- coef(gee.fit)
  
  for(i in 1:length(wald.stats)){
    current.contrast <- matrix(contrast.mat[i,], nrow=1)
    wald.outer <- (current.contrast %*% beta_hat)
    wald.inner <- solve(current.contrast %*% Sigma %*% t(current.contrast))
    wald.stats[i] <- t(wald.outer) %*% wald.inner %*% wald.outer
    wald.pvals[i] <- pchisq(q = wald.stats[i], df = 1, lower.tail = FALSE)
  }
  
  names(wald.pvals) <- rownames(contrast.mat)

  return(wald.pvals)
}
