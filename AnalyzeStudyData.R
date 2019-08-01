#################################################
#### Data Analysis
#################################################
require(nlme)
require(geepack)
#library(DynTxRegime)
#library(multcomp)

###########################################################
### Mixed Model
###########################################################

FitMixedModel <- function(observed.data, model.form = formula(Obsij ~ (Dwell + Music + Viz + Squeeze)^2)){
  mixed.fit <- gls(model = model.form, data = observed.data, 
                   correlation = corCompSymm(form = ~ 1 | ID),
                   na.action = na.omit)
  #Exclude intercept p-vals
  mixed.pvals <- summary(mixed.fit)$tTable[-1,4]
  return(mixed.pvals)
}

###########################################################
### GEE
###########################################################
FitGEEModel <- function(observed.data, model.form = formula(Obsij ~ (Dwell + Music + Viz + Squeeze)^2)){
  gee.fit <- geeglm(formula = model.form, id=ID, data = observed.data, 
                    corstr = "exchangeable", na.action="na.omit")
  #Exclude intercept pvals
  gee.pvals <- summary(gee.fit)$coefficients[-1,4]
  return(gee.pvals)
}

###########################################################
### Extract p-values of interest
###########################################################

#K <- diag(length(coef(m1)))[-1,]
#rownames(K) <- names(coef(m1))[-1]
#test1 <-glht(m1, linfct=K)




#gee.adjusted.pvals <- gee.pvals*10
