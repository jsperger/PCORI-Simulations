#################################################
#### Data Analysis
#################################################
require(nlme)
require(geepack)

###########################################################
### Fit Model
###########################################################

m1 <- gls(Obsij ~ (Dwell + Music + Viz + Squeeze)^2, data = observed.data, 
          correlation = corCompSymm(form = ~ 1 | ID),
          na.action = na.omit)

m3 <- geeglm(Obsij ~ (Dwell + Music + Viz + Squeeze)^2, id=ID, data = observed.data, 
          corstr = "exchangeable", na.action="na.omit")
###########################################################
### Calculate p-values of interest
###########################################################