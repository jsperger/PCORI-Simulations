#################################################
#### Proportion Assigned Best Treatment
#################################################

# For the no subgroup settings specifically
percent.best.treatment <- treatment.assignments[,10]/rowSums(treatment.assignments)
#################################################
#### Power
#################################################

gee.power <- sum(adjusted.gee.pvals < .05)/(nrow(adjusted.gee.pvals)*ncol(adjusted.gee.pvals))
mm.power <- sum(adjusted.mm.pvals < .05)/(nrow(adjusted.mm.pvals)*ncol(adjusted.mm.pvals))
# The lidocaine dwell has the lowest effect size, this is the power for that test vs Soc
gee.weakest.comp <- sum(adjusted.gee.pvals[,1] < .05)/nrow(adjusted.gee.pvals)
mm.weakest.comp <- sum(adjusted.mm.pvals[,1] < .05)/nrow(adjusted.mm.pvals)

