################################################################################
## Functions for analyzing Simulation runs
##
################################################################################

#' 
#' @param pval.mat A matrix of p-values. The rows should be simulation runs, and 
#' the columns unadjusted p-values for hypotheses.
AdjustPvals <- function(pval.mat, adj.method = "holm"){
  adjusted_pval_mat <- t(apply(pval.mat, 1, FUN = p.adjust, method = adj.method))
  colnames(adjusted_pval_mat) <- colnames(pval.mat)
  return(as_tibble(adjusted_pval_mat))
}

PowerPlot <- function(){
  cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
            "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

  cutoffHelper <- function(indf, cutoff.val){
    to_return <- indf %>% group_by(N) %>% 
      summarise(Power = sum(PercOracle > cutoff.val)/n(),
                .groups = "drop_last") %>% 
      mutate(Cutoff = cutoff.val)
    
    return(to_return)
  }
  
  perc_oracle_power <- map_dfr(c(seq(0, .5, by = .1), seq(0.5, 1, by = .005)), ~cutoffHelper(indf = sim_summaries, cutoff.val = .))
  
  ggplot(data = perc_oracle_power, aes(x = Cutoff, y = Power, color = N)) + 
    geom_line() + 
    xlim(0, 1)+
    ylim(0,1)+
    scale_colour_manual(values = cbp1) + 
    labs(x = "Percentage of Oracle Value Cutoff",
         y = "Probability of Attaining Value Above Cutoff")
}