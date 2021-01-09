


#' Solves for the total sample size given the proportion for each visit and a fixed
#' number of visits. Assumes that the visit proportions are ordered, so element one is the proportion
#' of patients with exactly one visit, element two is the proportion with exactly
#' two visits, and so on. 
#' 
#' @param total.visits the fixed sample size in terms of the number of visits/procedures
#' @param visit.props vector of proportions for each visit number. Must sum to one
#' 
#' @return an integer with the the sample size in terms of patients neccessary 
#' to achieve the sample size in visits given the proportion of patients with each
#' number of visits. 
SolveForNGivenVisitProps <- function(total.visits, visit.props){
  stopifnot(sum(visit.props) == 1)
  
  n <- ceiling(total.visits / sum(visit.props * 1:length(visit.props)))
  
  return(n)
}