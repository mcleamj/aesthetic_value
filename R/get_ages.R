#' Get Ages
#' 
#' @description 
#' Computes the age of all the species on the leaves of a phylgenetic tree.
#'
#' @param tree A phylogenetic object
#'
#' @return A vector of ages of the leaves of tree
#' 
#' @export

get_ages <- function(tree) {
  
  nsp <- length(tree$tip.label) # number of species
  
  tips <- which(tree$edge[ , 2] <= nsp) # get the starting and ending nodes of each edge
  # security to take only the ones inferior to the number of species ie "leaves" 
  
  ages <- tree$edge.length[tips] # the age of a species is actually the length of the edges from
  # the first node of the phylo tree to the species' leaf
  
  names(ages) <- tree$tip.label
  
  ages
}