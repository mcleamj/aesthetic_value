###################################################################################################
#' Functions used in the script analysis/01_esth_score.R
#' 
#' @author Matthew McLean, \email {mcleamj@@gmail.com},
#'         Juliette Langlois, \email{juliette.a.langlois@@gmail.com},
#'         Nicolas Mouquet, \email{nicolas.mouquet@@cnrs.fr},
#'         François Guilhaumon, \email{francois.guilhaumon@@ird.fr},
#'         Alienor Stahl, \email{a.stahl67@@gmail.com}
#'
#' @date 2021/08/11
##################################################################################################

#' total_value
#' Compute the total aesthetic value of an assemblage of species
#'
#' @param commu_id id of the assemblage
#' @param matrix_diversity dataframe of minimum two columns:
#'  the assemblage id
#'  the number of species in each assemblage.
#' @param matrix_present_missing dataframe of presence absence. 
#'   Species as columns, 
#'   Communities as rows
#' @param matrix_effect_species dataframe of at least two columns:
#'   the species name
#'   the effect of each species on the aesthetic value of an assemblage
#'
#' @return a number, the total aesthetic value of the assemblage
#' @export
#'
total_value <- 
  function(commu_id, matrix_diversity, matrix_present_missing, matrix_effect_species){
    
  # commu_id = survey_id; matrix_diversity = dataframe_survey 
  # matrix_present_missing = data_presab; matrix_effect_species  = effect_sp
    
  # we want the presence/absence of the species of the station
  vector_abs_pres <- matrix_present_missing[which(rownames(matrix_present_missing) == commu_id),]
  
  # number of species of the station
  nb_species <- matrix_diversity$nb_species[matrix_diversity$survey_id == commu_id]
  
  # calculating the log(score)
  E     <-  intercept_sr + slope_sr * log(nb_species) + 
    sum(matrix_effect_species$effect * vector_abs_pres)
  score <-  exp(E) 
  
  return(score)
  
} # eo function total_value


#' esth_sprichness
#' Compute the aesthetic value of an assemblage of species considering only the species richness
#'
#' @param commu_id id of the assemblage
#' @param matrix_diversity dataframe of minimum two columns:
#'  the assemblage id
#'  the number of species in each assemblage.
#'
#' @return a number, the aesthetic value of the assemblage considering only the species richness
#' @export
#'
esth_sprichness <- function(commu_id, matrix_diversity) {
  
  nb_species <- matrix_diversity$nb_species[matrix_diversity$survey_id == commu_id]
  E          <-  intercept_sr + slope_sr * log(nb_species)
  score      <-  exp(E) 
  
  return(score)
  
} # eo esth_sprichness

#' esth_sprichness
#' For an assemblage compute the number of species with respectively a positive or negative effect
#'  on the aesthetic value of the assemblage.
#'  
#' @param commu_id id of the assemblage
#' @param matrix_present_missing dataframe of presence absence. 
#'   Species as columns, 
#'   Communities as rows
#' @param matrix_effect_species dataframe of at least two columns:
#'   the species name
#'   the effect of each species on the aesthetic value of an assemblage
#'
#' @return two numbers: the number of species with a positive effect and the number of species with a negative effect
#' @export
#'
nb_species_sign <- function(commu_id, matrix_present_missing, matrix_effect_species) {
  
  # presence/absence of the species for the assemblage
  vector_abs_pres <- matrix_present_missing[commu_id,]
  
  vect            <-  (matrix_effect_species$effect * vector_abs_pres)
  pos             <-  length(which(vect>0))
  neg             <-  length(which(vect<0))
  
  return(c(pos, neg))
  
} # eo function nb_species_sign

