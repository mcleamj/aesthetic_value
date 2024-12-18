###################################################################################################
#'  Set up the main datafile with aesthe scores and classification from taxsize
#'
#' @author Nicolas Mouquet , \email{nicolas.mouquet@@cnrs.fr}
#' @author Matthew McLean, \email {mcleamj@@gmail.com},
#' @author Juliette Langlois, \email{juliette.a.langlois@@gmail.com}
#'
#'  Produce the file aesthe_species.csv with all species aesthe scores and classification taxo
#'         
#' @date 2022/04/20
##################################################################################################


#Get species list with esthe scores from Langlois et al. 2022 and update species names 

  aesthe_species <- read.csv2(here::here("data", "aesthe_langlois_2022.csv"))
  name_changes <- read.csv(here::here("data", "name_changes.csv"))

  aesthe_species$sp_name[aesthe_species$sp_name%in%name_changes$old]=name_changes$new
  
  aesthe_species$sp_name <- as.character(gsub("_"," ",aesthe_species$sp_name))
  aesthe_species$aesthe_score <- as.numeric(aesthe_species$aesthe_score)
  
  ## ARE THERE ANY DUPLICATES NOW?
  aesthe_species$sp_name[which(duplicated(aesthe_species$sp_name))]
  
  # THIS MEANS THE ORIGINAL DATA CONTAINED BOTH THE OLD AND NEW NAMES
  # FOR THESE 2 SPECIES
  # WE KEEP THE VERSION WITH THE HIGHER AESTHETIC SCORE 
  # AND DELETE THE DUPLICATE
  
  aesthe_species <- aesthe_species %>%
    group_by(sp_name) %>%
    slice_max(aesthe_score, with_ties = FALSE) %>%
    ungroup()

# Compute the aesthe contribution of each species
  # with parameters from Tribot, A.S, Deter, J., Claverie, T., Guillhaumon, F., Villeger, S., & Mouquet, N. (2019). Species diversity and composition drive the aesthetic value of coral reef fish assemblages. Biology letters, 15, 20190703, doi:10.1098/rsbl.2019.0703
  # positive and negative effect are relative to the espected effect computed with the species 
    
  ## Computing the aesthe_effect
  aesthe_species$aesthe_effect <- (log(aesthe_species$aesthe_score) - 7.3468679)/7.937672
  
#merge & save 
  
  write.csv2(x = aesthe_species, file = here::here("outputs", "aesthe_species.csv"),
            row.names = FALSE)







