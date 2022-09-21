###################################################################################################
#'  Compute metrics at the survey level
#'
#' @author Matthew McLean, \email {mcleamj@@gmail.com},
#'         Juliette Langlois, \email{juliette.a.langlois@@gmail.com},
#'         Nicolas Mouquet , \email{nicolas.mouquet@@cnrs.fr},
#'         
#'        
#' @date 2022/06/22
##################################################################################################


rm(list = ls())


# INIT ----

  # coefficient a and b from Tribot, A.S, Deter, J., Claverie, T., Guillhaumon, F., Villeger, S., & Mouquet, N. (2019). Species diversity and composition drive the aesthetic value of coral reef fish assemblages. Biology letters, 15, 20190703, doi:10.1098/rsbl.2019.0703
    intercept_sr <- 7.0772149
    slope_sr     <- 0.20439752
    
  #test branch
    
  #load the files 
    
    all_species <- read.csv(here::here("outputs", "all_species.csv"))
    
    sp_pres_matrix <- readRDS(here::here("outputs", "sp_pres_matrix.rds"))
    
    
  #load functions 
    
    source(here::here("R", "functions_esth.R"))
    
    
# ----
    
    
# AESTHETICS ----
    
# Compute the effect of each species
    
    # Computing the aesthe_effect
      all_species$aesthe_effect <- (log(all_species$esthe_score) - 7.3468679)/7.937672
    
    # plot
      plot(all_species$aesthe_effect[order(all_species$esthe_score)], type="l")
      abline(v = mean(all_species$esthe_score))
      abline(h = 0)
      length(all_species$aesthe_effect)
    
    # checking how many species with negative/positive effect on the score
      sum(all_species$aesthe_effect < 0) # 1805 species have a negative effect on the score of a community
      sum(all_species$aesthe_effect > 0) # only 656 species have a positive effect on the score
      
      
# Compute the survey aesthe scores 
      
      load(here::here("results", "00_presabs.RData"))
      
      str(sp_pres_matrix)
      
      surveyid_vect <- as.character(sp_pres_matrix$SurveyID)  
      

      sp_pres_matrix$SurveyID
      
      ptm       <- proc.time()
      survey_aesth <- do.call(rbind, 
                             parallel::mclapply((1:length(surveyid_vect)), function(i){   #

                              # i=1
                              
                               # presence_absence of the species of the survey
                                
                                vector_abs_pres <- sp_pres_matrix[which(sp_pres_matrix$SurveyID == surveyid_vect[i]),]
                                
                                vector_abs_pres <- vector_abs_pres[,-1]
                               
                               # species present
                                sp_survey <- colnames(vector_abs_pres)[vector_abs_pres[1,]>0]

                               # number of species of the survey
                                nb_species <- length(sp_survey)
                                
                               # aesthe of the survey
                                E <- intercept_sr + slope_sr * log(nb_species) + sum(all_species$aesthe_effect[all_species$sp_name %in%  gsub("_"," ",sp_survey)])
                                score <-  exp(E)

                                E          <-  intercept_sr + slope_sr * log(nb_species)
                                score_SR      <-  exp(E) 
                                
                                cbind.data.frame(SurveyID=surveyid_vect[i],nb_species=nb_species,aesthe_survey=score,aesthe_SR_survey=score_SR)

                             }, mc.cores = 7))
      proc.time() - ptm
      
      survey_aesth <- as.data.frame(survey_esth)
      
      write.csv(survey_aesth, here::here("outputs", "survey_aesth.csv"), row.names = FALSE)
      
      plot(survey_aesth$nb_species,survey_esth$aesthe_survey)
    

# ---- 