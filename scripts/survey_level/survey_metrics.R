###################################################################################################
#'  Compute community aesthe values at the survey level
#'
#' @author Nicolas Mouquet , \email{nicolas.mouquet@@cnrs.fr}
#' @author Matthew McLean, \email {mcleamj@@gmail.com},
#' @author Juliette Langlois, \email{juliette.a.langlois@@gmail.com}
#'         
#' 
#' Produce the file survey_aesth.csv which contains predicted aesthetic values for each survey 
#' Produce fig_1 b,c  
#'        
#' @date 2022/06/22
##################################################################################################

rm(list = ls())
library(ggplot2)
library(tidyverse)
library(lme4)

# INIT ----

  # coefficient a and b from Tribot, A.S, Deter, J., Claverie, T., Guillhaumon, F., Villeger, S., & Mouquet, N. (2019). Species diversity and composition drive the aesthetic value of coral reef fish assemblages. Biology letters, 15, 20190703, doi:10.1098/rsbl.2019.0703
  # will be used to compute the aesthetic scores of assemblages
    intercept_sr <- 7.0772149
    slope_sr     <- 0.20439752
    
  #load the files 
    
    aesthe_species <- read.csv2(here::here("outputs", "aesthe_species.csv"))
    sp_pres_matrix <- readRDS(here::here("outputs", "sp_pres_matrix.rds"))
    abund_matrix <- readRDS(here::here("outputs", "sp_abund_matrix.rds"))

# ----
    
# COMPUTE AESTHETICS ----
  #aesthe_survey_abund : predicted aesthe 
  #aesthe_SR_survey : predicted aesthe based only on species richness
      
      surveyid_vect <- as.character(sp_pres_matrix$SurveyID)  
      
      survey_aesth <- do.call(rbind, pbmcapply::pbmclapply((1:length(surveyid_vect)), function(i){   #

                                #i=1
                                survey_id <-  surveyid_vect[i]
                               # presence_absence of the species of the survey
                                
                                vector_abs_pres <- sp_pres_matrix[which(sp_pres_matrix$SurveyID == survey_id),]
                                vector_abs_pres <- vector_abs_pres[,-1]
                                
                                vector_abund <- abund_matrix[which(abund_matrix$SurveyID == survey_id),]
                                vector_abund <- vector_abund[,-1]
                                
                                vector_abund <- t(vector_abund) %>%
                                  as.data.frame() %>%
                                  rename(abundance = V1) %>%
                                  filter(abundance > 0) %>%
                                  rownames_to_column("sp_name") %>%
                                  mutate(sp_name = str_replace(sp_name, "_", " ")) %>%
                                  left_join(aesthe_species, by="sp_name") %>%
                                  mutate(abund_wtd_effect = log(abundance) * aesthe_effect)
                                
                               
                               # species present
                                sp_survey <- colnames(vector_abs_pres)[vector_abs_pres[1,]>0]

                               # number of species of the survey
                                nb_species <- length(sp_survey)
                                
                               # aesthe of the survey
                                E_presence <- intercept_sr + slope_sr * log(nb_species) + 
                                  sum(aesthe_species$aesthe_effect[aesthe_species$sp_name %in%  gsub("_"," ",sp_survey)])
                                score_presence <-  exp(E_presence)

                                E_abundance <- intercept_sr + slope_sr * log(nb_species) + 
                                  sum(vector_abund$abund_wtd_effect)
                                score_abundance <- exp(E_abundance)
                                
                                E <-  intercept_sr + slope_sr * log(nb_species)
                                score_SR  <-  exp(E) 
                                
                              # compute the number of species with positive and negative effect in each survey
                                
                                vect  <-  (aesthe_species$aesthe_effect * vector_abs_pres)
                                nb_sp_pos_survey   <-  length(which(vect>0))
                                nb_sp_neg_survey   <-  length(which(vect<0))
                                
                                cbind.data.frame(SurveyID=surveyid_vect[i],
                                                 nb_species=nb_species,
                                                 aesthe_survey_abund_pres=score_presence,
                                                 aesthe_survey_abund=score_abundance,
                                                 aesthe_SR_survey=score_SR,
                                                 nb_sp_pos_survey=nb_sp_pos_survey, 
                                                 nb_sp_neg_survey=nb_sp_neg_survey)

                             }, mc.cores = parallel::detectCores()-1))
      
      par(mfrow=c(1,3))
      plot(survey_aesth$nb_species,survey_aesth$aesthe_survey_abund_pres,
           xlab="Nb Species", ylab="Aesthetic Value (Presence)")
      points(survey_aesth$nb_species, survey_aesth$aesthe_SR_survey,col=2,pch=19)
      
      plot(survey_aesth$nb_species,survey_aesth$aesthe_survey_abund,
           xlab="Nb Species", ylab="Aesthetic Value (Abundance)")
      points(survey_aesth$nb_species, survey_aesth$aesthe_SR_survey,col=2,pch=19)
      
      plot(survey_aesth$aesthe_survey_abund_pres,survey_aesth$aesthe_survey_abund,
           xlab="Aesthetic Value (Presence)", ylab="Aesthetic Value (Abundance)")
      
      
      
      write.csv(survey_aesth, here::here("outputs", "survey_aesth.csv"), row.names = FALSE)
      
      
      # PLOT ABUNDANCE VS AESTHETIC SCORE
      
      abund_long <- abund_matrix %>%
        column_to_rownames("SurveyID") %>%
        pivot_longer(cols = everything(),
                       names_to = "sp_name",
                     values_to = "abundance") %>%
        mutate(sp_name = str_replace(sp_name, "_", " ")) %>%
        filter(abundance > 0) %>%
        left_join(aesthe_species, by="sp_name") %>%
        group_by(sp_name) %>%  # Group by species name
        mutate(aesthe_decile = ntile(aesthe_score, 10)) %>%  # Bin abundance into deciles
        ungroup()  # Ungroup after mutation
      
      hist(log1p(abund_long$abundance))
      
      ggplot(abund_long, aes(x=aesthe_decile, group=aesthe_decile, y=log(abundance))) +
        geom_boxplot()
      
      test <- abund_long %>%
        select(sp_name, aesthe_score, abundance) %>%
        group_by(sp_name) %>% 
        summarise_all(.funs=mean) %>%
        mutate(log_abund = log(abundance)) 
        
      ggplot(test, aes(x=aesthe_score,y=log_abund)) +
        geom_point() +
        geom_smooth(method="lm")
        
      
      ###################################################################
      ## RE COMPUTE AESTETHETICS WHILE INCORPORATING MODEL UNCERTAINTY ##
      ###################################################################
      
      # LOAD THE POSTERIOR DISTRIBUTIONS
      tribot_effects <- read_rds("data/tribot_effects.rds")
      
      n_samples <- nrow(tribot_effects)
      n_iterations <- length(surveyid_vect)
      
      # Generate predictions for each posterior sample
        
      results_list <- vector("list", n_iterations)
      
        for(i in 1:n_iterations){
          
          survey_id <-  surveyid_vect[i]
          # presence_absence of the species of the survey
          
          vector_abs_pres <- sp_pres_matrix[which(sp_pres_matrix$SurveyID == survey_id),]
          vector_abs_pres <- vector_abs_pres[,-1]
          
          vector_abund <- abund_matrix[which(abund_matrix$SurveyID == survey_id),]
          vector_abund <- vector_abund[,-1]
          
          vector_abund <- t(vector_abund) %>%
            as.data.frame() %>%
            rename(abundance = V1) %>%
            filter(abundance > 0) %>%
            rownames_to_column("sp_name") %>%
            mutate(sp_name = str_replace(sp_name, "_", " ")) %>%
            left_join(aesthe_species, by="sp_name") %>%
            mutate(abund_wtd_effect = log(abundance) * aesthe_effect)
          
          
          # species present
          sp_survey <- colnames(vector_abs_pres)[vector_abs_pres[1,]>0]
          
          # number of species of the survey
          nb_species <- length(sp_survey)
          
          score_holder <- matrix(NA, nrow = n_samples, ncol = 3)
          colnames(score_holder) <- c("score_presence", "score_abundance", "score_SR")
          
          for(j in 1:n_samples){
          
          # aesthe of the survey
          E_presence <- tribot_effects$b_Intercept[j] + tribot_effects$b_logsr[j] * log(nb_species) + 
            sum(aesthe_species$aesthe_effect[aesthe_species$sp_name %in%  gsub("_"," ",sp_survey)])
          score_presence <-  exp(E_presence)
          
          E_abundance <- tribot_effects$b_Intercept[j] + tribot_effects$b_logsr[j] * log(nb_species) + 
            sum(vector_abund$abund_wtd_effect)
          score_abundance <- exp(E_abundance)
          
          E <-  tribot_effects$b_Intercept[j] + tribot_effects$b_logsr[j] * log(nb_species)
          score_SR  <-  exp(E)
          
          score_holder[j,] <- c(score_presence, score_abundance, score_SR)
          
          }
          
          score_holder <- as.data.frame(score_holder)
          
          score_presence <- score_holder %>%
            select(score_presence) %>%
            summarise(median(score_presence)) %>%
            as.numeric()
          
          score_abundance <- score_holder %>%
            select(score_abundance) %>%
            summarise(median(score_abundance)) %>%
            as.numeric()
          
          score_SR <- score_holder %>%
            select(score_SR) %>%
            summarise(median(score_SR)) %>%
            as.numeric()
          
          # compute the number of species with positive and negative effect in each survey
          vect  <-  (aesthe_species$aesthe_effect * vector_abs_pres)
          nb_sp_pos_survey   <-  length(which(vect>0))
          nb_sp_neg_survey   <-  length(which(vect<0))
          
          results_list[[i]]  <- data.frame(SurveyID=surveyid_vect[i],
                           nb_species=nb_species,
                           aesthe_survey_abund_pres=score_presence,
                           aesthe_survey_abund=score_abundance,
                           aesthe_SR_survey=score_SR,
                           nb_sp_pos_survey=nb_sp_pos_survey, 
                           nb_sp_neg_survey=nb_sp_neg_survey)

        
        
        }
      
      results <- do.call(rbind, results_list)
      
      
      ## 0.999 correlation between original and version with uncertainty
      

# ---- 
    
#FIGURE 1b,c ----
#Fig1a was produced by other means     

    aesthe_species <- read.csv2(here::here("outputs", "aesthe_species.csv"))
    survey_aesth <- read.csv(here::here("outputs", "survey_aesth.csv"))
    sp_pres_matrix <- readRDS(here::here("outputs", "sp_pres_matrix.rds"))
    
    ##subset only species present in sp_pres_matrix.rds
    
    sp_survey <- gsub("_"," ",names(sp_pres_matrix)[-1])
    aesthe_species <- aesthe_species[aesthe_species$sp_name%in%sp_survey,]

    ##isolate two survey with same number of species but high and low aesthetic scores
    ##will be used in fig.1b and c
    
      temp  <- survey_aesth[which(survey_aesth$aesthe_survey_abund > 5000),]
      temp  <- temp[which(temp$aesthe_survey_abund < 5500),]
      temp  <- temp[which(temp$nb_species ==59),]
      n_sp <- temp$nb_species
      upsc <- temp[, c("SurveyID", "aesthe_survey_abund", "nb_species")]
      upsc <- upsc[1,]
      
      temp <- survey_aesth[which(survey_aesth$nb_species == upsc$nb_species),]
      lowsc <- temp[which(temp$aesthe_survey_abund == min(temp$aesthe_survey_abund)), c("SurveyID", "aesthe_survey_abund", "nb_species")]
    
      ###Fig.1b
    
        b <- ggplot(survey_aesth, ggplot2::aes(y=aesthe_survey_abund,x = nb_species)) +
          geom_point(col="#A0A0FA",alpha=0.5) +
          theme_bw()+
          theme(axis.text.x = element_text(size = 10),
                axis.title.x = element_text(size = 14),
                axis.text.y = element_text(size = 10),
                axis.title.y = element_text(size = 14))+
          geom_line(aes(y = aesthe_SR_survey, x = nb_species),col = "#7D7D7C",size=1.5)+
          geom_line(aes(y = aesthe_SR_survey, x = nb_species),col = "white",size=0.5)+
          labs(x = "Surveys species richness",
               y = "Survey aesthetic values")+
          geom_point(aes(x=lowsc$nb_species,y=lowsc$aesthe_survey_abund),colour="tomato",size=3)+
          geom_point(aes(x=upsc$nb_species,y=upsc$aesthe_survey_abund),colour="tomato",size=3)+
          geom_text(x=n_sp-9, y=lowsc$aesthe_survey_abund, label="Low",size=4,col="#7D7D7C")+
          geom_text(x=n_sp-9, y=upsc$aesthe_survey_abund, label="High",size=4,col="#7D7D7C")
        
      ###Fig.1c 
    
        sp_pres <- as.data.frame(sp_pres_matrix)
        rownames(sp_pres)=sp_pres$SurveyID
        sp_pres <- sp_pres[,-1]
        
        sp_lowsc <- colnames(sp_pres)[sp_pres[which(rownames(sp_pres)%in%lowsc$SurveyID),]>0]
        sp_upsc <- colnames(sp_pres)[sp_pres[which(rownames(sp_pres)%in%upsc$SurveyID),]>0]
        
        low_sc_effect  <-  cbind.data.frame(effect=aesthe_species$aesthe_effect[which(aesthe_species$sp_name %in% gsub("_"," ",sp_lowsc))])
        high_sc_effect <-  cbind.data.frame(effect=aesthe_species$aesthe_effect[which(aesthe_species$sp_name %in% gsub("_"," ",sp_upsc))])
        
        low_sc_effect$score_type  <- "Low"
        high_sc_effect$score_type <- "High"
        
        comp_surv <- rbind(low_sc_effect, high_sc_effect)
        
        c <- ggplot(comp_surv, aes(x = score_type, y = effect))+
          geom_boxplot(fill="white",width = .5, outlier.shape = NA) + 
          ylab("Species aesthetic effects") +
          xlab("Surveys")+
          theme_bw()+
          theme(axis.text.x = element_text(size = 10),
                axis.title.x = element_text(size = 14),
                axis.text.y = element_text(size = 10),
                axis.title.y = element_text(size = 14))+
          ylim(-0.04,0.04)+
          ggdist::stat_halfeye(adjust = .5,width = .6,.width = 0,
                               justification = -.3,alpha = .2,fill='blue',point_colour = NA)+
          geom_point(size = 1,alpha = .2,col="blue",
                     position = position_jitter(seed = 1, width = .1)) +
          geom_hline(yintercept=0, linetype="dashed", 
                     color = "gray", size=1)+
          scale_x_discrete(labels = c("Low", "High"))+
          ggplot2::theme(axis.line.x  = ggplot2::element_line(linetype = "blank"),
                         axis.ticks.x = ggplot2::element_blank())
        
        ###Assemble Fig 1 b and c and save 
  
          library(gridExtra)
          fig_1 <- gridExtra::arrangeGrob(b,c,ncol=2)
          fig_1
          ggsave(file=here::here("figures_tables","fig_1.tiff"), fig_1,width = 24, height = 10, dpi = 300, units = "cm", device='tiff') 
          ggsave(file=here::here("figures_tables","fig_1.pdf"), fig_1,width = 24, height = 10, dpi = 300, units = "cm", device='pdf') 
          
# ----
 
  
    
        
  