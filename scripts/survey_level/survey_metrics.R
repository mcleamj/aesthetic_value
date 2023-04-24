###################################################################################################
#'  Compute aesthe at the survey level
#'
#' @author Matthew McLean, \email {mcleamj@@gmail.com},
#'         Juliette Langlois, \email{juliette.a.langlois@@gmail.com},
#'         Nicolas Mouquet , \email{nicolas.mouquet@@cnrs.fr},
#' 
#' Produce the file survey_aesth.csv which contains predicted aesthetic values for each survey 
#' Produce FIG_2.png  
#'        
#' @date 2022/06/22
##################################################################################################

rm(list = ls())

# INIT ----

  # coefficient a and b from Tribot, A.S, Deter, J., Claverie, T., Guillhaumon, F., Villeger, S., & Mouquet, N. (2019). Species diversity and composition drive the aesthetic value of coral reef fish assemblages. Biology letters, 15, 20190703, doi:10.1098/rsbl.2019.0703
  # will be used to compute the aesthetic scores of assemblages
    intercept_sr <- 7.0772149
    slope_sr     <- 0.20439752
    
  #load the files 
    
    all_species <- read.csv(here::here("outputs", "all_species.csv"))
    sp_pres_matrix <- readRDS(here::here("outputs", "sp_pres_matrix.rds"))
    
# ----
    
# COMPUTE AESTHETICS ----
  #aesthe_survey : predicted aesthe 
  #aesthe_SR_survey : predicted aesthe based only on species richness
      
      surveyid_vect <- as.character(sp_pres_matrix$SurveyID)  
      
      survey_aesth <- do.call(rbind, 
                              pbmcapply::pbmclapply((1:length(surveyid_vect)), function(i){   #

                                #i=1
                                survey_id <-  surveyid_vect[i]
                               # presence_absence of the species of the survey
                                
                                vector_abs_pres <- sp_pres_matrix[which(sp_pres_matrix$SurveyID == survey_id),]
                                
                                vector_abs_pres <- vector_abs_pres[,-1]
                               
                               # species present
                                sp_survey <- colnames(vector_abs_pres)[vector_abs_pres[1,]>0]

                               # number of species of the survey
                                nb_species <- length(sp_survey)
                                
                               # aesthe of the survey
                                E <- intercept_sr + slope_sr * log(nb_species) + sum(all_species$aesthe_effect[all_species$sp_name %in%  gsub("_"," ",sp_survey)])
                                score <-  exp(E)

                                E <-  intercept_sr + slope_sr * log(nb_species)
                                score_SR  <-  exp(E) 
                                
                              # compute the number of species with positive and negative effect in each survey
                                
                                vect  <-  (all_species$aesthe_effect * vector_abs_pres)
                                nb_sp_pos_survey   <-  length(which(vect>0))
                                nb_sp_neg_survey   <-  length(which(vect<0))
                                
                                cbind.data.frame(SurveyID=surveyid_vect[i],nb_species=nb_species,aesthe_survey=score,aesthe_SR_survey=score_SR,
                                                 nb_sp_pos_survey=nb_sp_pos_survey, nb_sp_neg_survey=nb_sp_neg_survey)

                             }, mc.cores = 10))
      plot(survey_aesth$nb_species,survey_aesth$aesthe_survey)
      
      write.csv(survey_aesth, here::here("outputs", "survey_aesth.csv"), row.names = FALSE)

# ---- 
      
# PLOT THE FIGURE 1 ----
 
  all_species <- read.csv(here::here("outputs", "all_species.csv"))
  survey_aesth <- read.csv(here::here("outputs", "survey_aesth.csv"))
      
  ##add the ranks and the diff in aesthetic (between aesthe_survey and aesthe_SR_survey)
    survey_aesth$rank <- rank(survey_aesth$aesthe_survey)
    #survey_aesth$ecart_richness <- survey_aesth$aesthe_survey - survey_aesth$aesthe_SR_survey
    
  ##isolate two survey with same number of species but high and low aesthetic scores
  ##will be used in fig.2a and b

    temp  <- survey_aesth[which(survey_aesth$aesthe_survey > 3500),]
    temp  <- temp[which(temp$aesthe_survey < 3600),]
    temp  <- temp[which(temp$nb_species ==40),]
    upsc <- temp[, c("SurveyID", "aesthe_survey", "nb_species")]
    upsc <- upsc[1,]
    
    temp <- survey_aesth[which(survey_aesth$nb_species == upsc$nb_species),]
    lowsc <- temp[which(temp$aesthe_survey == min(temp$aesthe_survey)), c("SurveyID", "aesthe_survey", "nb_species")]
   
  ##FIG_2a
      
      a <- ggplot(survey_aesth, ggplot2::aes(y=aesthe_survey,x = nb_species)) +
        geom_point(col="royalblue1",alpha=0.5) +
        theme_bw()+
        geom_line(aes(y = aesthe_SR_survey, x = nb_species),col = "tomato")+
        labs(x = "Number of species in the survey",
             y = "Predicted aesthetic score")+
        geom_point(aes(x=lowsc$nb_species,y=lowsc$aesthe_survey),colour="tomato",size=3)+
        geom_point(aes(x=upsc$nb_species,y=upsc$aesthe_survey),colour="tomato",size=3)

  ##FIG_2b 
          
      sp_pres <- as.data.frame(sp_pres_matrix)
      rownames(sp_pres)=sp_pres$SurveyID
      sp_pres <- sp_pres[,-1]
  
      sp_lowsc <- colnames(sp_pres)[sp_pres[which(rownames(sp_pres)%in%lowsc$SurveyID),]>0]
      sp_upsc <- colnames(sp_pres)[sp_pres[which(rownames(sp_pres)%in%upsc$SurveyID),]>0]
      
      low_sc_effect  <-  cbind.data.frame(effect=all_species$aesthe_effect[which(all_species$sp_name %in% gsub("_"," ",sp_lowsc))])
      high_sc_effect <-  cbind.data.frame(effect=all_species$aesthe_effect[which(all_species$sp_name %in% gsub("_"," ",sp_upsc))])
          
      low_sc_effect$score_type  <- "low"
      high_sc_effect$score_type <- "high"
          
      comp_surv <- rbind(low_sc_effect, high_sc_effect)
      
      b <- ggplot(comp_surv, aes(x = score_type, y = effect))+
        geom_boxplot(fill="white",width = .5, outlier.shape = NA) + 
        ylab("Species aesthetic effects") +
        xlab("Survey")+
        theme_bw()+
        ggdist::stat_halfeye(adjust = .5,width = .6,.width = 0,
          justification = -.3,alpha = .2,fill='blue',point_colour = NA)+
        geom_point(size = 1,alpha = .2,col="blue",
          position = position_jitter(seed = 1, width = .1)) +
        geom_hline(yintercept=0, linetype="dashed", 
                   color = "gray", size=1)+
        scale_x_discrete(labels = c("High", "Low"))+
        ggplot2::theme(axis.line.x  = ggplot2::element_line(linetype = "blank"),
                       axis.ticks.x = ggplot2::element_blank())

    ##Assemble and save 
      
    fig_2 <- gridExtra::arrangeGrob(a,b,ncol=2)
    ggsave(here::here("figures_tables","fig_2.png"), plot = fig_2,
             width = 10, height = 4.5, dpi = 300, units = "in", device='png')
      
# ----     

      
      
      
      
      
      # -----
      
      # Violin plot of two stations with same nb sp ----
      
    
      ggplot2::ggsave(plot = survey_violin,
                      filename = here::here("results", "01_violin_esthval_sprich.pdf"), 
                      height = 18, width = 18, units = "cm", dpi = 320)
      
      # ----
      
      # Distribution of the aesthetic value of the surveys ----
      survey_esth <- read.csv(here::here("results", "01_survey_table_esth.csv"))
      densplot    <- 
        ggplot2::ggplot(survey_esth, ggplot2::aes(x = esth)) +
        ggplot2::geom_density(alpha = 0.7,
                              fill = viridis::viridis(1, alpha = 0.5, begin = 0.4, end = 0.6),
                              color = viridis::viridis(1, alpha = 0.5, begin = 0.4, end = 0.6)) + 
        ggplot2::theme_light() + 
        ggplot2::xlab("Aesthetic values") + 
        ggplot2::ylab("Density")  +
        ggplot2::scale_x_continuous(breaks = c(1000, 2000, 3000, 4000)) + 
        ggplot2::theme(axis.title      = ggplot2::element_text(size = 14, family = "serif"),
                       axis.text       = ggplot2::element_text(size = 9, family = "serif"), 
                       panel.grid      = ggplot2::element_blank(),
                       legend.position = "none")
      # save
      ggplot2::ggsave(plot = densplot,
                      filename = here::here("results", "01_survey_density_esth.pdf"), 
                      height = 18, width = 18, units = "cm", dpi = 320)
      
      # ----
      
      rm(comp_surv, data_commu, densplot, effect_sp, high_sc_effect, low_sc_effect, lowsc, plot,
         species_effect, species_table, survey_compo, survey_esth, survey_table, 
         survey_violin, upsc, intercept_sr, slope_sr, sp_lowsc, sp_upsc)
      
      
      