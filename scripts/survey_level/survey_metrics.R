###################################################################################################
#'  Compute aesthe at the survey level
#'
#' @author Nicolas Mouquet , \email{nicolas.mouquet@@cnrs.fr}
#' @author Matthew McLean, \email {mcleamj@@gmail.com},
#' @author Juliette Langlois, \email{juliette.a.langlois@@gmail.com}
#'         
#' 
#' Produce the file survey_aesth.csv which contains predicted aesthetic values for each survey 
#' Produce fig_1.png  
#'        
#' @date 2022/06/22
##################################################################################################

rm(list = ls())
library(ggplot2)

# INIT ----

  # coefficient a and b from Tribot, A.S, Deter, J., Claverie, T., Guillhaumon, F., Villeger, S., & Mouquet, N. (2019). Species diversity and composition drive the aesthetic value of coral reef fish assemblages. Biology letters, 15, 20190703, doi:10.1098/rsbl.2019.0703
  # will be used to compute the aesthetic scores of assemblages
    intercept_sr <- 7.0772149
    slope_sr     <- 0.20439752
    
  #load the files 
    
    aesthe_species <- read.csv2(here::here("outputs", "aesthe_species.csv"))
    sp_pres_matrix <- readRDS(here::here("outputs", "sp_pres_matrix.rds"))
    
# ----
    
# COMPUTE AESTHETICS ----
  #aesthe_survey : predicted aesthe 
  #aesthe_SR_survey : predicted aesthe based only on species richness
      
      surveyid_vect <- as.character(sp_pres_matrix$SurveyID)  
      
      survey_aesth <- do.call(rbind, pbmcapply::pbmclapply((1:length(surveyid_vect)), function(i){   #

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
                                E <- intercept_sr + slope_sr * log(nb_species) + sum(aesthe_species$aesthe_effect[aesthe_species$sp_name %in%  gsub("_"," ",sp_survey)])
                                score <-  exp(E)

                                E <-  intercept_sr + slope_sr * log(nb_species)
                                score_SR  <-  exp(E) 
                                
                              # compute the number of species with positive and negative effect in each survey
                                
                                vect  <-  (aesthe_species$aesthe_effect * vector_abs_pres)
                                nb_sp_pos_survey   <-  length(which(vect>0))
                                nb_sp_neg_survey   <-  length(which(vect<0))
                                
                                cbind.data.frame(SurveyID=surveyid_vect[i],nb_species=nb_species,aesthe_survey=score,aesthe_SR_survey=score_SR,
                                                 nb_sp_pos_survey=nb_sp_pos_survey, nb_sp_neg_survey=nb_sp_neg_survey)

                             }, mc.cores = parallel::detectCores()-1))
      
      plot(survey_aesth$nb_species,survey_aesth$aesthe_survey)
      
      write.csv(survey_aesth, here::here("outputs", "survey_aesth.csv"), row.names = FALSE)

# ---- 
      
#FIGURE 1 ----
 
  aesthe_species <- read.csv(here::here("outputs", "aesthe_species.csv"))
  survey_aesth <- read.csv(here::here("outputs", "survey_aesth.csv"))
  sp_pres_matrix <- readRDS(here::here("outputs", "sp_pres_matrix.rds"))
  
  #FIG 1 a 
  
  aesthe_species$rank <- rank(aesthe_species$aesthe_effect)
  
    nneg <- sum(aesthe_species$aesthe_effect < 0) # 1763 species have a negative effect
    npos <- sum(aesthe_species$aesthe_effect > 0) # 654 species have a positive effect
    
    ex <- c("Chaetodontoplus septentrionalis","Calloplesiops altivelis","Paracanthurus hepatus",
            "Chelmonops curiosus","Acanthurus blochii","Scarus rivulatus","Pseudolabrus guentheri",
            "Scolopsis ghanam","Sarpa salpa","Salmo salar","Arripis truttacea","Sardina pilchardus")
    
    colsp <- "#40403F"
    size_cr <- 6
    font_cr <- 3
    
    ex_fl <- c("Chaetodontoplus_septentrionalis_A_1.png","Calloplesiops_altivelis_A_2.png","Paracanthurus_hepatus_A_3.png",
               "Chelmonops_curiosus_A_2.png","Acanthurus_blochii_A_1.png","Scarus_rivulatus_M_2.png",
               "Pseudolabrus_guentheri_A_2.png","Scolopsis_ghanam_A_3.png","Sarpa_salpa_A_2.png",
               "Salmo_salar_A_2.png","Arripis_truttacea_A_1.png","Sardina_pilchardus_A_1.png")
  
    a <- ggplot(aesthe_species, ggplot2::aes(y=aesthe_effect,x = rank)) +
      geom_point(col="royalblue1",alpha=0.5)+
      theme_bw()+
      xlab("Ranks") + ylab("Net species aesthetic effects")+
      theme(axis.text.x = element_text(size = 10),
            axis.title.x = element_text(size = 14),
            axis.text.y = element_text(size = 10),
            axis.title.y = element_text(size = 14))+
      ylim(-0.04,0.04)+
      geom_hline(yintercept=0, linetype="dashed", 
                 color = "gray", size=0.8)+
      geom_vline(xintercept=nneg, linetype="dashed", 
                 color = "gray", size=0.8)+
      geom_text(x=100, y=-0.003, label=paste("n=",nneg),size=4,col="#7D7D7C")+
      geom_text(x=2400, y=+0.003, label=paste("n=",npos),size=4,col="#7D7D7C")+
      geom_point(aes(x=rank[sp_name%in%ex[1]],
                     y=aesthe_effect[sp_name%in%ex[1]]),
                 shape = 21, colour = colsp, fill = "white",size=size_cr)+
      geom_text(x=aesthe_species$rank[aesthe_species$sp_name%in%ex[1]], 
                y=aesthe_species$aesthe_effect[aesthe_species$sp_name%in%ex[1]], 
                label="1",size=font_cr,col=colsp)+
      geom_point(aes(x=rank[sp_name%in%ex[2]],
                     y=aesthe_effect[sp_name%in%ex[2]]),
                 shape = 21, colour = colsp, fill = "white",size=size_cr)+
      geom_text(x=aesthe_species$rank[aesthe_species$sp_name%in%ex[2]], 
                y=aesthe_species$aesthe_effect[aesthe_species$sp_name%in%ex[2]], 
                label="2",size=font_cr,col=colsp)+
      geom_point(aes(x=rank[sp_name%in%ex[3]],
                     y=aesthe_effect[sp_name%in%ex[3]]),
                 shape = 21, colour = colsp, fill = "white",size=size_cr)+
      geom_text(x=aesthe_species$rank[aesthe_species$sp_name%in%ex[3]], 
                y=aesthe_species$aesthe_effect[aesthe_species$sp_name%in%ex[3]], 
                label="3",size=font_cr,col=colsp)+
      geom_point(aes(x=rank[sp_name%in%ex[4]],
                     y=aesthe_effect[sp_name%in%ex[4]]),
                 shape = 21, colour = colsp, fill = "white",size=size_cr)+
      geom_text(x=aesthe_species$rank[aesthe_species$sp_name%in%ex[4]], 
                y=aesthe_species$aesthe_effect[aesthe_species$sp_name%in%ex[4]], 
                label="4",size=font_cr,col=colsp)+
      geom_point(aes(x=rank[sp_name%in%ex[5]],
                     y=aesthe_effect[sp_name%in%ex[5]]),
                 shape = 21, colour = colsp, fill = "white",size=size_cr)+
      geom_text(x=aesthe_species$rank[aesthe_species$sp_name%in%ex[5]], 
                y=aesthe_species$aesthe_effect[aesthe_species$sp_name%in%ex[5]], 
                label="5",size=font_cr,col=colsp)+
      geom_point(aes(x=rank[sp_name%in%ex[6]],
                     y=aesthe_effect[sp_name%in%ex[6]]),
                 shape = 21, colour = colsp, fill = "white",size=size_cr)+
      geom_text(x=aesthe_species$rank[aesthe_species$sp_name%in%ex[6]], 
                y=aesthe_species$aesthe_effect[aesthe_species$sp_name%in%ex[6]], 
                label="6",size=font_cr,col=colsp)+
      geom_point(aes(x=rank[sp_name%in%ex[7]],
                     y=aesthe_effect[sp_name%in%ex[7]]),
                 shape = 21, colour = colsp, fill = "white",size=size_cr)+
      geom_text(x=aesthe_species$rank[aesthe_species$sp_name%in%ex[7]], 
                y=aesthe_species$aesthe_effect[aesthe_species$sp_name%in%ex[7]], 
                label="7",size=font_cr,col=colsp)+
      geom_point(aes(x=rank[sp_name%in%ex[8]],
                     y=aesthe_effect[sp_name%in%ex[8]]),
                 shape = 21, colour = colsp, fill = "white",size=size_cr)+
      geom_text(x=aesthe_species$rank[aesthe_species$sp_name%in%ex[8]], 
                y=aesthe_species$aesthe_effect[aesthe_species$sp_name%in%ex[8]], 
                label="8",size=font_cr,col=colsp)+
      geom_point(aes(x=rank[sp_name%in%ex[9]],
                     y=aesthe_effect[sp_name%in%ex[9]]),
                 shape = 21, colour = colsp, fill = "white",size=size_cr)+
      geom_text(x=aesthe_species$rank[aesthe_species$sp_name%in%ex[9]], 
                y=aesthe_species$aesthe_effect[aesthe_species$sp_name%in%ex[9]], 
                label="9",size=font_cr,col=colsp)+
      geom_point(aes(x=rank[sp_name%in%ex[10]],
                     y=aesthe_effect[sp_name%in%ex[10]]),
                 shape = 21, colour = colsp, fill = "white",size=size_cr)+
      geom_text(x=aesthe_species$rank[aesthe_species$sp_name%in%ex[10]], 
                y=aesthe_species$aesthe_effect[aesthe_species$sp_name%in%ex[10]], 
                label="10",size=font_cr,col=colsp)+
      geom_point(aes(x=rank[sp_name%in%ex[11]],
                     y=aesthe_effect[sp_name%in%ex[11]]),
                 shape = 21, colour = colsp, fill = "white",size=size_cr)+
      geom_text(x=aesthe_species$rank[aesthe_species$sp_name%in%ex[11]], 
                y=aesthe_species$aesthe_effect[aesthe_species$sp_name%in%ex[11]], 
                label="11",size=font_cr,col=colsp)+
      geom_point(aes(x=rank[sp_name%in%ex[12]],
                     y=aesthe_effect[sp_name%in%ex[12]]),
                 shape = 21, colour = colsp, fill = "white",size=size_cr)+
      geom_text(x=aesthe_species$rank[aesthe_species$sp_name%in%ex[12]], 
                y=aesthe_species$aesthe_effect[aesthe_species$sp_name%in%ex[12]], 
                label="12",size=font_cr,col=colsp)
    
  
  #FIG 1 b,c 
  
  ##add the ranks and the diff in aesthetic (between aesthe_survey and aesthe_SR_survey)
    survey_aesth$rank <- rank(survey_aesth$aesthe_survey)
    #survey_aesth$ecart_richness <- survey_aesth$aesthe_survey - survey_aesth$aesthe_SR_survey
    
  ##isolate two survey with same number of species but high and low aesthetic scores
  ##will be used in fig.1b and c
    
  ##Fig.1b

    temp  <- survey_aesth[which(survey_aesth$aesthe_survey > 3500),]
    temp  <- temp[which(temp$aesthe_survey < 3600),]
    temp  <- temp[which(temp$nb_species ==40),]
    upsc <- temp[, c("SurveyID", "aesthe_survey", "nb_species")]
    upsc <- upsc[1,]
    
    temp <- survey_aesth[which(survey_aesth$nb_species == upsc$nb_species),]
    lowsc <- temp[which(temp$aesthe_survey == min(temp$aesthe_survey)), c("SurveyID", "aesthe_survey", "nb_species")]
      
      b <- ggplot(survey_aesth, ggplot2::aes(y=aesthe_survey,x = nb_species)) +
        geom_point(col="royalblue1",alpha=0.5) +
        theme_bw()+
        theme(axis.text.x = element_text(size = 10),
              axis.title.x = element_text(size = 14),
              axis.text.y = element_text(size = 10),
              axis.title.y = element_text(size = 14))+
        geom_line(aes(y = aesthe_SR_survey, x = nb_species),col = "#7D7D7C",size=1.5)+
        geom_line(aes(y = aesthe_SR_survey, x = nb_species),col = "white",size=0.5)+
        labs(x = "Number of species in the survey",
             y = "Predicted aesthetic score")+
        geom_point(aes(x=lowsc$nb_species,y=lowsc$aesthe_survey),colour="tomato",size=3)+
        geom_point(aes(x=upsc$nb_species,y=upsc$aesthe_survey),colour="tomato",size=3)+
        geom_text(x=48.5, y=lowsc$aesthe_survey, label="low",size=4,col="#7D7D7C")+
        geom_text(x=31.5, y=upsc$aesthe_survey, label="high",size=4,col="#7D7D7C")
        

    ##Fig.1c 
          
      sp_pres <- as.data.frame(sp_pres_matrix)
      rownames(sp_pres)=sp_pres$SurveyID
      sp_pres <- sp_pres[,-1]
  
      sp_lowsc <- colnames(sp_pres)[sp_pres[which(rownames(sp_pres)%in%lowsc$SurveyID),]>0]
      sp_upsc <- colnames(sp_pres)[sp_pres[which(rownames(sp_pres)%in%upsc$SurveyID),]>0]
      
      low_sc_effect  <-  cbind.data.frame(effect=aesthe_species$aesthe_effect[which(aesthe_species$sp_name %in% gsub("_"," ",sp_lowsc))])
      high_sc_effect <-  cbind.data.frame(effect=aesthe_species$aesthe_effect[which(aesthe_species$sp_name %in% gsub("_"," ",sp_upsc))])
          
      low_sc_effect$score_type  <- "low"
      high_sc_effect$score_type <- "high"
          
      comp_surv <- rbind(low_sc_effect, high_sc_effect)
      
      c <- ggplot(comp_surv, aes(x = score_type, y = effect))+
        geom_boxplot(fill="white",width = .5, outlier.shape = NA) + 
        ylab("Species aesthetic effects") +
        xlab("Survey")+
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
        scale_x_discrete(labels = c("High", "Low"))+
        ggplot2::theme(axis.line.x  = ggplot2::element_line(linetype = "blank"),
                       axis.ticks.x = ggplot2::element_blank())

    ##Assemble Fig 1 and save 
      
    library("cowplot")
      fig1 <- ggdraw() +
        draw_plot(a, x = 0, y = .4, width = 1, height = .6) +
        draw_plot(b, x = 0, y = 0, width = .5, height = .4) +
        draw_plot(c, x = 0.5, y = 0, width = 0.5, height = 0.4)
      
    ggsave(here::here("figures_tables","fig_1.png"), plot = fig1,
             width = 11, height = 8, dpi = 300, units = "in", device='png')
      
# ----     
 
  
    
        
  