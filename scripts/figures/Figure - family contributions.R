
#######################################################################################
#'  CODE TO PLOT THE FAMILY MPA CONTRIBUTION FIGURE
#'
#' @author Matthew McLean, \email {mcleamj@gmail.com},
#'         
#' @date AUGUST 9, 2023
########################################################################################

######################
## LIBRARY PACKAGES ##
######################

if(!require(plot3D)){install.packages("plot3D"); library(plot3D)}
if(!require(pals)){install.packages("pals"); library(pals)}
if(!require(brms)){install.packages("brms"); library(brms)}
if(!require(bayesplot)){install.packages("bayesplot"); library(bayesplot)}
if(!require(parallel)){install.packages("parallel"); library(parallel)}
if(!require(rstan)){install.packages("rstan"); library(rstan)}
if(!require(bayestestR)){install.packages("bayestestR"); library(bayestestR)}
if(!require(readr)){install.packages("readr"); library(readr)}
if(!require(tibble)){install.packages("tibble"); library(tibble)}
if(!require(dplyr)){install.packages("dplyr"); library(dplyr)}

##############################
## IMPORT DATA FROM OUTPUTS ##
##############################

fit_list <- read_rds("outputs/abundance_fit_list.rds")
family_info <- read_rds("outputs/family_info.rds")

##############################################
## WHICH FAMILIES HAVE HIGHEST MPA EFFECTS? ##
##############################################

family_MPA <- data.frame(family=colnames(fit_list),
                         MPA_effect_mean = colMeans(fit_list)) %>%
  arrange(desc(MPA_effect_mean))

##################################################
## CALCULATE AVERAGE AESTHETIC VALUE PER FAMILY ##
##################################################

aesthe_species <- read.csv2(here::here("data", "aesthe_langlois_2022.csv"))
aesthe_species$sp_name <- as.character(gsub("_"," ",aesthe_species$sp_name))
aesthe_species$aesthe_score <- as.numeric(aesthe_species$aesthe_score)

aesthe_species <- aesthe_species %>% rename(scientificname = sp_name)

aesthe_species <- merge(aesthe_species, family_info[,c("scientificname", "family")],
                        by="scientificname")
family_beauty <- aesthe_species %>%
  select(-scientificname) %>%
  group_by(family) %>%
  summarise_all(.funs=mean)


# SUMMARIZE MODEL OUPTUTS #

MPA_effect_summary <- as.data.frame(brms::posterior_summary(fit_list,
                                                            probs=c(0.10,0.25,0.75,0.90)))

MPA_effect_summary$family <- rownames(MPA_effect_summary)

MPA_effect_summary <- merge(MPA_effect_summary, family_beauty, by="family")
MPA_effect_summary <- MPA_effect_summary %>%
  arrange((Estimate))
MPA_effect_summary$log_beauty <- log(MPA_effect_summary$aesthe_score)

# # TAKE ONLY THE FAMILIES WITH HIGH MPA EFFECT (TOP X%)
# top_MPA_families <- MPA_effect_summary %>%
#   filter(Estimate >= quantile(MPA_effect_summary$Estimate, prob=0.75))
# 
# bottom_MPA_families <- MPA_effect_summary %>%
#   filter(Estimate <= quantile(MPA_effect_summary$Estimate, prob=0.25))

# TAKE ONLY THE FAMILIES WITH HIGH POS OR NEG MPA EFFECT (TOP X%)
top_MPA_families <- MPA_effect_summary %>%
  filter(abs_Estimate >= quantile(MPA_effect_summary$abs_Estimate, prob=0.75))

#########################################################
## FOREST PLOT OF MPA EFFECT COLORED BY AVERAGE BEAUTY ##
#########################################################

top_MPA_families <- top_MPA_families %>%
  arrange(aesthe_score)

source("R/variablecol.R")

plot_colors <- variablecol(colvar = top_MPA_families$aesthe_score, 
                           col = jet(n=nrow(top_MPA_families)), 
                           clim=range(family_beauty$aesthe_score))

graphics.off()
par(mar=c(4,12,4,4))
scatter2D(top_MPA_families$Estimate, seq(1:nrow(top_MPA_families)), xlim=c(min(top_MPA_families$Q10,na.rm = TRUE),max(top_MPA_families$Q90,na.rm = TRUE)),
          ylim=c(min(seq(1:nrow(top_MPA_families))-0.25),max(seq(1:nrow(top_MPA_families))+0.25)),cex=0,
          xlab="MPA Effect Size", ylab=NA, yaxt = "n",
          cex.lab=1.25, colvar = top_MPA_families$aesthe_score, col=jet(n=nrow(top_MPA_families)))
mtext(side=4, "Average Aesthetic Score", line=2, cex=1.1)
title("", line=1,
      font.main=1, cex.main=1.5)

par(lend=1)
x0 <- top_MPA_families$Q10
x1 <- top_MPA_families$Q90
y0 <- seq(1:nrow(top_MPA_families))
y1 <- seq(1:nrow(top_MPA_families))
segments(x0,y0,x1,y1, lwd=2, col=plot_colors)
x0 <- top_MPA_families$Q25
x1 <- top_MPA_families$Q75
y0 <- seq(1:nrow(top_MPA_families))
y1 <- seq(1:nrow(top_MPA_families))
segments(x0,y0,x1,y1, lwd=5, col=plot_colors)

points(top_MPA_families$Estimate, seq(1:nrow(top_MPA_families)),pch=21,col=1,bg=plot_colors, cex=1.5)
abline(v=0, lwd=1.5, col=adjustcolor("black",alpha.f = 0.75), lty=2)

axis(2, at = seq(1:nrow(top_MPA_families)), 
     labels = top_MPA_families$family,
     las=2, cex.axis=1)

title("Family Contributions to MPA Effect",
      line=1, font.main=1, cex.main=1.5)






