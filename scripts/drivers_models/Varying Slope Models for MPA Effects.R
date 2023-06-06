
####################################
## TEST FOR SPATIAL EFFECT OF MPA 
## FIRST USING VARYING SLOPES MODELS
## VARYING SLOPE BY COUNTRY AND
## BY ECOREGION 
####################################

####################################
## DEFINE VARIABLE COLOR FUNCTION ##
####################################

variablecol <- function(colvar, col, clim) {
  ncol <- length(col)
  colvar[colvar < min(clim)] <- NA
  colvar[colvar > max(clim)] <- NA
  rn <- clim[2] - clim[1]
  ifelse (rn != 0, Col <- col[1 + trunc((colvar - clim[1])/rn *
                                          (ncol - 1)+1e-15)], Col <- rep(col[1], ncol))               
  return(Col)
}

######################
## LIBRARY PACKAGES ##
######################

if(!require(data.table)){install.packages("data.table"); library(data.table)}
if(!require(rgdal)){install.packages("rgdal"); library(rgdal)}
if(!require(mapproj)){install.packages("mapproj"); library(mapproj)}
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
## IMPORT DATA "MODEL DATA" ##
##############################

model_data <- readRDS("data/model_data.rds")

names(model_data)

model_data$MPA <- as.factor(model_data$MPA)

##################################
## ADD NEW IMPUTED BENTHIC DATA ##
##################################

benthic_PCA <- read_rds("outputs/RLS_benthic_PCA_imputed.rds")
benthic_PCA <- benthic_PCA %>%
  select(SurveyID, PC1_imputation, PC2_imputation) %>%
  rename(PC1_imputed = PC1_imputation, PC2_imputed = PC2_imputation)

model_data <- merge(model_data, benthic_PCA, by = "SurveyID", all=T)


#####################################
## IMPORT AND ATTACH AESTHETIC DATA #
#####################################

aaesthe_surveyetic_data <- read.csv("outputs/survey_aesth.csv")

model_data <- merge(model_data, aaesthe_surveyetic_data, by="SurveyID")

##################################
## PARALLEL SETTINGS FOR MODELS ##
##################################

ncores = detectCores()
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#####################################################
## SCALE ALL THE NUMERIC PREDICTORS TO MEAN 0 SD 2 ##
#####################################################

z_score_2sd <- function(x){ (x - mean(x,na.rm=T)) / (2*sd(x,na.rm=T))}

z_vars <- model_data %>%
  select_if(is.numeric) %>%
  select(-any_of(c("SurveyID","SiteLongitude","SiteLatitude","aesthe_survey"))) %>%
  colnames()

standardized_data <- model_data %>%
  mutate_if(colnames(model_data) %in% z_vars, z_score_2sd)


###############################################
## MODEL USING SAME FORMAT FOR EFFECT OF MPA
## BUT WITH VARYING SLOPE FOR COUNTRY
###############################################

MPA_model_formula  <- bf(log(aesthe_survey) ~ MPA +
                           HDI2017 +
                           fshD +
                           gravtot2 +
                           as.factor(Temperature_Zone) +
                           #(1 | SiteCode) +
                           (1 + MPA | Country),
                           
                         family=gaussian())

MPA_model <- brm(MPA_model_formula,
                 data=standardized_data,
                 chains=3, iter=3000, cores=ncores,
                 c(set_prior("normal(0,3)", class = "b"),
                   set_prior("normal(0,3)", class="Intercept")))

MPA_post <- as.data.frame(as.matrix(MPA_model)) %>%
  select('b_MPANotake','b_MPARestrictedtake')



##############################################
## MODEL USING SAME FORMAT FOR EFFECT OF MPA
## BUT WITH VARYING SLOPE FOR ECOREGION
##############################################

ecoregion_slope_model_formula  <- bf(log(aesthe_survey) ~ MPA +
                           HDI2017 +
                           fshD +
                           gravtot2 +
                           as.factor(Temperature_Zone) +
                           (1 | SiteCode) +
                           (1 + MPA | Ecoregion),
                         
                         family=gaussian())

ecoregion_slope_model <- brm(ecoregion_slope_model_formula,
                 data=standardized_data,
                 chains=1, iter=2000, cores=ncores,
                 c(set_prior("normal(0,3)", class = "b"),
                   set_prior("normal(0,3)", class="Intercept")))
saveRDS(ecoregion_slope_model, "ecoregion_slope_model.rds")

eco_model_post <- as.data.frame(as.matrix(ecoregion_slope_model)) %>%
  #select('b_MPANotake','b_MPARestrictedtake')
  select('b_MPANotake':'b_gravtot2')

mcmc_intervals(eco_model_post)

eco_fix_ef <- fixef(ecoregion_slope_model)

eco_ran_ef <- ranef(ecoregion_slope_model)

eco_no_take <- as.data.frame(eco_ran_ef$Ecoregion[,,"MPANotake"])
eco_no_take$ECOREGION <- rownames(eco_no_take)

hist(eco_no_take$Estimate)




##########################################
## MAP OF MPA EFFECT FOR EACH ECOREGION ##
##########################################

ecoregion_info <- model_data %>%
  select(Ecoregion, SiteLongitude, SiteLatitude, PC1_imputed, PC2_imputed, nb_species, sst_mean) %>%
  group_by(Ecoregion) %>%
  summarise_all(.funs=mean) %>%
  rename(ECOREGION = Ecoregion)

regions <- readOGR("C:/Users/Matt/Documents/Biogeographic Trait Convergence Paper/Data & Code/MEOW-TNC", "meow_ecos")
regions <- regions[regions$ECOREGION %in% eco_no_take$ECOREGION,]

regions <- merge(regions, eco_no_take, by="ECOREGION")
regions <- merge(regions, ecoregion_info, by="ECOREGION")

col_range <- c( max(abs(regions$Estimate))*-1 , max(abs(regions$Estimate)) )

regions$color <- variablecol(regions$Estimate, col=(jet(nrow(regions))),
                             clim=col_range)

scatter2D(regions$SiteLongitude, regions$SiteLatitude, pch=19,
          colvar = regions$Estimate, clim=col_range,
          cex=0,xlim=c(-180,180),ylim=c(-60,80),
          xlab="Longitude (°)",ylab="Latitude (°)",
          cex.axis=1.25,cex.lab=1.25)
plot(regions, border="black", col=regions$color,add=TRUE)

map(database = "world",fill = TRUE, border="grey80", 
    col="grey80", ylim=c(-70,70),add=TRUE)
eco_corners <- par("usr")
polygon(c(eco_corners[1], eco_corners[1], eco_corners[2], eco_corners[2]),
        c(eco_corners[3], eco_corners[4], eco_corners[4], eco_corners[3]),lwd=1.35)
title("Ecoregion-Level MPA Effect")

plot(regions$Estimate ~ regions$SiteLatitude, pch=19,
     xlab="Latitude", ylab="MPA Effect Size")
title("Ecoregion-Level MPA Effect vs Latitude")

plot(regions$Estimate ~ regions$nb_species, pch=19,
     xlab="Species Richness", ylab="MPA Effect Size")
title("Ecoregion-Level MPA Effect vs Species Richness")

plot(regions$Estimate ~ regions$PC1_imputed, pch=19,
     xlab="x", ylab="MPA Effect Size")
title("Ecoregion-level MPA Effect vs. Benthic Composition")

plot(regions$Estimate ~ regions$PC2_imputed, pch=19,
     xlab="x", ylab="MPA Effect Size")
title("Ecoregion-level MPA Effect vs. Benthic Composition")


