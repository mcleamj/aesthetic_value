
#######################################################################################
#'  BAYSESIAN HIERARCHICAL MODEL OF THE EFFECT OF MPA ON AESTHETIC VALUE
#'  WITH VARYING SLOPES FOR COUNTRY OR ECOREGION
#'  TO ASSESS WHETHER THE EFFECT OF MPA VARIES ACROSS SPACE (LATITUDE)
#'
#' @author Matthew McLean, \email {mcleamj@gmail.com},
#'         
#' @date JUNE 9, 2023
########################################################################################

######################
## LIBRARY PACKAGES ##
######################

if(!require(data.table)){install.packages("data.table"); library(data.table)}
if(!require(mapproj)){install.packages("mapproj"); library(mapproj)}
if(!require(plot3D)){install.packages("plot3D"); library(plot3D)}
if(!require(ggplot2)){install.packages("ggplot2"); library(ggplot2)}
if(!require(pals)){install.packages("pals"); library(pals)}
if(!require(brms)){install.packages("brms"); library(brms)}
if(!require(bayesplot)){install.packages("bayesplot"); library(bayesplot)}
if(!require(parallel)){install.packages("parallel"); library(parallel)}
if(!require(rstan)){install.packages("rstan"); library(rstan)}
if(!require(cmdstanr)){install.packages("cmdstanr"); library(cmdstanr)}
if(!require(bayestestR)){install.packages("bayestestR"); library(bayestestR)}
if(!require(performance)){install.packages("performance"); library(performance)}
if(!require(readr)){install.packages("readr"); library(readr)}
if(!require(tibble)){install.packages("tibble"); library(tibble)}
if(!require(stringr)){install.packages("stringr"); library(stringr)}
if(!require(dplyr)){install.packages("dplyr"); library(dplyr)}

#####################################
## IMPORT DATA "STANDARDIZED DATA" ##
#####################################

standardized_data <- read_rds("outputs/standardized_data.rds")

##################################
## PARALLEL SETTINGS FOR MODELS ##
##################################

ncores = detectCores()
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#########################
## VARYING SLOPE MODEL ##
#########################

ecoregion_model_formula  <- bf(log(aesthe_survey_abund) ~ MPA + # TEST VARIABLE
                           HDI2017 + # CONTROL VARIABLE
                           gravtot2 + # CONTROL VARIABLE
                           as.factor(Temperature_Zone) + # CONTROL VARIABLE
                            # s(SiteLongitude, SiteLatitude, bs = "gp") + # SPATIAL GUASSIAN PROCESS
                           (1 | Country/SiteCode) + # RANDOM INTERCEPTS FOR SITES
                           (1 + MPA | Ecoregion), # RANDOM INTERCEPT AND SLOPE FOR ECOREGION
                         
                         family=gaussian())

ecoregion_model <- brm(ecoregion_model_formula,
                 data=standardized_data,
                 chains=4, iter=4000, cores=ncores,
                 backend = "cmdstanr",
                 #control = list(max_treedepth = 15),
                 c(set_prior("normal(0,3)", class = "b"),
                   set_prior("normal(0,3)", class="Intercept")))

saveRDS(ecoregion_model, "outputs/BIG_FILES/ecoregion_varying_slopes_model.rds")
ecoregion_model <- read_rds("outputs/BIG_FILES/ecoregion_varying_slopes_model.rds")

eco_model_post <- as.data.frame(as.matrix(ecoregion_model)) %>%
  #select('b_MPANotake','b_MPARestrictedtake')
  dplyr::select('b_MPANotake':'b_gravtot2')

mcmc_intervals(eco_model_post) # QUICK LOOK AT MODEL OUTPUT

################################
## EXTRACT MODEL COEFFICIENTS ##
################################

eco_no_take <- as.data.frame(coef(ecoregion_model)$Ecoregion[,,"MPANotake"]) %>%
  rename(Slope = "Estimate")
eco_no_take$ECOREGION <- rownames(eco_no_take)

# HISTOGRAM OF EFFECT SIZES
hist(eco_no_take$Slope)
abline(v=0, lty=2, lwd=2)

##############################
## CALCULATE PROJECTED GAIN ##
##############################

gain <- data.frame(Ecoregion = eco_intercept$ECOREGION,
                   fished = exp(eco_intercept$Intercept), 
                   no_take = exp(eco_no_take$Slope + eco_intercept$Intercept))
gain$diff <- gain$no_take - gain$fished
gain$perc_diff <- ((gain$no_take - gain$fished) / gain$fished)*100

mean(gain$perc_diff[gain$perc_diff>0])
