
##########################################
## SENSITIVITY TESTS FOR EXTREME PRIORS ##
##########################################

######################
## LIBRARY PACKAGES ##
######################

if(!require(data.table)){install.packages("data.table"); library(data.table)}
if(!require(mapproj)){install.packages("mapproj"); library(mapproj)}
if(!require(ggridges)){install.packages("ggridges"); library(ggridges)}
if(!require(ggplot2)){install.packages("ggplot2"); library(ggplot2)}
if(!require(plot3D)){install.packages("plot3D"); library(plot3D)}
if(!require(pals)){install.packages("pals"); library(pals)}
if(!require(brms)){install.packages("brms"); library(brms)}
if(!require(bayesplot)){install.packages("bayesplot"); library(bayesplot)}
if(!require(parallel)){install.packages("parallel"); library(parallel)}
if(!require(rstan)){install.packages("rstan"); library(rstan)}
if(!require(bayestestR)){install.packages("bayestestR"); library(bayestestR)}
if(!require(readr)){install.packages("readr"); library(readr)}
if(!require(performance)){install.packages("performance"); library(performance)}
if(!require(glmmTMB)){install.packages("glmmTMB"); library(glmmTMB)}
if(!require(tibble)){install.packages("tibble"); library(tibble)}
if(!require(dplyr)){install.packages("dplyr"); library(dplyr)}

#############################################
## IMPORT PREPARED STANDARDIZED MODEL DATA ##
#############################################

standardized_data <-  read_rds("outputs/standardized_data.rds")

##################################
## PARALLEL SETTINGS FOR MODELS ##
##################################

ncores = detectCores()
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

######################################
## SET GLOBAL PRIORS FOR ALL MODELS ##
######################################

prior <- c(set_prior("normal(0,0.1)", class = "b"),
           set_prior("normal(0,0.1)", class = "Intercept"))

control = list(
  adapt_delta = 0.95,            # Higher adapt_delta for stability
  max_treedepth = 15)            # Increased tree depth


init = "0"                       # Initialize at zero


#'#######################################################
#' NOW RUN THE INDIVIDUAL MODELS FOR EACH PREDICTOR
#' ACCORDING TO THE DAG (USING MINIMAL ADJUSTMENT SET)
#'#######################################################

################
## SST  MODEL ##
################

sst_model_formula <- bf(log(aesthe_survey_abund) ~ sst_mean +
                          abs_latitude +
                          as.factor(Temperature_Zone) +
                          (1 | Country/SiteCode),
                        
                        
                        family=gaussian())

sst_model <- brm(sst_model_formula,
                 data=standardized_data,
                 chains=4, iter=4000, cores=ncores,
                 prior=prior, control = control, init = init)

saveRDS(sst_model, "outputs/BIG_FILES/sst_model_strict.rds")
sst_model <- read_rds("outputs/BIG_FILES/sst_model_strict.rds")

sst_post <- as.data.frame(as.matrix(sst_model)) %>%
  select('b_sst_mean')

###################
## GRAVITY MODEL ##
###################

gravity_model_formula <- bf(log(aesthe_survey_abund) ~ gravtot2 +
                              HDI2017 +
                              abs_latitude + # ABSOLUTE VALUE FOR LINEARITY
                              as.factor(Temperature_Zone) +
                              (1 | Country/SiteCode),                           
                            
                            family=gaussian())

gravity_model <- brm(gravity_model_formula,
                     data=standardized_data,
                     chains=4, iter=4000, cores=ncores,
                     prior=prior, control = list(max_treedepth = 12))

saveRDS(gravity_model, "outputs/BIG_FILES/gravity_model_strict.rds")
gravity_model <- read_rds("outputs/BIG_FILES/gravity_model_strict.rds")

gravity_post <- as.data.frame(as.matrix(gravity_model)) %>%
  select('b_gravtot2')

###############
## NPP MODEL ##
###############

NPP_model_formula <- bf(log(aesthe_survey_abund) ~ NPP_mean +
                          abs_latitude +
                          BO_nitrate +
                          BO_phosphate +
                          sst_mean +
                          wave_energy +
                          as.factor(Temperature_Zone) +
                          (1 | Country/SiteCode),                           
                        
                        
                        family=gaussian())

NPP_model <- brm(NPP_model_formula,
                 data=standardized_data,
                 chains=4, iter=4000, cores=ncores,
                 prior=prior, control = list(max_treedepth = 12))

saveRDS(NPP_model, "outputs/BIG_FILES/NPP_model_strict.rds")
NPP_model <- read_rds("outputs/BIG_FILES/NPP_model_strict.rds")

NPP_post <- as.data.frame(as.matrix(NPP_model)) %>%
  select('b_NPP_mean')

#################
## DEPTH MODEL ##
#################

depth_model_formula <- bf(log(aesthe_survey_abund) ~ Depth +
                            as.factor(Temperature_Zone) +
                            (1 | Country/SiteCode),                          
                          
                          
                          family=gaussian())

depth_model <- brm(depth_model_formula,
                   data=standardized_data,
                   chains=4, iter=4000, cores=ncores,
                   prior=prior, control = list(max_treedepth = 12))

saveRDS(depth_model, "outputs/BIG_FILES/depth_model_strict.rds")
depth_model <- read_rds("outputs/BIG_FILES/depth_model_strict.rds")

depth_post <- as.data.frame(as.matrix(depth_model)) %>%
  select('b_Depth')

#######################
## FISH DEPEND MODEL ##
#######################

fshd_model_formula <- bf(log(aesthe_survey_abund) ~ fshD +
                           HDI2017 +
                           as.factor(Temperature_Zone) +
                           (1 | Country/SiteCode),                          
                         
                         
                         family=gaussian())

fshd_model <- brm(fshd_model_formula,
                  data=standardized_data,
                  chains=4, iter=4000, cores=ncores,
                  prior=prior, control = list(max_treedepth = 12))

saveRDS(fshd_model, "outputs/BIG_FILES/fshd_model_strict.rds")
fshd_model <- read_rds("outputs/BIG_FILES/fshd_model_strict.rds")

fshd_post <- as.data.frame(as.matrix(fshd_model)) %>%
  select('b_fshD')

#################
## HDI   MODEL ##
#################

HDI_model_formula <- bf(log(aesthe_survey_abund) ~ HDI2017 +
                          abs_latitude + 
                          sst_mean +
                          as.factor(Temperature_Zone) +
                          (1 | Country/SiteCode),                           
                        
                        
                        family=gaussian())

HDI_model <- brm(HDI_model_formula,
                 data=standardized_data,
                 chains=4, iter=4000, cores=ncores,
                 prior=prior, control = list(max_treedepth = 12))

saveRDS(HDI_model, "outputs/BIG_FILES/HDI_model_strict.rds")
HDI_model <- read_rds("outputs/BIG_FILES/HDI_model_strict.rds")

HDI_post <- as.data.frame(as.matrix(HDI_model)) %>%
  select('b_HDI2017')

###############
## MPA MODEL ##
###############

MPA_model_formula  <- bf(log(aesthe_survey_abund) ~ MPA +
                           HDI2017 +
                           gravtot2 +
                           as.factor(Temperature_Zone) +
                           (1 | Country/SiteCode),                          
                         
                         
                         family=gaussian())

MPA_model <- brm(MPA_model_formula,
                 data=standardized_data,
                 chains=4, iter=4000, cores=ncores,
                 prior=prior, control = list(max_treedepth = 12))

saveRDS(MPA_model, "outputs/BIG_FILES/MPA_model_strict.rds")
MPA_model <- read_rds("outputs/BIG_FILES/MPA_model_strict.rds")

MPA_post <- as.data.frame(as.matrix(MPA_model)) %>%
  select('b_MPANotake','b_MPARestrictedtake')

#################################
## BENTHIC COMPOSITION MODEL 1 ##
#################################

benthic_PC1_model_formula <- 
  bf(log(aesthe_survey_abund) ~ PC1_imputed +
       Depth + 
       gravtot2 +
       MPA + 
       NPP_mean +
       sst_mean +
       as.factor(Temperature_Zone) +
       (1 | Country/SiteCode),                           
     
     
     family=gaussian()) 

benthic_PC1_model <- brm(benthic_PC1_model_formula,
                         data=standardized_data,
                         chains=4, iter=4000, cores=ncores,
                         prior=prior, control = list(max_treedepth = 12))

saveRDS(benthic_PC1_model, "outputs/BIG_FILES/benthic_PC1_model_strict.rds")
benthic_PC1_model <- read_rds("outputs/BIG_FILES/benthic_PC1_model_strict.rds")

benthic_PC1_post <- as.data.frame(as.matrix(benthic_PC1_model))  %>%
  select('b_PC1_imputed')

#################################
## BENTHIC COMPOSITION MODEL 2 ##
#################################

benthic_PC2_model_formula <- 
  bf(log(aesthe_survey_abund) ~ PC2_imputed +
       Depth + 
       gravtot2 +
       MPA + 
       NPP_mean +
       sst_mean +
       as.factor(Temperature_Zone) +
       (1 | Country/SiteCode),                          
     
     
     family=gaussian()) 

benthic_PC2_model <- brm(benthic_PC2_model_formula,
                         data=standardized_data,
                         chains=4, iter=4000, cores=ncores,
                         prior=prior, control = list(max_treedepth = 12))

saveRDS(benthic_PC2_model, "outputs/BIG_FILES/benthic_PC2_model_strict.rds")
benthic_PC2_model <- read_rds("outputs/BIG_FILES/benthic_PC2_model_strict.rds")

benthic_PC2_post <- as.data.frame(as.matrix(benthic_PC2_model)) %>%
  select('b_PC2_imputed')

###############
## DHW MODEL ##
###############

DHW_model_formula <- bf(log(aesthe_survey_abund) ~ dhw_mean +
                          abs_latitude +
                          sst_mean +
                          as.factor(Temperature_Zone) +
                          (1 | Country/SiteCode),                           
                        
                        
                        family=gaussian())

DHW_model <- brm(DHW_model_formula,
                 data=standardized_data,
                 chains=4, iter=4000, cores=ncores,
                 prior=prior, control = list(max_treedepth = 12))

saveRDS(DHW_model, "outputs/BIG_FILES/DHW_model_strict.rds")
DHW_model <- read_rds("outputs/BIG_FILES/DHW_model_strict.rds")

DHW_post <- as.data.frame(as.matrix(DHW_model)) %>%
  select('b_dhw_mean')

###############################################
## RECOMBINE DAG MODELS AND RENAME VARIABLES ##
###############################################

dag_output <- data.frame(sst_post, gravity_post, NPP_post, depth_post,
                         fshd_post, HDI_post, MPA_post,
                         benthic_PC1_post, benthic_PC2_post, DHW_post)

names(dag_output) <- gsub("b_", "", names(dag_output))

dag_output <- dag_output %>%
  rename("Log Human Gravity" = gravtot2,
         "No Take MPA" = MPANotake,
         "Restricted Take MPA" = MPARestrictedtake,
         "Sea Surface Temperature" = sst_mean,
         "Net Primary Productivity" = NPP_mean,
         "Degree Heating Weeks" = dhw_mean,
         "Depth" =  Depth,
         "Fisheries Dependency" = fshD,
         "Human Dev. Index" = HDI2017,
         "Benthic Composition (PC1)" = PC1_imputed,
         "Benthic Composition (PC2)" = PC2_imputed)

saveRDS(dag_output, "outputs/dag_output_strict.rds")


###################################################################################
###################################################################################
## NOW THE FLAT PRIORS
###################################################################################
###################################################################################


######################################
## SET GLOBAL PRIORS FOR ALL MODELS ##
## ###################################

prior <- c(set_prior("normal(0,100)", class = "b"),
           set_prior("normal(0,100)", class = "Intercept"))

#'#######################################################
#' NOW RUN THE INDIVIDUAL MODELS FOR EACH PREDICTOR
#' ACCORDING TO THE DAG (USING MINIMAL ADJUSTMENT SET)
#'#######################################################

################
## SST  MODEL ##
################

sst_model_formula <- bf(log(aesthe_survey_abund) ~ sst_mean +
                          abs_latitude +
                          as.factor(Temperature_Zone) +
                          (1 | Country/SiteCode),
                        
                        
                        family=gaussian())

sst_model <- brm(sst_model_formula,
                 data=standardized_data,
                 chains=4, iter=4000, cores=ncores,
                 prior=prior, control = list(max_treedepth = 12))

saveRDS(sst_model, "outputs/BIG_FILES/sst_model_flat.rds")
sst_model <- read_rds("outputs/BIG_FILES/sst_model_flat.rds")

sst_post <- as.data.frame(as.matrix(sst_model)) %>%
  select('b_sst_mean')

###################
## GRAVITY MODEL ##
###################

gravity_model_formula <- bf(log(aesthe_survey_abund) ~ gravtot2 +
                              HDI2017 +
                              abs_latitude + # ABSOLUTE VALUE FOR LINEARITY
                              as.factor(Temperature_Zone) +
                              (1 | Country/SiteCode),                           
                            
                            family=gaussian())

gravity_model <- brm(gravity_model_formula,
                     data=standardized_data,
                     chains=4, iter=4000, cores=ncores,
                     prior=prior, control = list(max_treedepth = 12))

saveRDS(gravity_model, "outputs/BIG_FILES/gravity_model_flat.rds")
gravity_model <- read_rds("outputs/BIG_FILES/gravity_model_flat.rds")

gravity_post <- as.data.frame(as.matrix(gravity_model)) %>%
  select('b_gravtot2')

###############
## NPP MODEL ##
###############

NPP_model_formula <- bf(log(aesthe_survey_abund) ~ NPP_mean +
                          abs_latitude +
                          BO_nitrate +
                          BO_phosphate +
                          sst_mean +
                          wave_energy +
                          as.factor(Temperature_Zone) +
                          (1 | Country/SiteCode),                           
                        
                        
                        family=gaussian())

NPP_model <- brm(NPP_model_formula,
                 data=standardized_data,
                 chains=4, iter=4000, cores=ncores,
                 prior=prior, control = list(max_treedepth = 12))

saveRDS(NPP_model, "outputs/BIG_FILES/NPP_model_flat.rds")
NPP_model <- read_rds("outputs/BIG_FILES/NPP_model_flat.rds")

NPP_post <- as.data.frame(as.matrix(NPP_model)) %>%
  select('b_NPP_mean')

#################
## DEPTH MODEL ##
#################

depth_model_formula <- bf(log(aesthe_survey_abund) ~ Depth +
                            as.factor(Temperature_Zone) +
                            (1 | Country/SiteCode),                          
                          
                          
                          family=gaussian())

depth_model <- brm(depth_model_formula,
                   data=standardized_data,
                   chains=4, iter=4000, cores=ncores,
                   prior=prior, control = list(max_treedepth = 12))

saveRDS(depth_model, "outputs/BIG_FILES/depth_model_flat.rds")
depth_model <- read_rds("outputs/BIG_FILES/depth_model_flat.rds")

depth_post <- as.data.frame(as.matrix(depth_model)) %>%
  select('b_Depth')

#######################
## FISH DEPEND MODEL ##
#######################

fshd_model_formula <- bf(log(aesthe_survey_abund) ~ fshD +
                           HDI2017 +
                           as.factor(Temperature_Zone) +
                           (1 | Country/SiteCode),                          
                         
                         
                         family=gaussian())

fshd_model <- brm(fshd_model_formula,
                  data=standardized_data,
                  chains=4, iter=4000, cores=ncores,
                  prior=prior, control = list(max_treedepth = 12))

saveRDS(fshd_model, "outputs/BIG_FILES/fshd_model_flat.rds")
fshd_model <- read_rds("outputs/BIG_FILES/fshd_model_flat.rds")

fshd_post <- as.data.frame(as.matrix(fshd_model)) %>%
  select('b_fshD')

#################
## HDI   MODEL ##
#################

HDI_model_formula <- bf(log(aesthe_survey_abund) ~ HDI2017 +
                          abs_latitude + 
                          sst_mean +
                          as.factor(Temperature_Zone) +
                          (1 | Country/SiteCode),                           
                        
                        
                        family=gaussian())

HDI_model <- brm(HDI_model_formula,
                 data=standardized_data,
                 chains=4, iter=4000, cores=ncores,
                 prior=prior, control = list(max_treedepth = 12))

saveRDS(HDI_model, "outputs/BIG_FILES/HDI_model_flat.rds")
HDI_model <- read_rds("outputs/BIG_FILES/HDI_model_flat.rds")

HDI_post <- as.data.frame(as.matrix(HDI_model)) %>%
  select('b_HDI2017')

###############
## MPA MODEL ##
###############

MPA_model_formula  <- bf(log(aesthe_survey_abund) ~ MPA +
                           HDI2017 +
                           gravtot2 +
                           as.factor(Temperature_Zone) +
                           (1 | Country/SiteCode),                          
                         
                         
                         family=gaussian())

MPA_model <- brm(MPA_model_formula,
                 data=standardized_data,
                 chains=4, iter=4000, cores=ncores,
                 prior=prior, control = list(max_treedepth = 12))

saveRDS(MPA_model, "outputs/BIG_FILES/MPA_model_flat.rds")
MPA_model <- read_rds("outputs/BIG_FILES/MPA_model_flat.rds")

MPA_post <- as.data.frame(as.matrix(MPA_model)) %>%
  select('b_MPANotake','b_MPAReflatedtake')

#################################
## BENTHIC COMPOSITION MODEL 1 ##
#################################

benthic_PC1_model_formula <- 
  bf(log(aesthe_survey_abund) ~ PC1_imputed +
       Depth + 
       gravtot2 +
       MPA + 
       NPP_mean +
       sst_mean +
       as.factor(Temperature_Zone) +
       (1 | Country/SiteCode),                           
     
     
     family=gaussian()) 

benthic_PC1_model <- brm(benthic_PC1_model_formula,
                         data=standardized_data,
                         chains=4, iter=4000, cores=ncores,
                         prior=prior, control = list(max_treedepth = 12))

saveRDS(benthic_PC1_model, "outputs/BIG_FILES/benthic_PC1_model_flat.rds")
benthic_PC1_model <- read_rds("outputs/BIG_FILES/benthic_PC1_model_flat.rds")

benthic_PC1_post <- as.data.frame(as.matrix(benthic_PC1_model))  %>%
  select('b_PC1_imputed')

#################################
## BENTHIC COMPOSITION MODEL 2 ##
#################################

benthic_PC2_model_formula <- 
  bf(log(aesthe_survey_abund) ~ PC2_imputed +
       Depth + 
       gravtot2 +
       MPA + 
       NPP_mean +
       sst_mean +
       as.factor(Temperature_Zone) +
       (1 | Country/SiteCode),                          
     
     
     family=gaussian()) 

benthic_PC2_model <- brm(benthic_PC2_model_formula,
                         data=standardized_data,
                         chains=4, iter=4000, cores=ncores,
                         prior=prior, control = list(max_treedepth = 12))

saveRDS(benthic_PC2_model, "outputs/BIG_FILES/benthic_PC2_model_flat.rds")
benthic_PC2_model <- read_rds("outputs/BIG_FILES/benthic_PC2_model_flat.rds")

benthic_PC2_post <- as.data.frame(as.matrix(benthic_PC2_model)) %>%
  select('b_PC2_imputed')

###############
## DHW MODEL ##
###############

DHW_model_formula <- bf(log(aesthe_survey_abund) ~ dhw_mean +
                          abs_latitude +
                          sst_mean +
                          as.factor(Temperature_Zone) +
                          (1 | Country/SiteCode),                           
                        
                        
                        family=gaussian())

DHW_model <- brm(DHW_model_formula,
                 data=standardized_data,
                 chains=4, iter=4000, cores=ncores,
                 prior=prior, control = list(max_treedepth = 12))

saveRDS(DHW_model, "outputs/BIG_FILES/DHW_model_flat.rds")
DHW_model <- read_rds("outputs/BIG_FILES/DHW_model_flat.rds")

DHW_post <- as.data.frame(as.matrix(DHW_model)) %>%
  select('b_dhw_mean')

###############################################
## RECOMBINE DAG MODELS AND RENAME VARIABLES ##
###############################################

dag_output <- data.frame(sst_post, gravity_post, NPP_post, depth_post,
                         fshd_post, HDI_post, MPA_post,
                         benthic_PC1_post, benthic_PC2_post, DHW_post)

names(dag_output) <- gsub("b_", "", names(dag_output))

dag_output <- dag_output %>%
  rename("Log Human Gravity" = gravtot2,
         "No Take MPA" = MPANotake,
         "Reflated Take MPA" = MPAReflatedtake,
         "Sea Surface Temperature" = sst_mean,
         "Net Primary Productivity" = NPP_mean,
         "Degree Heating Weeks" = dhw_mean,
         "Depth" =  Depth,
         "Fisheries Dependency" = fshD,
         "Human Dev. Index" = HDI2017,
         "Benthic Composition (PC1)" = PC1_imputed,
         "Benthic Composition (PC2)" = PC2_imputed)

saveRDS(dag_output, "outputs/dag_output_flat.rds")
