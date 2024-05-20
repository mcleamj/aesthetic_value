#######################################################################################
#'  BAYSESIAN HIERARCHICAL MODELS TO ASSESS DRIVERS OF AESTHETIC VALUE
#'  BASED ON DIRECTED ACYCLIC GRAPHS (DAG)
#'  THE CORRESPONDING DAG CAN BE FOUND AT:
#'  
#'  http://dagitty.net/mSNlypU
#'  
#' @author Matthew McLean, \email {mcleamj@gmail.com},
#'         
#' @date
########################################################################################

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
if(!require(tibble)){install.packages("tibble"); library(tibble)}
if(!require(dplyr)){install.packages("dplyr"); library(dplyr)}

##############################
## IMPORT DATA "MODEL DATA" ##
##############################

model_data <- readRDS("data/model_data.rds")

names(model_data)

model_data$MPA <- as.factor(model_data$MPA)

######################################
## IMPORT AND ADD DIVERSITY METRICS ##
######################################

biodiv <- read_rds("outputs/survey_biodiversity.rds")

model_data <- merge(model_data, biodiv, by="SurveyID")

##################################
## ADD NEW IMPUTED BENTHIC DATA ##
##################################

benthic_PCA <- read_rds("outputs/RLS_benthic_PCA_imputed.rds")
benthic_PCA <- benthic_PCA %>%
  select(SurveyID, PC1_imputation, PC2_imputation) %>%
  dplyr::rename(PC1_imputed = PC1_imputation, PC2_imputed = PC2_imputation)

model_data <- merge(model_data, benthic_PCA, by = "SurveyID", all=T)

#############################################
## REPLACE TEMPERATE BENTHIC 'DATA' WITH 0 ##
#############################################

model_data$PC1_imputed[model_data$Temperature_Zone=="Temperate"] <- 0
model_data$PC2_imputed[model_data$Temperature_Zone=="Temperate"] <- 0

#####################################
## IMPORT AND ATTACH AESTHETIC DATA #
#####################################

aesthetic_survey_data <- read.csv("outputs/survey_aesth.csv")

model_data <- merge(model_data, aesthetic_survey_data, by="SurveyID")

############################
## HOW MUCH MISSING DATA? ##
############################

sapply(model_data, function(x) paste(round(sum(is.na(x))/length(x),2)*100,"%",sep=""))

###############################################
## DO ANY VARIABLES NEED TO BE TRANSFORMED ? ##
###############################################

model_data <- do.call(data.frame,lapply(model_data, function(x) replace(x, is.infinite(x),NA)))

num_vars <- select_if(model_data, is.numeric)
num_vars$SurveyID <- NULL

# graphics.off()
# par(mfrow=c(4,4))
# for(i in 1:ncol(num_vars)){
#   hist(num_vars[,i], main=colnames(num_vars)[i])
# }

graphics.off()
par(mfrow=c(3,3))
for(i in 1:ncol(num_vars)){
  
  num_var_min <- min(num_vars[,i], na.rm=TRUE)
  
  num_var_log <- if(num_var_min > 0) { log(num_vars[,i])
  } else {
    log(num_vars[,i]+ceiling(abs(num_var_min))+1)
  }
  
  hist(num_vars[,i], main=NA)
  title(colnames(num_vars)[i])
  
  hist(num_var_log, main=NA)
  title(paste("log", colnames(num_vars)[i]))
  
  hist((num_vars[,i]^2), main=NA)
  title(paste("square", colnames(num_vars)[i]))
  
}

# NB_SPECIES, TAXO_ENTROPY, GRAVITY, WAVE ENERGY, PHOSPHATE, NITRATE, NPP, BIOMASS ALL RIGHT-SKEWED, SHOULD BE LOG-TRANSFORMED
# DEPTH AND DHW COULD POTENTIALLY BE TRANSFORMED AS WELL

# Check minimum values
dplyr::summarise(model_data,
                 dplyr::across(c(nb_species, taxo_entropy, 
                                 gravtot2, wave_energy, 
                                 BO_phosphate, BO_nitrate,
                                 NPP_mean, Biomass,
                                 Depth, dhw_mean),
                               min, na.rm=TRUE))
# Transform data
model_data <- dplyr::mutate(model_data,
                             dplyr::across(c(nb_species, taxo_entropy, 
                                             gravtot2, wave_energy, 
                                             BO_phosphate, BO_nitrate,
                                             NPP_mean, Biomass),
                                           log))

# CREATE ABSOLUTE VALUE OF LATITUDE
model_data <- dplyr::mutate(model_data,
                             abs_latitude = abs(SiteLatitude))


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


#'########################################################
#' NOW RUN THE INDIVIDUAL MODELS FOR EACH PREDICTOR
#' ACCORDING TO THE DAG (USING MINIMAL ADJUSTMENT SET)
#'#######################################################

################
## SST  MODEL ##
################

sst_model_formula <- bf(log(aesthe_survey) ~ sst_mean +
                          abs_latitude +
                          as.factor(Temperature_Zone) +
                          (1 | Country/SiteCode),
  
                        family=gaussian())

sst_model <- brm(sst_model_formula,
                 data=standardized_data,
                 chains=4, iter=4000, cores=ncores,
                 c(set_prior("normal(0,3)", class = "b"),
                   set_prior("normal(0,3)", class="Intercept")))

saveRDS(sst_model, "outputs/BIG_FILES/sst_model.rds")
sst_model <- read_rds("outputs/BIG_FILES/sst_model.rds")

sst_post <- as.data.frame(as.matrix(sst_model)) %>%
  select('b_sst_mean')

###################
## GRAVITY MODEL ##
###################

gravity_model_formula <- bf(log(aesthe_survey) ~ gravtot2 +
                              HDI2017 +
                              abs_latitude + # ABSOLUTE VALUE FOR LINEARITY
                              as.factor(Temperature_Zone) +
                              (1 | Country/SiteCode),
                            family=gaussian())

gravity_model <- brm(gravity_model_formula,
                     data=standardized_data,
                     chains=4, iter=4000, cores=ncores,
                     c(set_prior("normal(0,3)", class = "b"),
                       set_prior("normal(0,3)", class="Intercept")))

saveRDS(gravity_model, "outputs/BIG_FILES/gravity_model.rds")
gravity_model <- read_rds("outputs/BIG_FILES/gravity_model.rds")

gravity_post <- as.data.frame(as.matrix(gravity_model)) %>%
  select('b_gravtot2')

###############
## NPP MODEL ##
###############

NPP_model_formula <- bf(log(aesthe_survey) ~ NPP_mean +
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
                 c(set_prior("normal(0,3)", class = "b"),
                   set_prior("normal(0,3)", class="Intercept")))

saveRDS(NPP_model, "outputs/BIG_FILES/NPP_model.rds")
NPP_model <- read_rds("outputs/BIG_FILES/NPP_model.rds")

NPP_post <- as.data.frame(as.matrix(NPP_model)) %>%
  select('b_NPP_mean')

#################
## DEPTH MODEL ##
#################

depth_model_formula <- bf(log(aesthe_survey) ~ Depth +
                            as.factor(Temperature_Zone) +
                            (1 | Country/SiteCode),
                            
                          family=gaussian())

depth_model <- brm(depth_model_formula,
                   data=standardized_data,
                   chains=4, iter=4000, cores=ncores,
                   c(set_prior("normal(0,3)", class = "b"),
                     set_prior("normal(0,3)", class="Intercept")))

saveRDS(depth_model, "outputs/BIG_FILES/depth_model.rds")
depth_model <- read_rds("outputs/BIG_FILES/depth_model.rds")

depth_post <- as.data.frame(as.matrix(depth_model)) %>%
  select('b_Depth')

#######################
## FISH DEPEND MODEL ##
#######################

fshd_model_formula <- bf(log(aesthe_survey) ~ fshD +
                           HDI2017 +
                           as.factor(Temperature_Zone) +
                           (1 | Country/SiteCode),
                         
                         family=gaussian())

fshd_model <- brm(fshd_model_formula,
                  data=standardized_data,
                  chains=4, iter=4000, cores=ncores,
                  c(set_prior("normal(0,3)", class = "b"),
                    set_prior("normal(0,3)", class="Intercept")))

saveRDS(fshd_model, "outputs/BIG_FILES/fshd_model.rds")
fshd_model <- read_rds("outputs/BIG_FILES/fshd_model.rds")

fshd_post <- as.data.frame(as.matrix(fshd_model)) %>%
  select('b_fshD')

#################
## HDI   MODEL ##
#################

HDI_model_formula <- bf(log(aesthe_survey) ~ HDI2017 +
                          abs_latitude + 
                          sst_mean +
                          as.factor(Temperature_Zone) +
                          (1 | Country/SiteCode),
                          
                        family=gaussian())

HDI_model <- brm(HDI_model_formula,
                 data=standardized_data,
                 chains=4, iter=4000, cores=ncores,
                 c(set_prior("normal(0,3)", class = "b"),
                   set_prior("normal(0,3)", class="Intercept")))

saveRDS(HDI_model, "outputs/BIG_FILES/HDI_model.rds")
HDI_model <- read_rds("outputs/BIG_FILES/HDI_model.rds")

HDI_post <- as.data.frame(as.matrix(HDI_model)) %>%
  select('b_HDI2017')

###############
## MPA MODEL ##
###############

MPA_model_formula  <- bf(log(aesthe_survey) ~ MPA +
                           HDI2017 +
                           gravtot2 +
                           as.factor(Temperature_Zone) +
                           (1 | Country/SiteCode),
                           
                         family=gaussian())

MPA_model <- brm(MPA_model_formula,
                 data=standardized_data,
                 chains=4, iter=4000, cores=ncores,
                 c(set_prior("normal(0,3)", class = "b"),
                   set_prior("normal(0,3)", class="Intercept")))

saveRDS(MPA_model, "outputs/BIG_FILES/MPA_model.rds")
MPA_model <- read_rds("outputs/BIG_FILES/MPA_model.rds")

MPA_post <- as.data.frame(as.matrix(MPA_model)) %>%
  select('b_MPANotake','b_MPARestrictedtake')

#################################
## BENTHIC COMPOSITION MODEL 1 ##
#################################

benthic_PC1_model_formula <- 
  bf(log(aesthe_survey) ~ PC1_imputed +
       Depth +
       abs_latitude +
       MPA + 
       NPP_mean +
       sst_mean +
       as.factor(Temperature_Zone) +
       (1 | Country/SiteCode),
       
     family=gaussian()) 

benthic_PC1_model <- brm(benthic_PC1_model_formula,
                         data=standardized_data,
                         chains=4, iter=4000, cores=ncores,
                         c(set_prior("normal(0,3)", class = "b"),
                           set_prior("normal(0,3)", class="Intercept")))

saveRDS(benthic_PC1_model, "outputs/BIG_FILES/benthic_PC1_model.rds")
benthic_PC1_model <- read_rds("outputs/BIG_FILES/benthic_PC1_model.rds")

benthic_PC1_post <- as.data.frame(as.matrix(benthic_PC1_model))  %>%
  select('b_PC1_imputed')

#################################
## BENTHIC COMPOSITION MODEL 2 ##
#################################

benthic_PC2_model_formula <- 
  bf(log(aesthe_survey) ~ PC2_imputed +
       Depth +
       abs_latitude +
       MPA + 
       NPP_mean +
       sst_mean +
       as.factor(Temperature_Zone) +
       (1 | Country/SiteCode),
       
     family=gaussian()) 

benthic_PC2_model <- brm(benthic_PC2_model_formula,
                         data=standardized_data,
                         chains=4, iter=4000, cores=ncores,
                         c(set_prior("normal(0,3)", class = "b"),
                           set_prior("normal(0,3)", class="Intercept")))

saveRDS(benthic_PC2_model, "outputs/BIG_FILES/benthic_PC2_model.rds")
benthic_PC2_model <- read_rds("outputs/BIG_FILES/benthic_PC2_model.rds")

benthic_PC2_post <- as.data.frame(as.matrix(benthic_PC2_model)) %>%
  select('b_PC2_imputed')

###############
## DHW MODEL ##
###############

DHW_model_formula <- bf(log(aesthe_survey) ~ dhw_mean +
                          abs_latitude +
                          sst_mean +
                          as.factor(Temperature_Zone) +
                          (1 | Country/SiteCode),
                          
                        family=gaussian())

DHW_model <- brm(DHW_model_formula,
                 data=standardized_data,
                 chains=4, iter=4000, cores=ncores,
                 c(set_prior("normal(0,3)", class = "b"),
                   set_prior("normal(0,3)", class="Intercept")))

saveRDS(DHW_model, "outputs/BIG_FILES/DHW_model.rds")
DHW_model <- read_rds("outputs/BIG_FILES/DHW_model.rds")

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

saveRDS(dag_output, "outputs/dag_output.rds")

dag_output <- read_rds("outputs/dag_output.rds")

dag_estimates <- data.frame(median=apply(dag_output, 2, median))
dag_estimates$abs_effect <- abs(dag_estimates$median)
dag_estimates <- dag_estimates %>%
  arrange(desc(abs_effect))

####################################
## RANK DAG ORDER BY EFFECT SIZE? ##
####################################

dag_output <- dag_output[,order(match(colnames(dag_output), rownames(dag_estimates)))]

#########################
## SIMPLE FOREST PLOTS ##
#########################

mcmc_areas_ridges(dag_output) +
  ggtitle("DAG-Based Models") +
  xlim(range(dag_output))

mcmc_intervals(dag_output) +
  ggtitle("DAG-Based Models") +
  xlim(range(dag_output))

#########################
## MODEL SUMMARY TABLE ##
#########################

posterior_summary(dag_output, probs=c(0.10,0.90))

####################
## MODEL CHECKING ##
####################

theme_set(theme_classic(base_size = 6))

# SST MODEL
performance::check_model(sst_model)#, check = "vif") # GOOD, POSTERIOR COULD BE BETTER, VIF HIGH DUE TO LATITUDE AND SST
graphics.off()
plot(sst_model, variable="b_sst_mean")
r2_bayes(sst_model)

# NPP MODEL
performance::check_model(NPP_model,check = "vif") # GOOD, POSTERIOR COULD BE BETTER, VIF HIGH DUE TO LATIUDE AND SST
graphics.off()
plot(NPP_model,variable="b_NPP_mean")
r2_bayes(NPP_model)

# DEPTH MODEL
performance::check_model(depth_model) # GOOD
graphics.off()
plot(depth_model, variable="b_Depth")
r2_bayes(depth_model)

# GRAVITY MODEL
performance::check_model(gravity_model) # GOOD
graphics.off()
plot(gravity_model, variable="b_gravtot2")
r2_bayes(gravity_model)

# MPA MODEL
performance::check_model(MPA_model) # GOOD
graphics.off()
plot(MPA_model, variable=c("b_MPANotake","b_MPARestrictedtake"))
r2_bayes(MPA_model)

# BENTHIC MODEL 1
performance::check_model(benthic_PC1_model) # GOOD, HIGH VIF FROM LATITUDE AND SST
graphics.off()
plot(benthic_PC1_model,variable="b_PC1_imputed")
r2_bayes(benthic_PC1_model)

# BENTHIC MODEL 2
performance::check_model(benthic_PC2_model) # GOOD, HIGH VIF FROM LAT AND SST
graphics.off()
plot(benthic_PC2_model, variable="b_PC2_imputed")
r2_bayes(benthic_PC2_model)

# DHW MODEL
performance::check_model(DHW_model) # GOOD, HIGH VIF FROM LAT AND SST
graphics.off()
plot(DHW_model, variable="b_dhw_mean")
r2_bayes(DHW_model)

# HDI MODEL
performance::check_model(HDI_model) # GOOD, HIGH VIF FROM LAT AND SST
graphics.off()
plot(HDI_model,variable="b_HDI2017")
r2_bayes(HDI_model)

# FISHERIES DEPENDENCY MODEL
performance::check_model(fshd_model) # GOOD
graphics.off()
plot(fshd_model, variable="b_fshD")
r2_bayes(fshd_model)


##############################################
## SENSITIVITY TESTS FOR SST AND NPP MODELS 
## THESE MODELS HAVE HIGH VIF
## DUE TO COLLINEARITY WITH LATITUDE
## CHECK MODEL RESULT WITH LATITUDE REMOVED
##############################################

sst_sens_formula <- bf(log(aesthe_survey) ~ sst_mean +
                         as.factor(Temperature_Zone) +
                         (1 | Country/SiteCode),
                       
                       family=gaussian())

sst_sense_model <- brm(sst_sens_formula,
                       data=standardized_data,
                       chains=4, iter=4000, cores=ncores,
                       c(set_prior("normal(0,3)", class = "b"),
                         set_prior("normal(0,3)", class="Intercept")))

saveRDS(sst_sense_model, "outputs/BIG_FILES/sst_sense_model.rds")

sst_sens_post <- as.data.frame(as.matrix(sst_sense_model)) %>%
  select('b_sst_mean')

performance::check_model(sst_sense_model, check="vif")

# RIDGE PLOT THE TWO POSTERIORS

sst_ridge <- data.frame(draws=c(sst_post$b_sst_mean, sst_sens_post$b_sst_mean),
                        model=rep(c("original","sensitivity"),each=nrow(sst_post)))

ggplot(sst_ridge, aes(x = draws, y = model)) +
  geom_density_ridges() 

###########################
## NPP SENSITIVITY MODEL ##
###########################

NPP_sens_formula <- bf(log(aesthe_survey) ~ NPP_mean +
                         BO_nitrate +
                         BO_phosphate +
                         sst_mean +
                         wave_energy +
                         as.factor(Temperature_Zone) +
                         (1 | Country/SiteCode),
                       
                       family=gaussian())

NPP_sens_model <- brm(NPP_sens_formula,
                      data=standardized_data,
                      chains=4, iter=4000, cores=ncores,
                      c(set_prior("normal(0,3)", class = "b"),
                        set_prior("normal(0,3)", class="Intercept")))

saveRDS(NPP_sens_model, "outputs/BIG_FILES/NPP_sens_model.rds")

NPP_sens_post <- as.data.frame(as.matrix(NPP_sens_model)) %>%
  select('b_NPP_mean')

performance::check_model(NPP_sens_model, check="vif")

# RIDGE PLOT THE TWO POSTERIORS

NPP_ridge <- data.frame(draws=c(NPP_post$b_NPP_mean, NPP_sens_post$b_NPP_mean),
                        model=rep(c("original","sensitivity"),each=nrow(NPP_post)))

ggplot(NPP_ridge, aes(x = draws, y = model)) +
  geom_density_ridges() 


########################################################
########################################################
## NOW THE MODELS WITH SPECIES RICHNESS CONTROLED FOR ##----------------------------------------
########################################################
########################################################

################
## SST  MODEL ##
################

sst_richness_formula <- bf(log(aesthe_survey) ~ sst_mean + nb_species + 
                          abs_latitude +
                          as.factor(Temperature_Zone) +
                          (1 | Country/SiteCode),
                        
                        family=gaussian())

sst_richness_model <- brm(sst_richness_formula,
                 data=standardized_data,
                 chains=4, iter=4000, cores=ncores,
                 c(set_prior("normal(0,3)", class = "b"),
                   set_prior("normal(0,3)", class="Intercept")))

saveRDS(sst_richness_model, "outputs/BIG_FILES/sst_richness_model.rds")
sst_richness_model <- read_rds("outputs/BIG_FILES/sst_richness_model.rds")

sst_richness_post <- as.data.frame(as.matrix(sst_richness_model)) %>%
  select('b_sst_mean')


###################
## GRAVITY MODEL ##
###################

gravity_richness_formula <- bf(log(aesthe_survey) ~ gravtot2 + nb_species +
                              HDI2017 +
                              abs_latitude + # ABSOLUTE VALUE FOR LINEARITY
                              as.factor(Temperature_Zone) +
                              (1 | Country/SiteCode),
                            family=gaussian())

gravity_richness_model <- brm(gravity_richness_formula,
                     data=standardized_data,
                     chains=4, iter=4000, cores=ncores,
                     c(set_prior("normal(0,3)", class = "b"),
                       set_prior("normal(0,3)", class="Intercept")))

saveRDS(gravity_richness_model, "outputs/BIG_FILES/gravity_richness_model.rds")
gravity_richness_model <- read_rds("outputs/BIG_FILES/gravity_richness_model.rds")

gravity_richness_post <- as.data.frame(as.matrix(gravity_richness_model)) %>%
  select('b_gravtot2')


###############
## NPP MODEL ##
###############

NPP_richness_formula <- bf(log(aesthe_survey) ~ NPP_mean + nb_species +
                          abs_latitude +
                          BO_nitrate +
                          BO_phosphate +
                          sst_mean +
                          wave_energy +
                          as.factor(Temperature_Zone) +
                          (1 | Country/SiteCode),
                        
                        family=gaussian())

NPP_richness_model <- brm(NPP_richness_formula,
                 data=standardized_data,
                 chains=4, iter=4000, cores=ncores,
                 c(set_prior("normal(0,3)", class = "b"),
                   set_prior("normal(0,3)", class="Intercept")))

saveRDS(NPP_richness_model, "outputs/BIG_FILES/NPP_richness_model.rds")
NPP_richness_model <- read_rds("outputs/BIG_FILES/NPP_richness_model.rds")

NPP_richness_post <- as.data.frame(as.matrix(NPP_richness_model)) %>%
  select('b_NPP_mean')

#################
## DEPTH MODEL ##
#################

depth_richness_formula <- bf(log(aesthe_survey) ~ Depth + nb_species +
                            as.factor(Temperature_Zone) +
                            (1 | Country/SiteCode),
                          
                          family=gaussian())

depth_richness_model <- brm(depth_richness_formula,
                   data=standardized_data,
                   chains=4, iter=4000, cores=ncores,
                   c(set_prior("normal(0,3)", class = "b"),
                     set_prior("normal(0,3)", class="Intercept")))

saveRDS(depth_richness_model, "outputs/BIG_FILES/depth_richness_model.rds")
depth_richness_model <- read_rds("outputs/BIG_FILES/depth_richness_model.rds")

depth_richness_post <- as.data.frame(as.matrix(depth_richness_model)) %>%
  select('b_Depth')

#######################
## FISH DEPEND MODEL ##
#######################

fshd_richness_formula <- bf(log(aesthe_survey) ~ fshD + nb_species +
                           HDI2017 +
                           #abs_latitude +
                           as.factor(Temperature_Zone) +
                           (1 | Country/SiteCode),
                         
                         family=gaussian())

fshd_richness_model <- brm(fshd_richness_formula,
                  data=standardized_data,
                  chains=4, iter=4000, cores=ncores,
                  c(set_prior("normal(0,3)", class = "b"),
                    set_prior("normal(0,3)", class="Intercept")))

saveRDS(fshd_richness_model, "outputs/BIG_FILES/fshd_richness_model.rds")
fshd_richness_model <- read_rds("outputs/BIG_FILES/fshd_richness_model.rds")

fshd_richness_post <- as.data.frame(as.matrix(fshd_richness_model)) %>%
  select('b_fshD')

#################
## HDI   MODEL ##
#################

HDI_richness_formula <- bf(log(aesthe_survey) ~ HDI2017 + nb_species +
                          abs_latitude + 
                            sst_mean +
                          as.factor(Temperature_Zone) +
                          (1 | Country/SiteCode),
                        
                        family=gaussian())

HDI_richness_model <- brm(HDI_richness_formula,
                 data=standardized_data,
                 chains=4, iter=4000, cores=ncores,
                 c(set_prior("normal(0,3)", class = "b"),
                   set_prior("normal(0,3)", class="Intercept")))

saveRDS(HDI_richness_model, "outputs/BIG_FILES/HDI_richness_model.rds")
HDI_richness_model <- read_rds("outputs/BIG_FILES/HDI_richness_model.rds")

HDI_richness_post <- as.data.frame(as.matrix(HDI_richness_model)) %>%
  select('b_HDI2017')

###############
## MPA MODEL ##
###############

MPA_richness_formula  <- bf(log(aesthe_survey) ~ MPA + nb_species +
                           HDI2017 +
                           gravtot2 +
                           as.factor(Temperature_Zone) +
                           (1 | Country/SiteCode),
                         
                         family=gaussian())

MPA_richness_model <- brm(MPA_richness_formula,
                 data=standardized_data,
                 chains=4, iter=4000, cores=ncores,
                 c(set_prior("normal(0,3)", class = "b"),
                   set_prior("normal(0,3)", class="Intercept")))

saveRDS(MPA_richness_model, "outputs/BIG_FILES/MPA_richness_model.rds")
MPA_richness_model <- read_rds("outputs/BIG_FILES/MPA_richness_model.rds")

MPA_richness_post <- as.data.frame(as.matrix(MPA_richness_model)) %>%
  select('b_MPANotake','b_MPARestrictedtake')

#################################
## BENTHIC COMPOSITION MODEL 1 ##
#################################

benthic_PC1_richness_formula <- 
  bf(log(aesthe_survey) ~ PC1_imputed + nb_species +
       Depth +
       abs_latitude +
       MPA + 
       NPP_mean +
       sst_mean +
       as.factor(Temperature_Zone) +
       (1 | Country/SiteCode),
     
     family=gaussian()) 

benthic_PC1_richness_model <- brm(benthic_PC1_richness_formula,
                         data=standardized_data,
                         chains=4, iter=4000, cores=ncores,
                         c(set_prior("normal(0,3)", class = "b"),
                           set_prior("normal(0,3)", class="Intercept")))

saveRDS(benthic_PC1_richness_model, "outputs/BIG_FILES/benthic_PC1_richness_model.rds")
benthic_PC1_richness_model <- read_rds("outputs/BIG_FILES/benthic_PC1_richness_model.rds")

benthic_PC1_richness_post <- as.data.frame(as.matrix(benthic_PC1_richness_model))  %>%
  select('b_PC1_imputed')

#################################
## BENTHIC COMPOSITION MODEL 2 ##
#################################

benthic_PC2_richness_formula <- 
  bf(log(aesthe_survey) ~ PC2_imputed + nb_species +
       Depth +
       abs_latitude +
       MPA + 
       NPP_mean +
       sst_mean +
       as.factor(Temperature_Zone) +
       (1 | Country/SiteCode),
     
     family=gaussian()) 

benthic_PC2_richness_model <- brm(benthic_PC2_richness_formula,
                         data=standardized_data,
                         chains=4, iter=4000, cores=ncores,
                         c(set_prior("normal(0,3)", class = "b"),
                           set_prior("normal(0,3)", class="Intercept")))

saveRDS(benthic_PC2_richness_model, "outputs/BIG_FILES/benthic_PC2_richness_model.rds")
benthic_PC2_richness_model <- read_rds("outputs/BIG_FILES/benthic_PC2_richness_model.rds")

benthic_PC2_richness_post <- as.data.frame(as.matrix(benthic_PC2_richness_model)) %>%
  select('b_PC2_imputed')

###############
## DHW MODEL ##
###############

DHW_richness_formula <- bf(log(aesthe_survey) ~ dhw_mean + nb_species +
                          abs_latitude +
                          sst_mean +
                          as.factor(Temperature_Zone) +
                          (1 | Country/SiteCode),
                        
                        family=gaussian())

DHW_richness_model <- brm(DHW_richness_formula,
                 data=standardized_data,
                 chains=4, iter=4000, cores=ncores,
                 c(set_prior("normal(0,3)", class = "b"),
                   set_prior("normal(0,3)", class="Intercept")))

saveRDS(DHW_richness_model, "outputs/BIG_FILES/DHW_richness_model.rds")
DHW_richness_model <- read_rds("outputs/BIG_FILES/DHW_richness_model.rds")

DHW_richness_post <- as.data.frame(as.matrix(DHW_richness_model)) %>%
  select('b_dhw_mean')


###############################################
## RECOMBINE DAG MODELS AND RENAME VARIABLES ##
###############################################

dag_output_richness <- data.frame(
  sst_richness_post, gravity_richness_post, NPP_richness_post, depth_richness_post,
  fshd_richness_post, HDI_richness_post, MPA_richness_post,
  benthic_PC1_richness_post, benthic_PC2_richness_post, DHW_richness_post)

names(dag_output_richness) <- gsub("b_", "", names(dag_output_richness))

dag_output_richness <- dag_output_richness %>%
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

saveRDS(dag_output_richness, "outputs/dag_output_richness.rds")
dag_output_richness <- read_rds("outputs/dag_output_richness.rds")

dag_estimates_richness <- data.frame(median=apply(dag_output_richness, 2, median))
dag_estimates_richness$abs_effect <- abs(dag_estimates_richness$median)
dag_estimates_richness <- dag_estimates_richness %>%
  arrange(desc(abs_effect))

dag_output_richness <- dag_output_richness[,order(match(colnames(dag_output_richness), rownames(dag_estimates_richness)))]

## MATCH ORDER TO AESTH VALUE MODELS?
dag_output_richness <- dag_output_richness[,order(match(colnames(dag_output_richness), rownames(dag_estimates)))]

#########################
## SIMPLE FOREST PLOTS ##
#########################

color_scheme_set("darkgray")
mcmc_areas_ridges(dag_output_richness) +
  ggtitle("X") +
  xlim(range(dag_output_richness))

#########################
## SIMPLE FOREST PLOTS ##
#########################

ggpubr::ggarrange(
  
  mcmc_intervals(dag_output) +
    ggtitle("Total Causal Effect") +
    xlim(range(dag_output)),
  
  mcmc_intervals(dag_output_richness) +
    ggtitle("Causal Effect After Richness") +
    xlim(range(dag_output_richness)),
  
  ncol = 2)

###########################################
## SIMPLE FOREST PLOTS SAME X AXIS SCALE ##
###########################################

ggpubr::ggarrange(
  
  mcmc_intervals(dag_output) +
    ggtitle("Total Causal Effect") +
    xlim(range(dag_output)),
  
  mcmc_intervals(dag_output_richness) +
    ggtitle("Causal Effect After Richness") +
    xlim(range(dag_output)),
  
  ncol = 2)

####################
## MODEL CHECKING ##
####################

# SST MODEL
performance::check_model(sst_richness_model) # 
graphics.off()
plot(sst_richness_model)
r2_bayes(sst_richness_model)

# NPP MODEL
performance::check_model(NPP_richness__model) # 
graphics.off()
plot(NPP_richness_model)
r2_bayes(NPP_richness_model)

# DEPTH MODEL
performance::check_model(depth_richness_model) # 
graphics.off()
plot(depth_richness_model)
r2_bayes(depth_richness_model)

# BIOMASS MODEL
performance::check_model(Biom_richness_model) # 
graphics.off()
plot(Biom_richness_model)
r2_bayes(Biom_richness_model)

# GRAVITY MODEL
performance::check_model(gravity_richness_model) # 
graphics.off()
plot(gravity_richness_model)
r2_bayes(gravity_richness_model)

# MPA MODEL
performance::check_model(MPA_richness_model) # 
graphics.off()
plot(MPA_richness_model)
r2_bayes(MPA_richness_model)

# BENTHIC MODEL 1
performance::check_model(benthic_PC1_richness_model) # 
graphics.off()
plot(benthic_PC1_richness_model)
r2_bayes(benthic_PC1_richness_model)

# BENTHIC MODEL 2
performance::check_model(benthic_PC2_richness_model) # 
graphics.off()
plot(benthic_PC2_richness_model)
r2_bayes(benthic_PC2_richness_model)

# DHW MODEL
performance::check_model(DHW_richness_model) # 
graphics.off()
plot(DHW_richness_model)
r2_bayes(DHW_richness_model)

# HDI MODEL
performance::check_model(HDI_richness_model) # 
graphics.off()
plot(HDI_richness_model)
r2_bayes(HDI_richness_model)

# FISHERIES DEPENDENCY MODEL
performance::check_model(fshd_richness_model) # 
graphics.off()
plot(fshd_richness_model)
r2_bayes(fshd_richness_model)

