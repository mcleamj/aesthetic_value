
#######################################################################################
#'
#'  CODE TO PREDICT THE OUTCOMES OF CHANGING FISHED SITES TO MPA SITES
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

model_data$sst_range <- model_data$sst_max - model_data$sst_min

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
  rename(PC1_imputed = PC1_imputation, PC2_imputed = PC2_imputation)

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

model_data <- model_data %>%
  mutate_if(colnames(model_data) %in% z_vars, z_score_2sd)



##############################################
## IMPORT THE SAVED FULL CAUSAL SALAD MODEL ##
##############################################

## EXTRACT MODEL ESTIMATES AND INFORMATION

full_model <- readRDS("outputs/full_model.RDS")

full_post <- as.data.frame(as.matrix(full_model)) %>%
  select('b_sst_mean':'b_MPARestrictedtake')

full_estimates <- data.frame(median=apply(full_post, 2, median))
full_estimates$abs_effect <- abs(full_estimates$median)
full_estimates <- full_estimates %>%
  arrange(desc(abs_effect))

full_post <- full_post[,order(match(colnames(full_post), rownames(full_estimates)))]
mcmc_intervals(full_post)

#######################################################
## SELECT ONLY THE RELEVANT VARIABLES FROM THE MODEL 
## AND SELECT ONLY FISHED SITES 
####################################################

model_variables <- names(full_model$data)

new_data <- model_data %>%
  select(any_of(model_variables), "SiteLongitude", "SiteLatitude") %>%
  filter(MPA == "Fishing")

table(new_data$MPA)

###############################################
## CONVERT FISHED SITES TO MPA NO TAKE SITES ##
###############################################

new_data$MPA[new_data$MPA=="Fishing"] <- "No take"
new_data$MPA <- droplevels(new_data$MPA)
levels(new_data$MPA)

###############################################
## CHECK AND REMOVE MISSING VALUES IF NEEDED ##
###############################################

sapply(new_data, function(x) paste(round(sum(is.na(x))/length(x),2)*100,"%",sep=""))

new_data <- na.omit(new_data)

#####################################################
## MAKE PREDICTIONS FROM FULL MODEL USING NEW DATA ##
#####################################################

MPA_scenario <- as.data.frame(predict(full_model, newdata = new_data))

###############################
## CONVERT LOG OUTPUT TO RAW ##
###############################

graphics.off()

MPA_test <- exp(MPA_scenario$Estimate) # NEED TO ADD INTERCEPT?

MPA_delta <- MPA_test - new_data$aesthe_survey
hist(MPA_delta)

MPA_percent_delta <- ((MPA_test - new_data$aesthe_survey)/new_data$aesthe_survey)*100
hist(MPA_percent_delta)

scatter2D(new_data$SiteLongitude, new_data$SiteLatitude, pch=19, 
          colvar = MPA_percent_delta)

plot(MPA_percent_delta ~ new_data$SiteLatitude)
plot(MPA_percent_delta ~ new_data$SiteLongitude)
plot(MPA_percent_delta ~ new_data$aesthe_survey)


###################################
###################################
## MAKE PREDICTIONS BASED ON DHW ##
###################################
###################################

# DHW INCREASES BY 2 FOLD

##############################################
## NEED TO CHANGE MPA STATUS BACK TO FISHED ##
##############################################

new_data$MPA <- as.character(new_data$MPA)
new_data$MPA[new_data$MPA=="No take"] <- "Fishing"
new_data$MPA <- as.factor(new_data$MPA)
levels(new_data$MPA)

####################
## DOUBLE THE DHW ##
####################

new_data$dhw_mean <- new_data$dhw_mean * 2

###########################################
## PREDCIT THE CHANGE USING THE NEW DATA ##
###########################################

DHW_scenario <-  as.data.frame(predict(full_model, newdata = new_data))

###############################
## CONVERT LOG OUTPUT TO RAW ##
###############################

DHW_test <- exp(DHW_scenario$Estimate) # NEED TO ADD INTERCEPT?

DHW_delta <- DHW_test - new_data$aesthe_survey
hist(DHW_delta)

DHW_percent_delta <- ((DHW_test - new_data$aesthe_survey)/new_data$aesthe_survey)*100
hist(DHW_percent_delta)

scatter2D(new_data$SiteLongitude, new_data$SiteLatitude, pch=19, 
          colvar = DHW_percent_delta)

plot(DHW_percent_delta ~ new_data$SiteLatitude)
plot(DHW_percent_delta ~ new_data$SiteLongitude)
plot(DHW_percent_delta ~ new_data$aesthe_survey)



###################################################
###################################################
## MAKE PREDICTIONS BASED ON BENTHIC COMPOSITION ##
###################################################
###################################################

# HOW MUCH TO CHANGE PC1 BY?

####################################
## CAN ONLY DO FOR TROPICAL SITES ##
####################################

new_benthic_data <- new_data %>%
  filter(Temperature_Zone=="Tropical")

####################################
## MAKE SURE MPA STATUS IS FISHED ##
####################################

levels(new_data$MPA)

##########################################################################################
## SHIFT ALL PC1 VALUES TO POSITIVE SIDE BY XX% (MORE ALGAE, LESS CORAL) 
## FIGURE OUT THE CORRELATION BETWEEN RAW CORAL COVER AND PC 1 VALUE
## USE THAT TO TRY SCENARIOS, E.G., DECREASE CORAL COVER BY 10% = HOW MUCH CHANGE IN PC1 
##########################################################################################

benthic_coral_data <- data.frame(coral_cover = benthic_raw$coral_imputation, 
                                 PC1_imputed = benthic_PCA$PC1_imputed)
plot(PC1_imputed ~ coral_cover, data=benthic_coral_data)
coral_PC1_regression <- lm(PC1_imputed ~ coral_cover, data=benthic_coral_data)
summary(coral_PC1_regression)
# LETS SAY ALL SITES LOSE 25% OF THEIR EXISTING COVER 
new_coral_data <- as.data.frame(benthic_coral_data$coral_cover - (benthic_coral_data$coral_cover*0.25))
colnames(new_coral_data) <- "coral_cover"
hist(new_coral_data$coral_cover)
# PREDICT CHANGE IN PC1 FROM SIMPLE LINEAR MODEL
new_PC1 <- predict(coral_PC1_regression, newdata = new_coral_data)

hist(benthic_coral_data$PC1_imputed)
hist(new_PC1)

############################################
## REPLACE PC1 VALUES IN NEW BENTHIC DATA ##
############################################

new_benthic_data$PC1_imputed <- new_PC1

###########################################
## PREDCIT THE CHANGE USING THE NEW DATA ##
###########################################

DHW_scenario <-  as.data.frame(predict(full_model, newdata = new_data))

###############################
## CONVERT LOG OUTPUT TO RAW ##
###############################

DHW_test <- exp(DHW_scenario$Estimate) # NEED TO ADD INTERCEPT?

DHW_delta <- DHW_test - new_data$aesthe_survey
hist(DHW_delta)

DHW_percent_delta <- ((DHW_test - new_data$aesthe_survey)/new_data$aesthe_survey)*100
hist(DHW_percent_delta)

scatter2D(new_data$SiteLongitude, new_data$SiteLatitude, pch=19, 
          colvar = DHW_percent_delta)

plot(DHW_percent_delta ~ new_data$SiteLatitude)
plot(DHW_percent_delta ~ new_data$SiteLongitude)
plot(DHW_percent_delta ~ new_data$aesthe_survey)








