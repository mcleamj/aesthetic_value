
#######################################################################################
#' CODE FOR ORGANIZING AND PREPARING COVARIATE DATA PRIOR TO RUNNING BAYESIAN MODELS
#'  
#' @author Matthew McLean, \email {mcleamj@gmail.com},
#'         
#' @date
########################################################################################

######################
## LIBRARY PACKAGES ##
######################

if(!require(data.table)){install.packages("data.table"); library(data.table)}
if(!require(readr)){install.packages("readr"); library(readr)}
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

##################################################
## ADD BENTHIC DATA WITH IMPUTED MISSING VALUES ##
##################################################

benthic_PCA <- read_rds("outputs/RLS_benthic_PCA_imputed.rds")
benthic_PCA <- benthic_PCA %>%
  select(SurveyID, PC1_imputation, PC2_imputation,
         PC3_imputation, PC4_imputation) %>%
  dplyr::rename(PC1_imputed = PC1_imputation, PC2_imputed = PC2_imputation,
                PC3_imputed = PC3_imputation, PC4_imputed = PC4_imputation)

model_data <- merge(model_data, benthic_PCA, by = "SurveyID", all.x=T)

benthic_categories <- read_rds("outputs/RLS_benthic_data_imputed.rds") %>%
  select(-c(SiteCode, SiteLongitude, SiteLatitude))

model_data <- merge(model_data, benthic_categories, by = "SurveyID", all.x=T)

##########################################################
## IMPORT AND ATTACH COMMUNITY WEIGHTED MEAN AESTHETICS ##
##########################################################

cwm_aes <- read_rds("outputs/cwm_aes.rds") %>%
  rename(cwm_aes = aesthe_score)
com_sum_aes <- read_rds("outputs/com_sum_aes.rds") %>%
  rename(com_sum_aes = aesthe_score)

model_data <- merge(model_data, cwm_aes, by="SurveyID")
model_data <- merge(model_data, com_sum_aes, by="SurveyID")

##############################################
## IMPORT AND ATTACH TROPHIC STRUCTURE DATA ##
##############################################

trophic_composition <- readRDS("outputs/trophic_composition.rds")

model_data <- merge(model_data, trophic_composition, by="SurveyID")

##################################################
## IMPORT AND ATTACH TAXONOMIC COMPOSITION DATA ##
##################################################

taxo_stucture <- readRDS("outputs/taxo_structure.rds")

model_data <- merge(model_data, taxo_stucture, by="SiteCode", all.x=TRUE)

###################################################
## IMPORT AND ATTACH FAMILY AND ORDER ABUNDANCES ##
###################################################

family_abundances <- readRDS("outputs/family_abundances.rds")

order_abundances <- readRDS("outputs/order_abundances.rds")

model_data <- merge(model_data, family_abundances, by="SiteCode", all.x=TRUE)

model_data <- merge(model_data, order_abundances, by="SiteCode", all.x=TRUE)

#######################################
## CALCULATE AND ADD TOTAL ABUNDANCE ##
#######################################

sp_abund <- readRDS("outputs/sp_abund_matrix.rds")
sp_abund_survey <- sp_abund$SurveyID
sp_abund <- sp_abund %>%
  column_to_rownames("SurveyID") 
sp_abund <- data.frame(SurveyID=sp_abund_survey, abundance = rowSums(sp_abund))

model_data <- merge(model_data, sp_abund, by="SurveyID", all.x=TRUE)


#############################################
## IMPORT AND ATTACH PHYLOGENETIC AGE DATA ##
#############################################

phylo_age_occurence <- read_rds("outputs/cwm_age_occurence.rds")
phylo_age_abundance <- read_rds("outputs/cwm_age_abundance.rds")

model_data <- merge(model_data, phylo_age_occurence, by="SurveyID")
model_data <- merge(model_data, phylo_age_abundance, by="SurveyID")


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

# NB_SPECIES, TAXO_ENTROPY, PHYLO RICHNESS, GRAVITY, WAVE ENERGY, PHOSPHATE, NITRATE, NPP, 
# BIOMASS, PHYLOGENETIC AGE, ALL RIGHT-SKEWED, SHOULD BE LOG-TRANSFORMED
# DEPTH AND DHW COULD POTENTIALLY BE TRANSFORMED AS WELL

# Check minimum values
dplyr::summarise(model_data,
                 dplyr::across(c(nb_species, taxo_entropy, 
                                 abundance,
                                 phylo_richness,
                                 gravtot2, wave_energy, 
                                 BO_phosphate, BO_nitrate,
                                 NPP_mean, Biomass,
                                 Depth, dhw_mean,
                                 phylo_age_abundance, phylo_age_occurence),
                               min, na.rm=TRUE))
# Transform data
model_data <- dplyr::mutate(model_data,
                            dplyr::across(c(nb_species, taxo_entropy, abundance,
                                            gravtot2, wave_energy, 
                                            BO_phosphate, BO_nitrate,
                                            NPP_mean, Biomass,
                                            phylo_age_abundance, phylo_age_occurence),
                                          log))

# CREATE ABSOLUTE VALUE OF LATITUDE
model_data <- dplyr::mutate(model_data,
                            abs_latitude = abs(SiteLatitude))


# ADD A DUPLICATE OF LATITUDE AND LONGITUDE 
# TO CREATE A SCALED AND ORIGINAL VERSION

model_data$Latitude_scaled <- model_data$SiteLatitude
model_data$Longitude_scaled <- model_data$SiteLongitude

#####################################################
## SCALE ALL THE NUMERIC PREDICTORS TO MEAN 0 SD 2 ##
#####################################################

z_score_2sd <- function(x){ (x - mean(x,na.rm=T)) / (2*sd(x,na.rm=T))}

z_vars <- model_data %>%
  select_if(is.numeric) %>%
  select(-any_of(c("SurveyID","SiteLongitude","SiteLatitude",
                   "aesthe_survey_pres", "aesthe_survey_abund"))) %>%
  colnames()

standardized_data <- model_data %>%
  mutate_if(colnames(model_data) %in% z_vars, z_score_2sd)

#############################################
## REPLACE TEMPERATE BENTHIC 'DATA' WITH 0 ##
#############################################

standardized_data$PC1_imputed[standardized_data$Temperature_Zone=="Temperate"] <- 0
standardized_data$PC2_imputed[standardized_data$Temperature_Zone=="Temperate"] <- 0
standardized_data$PC3_imputed[standardized_data$Temperature_Zone=="Temperate"] <- 0
standardized_data$PC4_imputed[standardized_data$Temperature_Zone=="Temperate"] <- 0

standardized_data$coral_imputation[standardized_data$Temperature_Zone=="Temperate"] <- 0
standardized_data$algae_imputation[standardized_data$Temperature_Zone=="Temperate"] <- 0
standardized_data$CCA_imputation[standardized_data$Temperature_Zone=="Temperate"] <- 0

#################################
## SAVE SURVEY ID AS CHARACTER ##
#################################

standardized_data$SurveyID <- as.character(standardized_data$SurveyID)

##################################################################
## ERROR FOUND - MULTIPLE SITE NAMES WITH IDENTICAL COORDINATES ##
## THESE SHOULD BE MERGED AT THE SITE LEVEL 
##################################################################

# Step 3: Merge duplicates under the first SiteCode for each coordinate pair
standardized_data <- standardized_data %>%
  group_by(SiteLongitude, SiteLatitude) %>%
  mutate(SiteCode = first(SiteCode)) %>%
  ungroup()

###################
## SAVE THE DATA ##
###################

saveRDS(standardized_data, "outputs/standardized_data.rds")
standardized_data <-  read_rds("outputs/standardized_data.rds")

