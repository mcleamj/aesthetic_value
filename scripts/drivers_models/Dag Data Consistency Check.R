
#######################################################################################
#'  SCRIPT TO RUN DAG DATA CONSISTENCY CHECKS FOR THE DIRECTED ACYLIC GRAPH (DAG)
#'  THE DAG CAN BE FOUND AT:
#'  
#'  dagitty.net/mxttk25
#'
#' @author Matthew McLean, \email {mcleamj@gmail.com},
#'         
#' @date JULY 07, 2023
########################################################################################

######################
## LIBRARY PACKAGES ##
######################

if(!require(dagitty)){install.packages("dagitty"); library(dagitty)}
if(!require(base64enc)){install.packages("base64enc"); library(base64enc)}
if(!require(dplyr)){install.packages("dplyr"); library(dplyr)}
if(!require(readr)){install.packages("readr"); library(readr)}

####################
# DOWNLOAD THE DAG #
####################

#DAG <- downloadGraph("dagitty.net/mxttk25") # ORIGINAL, OUTDATED DAG
DAG <- downloadGraph("dagitty.net/m6WdXviAT") # UPDATED FEB 2024

names(DAG)

##############################
## IMPORT DATA "MODEL DATA" ##
##############################

model_data <- readRDS("data/model_data.rds")

model_data$Temperature_Zone <- as.factor(model_data$Temperature_Zone)

names(model_data)

model_data$MPA <- as.factor(model_data$MPA)

######################################
## IMPORT AND ADD DIVERSITY METRICS ##
######################################

biodiv <- read_rds("outputs/survey_biodiversity.rds")

model_data <- merge(model_data, biodiv, by="SurveyID")

##############################################
## IMPORT AND ADD TROPHIC STRUCTURE METRICS ##
##############################################

trophic <- read_rds("outputs/trophic_structure.rds")

model_data <- merge(model_data, trophic, by="SurveyID")

##################################################
## IMPORT AND ADD TAXONOMIC COMPOSITION METRICS ##
##################################################

taxo_structure <- read_rds("outputs/taxo_structure.rds")
taxo_structure$Taxo_PC1 <- taxo_structure$Taxo_PC1*-1 #ROTATE FOR POS RELATIONSHIP

model_data <- merge(model_data, taxo_structure, by="SiteCode")

##################################
## ADD NEW IMPUTED BENTHIC DATA ##
##################################

benthic_PCA <- readRDS("outputs/RLS_benthic_PCA_imputed.rds")
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

###############################################
## DO ANY VARIABLES NEED TO BE TRANSFORMED ? ##
###############################################

model_data <- do.call(data.frame,lapply(model_data, function(x) replace(x, is.infinite(x),NA)))

num_vars <- select_if(model_data, is.numeric)
num_vars$SurveyID <- NULL

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

###################################################
# MAKE SURE NAMES MATCH BETWEEN VARIABLES AND DAG #
###################################################

names(DAG)

names(model_data)[which(colnames(model_data)=="PC1_imputed")] <- "Benthic Composition"
names(model_data)[which(colnames(model_data)=="dhw_mean")] <- "Degree Heating Weeks"
names(model_data)[which(colnames(model_data)=="fshD")] <- "Fisheries Dependency"
names(model_data)[which(colnames(model_data)=="gravtot2")] <- "Human Gravity"
names(model_data)[which(colnames(model_data)=="MPA")] <- "Marine Protected Area"
names(model_data)[which(colnames(model_data)=="NPP_mean")] <- "Net Primary Productivity"
names(model_data)[which(colnames(model_data)=="Biomass")] <- "Reef Fish Biomass"
names(model_data)[which(colnames(model_data)=="sst_mean")] <- "Sea Surface Temperature"
names(model_data)[which(colnames(model_data)=="BO_parmean")] <- "Solar Irradiance"
names(model_data)[which(colnames(model_data)=="wave_energy")] <- "Wave Energy"
names(model_data)[which(colnames(model_data)=="fun_richness")] <- "Functional/Phylogenetic Diversity"
#names(model_data)[which(colnames(model_data)=="phylo_richness")] <- "Functional/Phylogenetic Diversity"
names(model_data)[which(colnames(model_data)=="BO_nitrate")] <- "Nitrate"
names(model_data)[which(colnames(model_data)=="BO_phosphate")] <- "Phosphate"
names(model_data)[which(colnames(model_data)=="taxo_richness")] <- "taxonomic_diversity"
names(model_data)[which(colnames(model_data)=="PC1_trophic")] <- "trophic_structure"
names(model_data)[which(colnames(model_data)=="HDI2017")] <- "Human Development Index"
names(model_data)[which(colnames(model_data)=="aesthe_survey")] <- "aesthetic_value"
names(model_data)[which(colnames(model_data)=="abs_latitude")] <- "Latitude"
names(model_data)[which(colnames(model_data)=="Taxo_PC1")] <- "taxonomic_composition"

#####################################
## SUBSET TO TROPICS OR TEMPERATE? ##
#####################################

# model_data <- model_data %>%
#   filter(Temperature_Zone=="Tropical")

######################
# MATCH DATA TO DAG #
######################

model_data <- model_data[,colnames(model_data) %in% names(DAG)]
str(model_data)

##############
# REMOVE NAs #
##############

model_data <- na.omit(model_data)

###########################################
# WHAT VARIABLES IN DAG BUT NOT IN DATA ? #
###########################################

names(DAG)[!names(DAG) %in% names(model_data)]

#############################################
# THE TEST WON'T WORK WITH MISSING VARIABLES
# SO CREATE A FAKE UPWELLING VARIABLE 
#############################################

model_data$Upwelling <- rnorm(n=nrow(model_data),0,1)

names(DAG)[!names(DAG) %in% names(model_data)]

###############################
##  STRUCTURE OF VARIABLES ? ##
###############################

str(model_data)
model_data$`Marine Protected Area` <- as.integer(model_data$`Marine Protected Area`)

#####################################################
# evaluate the d-separation implications of the DAG #
#####################################################

test <- localTests(DAG,model_data) 

test2 <- data.frame(test)

#####################################################
# SUBSET DATA BASED ON CORRELATION VALUE OR P VALUE #
#####################################################

testf <- subset(test2, estimate >= 0.3 | estimate <=-0.3)

#########################################
# show final independencies that failed #
#########################################

testf

nrow(testf)/nrow(test)

