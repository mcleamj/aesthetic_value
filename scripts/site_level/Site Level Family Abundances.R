
#######################################################################################
#'  CODE TO CALCULATE FAMILY ABUNDANCES FOR EACH SITE
#'
#' @author Matthew McLean, \email {mcleamj@gmail.com},
#'         
#' @date Updated May, 2023
########################################################################################

######################
## LIBRARY PACKAGES ##
######################

#install.packages("data/worms", repos = NULL, type="source")
library(worms)

#if(!require(data.table)){install.packages("data.table"); library(data.table)}
#if(!require(mapproj)){install.packages("mapproj"); library(mapproj)}
#if(!require(plot3D)){install.packages("plot3D"); library(plot3D)}
#if(!require(pals)){install.packages("pals"); library(pals)}
#if(!require(brms)){install.packages("brms"); library(brms)}
#if(!require(FD)){install.packages("FD"); library(FD)}
#if(!require(bayesplot)){install.packages("bayesplot"); library(bayesplot)}
#if(!require(parallel)){install.packages("parallel"); library(parallel)}
#if(!require(rstan)){install.packages("rstan"); library(rstan)}
#if(!require(bayestestR)){install.packages("bayestestR"); library(bayestestR)}
if(!require(readr)){install.packages("readr"); library(readr)}
if(!require(tibble)){install.packages("tibble"); library(tibble)}
if(!require(dplyr)){install.packages("dplyr"); library(dplyr)}

##############################
## IMPORT DATA "MODEL DATA" ##
##############################

model_data <- readRDS("data/model_data.rds")

names(model_data)

model_data$MPA <- as.factor(model_data$MPA)

#####################################
## IMPORT AND ATTACH AESTHETIC DATA #
#####################################

aesthe_survey_data <- read.csv("outputs/survey_aesth.csv")

model_data <- merge(model_data, aesthe_survey_data, by="SurveyID")

##'##################################
##' IMPORT ABUNDANCE INFORMATION 
##' AGGREGATE OCCURENCES TO SITE LEVEL
##' AND CALCULATE FAMILY ABUNDANCES 
##' FOR EACH SITE
##' #################################

survey_sp_abund <- readr::read_rds("outputs/sp_abund_matrix.rds")

survey_species <- survey_sp_abund %>%
  dplyr::select(where(is.numeric)) %>%
  colnames() %>%
  as.data.frame() %>%
  dplyr::rename("species" = ".")

survey_species$species <- gsub("_", " ", survey_species$species) 
survey_species <- survey_species %>%
  dplyr::rename("species_name" = "species")

survey_species$scientificname <- gsub(" spp.", "", survey_species$species_name)

family_info <- worms::wormsbynames(survey_species$scientificname) 
family_info <- merge(family_info, survey_species, by="scientificname")
saveRDS(family_info, "outputs/family_info.rds")
family_info <- read_rds("outputs/family_info.rds")

##############################################
## CALCULATE FAMILY ABUNDANCES AT EACH SITE ##
##############################################

site_sp_abund <- merge(model_data[,c("SiteCode","SurveyID")], survey_sp_abund,
                     by="SurveyID")

site_sp_abund <- site_sp_abund %>% # AVERAGE ABUNDACE OF EACH SPECIES PER SITE (AVG OF 1 OR 2 TRANSECTS)
  select(-SurveyID) %>%
  group_by(SiteCode) %>%
  summarise_all(.funs = mean, na.rm=T) 

family_trait <- family_info %>%
  select(species_name, family) %>%
  arrange(species_name) 

family_abundances <- site_sp_abund %>%
  pivot_longer(!SiteCode, names_to = "species_name", values_to = "abundance") %>%
  mutate("species_name" = gsub("_", " ", species_name )) %>%
  left_join(family_trait, by="species_name") %>%
  select(-species_name) %>%
  group_by(SiteCode, family) %>%
  summarise_all(.funs=sum) %>% # SUM OF INDIVIDUALS PER FAMILY (SUM OF THE AVERAGE ABUNDANCES PER SITE)
  pivot_wider(names_from="family", values_from = "abundance") 

saveRDS(family_abundances, "outputs/family_abundances.rds")

