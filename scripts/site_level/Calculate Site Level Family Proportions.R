
#######################################################################################
#'  CODE TO CALCULATE FAMILY PROPORTIONS FOR EACH SITE
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

aaesthe_surveyetic_data <- read.csv("outputs/survey_aesth.csv")

model_data <- merge(model_data, aaesthe_surveyetic_data, by="SurveyID")

##'##################################
##' IMPORT OCCURENCE INFORMATION 
##' AGGREGATE OCCURENCES TO SITE LEVEL
##' AND CALCULATE FAMILY PROPORTIONS 
##' FOR EACH SITE
##' #################################

survey_sp_occ <- readr::read_rds("outputs/sp_pres_matrix.rds")

survey_species <- survey_sp_occ %>%
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

## CALCULATE FAMILY PROPORTIONS AT EACH SITE

site_sp_occ <- merge(model_data[,c("SiteCode","SurveyID")], survey_sp_occ,
                     by="SurveyID")
site_sp_occ <- site_sp_occ %>%
  select(-SurveyID) %>%
  group_by(SiteCode) %>%
  summarise_all(.funs = mean, na.rm=T) %>%
  mutate_if(is.numeric, ~1 * (. > 0))

family_trait <- family_info %>%
  select(species_name, family) %>%
  arrange(species_name) %>%
  column_to_rownames("species_name")

site_sp_occ <- site_sp_occ %>%
  column_to_rownames("SiteCode")
colnames(site_sp_occ) <- gsub("_", " ", colnames(site_sp_occ) )
site_sp_occ <- site_sp_occ %>%
  select(order(colnames(site_sp_occ)))

identical(colnames(site_sp_occ), rownames(family_trait))

family_proportions <- functcomp(as.matrix(family_trait),
                              as.matrix(site_sp_occ),
                              CWM.type = "all")

colnames(family_proportions) <- gsub("family_", "", colnames(family_proportions))

family_proportions <- family_proportions %>%
  rownames_to_column("SiteCode")

saveRDS(family_proportions, "outputs/family_proportions.rds")

