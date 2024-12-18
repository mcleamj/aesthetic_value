
#######################################################################################
#'  CODE TO CALCULATE ORDER ABUNDANCES FOR EACH SITE
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

if(!require(tidyverse)){install.packages("dplyr"); library(dplyr)}

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
##' AND CALCULATE order ABUNDANCES 
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

order_info <- worms::wormsbynames(survey_species$scientificname) 
order_info <- merge(order_info, survey_species, by="scientificname")
saveRDS(order_info, "outputs/order_info.rds")
order_info <- read_rds("outputs/order_info.rds")

##############################################
## CALCULATE order ABUNDANCES AT EACH SITE ##
##############################################

site_sp_abund <- merge(model_data[,c("SiteCode","SurveyID")], survey_sp_abund,
                       by="SurveyID")

site_sp_abund <- site_sp_abund %>% # AVERAGE ABUNDACE OF EACH SPECIES PER SITE (AVG OF 1 OR 2 TRANSECTS)
  select(-SurveyID) %>%
  group_by(SiteCode) %>%
  summarise_all(.funs = mean, na.rm=T) 

order_trait <- order_info %>%
  select(species_name, order) %>%
  arrange(species_name) 

order_abundances <- site_sp_abund %>%
  pivot_longer(!SiteCode, names_to = "species_name", values_to = "abundance") %>%
  mutate("species_name" = gsub("_", " ", species_name )) %>%
  left_join(order_trait, by="species_name") %>%
  select(-species_name) %>%
  group_by(SiteCode, order) %>%
  summarise_all(.funs=sum) %>% # SUM OF INDIVIDUALS PER order (SUM OF THE AVERAGE ABUNDANCES PER SITE)
  pivot_wider(names_from="order", values_from = "abundance") 

saveRDS(order_abundances, "outputs/order_abundances.rds")

