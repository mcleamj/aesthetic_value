###################################################################################################
#'  Information at the survey scale
#'
#' @author Matthew McLean, \email {mcleamj@@gmail.com},
#'         Juliette Langlois, \email{juliette.a.langlois@@gmail.com},
#'         Nicolas Mouquet , \email{nicolas.mouquet@@cnrs.fr},
#'         
#'         
#'
#' @date 2021/08/11
##################################################################################################

library(tidyverse)
library(DataCombine)

rm(list = ls())

## Import list of species with aesthethic values
#esth_sp <- read.csv(here::here("data", "esthe_table.csv"))
# esth_sp <- read.csv(here::here("data", "aesthe_langlois_2022.csv"),
#                     sep=";")
esth_sp <- read.csv(here::here("outputs", "aesthe_species.csv"),
                    sep=";")


## Any duplicated species?
esth_sp$sp_name[which(duplicated(esth_sp$sp_name))]

# Remove duplicates

## Import cleaned RLS data - errors removed, species names corrected
survey_spcompo <- readRDS(here::here("data", "RLS_fishdata_cleaned.rds"))
length(unique(survey_spcompo$SurveyID)) # 7017 SURVYES
length(unique(survey_spcompo$SiteCode)) # 3529 SITES

# Replace spaces with underscores
names(survey_spcompo)
length(unique(survey_spcompo$species))
survey_spcompo$SurveyID <- as.character(survey_spcompo$SurveyID)
survey_spcompo$species <- gsub(" ", "_", survey_spcompo$species)
survey_spcompo$TAXONOMIC_NAME <- gsub(" ", "_", survey_spcompo$TAXONOMIC_NAME)
survey_spcompo$SPECIES_NAME <- gsub(" ", "_", survey_spcompo$SPECIES_NAME)

esth_sp$sp_name <- gsub(" ", "_", esth_sp$sp_name)

## How many species in each data set?
length(unique(survey_spcompo$species))  #2655 in RLS
length(unique(esth_sp$sp_name)) #2417 in aesthetic scores

## How many species in common between the two data sets?
length(which(esth_sp$sp_name %in% survey_spcompo$species)) #2254 species in common

# How many species have an aesthetic value score but not observed in RLS?
# Which species?
missing_esth_sp <- data.frame(missing=esth_sp$sp_name[which(!esth_sp$sp_name %in% survey_spcompo$species)])
length(missing_esth_sp$missing) #163 species

# How many species observed in RLS but don't have an aesthetic value score?
# Which species?
missing_survey_sp <- data.frame(missing=survey_spcompo$species[which(!survey_spcompo$species %in% esth_sp$sp_name)])
missing_survey_sp <- missing_survey_sp %>%
  distinct()
length(missing_survey_sp$missing) # 401 species

# Which are actual species and which are Genus spp?
missing_survey_genera <- missing_survey_sp %>%
  filter(grepl("spp",missing))
length(missing_survey_genera$missing) #164 Genus spp.

missing_survey_sp <- missing_survey_sp %>%
  filter(!grepl("spp",missing))
length(missing_survey_sp$missing) #237 actual species don't have an aesthetic score

## Are species with aesthetic scores but not found in RLS found in the old (outdated) RLS species names?
length(which(missing_esth_sp$missing %in% survey_spcompo$SPECIES_NAME)) # at least 18 of them are

# Change the outdated esth_sp name to the updated RLS name
sp_ref <- survey_spcompo %>%
  select(SPECIES_NAME, species) %>%
  distinct()

esth_sp <- FindReplace(esth_sp, "sp_name", sp_ref, from="SPECIES_NAME", to="species",
            exact = T, vector = F)

length(unique(esth_sp$sp_name)) - length(which(esth_sp$sp_name %in% survey_spcompo$species)) # 20 (163-143) have been replaced

#20 NAMES ARE OUTDATED IN THE AESTHETIC VALUE DATA

## Subset the two data sets to only matching species
survey_spcompo <- survey_spcompo[which(survey_spcompo$species %in% esth_sp$sp_name),]
esth_sp <- esth_sp[which(esth_sp$sp_name %in% survey_spcompo$species),]

length(unique(esth_sp$sp_name)) # How many aesthethic value species left?
length(unique(survey_spcompo$species)) # How many RLS species left?

identical(sort(unique(survey_spcompo$species)), sort(unique(esth_sp$sp_name))) # Are they identical?


## community matrix by species ##

#---- Biomass Matrix -----#

sp_biomass_matrix <- survey_spcompo %>%
  group_by(species, SurveyID) %>%
  summarise(sum(Biomass)) %>%
  rename(Biomass = 'sum(Biomass)')

sp_biomass_matrix <- spread(sp_biomass_matrix, species, value=Biomass)
sp_biomass_matrix[is.na(sp_biomass_matrix)] <- 0
saveRDS(sp_biomass_matrix, "outputs/sp_biomass_matrix.rds")

#---- Abundance Matrix -----#

sp_abund_matrix <- survey_spcompo %>%
  group_by(species, SurveyID) %>%
  summarise(sum(Num)) %>%
  rename(Abund = 'sum(Num)')

sp_abund_matrix <- spread(sp_abund_matrix, species, value=Abund)
sp_abund_matrix[is.na(sp_abund_matrix)] <- 0
saveRDS(sp_abund_matrix, "outputs/sp_abund_matrix.rds")

#---- Presence Matrix -----#

sp_pres_matrix <- sp_abund_matrix
sp_pres_matrix[,2:ncol(sp_abund_matrix)][sp_pres_matrix[,2:ncol(sp_abund_matrix)]>0] <- 1
saveRDS(sp_pres_matrix, "outputs/sp_pres_matrix.rds")


# ## community matrix by species & size class ##
# 
# #---- Biomass Matrix -----#
# 
# sizeclass_biomass_matrix <- survey_spcompo %>%
#   group_by(species, Sizeclass, SurveyID) %>%
#   summarise(sum(Biomass)) %>%
#   rename(Biomass = 'sum(Biomass)')
# 
# sizeclass_biomass_matrix$species_size <- paste(sizeclass_biomass_matrix$species, sizeclass_biomass_matrix$Sizeclass)
# 
# sizeclass_biomass_matrix <- spread(sizeclass_biomass_matrix, species_size, value=Biomass)
# sizeclass_biomass_matrix[is.na(sizeclass_biomass_matrix)] <- 0
# 
# #---- Abundance Matrix -----#
# 
# sizeclass_abund_matrix <- survey_spcompo %>%
#   group_by(species, SurveyID) %>%
#   summarise(sum(Num)) %>%
#   rename(Abund = 'sum(Num)')
# 
# sizeclass_abund_matrix <- spread(sizeclass_abund_matrix, species, value=Abund)
# sizeclass_abund_matrix[is.na(sizeclass_abund_matrix)] <- 0
# 
# #---- Presence Matrix -----#
# 
# sizeclass_pres_matrix <- sizeclass_abund_matrix
# sizeclass_pres_matrix[,2:ncol(sizeclass_pres_matrix)][sizeclass_pres_matrix[,2:ncol(sizeclass_pres_matrix)]>0] <- 1
