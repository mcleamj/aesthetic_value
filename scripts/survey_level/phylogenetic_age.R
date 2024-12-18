

###############################################################
## CODE TO CALUCLATE AVERAGE PHYLOGENETIC AGE OF ASSEMBLAGES ##
###############################################################

library(FD)
library(fishtree)
library(tidyverse)

#########################################################################################
# Load data

## UPDATED NOVEMBER 2024 TO ANALYZE WITH ABUNDANCE INSTEAD OF BIOMASS
## PER PEER REVIEW  

survey_sp_abund <- readr::read_rds("outputs/sp_abund_matrix.rds")
survey_sp_abund <- survey_sp_abund %>%
  column_to_rownames("SurveyID") 

survey_sp_occ <- readr::read_rds("outputs/sp_pres_matrix.rds")
survey_sp_occ <- survey_sp_occ %>%
  column_to_rownames("SurveyID") 

sum(is.na(survey_sp_abund))
sum(is.na(survey_sp_occ))

# LOG TRANSFORM ABUNDANCE DATA #
survey_sp_abund <- log10(survey_sp_abund + 1)

multi_tree <- fishtree_complete_phylogeny(
  species = colnames(survey_sp_abund)
  #mc.cores = getOption("mc.cores", 1L)
)

tree <- multi_tree[[1]]

phylo_occ <- survey_sp_occ %>%
  as.data.frame() %>%
  select_if(colnames(survey_sp_occ) %in% tree$tip.label) %>%
  as.matrix()

phylo_abundance <- survey_sp_abund %>%
  select_if(colnames(survey_sp_abund) %in% tree$tip.label) %>%
  as.matrix() 

source("R/get_ages.R")
phylo_ages <- as.data.frame(get_ages(tree)) %>%
  rename(phylo_age = "get_ages(tree)") %>%
  rownames_to_column("species") %>%
  arrange(species)

identical(phylo_ages$species, colnames(phylo_occ))

###########################################
## AVVERAGE AGE - COMMUNITY WEIGHTED MEAN 
## OF OCCURENCE AND ABUNDANCE 
###########################################

phylo_ages <- phylo_ages %>%
  column_to_rownames("species")

cwm_age_occurence <- functcomp(as.matrix(phylo_ages), phylo_occ) %>%
  rownames_to_column("SurveyID") %>%
  rename(phylo_age_occurence = phylo_age)

cwm_age_abundance <- functcomp(as.matrix(phylo_ages), phylo_abundance) %>%
  rownames_to_column("SurveyID") %>%
  rename(phylo_age_abundance = phylo_age)

saveRDS(cwm_age_occurence, "outputs/cwm_age_occurence.rds")
saveRDS(cwm_age_abundance, "outputs/cwm_age_abundance.rds")

