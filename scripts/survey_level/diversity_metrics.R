###################################################################################################
#'  SCRIPT TO COMPUTE TAXO/FUN/PHYLO DIVERSITY METRICS AT THE SURVEY LEVEL
#'  
#' @author Mattia Ghilardi \email {mattia.ghilardi@@leibniz-zmt.de}
#'         Matthew McLean, \email {mcleamj@@gmail.com}
#'         
#' @date updated May 2024

##################################################################################################

#########################################################################################
#### Clean environment
rm(list = ls())

#########################################################################################
#### Load packages

library(tidyverse)
#library(readr)
#library(dplyr)
#library(tibble)
library(mFD)
library(FD)
library(FishPhyloMaker)
library(hillR)
library(funrar)
library(fishtree)
library(geiger)

#########################################################################################
# Load data

survey_sp_biom <- readr::read_rds("outputs/sp_biomass_matrix.rds")
survey_sp_biom <- survey_sp_biom %>%
  column_to_rownames("SurveyID") 

survey_sp_occ <- readr::read_rds("outputs/sp_pres_matrix.rds")
survey_sp_occ <- survey_sp_occ %>%
  column_to_rownames("SurveyID") 

sp_traits <- read.csv("Data/All.RLS.traits.06222022.csv")

# Select only traits of interest
sp_traits <- sp_traits %>%
  select(Species, Gregariousness, Water.Column, Diet, Habitat, TrophicLevel, BodySize)

sp_traits <- sp_traits %>%
  mutate(Species = str_replace(Species, " ", "_")) %>%
  filter(Species %in% colnames(survey_sp_biom))

survey_sp_biom <- survey_sp_biom %>%
  select_if(colnames(survey_sp_biom) %in% sp_traits$Species)

survey_sp_occ <- survey_sp_occ %>%
  select_if(colnames(survey_sp_occ) %in% sp_traits$Species) %>%
  as.matrix()

identical(sort(unique(sp_traits$Species)), sort(unique(colnames(survey_sp_biom))))

sum(is.na(survey_sp_biom))
sum(is.na(survey_sp_occ))

# SURVEYS WITH NO OBSERVATIONS?
survey_sp_biom <- survey_sp_biom[-which(rowSums(survey_sp_biom)==0),]
survey_sp_occ <- survey_sp_occ[-which(rowSums(survey_sp_occ)==0),]

# LOG TRANSFORM BIOMASS DATA #
survey_sp_biom <- log10(survey_sp_biom + 1)

# Relative biomass
survey_sp_biom_rel <- survey_sp_biom %>%
  as.matrix() %>%
  make_relative()

table(rowSums(survey_sp_biom_rel))

sum(is.na(survey_sp_biom_rel))

## Taxonomic diversity

# Species richness and Shannon-like entropy (eq number of species, exp(H'))
survey_biodiversity <- tibble::tibble(
  SurveyID = rownames(survey_sp_biom),
  taxo_richness = apply(survey_sp_occ, 1,  sum),
  taxo_entropy = apply(survey_sp_biom_rel, 1,
                       function(x) {exp((-1) * sum(x[x != 0] * log(x[x != 0])))})
  )

## Functional diversity

# Rearrange trait data
sp_traits <- sp_traits %>%
  column_to_rownames("Species") %>%
  mutate(across(c(Diet, Habitat), as.factor)) %>%
  mutate(across(c(Gregariousness, Water.Column, TrophicLevel, BodySize), ordered))
str(sp_traits)
  

# Type of traits for mFD::funct.dist()
traits_cat <- tibble::tibble(
  trait_name = c("Gregariousness", "Water.Column", "Diet","Habitat","TrophicLevel","BodySize"),
  trait_type = c("O", "O", "N", "N", "O", "O")
  )

# Gower distance
sp_dist <- mFD::funct.dist(sp_tr = sp_traits,
                           tr_cat = traits_cat,
                           metric = "gower")
# there is a warning due to species having same trait values
# not an issue for computation of Chao et al 2019 indices

# Compute functional richness and functional entropy based on Chao 2019

survey_biodiversity <- survey_biodiversity %>%
  # Richness on species occurrences with q=0
  bind_cols(fun_richness = mFD::alpha.fd.hill(asb_sp_w = survey_sp_occ,
                                              sp_dist = sp_dist,
                                              q = 0,
                                              tau = "mean",
                                              details_returned = FALSE)[, 1],
            # Shannon-like entropy on species relative biomass with q=1
            fun_entropy = mFD::alpha.fd.hill(asb_sp_w = survey_sp_biom_rel,
                                             sp_dist = sp_dist,
                                             q = 1,
                                             tau = "mean",
                                             details_returned = FALSE)[,1]
            )

## Phylogenetic diversity

# Following Chao et al. 2010

multi_tree <- fishtree_complete_phylogeny(
  species = colnames(survey_sp_biom),
  mc.cores = getOption("mc.cores", 1L)
)

tree <- multi_tree[[1]]

phylo_occ <- survey_sp_occ %>%
  as.data.frame() %>%
  select_if(colnames(survey_sp_occ) %in% tree$tip.label) %>%
  as.matrix()

phylo_biomass <- survey_sp_biom %>%
  select_if(colnames(survey_sp_biom) %in% tree$tip.label) %>%
  as.matrix() %>%
  make_relative()

survey_biodiversity <- survey_biodiversity %>%
  # Richness on species occurrences with q=0
  bind_cols(phylo_richness = hillR::hill_phylo(comm = phylo_occ,
                                               tree = tree,
                                               q = 0),
            # Shannon-like entropy on species relative biomass with q=1
            phylo_entropy = hillR::hill_phylo(comm = phylo_biomass,
                                              tree = tree,
                                              q = 1)
  )

# Save results
readr::write_rds(survey_biodiversity,
                 here::here("outputs/survey_biodiversity.rds"),
                 compress = "gz")
