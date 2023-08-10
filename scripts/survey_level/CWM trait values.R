
###################################################################################################
#' Script to calculate community weighted mean trait values
#'
#' @author Matthew McLean, \email {mcleamj@@gmail.com},
#'         Nicolas Mouquet , \email{nicolas.mouquet@@cnrs.fr},
#'         
#'    
#' @date 2021/08/11
##################################################################################################

library(stringr)
library(tidyverse)
library(DataCombine)
library(FD)
library(factoextra)
library(plot3D)

rm(list = ls())

#########################################################################################
# Load data

site_info <- read_rds("data/RLS_sitesInfos_tropical.rds")

sp_ref <- read.table("data/sp_name_ref.txt")
species_traits <- read.csv("data/species_traits.csv")
survey_sp_biom <- readr::read_rds("outputs/sp_biomass_matrix.rds")
survey_sp_occ <- readr::read_rds("outputs/sp_pres_matrix.rds")

survey_sp_biom <- survey_sp_biom %>%
  column_to_rownames("SurveyID") 

survey_sp_occ <- survey_sp_occ %>%
  column_to_rownames("SurveyID") 

species_traits <- species_traits %>%
  mutate(Species = str_replace(Species, " ", "_")) 

species_traits <- FindReplace(species_traits, "Species", sp_ref, from="SPECIES_NAME", to="species",
                       exact = T, vector = F)

species_traits <- species_traits %>%
  distinct(Species, .keep_all = TRUE) %>%
  filter(Species %in% colnames(survey_sp_biom))

survey_sp_biom <- survey_sp_biom %>%
  select_if(colnames(survey_sp_biom) %in% species_traits$Species)

survey_sp_occ <- survey_sp_occ %>%
  select_if(colnames(survey_sp_occ) %in% species_traits$Species) %>%
  as.matrix()

identical(sort(unique(species_traits$Species)), sort(unique(colnames(survey_sp_biom))))

sum(is.na(survey_sp_biom))
sum(is.na(survey_sp_occ))

# SURVEYS WITH NO OBSERVATIONS?
survey_sp_biom <- survey_sp_biom[-which(rowSums(survey_sp_biom)==0),]
survey_sp_occ <- survey_sp_occ[-which(rowSums(survey_sp_occ)==0),]

# LOG TRANSFORM BIOMASS DATA #
survey_sp_biom <- log10(survey_sp_biom + 1)

# CWM TRAITS USING OCCURENENCE DATA
species_traits <- species_traits %>%
  arrange(Species) %>%
  column_to_rownames("Species") %>%
  select(Diet, Trophic.Level, BodySize)
str(species_traits)
hist(species_traits$Trophic.Level)
hist(species_traits$BodySize)
species_traits <- species_traits %>%
  mutate(BodySize = log(BodySize)) %>%
  mutate(Diet = as.factor(Diet))

sapply(species_traits, function(x) paste(round(sum(is.na(x))/length(x),2)*100,"%",sep=""))
str(species_traits)

cwm_traits_cont <- functcomp(as.matrix(select(species_traits, c(BodySize, Trophic.Level))), as.matrix(survey_sp_occ))
cwm_traits_cat <-  functcomp(as.matrix(select(species_traits, Diet)), as.matrix(survey_sp_occ),CWM.type = "all")

cwm_traits <- data.frame(cwm_traits_cont, cwm_traits_cat)                       


# PCA FOR TROPHIC GUILD STRUCTURE

trophic_PCA <- prcomp(cwm_traits_cat, scale. = FALSE)

fviz_pca_var(trophic_PCA, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)
## GROUP DATA

trophic_structure <- data.frame(cwm_traits,
                            trophic_PCA$x[,1:4]) %>%
  dplyr::rename(PC1_trophic = PC1, PC2_trophic = PC2,
         PC3_trophic = PC3, PC4_trophic = PC4)
trophic_structure$SurveyID <- rownames(trophic_structure)

saveRDS(trophic_structure, "outputs/trophic_structure.rds")

