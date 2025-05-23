
###################################################################################################
#' Script to calculate community weighted mean trait values
#'
#' @author Matthew McLean, \email {mcleamj@@gmail.com},
#'         Nicolas Mouquet , \email{nicolas.mouquet@@cnrs.fr},
#'         
#'    
#' @date updated May 2024
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
#survey_sp_biom <- readr::read_rds("outputs/sp_biomass_matrix.rds")
survey_sp_abund <- readr::read_rds("outputs/sp_abund_matrix.rds")
survey_sp_occ <- readr::read_rds("outputs/sp_pres_matrix.rds")

survey_sp_abund <- survey_sp_abund %>%
  column_to_rownames("SurveyID") 

survey_sp_occ <- survey_sp_occ %>%
  column_to_rownames("SurveyID") 

species_traits <- species_traits %>%
  mutate(Species = str_replace(Species, " ", "_")) 

species_traits <- FindReplace(species_traits, "Species", sp_ref, from="SPECIES_NAME", to="species",
                       exact = T, vector = F)

species_traits <- species_traits %>%
  distinct(Species, .keep_all = TRUE) %>%
  filter(Species %in% colnames(survey_sp_abund))

survey_sp_abund <- survey_sp_abund %>%
  select_if(colnames(survey_sp_abund) %in% species_traits$Species)

survey_sp_occ <- survey_sp_occ %>%
  select_if(colnames(survey_sp_occ) %in% species_traits$Species) %>%
  as.matrix()

identical(sort(unique(species_traits$Species)), sort(unique(colnames(survey_sp_abund))))

sum(is.na(survey_sp_abund))
sum(is.na(survey_sp_occ))

# SURVEYS WITH NO OBSERVATIONS?
survey_sp_abund <- survey_sp_abund[-which(rowSums(survey_sp_abund)==0),]
survey_sp_occ <- survey_sp_occ[-which(rowSums(survey_sp_occ)==0),]

# LOG TRANSFORM ABUNDANCE DATA #
survey_log_abund <- log10(survey_sp_abund + 1)

# CWM TRAITS USING ABUNDANCE DATA
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

##############################################################
# SPECIFY WHETHER TO USE OCCURENCE OR ABUNDANCE DATA BELOW
# UPDATED IN THE MANUSCRIPT - ALL METRICS NOW WITH ABUNDANCE
##############################################################

cwm_traits_cont <- functcomp(as.matrix(select(species_traits, c(BodySize, Trophic.Level))), as.matrix(survey_log_abund))
cwm_traits_cat <-  functcomp(as.matrix(select(species_traits, Diet)), as.matrix(survey_log_abund),CWM.type = "all")

# NEED TO TRANSFORM THE CATEGORICAL TRAITS
cwm_traits_cat <- asin(sqrt(cwm_traits_cat))

cwm_traits <- data.frame(cwm_traits_cont, cwm_traits_cat)                       

# PCA FOR TROPHIC GUILD COMPOSITION

trophic_PCA <- prcomp(cwm_traits_cat, scale. = FALSE)

fviz_pca_var(trophic_PCA, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = FALSE)

# ROTATE THE ORDINATION AXES BASED ON EXPLORATORY ANALYSIS
# THIS WILL PROVIDE A POSITIVE RELATIONSHIP BETWEEN 
# PC AXES AND AESTHETIC VALUE
# WHICH WILL MAKE THE RESULTS MORE INTUITIVE AND CLEAR
# BUT WILL HAVE NO EFFECT ON THE OVERALL RESULT 
# OR THE ORDINATION ITSELF

# FLIP PC1 AND PC2

trophic_PCA$x[,1] <- trophic_PCA$x[,1]*-1
trophic_PCA$rotation[,1] <- trophic_PCA$rotation[,1] * -1

trophic_PCA$x[,2] <- trophic_PCA$x[,2]*-1
trophic_PCA$rotation[,2] <- trophic_PCA$rotation[,2] * -1

# RE PLOT PCA
fviz_pca_var(trophic_PCA, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = FALSE)


## GROUP DATA
trophic_composition <- data.frame(cwm_traits,
                            trophic_PCA$x[,1:4]) %>%
  dplyr::rename(PC1_trophic = PC1, PC2_trophic = PC2,
         PC3_trophic = PC3, PC4_trophic = PC4)
trophic_composition$SurveyID <- rownames(trophic_composition)

saveRDS(trophic_composition, "outputs/trophic_composition.rds")


############################################
## ALTERNATIVE 
## CALCUALTE TOTAL ABUNDANCE OF EACH GROUP
## RATHER THAN CWM PROPORTION
############################################

trophic_abundance <- survey_sp_abund %>%
  rownames_to_column("SurveyID") %>%
  pivot_longer(cols = -SurveyID, names_to = "Species", values_to = "abundance") %>%
  filter(abundance != 0) %>%
  left_join(species_traits, by="Species") %>%
  select(SurveyID, abundance, Diet) %>%
  group_by(SurveyID, Diet) %>%
  summarise_all(.funs=sum) %>%
  pivot_wider(id_cols = SurveyID, names_from = Diet, values_from = abundance ) %>%
  mutate_all(~ replace(., is.na(.), 0))

log_trophic_abund <- trophic_abundance %>%
  mutate_all(~ log10(. + 1)) %>%
  column_to_rownames("SurveyID")

trophic_abund_PCA <- prcomp(log_trophic_abund, scale. = FALSE)

fviz_pca_var(trophic_abund_PCA, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = FALSE)

