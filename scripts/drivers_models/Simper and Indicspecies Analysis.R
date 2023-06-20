


#######################################################################################
#'  CODE TO EXAME SPECIES AND FAMILY DIFFERENCES BETWEEN MPAS AND FISHED SITES
#'  
#'  dagitty.net/mxttk25
#'
#' @author Matthew McLean, \email {mcleamj@gmail.com},
#'         
#' @date JUNE 9, 2023
########################################################################################


######################
## LIBRARY PACKAGES ##
######################

if(!require(data.table)){install.packages("data.table"); library(data.table)}
if(!require(DataCombine)){install.packages("DataCombine"); library(DataCombine)}
if(!require(mapproj)){install.packages("mapproj"); library(mapproj)}
if(!require(plot3D)){install.packages("plot3D"); library(plot3D)}
if(!require(pals)){install.packages("pals"); library(pals)}
if(!require(ggplot2)){install.packages("ggplot2"); library(ggplot2)}
if(!require(viridis)){install.packages("viridis"); library(viridis)}
if(!require(brms)){install.packages("brms"); library(brms)}
if(!require(bayesplot)){install.packages("bayesplot"); library(bayesplot)}
if(!require(parallel)){install.packages("parallel"); library(parallel)}
if(!require(rstan)){install.packages("rstan"); library(rstan)}
if(!require(bayestestR)){install.packages("bayestestR"); library(bayestestR)}
if(!require(indicspecies)){install.packages("indicspecies"); library(indicspecies)}
if(!require(vegan)){install.packages("vegan"); library(vegan)}
if(!require(readr)){install.packages("readr"); library(readr)}
if(!require(tibble)){install.packages("tibble"); library(tibble)}
if(!require(dplyr)){install.packages("dplyr"); library(dplyr)}

library(worms)


#################################################################
## SIMPER ANALYSIS AND INDISPECIES ANALYSIS BY ECOREGION
## TAKE MEANS AND MAKE SOME KIND OF RANKED BOXPLOT BY FAMILY 
##################################################################

## AGGREGATE OCCURENCES BY FAMILY ? GENUS ? 
## SPECIES COUNTS PER FAMILY - GET A MATRIX WITH NUM OF SPECIES PER FAMILY
## OR JUST AN OCCURNECE PER FAMILY?

# GET FAMILY INFO FOR SPECIES LIST

# GET GLOBAL SIMPER AND INDISPECIES, THEN BY ECOREGION
# ONLY ECOREGIONS WITH BOTH FISHED AND NOT FISHED
# ONLY FOR ECOREGIONS WHERE MPA HAS HIGHER AESTHETIC VALUE OR RICHNESS?
# FOR NOW, ALL ECOREGIONS
# OR** - RUN IT ON THE SPECIES (OR GENUS) LEVEL FOR EACH ECOREGION - I TINK THIS IS THE WINNER
# SAVE THE CONTRIBUTION SCORE FOR EACH SPECIES
# THEN CALCULATE AVERAGE SCORE PER FAMILY AFTERWARDS

# IMPORT SPECIES OCCURENCE MATRIX
# COMBINE OCCURENCE MATRIX WITH SURVEY INFORMATION
# FILTER TO ONLY SELECTED ECOREGION
# SUM BY FAMILY PER SITE PER ECOREGION

survey_sp_occ <- readr::read_rds("outputs/sp_pres_matrix.rds")

survey_info <- model_data %>%
  select(SurveyID, SiteCode, Ecoregion, MPA)

survey_sp_occ <- merge(survey_info, survey_sp_occ, by="SurveyID")

# AGGREGATE TO SITE LEVEL AND GET OCCURENCE PER SITE FOR EACH SPECIES
# NOT HAPPY WITH THIS CHUNK OF CODE
survey_sp_occ <- survey_sp_occ %>% 
  filter(SiteCode %in% no_take_only$SiteCode) %>% # ONLY NO TAKE MPA, ONLY ECOREGIONS WITH FISHED AND MPA SITES
  select(-SurveyID) %>%
  group_by(SiteCode, Ecoregion, MPA) %>%
  summarise_all(.funs = mean) %>% # TAKE MEAN PER SPECIES BY SITE (MEAN OF 0S AND 1S)
  mutate_if(is.numeric, ~1 * (. > 0)) %>% # CONVERT ALL VALUES GREATER THAN 0 TO 1
  as.data.frame()

survey_species <- survey_sp_occ %>%
  ungroup() %>%
  dplyr::select(where(is.numeric)) %>%
  colnames() %>%
  as.data.frame() %>%
  rename("species"= ".")

survey_species$species <- gsub("_", " ", survey_species$species)
survey_species$scientificname <- gsub(" spp.", "", survey_species$species)

family_info <- worms::wormsbynames(survey_species$scientificname) 

#################################################
## CALCULATE PHYLOGENETIC AGES 
## USING CODE FROM LANGLOIS ET AL PLOS BIOLOGY   
#################################################

survey_species$species_underscore <- gsub(" ", "_",  survey_species$species)

# Download the phylogenetic tree for all fishes to compute the age
set100_all  <- fishtree::fishtree_complete_phylogeny(survey_species$species_underscore, mc.cores = 1)
# Dropping names not in list
set100 <- parallel::mclapply(set100_all,function(x){
  ape::drop.tip(x, x$tip.label[!is.element(x$tip.label,as.character(survey_species$species_underscore))])
}, mc.cores = 1)

# Compute the age of all species found by fishtree (for each of the 100 trees)
ages      <- do.call(cbind, lapply(set100, get_ages)) 
ages_mean <- apply(ages, 1, mean, na.rm = T) # mean across the trees
ages_mean <- ages_mean[match(survey_species$species_underscore, names(ages_mean))]
species_ages <- data.frame(scientificname=survey_species$scientificname, ages_mean)
family_ages <- family_info %>%
  select(scientificname, family)
family_ages <- merge(family_ages, species_ages, by="scientificname") %>%
  select(family, ages_mean) %>%
  group_by(family) %>%
  summarise_all(.funs=mean, na.rm=T)
family_ages$log_age <- log(family_ages$ages_mean)


##############################################
## LOOP SIMPER AND INDISPECIES BY ECOREGION ##
##############################################

simper_table <- NULL
indicator_table <- NULL

for (i in unique(survey_sp_occ$Ecoregion)){
  
  sub_eco <- survey_sp_occ %>%
    filter(Ecoregion==i) %>%
    arrange(MPA)
  
  sub_group <- sub_eco$MPA
  sub_group <- droplevels(sub_group)
  
  sub_fish <- sub_eco %>% # ONLY TAKE FISH WITH ACTUAL OCCURENCE
    ungroup() %>%
    select(where( ~ is.numeric(.x) && sum(.x) != 0))
  
  indval = multipatt(sub_fish, sub_group, control = how(nperm=999), func = "r.g")
  #summary(indval, alpha=1) # SIG P VALUE SHORTENS THE LIST
  ind_species <- indval$sign
  ind_species$species <- rownames(ind_species)
  
  indicator_table <- rbind(indicator_table, data.frame(Ecoregion=i,
                                                       ind_species))
  
  simp = simper(sub_fish, sub_group, permutations = 999)
  simp <- as.data.frame(simp$`Fishing_No take`)
  simp$MPA_species <- ifelse(simp$avb > simp$ava, "yes", "no")
  #simp <- simp %>%
  #  filter(MPA_species == "yes")
  
  simper_table <- rbind(simper_table, data.frame(Ecoregion=i,
                                                 simp[,c("species","average","ratio","MPA_species")]))
  
}


simper_table$species <- gsub("_"," ", simper_table$species)
simper_table$scientificname <- gsub(" spp.", "", simper_table$species)
simper_table$family <- simper_table$scientificname
simper_table <- FindReplace(simper_table, "family", family_info, from="scientificname", to="family",
                            exact = T, vector = F)
sum(is.na(simper_table$family))
simper_table <- merge(simper_table, species_ages, by="scientificname")


indicator_table$species <- gsub("_"," ", indicator_table$species)
indicator_table$scientificname <- gsub(" spp.", "", indicator_table$species)
indicator_table$family <- indicator_table$scientificname
indicator_table <- FindReplace(indicator_table, "family", family_info, from="scientificname", to="family",
                               exact = T, vector = F)
sum(is.na(indicator_table$family))
indicator_table$MPA_species <- ifelse(indicator_table$s.No.take==1, "yes", "no")
indicator_table <- merge(indicator_table, species_ages, by="scientificname")



#' #########################################################
#' CALCULATE AVERAGE VALUE PER FAMILY PER ECOREGION 
#' THEN MAKE BOXPLOT
#' HOW DOES THIS COMPARE TO ONE GLOBAL BOXPLOT
#' #########################################################

simp_family <- simper_table %>%
  select(Ecoregion, family, MPA_species, average) %>%
  group_by(Ecoregion, MPA_species, family) %>%
  summarise_all(.funs=mean)
length(unique(simp_family$family))
simp_family <- merge(simp_family, family_ages, by="family")

#   WHICH FAMILIES HAVE THE BIGGEST DISPARITY?
# simp_family_mean <- simper_table %>%
#   select(family, MPA_species, average) %>%
#   group_by(family, MPA_species) %>%
#   summarise_all(.funs=mean)

indicator_family <- indicator_table %>%
  select(Ecoregion, family, MPA_species, stat) %>%
  group_by(Ecoregion, MPA_species, family) %>%
  summarise_all(.funs=mean)
length(unique(indicator_family$family))
indicator_family <- merge(indicator_family, family_ages, by="family")

## HOW TO DO THESE PLOTS - COLOR BY PHYLOGENY SOMEHOW?
## CUT TO ONLY A CERTAIN NUMBER OF FAMILIES?


## FIRST SIMPER ANALYSIS

## MPA CONTRIBUTORS ##
ggplot(data=
         filter(simp_family, MPA_species=="yes"),
       aes(x = reorder(family, log(average), fun = median, .desc =TRUE), 
           y = log(average))) + 
  geom_boxplot(aes(fill=log_age)) + 
  #geom_jitter(position=position_jitter(0.2)) +
  theme_bw(base_size = 10) +
  xlab("Family") +
  ylab("") +
  scale_fill_gradient(low = "red",
                      high = "blue") +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  ggtitle("Family Contributions to MPA Sites")

## FISHED SITE CONTRIBUTORS ##
ggplot(data=
         filter(simp_family, MPA_species=="no"),
       aes(x = reorder(family, ages_mean, fun = median, .desc =TRUE), y = average)) + 
  geom_boxplot(aes(fill = reorder(family, ages_mean, fun = median, .desc =TRUE))) + 
  #geom_jitter(position=position_jitter(0.2)) +
  theme_bw(base_size = 10) +
  xlab("Family") +
  ylab("") +
  scale_fill_discrete(guide = guide_legend(title = "Family")) + 
  theme(legend.position = "none")  +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  ggtitle("Family Contributions to Fished Sites")



### PAIRED BOXPLOT
ggplot(data=simp_family, 
       aes(x=family, y=log(average), fill=MPA_species)) +
  geom_boxplot() +
  scale_x_discrete(guide = guide_axis(angle = 90))




## INDICATOR ANALYSIS

## MPA CONTRIBUTORS ##
ggplot(data=
         filter(indicator_family, MPA_species=="yes"),
       aes(x = reorder(family, -stat, fun = median, .desc =TRUE), y = stat)) + 
  geom_boxplot(aes(fill = reorder(family, stat, fun = median, .desc =TRUE))) + 
  #geom_jitter(position=position_jitter(0.2)) +
  theme_bw(base_size = 10) +
  xlab("Family") +
  ylab("") +
  #scale_fill_discrete(guide = guide_legend(title = "Family")) + 
  theme(legend.position = "none")  +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  ggtitle("Family Contributions to MPA Sites")


## FISHED SITE CONTRIBUTORS ##
ggplot(data=
         filter(indicator_family, MPA_species=="no"),
       aes(x = reorder(family, -stat, fun = median, .desc =TRUE), y = stat)) + 
  geom_boxplot(aes(fill = reorder(family, stat, fun = median, .desc =TRUE))) +
  theme_bw(base_size = 10) +
  xlab("Family") +
  ylab("") +
  #scale_fill_discrete(guide = guide_legend(title = "Family")) + 
  theme(legend.position = "none")  +
  scale_x_discrete(guide = guide_axis(angle = 90)) + 
  ggtitle("Family Contributions to Fished Sites")





