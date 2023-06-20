
#######################################################################################
#'  CODE TO DISENTANGLE THE MECHANISMS BY WHICH MPA AFFECTS AESTHETIC VALUE
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

##############################
## IMPORT DATA "MODEL DATA" ##
##############################

model_data <- readRDS("data/model_data.rds")

model_data$Temperature_Zone <- as.factor(model_data$Temperature_Zone)

names(model_data)

model_data$MPA <- as.factor(model_data$MPA)

model_data$sst_range <- model_data$sst_max - model_data$sst_min

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

##################################
## ADD NEW IMPUTED BENTHIC DATA ##
##################################

benthic_PCA <- readRDS("outputs/RLS_benthic_PCA_imputed.rds")
benthic_PCA <- benthic_PCA %>%
  select(SurveyID, PC1_imputation, PC2_imputation) %>%
  rename(PC1_imputed = PC1_imputation, PC2_imputed = PC2_imputation)

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

############################
## HOW MUCH MISSING DATA? ##
############################

sapply(model_data, function(x) paste(round(sum(is.na(x))/length(x),2)*100,"%",sep=""))


##############################################
## FIRST EXAMINE WHETHER RICHNESS IS HIGHER
## IN MPAS THAN FISHED SITES 
##############################################

# CALCULATE AVERAGE SP RICHNESS PER SITE
site_info <- model_data %>%
  select(SiteCode, SiteLongitude, SiteLatitude, Country, Ecoregion, MPA)
site_info <- site_info[!duplicated(site_info$SiteCode),]

site_richness <- model_data %>%
  select(SiteCode, nb_species) %>%
  group_by(SiteCode) %>%
  summarise_all(.funs=mean)

site_richness <- merge(site_richness, site_info, by="SiteCode")

# MAKE A DATA FRAME TO FOCUS ONLY ON NO TAKE MPAS
# ONLY KEEP ECOREGIONS WITH BOTH FISHED AND NO-TAKE SITES
no_take_only <- site_richness %>%
  filter(MPA != "Restricted take")
no_take_only$MPA <- droplevels(no_take_only$MPA)

selected_regions <- no_take_only %>%
  group_by(Ecoregion, MPA) %>%
  summarise(n()) %>%
  group_by(Ecoregion) %>%
  filter(n()>1)

length(unique(model_data$Ecoregion))
length(unique(selected_regions$Ecoregion)) 
# ONLY 43/96 REGIONS HAVE BOTH TYPES OF SITES

no_take_only <- no_take_only %>%
  filter(Ecoregion %in% selected_regions$Ecoregion)
length(unique(no_take_only$Ecoregion))

# SIMPLE BOXPLOTS
graphics.off()

# ALL DATA
boxplot(site_richness$nb_species ~ site_richness$MPA)
boxplot(log(site_richness$nb_species) ~ site_richness$MPA)

# NO TAKE, FILTERED REGIONS
boxplot(no_take_only$nb_species ~ no_take_only$MPA)
boxplot(log(no_take_only$nb_species) ~ no_take_only$MPA)

# ON GLOBAL POOLED SCALE - YES, MPAS ARE HIGHER
# FOR FILTERED REGIONS - YES

# ACROSS ECOREGIONS
ggplot(data = no_take_only,
       aes(x = MPA,
           y = nb_species,
           group = Ecoregion,
           color = factor(Ecoregion))) +
  geom_line(size = 2,
            alpha = 0.5) +
  geom_point(size = 3) +
  theme(legend.position = "none")  

graphics.off()
par(mfrow=c(3,3))
rich_count <- NULL
for(i in unique(no_take_only$Ecoregion)){
  sub_eco <- no_take_only %>%
    filter(Ecoregion==i)
  boxplot(sub_eco$nb_species ~ sub_eco$MPA)
  eco_mean <- sub_eco %>%
    select(Ecoregion, MPA, nb_species) %>%
    group_by(Ecoregion, MPA) %>%
    summarise_all(.funs=mean)
  eco_mean$MPA <- as.character(eco_mean$MPA)
  rich_count <- rbind(rich_count, 
                      which_higher=print(eco_mean$MPA[which.max(eco_mean$nb_species)]))
}

table(rich_count)

# 30 ECOREGION HAVE HIGHER AVERAGE RICHNESS IN MPA 
# 13 HAVE HIGHER AVERAGE RICHNESS IN FISHED





#############################
## PREPARE DATA FOR MODELS ##
#############################
###############################################
## DO ANY VARIABLES NEED TO BE TRANSFORMED ? ##
###############################################

model_data <- do.call(data.frame,lapply(model_data, function(x) replace(x, is.infinite(x),NA)))

num_vars <- select_if(model_data, is.numeric)
num_vars$SurveyID <- NULL

# graphics.off()
# par(mfrow=c(4,4))
# for(i in 1:ncol(num_vars)){
#   hist(num_vars[,i], main=colnames(num_vars)[i])
# }

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

##################################
## PARALLEL SETTINGS FOR MODELS ##
##################################

ncores = detectCores()
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#####################################################
## SCALE ALL THE NUMERIC PREDICTORS TO MEAN 0 SD 2 ##
#####################################################

z_score_2sd <- function(x){ (x - mean(x,na.rm=T)) / (2*sd(x,na.rm=T))}

# WHAT ABOUT THE TROPHIC TRAITS?
z_vars <- model_data %>%
  select_if(is.numeric) %>%
  select(-any_of(c("SurveyID","SiteLongitude","SiteLatitude","aesthe_survey"))) %>%
  colnames()

standardized_data <- model_data %>%
  mutate_if(colnames(model_data) %in% z_vars, z_score_2sd)


########################################################
## MODEL TO ASSESS WHEHTER RICHNESS IS HIGHER IN MPAS ##
## WITH SITE NESTED IN COUNTRY AS RANDOM EFFECTS
########################################################

richness_model_formula  <- bf(nb_species ~ MPA +
                                Temperature_Zone +
                                (1 | Country/SiteCode),
                              
                              family=gaussian())

richness_model <- brm(richness_model_formula,
                     
                       data=model_data,
                      
                      chains=1, iter=3000, cores=ncores,
                      c(set_prior("normal(0,3)", class = "b"),
                        set_prior("normal(0,3)", class="Intercept")))

conditional_effects(richness_model, "MPA")
summary(richness_model)

# POSITIVE EFFECT OF MPA ON RICHNESS

########################################################
## MODEL TO ASSESS WHEHTER RICHNESS IS HIGHER IN MPAS 
## WITH SITE NESTED IN ECOREGION AS RANDOM EFFECTS
########################################################

richness_model_formula  <- bf(nb_species ~ MPA +
                                Temperature_Zone +
                                (1 | Ecoregion/SiteCode),
                              
                              family=gaussian())

richness_model <- brm(richness_model_formula,
                      
                      data=model_data,
                     
                       chains=1, iter=3000, cores=ncores,
                      c(set_prior("normal(0,3)", class = "b"),
                        set_prior("normal(0,3)", class="Intercept")))

summary(richness_model)
conditional_effects(richness_model)

# POSITIVE EFFECT OF MPA ON RICHNESS


########################################################
## MODEL TO ASSESS WHEHTER RICHNESS IS HIGHER IN MPAS 
## USING SITE NESTED IN ECOREGION AS RANDOM EFFECTS
## FOR ONLY ECOREGIONS WITH BOTH TYPES OF SITES
########################################################

sub_eco_data <- model_data %>%
  filter(Ecoregion %in% selected_regions$Ecoregion)

richness_model_formula  <- bf(nb_species ~ MPA +
                                Temperature_Zone +
                                (1 | Ecoregion/SiteCode),
                              
                              family=gaussian())

richness_model <- brm(richness_model_formula,
                      
                      data=sub_eco_data,
                     
                       chains=1, iter=3000, cores=ncores,
                      c(set_prior("normal(0,3)", class = "b"),
                        set_prior("normal(0,3)", class="Intercept")))

summary(richness_model)
conditional_effects(richness_model)

# POSITIVE EFFECT OF MPA ON RICHNESS


###################################################################
#' TO ASSESS THE EFFECT OF EACH PATHWAY OF THE
#' TOTAL CAUSAL EFFECT OF MPA, PERFORM THE FOLLOWING
#' CALCULATE THE TOTAL CAUSAL EFFECT OF MPA ACCORDING TO THE DAG
#' THEN TEST THE INDEPENDENT CONTRIBUTION OF EACH PATHWAY
#' BY BLOCKING THE OTHER PATHWAYS
#' THE PATHWAYS OF MPA EFFECT ARE:
#' TAXO DIVERSITY, FUN DIVERSITY, PHYLO DIVERSITY,
#' BENTHIC COMPOSITION, BIOMASS, TROPHIC STRUCTURE
####################################################################


###################################################
## MPA TOTAL CAUSAL EFFECT FROM DAG-BASED MODELS ##
###################################################

MPA_model_formula  <- bf(log(aesthe_survey) ~ MPA +
                           HDI2017 +
                           fshD +
                           gravtot2 +
                           as.factor(Temperature_Zone) +
                           (1 | Country/SiteCode),
                         
                         family=gaussian())

MPA_model <- brm(MPA_model_formula,
                 data=model_data,
                 chains=3, iter=3000, cores=ncores,
                 c(set_prior("normal(0,3)", class = "b"),
                   set_prior("normal(0,3)", class="Intercept")))

MPA_post <- as.data.frame(as.matrix(MPA_model)) %>%
  select('b_MPANotake') %>%
  rename("MPA Total Causal Effect" = b_MPANotake)


#############################
## MPA "PURE" CONTRIBUTION ##
#############################

MPA_pure_model_formula  <- bf(log(aesthe_survey) ~ MPA +
                               HDI2017 + # CONTROL
                               fshD + # CONTROL
                               gravtot2 + # CONTROL
                               Temperature_Zone + # CONTROL
                               
                                taxo_entropy + # BLOCKED
                               fun_entropy + # BLOCKED
                               phylo_entropy + # BLOCKED
                               Biomass + # BLOCKED
                               PC1_imputed + # BLOCKED
                               PC2_imputed + # BLOCKED
                               PC1_trophic + # BLOCKED
                               PC2_trophic + # BLOCKED
                               
                               (1 | Country/SiteCode),
                             
                             family=gaussian())

MPA_pure_model <- brm(MPA_pure_model_formula,
                     
                     data=standardized_data,
                     
                     chains=3, iter=3000, cores=ncores,
                     c(set_prior("normal(0,3)", class = "b"),
                       set_prior("normal(0,3)", class="Intercept")))

MPA_pure_post <- as.data.frame(as.matrix(MPA_pure_model)) %>%
  select('b_MPANotake') %>% 
  rename("MPA Pure Direct Effect" = b_MPANotake)



#################################
## TAXO DIVERSITY CONTRIBUTION ##
#################################

tax_div_model_formula  <- bf(log(aesthe_survey) ~ MPA +
                           HDI2017 + # CONTROL
                           fshD + # CONTROL
                           gravtot2 + # CONTROL
                           Temperature_Zone + # CONTROL
                             
                             fun_entropy + # BLOCKED
                             phylo_entropy + # BLOCKED
                             Biomass + # BLOCKED
                             PC1_imputed + # BLOCKED
                             PC2_imputed + # BLOCKED
                             PC1_trophic + # BLOCKED
                             PC2_trophic + # BLOCKED
                             
                           (1 | Country/SiteCode),
                         
                         family=gaussian())

tax_div_model <- brm(tax_div_model_formula,
                     
                 data=standardized_data,
                 
                 chains=3, iter=3000, cores=ncores,
                 c(set_prior("normal(0,3)", class = "b"),
                   set_prior("normal(0,3)", class="Intercept")))

tax_div_post <- as.data.frame(as.matrix(tax_div_model)) %>%
  select('b_MPANotake') %>% 
  rename("Taxonomic Diversity" = b_MPANotake)



################################
## FUN DIVERSITY CONTRIBUTION ##
################################

fun_div_model_formula  <- bf(log(aesthe_survey) ~ MPA +
                               HDI2017 + # CONTROL
                               fshD + # CONTROL
                               gravtot2 + # CONTROL
                               Temperature_Zone + # CONTROL
                               
                               taxo_entropy + # BLOCKED
                               phylo_entropy + # BLOCKED
                               Biomass + # BLOCKED
                               PC1_imputed + # BLOCKED
                               PC2_imputed + # BLOCKED
                               PC1_trophic + # BLOCKED
                               PC2_trophic + # BLOCKED
                               
                               (1 | Country/SiteCode),
                             
                             family=gaussian())

fun_div_model <- brm(fun_div_model_formula,
                     
                     data=standardized_data,
                     
                     chains=3, iter=3000, cores=ncores,
                     c(set_prior("normal(0,3)", class = "b"),
                       set_prior("normal(0,3)", class="Intercept")))

fun_div_post <- as.data.frame(as.matrix(fun_div_model)) %>%
  select('b_MPANotake') %>%
  rename("Functional Diversity" = b_MPANotake)



##################################
## PHYLO DIVERSITY CONTRIBUTION ##
##################################

phylo_div_model_formula  <- bf(log(aesthe_survey) ~ MPA +
                               HDI2017 + # CONTROL
                               fshD + # CONTROL
                               gravtot2 + # CONTROL
                               Temperature_Zone + # CONTROL
                               
                               taxo_entropy + # BLOCKED
                               fun_entropy + # BLOCKED
                               Biomass + # BLOCKED
                               PC1_imputed + # BLOCKED
                               PC2_imputed + # BLOCKED
                                 PC1_trophic + # BLOCKED
                                 PC2_trophic + # BLOCKED
                                 
                               (1 | Country/SiteCode),
                             
                             family=gaussian())

phylo_div_model <- brm(phylo_div_model_formula,
                       
                     data=standardized_data,
                     
                     chains=3, iter=3000, cores=ncores,
                     c(set_prior("normal(0,3)", class = "b"),
                       set_prior("normal(0,3)", class="Intercept")))

phylo_div_post <- as.data.frame(as.matrix(phylo_div_model)) %>%
  select('b_MPANotake') %>%
  rename("Phylogenetic Diversity" = b_MPANotake)




##########################
## BIOMASS CONTRIBUTION ##
##########################

biomass_model_formula  <- bf(log(aesthe_survey) ~ MPA +
                                 HDI2017 + # CONTROL
                                 fshD + # CONTROL
                                 gravtot2 + # CONTROL
                                 Temperature_Zone + # CONTROL
                                 
                                 taxo_entropy + # BLOCKED
                                 fun_entropy + # BLOCKED
                                 phylo_entropy + # BLOCKED
                                 PC1_imputed + # BLOCKED
                                 PC2_imputed + # BLOCKED
                               PC1_trophic + # BLOCKED
                               PC2_trophic + # BLOCKED
                               
                                 (1 | Country/SiteCode),
                               
                               family=gaussian())

biomass_model <- brm(biomass_model_formula,
                     
                       data=standardized_data,
                     
                       chains=3, iter=3000, cores=ncores,
                       c(set_prior("normal(0,3)", class = "b"),
                         set_prior("normal(0,3)", class="Intercept")))

biomass_post <- as.data.frame(as.matrix(biomass_model)) %>%
  select('b_MPANotake') %>%
  rename("Biomass" = b_MPANotake)




##########################
## BENTHIC CONTRIBUTION ##
##########################

## DOES THIS NEED TO BE SPLIT BY PC1 AND PC2?

benthic_model_formula  <- bf(log(aesthe_survey) ~ MPA +
                               HDI2017 + # CONTROL
                               fshD + # CONTROL
                               gravtot2 + # CONTROL
                               Temperature_Zone + # CONTROL
                               
                               taxo_entropy + # BLOCKED
                               fun_entropy + # BLOCKED
                               phylo_entropy + # BLOCKED
                               Biomass + # BLOCKED
                                 PC1_trophic + # BLOCKED
                                 PC2_trophic + # BLOCKED
                                 
                               (1 | Country/SiteCode),
                             
                             family=gaussian())

benthic_model <- brm(benthic_model_formula,
                     
                     data=standardized_data,
                     
                     chains=3, iter=3000, cores=ncores,
                     c(set_prior("normal(0,3)", class = "b"),
                       set_prior("normal(0,3)", class="Intercept")))

benthic_post <- as.data.frame(as.matrix(benthic_model)) %>%
  select('b_MPANotake') %>%
  rename("Benthic Composition" = b_MPANotake)


##########################
## TROPHIC CONTRIBUTION ##
##########################

trophic_model_formula  <- bf(log(aesthe_survey) ~ MPA +
                               HDI2017 + # CONTROL
                               fshD + # CONTROL
                               gravtot2 + # CONTROL
                               Temperature_Zone + # CONTROL
                               
                               taxo_entropy + # BLOCKED
                               fun_entropy + # BLOCKED
                               phylo_entropy + # BLOCKED
                               Biomass + # BLOCKED
                               PC1_imputed + # BLOCKED
                               PC2_imputed + # BLOCKED
                               
                               (1 | Country/SiteCode),
                             
                             family=gaussian())

trophic_model <- brm(trophic_model_formula,
                     
                     data=standardized_data,
                     
                     chains=3, iter=3000, cores=ncores,
                     c(set_prior("normal(0,3)", class = "b"),
                       set_prior("normal(0,3)", class="Intercept")))

trophic_post <- as.data.frame(as.matrix(trophic_model)) %>%
  select('b_MPANotake') %>%
  rename("Trophic Structure" = b_MPANotake)




###########################
## GROUP MODEL ESTIMATES ##
###########################

model_outputs <- data.frame(MPA_pure_post, tax_div_post, fun_div_post, phylo_div_post,
                            biomass_post, 
                            benthic_post, trophic_post)

saveRDS(model_outputs, "MPA_mechanism_effects.rds")

model_estimates <- data.frame(median=apply(model_outputs, 2, median))
model_estimates$abs_effect <- abs(model_estimates$median)
model_estimates <- model_estimates %>%
  arrange(desc(abs_effect))
model_outputs <- model_outputs[,order(match(colnames(model_outputs), rownames(model_estimates)))]

mcmc_intervals(model_outputs)

sum(model_estimates$median)

colMeans(MPA_post)

# SUM OF INDIVIDUAL PATHWAYS IS 2X HIGHER THAN TOTAL MPA EFFECT - WHY? OVERLAP??
