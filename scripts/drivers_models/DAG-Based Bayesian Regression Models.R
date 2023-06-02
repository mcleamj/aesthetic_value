#######################################################################################
#'  BAYSESIAN HIERARCHICAL MODELS TO ASSESS DRIVERS OF AESTHETIC VALUE
#'  BASED ON DIRECTED ACYCLIC GRAPHS (DAG)
#'  THE CORRESPONDING DAG CAN BE FOUND AT:
#'  
#'  dagitty.net/mxttk25
#'
#' @author Matthew McLean, \email {mcleamj@gmail.com},
#'         
#' @date
########################################################################################

######################
## LIBRARY PACKAGES ##
######################

if(!require(data.table)){install.packages("data.table"); library(data.table)}
if(!require(mapproj)){install.packages("mapproj"); library(mapproj)}
if(!require(plot3D)){install.packages("plot3D"); library(plot3D)}
if(!require(pals)){install.packages("pals"); library(pals)}
if(!require(brms)){install.packages("brms"); library(brms)}
if(!require(bayesplot)){install.packages("bayesplot"); library(bayesplot)}
if(!require(parallel)){install.packages("parallel"); library(parallel)}
if(!require(rstan)){install.packages("rstan"); library(rstan)}
if(!require(bayestestR)){install.packages("bayestestR"); library(bayestestR)}
if(!require(readr)){install.packages("readr"); library(readr)}
if(!require(tibble)){install.packages("tibble"); library(tibble)}
if(!require(dplyr)){install.packages("dplyr"); library(dplyr)}

##############################
## IMPORT DATA "MODEL DATA" ##
##############################

model_data <- readRDS("data/model_data.rds")

names(model_data)

model_data$MPA <- as.factor(model_data$MPA)

model_data$sst_range <- model_data$sst_max - model_data$sst_min

######################################
## IMPORT AND ADD DIVERSITY METRICS ##
######################################

biodiv <- read_rds("outputs/survey_biodiversity.rds")

model_data <- merge(model_data, biodiv, by="SurveyID")

##################################
## ADD NEW IMPUTED BENTHIC DATA ##
##################################

benthic_PCA <- read.table("outputs/RLS_benthic_PCA_imputed.txt")
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

aaesthe_surveyetic_data <- read.csv("outputs/survey_aesth.csv")

model_data <- merge(model_data, aaesthe_surveyetic_data, by="SurveyID")

############################
## HOW MUCH MISSING DATA? ##
############################

sapply(model_data, function(x) paste(round(sum(is.na(x))/length(x),2)*100,"%",sep=""))

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

z_vars <- model_data %>%
  select_if(is.numeric) %>%
  select(-any_of(c("SurveyID","SiteLongitude","SiteLatitude","aesthe_survey"))) %>%
  colnames()

model_data <- model_data %>%
  mutate_if(colnames(model_data) %in% z_vars, z_score_2sd)


##########################
## INTERCEPT ONLY MODEL ##
##########################

int_only_formula <- bf(log(aesthe_survey) ~ 
                         (1 | Country/SiteCode),
                       
                       family = gaussian())

int_only_model <- brm(int_only_formula,
                      data=model_data,
                      chains=3, iter=4000, cores=ncores,
                      refresh=500,
                      set_prior("normal(0,3)", class="Intercept"))

int_output <- data.frame(as.matrix(int_only_model)) ## NEED TO BACK CALCULATE TO RAW VALUE

site_estimates <- int_output %>%
  select(starts_with("r_Country.SiteCode"))

SiteCode = sapply(strsplit(names(site_estimates),split="[_]"), function(x) x[3])
SiteCode <- gsub(".Intercept.","",SiteCode)
SiteCode <- gsub("\\.","-",SiteCode)

names(site_estimates) <- SiteCode

site_estimates <- exp(site_estimates + int_output[,1]) ## BACK CALCULATE TO RAW VALUE

site_estimates <- as.data.frame(posterior_summary(site_estimates))
site_estimates$SiteCode <- rownames(site_estimates)

site_info <- model_data %>%
  select(SiteCode, SiteLongitude, SiteLatitude)
site_info <- site_info[!duplicated(site_info$SiteCode),]

site_estimates <- merge(site_info, site_estimates, by="SiteCode")
write.table(site_estimates,"int_only_site_estimates.txt")

raw_site_mean <- model_data %>%
  select(SiteCode, SiteLongitude, SiteLatitude, aesthe_survey) %>%
  group_by(SiteCode) %>%
  summarise_all(mean)
write.table(raw_site_mean, "raw_site_mean.txt")

compare_site <- merge(site_estimates, raw_site_mean, by="SiteCode")
plot(log(compare_site$aesthe_survey), log(compare_site$Estimate))

########################
## CAUSAL SALAD MODEL ##
########################

full_model_formula <- 
  bf(log(aesthe_survey) ~
       as.factor(Temperature_Zone) +
       sst_mean +
       NPP_mean +
       Depth +
       fshD +
       dhw_mean +
       PC1_imputed +
       PC2_imputed +
       HDI2017 +
       gravtot2 +
       MPA +
       (1 | Country/SiteCode),
     family = gaussian()) 

full_model <- brm(full_model_formula,
                  data=model_data,
                  chains=3, iter=3000, cores=ncores,
                  refresh=500,
                  c(set_prior("normal(0,3)", class = "b"),
                    set_prior("normal(0,3)", class="Intercept")))

r2_bayes(full_model)
full_post <- as.data.frame(as.matrix(full_model)) %>%
  select('b_sst_mean':'b_MPARestrictedtake')

full_estimates <- data.frame(median=apply(full_post, 2, median))
full_estimates$abs_effect <- abs(full_estimates$median)
full_estimates <- full_estimates %>%
  arrange(desc(abs_effect))

full_post <- full_post[,order(match(colnames(full_post), rownames(full_estimates)))]
mcmc_intervals(full_post)

################
## SST  MODEL ##
################

sst_model_formula <- bf(log(aesthe_survey) ~ sst_mean +
                          abs_latitude +
                          as.factor(Temperature_Zone) +
                          (1 | Country/SiteCode),
  
                        family=gaussian())

sst_model <- brm(sst_model_formula,
                 data=model_data,
                 chains=3, iter=3000, cores=ncores,
                 c(set_prior("normal(0,3)", class = "b"),
                   set_prior("normal(0,3)", class="Intercept")))

sst_post <- as.data.frame(as.matrix(sst_model)) %>%
  select('b_sst_mean')

###################
## GRAVITY MODEL ##
###################

gravity_model_formula <- bf(log(aesthe_survey) ~ gravtot2 +
                              HDI2017 +
                              abs_latitude + # ABSOLUTE VALUE FOR LINEARITY
                              as.factor(Temperature_Zone) +
                              (1 | Country/SiteCode),
                            family=gaussian())

gravity_model <- brm(gravity_model_formula,
                     data=model_data,
                     chains=3, iter=3000, cores=ncores,
                     c(set_prior("normal(0,3)", class = "b"),
                       set_prior("normal(0,3)", class="Intercept")))

gravity_post <- as.data.frame(as.matrix(gravity_model)) %>%
  select('b_gravtot2')

###############
## NPP MODEL ##
###############

NPP_model_formula <- bf(log(aesthe_survey) ~ NPP_mean +
                          abs_latitude +
                          BO_nitrate +
                          BO_phosphate +
                          sst_mean +
                          wave_energy +
                          as.factor(Temperature_Zone) +
                          (1 | Country/SiteCode),
                          
                        family=gaussian())

NPP_model <- brm(NPP_model_formula,
                 data=model_data,
                 chains=3, iter=3000, cores=ncores,
                 c(set_prior("normal(0,3)", class = "b"),
                   set_prior("normal(0,3)", class="Intercept")))

NPP_post <- as.data.frame(as.matrix(NPP_model)) %>%
  select('b_NPP_mean')

#################
## DEPTH MODEL ##
#################

depth_model_formula <- bf(log(aesthe_survey) ~ Depth +
                            as.factor(Temperature_Zone) +
                            (1 | Country/SiteCode),
                            
                          family=gaussian())

depth_model <- brm(depth_model_formula,
                   data=model_data,
                   chains=3, iter=3000, cores=ncores,
                   c(set_prior("normal(0,3)", class = "b"),
                     set_prior("normal(0,3)", class="Intercept")))

depth_post <- as.data.frame(as.matrix(depth_model)) %>%
  select('b_Depth')

#######################
## FISH DEPEND MODEL ##
#######################

fshd_model_formula <- bf(log(aesthe_survey) ~ fshD +
                           HDI2017 +
                           as.factor(Temperature_Zone) +
                           (1 | Country/SiteCode),
                         
                         family=gaussian())

fshd_model <- brm(fshd_model_formula,
                  data=model_data,
                  chains=3, iter=3000, cores=ncores,
                  c(set_prior("normal(0,3)", class = "b"),
                    set_prior("normal(0,3)", class="Intercept")))

fshd_post <- as.data.frame(as.matrix(fshd_model)) %>%
  select('b_fshD')

#################
## HDI   MODEL ##
#################

HDI_model_formula <- bf(log(aesthe_survey) ~ HDI2017 +
                          abs_latitude + 
                          as.factor(Temperature_Zone) +
                          (1 | Country/SiteCode),
                          
                        family=gaussian())

HDI_model <- brm(HDI_model_formula,
                 data=model_data,
                 chains=3, iter=3000, cores=ncores,
                 c(set_prior("normal(0,3)", class = "b"),
                   set_prior("normal(0,3)", class="Intercept")))

HDI_post <- as.data.frame(as.matrix(HDI_model)) %>%
  select('b_HDI2017')

###############
## MPA MODEL ##
###############

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
  select('b_MPANotake','b_MPARestrictedtake')

#################################
## BENTHIC COMPOSITION MODEL 1 ##
#################################

benthic_PC1_model_formula <- 
  bf(log(aesthe_survey) ~ PC1_imputed +
       Depth +
       abs_latitude +
       MPA + 
       NPP_mean +
       sst_mean +
       as.factor(Temperature_Zone) +
       (1 | Country/SiteCode),
       
     family=gaussian()) 

benthic_PC1_model <- brm(benthic_PC1_model_formula,
                         data=model_data,
                         chains=3, iter=3000, cores=ncores,
                         c(set_prior("normal(0,3)", class = "b"),
                           set_prior("normal(0,3)", class="Intercept")))

benthic_PC1_post <- as.data.frame(as.matrix(benthic_PC1_model))  %>%
  select('b_PC1_imputed')

#################################
## BENTHIC COMPOSITION MODEL 2 ##
#################################

benthic_PC2_model_formula <- 
  bf(log(aesthe_survey) ~ PC2_imputed +
       Depth +
       abs_latitude +
       MPA + 
       NPP_mean +
       sst_mean +
       as.factor(Temperature_Zone) +
       (1 | Country/SiteCode),
       
     family=gaussian()) 

benthic_PC2_model <- brm(benthic_PC2_model_formula,
                         data=model_data,
                         chains=3, iter=3000, cores=ncores,
                         c(set_prior("normal(0,3)", class = "b"),
                           set_prior("normal(0,3)", class="Intercept")))

benthic_PC2_post <- as.data.frame(as.matrix(benthic_PC2_model)) %>%
  select('b_PC2_imputed')

###############
## DHW MODEL ##
###############

DHW_model_formula <- bf(log(aesthe_survey) ~ dhw_mean +
                          abs_latitude +
                          sst_mean +
                          as.factor(Temperature_Zone) +
                          (1 | Country/SiteCode),
                          
                        family=gaussian())

DHW_model <- brm(DHW_model_formula,
                 data=model_data,
                 chains=3, iter=3000, cores=ncores,
                 c(set_prior("normal(0,3)", class = "b"),
                   set_prior("normal(0,3)", class="Intercept")))

DHW_post <- as.data.frame(as.matrix(DHW_model)) %>%
  select('b_dhw_mean')

###############################################
## RECOMBINE DAG MODELS AND RENAME VARIABLES ##
###############################################

dag_output <- data.frame(sst_post, gravity_post, NPP_post, depth_post,
                         fshd_post, HDI_post, MPA_post,
                         benthic_PC1_post, benthic_PC2_post, DHW_post)

write.table(dag_output, "dag_output.txt")
write.table(full_post, "full_post.txt")

names(full_post) <- gsub("b_", "", names(full_post))
names(dag_output) <- gsub("b_", "", names(dag_output))

dag_output <- dag_output %>%
  rename("Log Human Gravity" = gravtot2,
         "No Take MPA" = MPANotake,
         "Restricted Take MPA" = MPARestrictedtake,
         "Sea Surface Temperature" = sst_mean,
         "Net Primary Productivity" = NPP_mean,
         "Degree Heating Weeks" = dhw_mean,
         "Depth" =  Depth,
         "Fisheries Dependency" = fshD,
         "Human Dev. Index" = HDI2017,
         "Benthic Composition (PC1)" = PC1_imputed,
         "Benthic Composition (PC2)" = PC2_imputed)


full_post <- full_post %>%
  rename("Log Human Gravity" = gravtot2,
    "No Take MPA" = MPANotake,
    "Restricted Take MPA" = MPARestrictedtake,
    "Sea Surface Temperature" = sst_mean,
    "Net Primary Productivity" = NPP_mean,
    "Degree Heating Weeks" = dhw_mean,
    "Depth" =  Depth,
    "Fisheries Dependency" = fshD,
    "Human Dev. Index" = HDI2017,
    "Benthic Composition (PC1)" = PC1_imputed,
    "Benthic Composition (PC2)" = PC2_imputed)

dag_estimates <- data.frame(median=apply(dag_output, 2, median))
dag_estimates$abs_effect <- abs(dag_estimates$median)
dag_estimates <- dag_estimates %>%
  arrange(desc(abs_effect))

# MATCH DAG ORDER TO CAUSAL SALAD ORDER?
dag_output <- dag_output[,order(match(colnames(dag_output), colnames(full_post)))]

#########################
## SIMPLE FOREST PLOTS ##
#########################

color_scheme_set("darkgray")
mcmc_areas_ridges(full_post) +
  ggtitle("Causal Salad Model") +
  xlim(range(full_post, dag_output))

mcmc_areas_ridges(dag_output) +
  ggtitle("DAG-Based Models") +
  xlim(range(full_post, dag_output))

#########################
## SIMPLE FOREST PLOTS ##
#########################

color_scheme_set("darkgray")
mcmc_intervals(full_post) +
  ggtitle("Causal Salad Model") +
  xlim(range(full_post, dag_output))

mcmc_intervals(dag_output) +
  ggtitle("DAG-Based Models") +
  xlim(range(full_post, dag_output))

ggpubr::ggarrange(mcmc_intervals(full_post) +
                    ggtitle("Causal Salad Model") +
                    xlim(range(full_post, dag_output)),
                  
                  mcmc_intervals(dag_output) +
                    ggtitle("DAG-Based Models") +
                    xlim(range(full_post, dag_output)),
                  
                  ncol = 2)

####################################
## RANK DAG ORDER BY EFFECT SIZE? ##
####################################

dag_output <- dag_output[,order(match(colnames(dag_output), rownames(dag_estimates)))]
full_post <- full_post[,order(match(colnames(full_post), colnames(dag_output)))]

#########################
## SIMPLE FOREST PLOTS ##
#########################

ggpubr::ggarrange(
  
  mcmc_intervals(dag_output) +
    ggtitle("DAG-Based Models") +
    xlim(range(full_post, dag_output)),
  
  mcmc_intervals(full_post) +
    ggtitle("Causal Salad Model") +
    xlim(range(full_post, dag_output)),
  
  ncol = 2)

####################
## MODEL CHECKING ##
####################

# NEEDS TO BE UPDATED

########################################################
########################################################
## NOW THE MODELS WITH SPECIES RICHNESS CONTROLED FOR ##----------------------------------------
########################################################
########################################################

################
## SST  MODEL ##
################

sst_richness_formula <- bf(log(aesthe_survey) ~ sst_mean + nb_species + 
                          abs_latitude +
                          as.factor(Temperature_Zone) +
                          (1 | Country/SiteCode),
                        
                        family=gaussian())

sst_richness_model <- brm(sst_richness_formula,
                 data=model_data,
                 chains=3, iter=3000, cores=ncores,
                 c(set_prior("normal(0,3)", class = "b"),
                   set_prior("normal(0,3)", class="Intercept")))

sst_richness_post <- as.data.frame(as.matrix(sst_richness_model)) %>%
  select('b_sst_mean')


###################
## GRAVITY MODEL ##
###################

gravity_richness_formula <- bf(log(aesthe_survey) ~ gravtot2 + nb_species +
                              HDI2017 +
                              abs_latitude + # ABSOLUTE VALUE FOR LINEARITY
                              as.factor(Temperature_Zone) +
                              (1 | Country/SiteCode),
                            family=gaussian())

gravity_richness_model <- brm(gravity_richness_formula,
                     data=model_data,
                     chains=3, iter=3000, cores=ncores,
                     c(set_prior("normal(0,3)", class = "b"),
                       set_prior("normal(0,3)", class="Intercept")))

gravity_richness_post <- as.data.frame(as.matrix(gravity_richness_model)) %>%
  select('b_gravtot2')


###############
## NPP MODEL ##
###############

NPP_richness_formula <- bf(log(aesthe_survey) ~ NPP_mean + nb_species +
                          abs_latitude +
                          BO_nitrate +
                          BO_phosphate +
                          sst_mean +
                          wave_energy +
                          as.factor(Temperature_Zone) +
                          (1 | Country/SiteCode),
                        
                        family=gaussian())

NPP_richness_model <- brm(NPP_richness_formula,
                 data=model_data,
                 chains=3, iter=3000, cores=ncores,
                 c(set_prior("normal(0,3)", class = "b"),
                   set_prior("normal(0,3)", class="Intercept")))

NPP_richness_post <- as.data.frame(as.matrix(NPP_richness_model)) %>%
  select('b_NPP_mean')

#################
## DEPTH MODEL ##
#################

depth_richness_formula <- bf(log(aesthe_survey) ~ Depth + nb_species +
                            as.factor(Temperature_Zone) +
                            (1 | Country/SiteCode),
                          
                          family=gaussian())

depth_richness_model <- brm(depth_richness_formula,
                   data=model_data,
                   chains=3, iter=3000, cores=ncores,
                   c(set_prior("normal(0,3)", class = "b"),
                     set_prior("normal(0,3)", class="Intercept")))

depth_richness_post <- as.data.frame(as.matrix(depth_richness_model)) %>%
  select('b_Depth')

#######################
## FISH DEPEND MODEL ##
#######################

fshd_richness_formula <- bf(log(aesthe_survey) ~ fshD + nb_species +
                           HDI2017 +
                           #abs_latitude +
                           as.factor(Temperature_Zone) +
                           (1 | Country/SiteCode),
                         
                         family=gaussian())

fshd_richness_model <- brm(fshd_richness_formula,
                  data=model_data,
                  chains=3, iter=3000, cores=ncores,
                  c(set_prior("normal(0,3)", class = "b"),
                    set_prior("normal(0,3)", class="Intercept")))

fshd_richness_post <- as.data.frame(as.matrix(fshd_richness_model)) %>%
  select('b_fshD')

#################
## HDI   MODEL ##
#################

HDI_richness_formula <- bf(log(aesthe_survey) ~ HDI2017 + nb_species +
                          abs_latitude + 
                          as.factor(Temperature_Zone) +
                          (1 | Country/SiteCode),
                        
                        family=gaussian())

HDI_richness_model <- brm(HDI_richness_formula,
                 data=model_data,
                 chains=3, iter=3000, cores=ncores,
                 c(set_prior("normal(0,3)", class = "b"),
                   set_prior("normal(0,3)", class="Intercept")))

HDI_richness_post <- as.data.frame(as.matrix(HDI_richness_model)) %>%
  select('b_HDI2017')

###############
## MPA MODEL ##
###############

MPA_richness_formula  <- bf(log(aesthe_survey) ~ MPA + nb_species +
                           HDI2017 +
                           fshD +
                           gravtot2 +
                           as.factor(Temperature_Zone) +
                           (1 | Country/SiteCode),
                         
                         family=gaussian())

MPA_richness_model <- brm(MPA_richness_formula,
                 data=model_data,
                 chains=3, iter=3000, cores=ncores,
                 c(set_prior("normal(0,3)", class = "b"),
                   set_prior("normal(0,3)", class="Intercept")))

MPA_richness_post <- as.data.frame(as.matrix(MPA_richness_model)) %>%
  select('b_MPANotake','b_MPARestrictedtake')

#################################
## BENTHIC COMPOSITION MODEL 1 ##
#################################

benthic_PC1_richness_formula <- 
  bf(log(aesthe_survey) ~ PC1_imputed + nb_species +
       Depth +
       abs_latitude +
       MPA + 
       NPP_mean +
       sst_mean +
       as.factor(Temperature_Zone) +
       (1 | Country/SiteCode),
     
     family=gaussian()) 

benthic_PC1_ricness_model <- brm(benthic_PC1_richness_formula,
                         data=model_data,
                         chains=3, iter=3000, cores=ncores,
                         c(set_prior("normal(0,3)", class = "b"),
                           set_prior("normal(0,3)", class="Intercept")))

benthic_PC1_richness_post <- as.data.frame(as.matrix(benthic_PC1_ricness_model))  %>%
  select('b_PC1_imputed')

#################################
## BENTHIC COMPOSITION MODEL 2 ##
#################################

benthic_PC2_richness_formula <- 
  bf(log(aesthe_survey) ~ PC2_imputed + nb_species +
       Depth +
       abs_latitude +
       MPA + 
       NPP_mean +
       sst_mean +
       as.factor(Temperature_Zone) +
       (1 | Country/SiteCode),
     
     family=gaussian()) 

benthic_PC2_richness_model <- brm(benthic_PC2_richness_formula,
                         data=model_data,
                         chains=3, iter=3000, cores=ncores,
                         c(set_prior("normal(0,3)", class = "b"),
                           set_prior("normal(0,3)", class="Intercept")))

benthic_PC2_richness_post <- as.data.frame(as.matrix(benthic_PC2_richness_model)) %>%
  select('b_PC2_imputed')

###############
## DHW MODEL ##
###############

DHW_richness_formula <- bf(log(aesthe_survey) ~ dhw_mean + nb_species +
                          abs_latitude +
                          sst_mean +
                          as.factor(Temperature_Zone) +
                          (1 | Country/SiteCode),
                        
                        family=gaussian())

DHW_richness_model <- brm(DHW_richness_formula,
                 data=model_data,
                 chains=3, iter=3000, cores=ncores,
                 c(set_prior("normal(0,3)", class = "b"),
                   set_prior("normal(0,3)", class="Intercept")))

DHW_richness_post <- as.data.frame(as.matrix(DHW_richness_model)) %>%
  select('b_dhw_mean')


###############################################
## RECOMBINE DAG MODELS AND RENAME VARIABLES ##
###############################################

dag_output_richness <- data.frame(
  sst_richness_post, gravity_richness_post, NPP_richness_post, depth_richness_post,
  fshd_richness_post, HDI_richness_post, MPA_richness_post,
  benthic_PC1_richness_post, benthic_PC2_richness_post, DHW_richness_post)

write.table(dag_output_richness, "dag_output_richness.txt")

names(dag_output_richness) <- gsub("b_", "", names(dag_output_richness))

dag_output_richness <- dag_output_richness %>%
  rename("Log Human Gravity" = gravtot2,
    "No Take MPA" = MPANotake,
    "Restricted Take MPA" = MPARestrictedtake,
    "Sea Surface Temperature" = sst_mean,
    "Net Primary Productivity" = NPP_mean,
    "Degree Heating Weeks" = dhw_mean,
    "Depth" =  Depth,
    "Fisheries Dependency" = fshD,
    "Human Dev. Index" = HDI2017,
    "Benthic Composition (PC1)" = PC1_imputed,
    "Benthic Composition (PC2)" = PC2_imputed)

dag_estimates_richness <- data.frame(median=apply(dag_output_richness, 2, median))
dag_estimates_richness$abs_effect <- abs(dag_estimates_richness$median)
dag_estimates_richness <- dag_estimates_richness %>%
  arrange(desc(abs_effect))

dag_output_richness <- dag_output_richness[,order(match(colnames(dag_output_richness), rownames(dag_estimates_richness)))]

#########################
## SIMPLE FOREST PLOTS ##
#########################

color_scheme_set("darkgray")
mcmc_areas_ridges(dag_output_richness) +
  ggtitle("X") +
  xlim(range(dag_output_richness))

#########################
## SIMPLE FOREST PLOTS ##
#########################

dag_output <- read.table("dag_output.txt")

dag_output <- dag_output[,order(match(colnames(dag_output), rownames(dag_estimates_richness)))]

ggpubr::ggarrange(
  
  mcmc_intervals(dag_output) +
    ggtitle("Total Causal Effect") +
    xlim(range(dag_output)),
  
  mcmc_intervals(dag_output_richness) +
    ggtitle("Causal Effect After Richness") +
    xlim(range(dag_output_richness)),
  
  ncol = 2)
