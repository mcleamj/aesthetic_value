#######################################################################################
#'  BAYSESIAN HIERARCHICAL MODELS TO ASSESS DRIVERS OF AESTHETIC VALUE
#'  BASED ON DIRECTED ACYCLIC GRAPHS (DAG)
#'  THE CORRESPONDING DAG CAN BE FOUND AT:
#'  
#'  http://dagitty.net/m6WdXviAT
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
if(!require(ggridges)){install.packages("ggridges"); library(ggridges)}
if(!require(ggplot2)){install.packages("ggplot2"); library(ggplot2)}
if(!require(plot3D)){install.packages("plot3D"); library(plot3D)}
if(!require(pals)){install.packages("pals"); library(pals)}
if(!require(brms)){install.packages("brms"); library(brms)}
if(!require(bayesplot)){install.packages("bayesplot"); library(bayesplot)}
if(!require(parallel)){install.packages("parallel"); library(parallel)}
if(!require(rstan)){install.packages("rstan"); library(rstan)}
if(!require(bayestestR)){install.packages("bayestestR"); library(bayestestR)}
if(!require(readr)){install.packages("readr"); library(readr)}
if(!require(performance)){install.packages("performance"); library(performance)}
if(!require(glmmTMB)){install.packages("glmmTMB"); library(glmmTMB)}
if(!require(tibble)){install.packages("tibble"); library(tibble)}
if(!require(dplyr)){install.packages("dplyr"); library(dplyr)}

################################
## IMPORT PREPARED MODEL DATA ##
################################

standardized_data <-  read_rds("outputs/standardized_data.rds")

##################################
## PARALLEL SETTINGS FOR MODELS ##
##################################

ncores = detectCores()
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

######################################
## SET GLOBAL PRIORS FOR ALL MODELS ##
## ###################################

prior <- c(set_prior("normal(0,3)", class = "b"),
           set_prior("normal(0,3)", class = "Intercept"),
           set_prior("exp)", class = "sd"))

#'#######################################################
#' NOW RUN THE INDIVIDUAL MODELS FOR EACH PREDICTOR
#' ACCORDING TO THE DAG (USING MINIMAL ADJUSTMENT SET)
#'#######################################################

################
## SST  MODEL ##
################

sst_model_formula <- bf(log(aesthe_survey) ~ sst_mean +
                          abs_latitude +
                          as.factor(Temperature_Zone) +
                          (1 | Country/SiteCode),
  
                        family=gaussian())

sst_model <- brm(sst_model_formula,
                 data=standardized_data,
                 prior=prior)

saveRDS(sst_model, "outputs/BIG_FILES/sst_model.rds")
sst_model <- read_rds("outputs/BIG_FILES/sst_model.rds")

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
                     data=standardized_data,
                     chains=4, iter=4000, cores=ncores,
                     prior=prior)

saveRDS(gravity_model, "outputs/BIG_FILES/gravity_model.rds")
gravity_model <- read_rds("outputs/BIG_FILES/gravity_model.rds")

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
                 data=standardized_data,
                 chains=4, iter=4000, cores=ncores,
                 prior=prior)

saveRDS(NPP_model, "outputs/BIG_FILES/NPP_model.rds")
NPP_model <- read_rds("outputs/BIG_FILES/NPP_model.rds")

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
                   data=standardized_data,
                   chains=4, iter=4000, cores=ncores,
                   prior=prior)

saveRDS(depth_model, "outputs/BIG_FILES/depth_model.rds")
depth_model <- read_rds("outputs/BIG_FILES/depth_model.rds")

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
                  data=standardized_data,
                  chains=4, iter=4000, cores=ncores,
                  prior=prior)

saveRDS(fshd_model, "outputs/BIG_FILES/fshd_model.rds")
fshd_model <- read_rds("outputs/BIG_FILES/fshd_model.rds")

fshd_post <- as.data.frame(as.matrix(fshd_model)) %>%
  select('b_fshD')

#################
## HDI   MODEL ##
#################

HDI_model_formula <- bf(log(aesthe_survey) ~ HDI2017 +
                          abs_latitude + 
                          sst_mean +
                          as.factor(Temperature_Zone) +
                          (1 | Country/SiteCode),
                          
                        family=gaussian())

HDI_model <- brm(HDI_model_formula,
                 data=standardized_data,
                 chains=4, iter=4000, cores=ncores,
                 prior=prior)

saveRDS(HDI_model, "outputs/BIG_FILES/HDI_model.rds")
HDI_model <- read_rds("outputs/BIG_FILES/HDI_model.rds")

HDI_post <- as.data.frame(as.matrix(HDI_model)) %>%
  select('b_HDI2017')

###############
## MPA MODEL ##
###############

MPA_model_formula  <- bf(log(aesthe_survey) ~ MPA +
                           HDI2017 +
                           gravtot2 +
                           as.factor(Temperature_Zone) +
                           (1 | Country/SiteCode),
                           
                         family=gaussian())

MPA_model <- brm(MPA_model_formula,
                 data=standardized_data,
                 chains=4, iter=4000, cores=ncores,
                 prior=prior)

saveRDS(MPA_model, "outputs/BIG_FILES/MPA_model.rds")
MPA_model <- read_rds("outputs/BIG_FILES/MPA_model.rds")

MPA_post <- as.data.frame(as.matrix(MPA_model)) %>%
  select('b_MPANotake','b_MPARestrictedtake')

#################################
## BENTHIC COMPOSITION MODEL 1 ##
#################################

benthic_PC1_model_formula <- 
  bf(log(aesthe_survey) ~ PC1_imputed +
       Depth + 
       gravtot2 +
       MPA + 
       NPP_mean +
       sst_mean +
       as.factor(Temperature_Zone) +
       (1 | Country/SiteCode),
       
     family=gaussian()) 

benthic_PC1_model <- brm(benthic_PC1_model_formula,
                         data=standardized_data,
                         chains=4, iter=4000, cores=ncores,
                         prior=prior)

saveRDS(benthic_PC1_model, "outputs/BIG_FILES/benthic_PC1_model.rds")
benthic_PC1_model <- read_rds("outputs/BIG_FILES/benthic_PC1_model.rds")

benthic_PC1_post <- as.data.frame(as.matrix(benthic_PC1_model))  %>%
  select('b_PC1_imputed')

#################################
## BENTHIC COMPOSITION MODEL 2 ##
#################################

benthic_PC2_model_formula <- 
  bf(log(aesthe_survey) ~ PC2_imputed +
       Depth + 
       gravtot2 +
       MPA + 
       NPP_mean +
       sst_mean +
       as.factor(Temperature_Zone) +
       (1 | Country/SiteCode),
       
     family=gaussian()) 

benthic_PC2_model <- brm(benthic_PC2_model_formula,
                         data=standardized_data,
                         chains=4, iter=4000, cores=ncores,
                         prior=prior)

saveRDS(benthic_PC2_model, "outputs/BIG_FILES/benthic_PC2_model.rds")
benthic_PC2_model <- read_rds("outputs/BIG_FILES/benthic_PC2_model.rds")

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
                 data=standardized_data,
                 chains=4, iter=4000, cores=ncores,
                 prior=prior)

saveRDS(DHW_model, "outputs/BIG_FILES/DHW_model.rds")
DHW_model <- read_rds("outputs/BIG_FILES/DHW_model.rds")

DHW_post <- as.data.frame(as.matrix(DHW_model)) %>%
  select('b_dhw_mean')

###############################################
## RECOMBINE DAG MODELS AND RENAME VARIABLES ##
###############################################

dag_output <- data.frame(sst_post, gravity_post, NPP_post, depth_post,
                         fshd_post, HDI_post, MPA_post,
                         benthic_PC1_post, benthic_PC2_post, DHW_post)

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

saveRDS(dag_output, "outputs/dag_output.rds")

dag_output <- read_rds("outputs/dag_output.rds")

dag_estimates <- data.frame(median=apply(dag_output, 2, median))
dag_estimates$abs_effect <- abs(dag_estimates$median)
dag_estimates <- dag_estimates %>%
  arrange(desc(abs_effect))

####################################
## RANK DAG ORDER BY EFFECT SIZE? ##
####################################

dag_output <- dag_output[,order(match(colnames(dag_output), rownames(dag_estimates)))]

#########################
## SIMPLE FOREST PLOTS ##
#########################

mcmc_areas_ridges(dag_output) +
  ggtitle("DAG-Based Models") +
  xlim(range(dag_output))

mcmc_intervals(dag_output) +
  ggtitle("DAG-Based Models") +
  xlim(range(dag_output))

#########################
## MODEL SUMMARY TABLE ##
#########################

posterior_summary(dag_output, probs=c(0.10,0.90))

####################
## MODEL CHECKING ##
####################

# READ MODELS IN TO START HERE
sst_model <- read_rds("outputs/BIG_FILES/sst_model.rds")
gravity_model <- read_rds("outputs/BIG_FILES/gravity_model.rds")
NPP_model <- read_rds("outputs/BIG_FILES/NPP_model.rds")
depth_model <- read_rds("outputs/BIG_FILES/depth_model.rds")
fshd_model <- read_rds("outputs/BIG_FILES/fshd_model.rds")
HDI_model <- read_rds("outputs/BIG_FILES/HDI_model.rds")
MPA_model <- read_rds("outputs/BIG_FILES/MPA_model.rds")
benthic_PC1_model <- read_rds("outputs/BIG_FILES/benthic_PC1_model.rds")
benthic_PC2_model <- read_rds("outputs/BIG_FILES/benthic_PC2_model.rds")
DHW_model <- read_rds("outputs/BIG_FILES/DHW_model.rds")

dev.size()
dev.size(units = "in")
dev.size(units = "px")
dev.size(units= "cm")

dev.off()

# SST MODEL
performance::check_model(sst_model, panel = TRUE) 
  #, check = "vif") # GOOD, POSTERIOR COULD BE BETTER, VIF HIGH DUE TO LATITUDE AND SST
#graphics.off()
#plot(sst_model, variable="b_sst_mean")
#r2_bayes(sst_model)

ggsave("figures_tables/model_checks/sst_model_check.png",
       dpi = 300,
       height = 10, width= 8, 
       units = "in",
       device = "png"
) 


# NPP MODEL
performance::check_model(NPP_model)
# GOOD, POSTERIOR COULD BE BETTER, VIF HIGH DUE TO LATIUDE AND SST
#graphics.off()
#plot(NPP_model,variable="b_NPP_mean")
#r2_bayes(NPP_model)

ggsave("figures_tables/model_checks/NPP_model_check.png",
       dpi = 300,
       height = 10, width= 8, 
       units = "in",
       device = "png"
) 

# DEPTH MODEL
performance::check_model(depth_model) # GOOD
#graphics.off()
#plot(depth_model, variable="b_Depth")
#r2_bayes(depth_model)

ggsave("figures_tables/model_checks/depth_model_check.png",
       dpi = 300,
       height = 10, width= 8, 
       units = "in",
       device = "png"
)

# GRAVITY MODEL
performance::check_model(gravity_model) # GOOD
#graphics.off()
#plot(gravity_model, variable="b_gravtot2")
#r2_bayes(gravity_model)

ggsave("figures_tables/model_checks/gravity_model_check.png",
       dpi = 300,
       height = 10, width= 8, 
       units = "in",
       device = "png"
)

# MPA MODEL
performance::check_model(MPA_model) # GOOD
#graphics.off()
#plot(MPA_model, variable=c("b_MPANotake","b_MPARestrictedtake"))
#r2_bayes(MPA_model)

ggsave("figures_tables/model_checks/MPA_model_check.png",
       dpi = 300,
       height = 10, width= 8, 
       units = "in",
       device = "png"
)

# BENTHIC MODEL 1
performance::check_model(benthic_PC1_model) # GOOD, HIGH VIF FROM LATITUDE AND SST
#graphics.off()
#plot(benthic_PC1_model,variable="b_PC1_imputed")
#r2_bayes(benthic_PC1_model)

ggsave("figures_tables/model_checks/benthic_pc1_model_check.png",
       dpi = 300,
       height = 10, width= 8, 
       units = "in",
       device = "png"
)

# BENTHIC MODEL 2
performance::check_model(benthic_PC2_model) # GOOD, HIGH VIF FROM LAT AND SST
#graphics.off()
#plot(benthic_PC2_model, variable="b_PC2_imputed")
#r2_bayes(benthic_PC2_model)

ggsave("figures_tables/model_checks/benthic_pc2_model_check.png",
       dpi = 300,
       height = 10, width= 8, 
       units = "in",
       device = "png"
)

# DHW MODEL
performance::check_model(DHW_model) # GOOD, HIGH VIF FROM LAT AND SST
#graphics.off()
#plot(DHW_model, variable="b_dhw_mean")
#r2_bayes(DHW_model)

ggsave("figures_tables/model_checks/DHW_model_check.png",
       dpi = 300,
       height = 10, width= 8, 
       units = "in",
       device = "png"
)

# HDI MODEL
performance::check_model(HDI_model) # GOOD, HIGH VIF FROM LAT AND SST
#graphics.off()
#plot(HDI_model,variable="b_HDI2017")
#r2_bayes(HDI_model)

ggsave("figures_tables/model_checks/HDI_model_check.png",
       dpi = 300,
       height = 10, width= 8, 
       units = "in",
       device = "png"
)

# FISHERIES DEPENDENCY MODEL
performance::check_model(fshd_model) # GOOD
#graphics.off()
#plot(fshd_model, variable="b_fshD")
#r2_bayes(fshd_model)

ggsave("figures_tables/model_checks/fshd_model_model_check.png",
       dpi = 300,
       height = 10, width= 8, 
       units = "in",
       device = "png"
)


##############################################
## SENSITIVITY TESTS FOR MODELS WITH HIGH VIF
## THESE MODELS HAVE HIGH VIF
## DUE TO COLLINEARITY WITH LATITUDE
## CHECK MODEL RESULT WITH LATITUDE REMOVED
## SST, NPP, DHW, HDI
##############################################

sst_sens_formula <- bf(log(aesthe_survey) ~ sst_mean +
                         as.factor(Temperature_Zone) +
                         (1 | Country/SiteCode),
                       
                       family=gaussian())

sst_sense_model <- brm(sst_sens_formula,
                       data=standardized_data,
                       chains=4, iter=4000, cores=ncores,
                       prior=prior)

saveRDS(sst_sense_model, "outputs/BIG_FILES/sst_sense_model.rds")
sst_sense_model <- read_rds("outputs/BIG_FILES/sst_sense_model.rds")

sst_sens_post <- as.data.frame(as.matrix(sst_sense_model)) %>%
  select('b_sst_mean')

performance::check_model(sst_sense_model, check="vif")

# RIDGE PLOT THE TWO POSTERIORS

sst_ridge <- data.frame(draws=c(sst_post$b_sst_mean, sst_sens_post$b_sst_mean),
                        model=rep(c("original","sensitivity"),each=nrow(sst_post)))

ggplot(sst_ridge, aes(x = draws, y = model, fill=model)) +
  geom_density_ridges() +
  ggtitle("SST Effect Sizes")


###########################
## NPP SENSITIVITY MODEL ##
###########################

NPP_sens_formula <- bf(log(aesthe_survey) ~ NPP_mean +
                         BO_nitrate +
                         BO_phosphate +
                         sst_mean +
                         wave_energy +
                         as.factor(Temperature_Zone) +
                         (1 | Country/SiteCode),
                       
                       family=gaussian())

NPP_sens_model <- brm(NPP_sens_formula,
                      data=standardized_data,
                      chains=4, iter=4000, cores=ncores,
                      prior=prior)

saveRDS(NPP_sens_model, "outputs/BIG_FILES/NPP_sens_model.rds")
NPP_sens_model <- read_rds("outputs/BIG_FILES/NPP_sens_model.rds")

NPP_sens_post <- as.data.frame(as.matrix(NPP_sens_model)) %>%
  select('b_NPP_mean')

performance::check_model(NPP_sens_model, check="vif")

# RIDGE PLOT THE TWO POSTERIORS

NPP_ridge <- data.frame(draws=c(NPP_post$b_NPP_mean, NPP_sens_post$b_NPP_mean),
                        model=rep(c("original","sensitivity"),each=nrow(NPP_post)))

ggplot(NPP_ridge, aes(x = draws, y = model, fill=model)) +
  geom_density_ridges() +
  ggtitle("NPP Effect Sizes")


###########################
## DHW SENSITIVITY MODEL ##
###########################

DHW_sens_formula <- bf(log(aesthe_survey) ~ dhw_mean +
                         sst_mean +
                         as.factor(Temperature_Zone) +
                         (1 | Country/SiteCode),
                       
                       family=gaussian())

DHW_sens_model <- brm(DHW_sens_formula,
                      data=standardized_data,
                      chains=4, iter=4000, cores=ncores,
                      prior=prior)

saveRDS(DHW_sens_model, "outputs/BIG_FILES/DHW_sens_model.rds")
DHW_sens_model <- read_rds("outputs/BIG_FILES/DHW_sens_model.rds")

DHW_sens_post <- as.data.frame(as.matrix(DHW_sens_model)) %>%
  select('b_dhw_mean')

performance::check_model(DHW_sens_model, check="vif")

# RIDGE PLOT THE TWO POSTERIORS

DHW_ridge <- data.frame(draws=c(DHW_post$b_dhw_mean, DHW_sens_post$b_dhw_mean),
                        model=rep(c("original","sensitivity"),each=nrow(DHW_post)))

ggplot(DHW_ridge, aes(x = draws, y = model, fill=model)) +
  geom_density_ridges() +
  ggtitle("DHW Effect Sizes")


###########################
## HDI SENSITIVITY MODEL ##
###########################

HDI_sens_formula <- bf(log(aesthe_survey) ~ HDI2017 +
                         sst_mean +
                         as.factor(Temperature_Zone) +
                         (1 | Country/SiteCode),
                       
                       family=gaussian())

HDI_sens_model <- brm(HDI_sens_formula,
                      data=standardized_data,
                      chains=4, iter=4000, cores=ncores,
                      prior=prior)

saveRDS(HDI_sens_model, "outputs/BIG_FILES/HDI_sens_model.rds")
HDI_sens_model <- read_rds("outputs/BIG_FILES/HDI_sens_model.rds")

HDI_sens_post <- as.data.frame(as.matrix(HDI_sens_model)) %>%
  select('b_HDI2017')

performance::check_model(HDI_sens_model, check="vif")

# RIDGE PLOT THE TWO POSTERIORS

HDI_ridge <- data.frame(draws=c(HDI_post$b_HDI2017, HDI_sens_post$b_HDI2017),
                        model=rep(c("original","sensitivity"),each=nrow(HDI_post)))

ggplot(HDI_ridge, aes(x = draws, y = model, fill=model)) +
  geom_density_ridges() +
  ggtitle("HDI Effect Sizes")



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
                 data=standardized_data,
                 chains=4, iter=4000, cores=ncores,
                 prior=prior)

saveRDS(sst_richness_model, "outputs/BIG_FILES/sst_richness_model.rds")
sst_richness_model <- read_rds("outputs/BIG_FILES/sst_richness_model.rds")

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
                     data=standardized_data,
                     chains=4, iter=4000, cores=ncores,
                     prior=prior)

saveRDS(gravity_richness_model, "outputs/BIG_FILES/gravity_richness_model.rds")
gravity_richness_model <- read_rds("outputs/BIG_FILES/gravity_richness_model.rds")

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
                 data=standardized_data,
                 chains=4, iter=4000, cores=ncores,
                 prior=prior)

saveRDS(NPP_richness_model, "outputs/BIG_FILES/NPP_richness_model.rds")
NPP_richness_model <- read_rds("outputs/BIG_FILES/NPP_richness_model.rds")

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
                   data=standardized_data,
                   chains=4, iter=4000, cores=ncores,
                   prior=prior)

saveRDS(depth_richness_model, "outputs/BIG_FILES/depth_richness_model.rds")
depth_richness_model <- read_rds("outputs/BIG_FILES/depth_richness_model.rds")

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
                  data=standardized_data,
                  chains=4, iter=4000, cores=ncores,
                  prior=prior)

saveRDS(fshd_richness_model, "outputs/BIG_FILES/fshd_richness_model.rds")
fshd_richness_model <- read_rds("outputs/BIG_FILES/fshd_richness_model.rds")

fshd_richness_post <- as.data.frame(as.matrix(fshd_richness_model)) %>%
  select('b_fshD')

#################
## HDI   MODEL ##
#################

HDI_richness_formula <- bf(log(aesthe_survey) ~ HDI2017 + nb_species +
                          abs_latitude + 
                            sst_mean +
                          as.factor(Temperature_Zone) +
                          (1 | Country/SiteCode),
                        
                        family=gaussian())

HDI_richness_model <- brm(HDI_richness_formula,
                 data=standardized_data,
                 chains=4, iter=4000, cores=ncores,
                 prior=prior)

saveRDS(HDI_richness_model, "outputs/BIG_FILES/HDI_richness_model.rds")
HDI_richness_model <- read_rds("outputs/BIG_FILES/HDI_richness_model.rds")

HDI_richness_post <- as.data.frame(as.matrix(HDI_richness_model)) %>%
  select('b_HDI2017')

###############
## MPA MODEL ##
###############

MPA_richness_formula  <- bf(log(aesthe_survey) ~ MPA + nb_species +
                           HDI2017 +
                           gravtot2 +
                           as.factor(Temperature_Zone) +
                           (1 | Country/SiteCode),
                         
                         family=gaussian())

MPA_richness_model <- brm(MPA_richness_formula,
                 data=standardized_data,
                 chains=4, iter=4000, cores=ncores,
                 prior=prior)

saveRDS(MPA_richness_model, "outputs/BIG_FILES/MPA_richness_model.rds")
MPA_richness_model <- read_rds("outputs/BIG_FILES/MPA_richness_model.rds")

MPA_richness_post <- as.data.frame(as.matrix(MPA_richness_model)) %>%
  select('b_MPANotake','b_MPARestrictedtake')

#################################
## BENTHIC COMPOSITION MODEL 1 ##
#################################

benthic_PC1_richness_formula <- 
  bf(log(aesthe_survey) ~ PC1_imputed + nb_species +
       Depth + 
       gravtot2 +
       MPA + 
       NPP_mean +
       sst_mean +
       as.factor(Temperature_Zone) +
       (1 | Country/SiteCode),
     
     family=gaussian()) 

benthic_PC1_richness_model <- brm(benthic_PC1_richness_formula,
                         data=standardized_data,
                         chains=4, iter=4000, cores=ncores,
                         prior=prior)

saveRDS(benthic_PC1_richness_model, "outputs/BIG_FILES/benthic_PC1_richness_model.rds")
benthic_PC1_richness_model <- read_rds("outputs/BIG_FILES/benthic_PC1_richness_model.rds")

benthic_PC1_richness_post <- as.data.frame(as.matrix(benthic_PC1_richness_model))  %>%
  select('b_PC1_imputed')

#################################
## BENTHIC COMPOSITION MODEL 2 ##
#################################

benthic_PC2_richness_formula <- 
  bf(log(aesthe_survey) ~ PC2_imputed + nb_species +
       Depth + 
       gravtot2 +
       MPA + 
       NPP_mean +
       sst_mean +
       as.factor(Temperature_Zone) +
       (1 | Country/SiteCode),
     
     family=gaussian()) 

benthic_PC2_richness_model <- brm(benthic_PC2_richness_formula,
                         data=standardized_data,
                         chains=4, iter=4000, cores=ncores,
                         prior=prior)

saveRDS(benthic_PC2_richness_model, "outputs/BIG_FILES/benthic_PC2_richness_model.rds")
benthic_PC2_richness_model <- read_rds("outputs/BIG_FILES/benthic_PC2_richness_model.rds")

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
                 data=standardized_data,
                 chains=4, iter=4000, cores=ncores,
                 prior=prior)

saveRDS(DHW_richness_model, "outputs/BIG_FILES/DHW_richness_model.rds")
DHW_richness_model <- read_rds("outputs/BIG_FILES/DHW_richness_model.rds")

DHW_richness_post <- as.data.frame(as.matrix(DHW_richness_model)) %>%
  select('b_dhw_mean')


###############################################
## RECOMBINE DAG MODELS AND RENAME VARIABLES ##
###############################################

dag_output_richness <- data.frame(
  sst_richness_post, gravity_richness_post, NPP_richness_post, depth_richness_post,
  fshd_richness_post, HDI_richness_post, MPA_richness_post,
  benthic_PC1_richness_post, benthic_PC2_richness_post, DHW_richness_post)

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

saveRDS(dag_output_richness, "outputs/dag_output_richness.rds")
dag_output_richness <- read_rds("outputs/dag_output_richness.rds")

dag_estimates_richness <- data.frame(median=apply(dag_output_richness, 2, median))
dag_estimates_richness$abs_effect <- abs(dag_estimates_richness$median)
dag_estimates_richness <- dag_estimates_richness %>%
  arrange(desc(abs_effect))

dag_output_richness <- dag_output_richness[,order(match(colnames(dag_output_richness), rownames(dag_estimates_richness)))]

## MATCH ORDER TO AESTH VALUE MODELS?
dag_output_richness <- dag_output_richness[,order(match(colnames(dag_output_richness), rownames(dag_estimates)))]

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

ggpubr::ggarrange(
  
  mcmc_intervals(dag_output) +
    ggtitle("Total Causal Effect") +
    xlim(range(dag_output)),
  
  mcmc_intervals(dag_output_richness) +
    ggtitle("Causal Effect After Richness") +
    xlim(range(dag_output_richness)),
  
  ncol = 2)

###########################################
## SIMPLE FOREST PLOTS SAME X AXIS SCALE ##
###########################################

ggpubr::ggarrange(
  
  mcmc_intervals(dag_output) +
    ggtitle("Total Causal Effect") +
    xlim(range(dag_output)),
  
  mcmc_intervals(dag_output_richness) +
    ggtitle("Causal Effect After Richness") +
    xlim(range(dag_output)),
  
  ncol = 2)

####################
## MODEL CHECKING ##
####################

# SST MODEL
performance::check_model(sst_richness_model) # 
# graphics.off()
# plot(sst_richness_model)
# r2_bayes(sst_richness_model)

# NPP MODEL
performance::check_model(NPP_richness__model) # 
# graphics.off()
# plot(NPP_richness_model)
# r2_bayes(NPP_richness_model)

# DEPTH MODEL
performance::check_model(depth_richness_model) # 
# graphics.off()
# plot(depth_richness_model)
# r2_bayes(depth_richness_model)

# GRAVITY MODEL
performance::check_model(gravity_richness_model) # 
# graphics.off()
# plot(gravity_richness_model)
# r2_bayes(gravity_richness_model)

# MPA MODEL
performance::check_model(MPA_richness_model) # 
# graphics.off()
# plot(MPA_richness_model)
# r2_bayes(MPA_richness_model)

# BENTHIC MODEL 1
performance::check_model(benthic_PC1_richness_model) # 
# graphics.off()
# plot(benthic_PC1_richness_model)
# r2_bayes(benthic_PC1_richness_model)

# BENTHIC MODEL 2
performance::check_model(benthic_PC2_richness_model) # 
# graphics.off()
# plot(benthic_PC2_richness_model)
# r2_bayes(benthic_PC2_richness_model)

# DHW MODEL
performance::check_model(DHW_richness_model) # 
# graphics.off()
# plot(DHW_richness_model)
# r2_bayes(DHW_richness_model)

# HDI MODEL
performance::check_model(HDI_richness_model) # 
# graphics.off()
# plot(HDI_richness_model)
# r2_bayes(HDI_richness_model)

# FISHERIES DEPENDENCY MODEL
performance::check_model(fshd_richness_model) # 
# graphics.off()
# plot(fshd_richness_model)
# r2_bayes(fshd_richness_model)

