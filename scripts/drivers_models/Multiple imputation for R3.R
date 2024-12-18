

library(brms)
library(bayesplot)
library(parallel)
library(rstan)
library(ggplot2)
library(ggridges)
library(dplyr)

###################################################
## SCRIPT FOR MULTIPLE IMPUTATION FOR REVIEWER 3 ##
###################################################

# RUN A MISSFOREST IMPUTATION ON BENTHIC DATA (ALREADY COMPLETED IN BENTHIC IMPUTATION CODE)
# EXTRACT THE BENTHIC DATA FROM EACH IMPUTATION
# RUN A BAYESIAN MODEL FOR THAT EXTRACTION
# SAVE THE POSTERIOR DISTRIBUTION
# REPEAT 10 TIMES
# COMBINE POSTERIOR DISTRIBUTIONS 


################################################
## UPLOAD THE 10 SETS OF IMPUTED BENTHIC DATA ##
################################################

PC1_multiple <- readRDS("outputs/PC1_multiple.rds")
PC2_multiple <- readRDS("outputs/PC2_multiple.rds")

######################################
## Z SCORE THE MULTIPLE IMPUTATIONS ##
######################################

z_score_2sd <- function(x){ (x - mean(x,na.rm=T)) / (2*sd(x,na.rm=T))}

z_vars <- PC1_multiple %>%
  select_if(is.numeric) %>%
  select(-any_of(c("SurveyID","SiteLongitude","SiteLatitude","SiteCode"))) %>%
  colnames()

PC1_multiple <- PC1_multiple %>%
  mutate_if(colnames(PC1_multiple) %in% z_vars, z_score_2sd)

PC2_multiple <- PC2_multiple %>%
  mutate_if(colnames(PC2_multiple) %in% z_vars, z_score_2sd)

################################
## IMPORT PREPARED MODEL DATA ##
################################

standardized_data <-  readRDS("outputs/standardized_data.rds")

# DELETE THE BENTHIC DATA ALREADY IN THERE

standardized_data$PC1_imputed <- NULL
standardized_data$PC2_imputed <- NULL

##################################
## PARALLEL SETTINGS FOR MODELS ##
##################################

ncores = detectCores()
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

######################################
## SET GLOBAL PRIORS FOR ALL MODELS ##
######################################

prior <- c(set_prior("normal(0,3)", class = "b"),
           set_prior("normal(0,3)", class = "Intercept"))
#set_prior("exp(1)", class = "sd"))

###########################################################
## RUN A GIANT LOOP TO RUN 10 MODELS AND SAVE THE OUTPUT ##
###########################################################


PC1_posteriors <- NULL

for( i in 1:10){
  
  col_name <- names(PC1_multiple)[i+4]
  
  PC1_iteration <- PC1_multiple %>%
    select(SurveyID, col_name)
  
  colnames(PC1_iteration)[2] <- "PC1_imputed"
  
  model_data <- merge(standardized_data, PC1_iteration, by="SurveyID", all.x=TRUE)
  
  model_data$PC1_imputed[model_data$Temperature_Zone=="Temperate"] <- 0
  
  #################################
  ## BENTHIC COMPOSITION MODEL 1 ##
  #################################
  
  benthic_PC1_model_formula <- 
    bf(log(aesthe_survey_abund) ~ PC1_imputed +
         Depth + 
         gravtot2 +
         MPA + 
         NPP_mean +
         sst_mean +
         as.factor(Temperature_Zone) +
         (1 | Country/SiteCode),
       
       family=gaussian()) 
  
  benthic_PC1_model <- brm(benthic_PC1_model_formula,
                           data=model_data,
                           chains=4, iter=4000, cores=ncores,
                           prior=prior)
  
  benthic_PC1_post <- as.data.frame(as.matrix(benthic_PC1_model))  %>%
    select('b_PC1_imputed')
  
  # SAVE POSTERIOR OUTPUTS
  PC1_posteriors <- rbind(PC1_posteriors, data.frame(benthic_PC1_post, iteration=i))
  
}

saveRDS(PC1_posteriors, "outputs/PC1_posteriors.rds")

## RIDGE PLOT THE POSTERIORS FROM EACH ITERATION

ggplot(PC1_posteriors, aes(x = b_PC1_imputed, y = factor(iteration), fill = factor(iteration))) +
  geom_density_ridges(alpha = 0.5, scale = 1.5) +  # alpha controls transparency
  scale_fill_viridis_d(option = "C") +  # Use a color scale with good contrast
  theme_ridges() +  # Nice theme for ridgeline plots
  labs(x = "Values", y = "Iteration", title = "Posterior Distributions") +
  theme(legend.position = "none")  # Optionally remove the legend

############################################
## UPLOAD THE ORIGINAL BENTHIC DATA MODEL ##
############################################

original_model <- readRDS("outputs/BIG_FILES/benthic_PC1_model.rds")
original_posterior <- as.data.frame(as.matrix(original_model))  %>%
  select('b_PC1_imputed')

##########################################################################
## COMBINE THE MULTIPLE IMPUTATION POSTERIOR AND THE ORIGINAL POSTERIOR ##
##########################################################################

combined_data <- data.frame(values=c(PC1_posteriors$b_PC1_imputed, original_posterior$b_PC1_imputed))
combined_data$identity <- c(rep("multiple", nrow(PC1_posteriors)), rep("original", nrow(original_posterior)))

ggplot(combined_data, aes(x = values, y = factor(identity), fill = factor(identity))) +
  geom_density_ridges(alpha = 0.5, scale = 1.5) +  # alpha controls transparency
  scale_fill_viridis_d(option = "C") +  # Use a color scale with good contrast
  theme_ridges() +  # Nice theme for ridgeline plots
  labs(x = "Values", y = "Identity", title = "Posterior Distributions") +
  theme(legend.position = "none")  # Optionally remove the legend





###########################################################
## RUN A GIANT LOOP TO RUN 10 MODELS AND SAVE THE OUTPUT ##
###########################################################


PC2_posteriors <- NULL

for( i in 1:10){
  
  col_name <- names(PC1_multiple)[i+4]
  
  PC2_iteration <- PC2_multiple %>%
    select(SurveyID, col_name)
  
  colnames(PC2_iteration)[2] <- "PC2_imputed"
  
  model_data <- merge(standardized_data, PC2_iteration, by="SurveyID", all.x=TRUE)
  
  model_data$PC2_imputed[model_data$Temperature_Zone=="Temperate"] <- 0
  
  #################################
  ## BENTHIC COMPOSITION MODEL 2 ##
  #################################
  
  benthic_PC2_model_formula <- 
    bf(log(aesthe_survey_abund) ~ PC2_imputed +
         Depth + 
         gravtot2 +
         MPA + 
         NPP_mean +
         sst_mean +
         as.factor(Temperature_Zone) +
         (1 | Country/SiteCode),
       
       family=gaussian()) 
  
  benthic_PC2_model <- brm(benthic_PC2_model_formula,
                           data=model_data,
                           chains=4, iter=4000, cores=ncores,
                           prior=prior)
  
  benthic_PC2_post <- as.data.frame(as.matrix(benthic_PC2_model))  %>%
    select('b_PC2_imputed')
  
  # SAVE POSTERIOR OUTPUTS
  PC2_posteriors <- rbind(PC2_posteriors, data.frame(benthic_PC2_post, iteration=i))
  
}

saveRDS(PC2_posteriors, "outputs/PC2_posteriors.rds")

## RIDGE PLOT THE POSTERIORS FROM EACH ITERATION

ggplot(PC2_posteriors, aes(x = b_PC2_imputed, y = factor(iteration), fill = factor(iteration))) +
  geom_density_ridges(alpha = 0.5, scale = 1.5) +  # alpha controls transparency
  scale_fill_viridis_d(option = "C") +  # Use a color scale with good contrast
  theme_ridges() +  # Nice theme for ridgeline plots
  labs(x = "Values", y = "Iteration", title = "Posterior Distributions") +
  theme(legend.position = "none")  # Optionally remove the legend

############################################
## UPLOAD THE ORIGINAL BENTHIC DATA MODEL ##
############################################

original_model <- readRDS("outputs/BIG_FILES/benthic_PC2_model.rds")
original_posterior <- as.data.frame(as.matrix(original_model))  %>%
  select('b_PC2_imputed')

##########################################################################
## COMBINE THE MULTIPLE IMPUTATION POSTERIOR AND THE ORIGINAL POSTERIOR ##
##########################################################################

combined_data <- data.frame(values=c(PC2_posteriors$b_PC2_imputed, original_posterior$b_PC2_imputed))
combined_data$identity <- c(rep("multiple", nrow(PC2_posteriors)), rep("original", nrow(original_posterior)))

ggplot(combined_data, aes(x = values, y = factor(identity), fill = factor(identity))) +
  geom_density_ridges(alpha = 0.5, scale = 1.5) +  # alpha controls transparency
  scale_fill_viridis_d(option = "C") +  # Use a color scale with good contrast
  theme_ridges() +  # Nice theme for ridgeline plots
  labs(x = "Values", y = "Identity", title = "Posterior Distributions") +
  theme(legend.position = "none")  # Optionally remove the legend

