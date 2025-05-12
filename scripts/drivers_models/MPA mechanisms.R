
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

if(!require(brms)){install.packages("brms"); library(brms)}
if(!require(bayesplot)){install.packages("bayesplot"); library(bayesplot)}
if(!require(parallel)){install.packages("parallel"); library(parallel)}
if(!require(rstan)){install.packages("rstan"); library(rstan)}
library(cmdstanr)
if(!require(readr)){install.packages("readr"); library(readr)}
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

#################################
## SET CMDSTANR AS THE BACKEND ##
#################################

options(brms.backend = "cmdstanr")

###################################################################
#' TO ASSESS THE EFFECT OF EACH PATHWAY OF THE
#' TOTAL CAUSAL EFFECT OF MPA, PERFORM THE FOLLOWING
#' USE THE MODEL FORMULA TO ESTIAMTE THE TOTAL CAUSAL EFFECT OF MPA
#' AND BLOCK EACH PATHWAY ONE BY ONE BY ADDING EACH PATHWAY
#' INDIVIDUALLY TO THE MODEL
#' THE CONTRIBUTION FOR EACH VARIABLE SHOULD BE THE AMOUNT
#' BY WHICH THE TOTAL EFFECT DECREASES WHEN THE PATHWAY IS BLOCKED 
#' THE PATHWAYS OF MPA EFFECT ARE:
#' TAXO DIVERSITY, FUN DIVERSITY, PHYLO DIVERSITY,
#' BENTHIC COMPOSITION, BIOMASS, TROPHIC STRUCTURE
####################################################################


###################################################
## MPA TOTAL CAUSAL EFFECT FROM DAG-BASED MODELS ##
###################################################

MPA_model_formula  <- bf(log(aesthe_survey_abund) ~ MPA +
                           HDI2017 + # CONTROL
                           gravtot2 + # CONTROL
                           as.factor(Temperature_Zone) + # CONTROL
                         
                           (1 | Country/SiteCode),
                         
                         family=gaussian())

MPA_model <- brm(MPA_model_formula,
                 data=standardized_data,
                 control = list(max_treedepth = 15),
                 chains=4, iter=4000, cores=ncores,
                 threads = threading(8),
                 c(set_prior("normal(0,3)", class = "b"),
                   set_prior("normal(0,3)", class="Intercept")))

MPA_model <- read_rds("outputs/BIG_FILES/MPA_model.rds")

MPA_post <- as.data.frame(as.matrix(MPA_model)) %>%
  select('b_MPANotake') %>%
  dplyr::rename("MPA Total Causal Effect" = b_MPANotake)

#################################
## TAXO DIVERSITY CONTRIBUTION ##
#################################

tax_div_formula  <- bf(log(aesthe_survey_abund) ~ MPA +
                               HDI2017 + # CONTROL
                               gravtot2 + # CONTROL
                               Temperature_Zone + # CONTROL
                               
                               taxo_entropy + # BLOCKED
                       
                               (1 | Country/SiteCode),
                             
                             family=gaussian())

tax_div_model <- brm(tax_div_formula,
                     
                     data=standardized_data,
                     
                     control = list(max_treedepth = 15),
                     chains=4, iter=4000, cores=ncores,
                     threads = threading(8),
                     c(set_prior("normal(0,3)", class = "b"),
                       set_prior("normal(0,3)", class="Intercept")))

saveRDS(tax_div_model, "outputs/BIG_FILES/MPA_mechanisms_tax_div.rds")
tax_div_model <- read_rds("outputs/BIG_FILES/MPA_mechanisms_tax_div.rds")

tax_div_post <- as.data.frame(as.matrix(tax_div_model)) %>%
  select('b_MPANotake') 

tax_div_post <- data.frame(MPA_post$`MPA Total Causal Effect` - tax_div_post$b_MPANotake)
colnames(tax_div_post) <- "tax_entropy"

################################
## FUN DIVERSITY CONTRIBUTION ##
################################

fun_div_formula  <- bf(log(aesthe_survey_abund) ~ MPA +
                               HDI2017 + # CONTROL
                               gravtot2 + # CONTROL
                               Temperature_Zone + # CONTROL
                               
                               fun_entropy + # BLOCKED
                         
                               (1 | Country/SiteCode),
                             
                             family=gaussian())

fun_div_model <- brm(fun_div_formula,
                     
                     data=standardized_data,
                     
                     control = list(max_treedepth = 15),
                     chains=4, iter=4000, cores=ncores,
                     threads = threading(8),
                     c(set_prior("normal(0,3)", class = "b"),
                       set_prior("normal(0,3)", class="Intercept")))

saveRDS(fun_div_model, "outputs/BIG_FILES/MPA_mechanisms_fun_div.rds")
fun_div_model <- read_rds("outputs/BIG_FILES/MPA_mechanisms_fun_div.rds")

fun_div_post <- as.data.frame(as.matrix(fun_div_model)) %>%
  select('b_MPANotake') 

fun_div_post <- data.frame(MPA_post$`MPA Total Causal Effect` - fun_div_post$b_MPANotake)
colnames(fun_div_post) <- "fun_entropy"



##################################
## PHYLO DIVERSITY CONTRIBUTION ##
##################################

phylo_div_formula  <- bf(log(aesthe_survey_abund) ~ MPA +
                                 HDI2017 + # CONTROL
                                 gravtot2 + # CONTROL
                                 Temperature_Zone + # CONTROL
   
                                 phylo_entropy + # BLOCKED
                                 
                                 (1 | Country/SiteCode),
                               
                               family=gaussian())

phylo_div_model <- brm(phylo_div_formula,
                       
                       data=standardized_data,
                       
                       control = list(max_treedepth = 15),
                       chains=4, iter=4000, cores=ncores,
                       threads = threading(8),
                       c(set_prior("normal(0,3)", class = "b"),
                         set_prior("normal(0,3)", class="Intercept")))

saveRDS(phylo_div_model, "outputs/BIG_FILES/MPA_mechanisms_phylo_div.rds")
phylo_div_model <- read_rds("outputs/BIG_FILES/MPA_mechanisms_phylo_div.rds")

phylo_div_post <- as.data.frame(as.matrix(phylo_div_model)) %>%
  select('b_MPANotake')

phylo_div_post <- data.frame(MPA_post$`MPA Total Causal Effect` - phylo_div_post$b_MPANotake)
colnames(phylo_div_post) <- "phylo_entropy"


##########################
## BENTHIC CONTRIBUTION ##
##########################

benthic_formula  <- bf(log(aesthe_survey_abund) ~ MPA +
                                   HDI2017 + # CONTROL
                                   gravtot2 + # CONTROL
                                   Temperature_Zone + # CONTROL
                                   
                                   PC1_imputed + # BLOCKED
                                   PC2_imputed + # BLOCKED
                                   PC3_imputed + # BLOCKED  
                                   PC4_imputed + # BLOCKED
                         
                                   (1 | Country/SiteCode),
                                 
                                 family=gaussian())

benthic_model <- brm(benthic_formula,
                         
                         data=standardized_data,
                         
                     control = list(max_treedepth = 15),
                     chains=4, iter=4000, cores=ncores,
                     threads = threading(8),
                         c(set_prior("normal(0,3)", class = "b"),
                           set_prior("normal(0,3)", class="Intercept")))

saveRDS(benthic_model, "outputs/BIG_FILES/MPA_mechanisms_benthic.rds")
benthic_model <- read_rds("outputs/BIG_FILES/MPA_mechanisms_benthic.rds")

benthic_post <- as.data.frame(as.matrix(benthic_model)) %>%
  select('b_MPANotake') 

benthic_post <- data.frame(MPA_post$`MPA Total Causal Effect` - benthic_post$b_MPANotake)
colnames(benthic_post) <- "Benthic_Composition"


##########################
## TROPHIC CONTRIBUTION ##
##########################

# WHAT ABOUT DIRECTIONALITY OF PC AXES? DOES NOT MATTER, ITS JUST VARIATION EXPLAINED

trophic_formula  <- bf(log(aesthe_survey_abund) ~ MPA +
                                   HDI2017 + # CONTROL
                                   gravtot2 + # CONTROL
                                   Temperature_Zone + # CONTROL
                                   
                                   PC1_trophic + # BLOCKED
                                   PC2_trophic + # BLOCKED
                                   PC3_trophic + # BLOCKED
                                   PC4_trophic + # BLOCKED
                                   
                       (1 | Country/SiteCode),
                                 
                                 family=gaussian())

trophic_model <- brm(trophic_formula,
                         
                         data=standardized_data,
                         
                     control = list(max_treedepth = 15),
                     chains=4, iter=4000, cores=ncores,
                     threads = threading(8),
                         c(set_prior("normal(0,3)", class = "b"),
                           set_prior("normal(0,3)", class="Intercept")))

saveRDS(trophic_model, "outputs/BIG_FILES/MPA_mechanisms_trophic.rds")
trophic_model <- read_rds("outputs/BIG_FILES/MPA_mechanisms_trophic.rds")

trophic_post <- as.data.frame(as.matrix(trophic_model)) %>%
  select('b_MPANotake')
trophic_post
trophic_post <- data.frame(MPA_post$`MPA Total Causal Effect` - trophic_post$b_MPANotake)
colnames(trophic_post) <- "Trophic_Composition"


########################################
## TAXONOMIC COMPOSITION CONTRIBUTION ##
########################################

taxo_comp_formula  <- bf(log(aesthe_survey_abund) ~ MPA +
                                HDI2017 + # CONTROL
                                gravtot2 + # CONTROL
                                Temperature_Zone + # CONTROL
                                
                                Taxo_PC1 + # BLOCKED
                                Taxo_PC2 + # BLOCKED
                                Taxo_PC3 + # BLOCKED
                                Taxo_PC4 + # BLOCKED
                         
                         (1 | Country/SiteCode),
                              
                              family=gaussian())

taxo_comp_model <- brm(taxo_comp_formula,
                      
                      data=standardized_data,
                     
                       control = list(max_treedepth = 15),
                      chains=4, iter=4000, cores=ncores,
                      threads = threading(8),
                      c(set_prior("normal(0,3)", class = "b"),
                        set_prior("normal(0,3)", class="Intercept")))

saveRDS(taxo_comp_model, "outputs/BIG_FILES/MPA_mechanisms_taxonomic.rds")
taxo_comp_model <- read_rds("outputs/BIG_FILES/MPA_mechanisms_taxonomic.rds")

taxo_comp_post <- as.data.frame(as.matrix(taxo_comp_model)) %>%
  select('b_MPANotake')

taxo_comp_post <- data.frame(MPA_post$`MPA Total Causal Effect` - taxo_comp_post$b_MPANotake)
colnames(taxo_comp_post) <- "Taxonomic_Composition"


###########################
## GROUP MODEL ESTIMATES ##
###########################

model_outputs <- data.frame(MPA_post, tax_div_post, fun_div_post, phylo_div_post,
                            benthic_post, trophic_post, taxo_comp_post)

model_estimates <- data.frame(median=apply(model_outputs, 2, median))
model_estimates$variable <- rownames(model_estimates)
model_estimates$abs_effect <- abs(model_estimates$median)
model_estimates <- model_estimates %>%
  arrange(desc(abs_effect))
model_outputs <- model_outputs[,order(match(colnames(model_outputs), model_estimates$variable))]

saveRDS(model_outputs, "outputs/MPA_mechanism_effects.rds")
model_outputs <- read_rds("outputs/MPA_mechanism_effects.rds")

mcmc_intervals(model_outputs)

model_estimates %>%
  filter(variable != "MPA.Total.Causal.Effect") %>%
  summarise(sum((median)))

model_estimates %>%
  filter(variable == "MPA.Total.Causal.Effect") %>%
  summarise(sum((median)))





