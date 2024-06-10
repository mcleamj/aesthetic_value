
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
if(!require(readr)){install.packages("readr"); library(readr)}
if(!require(tibble)){install.packages("tibble"); library(tibble)}
if(!require(dplyr)){install.packages("dplyr"); library(dplyr)}

################################
## IMPORT PREPARED MODEL DATA ##
################################

standardized_data <-  read_rds("outputs/standardized_data.rds")

#############################
## ADD TROPHIC COMPOSITION ##
#############################

trophic <- read_rds("outputs/trophic_composition.rds")

standardized_data <- merge(standardized_data, trophic, 
                           by="SurveyID")

##################################
## PARALLEL SETTINGS FOR MODELS ##
##################################

ncores = detectCores()
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())






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

MPA_model_formula  <- bf(log(aesthe_survey) ~ MPA +
                           HDI2017 + # CONTROL
                           fshD + # CONTROL
                           gravtot2 + # CONTROL
                           as.factor(Temperature_Zone) + # CONTROL
                           
                           (1 | Country/SiteCode),
                         
                         family=gaussian())

MPA_model <- brm(MPA_model_formula,
                 data=standardized_data,
                 chains=4, iter=4000, cores=ncores,
                 c(set_prior("normal(0,3)", class = "b"),
                   set_prior("normal(0,3)", class="Intercept")))
MPA_model <- read_rds("outputs/BIG_FILES/MPA_model.rds")

MPA_post <- as.data.frame(as.matrix(MPA_model)) %>%
  select('b_MPANotake') %>%
  dplyr::rename("MPA Total Causal Effect" = b_MPANotake)

#################################
## TAXO DIVERSITY CONTRIBUTION ##
#################################

tax_div_model_formula  <- bf(log(aesthe_survey) ~ MPA +
                               HDI2017 + # CONTROL
                               fshD + # CONTROL
                               gravtot2 + # CONTROL
                               Temperature_Zone + # CONTROL
                               
                               taxo_richness + # BLOCKED
                               
                               
                               (1 | Country/SiteCode),
                             
                             family=gaussian())

tax_div_model <- brm(tax_div_model_formula,
                     
                     data=standardized_data,
                     
                     chains=4, iter=4000, cores=ncores,
                     c(set_prior("normal(0,3)", class = "b"),
                       set_prior("normal(0,3)", class="Intercept")))

saveRDS(tax_div_model, "outputs/BIG_FILES/MPA_mechanisms_tax_div.rds")
tax_div_model <- read_rds("outputs/BIG_FILES/MPA_mechanisms_tax_div.rds")

tax_div_post <- as.data.frame(as.matrix(tax_div_model)) %>%
  select('b_MPANotake') 

tax_div_post <- data.frame(MPA_post$`MPA Total Causal Effect` - tax_div_post$b_MPANotake)
colnames(tax_div_post) <- "tax_richness"

################################
## FUN DIVERSITY CONTRIBUTION ##
################################

fun_div_model_formula  <- bf(log(aesthe_survey) ~ MPA +
                               HDI2017 + # CONTROL
                               fshD + # CONTROL
                               gravtot2 + # CONTROL
                               Temperature_Zone + # CONTROL
                               
                               fun_richness + # BLOCKED
                               
                               (1 | Country/SiteCode),
                             
                             family=gaussian())

fun_div_model <- brm(fun_div_model_formula,
                     
                     data=standardized_data,
                     
                     chains=4, iter=4000, cores=ncores,
                     c(set_prior("normal(0,3)", class = "b"),
                       set_prior("normal(0,3)", class="Intercept")))

saveRDS(fun_div_model, "outputs/BIG_FILES/MPA_mechanisms_fun_div.rds")
fun_div_model <- read_rds("outputs/BIG_FILES/MPA_mechanisms_fun_div.rds")

fun_div_post <- as.data.frame(as.matrix(fun_div_model)) %>%
  select('b_MPANotake') 

fun_div_post <- data.frame(MPA_post$`MPA Total Causal Effect` - fun_div_post$b_MPANotake)
colnames(fun_div_post) <- "fun_richness"



##################################
## PHYLO DIVERSITY CONTRIBUTION ##
##################################

phylo_div_model_formula  <- bf(log(aesthe_survey) ~ MPA +
                                 HDI2017 + # CONTROL
                                 fshD + # CONTROL
                                 gravtot2 + # CONTROL
                                 Temperature_Zone + # CONTROL
                                 
                                 phylo_richness + # BLOCKED
                                 
                                 (1 | Country/SiteCode),
                               
                               family=gaussian())

phylo_div_model <- brm(phylo_div_model_formula,
                       
                       data=standardized_data,
                       
                       chains=4, iter=4000, cores=ncores,
                       c(set_prior("normal(0,3)", class = "b"),
                         set_prior("normal(0,3)", class="Intercept")))

saveRDS(phylo_div_model, "outputs/BIG_FILES/MPA_mechanisms_phylo_div.rds")
phylo_div_model <- read_rds("outputs/BIG_FILES/MPA_mechanisms_phylo_div.rds")

phylo_div_post <- as.data.frame(as.matrix(phylo_div_model)) %>%
  select('b_MPANotake')

phylo_div_post <- data.frame(MPA_post$`MPA Total Causal Effect` - phylo_div_post$b_MPANotake)
colnames(phylo_div_post) <- "phylo_richness"


##############################
## BENTHIC PC1 CONTRIBUTION ##
##############################

benthic_PC1_model_formula  <- bf(log(aesthe_survey) ~ MPA +
                                   HDI2017 + # CONTROL
                                   fshD + # CONTROL
                                   gravtot2 + # CONTROL
                                   Temperature_Zone + # CONTROL
                                   
                                   PC1_imputed + # BLOCKED
                                   
                                   (1 | Country/SiteCode),
                                 
                                 family=gaussian())

benthic_PC1_model <- brm(benthic_PC1_model_formula,
                         
                         data=standardized_data,
                         
                         chains=4, iter=4000, cores=ncores,
                         c(set_prior("normal(0,3)", class = "b"),
                           set_prior("normal(0,3)", class="Intercept")))

saveRDS(benthic_PC1_model, "outputs/BIG_FILES/MPA_mechanisms_benthic_1.rds")
benthic_PC1_model <- read_rds("outputs/BIG_FILES/MPA_mechanisms_benthic_1.rds")

benthic_PC1_post <- as.data.frame(as.matrix(benthic_PC1_model)) %>%
  select('b_MPANotake') 

benthic_PC1_post <- data.frame(MPA_post$`MPA Total Causal Effect` - benthic_PC1_post$b_MPANotake)
colnames(benthic_PC1_post) <- "PC1_Benthic"


##############################
## BENTHIC PC2 CONTRIBUTION ##
##############################

benthic_PC2_model_formula  <- bf(log(aesthe_survey) ~ MPA +
                                   HDI2017 + # CONTROL
                                   fshD + # CONTROL
                                   gravtot2 + # CONTROL
                                   Temperature_Zone + # CONTROL
                                   
                                   PC2_imputed + # BLOCKED
                                   
                                   (1 | Country/SiteCode),
                                 
                                 family=gaussian())

benthic_PC2_model <- brm(benthic_PC2_model_formula,
                         
                         data=standardized_data,
                         
                         chains=4, iter=4000, cores=ncores,
                         c(set_prior("normal(0,3)", class = "b"),
                           set_prior("normal(0,3)", class="Intercept")))

saveRDS(benthic_PC2_model, "outputs/BIG_FILES/MPA_mechanisms_benthic_2.rds")
benthic_PC2_model <- read_rds("outputs/BIG_FILES/MPA_mechanisms_benthic_2.rds")

benthic_PC2_post <- as.data.frame(as.matrix(benthic_PC2_model)) %>%
  select('b_MPANotake') 

benthic_PC2_post <- data.frame(MPA_post$`MPA Total Causal Effect` - benthic_PC2_post$b_MPANotake)
colnames(benthic_PC2_post) <- "PC2_Benthic"


##############################
## TROPHIC PC1 CONTRIBUTION ##
##############################

# WHAT ABOUT DIRECTIONALITY OF PC AXES? - GOOD

trophic_PC1_model_formula  <- bf(log(aesthe_survey) ~ MPA +
                                   HDI2017 + # CONTROL
                                   fshD + # CONTROL
                                   gravtot2 + # CONTROL
                                   Temperature_Zone + # CONTROL
                                   
                                   PC1_trophic + # BLOCKED
                                   
                                   (1 | Country/SiteCode),
                                 
                                 family=gaussian())

trophic_PC1_model <- brm(trophic_PC1_model_formula,
                         
                         data=standardized_data,
                         
                         chains=4, iter=4000, cores=ncores,
                         c(set_prior("normal(0,3)", class = "b"),
                           set_prior("normal(0,3)", class="Intercept")))

saveRDS(trophic_PC1_model, "outputs/BIG_FILES/MPA_mechanisms_trophic_1.rds")
trophic_PC1_model <- read_rds("outputs/BIG_FILES/MPA_mechanisms_trophic_1.rds")

trophic_PC1_post <- as.data.frame(as.matrix(trophic_PC1_model)) %>%
  select('b_MPANotake')

trophic_PC1_post <- data.frame(MPA_post$`MPA Total Causal Effect` - trophic_PC1_post$b_MPANotake)
colnames(trophic_PC1_post) <- "PC1_Trophic"


##############################
## TROPHIC PC2 CONTRIBUTION ##
##############################

trophic_PC2_model_formula  <- bf(log(aesthe_survey) ~ MPA +
                                   HDI2017 + # CONTROL
                                   fshD + # CONTROL
                                   gravtot2 + # CONTROL
                                   Temperature_Zone + # CONTROL
                                   
                                   PC2_trophic + # BLOCKED
                                   
                                   (1 | Country/SiteCode),
                                 
                                 family=gaussian())

trophic_PC2_model <- brm(trophic_PC2_model_formula,
                         
                         data=standardized_data,
                         
                         chains=4, iter=4000, cores=ncores,
                         c(set_prior("normal(0,3)", class = "b"),
                           set_prior("normal(0,3)", class="Intercept")))

saveRDS(trophic_PC2_model, "outputs/BIG_FILES/MPA_mechanisms_trophic_2.rds")
trophic_PC2_model <- read_rds("outputs/BIG_FILES/MPA_mechanisms_trophic_2.rds")

trophic_PC2_post <- as.data.frame(as.matrix(trophic_PC2_model)) %>%
  select('b_MPANotake')

trophic_PC2_post <- data.frame(MPA_post$`MPA Total Causal Effect` - trophic_PC2_post$b_MPANotake)
colnames(trophic_PC2_post) <- "PC2_Trophic"


###########################
## Taxo PC1 CONTRIBUTION ##
###########################

# WHAT ABOUT DIRECTIONALITY OF PC AXES? - GOOD

Taxo_PC1_model_formula  <- bf(log(aesthe_survey) ~ MPA +
                                HDI2017 + # CONTROL
                                fshD + # CONTROL
                                gravtot2 + # CONTROL
                                Temperature_Zone + # CONTROL
                                
                                Taxo_PC1 + # BLOCKED
                                
                                (1 | Country/SiteCode),
                              
                              family=gaussian())

Taxo_PC1_model <- brm(Taxo_PC1_model_formula,
                      
                      data=standardized_data,
                      
                      chains=4, iter=4000, cores=ncores,
                      c(set_prior("normal(0,3)", class = "b"),
                        set_prior("normal(0,3)", class="Intercept")))

saveRDS(Taxo_PC1_model, "outputs/BIG_FILES/MPA_mechanisms_taxonomic_1.rds")
Taxo_PC1_model <- read_rds("outputs/BIG_FILES/MPA_mechanisms_taxonomic_1.rds")

Taxo_PC1_post <- as.data.frame(as.matrix(Taxo_PC1_model)) %>%
  select('b_MPANotake')

Taxo_PC1_post <- data.frame(MPA_post$`MPA Total Causal Effect` - Taxo_PC1_post$b_MPANotake)
colnames(Taxo_PC1_post) <- "Taxo_PC1"


##############################
## Taxo PC2 CONTRIBUTION ##
##############################

Taxo_PC2_model_formula  <- bf(log(aesthe_survey) ~ MPA +
                                HDI2017 + # CONTROL
                                fshD + # CONTROL
                                gravtot2 + # CONTROL
                                Temperature_Zone + # CONTROL
                                
                                Taxo_PC2 + # BLOCKED
                                
                                (1 | Country/SiteCode),
                              
                              family=gaussian())

Taxo_PC2_model <- brm(Taxo_PC2_model_formula,
                      
                      data=standardized_data,
                      
                      chains=4, iter=4000, cores=ncores,
                      c(set_prior("normal(0,3)", class = "b"),
                        set_prior("normal(0,3)", class="Intercept")))

saveRDS(Taxo_PC2_model, "outputs/BIG_FILES/MPA_mechanisms_taxonomic_2.rds")
Taxo_PC2_model <- read_rds("outputs/BIG_FILES/MPA_mechanisms_taxonomic_2.rds")

Taxo_PC2_post <- as.data.frame(as.matrix(Taxo_PC2_model)) %>%
  select('b_MPANotake')

Taxo_PC2_post <- data.frame(MPA_post$`MPA Total Causal Effect` - Taxo_PC2_post$b_MPANotake)
colnames(Taxo_PC2_post) <- "Taxo_PC2"

###########################
## GROUP MODEL ESTIMATES ##
###########################

model_outputs <- data.frame(MPA_post, tax_div_post, fun_div_post, phylo_div_post,
                            benthic_PC1_post, benthic_PC2_post,
                            trophic_PC1_post, trophic_PC2_post,
                            Taxo_PC1_post, Taxo_PC2_post)

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





