
#######################################################################################
#'  BAYSESIAN HIERARCHICAL MODEL OF THE EFFECT OF MPA ON AESTHETIC VALUE
#'  WITH VARYING SLOPES FOR COUNTRY OR ECOREGION
#'  TO ASSESS WHETHER THE EFFECT OF MPA VARIES ACROSS SPACE (LATITUDE)
#'
#' @author Matthew McLean, \email {mcleamj@gmail.com},
#'         
#' @date JUNE 9, 2023
########################################################################################

######################
## LIBRARY PACKAGES ##
######################

if(!require(data.table)){install.packages("data.table"); library(data.table)}
if(!require(mapproj)){install.packages("mapproj"); library(mapproj)}
if(!require(plot3D)){install.packages("plot3D"); library(plot3D)}
if(!require(ggplot2)){install.packages("ggplot2"); library(ggplot2)}
if(!require(pals)){install.packages("pals"); library(pals)}
if(!require(brms)){install.packages("brms"); library(brms)}
if(!require(bayesplot)){install.packages("bayesplot"); library(bayesplot)}
if(!require(parallel)){install.packages("parallel"); library(parallel)}
if(!require(rstan)){install.packages("rstan"); library(rstan)}
if(!require(bayestestR)){install.packages("bayestestR"); library(bayestestR)}
if(!require(performance)){install.packages("performance"); library(performance)}
if(!require(readr)){install.packages("readr"); library(readr)}
if(!require(tibble)){install.packages("tibble"); library(tibble)}
if(!require(stringr)){install.packages("stringr"); library(stringr)}
if(!require(dplyr)){install.packages("dplyr"); library(dplyr)}

#####################################
## IMPORT DATA "STANDARDIZED DATA" ##
#####################################

standardized_data <- read_rds("outputs/standardized_data.rds")

##################################
## PARALLEL SETTINGS FOR MODELS ##
##################################

ncores = detectCores()
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#########################
## VARYING SLOPE MODEL ##
#########################

ecoregion_model_formula  <- bf(log(aesthe_survey) ~ MPA + # TEST VARIABLE
                           HDI2017 + # CONTROL VARIABLE
                           fshD + # CONTROL VARIABLE
                           gravtot2 + # CONTROL VARIABLE
                           as.factor(Temperature_Zone) + # CONTROL VARIABLE
                           (1 | Country/SiteCode) + # RANDOM INTERCEPTS FOR SITES
                           (1 + MPA | Ecoregion), # RANDOM INTERCEPT AND SLOPE FOR ECOREGION
                         
                         family=gaussian())

ecoregion_model <- brm(ecoregion_model_formula,
                 data=standardized_data,
                 chains=4, iter=4000, cores=ncores,
                 c(set_prior("normal(0,3)", class = "b"),
                   set_prior("normal(0,3)", class="Intercept")))

saveRDS(ecoregion_model, "outputs/BIG_FILES/ecoregion_varying_slopes_model.rds")
ecoregion_model <- read_rds("outputs/BIG_FILES/ecoregion_varying_slopes_model.rds")

eco_model_post <- as.data.frame(as.matrix(ecoregion_model)) %>%
  #select('b_MPANotake','b_MPARestrictedtake')
  select('b_MPANotake':'b_gravtot2')

mcmc_intervals(eco_model_post) # QUICK LOOK AT MODEL OUTPUT

################################
## EXTRACT MODEL COEFFICIENTS ##
################################

eco_intercept <- as.data.frame(coef(ecoregion_model)$Ecoregion[,,"Intercept"]) %>%
  rename(Intercept = "Estimate")
eco_intercept$ECOREGION <- rownames(eco_intercept)

eco_no_take <- as.data.frame(coef(ecoregion_model)$Ecoregion[,,"MPANotake"]) %>%
  rename(Slope = "Estimate")
eco_no_take$ECOREGION <- rownames(eco_no_take)

# HISTOGRAM OF EFFECT SIZES
hist(eco_no_take$Slope)
abline(v=0, lty=2, lwd=2)

########################
## ADD ECOREGION INFO ##
########################

ecoregion_info <- model_data %>%
  select(Ecoregion, SiteLongitude, SiteLatitude, PC1_imputed, PC2_imputed, nb_species, sst_mean) %>%
  group_by(Ecoregion) %>%
  summarise_all(.funs=mean) %>%
  rename(ECOREGION = Ecoregion)

###################################
## IS THERE A LATITUDINAL BIAS IN 
## MPA/FISHED SITE OCCURENCE?
###################################

MPA_by_region <- model_data %>%
  select(SiteCode, Ecoregion, MPA) %>%
  group_by(Ecoregion, MPA) %>%
  summarise(n()) %>%
  rename(ECOREGION = Ecoregion, count = "n()")
MPA_by_region <- merge(MPA_by_region, ecoregion_info[,c("ECOREGION","SiteLongitude","SiteLatitude","nb_species")],
                       by="ECOREGION")

ggplot(MPA_by_region, 
       aes(x = SiteLatitude, y = count)) + 
  geom_point() + 
  facet_grid(~MPA)

ggplot(MPA_by_region, 
       aes(x = nb_species, y = count)) + 
  geom_point() + 
  facet_grid(~MPA)

# NO OBVIOUS BIAS BY LATITUDE OR SPECIES RICHNESS

# MAKE A DATA FRAME TO FOCUS ONLY ON NO TAKE MPAS
# ONLY KEEP ECOREGIONS WITH BOTH FISHED AND NO-TAKE SITES
no_take_only <- model_data %>%
  select(SiteCode, SiteLongitude, SiteLatitude, 
         Country, Ecoregion, MPA) %>%
  filter(MPA != "Restricted take")
no_take_only$MPA <- droplevels(no_take_only$MPA)
no_take_only <- no_take_only[!duplicated(no_take_only$SiteCode),]

selected_regions <- no_take_only %>%
  select(SiteCode, Ecoregion, MPA) %>%
  group_by(Ecoregion, MPA) %>%
  summarise(n()) %>%
  group_by(Ecoregion) %>%
  filter(n()>1)

length(unique(model_data$Ecoregion))
length(unique(selected_regions$Ecoregion)) 
# ONLY 43/96 REGIONS HAVE BOTH TYPES OF SITES

graphics.off()
ecoregion_info$MPA_and_fished <- ecoregion_info$ECOREGION %in% selected_regions$Ecoregion
boxplot(ecoregion_info$SiteLatitude ~ ecoregion_info$MPA_and_fished)

#' LATITUDINAL GRADIENT IN MPA EFFECT IS NOT DUE TO
#' THERE ONLY BEING ECOREGIONS WITH BOTH FISHED AND NO TAKE SITES
#' NEAR THE EQUATOR, IN FACT THEY ARE MORE PREVELANT IN THE SOUTHERN HEMISPHERE

