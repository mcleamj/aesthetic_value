
#######################################################################################
#'  CODE TO EXAMINE FAMILY LEVEL CONTRIBUTION TO HIGHER AESTHETIC VALUE
#'  IN MPA SITES, USING BAYESIAN REGRESSION MODELS
#'
#' @author Matthew McLean, \email {mcleamj@gmail.com},
#'         
#' @date Updated May, 2023
########################################################################################

######################
## LIBRARY PACKAGES ##
######################

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

####################################
## IMPORT STANDARDIZED MODEL DATA ##
####################################

standardized_data <- readRDS("outputs/standardized_data.rds")

##################################
## PARALLEL SETTINGS FOR MODELS ##
##################################

ncores = detectCores()
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

################################
# IMPORT FAMILY ABUNDANCE DATA #
################################

family_abundances <- read_rds("outputs/family_abundances.rds")

#family_abundances <- family_abundances %>%
#  rownames_to_column("SiteCode")


#'###########################################
#' FILTER TO ONLY FISHED AND NO TAKE SITES 
#' FITLER TO ONLY ECOREGIONS 
#' WITH BOTH FISHED AND NO TAKE
#'############################################

family_model_data <- standardized_data %>%
  filter(MPA != "Restricted take") 
family_model_data$MPA <- droplevels(family_model_data$MPA)

selected_regions <- family_model_data %>%
  group_by(Ecoregion, MPA) %>%
  dplyr::summarise(n()) %>%
  group_by(Ecoregion) %>%
  dplyr::filter(n()>1)

family_model_data <- family_model_data %>%
  filter(Ecoregion %in% selected_regions$Ecoregion)

##########################################
# WHICH ARE THE MOST DOMINANT FAMILIES ? #
##########################################

family_dom <- family_abundances %>%
  ungroup() %>%
  select(-SiteCode) %>%
  colMeans() %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  dplyr::rename("mean" = ".",
                "family" = rowname) %>%
  arrange(desc(mean))
family_dom$log_mean <- log(family_dom$mean)
family_dom$mean <- round(family_dom$mean, digits = 5)
family_dom$family <- gsub("family_","",family_dom$family)


# HOW MANY SITES AND ECOREGIONS PER FAMILY
# HOW MANY FAMILIES PER SITE?

family_per_site <- family_abundances %>%
  ungroup() %>%
  select(-SiteCode) %>%
  mutate_if(is.numeric, ~1 * (. > 0)) %>%
  rowSums()
family_per_site <- data.frame(SiteCode=family_abundances$SiteCode,
                              family_per_site)

site_per_family <-  family_abundances %>%
  ungroup() %>%
  select(-SiteCode) %>%
  mutate_if(is.numeric, ~1 * (. > 0)) %>%
  colSums() %>%
  as.data.frame() %>%
  dplyr::rename(n_sites = ".")
site_per_family$family <- rownames(site_per_family)

family_dom <- merge(family_dom, site_per_family, by="family")

dom_families <- family_dom %>%
  filter(mean >= quantile(family_dom$mean, probs = 0.75))
##filter(n_sites >=10)

#######################################
## PREPARE AND MERGE DATA FOR MODELS ##
#######################################

## COMBINE THE FAMILY ABUNDANCES AND MODEL DATA
## AT THE SITE LEVEL FOR MPA MODELS

MPA_standardized_data <- family_model_data %>%
  select(SiteCode, SiteLongitude, SiteLatitude, Country, Ecoregion, MPA)

MPA_standardized_data <- MPA_standardized_data[!duplicated(MPA_standardized_data$SiteCode),]  

abundance_standardized_data <- merge(MPA_standardized_data, family_abundances, by="SiteCode")

################################################
## HOW MANY SITES PER ECOREGION AND COUNTRY ? ##
################################################

site_per_eco <- family_model_data %>%
  select(SiteCode, Ecoregion) %>%
  group_by(Ecoregion) %>%
  dplyr::summarise(n())

site_per_country <- family_model_data %>%
  select(SiteCode, Country) %>%
  group_by(Country) %>%
  dplyr::summarise(n())


###########################################
## LOOP TO RUN ONE MODEL FOR EACH FAMILY 
## FAMILY ABUNDANCE IN FUNCTION OF MPA 
## WITH A RANDOM EFFECT FOR ECOREGION 
###########################################

abundance_fit_list <- as.data.frame(matrix(ncol=nrow(family_dom), nrow=8000))
colnames(abundance_fit_list) <- family_dom$family

start_time <- Sys.time()

for(i in 1:nrow(family_dom)){
  sub_model <- brm(formula=brmsformula(paste("log10(",family_dom$family[i],"+1)", 
  "~ MPA + s(SiteLongitude, SiteLatitude) + (1 | Ecoregion)")),
                   family=gaussian(),
                   data=abundance_standardized_data,
                   chains=4,
                   cores=ncores, 
                   iter=4000,
                   control = list(max_treedepth = 12))
  
  sub_posterior <- as.data.frame(as.matrix(sub_model)) %>%
    select(b_Intercept, b_MPANotake)
  
  abundance_fit_list[,i] <- sub_posterior$b_MPANotake
  
}

end_time <- Sys.time()
end_time - start_time

mcmc_intervals(abundance_fit_list)

saveRDS(abundance_fit_list,"outputs/abundance_fit_list_spatial.rds")

