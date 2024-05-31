
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

#install.packages("data/worms", repos = NULL, type="source")
library(worms)

if(!require(data.table)){install.packages("data.table"); library(data.table)}
if(!require(mapproj)){install.packages("mapproj"); library(mapproj)}
if(!require(plot3D)){install.packages("plot3D"); library(plot3D)}
if(!require(pals)){install.packages("pals"); library(pals)}
if(!require(brms)){install.packages("brms"); library(brms)}
if(!require(FD)){install.packages("FD"); library(FD)}
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

#####################################
## IMPORT AND ATTACH AESTHETIC DATA #
#####################################

aaesthe_surveyetic_data <- read.csv("outputs/survey_aesth.csv")

model_data <- merge(model_data, aaesthe_surveyetic_data, by="SurveyID")

##################################
## PARALLEL SETTINGS FOR MODELS ##
##################################

ncores = detectCores()
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#################################
# IMPORT FAMILY PROPORTION DATA #
#################################

family_proportions <- read_rds("outputs/family_proportions.rds")

#family_proportions <- family_proportions %>%
#  rownames_to_column("SiteCode")

#############################################
## FAMILY COUNTS PER SITE (N SP IN FAMILY) ##
#############################################

family_count <- data.frame(SiteCode=rownames(site_sp_occ), site_sp_occ)

family_count <- tidyr::gather(family_count,
                              scientificname,
                              presence,
                              "Abalistes.stellatus":"Zoramia.viridiventer")

family_count$scientificname <- gsub("\\."," ", family_count$scientificname)

family_count <- merge(family_count, family_info[,c("scientificname","family")], by="scientificname")
head(family_count)

family_count <- family_count %>%
  select(SiteCode, family, presence) %>%
  group_by(SiteCode, family) %>%
  summarise_all(.funs=sum)

family_count <- tidyr::spread(family_count,
                      family,
                      presence)

#'###########################################
#' FILTER TO ONLY FISHED AND NO TAKE SITES 
#' FITLER TO ONLY ECOREGIONS 
#' WITH BOTH FISHED AND NO TAKE
#'############################################

family_model_data <- model_data %>%
  filter(MPA != "Restricted take") 
family_model_data$MPA <- droplevels(family_model_data$MPA)

selected_regions <- family_model_data %>%
  group_by(Ecoregion, MPA) %>%
  dplyr::summarise(n()) %>%
  group_by(Ecoregion) %>%
  dplyr::filter(n()>1)

family_model_data <- family_model_data %>%
  filter(Ecoregion %in% selected_regions$Ecoregion)

##########################################################
# WHICH ARE THE MOST DOMINANT FAMILIES (BY PROPORTIION)?
##########################################################

family_abund <- family_proportions %>%
  select(-SiteCode) %>%
  colMeans() %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  dplyr::rename("mean" = ".",
                "family" = rowname) %>%
  arrange(desc(mean))
family_abund$log_mean <- log(family_abund$mean)
family_abund$mean <- round(family_abund$mean, digits = 5)
family_abund$family <- gsub("family_","",family_abund$family)


# HOW MANY SITES AND ECOREGIONS PER FAMILY
# HOW MANY FAMILIES PER SITE?

family_per_site <- family_proportions %>%
  select(-SiteCode) %>%
  mutate_if(is.numeric, ~1 * (. > 0)) %>%
  rowSums()
family_per_site <- data.frame(SiteCode=family_proportions$SiteCode,
                              family_per_site)

site_per_family <-  family_proportions %>%
  select(-SiteCode) %>%
  mutate_if(is.numeric, ~1 * (. > 0)) %>%
  colSums() %>%
  as.data.frame() %>%
  dplyr::rename(n_sites = ".")
site_per_family$family <- rownames(site_per_family)

family_abund <- merge(family_abund, site_per_family, by="family")

dom_families <- family_abund %>%
  filter(mean >= quantile(family_abund$mean, probs = 0.75))
  ##filter(n_sites >=10)

#######################################
## PREPARE AND MERGE DATA FOR MODELS ##
#######################################

## ADD OR SUBTRACT 0.001 TO BOUND FOR BETA MODELS

family_proportions[family_proportions==0] <- 0.001
family_proportions[family_proportions==1] <- 0.999

## COMBINE THE FAMILY PROPORTIONS AND MODEL DATA
## AT THE SITE LEVEL FOR MPA MODELS

MPA_model_data <- family_model_data %>%
  select(SiteCode, Country, Ecoregion, MPA)

MPA_model_data <- MPA_model_data[!duplicated(MPA_model_data$SiteCode),]  

proportion_model_data <- merge(MPA_model_data, family_proportions, by="SiteCode")

count_model_data <- merge(MPA_model_data, family_count, by="SiteCode")


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
## FAMILY PROPORTION IN FUNCTION OF MPA 
## WITH A RANDOM EFFECT FOR ECOREGION 
###########################################

proportion_fit_list <- as.data.frame(matrix(ncol=nrow(family_abund), nrow=8000))
colnames(proportion_fit_list) <- family_abund$family

start_time <- Sys.time()

for(i in 1:nrow(family_abund)){
  sub_model <- brm(formula=brmsformula(paste(family_abund$family[i], "~ MPA + (1 | Ecoregion)")),
                   family=Beta(link = "logit", link_phi = "log"),
                   data=proportion_model_data,
                   chains=4,
                   cores=ncores, 
                   iter=4000)
                   #control = list(adapt_delta=0.95))
  sub_posterior <- as.data.frame(as.matrix(sub_model)) %>%
    select(b_Intercept, b_MPANotake)
  
  proportion_fit_list[,i] <- sub_posterior$b_MPANotake
  
}

end_time <- Sys.time()
end_time - start_time

mcmc_intervals(proportion_fit_list)

saveRDS(proportion_fit_list,"outputs/proportion_fit_list.rds")
proportion_fit_list <- read_rds("outputs/proportion_fit_list.rds")



###########################################
## LOOP TO RUN ONE MODEL FOR EACH FAMILY 
## FAMILY COUNT IN FUNCTION OF MPA 
## WITH A RANDOM EFFECT FOR ECOREGION 
###########################################

count_fit_list <- as.data.frame(matrix(ncol=nrow(family_abund), nrow=8000))
colnames(count_fit_list) <- family_abund$family

start_time <- Sys.time()

for(i in 1:nrow(family_abund)){
  sub_model <- brm(formula=brmsformula(paste(family_abund$family[i], "~ MPA + (1 | Ecoregion)")),
                   family=negbinomial(link="log", link_shape = "log"),
                   data=count_model_data,
                   chains=4,
                   cores=ncores, 
                   iter=4000)
  #control = list(adapt_delta=0.95))
  sub_posterior <- as.data.frame(as.matrix(sub_model)) %>%
    select(b_Intercept, b_MPANotake)
  
  count_fit_list[,i] <- sub_posterior$b_MPANotake
  
}

end_time <- Sys.time()
end_time - start_time

mcmc_intervals(count_fit_list)

saveRDS(count_fit_list,"outputs/count_fit_list.rds")
count_fit_list <- read_rds("outputs/count_fit_list.rds")

