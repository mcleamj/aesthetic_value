
#######################################################################################
#'  CODE TO EXAMINE FAMILY LEVEL CONTRIBUTION TO HIGHER AESTHETIC VALUE
#'  IN MPA SITES, USING BAYESIAN REGRESSION MODELS
#'
#' @author Matthew McLean, \email {mcleamj@gmail.com},
#'         
#' @date JUNE 9, 2023
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

##'##################################
##' IMPORT OCCURENCE INFORMATION 
##' AGGREGATE OCCURENCES TO SITE LEVEL
##' AND CALCULATE FAMILY PROPORTIONS 
##' FOR EACH SITE
##' #################################

survey_sp_occ <- readr::read_rds("outputs/sp_pres_matrix.rds")

survey_species <- survey_sp_occ %>%
  dplyr::select(where(is.numeric)) %>%
  colnames() %>%
  as.data.frame() %>%
  dplyr::rename("species"= ".")

survey_species$species <- gsub("_", " ", survey_species$species) 
survey_species <- survey_species %>%
  rename("species_name" = "species")

survey_species$scientificname <- gsub(" spp.", "", survey_species$species_name)

#family_info <- worms::wormsbynames(survey_species$scientificname) 
#family_info <- merge(family_info, survey_species, by="scientificname")
#saveRDS(family_info, "outputs/family_info.rds")
family_info <- read_rds("outputs/family_info.rds")

## CALCULATE FAMILY PROPORTIONS AT EACH SITE

site_sp_occ <- merge(model_data[,c("SiteCode","SurveyID")], survey_sp_occ,
                     by="SurveyID")
site_sp_occ <- site_sp_occ %>%
  select(-SurveyID) %>%
  group_by(SiteCode) %>%
  summarise_all(.funs = mean, na.rm=T) %>%
  mutate_if(is.numeric, ~1 * (. > 0))

family_trait <- family_info %>%
  select(species_name, family) %>%
  arrange(species_name) %>%
  column_to_rownames("species_name")

site_sp_occ <- site_sp_occ %>%
  column_to_rownames("SiteCode")
colnames(site_sp_occ) <- gsub("_", " ", colnames(site_sp_occ) )
site_sp_occ <- site_sp_occ %>%
 select(order(colnames(site_sp_occ)))

identical(colnames(site_sp_occ), rownames(family_trait))

#family_proportions <- functcomp(as.matrix(family_trait),
#                               as.matrix(site_sp_occ),
#                               CWM.type = "all")
#saveRDS(family_proportions, "outputs/family_proportions.rds")

family_proportions <- read_rds("outputs/family_proportions.rds")

colnames(family_proportions) <- gsub("family_", "", colnames(family_proportions))

family_proportions <- family_proportions %>%
  rownames_to_column("SiteCode")

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
  rename(n_sites = ".")
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
  summarise(n())

site_per_country <- family_model_data %>%
  select(SiteCode, Country) %>%
  group_by(Country) %>%
  summarise(n())


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


##############################################
## WHICH FAMILIES HAVE HIGHEST MPA EFFECTS? ##
##############################################

family_MPA <- data.frame(family=colnames(proportion_fit_list),
                            MPA_effect_mean = colMeans(proportion_fit_list)) %>%
  arrange(desc(MPA_effect_mean))


##################################################
## CALCULATE AVERAGE AESTHETIC VALUE PER FAMILY ##
##################################################

aesthe_species <- read.csv2(here::here("data", "aesthe_langlois_2022.csv"))
aesthe_species$sp_name <- as.character(gsub("_"," ",aesthe_species$sp_name))
aesthe_species$aesthe_score <- as.numeric(aesthe_species$aesthe_score)

aesthe_species <- aesthe_species %>% rename(scientificname = sp_name)

aesthe_species <- merge(aesthe_species, family_info[,c("scientificname", "family")],
                        by="scientificname")
family_beauty <- aesthe_species %>%
  select(-scientificname) %>%
  group_by(family) %>%
  summarise_all(.funs=mean)












#########################################################
## FOREST PLOT OF MPA EFFECT COLORED BY AVERAGE BEAUTY ##
#########################################################

MPA_effect_summary <- as.data.frame(brms::posterior_summary(fit_list,
                                              probs=c(0.10,0.25,0.75,0.90)))

MPA_effect_summary$family <- rownames(MPA_effect_summary)

MPA_effect_summary <- merge(MPA_effect_summary, family_beauty, by="family")
MPA_effect_summary <- MPA_effect_summary %>%
  arrange((Estimate))
MPA_effect_summary$log_beauty <- log(MPA_effect_summary$aesthe_score)

# TAKE ONLY THE FAMILIES WITH HIGH MPA EFFECT (TOP X%)
top_MPA_families <- MPA_effect_summary %>%
  filter(Estimate >= quantile(MPA_effect_summary$Estimate, prob=0.75))

plot_colors <- variablecol(colvar = top_MPA_families$aesthe_score, col = jet(n=nrow(top_MPA_families)), clim=range(family_beauty$aesthe_score))

graphics.off()
par(mar=c(4,12,4,4))
scatter2D(top_MPA_families$Estimate, seq(1:nrow(top_MPA_families)), xlim=c(min(top_MPA_families$Q10,na.rm = TRUE),max(top_MPA_families$Q90,na.rm = TRUE)),
          ylim=c(min(seq(1:nrow(top_MPA_families))-0.25),max(seq(1:nrow(top_MPA_families))+0.25)),cex=0,
          xlab="MPA Effect Size", ylab=NA, yaxt = "n",
          cex.lab=1.25, colvar = top_MPA_families$aesthe_score, col=jet(n=nrow(top_MPA_families)))
mtext(side=4, "Average Aesthetic Value", line=2, cex=1.1)
title("", line=1,
      font.main=1, cex.main=1.5)

par(lend=1)
x0 <- top_MPA_families$Q10
x1 <- top_MPA_families$Q90
y0 <- seq(1:nrow(top_MPA_families))
y1 <- seq(1:nrow(top_MPA_families))
segments(x0,y0,x1,y1, lwd=2, col=plot_colors)
x0 <- top_MPA_families$Q25
x1 <- top_MPA_families$Q75
y0 <- seq(1:nrow(top_MPA_families))
y1 <- seq(1:nrow(top_MPA_families))
segments(x0,y0,x1,y1, lwd=5, col=plot_colors)

points(top_MPA_families$Estimate, seq(1:nrow(top_MPA_families)),pch=21,col=1,bg=plot_colors, cex=1.5)
abline(v=0, lwd=1.5, col=adjustcolor("black",alpha.f = 0.75), lty=2)

axis(2, at = seq(1:nrow(top_MPA_families)), 
     labels = top_MPA_families$family,
     las=2, cex.axis=1)

