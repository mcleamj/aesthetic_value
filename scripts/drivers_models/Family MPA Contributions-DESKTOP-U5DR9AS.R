
#######################################################################################
#'  CODE TO EXAMINE FAMILY LEVEL CONTRIBUTION TO HIGHER AESTHETIC VALUE
#'  IN MPA SITES, USING A BAYESIAN MODEL WITH INTERACTIONS
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
if(!require(rgdal)){install.packages("rgdal"); library(rgdal)}
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

##################################
## ADD NEW IMPUTED BENTHIC DATA ##
##################################

# benthic_PCA <- read_rds("outputs/RLS_benthic_PCA_imputed.rds")
# benthic_PCA <- benthic_PCA %>%
#   select(SurveyID, PC1_imputation, PC2_imputation) %>%
#   rename(PC1_imputed = PC1_imputation, PC2_imputed = PC2_imputation)
# 
# model_data <- merge(model_data, benthic_PCA, by = "SurveyID", all=T)

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

#####################################################
## SCALE ALL THE NUMERIC PREDICTORS TO MEAN 0 SD 2 ##
#####################################################

# z_score_2sd <- function(x){ (x - mean(x,na.rm=T)) / (2*sd(x,na.rm=T))}
# 
# z_vars <- model_data %>%
#   select_if(is.numeric) %>%
#   select(-any_of(c("SurveyID","SiteLongitude","SiteLatitude","aesthe_survey"))) %>%
#   colnames()
# 
# standardized_data <- model_data %>%
#   mutate_if(colnames(model_data) %in% z_vars, z_score_2sd)


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

# family_proportions <- functcomp(as.matrix(family_trait),
#                                as.matrix(site_sp_occ),
#                                CWM.type = "all")
# saveRDS(family_proportions, "outputs/family_proportions.rds")
family_proportions <- read_rds("outputs/family_proportions.rds")

colnames(family_proportions) <- gsub("family_", "", colnames(family_proportions))

family_proportions <- family_proportions %>%
  rownames_to_column("SiteCode")


##################################################
## IMPORT ECOREGION VARYING SLOPES EFFECT SIZES ##
##################################################

# ecoregion_model <- read_rds("outputs/BIG_FILES/ecoregion_varying_slopes_model.rds")
# 
# eco_intercept <- as.data.frame(coef(ecoregion_model)$Ecoregion[,,"Intercept"]) %>%
#   rename(Intercept = "Estimate")
# eco_intercept$ECOREGION <- rownames(eco_intercept)
# 
# eco_no_take <- as.data.frame(coef(ecoregion_model)$Ecoregion[,,"MPANotake"]) %>%
#   rename(Slope = "Estimate")
# eco_no_take$ECOREGION <- rownames(eco_no_take)
# 
# # TOP 10% OF EFFECT SIZES
# top_eco <- eco_no_take %>%
#   filter(Slope >= quantile(eco_no_take$Slope, probs = 0.90))

#########################################
## SUBSET MODEL DATA TO TOP ECOREGIONS 
## AND NO TAKE AND FISHED SITES ONLY 
#########################################

# sub_data <- standardized_data %>%
#   filter(Ecoregion %in% top_eco$ECOREGION) %>%
#   filter(MPA != "Restricted take")
# sub_data$MPA <- droplevels(sub_data$MPA)
# 
# sub_data <- merge(sub_data, family_proportions, by="SurveyID")

##############################################
## SMALL TEST MODEL WITH FAMILY INTERACTION ##
##############################################

# test_model_formula  <- bf(log(aesthe_survey) ~ 
#                             MPA:Chaetodontidae + 
#                             MPA:Pomacanthidae +
#                             MPA:Carangidae +
#                             MPA:Scombridae +
#                                  HDI2017 + # CONTROL VARIABLE
#                                  fshD + # CONTROL VARIABLE
#                                  gravtot2 + # CONTROL VARIABLE
#                                  as.factor(Temperature_Zone) + # CONTROL VARIABLE
#                                  (1 | SiteCode) + # RANDOM INTERCEPTS FOR SITES
#                                  (1 | Ecoregion),
# 
#                                family=gaussian())
# 
# test_model <- brm(test_model_formula,
#                   
#                        data=sub_data,
#                   
#                        chains=3, iter=3000, cores=ncores,
#                        c(set_prior("normal(0,3)", class = "b"),
#                          set_prior("normal(0,3)", class="Intercept")))
# 
# test_model_post <- as.data.frame(as.matrix(test_model))
# 
# conditional_effects(test_model, "Chaetodontidae:MPA")
# conditional_effects(test_model, "MPA:Pomacanthidae")
# conditional_effects(test_model, "MPA:Carangidae")
# conditional_effects(test_model, "MPA:Scombridae")


#'######################################################
#' DO I MODEL EACH FAMILY IN FUNCTION OF MPA VS FISHED? 
#'


#'###########################################
#' FILTER TO ONLY FISHED AND NO TAKE SITES 
#' COULD ALSO FITLER TO ONLY ECOREGIONS 
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

#########################################
# WHICH ARE THE MOST DOMINANT FAMILIES?
#########################################

family_abund <- family_proportions %>%
  select(-SiteCode) %>%
  colMeans() %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  dplyr::rename("mean" = ".", 
                "family" = rowname) %>%
  arrange(desc(mean))
family_abund$log_mean <- log(family_abund$mean)
family_abund$family <- gsub("family_","",family_abund$family)

dom_families <- family_abund %>%
  filter(mean >= quantile(family_abund$mean, probs = 0.75))

## ADD OR SUBTRACT 0.001 TO BOUND FOR BETA MODELS

family_proportions[family_proportions==0] <- 0.001
family_proportions[family_proportions==1] <- 0.999

## COMBINE THE FAMILY PROPORTIONS AND MODEL DATA
## AT THE SITE LEVEL FOR MPA MODELS
## AT THE SURVEY LEVEL FOR BEAUTY MODELS

MPA_model_data <- family_model_data %>%
  select(SiteCode, Country, Ecoregion, MPA)
MPA_model_data <- MPA_model_data[!duplicated(MPA_model_data$SiteCode),]  
MPA_model_data <- merge(MPA_model_data, family_proportions, by="SiteCode")

beauty_model_data <- merge(family_model_data, 
                           family_proportions, by="SiteCode")


########################################################
## MODEL FAMILY PROPORTIONS IN FUNCTION OF MPA STATUS ##
########################################################

# paste(dom_families$family, sep="", collapse = ",")
# 
# family_MPA_formula  <- bf(Scaridae 
#                             ~ MPA:log(aesthe_survey) +
#                             
#                             # (1 | SiteCode)
#                               (1 | Ecoregion),
#                           
#                           family=Beta(link = "logit", link_phi = "log"))
# 
# family_MPA_model <- brm(family_MPA_formula,
#                   
#                   data=family_model_data,
#                   
#                   chains=1, iter=3000, cores=ncores,
#                   c(set_prior("normal(0,3)", class = "b"),
#                     set_prior("normal(0,3)", class="Intercept")))
# 
# conditional_effects(family_MPA_model, "Sc")

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


##############################################################################
## LOOP TO RUN ONE MODEL FOR EACH FAMILY 
## FAMILY PROPORTION IN FUNCTION OF MPA 
## WITH A RANDOM EFFECT FOR ECOREGION (WHAT ABOUT SITE NESTED IN ECOREGION)
##############################################################################

fit_list <- as.data.frame(matrix(ncol=nrow(dom_families), nrow=4000))
colnames(fit_list) <- dom_families$family

start_time <- Sys.time()

for(i in 1:nrow(dom_families)){
  sub_model <- brm(formula=brmsformula(paste(dom_families$family[i], "~ MPA + (1 | Ecoregion)")),
                   family=Beta(link = "logit", link_phi = "log"),
                   data=MPA_model_data,
                   chains=4,
                   cores=ncores, 
                   iter=2000)
                   #control = list(adapt_delta=0.95))
  sub_posterior <- as.data.frame(as.matrix(sub_model)) %>%
    select(b_MPANotake)
  
  fit_list[,i] <- sub_posterior$b_MPANotake
  
}

end_time <- Sys.time()
end_time - start_time

mcmc_intervals(fit_list)

saveRDS(fit_list,"outputs/fit_list.rds")
fit_list <- read_rds("outputs/fit_list.rds")


####################################################
## MODEL BEAUTY IN FUNCTION OF FAMILY PROPORTIONS ##
####################################################

paste(dom_families$family, sep="", collapse = "+")
paste(family_abund$family, sep="", collapse = "+")

beauty_family_formula <- bf(log(aesthe_survey) ~ 
                              Labridae+Pomacentridae+Serranidae+Chaetodontidae+Acanthuridae+Kyphosidae+Scaridae+Sparidae+Monacanthidae+Mullidae+Blenniidae+Odacidae+Gobiidae+Pomacanthidae+Tetraodontidae+Lutjanidae+Apogonidae+Tripterygiidae+Carangidae+Plesiopidae+Latridae+Pempheridae+Balistidae+Haemulidae+Cheilodactylidae+Holocentridae+Cirrhitidae+Pinguipedidae+Diodontidae+Aplodactylidae+Sebastidae+Siganidae+Nemipteridae+Caesionidae+Dinolestidae+Enoplosidae+Lethrinidae+Scorpaenidae+Zanclidae+Aracanidae+Muraenidae+Synodontidae+Nototheniidae+Cottidae+Pseudochromidae+Ostraciidae+Embiotocidae+Gadidae+Microdesmidae+Hexagrammidae+Chironemidae+Labrisomidae+Moridae+Monodactylidae+Gerreidae+Pholidae+Mugilidae+Sphyraenidae+Pentacerotidae+Atherinidae+Platycephalidae+Trachichthyidae+Clinidae+Arripidae+Sciaenidae+Stichaeidae+Neosebastidae+Aulopidae+Ephippidae+Priacanthidae+Plotosidae+Scombridae+Malacanthidae+Bovichtidae+Latidae+Gobiesocidae+Congiopodidae+Terapontidae+Grammatidae+Berycidae+Callionymidae+Agonidae+Echeneidae+Aulorhynchidae+Glaucosomatidae+Oplegnathidae+Tetrarogidae+Moronidae+Sillaginidae+Chaenopsidae+Triglidae+Channichthyidae+Spratelloididae+Belonidae+Engraulidae+Clupeidae+Opistognathidae+Ophichthidae+Pataecidae+Parascorpididae+Anguillidae+Aploactinidae+Chanidae+Callanthiidae+Alosidae+Ammodytidae+Congridae+Dorosomatidae+Dichistiidae+Gasterosteidae+Liparidae+Salmonidae+Pomatomidae+Pholidichthyidae+Zeidae+Rhycheridae+Ambassidae+Gempylidae+Batrachoididae+Anarhichadidae+Synanceiidae+Monocentridae+Gnathanacanthidae+Hemiramphidae+Creediidae+Lotidae+Scombropidae+Trachinidae+Molidae+Ophidiidae+Ogcocephalidae+Centrolophidae+Megalopidae+Bythitidae+Percophidae+Kuhliidae+Elopidae+Centrogenyidae+Rachycentridae+
                              (1 | Ecoregion),
                          family=gaussian())

beauty_family_model <- brm(beauty_family_formula,
                    
                    data=beauty_model_data,
                    
                    chains=4, iter=4000, cores=ncores)

beauty_family_post <- data.frame(as.matrix(beauty_family_model)) %>%
  select(b_Labridae:b_Caesionidae)

mcmc_intervals(beauty_family_post)                   

conditional_effects(beauty_family_model)


#'#######################################
#' FOREST PLOT, FAMILY CONTRIBUTION TO 
#'#######################################

fit_list <- read_rds("outputs/fit_list.rds")

model_coefs <- brms::posterior_summary(fit_list, probs=c(0.10,0.25,0.75,0.90)) %>%
  as.data.frame
model_coefs <- model_coefs %>%
  arrange(rownames(model_coefs))

beauty_values <- posterior_summary(beauty_family_post) %>%
  as.data.frame() %>%
  select(Estimate) %>%
  rename(beauty=Estimate)
rownames(beauty_values) <- gsub("b_","",rownames(beauty_values))
beauty_values <- beauty_values %>%
  arrange(rownames(beauty_values))

identical(rownames(model_coefs), rownames(beauty_values))

model_coefs <- cbind(model_coefs, beauty_values) 
model_coefs <- model_coefs %>%  
  arrange(Estimate)

model_coefs$label <- rownames(model_coefs)

model_coefs$top_families <- ifelse(model_coefs$Estimate >= quantile(model_coefs$Estimate, probs=0.75) &
                                     model_coefs$beauty >= quantile(model_coefs$beauty, probs=0.75),
                                   "yes","no")
                                     
                                     

plot_colors <- variablecol(colvar = model_coefs$beauty, col = jet(n=nrow(model_coefs)), clim=range(model_coefs$beauty))

graphics.off()
par(mar=c(4,12,4,4))
scatter2D(model_coefs$Estimate, seq(1:nrow(model_coefs)), xlim=c(min(model_coefs$Q10,na.rm = TRUE),max(model_coefs$Q90,na.rm = TRUE)),
     ylim=c(min(seq(1:nrow(model_coefs))-0.25),max(seq(1:nrow(model_coefs))+0.25)),cex=0,
     xlab="MPA Effect Size", ylab=NA, yaxt = "n",
     cex.lab=1.25, colvar = model_coefs$beauty, col=jet(n=nrow(model_coefs)))
mtext(side=4, "Aesthetic Value Effect", line=1.5, cex=1.25)
title("", line=1,
      font.main=1, cex.main=1.5)

par(lend=1)
x0 <- model_coefs$Q10
x1 <- model_coefs$Q90
y0 <- seq(1:nrow(model_coefs))
y1 <- seq(1:nrow(model_coefs))
segments(x0,y0,x1,y1, lwd=2, col=plot_colors)
x0 <- model_coefs$Q25
x1 <- model_coefs$Q75
y0 <- seq(1:nrow(model_coefs))
y1 <- seq(1:nrow(model_coefs))
segments(x0,y0,x1,y1, lwd=5, col=plot_colors)

points(model_coefs$Estimate, seq(1:nrow(model_coefs)),pch=21,col=1,bg=plot_colors, 
       cex=ifelse(model_coefs$top_families=="yes",2.75,2), 
       lwd=ifelse(model_coefs$top_families=="yes",3,1))
abline(v=0, lwd=1.5, col=adjustcolor("black",alpha.f = 0.75), lty=2)

axis(2, at = seq(1:nrow(model_coefs)), 
     labels = model_coefs$label,
     las=2, cex.axis=1.)







