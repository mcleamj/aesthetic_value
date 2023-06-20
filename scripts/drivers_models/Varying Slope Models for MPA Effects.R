
#######################################################################################
#'  BAYSESIAN HIERARCHICAL MODEL OF THE EFFECT OF MPA ON AESTHETIC VALUE
#'  WITH VARYING SLOPES FOR COUNTRY OR ECOREGION
#'  TO ASSESS WHETHER THE EFFECT OF MPA VARIES ACROSS SPACE (LATITUDE)
#'
#' @author Matthew McLean, \email {mcleamj@gmail.com},
#'         
#' @date JUNE 9, 2023
########################################################################################

####################################
## DEFINE VARIABLE COLOR FUNCTION ##
####################################

variablecol <- function(colvar, col, clim) {
  ncol <- length(col)
  colvar[colvar < min(clim)] <- NA
  colvar[colvar > max(clim)] <- NA
  rn <- clim[2] - clim[1]
  ifelse (rn != 0, Col <- col[1 + trunc((colvar - clim[1])/rn *
                                          (ncol - 1)+1e-15)], Col <- rep(col[1], ncol))               
  return(Col)
}

######################
## LIBRARY PACKAGES ##
######################

if(!require(data.table)){install.packages("data.table"); library(data.table)}
if(!require(rgdal)){install.packages("rgdal"); library(rgdal)}
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

##################################
## ADD NEW IMPUTED BENTHIC DATA ##
##################################

benthic_PCA <- read_rds("outputs/RLS_benthic_PCA_imputed.rds")
benthic_PCA <- benthic_PCA %>%
  select(SurveyID, PC1_imputation, PC2_imputation) %>%
  rename(PC1_imputed = PC1_imputation, PC2_imputed = PC2_imputation)

model_data <- merge(model_data, benthic_PCA, by = "SurveyID", all=T)

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

z_score_2sd <- function(x){ (x - mean(x,na.rm=T)) / (2*sd(x,na.rm=T))}

z_vars <- model_data %>%
  select_if(is.numeric) %>%
  select(-any_of(c("SurveyID","SiteLongitude","SiteLatitude","aesthe_survey"))) %>%
  colnames()

standardized_data <- model_data %>%
  mutate_if(colnames(model_data) %in% z_vars, z_score_2sd)

###############################################
## MODEL USING SAME FORMAT FOR EFFECT OF MPA
## BUT WITH VARYING SLOPE FOR COUNTRY
###############################################

country_model_formula  <- bf(log(aesthe_survey) ~ MPA +
                           HDI2017 +
                           fshD +
                           gravtot2 +
                           as.factor(Temperature_Zone) +
                           (1 | SiteCode) +
                           (1 + MPA | Country),
                           
                         family=gaussian())

country_model <- brm(country_model_formula,
                 data=standardized_data,
                 chains=3, iter=3000, cores=ncores,
                 c(set_prior("normal(0,3)", class = "b"),
                   set_prior("normal(0,3)", class="Intercept")))

country_model_post <- as.data.frame(as.matrix(country_model)) %>%
  select('b_MPANotake','b_MPARestrictedtake')

##############################################
## MODEL USING SAME FORMAT FOR EFFECT OF MPA
## BUT WITH VARYING SLOPE FOR ECOREGION
##############################################

# RE-DO THIS ANALYSIS FOR ONLY ECOREGIONS WHERE THERE
# ARE BOTH FISHED AND NO TAKE SITES

ecoregion_model_formula  <- bf(log(aesthe_survey) ~ MPA + # TEST VARIABLE
                           HDI2017 + # CONTROL VARIABLE
                           fshD + # CONTROL VARIABLE
                           gravtot2 + # CONTROL VARIABLE
                           as.factor(Temperature_Zone) + # CONTROL VARIABLE
                           (1 | SiteCode) + # RANDOM INTERCEPTS FOR SITES
                           (1 + MPA | Ecoregion), # RANDOM INTERCEPT AND SLOPE FOR ECOREGION
                         
                         family=gaussian())

ecoregion_model <- brm(ecoregion_model_formula,
                 data=standardized_data,
                 chains=1, iter=2000, cores=ncores,
                 c(set_prior("normal(0,3)", class = "b"),
                   set_prior("normal(0,3)", class="Intercept")))
saveRDS(ecoregion_model, "ecoregion_model")
ecoregion_model <- read_rds("ecoregion_model.rds")

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

##########################################
## MAP OF MPA EFFECT FOR EACH ECOREGION ##
##########################################

ecoregion_info <- model_data %>%
  select(Ecoregion, SiteLongitude, SiteLatitude, PC1_imputed, PC2_imputed, nb_species, sst_mean) %>%
  group_by(Ecoregion) %>%
  summarise_all(.funs=mean) %>%
  rename(ECOREGION = Ecoregion)

# IMPORT SHAPEFILE OF ECOREGIONS
regions <- readOGR("data/MEOW-TNC", "meow_ecos")
regions <- regions[regions$ECOREGION %in% eco_no_take$ECOREGION,]

# COMBINE ECOREGION EFFECT SIZES AND INFORMATION WITH SHAPEFILE INFO
regions <- merge(regions, eco_no_take, by="ECOREGION")
regions <- merge(regions, ecoregion_info, by="ECOREGION")

# ASSIGN A COLOR TO EACH ECOREGION BASED ON EFFECT SIZE - CAN CHANGE PALETTE AND RANGE HERE
col_range <- c( max(abs(regions$Slope))*-1 , max(abs(regions$Slope)) )
regions$color <- variablecol(regions$Slope, col=(jet(nrow(regions))),
                             clim=col_range)

# MAKE MAP WITH BASE R
scatter2D(regions$SiteLongitude, regions$SiteLatitude, pch=19,
          colvar = regions$Slope, clim=col_range,
          cex=0,xlim=c(-180,180),ylim=c(-60,80),
          xlab="Longitude (°)",ylab="Latitude (°)",
          cex.axis=1.25,cex.lab=1.25)
plot(regions, border="black", col=regions$color,add=TRUE)

map(database = "world",fill = TRUE, border="grey80", 
    col="grey80", ylim=c(-70,70),add=TRUE)
eco_corners <- par("usr")
polygon(c(eco_corners[1], eco_corners[1], eco_corners[2], eco_corners[2]),
        c(eco_corners[3], eco_corners[4], eco_corners[4], eco_corners[3]),lwd=1.35)
title("Ecoregion-Level MPA Effect")

# PLOT EFFECT SIZE IN FUNCTION OF LATITUDE AND SPECIES RICHNESS
par(mfrow=c(1,2))
plot(regions$Slope ~ regions$SiteLatitude, pch=19,
     xlab="Latitude", ylab="MPA Effect Size")
title("Ecoregion-Level MPA Effect vs Latitude")

plot(regions$Slope ~ regions$nb_species, pch=19,
     xlab="Species Richness", ylab="MPA Effect Size")
title("Ecoregion-Level MPA Effect vs Species Richness")



####################################
## MAKE FIGURE OF MAP AND EFFECTS ##
####################################

## FIGURE LAYOUT ##

graphics.off()
par(oma=c(0,0,0,0))
par(mar=c(4,4,4,4))

m2 <- matrix(c(1,1,1,1,1,1,1,1,1,2,2,2,
               1,1,1,1,1,1,1,1,1,2,2,2,
               1,1,1,1,1,1,1,1,1,2,2,2,
               1,1,1,1,1,1,1,1,1,2,2,2,
               1,1,1,1,1,1,1,1,1,2,2,2,
               1,1,1,1,1,1,1,1,1,2,2,2),
             nrow = 6, 
             ncol= 12,
             byrow = TRUE)

layout(m2, widths = c(1,1,1,1,1,1,1,1,1,1,1,1),
       heights = c(1,1,1,1,1,1))

layout.show(n=2)


###############
## VERSION 1 ##
###############

# MAP PANEL #
col_range <- c( max(abs(regions$Slope))*-1 , max(abs(regions$Slope)) )
#col_range <- range(regions$Slope)
regions$color <- variablecol(regions$Slope, 
                             col=rev((brewer.rdylbu(nrow(regions)))),
                             clim=col_range)
scatter2D(regions$SiteLongitude, regions$SiteLatitude, pch=19,
          colvar = regions$Slope, clim=col_range,
          col=rev(brewer.rdylbu(nrow(regions))),
          cex=0,xlim=c(-180,180),ylim=c(-60,80),
          xlab="Longitude (°)",ylab="Latitude (°)",
          cex.axis=1.25,cex.lab=1.25)
plot(regions, border="black", col=regions$color,add=TRUE)
mtext(text="Effect Size",
      side=4,line=1.5,cex=1)
map(database = "world",fill = TRUE, border="grey80", 
    col="grey80", ylim=c(-70,70),add=TRUE)
eco_corners <- par("usr")
polygon(c(eco_corners[1], eco_corners[1], eco_corners[2], eco_corners[2]),
        c(eco_corners[3], eco_corners[4], eco_corners[4], eco_corners[3]),lwd=1.35)
title("Ecoregion-Level MPA Effect",line=1,
      font.main=1, cex.main=1.25)


# SLOPES PANEL #
fig_data <- eco_intercept %>%
  select(ECOREGION, Intercept) %>%
  rename(Fished = Intercept)

fig_data$No_Take <- (eco_intercept$Intercept + eco_no_take$Slope)
fig_data <- merge(fig_data, ecoregion_info[,c("ECOREGION","SiteLatitude"),by="ECOREGION"])
fig_data$abs_lat <- abs(fig_data$SiteLatitude)
fig_data$x1 <- rep(1)
fig_data$x2 <- rep(2)

scatter2D(x = rep(1,nrow(fig_data)),
          y=fig_data$Fished,
          xlim=c(0.9,2.1),
          ylim=c(7.05,8.15),
          colvar = abs(fig_data$SiteLatitude),
          col=rev(ocean.balance(n=nrow(fig_data))),
          cex=1,pch=19,
          xlab=NA, ylab="Log (Aesthetic Value)",
          cex.lab=1.25,
          xaxt="n")
axis(side=1,at=c(1,2),labels=c("Fished", "No Take MPA"),
     cex.axis=1.25)
scatter2D(x = rep(2,nrow(fig_data)),
          y=fig_data$No_Take,
          xlim=c(0.75,2.25),
          colvar = abs(fig_data$SiteLatitude),
          col=rev(ocean.balance(n=nrow(fig_data))),
          cex=1,pch=19,
          add=TRUE,colkey=FALSE)
segments2D(x0=fig_data$x1, x1=fig_data$x2,
           y0=fig_data$Fished, y1=fig_data$No_Take,
           colvar = fig_data$abs_lat,lwd=2,
           col=rev(ocean.balance(n=nrow(fig_data))),
           add=TRUE)
mtext(text="Absolute value of Latitude",
      side=4,line=1,cex=1)
title("Ecoregion-Level MPA Effect", line=1,
      font.main=1, cex.main=1.25)




###############
## VERSION 2 ##
###############

# MAP PANEL #
col_range <- c( max(abs(regions$Slope))*-1 , max(abs(regions$Slope)) )
#col_range <- range(regions$Slope)
regions$color <- variablecol(regions$Slope, 
                             col=rev((brewer.rdylbu(nrow(regions)))),
                             clim=col_range)
scatter2D(regions$SiteLongitude, regions$SiteLatitude, pch=19,
          colvar = regions$Slope, clim=col_range,
          col=rev(brewer.rdylbu(nrow(regions))),
          cex=0,xlim=c(-180,180),ylim=c(-60,80),
          xlab="Longitude (°)",ylab="Latitude (°)",
          cex.axis=1.25,cex.lab=1.25)
plot(regions, border="black", col=regions$color,add=TRUE)
mtext(text="Effect Size",
      side=4,line=1.5,cex=1)
map(database = "world",fill = TRUE, border="grey80", 
    col="grey80", ylim=c(-70,70),add=TRUE)
eco_corners <- par("usr")
polygon(c(eco_corners[1], eco_corners[1], eco_corners[2], eco_corners[2]),
        c(eco_corners[3], eco_corners[4], eco_corners[4], eco_corners[3]),lwd=1.35)
title("Ecoregion-Level MPA Effect",line=1,
      font.main=1, cex.main=1.25)

# FOREST PLOT PANEL #
forest_plot <- eco_no_take
forest_plot <- merge(forest_plot, ecoregion_info[,c("ECOREGION","SiteLatitude")],by="ECOREGION")
forest_plot$abs_lat <- abs(forest_plot$SiteLatitude)
forest_plot <- forest_plot %>%
  arrange(Slope)
par(lend=1)
x0 <- forest_plot$Q2.5
x1 <- forest_plot$Q97.5
y0 <- seq(1:nrow(forest_plot))
y1 <- seq(1:nrow(forest_plot))
segments2D(x0,y0,x1,y1, lwd=2, 
           colvar = forest_plot$abs_lat,
           col=rev(ocean.balance(nrow(forest_plot))),
           xlim=c(min(forest_plot$Q2.5,na.rm = TRUE),max(forest_plot$Q97.5,na.rm = TRUE)),
           ylim=c(min(seq(1:nrow(forest_plot))-0.25),max(seq(1:nrow(forest_plot))+0.25)),
           xlab="Effect Size", ylab=NA, yaxt = "n",
           cex.lab=1.25)
mtext(text="Absolute value of Latitude",
      side=4,line=1,cex=1)
title("Ecoregion-Level No Take MPA Effect", line=1,
      font.main=1, cex.main=1.25)
scatter2D(forest_plot$Slope, seq(1:nrow(forest_plot)),
          cex=1.2,
          colvar = forest_plot$abs_lat,
          col=rev(ocean.balance(nrow(forest_plot))),
          add=TRUE)
abline(v=0, lwd=2, col=adjustcolor("black",alpha.f = 0.75), lty=2)





#################
## VERSION 3   ##
#################

col_range <- c( max(abs(regions$Slope))*-1 , max(abs(regions$Slope)) )
regions$color <- variablecol(regions$Slope, rev(brewer.rdylbu(nrow(regions))),
                             clim=col_range)
# MAP PANEL #
col_range <- c( max(abs(regions$Slope))*-1 , max(abs(regions$Slope)) )
#col_range <- range(regions$Slope)
regions$color <- variablecol(regions$Slope, 
                             col=rev((brewer.rdylbu(nrow(regions)))),
                             clim=col_range)
scatter2D(regions$SiteLongitude, regions$SiteLatitude, pch=19,
          colvar = regions$Slope, clim=col_range,
          col=rev(brewer.rdylbu(nrow(regions))),
          cex=0,xlim=c(-180,180),ylim=c(-60,80),
          xlab="Longitude (°)",ylab="Latitude (°)",
          cex.axis=1.25,cex.lab=1.25)
plot(regions, border="black", col=regions$color,add=TRUE)
mtext(text="Effect Size",
      side=4,line=1.5,cex=1)
map(database = "world",fill = TRUE, border="grey80", 
    col="grey80", ylim=c(-70,70),add=TRUE)
eco_corners <- par("usr")
polygon(c(eco_corners[1], eco_corners[1], eco_corners[2], eco_corners[2]),
        c(eco_corners[3], eco_corners[4], eco_corners[4], eco_corners[3]),lwd=1.35)
title("Ecoregion-Level MPA Effect",line=1,
      font.main=1, cex.main=1.25)


# LATITUDE PLOT PANEL #

scatter2D(regions$Slope, regions$SiteLatitude, pch=19,
     xlab="Effect Size", ylab="Latitude",
     colvar = regions$Slope, clim=col_range,
     col=rev(brewer.rdylbu(nrow(regions))),
     cex=0)

points(regions$Slope, regions$SiteLatitude, pch=21,
     cex=1.25, col="black", bg=regions$color,
     xlab="Effect Size", ylab="Latitude")

plot(regions$Slope, regions$SiteLatitude, pch=21,
     cex=1.25, col="black", bg="grey",
     xlab="Effect Size", ylab="Latitude")



###################################
## IS THERE A LATITUDINAL BIAS IN 
## MPA/FISHED SITE OCCURENCE
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
