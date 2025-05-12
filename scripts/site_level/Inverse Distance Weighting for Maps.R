
#########################################################################################################
#' SCRIPT TO COMPUTE INVERSE DISTANCE WEIGHTING INTERPOLATION OF SITE-LEVEL AESTEHTIC VALUE FOR MAPPING
#'
#' @author Ella Clausius \email {ella.clausius@@utas.edu.au}
#'         Matthew McLean, \email {mcleamj@@gmail.com}
#' 
#########################################################################################################

library(sf)
library(raster)
library(geosphere)
library(plot3D)
library(mapproj)
library(pals)
library(readr)
library(dplyr)

###########################################
## IMPORT AND MERGE SURVEY AND SITE DATA ##
###########################################

survey_esth <- read.csv("outputs/survey_aesth.csv")
site_info <- readRDS("data/RLS_sitesInfos_tropical.rds")

##################################################
## CALCULATE RESIDUAL (DEVIATION FROM EXPECTED) ##
##################################################

survey_esth$residual <- survey_esth$aesthe_survey_abund - survey_esth$aesthe_SR_survey

###################
## PLOT THE DATA ##
###################

scatter2D(survey_esth$nb_species, survey_esth$aesthe_survey_abund,
          pch=19, cex=1, colvar = survey_esth$residual,
          xlab="Species Richness", ylab="Observed Aesthetic Value")
points(survey_esth$nb_species, survey_esth$aesthe_SR_survey, pch=19, col=1)

#########################################
## MERGE DATA AND CALCULATE SITE MEANS ##
#########################################

survey_site <- site_info %>%
  select(SurveyID, SiteCode, SiteLatitude, SiteLongitude)

survey_esth <- merge(survey_esth, survey_site, by="SurveyID")

site_esth <- survey_esth %>%
  select(-SurveyID) %>%
  group_by(SiteCode) %>%
  summarise_all(mean)

####################################################
## CALCULATE LOG OF AESTHETIC VALUE FOR EACH SITE ##
####################################################

site_esth$log_esth <- log(site_esth$aesthe_survey_abund)

##########################################
## IMPORT AND ATTACH MODELED INTERCEPTS ##
##########################################

site_intercepts <- readRDS("outputs/site_intercepts.rds")

site_esth <- merge(site_esth, site_intercepts, by="SiteCode")

########################################################
## CALCULATE INVERSE DISTANCE WEIGHTING INTERPOLATION ##
########################################################

IDW_sites <- site_esth %>%
  select(SiteLongitude, SiteLatitude, log_esth, Estimate.Intercept, residual)

sites_xy <- cbind(IDW_sites$SiteLongitude, IDW_sites$SiteLatitude)

res <- c(1000, 1000) #resolution of map x, y
my_buffer <- 150*1000 # The size of the buffers around the points in meters
p <- 1.6 #power for IDW

#create raster with resolution defined above 
r <- raster(nrow = res[1], ncol = res[2]) 

#pulls raster cells that match coordinates from sites_xy
icell <- cellFromXY(r, sites_xy) 

#assigns grid cells present in icell a value of 1 
r[icell] <- 1 

#computes the distance for all cells that are NA to the nearest cell that is not NA (i.e., to the sites in sites_xy)
rdist <- distance(r) ##this may take a while to compute if working with a large number of sites 

#list of grid cells in rdist that fall within the defined buffer 
ifilter <- which(rdist[] < my_buffer)  

#generate a list of coordinates for each of the grid cells that fall within the buffer of each site
xyinterp <- xyFromCell(r, ifilter) 

#calculate the distance of each grid cell inside the buffer from each site
xydist <- distm(sites_xy, xyinterp)/1000 

#####################
## AESTHETIC VALUE ##
#####################

#pull values of interest
values_esth <- IDW_sites %>% 
  pull(log_esth)

#weights each grid cell value by it's distance from the site according to the power you assign. 
w <- 1/(xydist^p) 

indic <- matrix(values_esth, nrow = 1)
isreal <- which(!is.na(indic))
nvals <- length(isreal)

val <- (indic[isreal] %*% w[isreal,])/colSums(w[isreal,])
val <- as.numeric(val)

rval <- r
rval[ifilter] <- as.numeric(val)
esth_rast <- rval 

## CREATE DATA FRAME
log_esth_IDW <- data.frame(xyinterp, val)
names(log_esth_IDW) <- c("lon", "lat", "log_esth")

log_esth_IDW <- log_esth_IDW %>%
  arrange(log_esth)

log_esth_IDW$lon_2 <- ifelse(log_esth_IDW$lon <0, log_esth_IDW$lon+360, log_esth_IDW$lon)

## SAVE OUTPUT
saveRDS(log_esth_IDW, "outputs/log_esth_IDW.rds")



#######################
## MODELED INTERCEPT ##
#######################

#pull values of interest
values_intercept <- IDW_sites %>% 
  pull(Estimate.Intercept)

#weights each grid cell value by it's distance from the site according to the power you assign. 
w <- 1/(xydist^p) 

indic <- matrix(values_intercept, nrow = 1)
isreal <- which(!is.na(indic))
nvals <- length(isreal)

val <- (indic[isreal] %*% w[isreal,])/colSums(w[isreal,])
val <- as.numeric(val)

rval <- r
rval[ifilter] <- as.numeric(val)
esth_rast <- rval 

## CREATE DATA FRAME
intercept_IDW <- data.frame(xyinterp, val)
names(intercept_IDW) <- c("lon", "lat", "Estimate.Intercept")

intercept_IDW <- intercept_IDW %>%
  arrange(Estimate.Intercept)

intercept_IDW$lon_2 <- ifelse(intercept_IDW$lon <0, intercept_IDW$lon+360, intercept_IDW$lon)

## SAVE OUTPUT
saveRDS(intercept_IDW, "outputs/intercept_IDW.rds")


####################
## RESIDUAL VALUE ##
####################

#pull values of interest
values_resid <- IDW_sites %>% 
  pull(residual)

#weights each grid cell value by it's distance from the site according to the power you assign. 
w <- 1/(xydist^p) 

indic <- matrix(values_resid, nrow = 1)
isreal <- which(!is.na(indic))
nvals <- length(isreal)

val <- (indic[isreal] %*% w[isreal,])/colSums(w[isreal,])
val <- as.numeric(val)

rval <- r
rval[ifilter] <- as.numeric(val)
esth_rast <- rval 

## CREATE DATA FRAME
residual_IDW <- data.frame(xyinterp, val)
names(residual_IDW) <- c("lon", "lat", "residual")

residual_IDW <- residual_IDW %>%
  arrange(residual)

residual_IDW$lon_2 <- ifelse(residual_IDW$lon <0, residual_IDW$lon+360, residual_IDW$lon)

## SAVE OUTPUT
saveRDS(residual_IDW, "outputs/residual_IDW.rds")
