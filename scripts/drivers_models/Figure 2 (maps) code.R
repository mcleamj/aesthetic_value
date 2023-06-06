
#######################################################################################
#'  FIGURE 2
#'  GLOBAL MAPS OF AESTHETIC VALUE AND ITS DRIVERS
#'  
#'
#' @author Matthew McLean, \email {mcleamj@gmail.com},
#'         
#'
#' @date
########################################################################################

######################
## LIBRARY PACKAGES ##
######################

library(sf)
library(raster)
library(geosphere)
library(plot3D)
library(mapproj)
library(pals)
library(brms)
library(dplyr)

###########################################
## IMPORT AND MERGE SURVEY AND SITE DATA ##
###########################################

survey_esth <- read.csv("outputs/survey_aesth.csv")
site_info <- readRDS("data/RLS_sitesInfos_tropical.rds")

##################################################
## CALCULATE RESIDUAL (DEVIATION FROM EXPECTED) ##
##################################################

survey_esth$residual <- survey_esth$aesthe_survey - survey_esth$aesthe_SR_survey

###################
## PLOT THE DATA ##
###################

scatter2D(survey_esth$nb_species, survey_esth$aesthe_survey,
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

site_esth$log_esth <- log(site_esth$aesthe_survey)

##############################################################################################
## CALCULATE INVERSE DISTANCE WEIGHTING INTERPOLATION FOR LOG AESTHETIC VALUE AND RESIDUALS ##
##############################################################################################

IDW_sites <- site_esth %>%
  select(SiteLongitude, SiteLatitude, log_esth, residual)

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
saveRDS(log_esth_IDW, "log_esth_IDW.rds")

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
saveRDS(residual_IDW, "residual_IDW.rds")

##################################
## MAKE THE FIGURE USING BASE R ##
##################################

## FIGURE LAYOUT ##

graphics.off()
par(oma=c(0,0,0,0))
par(mar=c(4,4,4,4))

m2 <- matrix(c(1,1,1,1,1,1,1,0,0,0,
               1,1,1,1,1,1,1,0,2,2,
               1,1,1,1,1,1,1,0,2,2,
               1,1,1,1,1,1,1,0,2,2,
               1,1,1,1,1,1,1,0,2,2,
               1,1,1,1,1,1,1,0,2,2,
               3,3,3,3,3,3,3,0,2,2,
               3,3,3,3,3,3,3,0,2,2,
               3,3,3,3,3,3,3,0,2,2,
               3,3,3,3,3,3,3,0,2,2,
               3,3,3,3,3,3,3,0,2,2,
               3,3,3,3,3,3,3,0,0,0),
             nrow = 12, 
             ncol= 10,
             byrow = TRUE)

layout(m2, widths = c(1,1,1,1,1,1,1,0.8,1,1),
       heights = c(1,1,1,1,1,1,1,1,1,1,1,1))

layout.show(n=3)

##################################
## PLOT THE AESTHETIC VALUE MAP ##
##################################

log_esth_IDW <- readRDS("log_esth_IDW.rds")

scatter2D(log_esth_IDW$lon, log_esth_IDW$lat, pch=19,
          colvar = log_esth_IDW$log_esth, cex=0,
          col=(jet(n=nrow(IDW_sites))),
          xlab="Longitude", ylab="Latitude", 
          cex.lab=1.5)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = adjustcolor("grey80",alpha=0.2))
scatter2D(log_esth_IDW$lon, log_esth_IDW$lat, pch=19,
          colvar = log_esth_IDW$log_esth, cex=0.75,
          col=(jet(n=nrow(IDW_sites))),
          add=TRUE)
mtext(side=4, line=1, text="log(Aesthetic Value)")
mtext(side=3, adj=-0.05, text="A", line=1.5, font=2, cex=1.5)

map(database = "world", fill=T, col="grey85", border="grey85",add=T)

title("Global Distribution of Aesthetic Value", line=1,
      font.main=1, cex.main=1.5)

##########################
## PLOT THE FOREST PLOT ##
##########################

# model_output <- read.table("dag_output.txt")
# model_coefs <- data.frame(posterior_summary(model_output, probs = c(0.025,0.1,0.25,0.75,0.9,0.975)))
# 
# write.csv(model_coefs, "model_coefs.csv")

model_coefs <- read.csv("model_coefs.csv")
model_coefs$abs_effect <- abs(model_coefs$Estimate)

model_coefs <- model_coefs %>%
  arrange(type, abs_effect)
# model_coefs$variable <- rownames(model_coefs)
# model_coefs$variable <- gsub("b_","",model_coefs$variable)

plot(model_coefs$Estimate, seq(1:nrow(model_coefs)), xlim=c(min(model_coefs$Q2.5,na.rm = TRUE),max(model_coefs$Q97.5,na.rm = TRUE)),
     ylim=c(min(seq(1:nrow(model_coefs))-0.25),max(seq(1:nrow(model_coefs))+0.25)),cex=0,
     xlab="Standardized Effect Size", ylab=NA, yaxt = "n",
     cex.lab=1.5)
title("Drivers of Aesthetic Value", line=1,
      font.main=1, cex.main=1.5)

corners <- par("usr")

polygon(x=c(corners[2],corners[2],corners[1],corners[1]),
        y=c(corners[3],
            min(which(model_coefs$type=="social")-0.5),
            min(which(model_coefs$type=="social")-0.5),
            corners[3]),
        col=adjustcolor("lightgreen",alpha.f = 0.15),lwd=1)

polygon(x=c(corners[2],corners[2],corners[1],corners[1]),
        y=c(max(which(model_coefs$type=="env")+0.5),
            corners[4], corners[4],
            max(which(model_coefs$type=="env")+0.5)),
        col=adjustcolor("lightblue",alpha.f = 0.15),lwd=1)

par(lend=1)
x0 <- model_coefs$Q2.5
x1 <- model_coefs$Q97.5
y0 <- seq(1:nrow(model_coefs))
y1 <- seq(1:nrow(model_coefs))
segments(x0,y0,x1,y1, lwd=2)
x0 <- model_coefs$Q25
x1 <- model_coefs$Q75
y0 <- seq(1:nrow(model_coefs))
y1 <- seq(1:nrow(model_coefs))
segments(x0,y0,x1,y1, lwd=5)

points(model_coefs$Estimate, seq(1:nrow(model_coefs)),pch=21,col=1,bg="grey", cex=2.25)
abline(v=0, lwd=1.5, col=adjustcolor("black",alpha.f = 0.75), lty=2)

axis(2, at = seq(1:nrow(model_coefs)), 
     labels = model_coefs$label,
     las=2, cex.axis=1.25)

mtext(side=3, adj=-0.05, text="C", line=1.5, font=2, cex=1.5)
# mtext(side=4, line=1, text="Social", 0.75)
# mtext(side=4, line=1, text="Environmental", adj=0.25)

#legend("bottomright", legend="Social",bty="n",bg=NA, cex=1)
#legend("bottomright", legend="Environmental",bty="n",bg=NA)

text("Social", x = corners[2]-0.05, y = (max(which(model_coefs$type=="env")+0.5) + corners[4])/2 , 
     srt = 90, adj =c(0.5), cex=1.5)

text("Environmental", x = corners[2]-0.05, y = (max(which(model_coefs$type=="env")+0.5) + corners[3])/2 , 
     srt = 90, adj =c(0.5), cex=1.5)

###########################
## PLOT THE RESIDUAL MAP ##
###########################

residual_IDW <- readRDS("residual_IDW.rds")
resid_clim <- c(-max(abs(residual_IDW$residual)),max(abs(residual_IDW$residual)))
scatter2D(residual_IDW$lon, residual_IDW$lat, pch=19,
          colvar = residual_IDW$residual, 
          cex=0,
          col=(jet(n=nrow(IDW_sites))),
          clim=resid_clim,
          xlab="Longitude", ylab="Latitude", 
          cex.lab=1.5)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = adjustcolor("grey80",alpha=0.2))
scatter2D(residual_IDW$lon, residual_IDW$lat, pch=19,
          colvar = residual_IDW$residual, cex=0.75,
          col=(jet(n=nrow(IDW_sites))),
          clim=resid_clim,
          add=TRUE)
mtext(side=4, line=1, text="Residual Value")
mtext(side=3, adj=-0.05, text="B", line=1.5, font=2, cex=1.5)

map(database = "world", fill=T, col="grey85", border="grey85",add=T)

title("Deviation from Expected Given Species Richness",line=1,
      font.main=1, cex.main=1.5)


##################################
## MAKE THE FIGURE USING GGPLOT ##
##################################













#######################
## WITH ONLY ONE MAP ##
#######################

##################################
## MAKE THE FIGURE USING BASE R ##
##################################

## FIGURE LAYOUT ##

graphics.off()
par(oma=c(0,0,0,0))
par(mar=c(4,4,4,4))

m2 <- matrix(c(1,1,1,1,1,1,1,0,2,2,
               1,1,1,1,1,1,1,0,2,2,
               1,1,1,1,1,1,1,0,2,2,
               1,1,1,1,1,1,1,0,2,2,
               1,1,1,1,1,1,1,0,2,2),
             nrow = 5, 
             ncol= 10,
             byrow = TRUE)

layout(m2, widths = c(1,1,1,1,1,1,1,0.8,1,1),
       heights = c(1,1,1,1,1))

layout.show(n=2)

##################################
## PLOT THE AESTHETIC VALUE MAP ##
##################################

log_esth_IDW <- readRDS("log_esth_IDW.rds")

scatter2D(log_esth_IDW$lon, log_esth_IDW$lat, pch=19,
          colvar = log_esth_IDW$log_esth, cex=0,
          col=(jet(n=nrow(IDW_sites))),
          xlab="Longitude", ylab="Latitude", 
          cex.lab=1.5)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = adjustcolor("grey80",alpha=0.2))
scatter2D(log_esth_IDW$lon, log_esth_IDW$lat, pch=19,
          colvar = log_esth_IDW$log_esth, cex=0.75,
          col=(jet(n=nrow(IDW_sites))),
          add=TRUE)
mtext(side=4, line=1, text="log(Aesthetic Value)")
mtext(side=3, adj=-0.05, text="A", line=1.5, font=2, cex=1.5)

map(database = "world", fill=T, col="grey85", border="grey85",add=T)

title("Global Distribution of Aesthetic Value", line=1,
      font.main=1, cex.main=1.5)

##########################
## PLOT THE FOREST PLOT ##
##########################

# model_output <- read.table("dag_output.txt")
# model_coefs <- data.frame(posterior_summary(model_output, probs = c(0.025,0.1,0.25,0.75,0.9,0.975)))
# 
# write.csv(model_coefs, "model_coefs.csv")

model_coefs <- read.csv("model_coefs.csv")
model_coefs$abs_effect <- abs(model_coefs$Estimate)

model_coefs <- model_coefs %>%
  arrange(type, abs_effect)
# model_coefs$variable <- rownames(model_coefs)
# model_coefs$variable <- gsub("b_","",model_coefs$variable)

plot(model_coefs$Estimate, seq(1:nrow(model_coefs)), xlim=c(min(model_coefs$Q2.5,na.rm = TRUE),max(model_coefs$Q97.5,na.rm = TRUE)),
     ylim=c(min(seq(1:nrow(model_coefs))-0.25),max(seq(1:nrow(model_coefs))+0.25)),cex=0,
     xlab="Standardized Effect Size", ylab=NA, yaxt = "n",
     cex.lab=1.5)
title("Drivers of Aesthetic Value", line=1,
      font.main=1, cex.main=1.5)

corners <- par("usr")

polygon(x=c(corners[2],corners[2],corners[1],corners[1]),
        y=c(corners[3],
            min(which(model_coefs$type=="social")-0.5),
            min(which(model_coefs$type=="social")-0.5),
            corners[3]),
        col=adjustcolor("lightgreen",alpha.f = 0.15),lwd=1)

polygon(x=c(corners[2],corners[2],corners[1],corners[1]),
        y=c(max(which(model_coefs$type=="env")+0.5),
            corners[4], corners[4],
            max(which(model_coefs$type=="env")+0.5)),
        col=adjustcolor("lightblue",alpha.f = 0.15),lwd=1)

par(lend=1)
x0 <- model_coefs$Q2.5
x1 <- model_coefs$Q97.5
y0 <- seq(1:nrow(model_coefs))
y1 <- seq(1:nrow(model_coefs))
segments(x0,y0,x1,y1, lwd=2)
x0 <- model_coefs$Q25
x1 <- model_coefs$Q75
y0 <- seq(1:nrow(model_coefs))
y1 <- seq(1:nrow(model_coefs))
segments(x0,y0,x1,y1, lwd=5)

points(model_coefs$Estimate, seq(1:nrow(model_coefs)),pch=21,col=1,bg="grey", cex=2.25)
abline(v=0, lwd=1.5, col=adjustcolor("black",alpha.f = 0.75), lty=2)

axis(2, at = seq(1:nrow(model_coefs)), 
     labels = model_coefs$label,
     las=2, cex.axis=1.25)

mtext(side=3, adj=-0.05, text="C", line=1.5, font=2, cex=1.5)
# mtext(side=4, line=1, text="Social", 0.75)
# mtext(side=4, line=1, text="Environmental", adj=0.25)

#legend("bottomright", legend="Social",bty="n",bg=NA, cex=1)
#legend("bottomright", legend="Environmental",bty="n",bg=NA)

text("Social", x = corners[2]-0.05, y = (max(which(model_coefs$type=="env")+0.5) + corners[4])/2 , 
     srt = 90, adj =c(0.5), cex=1.5)

text("Environmental", x = corners[2]-0.05, y = (max(which(model_coefs$type=="env")+0.5) + corners[3])/2 , 
     srt = 90, adj =c(0.5), cex=1.5)

