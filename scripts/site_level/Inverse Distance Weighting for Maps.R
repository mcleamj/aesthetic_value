
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
library(dplyr)

################################################
## IMPORT SURVEY-LEVEL AESTHETIC VALUE SCORES
## AND CALCULATE RAW MEAN PER SITE
################################################

aesth_survey <- read.csv("outputs/survey_aesth.csv")

site_info <- read_rds("data/RLS_sitesInfos_tropical.rds") %>%
  select(SurveyID, SiteCode, SiteLatitude, SiteLongitude)

aesth_survey <- merge(aesth_survey, site_info, by="SurveyID")

aesth_site <- aesth_survey %>%
  select(SiteCode, SiteLatitude, SiteLongitude, aesthe_survey) %>%
  group_by(SiteCode) %>%
  summarise_all(.funs=mean)

aesth_site$log_aesth <- log(aesth_site$aesthe_survey) # CREATE LOG-TRANSFORMED VERSION OF DATA

scatter2D(aesth_site$SiteLongitude, aesth_site$SiteLatitude, pch=19,
          colvar = aesth_site$log_aesth)

res <- c(1000, 1000) #resolution of map x, y
my_buffer <- 150*1000 # The size of the buffers around the points in metres
p <- 1.6 #power for IDW 
#lower values mean its more influenced by distances further away
# higher values = less smoothing 

sites <- aesth_site

sites_xy <- cbind(sites$SiteLongitude, sites$SiteLatitude)

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

#pull values of interest, which in this case is b20
values <- sites %>% 
  pull(log_aesth)

#weights each grid cell value by it's distance from the site according to the power you assign. 
w <- 1/(xydist^p) 

indic <- matrix(values, nrow = 1)
isreal <- which(!is.na(indic))
nvals <- length(isreal)

val <- (indic[isreal] %*% w[isreal,])/colSums(w[isreal,])
val <- as.numeric(val)

rval <- r
rval[ifilter] <- as.numeric(val)
esth_rast <- rval 

#create dataframe of coordinates and values
dat1 <- data.frame(xyinterp, val)
names(dat1) <- c("lon", "lat", "log_aesth")

dat1 <- dat1 %>%
  arrange(log_aesth)

dat1$lon_2 <- ifelse(dat1$lon <0, dat1$lon+360, dat1$lon) #REARRANGE LAT/LON


graphics.off()
scatter2D(dat1$lon_2, dat1$lat, pch=19,
          colvar = dat1$log_aesth, cex=0.75,
          col=(jet(n=nrow(sub_esth))))
map(database = "world2", fill=T, col="grey90", border="grey90",add=T)

#################
## SAVE OUTPUT ##
#################

saveRDS(dat1, "outputs/IDW_raw_site_mean.rds")

graphics.off()
scatter2D(dat1$lon, dat1$lat, pch=19,
          colvar = dat1$log_aesth, cex=0.75,
          col=(jet(n=nrow(sub_esth))))

map(database = "world", fill=T, col="grey90", border="grey90",add=T)
map(database = "world", fill=T, col="grey60", border="grey60",add=T)




