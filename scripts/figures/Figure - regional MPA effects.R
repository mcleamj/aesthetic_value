

#######################################################################################
#'  CODE FOR MAP OF ECOREGION LEVEL MPA EFFECTS
#'
#' @author Matthew McLean, \email {mcleamj@gmail.com},
#'         
#' @date JUNE 9, 2023
########################################################################################

######################
## LIBRARY PACKAGES ##
######################

if(!require(sf)){install.packages("sf"); library(sf)}
if(!require(mapproj)){install.packages("mapproj"); library(mapproj)}
if(!require(plot3D)){install.packages("plot3D"); library(plot3D)}
if(!require(ggplot2)){install.packages("ggplot2"); library(ggplot2)}
if(!require(pals)){install.packages("pals"); library(pals)}
if(!require(readr)){install.packages("readr"); library(readr)}
if(!require(dplyr)){install.packages("dplyr"); library(dplyr)}

source("R/variablecol.R")

#######################
## IMPORT MODEL DATA ##
#######################

model_data <- read_rds("data/model_data.rds")

#######################################
## IMPORT VARYING SLOPE MODEL OUTPUT ##
#######################################

ecoregion_model <- read_rds("outputs/BIG_FILES/ecoregion_varying_slopes_model.rds")

################################
## EXTRACT MODEL COEFFICIENTS ##
################################

eco_no_take <- as.data.frame(coef(ecoregion_model)$Ecoregion[,,"MPANotake"]) %>%
  rename(Slope = "Estimate")
eco_no_take$ECOREGION <- rownames(eco_no_take)

##########################################
## MAP OF MPA EFFECT FOR EACH ECOREGION ##
##########################################

ecoregion_info <- model_data %>%
  dplyr::select(Ecoregion, SiteLongitude, SiteLatitude) %>%
  group_by(Ecoregion) %>%
  summarise_all(.funs=mean) %>%
  rename(ECOREGION = Ecoregion)

# IMPORT SHAPEFILE OF ECOREGIONS
#regions <- readOGR("data/MEOW-TNC", "meow_ecos")
regions <- read_sf("data/MEOW-TNC", "meow_ecos")
regions <- regions[regions$ECOREGION %in% eco_no_take$ECOREGION,]

# COMBINE ECOREGION EFFECT SIZES AND INFORMATION WITH SHAPEFILE INFO
regions <- merge(regions, eco_no_take, by="ECOREGION")
regions <- merge(regions, ecoregion_info, by="ECOREGION")

# ASSIGN A COLOR TO EACH ECOREGION BASED ON EFFECT SIZE - CAN CHANGE PALETTE AND RANGE HERE
col_range <- c( max(abs(regions$Slope))*-1 , max(abs(regions$Slope)) )
regions$color <- variablecol(regions$Slope, col=(jet(nrow(regions))),
                             clim=col_range)


###############################
## SAVE THE FILE AS A FIGURE ##
###############################

# graphics.off()
# tiff(file=here::here("figures_tables","MPA_Ecoregions.tiff"),
#      width = 12,
#      height =4.5,
#      units="in",
#      res = 300)

####################################
## MAKE FIGURE OF MAP AND EFFECTS ##
####################################

## FIGURE LAYOUT ##

par(oma=c(0,0,0,0))
par(mar=c(4,4,4,4))

m2 <- matrix(c(1,1,1,1,1,1,1,1,1,2,2,
               1,1,1,1,1,1,1,1,1,2,2,
               1,1,1,1,1,1,1,1,1,2,2,
               1,1,1,1,1,1,1,1,1,2,2,
               1,1,1,1,1,1,1,1,1,2,2,
               1,1,1,1,1,1,1,1,1,2,2),
             nrow = 6, 
             ncol= 11,
             byrow = TRUE)

layout(m2, widths = c(1,1,1,1,1,1,1,1,1,1.1,1.1),
       heights = c(1,1,1,1,1,1))

#layout.show(n=2)


#################
## FIGURE  ##
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
corners <- par("usr")
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = adjustcolor("grey80",alpha=0.2))
plot(regions, border="black", col=regions$color,add=TRUE)
mtext(text="Effect Size",
      side=4,line=1.5,cex=1)
map(database = "world",fill = TRUE, border="grey80", 
    col="grey80", ylim=c(-70,70),add=TRUE)
eco_corners <- par("usr")
polygon(c(eco_corners[1], eco_corners[1], eco_corners[2], eco_corners[2]),
        c(eco_corners[3], eco_corners[4], eco_corners[4], eco_corners[3]),lwd=1.35)
title("Predicted Effect of No-Take MPAs by Ecoregion",line=1,
      font.main=1, cex.main=1.25)


# LATITUDE PLOT PANEL #

scatter2D(regions$Slope, regions$SiteLatitude, pch=19,
          xlab="Effect Size", ylab=NA, yaxt="n",
          colvar = regions$Slope, clim=col_range,
          col=rev(brewer.rdylbu(nrow(regions))),
          cex=0,
          colkey = FALSE)

points(regions$Slope, regions$SiteLatitude, pch=21,
       cex=1.75, col="black", bg=regions$color,
       xlab="Effect Size", ylab=NA)

#dev.off() 