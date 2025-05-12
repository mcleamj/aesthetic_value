
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
library(readr)
library(dplyr)

###################################
## IMPORT NECESSARY DATA OUTPUTS ##
###################################

# IDW

aes_models <- read_rds("outputs/dag_output.rds")
aes_coefs <- as.data.frame(posterior_summary(aes_models, 
                                             probs=c(0.05,0.25,0.75,0.95)))

aes_coefs$variable <- rownames(aes_coefs)
aes_coefs$type <- aes_coefs$variable
aes_coefs <- aes_coefs %>%
  mutate(type=recode(type, 
                     "Sea Surface Temperature" = "Environmental",
                     "Log Human Gravity" = "Anthropogenic",
                     "Net Primary Productivity" = "Environmental",
                     "Depth" = "Environmental",
                     "Fisheries Dependency" = "Anthropogenic",
                     "Human Dev. Index" = "Anthropogenic",
                     "No Take MPA"  = "Anthropogenic",
                     "Restricted Take MPA" = "Anthropogenic",
                     "Benthic Composition (PC1)" = "Environmental",
                     "Benthic Composition (PC2)" = "Environmental",
                     "Degree Heating Weeks" = "Environmental"))

aes_coefs$abs_effect <- abs(aes_coefs$Estimate)

aes_coefs <- aes_coefs %>%
  arrange(type, abs_effect)

deviation_models <- read_rds("outputs/dag_output_richness.rds")
deviation_coefs <- as.data.frame(posterior_summary(deviation_models, probs=c(0.05,0.25,0.75,0.95)))

deviation_coefs$variable <- rownames(deviation_coefs)
deviation_coefs$type <- deviation_coefs$variable
deviation_coefs <- deviation_coefs %>%
  mutate(type=recode(type, 
                     "Sea Surface Temperature" = "Environmental",
                     "Log Human Gravity" = "Anthropogenic",
                     "Net Primary Productivity" = "Environmental",
                     "Depth" = "Environmental",
                     "Fisheries Dependency" = "Anthropogenic",
                     "Human Dev. Index" = "Anthropogenic",
                     "No Take MPA"  = "Anthropogenic",
                     "Restricted Take MPA" = "Anthropogenic",
                     "Benthic Composition (PC1)" = "Environmental",
                     "Benthic Composition (PC2)" = "Environmental",
                     "Degree Heating Weeks" = "Environmental"))

deviation_coefs$abs_effect <- abs(deviation_coefs$Estimate)

deviation_coefs <- deviation_coefs[order(match(deviation_coefs$variable, aes_coefs$variable)),]


log_esth_IDW <- readRDS("outputs/log_esth_IDW.rds")

residual_IDW <- readRDS("outputs/residual_IDW.rds")


#########################
## DEFINE POINT COLORS ##
#########################

anthro_color <- adjustcolor("grey50",alpha.f = 1)
env_color <- adjustcolor("grey80", alpha.f = 1)

# anthro_color <- adjustcolor("blue",alpha.f = 0.5)
# env_color <- adjustcolor("green", alpha.f = 0.5)


###############################
## SAVE THE FILE AS A FIGURE ##
###############################

# graphics.off()
# tiff(file=here::here("figures_tables","TEST_TEST.tiff"), 
#      # width = 10*3.5, 
#      # height = 6*3.5,
#      width = 17.8,
#      height= 10.7,
#      units="cm",
#      res = 600)

##################################
## MAKE THE FIGURE USING BASE R ##
##################################

## FIGURE LAYOUT ##

par(oma=c(0,0,0,0))
par(mar=c(4,4,4,4))

m2 <- matrix(c(1,1,1,1,1,1,1,1,0,0,2,2,2,
               1,1,1,1,1,1,1,1,0,0,2,2,2,
               1,1,1,1,1,1,1,1,0,0,2,2,2,
               1,1,1,1,1,1,1,1,0,0,2,2,2,
               1,1,1,1,1,1,1,1,0,0,2,2,2,
               1,1,1,1,1,1,1,1,0,0,2,2,2,
               3,3,3,3,3,3,3,3,0,0,4,4,4,
               3,3,3,3,3,3,3,3,0,0,4,4,4,
               3,3,3,3,3,3,3,3,0,0,4,4,4,
               3,3,3,3,3,3,3,3,0,0,4,4,4,
               3,3,3,3,3,3,3,3,0,0,4,4,4,
               3,3,3,3,3,3,3,3,0,0,4,4,4),
             nrow = 12, 
             ncol= 13,
             byrow = TRUE)

layout(m2, widths = c(1,1,1,1,1,1,1,1,1,0.5,1,1,1),
       heights = c(1,1,1,1,1,1,1,1,1,1,1,1))

layout.show(n=4)


##################################
## PLOT THE AESTHETIC VALUE MAP ##
##################################

scatter2D(log_esth_IDW$lon, log_esth_IDW$lat, pch=19,
          colvar = log_esth_IDW$log_esth, cex=0,
          col=(jet(n=nrow(log_esth_IDW))),
          xlab="Longitude", ylab="Latitude", 
          cex.lab=1.5)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = adjustcolor("grey80",alpha=0.2))
#rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = adjustcolor("grey20",alpha=0.2))
scatter2D(log_esth_IDW$lon, log_esth_IDW$lat, pch=19,
          colvar = log_esth_IDW$log_esth, cex=0.75,
          col=(jet(n=nrow(log_esth_IDW))),
          add=TRUE)
mtext(side=4, line=1.5, text="log(Aesthetic Value)")
mtext(side=3, adj=-0.05, text="A", line=1.5, font=2, cex=1.25)

map(database = "world", fill=T, col="grey80", border="grey80",add=T)
#map(database = "world", fill=T, col="grey70", border="grey70",add=T)
#map(database = "world", fill=T,col="grey80", border="black",add=T)

title("Global Distribution of Aesthetic Value", line=1,
      font.main=1, cex.main=1.5)

##########################################
## PLOT THE AESTHETIC VALUE FOREST PLOT ##
##########################################

plot(aes_coefs$Estimate, seq(1:nrow(aes_coefs)), xlim=c(min(aes_coefs$Q5,na.rm = TRUE),max(aes_coefs$Q95,na.rm = TRUE)),
     ylim=c(min(seq(1:nrow(aes_coefs))-0.25),max(seq(1:nrow(aes_coefs))+0.25)),cex=0,
     xlab="Standardized Effect Size", ylab=NA, yaxt = "n",
     cex.lab=1.5)
title("Drivers of Aesthetic Value", line=1,
      font.main=1, cex.main=1.5)

corners <- par("usr")

abline(v=0, lwd=1.5, col=adjustcolor("black",alpha.f = 0.75), lty=2)

par(lend=1)
x0 <- aes_coefs$Q5
x1 <- aes_coefs$Q95
y0 <- seq(1:nrow(aes_coefs))
y1 <- seq(1:nrow(aes_coefs))
segments(x0,y0,x1,y1, lwd=2)
x0 <- aes_coefs$Q25
x1 <- aes_coefs$Q75
y0 <- seq(1:nrow(aes_coefs))
y1 <- seq(1:nrow(aes_coefs))
segments(x0,y0,x1,y1, lwd=5)
points(aes_coefs$Estimate, seq(1:nrow(aes_coefs)),col="white",bg="white",cex=1.75,
       pch=ifelse(aes_coefs$type=="Anthropogenic",21,21))
points(aes_coefs$Estimate, seq(1:nrow(aes_coefs)),col=1,cex=1.75,
       pch=ifelse(aes_coefs$type=="Anthropogenic",21,21),
       bg=ifelse(aes_coefs$type=="Anthropogenic",anthro_color,env_color))

h_line <- which(aes_coefs$type=="Environmental")[1]-0.5
abline(h=h_line, col=adjustcolor("black",alpha.f = 0.75), lty=1)

axis(2, at = seq(1:nrow(aes_coefs)), 
     labels = aes_coefs$variable,
     las=2, cex.axis=1.25)

#text(max(aes_coefs$Estimate/1.5),which(aes_coefs$type=="Environmental")[1],"Environmental",cex=1.25)
#text(max(aes_coefs$Estimate/1.5),which(aes_coefs$type=="Anthropogenic")[1],"Anthropogenic",cex=1.25)

#mtext(side=4, adj=0, text="Anthropogenic", line=1, font=1, cex=1)
#mtext(side=4, adj=0.75, text="Environmental", line=1, font=1, cex=1)

mtext(side=3, adj=-0.25, text="B", line=1.5, font=2, cex=1.25)

legend("bottomright", legend=c("Environmental","Anthropogenic"),
       bty="n",bg=NA, cex=1, pch=19, col=c(env_color,anthro_color))


###########################
## PLOT THE RESIDUAL MAP ##
###########################

resid_clim <- c(-max(abs(residual_IDW$residual)),max(abs(residual_IDW$residual)))

scatter2D(residual_IDW$lon, residual_IDW$lat, pch=19,
          colvar = residual_IDW$residual, 
          cex=0,
          col=(jet(n=nrow(residual_IDW))),          clim=resid_clim,
          xlab="Longitude", ylab="Latitude", 
          cex.lab=1.5)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = adjustcolor("grey80",alpha=0.2))
scatter2D(residual_IDW$lon, residual_IDW$lat, pch=19,
          colvar = residual_IDW$residual, cex=0.75,
          col=(jet(n=nrow(residual_IDW))),
          clim=resid_clim,
          add=TRUE)
mtext(side=4, line=1.5, text="Residual Value")
mtext(side=3, adj=-0.05, text="C", line=1.5, font=2, cex=1.25)

map(database = "world", fill=T, col="grey80", border="grey80",add=T)
#map(database = "world", fill=T, col="grey70", border="grey70",add=T)

title("Deviation from Expected Given Species Richness",line=1,
      font.main=1, cex.main=1.5)



##########################################
## PLOT THE DEVIATIONA VALUE FOREST PLOT ##
##########################################

plot(deviation_coefs$Estimate, seq(1:nrow(deviation_coefs)), 
     xlim=c(min(aes_coefs$Q5,na.rm = TRUE),max(aes_coefs$Q95,na.rm = TRUE)),
     ylim=c(min(seq(1:nrow(deviation_coefs))-0.25),max(seq(1:nrow(deviation_coefs))+0.25)),cex=0,
     xlab="Standardized Effect Size", ylab=NA, yaxt = "n",
     cex.lab=1.5)
title("Drivers of Deviation Value", line=1,
      font.main=1, cex.main=1.5)

corners <- par("usr")

abline(v=0, lwd=1.5, col=adjustcolor("black",alpha.f = 0.75), lty=2)

par(lend=1)
x0 <- deviation_coefs$Q5
x1 <- deviation_coefs$Q95
y0 <- seq(1:nrow(deviation_coefs))
y1 <- seq(1:nrow(deviation_coefs))
segments(x0,y0,x1,y1, lwd=2)
x0 <- deviation_coefs$Q25
x1 <- deviation_coefs$Q75
y0 <- seq(1:nrow(deviation_coefs))
y1 <- seq(1:nrow(deviation_coefs))
segments(x0,y0,x1,y1, lwd=5)
points(deviation_coefs$Estimate, seq(1:nrow(deviation_coefs)),pch=21,col="white",bg="white",cex=1.75)
points(deviation_coefs$Estimate, seq(1:nrow(deviation_coefs)),pch=21,col=1,cex=1.75,
       bg=ifelse(deviation_coefs$type=="Anthropogenic",anthro_color,env_color))

h_line <- which(aes_coefs$type=="Environmental")[1]-0.5
abline(h=h_line, col=adjustcolor("black",alpha.f = 0.75), lty=1)

axis(2, at = seq(1:nrow(deviation_coefs)), 
     labels = deviation_coefs$variable,
     las=2, cex.axis=1.25)

#text(max(deviation_coefs$Estimate/1.5),which(deviation_coefs$type=="Environmental")[1],"Environmental",cex=1.25)
#text(max(deviation_coefs$Estimate/1.5),which(deviation_coefs$type=="Anthropogenic")[1],"Anthropogenic",cex=1.25)

#mtext(side=4, adj=0, text="Anthropogenic", line=1, font=1, cex=1)
#mtext(side=4, adj=0.75, text="Environmental", line=1, font=1, cex=1)

mtext(side=3, adj=-0.25, text="D", line=1.5, font=2, cex=1.25)

legend("bottomright", legend=c("Environmental","Anthropogenic"),
       bty="n",bg=NA, cex=1, pch=19, col=c(env_color,anthro_color))

#dev.off() 


