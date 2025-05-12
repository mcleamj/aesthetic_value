
library(brms)
library(bayesplot)
library(readr)
library(dplyr)

# LOAD MPA MECHANSIM EFFECT SIZE DATA 

MPA_models <- read_rds("outputs/MPA_mechanism_effects.rds")
MPA_coefs <- as.data.frame(posterior_summary(MPA_models, 
                                             probs=c(0.05,0.25,0.75,0.95)))

MPA_coefs$variable <- rownames(MPA_coefs)

MPA_coefs <- MPA_coefs %>%
  mutate(variable=recode(variable, 
                     "tax_entropy" = "Taxonomic Diversity",
                     "phylo_entropy" = "Phylogenetic Diversity",
                     "Taxonomic_Composition" = "Taxonomic Composition",
                     "fun_entropy" = "Functional Diversity",
                     "Benthic_Composition" = "Benthic Composition",
                     "Trophic_Composition" = "Trophic Composition",
                     "MPA.Total.Causal.Effect"="MPA Total Causal Effect"))

MPA_coefs$abs_effect <- abs(MPA_coefs$Estimate)

MPA_coefs <- MPA_coefs %>%
  arrange(abs_effect)

MPA_coefs$sign <- ifelse(MPA_coefs$Estimate>0,"pos","neg")

MPA_coefs$color <- ifelse(MPA_coefs$sign=="pos","blue","red")

MPA_coefs$color_alpha <- NA

for(i in 1:nrow(MPA_coefs)){
  MPA_coefs$color_alpha[i] <- adjustcolor(MPA_coefs$color[i], 
                                          alpha.f = abs(MPA_coefs$Estimate[i])*50)
}

MPA_coefs$size <- rep(2)

MPA_coefs$color[MPA_coefs$variable=="MPA Total Causal Effect"] <- "black"
MPA_coefs$color_alpha[MPA_coefs$variable=="MPA Total Causal Effect"] <- 
  adjustcolor('black', alpha.f = 1)

MPA_coefs$size[MPA_coefs$variable=="MPA Total Causal Effect"] <- 2.5


###########################
## PLOT THE  FOREST PLOT ##
###########################

graphics.off()

par(mar=c(4,16,2,4))

plot(MPA_coefs$Estimate, seq(1:nrow(MPA_coefs)), xlim=c(min(MPA_coefs$Q5,na.rm = TRUE),max(MPA_coefs$Q95,na.rm = TRUE)),
     ylim=c(min(seq(1:nrow(MPA_coefs))-0.25),max(seq(1:nrow(MPA_coefs))+0.25)),cex=0,
     xlab="Standardized Effect Size", ylab=NA, yaxt = "n",
     cex.lab=1.25)
title("Contributions to MPA Effect", line=1,
      font.main=1, cex.main=1.2)

corners <- par("usr")

abline(v=0, lwd=1.5, col=adjustcolor("black",alpha.f = 0.75), lty=2)

par(lend=1)
x0 <- MPA_coefs$Q5
x1 <- MPA_coefs$Q95
y0 <- seq(1:nrow(MPA_coefs))
y1 <- seq(1:nrow(MPA_coefs))
segments(x0,y0,x1,y1, lwd=MPA_coefs$size,col=MPA_coefs$color_alpha)
x0 <- MPA_coefs$Q25
x1 <- MPA_coefs$Q75
y0 <- seq(1:nrow(MPA_coefs))
y1 <- seq(1:nrow(MPA_coefs))
segments(x0,y0,x1,y1, lwd=MPA_coefs$size+3, col=MPA_coefs$color_alpha)
points(MPA_coefs$Estimate, seq(1:nrow(MPA_coefs)),pch=21,col="white",bg="white",cex=MPA_coefs$size)
points(MPA_coefs$Estimate, seq(1:nrow(MPA_coefs)),pch=21,col=1,cex=MPA_coefs$size,
       bg=MPA_coefs$color_alpha)

h_line <- which(MPA_coefs$type=="Environmental")[1]-0.5
abline(h=h_line, col=adjustcolor("black",alpha.f = 0.75), lty=1)

axis(2, at = seq(1:nrow(MPA_coefs)), 
     labels = MPA_coefs$variable,
     las=2, cex.axis=1)

abline(h=nrow(MPA_coefs)-0.5, lty=1, lwd=1)





