

library(DescTools)
library(brms)
library(bayesplot)
library(readr)
library(dplyr)

## IMPORT EFFECT SIZES

model_outputs <- read_rds("outputs/MPA_mechanism_effects.rds")
path_effects <- as.data.frame(posterior_summary(model_outputs))
path_effects$variable <- rownames(path_effects)

mcmc_intervals(model_outputs)

# SUM THE EFFECTS OF PC1 AND PC2 TOGETHER FOR 
# VARIABLES BASED ON ORDINATIONS

tax_div_effect <- path_effects %>%
  filter(variable=="tax_richness") %>%
  select(Estimate) %>%
  as.numeric()

phylo_div_effect <- path_effects %>%
  filter(variable=="phylo_richness") %>%
  select(Estimate) %>%
  as.numeric()

tax_comp_effect <- path_effects %>%
  filter(variable=="Taxo_PC1" |
           variable=="Taxo_PC2") %>%
  select(Estimate) %>%
  summarise(sum(Estimate)) %>%
  as.numeric()

trophic_comp_effect <- path_effects %>%
  filter(variable=="PC1_Trophic" |
           variable=="PC2_Trophic" ) %>%
  select(Estimate) %>%
  summarise(sum(Estimate)) %>%
  as.numeric()

benthic_comp_effect <- path_effects %>%
  filter(variable=="PC1_Benthic" |
           variable=="PC2_Benthic" ) %>%
  select(Estimate) %>%
  summarise(sum(Estimate)) %>%
  as.numeric()


# PLOT SETTINGS

x_range <- c(0,12)
y_range <- c(0,6)

x_rad_big <- max(x_range)/9
y_rad_big <- max(y_range)/9

x_rad_small <- x_rad_big*0.9
y_rad_small <- y_rad_big*0.9

x_center <- max(x_range)/2
y_center <- max(y_range)/2

x_left <- min(x_range)+1
x_right <- max(x_range)-1

y_bottom <- min(y_range)+0.5
y_top <- max(y_range)-0.5
y_bottom_mid <- (y_bottom + y_center)/2
y_top_mid <- (y_top + y_center)/2

## MAKE FIGURE

plot(NA, xlim=x_range,ylim=y_range, xlab=NA,ylab=NA,
     xaxt="n", yaxt="n")

# MPA
DrawEllipse(x = x_left, y = y_center, radius.x = x_rad_big, radius.y = y_rad_big, rot = 0,
            nv = 100, border = "black", col = "grey",
            lwd = 1.5)
text(x_left, y_center, "MPA",
     font=2, col="black", cex=1.25)

# TAX DIVERSITY
DrawEllipse(x = x_center, y = y_top, radius.x = x_rad_small, radius.y = y_rad_small, rot = 0,
            nv = 100, border = "black", col = "grey",
            lwd = 1.5)
text(x_center, y_top, "Taxonomic\nDiversity",
     font=2, col="black")

# PHYLO DIVERSITY
DrawEllipse(x = x_center, y = y_top_mid, radius.x = x_rad_small, radius.y = y_rad_small, rot = 0,
            nv = 100, border = "black", col = "grey",
            lwd = 1.5)
text(x_center, y_top_mid, "Phylogenetic\nDiversity",
     font=2, col="black")

# TAXO COMPOSITION
DrawEllipse(x = x_center, y = y_center, radius.x = x_rad_small, radius.y = y_rad_small, rot = 0,
            nv = 100, border = "black", col = "grey",
            lwd = 1.5)
text(x_center, y_center, "Taxonomic\nComposition",
     font=2, col="black")

# TROPHIC COMPOSITION
DrawEllipse(x = x_center, y = y_bottom_mid, radius.x = x_rad_small, radius.y = y_rad_small, rot = 0,
            nv = 100, border = "black", col = "grey",
            lwd = 1.5)
text(x_center, y_bottom_mid, "Trophic\nComposition",
     font=2, col="black")

# BENTHIC COMPOSITION
DrawEllipse(x = x_center, y = y_bottom, radius.x = x_rad_small, radius.y = y_rad_small, rot = 0,
            nv = 100, border = "black", col = "grey",
            lwd = 1.5)
text(x_center, y_bottom, "Benthic\nComposition",
     font=2, col="black")

# AESTHETIC
DrawEllipse(x = x_right, y = y_center, radius.x = x_rad_big, radius.y = y_rad_big, rot = 0,
            nv = 100, border = "black", col = "grey",
            lwd = 1.5)
text(x_right, y_center, "Aesthetic\nValue", col="black",
     font=2, cex=1.25)

## ARROWS

# MPA ARROWS
arrows(x0=x_left+x_rad_big, x1=x_center-x_rad_big,
       y0=y_center+0.2, y1=y_top,
       length=0.1,
       lwd=2, lty=2)

arrows(x0=x_left+x_rad_big, x1=x_center-x_rad_big,
       y0=y_center+0.1, y1=y_top_mid,
       length=0.1,
       lwd=2, lty=2)

arrows(x0=x_left+x_rad_big, x1=x_center-x_rad_big,
       y0=y_center, y1=y_center,
       length=0.1,
       lwd=2, lty=2)

arrows(x0=x_left+x_rad_big, x1=x_center-x_rad_big,
       y0=y_center-0.1, y1=y_bottom_mid,
       length=0.1,
       lwd=2, lty=2)

arrows(x0=x_left+x_rad_big, x1=x_center-x_rad_big,
       y0=y_center-0.2, y1=y_bottom,
       length=0.1,
       lwd=2, lty=2)

# EFFECT SIZE ARROWS

# MULTIPLE ALL EFFECTS BY A CONSTANT 
# TO INCREASE SIZE FOR PLOT

arrow_constant <- 100

# TAX DIV
arrows(x0=x_center+x_rad_big, x1=x_right-x_rad_big, 
       y0=y_top, y1=y_center+0.2,
       length=0.1,
       lwd=abs(tax_div_effect*arrow_constant), lty=1,
       col=ifelse(tax_div_effect>0,"blue","red"))

# PHYLO DIV
arrows(x0=x_center+x_rad_big, x1=x_right-x_rad_big, 
       y0=y_top_mid, y1=y_center+0.1,
       length=0.1,
       lwd=abs(phylo_div_effect*arrow_constant), lty=1,
       col=ifelse(phylo_div_effect>0,"blue","red"))

# TAX COMPOSITION
arrows(x0=x_center+x_rad_big, x1=x_right-x_rad_big,
       y0=y_center, y1=y_center,
       length=0.1,
       lwd=abs(tax_comp_effect*arrow_constant), lty=1,
       col=ifelse(tax_comp_effect>0,"blue","red"))


# TROPHIC COMP
arrows(x0=x_center+x_rad_big, x1=x_right-x_rad_big, 
       y0=y_bottom_mid, y1=y_center-0.1,
       length=0.1,
       lwd=abs(trophic_comp_effect*arrow_constant),  lty=1,
       col=ifelse(trophic_comp_effect>0,"blue","red"))

# BENTHIC COMP
arrows(x0=x_center+x_rad_big, x1=x_right-x_rad_big, 
       y0=y_bottom, y1=y_center-0.2,
       length=0.1,
       lwd=abs(benthic_comp_effect*arrow_constant), lty=1,
       col=ifelse(benthic_comp_effect>0,"blue","red"))

legend("bottomright",
       pch=19, col=c("blue","red"),
       legend=c("Positive Effect", "Negative Effect"),
       border=NA,
       cex=1)
