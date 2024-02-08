

MPA_models <- read_rds("outputs/MPA_mechanism_effects.rds")
MPA_coefs <- as.data.frame(posterior_summary(MPA_models, 
                                             probs=c(0.05,0.25,0.75,0.95)))

MPA_coefs$variable <- rownames(MPA_coefs)

MPA_coefs <- MPA_coefs %>%
  mutate(variable=recode(variable, 
                     "tax_richness" = "Taxonomic Richness",
                     "phylo_richness" = "Phylogenetic Richness",
                     "Taxo_PC2" = "Taxonomic Composition (PC2)",
                     "Taxo_PC1" = "Taxonomic Composition (PC1)",
                     "fun_richness" = "Functional Richness",
                     "PC2_Benthic" = "Benthic Composition (PC2)",
                     "PC1_Benthic" = "Benthic Composition (PC1)",
                     "PC1_Trophic" = "Trophic Composition (PC1)",
                     "PC2_Trophic" = "Trophic Composition (PC2)",
                     "MPA.Total.Causal.Effect"="MPA Total Causal Effect"))

MPA_coefs$abs_effect <- abs(MPA_coefs$Estimate)

MPA_coefs <- MPA_coefs %>%
  arrange(abs_effect)

MPA_coefs$sign <- ifelse(MPA_coefs$Estimate>0,"pos","neg")

MPA_coefs$color <- ifelse(MPA_coefs$sign=="pos","blue","red")

MPA_coefs$size <- rep(2)

MPA_coefs$color[MPA_coefs$variable=="MPA Total Causal Effect"] <- "black"

MPA_coefs$size[MPA_coefs$variable=="MPA Total Causal Effect"] <- 3


###########################
## PLOT THE  FOREST PLOT ##
###########################


par(mar=c(4,16,4,4))

plot(MPA_coefs$Estimate, seq(1:nrow(MPA_coefs)), xlim=c(min(MPA_coefs$Q5,na.rm = TRUE),max(MPA_coefs$Q95,na.rm = TRUE)),
     ylim=c(min(seq(1:nrow(MPA_coefs))-0.25),max(seq(1:nrow(MPA_coefs))+0.25)),cex=0,
     xlab="Standardized Effect Size", ylab=NA, yaxt = "n",
     cex.lab=1.5)
title("Contributions to MPA Effect", line=1,
      font.main=1, cex.main=1.5)

corners <- par("usr")

par(lend=1)
x0 <- MPA_coefs$Q5
x1 <- MPA_coefs$Q95
y0 <- seq(1:nrow(MPA_coefs))
y1 <- seq(1:nrow(MPA_coefs))
segments(x0,y0,x1,y1, lwd=MPA_coefs$size,col=MPA_coefs$color)
x0 <- MPA_coefs$Q25
x1 <- MPA_coefs$Q75
y0 <- seq(1:nrow(MPA_coefs))
y1 <- seq(1:nrow(MPA_coefs))
segments(x0,y0,x1,y1, lwd=MPA_coefs$size+3, col=MPA_coefs$color)
points(MPA_coefs$Estimate, seq(1:nrow(MPA_coefs)),pch=21,col="white",bg="white",cex=MPA_coefs$size)
points(MPA_coefs$Estimate, seq(1:nrow(MPA_coefs)),pch=21,col=1,cex=MPA_coefs$size,
       bg=MPA_coefs$color)

abline(v=0, lwd=1.5, col=adjustcolor("black",alpha.f = 0.75), lty=2)
h_line <- which(MPA_coefs$type=="Environmental")[1]-0.5
abline(h=h_line, col=adjustcolor("black",alpha.f = 0.75), lty=1)

axis(2, at = seq(1:nrow(MPA_coefs)), 
     labels = MPA_coefs$variable,
     las=2, cex.axis=1.25)

abline(h=nrow(MPA_coefs)-0.5, lty=2, lwd=2)

