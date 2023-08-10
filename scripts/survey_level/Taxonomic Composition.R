


#######################################################################################
#'  CODE TO CALCULATE METRICS OF TAXONOMIC COMPOSITION 
#'  AT THE SITE LEVEL
#'  
#'  dagitty.net/mxttk25
#'
#' @author Matthew McLean, \email {mcleamj@gmail.com},
#'         
#' @date JULY 6, 2023
########################################################################################

######################
## LIBRARY PACKAGES ##
######################

if(!require(factoextra)){install.packages("factoextra"); library(factoextra)}
if(!require(betapart)){install.packages("betapart"); library(betapart)}
if(!require(plot3D)){install.packages("plot3D"); library(plot3D)}
if(!require(pals)){install.packages("pals"); library(pals)}
if(!require(parallel)){install.packages("parallel"); library(parallel)}
if(!require(readr)){install.packages("readr"); library(readr)}
if(!require(tibble)){install.packages("tibble"); library(tibble)}
if(!require(dplyr)){install.packages("dplyr"); library(dplyr)}


family_proportions <- read_rds("outputs/family_proportions.rds")

colnames(family_proportions) <- gsub("family_", "", colnames(family_proportions))

family_proportions <- family_proportions %>%
  rownames_to_column("SiteCode")

model_data <- readRDS("data/model_data.rds") %>%
  select(SiteCode, MPA, SiteLatitude, SiteLongitude, Ecoregion)
model_data <- model_data[!duplicated(model_data$SiteCode),]

family_PCA_data <- merge(model_data, family_proportions,
                        by="SiteCode")

family_arc_sine <- family_PCA_data %>%
  select(-c(SiteCode, MPA, SiteLatitude,SiteLongitude,Ecoregion))
family_arc_sine <- asin(sqrt(family_arc_sine))

family_PCA <- prcomp(family_arc_sine)

taxo_structure <- data.frame(SiteCode = family_PCA_data$SiteCode,
                             family_PCA$x[,1:2]) %>%
  dplyr::rename("Taxo_PC1" = PC1,
                "Taxo_PC2" = PC2)

saveRDS(taxo_structure, "outputs/taxo_structure.rds")

fviz_pca_var(family_PCA, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)

plot(family_PCA$x[,1], family_PCA$x[,2], 
     col=as.integer(as.factor(family_PCA_data$MPA)))
boxplot(family_PCA$x[,1] ~ family_PCA_data$MPA)
boxplot(family_PCA$x[,2] ~ family_PCA_data$MPA)

graphics.off()
plot(family_PCA$x[,1] ~ family_PCA_data$SiteLatitude)

par(mfrow=c(2,2))
for(i in unique(dom_families$family)) {
  plot(family_arc_sine[,i] ~ family_PCA_data$SiteLatitude,
       xlab="Latitude", ylab="Proportion of Species (transformed)")
  title(i)
}



tax_dis <- betapart::beta.pair.abund(family_arc_sine, index.family = "bray")
saveRDS(tax_dis, "outputs/tax_dis.rds")
tax_dis_balanced <- tax_dis$beta.bray.bal
tax_dis_gradient <- tax_dis$beta.bray.gra

total_PCOA <- ape::pcoa(tax_dis$beta.bray)
balanced_PCOA <- ape::pcoa(tax_dis_balanced) #turnover
gradient_PCOA <- ape::pcoa(tax_dis_gradient) #nestedness/abundance driven

saveRDS(total_PCOA, "outputs/total_PCOA.rds")
saveRDS(balanced_PCOA, "outputs/balanced_PCOA.rds")
saveRDS(gradient_PCOA, "outputs/gradient_PCOA.rds")

plot(total_PCOA$vectors[,1] ~ family_PCA_data$SiteLatitude)
plot(balanced_PCOA$vectors[,1] ~ family_PCA_data$SiteLatitude)
plot(gradient_PCOA$vectors[,1] ~ family_PCA_data$SiteLatitude)

plot(family_PCA$x[,1], total_PCOA$vectors[,1])
plot(family_PCA$x[,1], balanced_PCOA$vectors[,1])
plot(family_PCA$x[,1], gradient_PCOA$vectors[,1])




