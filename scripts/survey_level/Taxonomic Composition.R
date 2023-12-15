


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

if(!require(ape)){install.packages("ape"); library(ape)}
if(!require(vegan)){install.packages("vegan"); library(vegan)}
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

model_data <- readRDS("data/model_data.rds")

site_survey <- model_data %>%
  select(SurveyID, SiteCode)

model_data <- model_data %>%
  select(SiteCode, MPA, SiteLatitude, SiteLongitude, Ecoregion)
model_data <- model_data[!duplicated(model_data$SiteCode),]

family_PCA_data <- merge(model_data, family_proportions,
                        by="SiteCode") %>%
  arrange(SiteCode)

###############################################
## ARC SINE TRANSOFMRATON OF PROPORTION DATA ##
###############################################

family_arc_sine <- family_PCA_data %>%
  select(-c(SiteCode, MPA, SiteLatitude,SiteLongitude,Ecoregion))
family_arc_sine <- asin(sqrt(family_arc_sine))

############################
# PCA ON THE ARC SINE DATA #
############################

family_PCA <- prcomp(family_arc_sine, scale. = FALSE)

taxo_structure <- data.frame(SiteCode = family_PCA_data$SiteCode,
                             family_PCA$x[,1:2]) %>%
  dplyr::rename("Taxo_PC1" = PC1,
                "Taxo_PC2" = PC2)

saveRDS(taxo_structure, "outputs/taxo_structure.rds")

fviz_pca_var(family_PCA, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = FALSE)

# ROTATE AXIS 1 SO THAT IT IS POSITIVELY
# ASSOCIATED WITH AESTHETIC VALUE
# TO BE MORE INTUITIVE
# RE PLOT THE PCA
family_PCA$x[,1] <- family_PCA$x[,1]*-1
family_PCA$rotation[,1] <- family_PCA$rotation[,1] * -1
fviz_pca_var(family_PCA, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = FALSE)



plot(family_PCA$x[,1], family_PCA$x[,2], 
     col=as.integer(as.factor(family_PCA_data$MPA)))
boxplot(family_PCA$x[,1] ~ family_PCA_data$MPA)
boxplot(family_PCA$x[,2] ~ family_PCA_data$MPA)

graphics.off()
plot(family_PCA$x[,1] ~ family_PCA_data$SiteLatitude)

##################################
## FAMILY LEVEL LATITUDE TRENDS ##
##################################

par(mfrow=c(2,2))
for(i in unique(dom_families$family)) {
  plot(family_arc_sine[,i] ~ family_PCA_data$SiteLatitude,
       xlab="Latitude", ylab="Proportion of Species (transformed)")
  title(i)
}

###############################################
## BRAY CURTIS AND PCOA ON THE ARC SINE DATA ##
###############################################

bray_family <- vegdist(family_arc_sine)

family_pco <- pcoa(bray_family)

biplot.pcoa(family_pco)

# RELATIONSHIP BETWEEN PCA AND PCOA ?

plot(family_PCA$x[,1], family_pco$vectors[,1])
cor.test(family_PCA$x[,1], family_pco$vectors[,1])
# 99%

plot(family_PCA$x[,2], family_pco$vectors[,2])
cor.test(family_PCA$x[,2], family_pco$vectors[,2])
# 93%

#########################################################
## BETA PART BRAY CURTIS DECOMPOSITION AT FAMILY LEVEL ##
#########################################################

# tax_dis <- betapart::beta.pair.abund(family_arc_sine, index.family = "bray")
# saveRDS(tax_dis, "outputs/tax_dis.rds")
# tax_dis_balanced <- tax_dis$beta.bray.bal
# tax_dis_gradient <- tax_dis$beta.bray.gra
# 
# total_PCOA <- ape::pcoa(tax_dis$beta.bray)
# balanced_PCOA <- ape::pcoa(tax_dis_balanced) #turnover
# gradient_PCOA <- ape::pcoa(tax_dis_gradient) #nestedness/abundance driven
# 
# saveRDS(total_PCOA, "outputs/total_PCOA.rds")
# saveRDS(balanced_PCOA, "outputs/balanced_PCOA.rds")
# saveRDS(gradient_PCOA, "outputs/gradient_PCOA.rds")
# 
# plot(total_PCOA$vectors[,1] ~ family_PCA_data$SiteLatitude)
# plot(balanced_PCOA$vectors[,1] ~ family_PCA_data$SiteLatitude)
# plot(gradient_PCOA$vectors[,1] ~ family_PCA_data$SiteLatitude)
# 
# plot(family_PCA$x[,1], total_PCOA$vectors[,1])
# plot(family_PCA$x[,1], balanced_PCOA$vectors[,1])
# plot(family_PCA$x[,1], gradient_PCOA$vectors[,1])

################################################
## WHAT ABOUT BETA PART JACCARD?
## APPLIED TO SPECIES PRESENCE ABSENCE DATA ?
################################################

species_data <- readRDS("outputs/sp_pres_matrix.rds") 

species_data <- merge(site_survey, species_data, by="SurveyID")
species_data <- species_data %>%
  select(-SurveyID) %>%
  group_by(SiteCode) %>%
  summarise_all(.funs=sum) %>%
  arrange(SiteCode)

identical(family_PCA_data$SiteCode, species_data$SiteCode)

species_matrix <- species_data %>%
  select(-SiteCode)
species_matrix[species_matrix>=1] <- 1

raw_jaccard <- vegdist(species_matrix, method = "jaccard")

raw_jaccard_pco <- pcoa(raw_jaccard)

plot(raw_jaccard_pco$vectors[,1] ~ family_PCA_data$SiteLatitude)

plot(family_PCA$x[,1], raw_jaccard_pco$vectors[,1])
cor.test(family_PCA$x[,1], raw_jaccard_pco$vectors[,1])

plot(family_PCA$x[,2], raw_jaccard_pco$vectors[,2])
cor.test(family_PCA$x[,2], raw_jaccard_pco$vectors[,2])

beta_jac <- betapart::beta.pair(species_matrix, index.family = "jaccard")
jac_total <- beta_jac$beta.jac
jac_turnover <- beta_jac$beta.jtu
jac_nested <- beta_jac$beta.jne

jac_pco <- pcoa(jac_nested)

plot(jac_pco$vectors[,1] ~ family_PCA_data$SiteLatitude)

plot(family_PCA$x[,1], jac_pco$vectors[,1])
cor.test(family_PCA$x[,1], jac_pco$vectors[,1])

biplot.pcoa(raw_jaccard_pco)

