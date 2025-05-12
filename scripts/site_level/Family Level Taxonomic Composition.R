


#######################################################################################
#'  CODE TO CALCULATE METRICS OF TAXONOMIC COMPOSITION 
#'  AT THE SITE LEVEL
#'  
#'  dagitty.net/mxttk25
#'
#' @author Matthew McLean, \email {mcleamj@gmail.com},
#'         
#' @date Updated May 2024
########################################################################################

######################
## LIBRARY PACKAGES ##
######################

if(!require(ggplot2)){install.packages("ggplot2"); library(ggplot2)}
if(!require(ggrepel)){install.packages("ggrepel"); library(ggrepel)}
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

##################################
## UPLOAD FAMILY ABUNDANCE DATA ##
##################################

family_abundances <- read_rds("outputs/family_abundances.rds") %>%
  select(where(~ any(. != 0)))

#################################################
# UPLOAD MODEL DATA AND MERGE WITH FAMILY DATA ##
#################################################

standardized_data <- readRDS("outputs/standardized_data.rds")

site_survey <- standardized_data %>%
  select(SurveyID, SiteCode)

standardized_data <- standardized_data %>%
  select(SiteCode, MPA, SiteLatitude, SiteLongitude, Ecoregion, aesthe_survey_abund) %>%
  mutate(log_aes = log(aesthe_survey_abund)) %>%
  na.omit()
standardized_data <- standardized_data[!duplicated(standardized_data$SiteCode),]

family_PCA_data <- merge(standardized_data, family_abundances,
                         by="SiteCode") %>%
  arrange(SiteCode)

##################################
## TRANSOFMR RAW ABUNDANCE DATA ##
##################################

family_transform <- family_PCA_data %>%
  select(-c(SiteCode, MPA, SiteLatitude,SiteLongitude,Ecoregion, aesthe_survey_abund, log_aes))

family_log <- log10(family_transform + 1)

family_log_hellinger <- decostand(family_log, method = "hellinger")

############################################
## HOW STRONG IS THE ABUNDANCE GRADIENT ? ##
############################################

plot(rowSums(family_log) ~ family_PCA_data$SiteLatitude) # STRONG

#########################
## PCOA ON BRAY CURTIS ##
#########################

# family_bray <- vegdist(family_log, method = "bray")
# 
# family_pcoa <- pco(family_bray)
# 
# site_scores <- data.frame(scores(family_pcoa, display = "sites"), select(family_PCA_data, SiteCode:Ecoregion)) %>%
#   dplyr::rename(PCoA1 = MDS1,
#                 PCoA2 = MDS2) 
# site_scores$sum_log_abundance <- rowSums(family_log)
# 
# species_fit <- envfit(family_pcoa, family_log, permutations = 999)
# 
# top_species <- as.data.frame(species_fit$vectors$arrows) %>%
#   dplyr::rename(PCoA1 = MDS1,
#                 PCoA2 = MDS2) 
# 
# top_species$R2 <- species_fit$vectors$r
# top_species <- top_species[order(-top_species$R2), ] # Order by R² values
# top_species <- top_species[1:25, ]  # Select top 5 species
# top_species$species <- rownames(top_species)
# 
# ggplot(site_scores, aes(x = PCoA1, y = PCoA2, color=sum_log_abundance)) +
#   geom_point(size = 3) +
#   geom_segment(data = top_species, aes(x = 0, y = 0, xend = PCoA1, yend = PCoA2),
#                arrow = arrow(length = unit(0.2, "cm")), color = "red") +
#   geom_text_repel(data = top_species, aes(x = PCoA1, y = PCoA2, label = species),
#                   color = "red", size = 4) +
#   labs(x = "PCoA Axis 1", y = "PCoA Axis 2",
#        title = "PCoA of Bray-Curtis Distance with Top Species Vectors") +
#   theme_minimal()
# 
# plot(site_scores$PCoA1 ~ site_scores$SiteLatitude)

############################
# PCA ON THE HELLINGER DATA #
############################

family_PCA <- prcomp(family_log_hellinger, scale. = FALSE)
fviz_eig(family_PCA, addlabels = TRUE, ylim = c(0, 100))

fviz_pca_var(family_PCA, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = FALSE)

plot(family_PCA$x[,1] ~ family_PCA_data$SiteLatitude)

# ROTATE THE ORDINATION AXES BASED ON EXPLORATORY ANALYSIS
# THIS WILL PROVIDE A POSITIVE RELATIONSHIP BETWEEN 
# PC AXES AND AESTHETIC VALUE
# WHICH WILL MAKE THE RESULTS MORE INTUITIVE AND CLEAR
# BUT WILL HAVE NO EFFECT ON THE OVERALL RESULT 
# OR THE ORDINATION ITSELF

# FLIP PC1 AND PC2

family_PCA$x[,1] <- family_PCA$x[,1]*-1
family_PCA$rotation[,1] <- family_PCA$rotation[,1] * -1

#family_PCA$x[,2] <- family_PCA$x[,2]*-1
#family_PCA$rotation[,2] <- family_PCA$rotation[,2] * -1

fviz_pca_var(family_PCA, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = FALSE,
             labelsize = 3) 

ggsave(file=here::here("figures_tables/revised_submission/supplementary_figures",
                       "Taxonomic_PCA.png"), 
       device = 'png',
       dpi = 300, 
       units = 'in', 
       width = 8, 
       height = 6, 
       bg = 'white')

# SAVE THE OUTPUT
taxo_structure <- data.frame(SiteCode = family_PCA_data$SiteCode,
                             family_PCA$x[,1:4]) %>%
  dplyr::rename("Taxo_PC1" = PC1,
                "Taxo_PC2" = PC2,
                "Taxo_PC3" = PC3,
                "Taxo_PC4" = PC4)

saveRDS(taxo_structure, "outputs/taxo_structure.rds")


##############################################
## REDUNDANCY ANALYSIS WITH AESTHETIC VALUE ## 
##############################################

# rda_aes <- rda(family_log_hellinger ~ family_PCA_data$log_aes)
# plot(rda_aes)
# 
# rda_scores <- as.data.frame(scores(rda_aes, display = "sites"))
# rda_scores$log_aes <- family_PCA_data$log_aes
# 
# ggplot(rda_scores, aes(x=RDA1, y=PC1, color=log_aes)) +
#   geom_point() +
#   scale_color_gradientn(colors = pals::jet(100))
# 
# taxo_structure <- data.frame(SiteCode = family_PCA_data$SiteCode,
#            Taxo_RDA1 = rda_scores$RDA1) 
# 
# saveRDS(taxo_structure, "outputs/taxo_structure.rds")


#########################################################
## BETA PART BRAY CURTIS DECOMPOSITION AT FAMILY LEVEL ##
#########################################################

# tax_dis <- betapart::beta.pair.abund(family_log, index.family = "bray")
# saveRDS(tax_dis, "outputs/tax_dis.rds")
# tax_dis_turnover <- tax_dis$beta.bray.bal
# tax_dis_gradient <- tax_dis$beta.bray.gra
# 
# total_PCOA <- ape::pcoa(tax_dis$beta.bray)
# turnover_PCOA <- ape::pcoa(tax_dis_turnover) #turnover
# gradient_PCOA <- ape::pcoa(tax_dis_gradient) #nestedness/abundance driven
# 
# saveRDS(total_PCOA, "outputs/total_PCOA.rds")
# saveRDS(turnover_PCOA, "outputs/turnover_PCOA")
# saveRDS(gradient_PCOA, "outputs/gradient_PCOA.rds")
# 
# plot(total_PCOA$vectors[,1] ~ family_PCA_data$SiteLatitude)
# plot(turnover_PCOA$vectors[,1] ~ family_PCA_data$SiteLatitude)
# plot(gradient_PCOA$vectors[,1] ~ family_PCA_data$SiteLatitude)
# 
# plot(total_PCOA$vectors[,1], turnover_PCOA$vectors[,1])
# plot(total_PCOA$vectors[,1], gradient_PCOA$vectors[,1])
# 
# plot(family_PCA$x[,1], total_PCOA$vectors[,1])
# cor.test(family_PCA$x[,1], total_PCOA$vectors[,1])
# 
# plot(family_PCA$x[,1], turnover_PCOA$vectors[,1])
# cor.test(family_PCA$x[,1], turnover_PCOA$vectors[,1])
# 
# plot(family_PCA$x[,1], gradient_PCOA$vectors[,1])
# cor.test(family_PCA$x[,1], gradient_PCOA$vectors[,1])

##
## 
## 

