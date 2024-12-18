

#################################################
## CODE TO EXAMINE LOCAL TAXONOMIC DIFFERENCES ##
#################################################

######################
## LIBRARY PACKAGES ##
######################

if(!require(pals)){install.packages("pals"); library(pals)}
if(!require(vegan)){install.packages("vegan"); library(vegan)}
if(!require(factoextra)){install.packages("factoextra"); library(factoextra)}
if(!require(FactoMineR)){install.packages("FactoMineR"); library(FactoMineR)}
if(!require(ggrepel)){install.packages("ggrepel"); library(ggrepel)}
if(!require(tidyverse)){install.packages("tidyverse"); library(tidyverse)}

#############################################
## IMPORT PREPARED STANDARDIZED MODEL DATA ##
#############################################

standardized_data <-  read_rds("outputs/standardized_data.rds")
standardized_data$aes_deviation <- standardized_data$aesthe_survey_abund - standardized_data$aesthe_SR_survey

# FILTER TO GREAT BARRIER REEF
# Central and Southern Great Barrier Reef
# Torres Strait Northern Great Barrier Reef

standardized_data <- standardized_data %>%
  filter(Ecoregion == "Torres Strait Northern Great Barrier Reef") %>%
  #filter(Ecoregion == "Central and Southern Great Barrier Reef")
  select(SiteCode, aesthe_survey_abund, aes_deviation) %>%
  mutate(log_aes = log(aesthe_survey_abund)) %>%
  group_by(SiteCode) %>%
  summarise_all(.funs=mean)

## FILTER TO MOST BEAUTIFUL AND MOST UGLY SURVEYS
standardized_data <- standardized_data %>%
  filter(log_aes >= quantile(log_aes, probs=0.75) |
           log_aes <= quantile(log_aes, probs=0.25))


##################################
## UPLOAD FAMILY ABUNDANCE DATA ##
##################################

family_abundances <- read_rds("outputs/family_abundances.rds")%>%
  filter(SiteCode %in% standardized_data$SiteCode)  %>%
  column_to_rownames("SiteCode") %>%
  select_if(colSums(.) != 0)

############################
## LOG HELLINGER THEN PCA ##
############################

abund_log <- log10(family_abundances + 1)
log_hellinger <- decostand(abund_log, method = "hellinger")

#############################
# PCA ON THE HELLINGER DATA #
#############################

# Run PCA using vegan's rda function
pca_result <- rda(log_hellinger, scale = FALSE)  # Scale data for PCA

# Extract PCA site scores (coordinates for points)
site_scores <- as.data.frame(scores(pca_result, display = "sites"))
site_scores$ColorVar <- standardized_data$log_aes
# Extract PCA species scores (coordinates for vectors)
species_scores <- as.data.frame(scores(pca_result, display = "species"))
species_scores$Species <- rownames(species_scores)

# Create a ggplot
ggplot(site_scores, aes(x = PC1, y = PC2, color = ColorVar)) +
  geom_point(size = 3) +  # Add points
  geom_segment(data = species_scores, 
               aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +  # Add species vectors
  geom_text_repel(data = species_scores, 
                  aes(x = PC1, y = PC2, label = Species), 
                  size = 4, color = "black") +  # Label species vectors
  labs(x = "PC1", y = "PC2", color = "Log(Aesthetic Value)") +  # Axis labels
  scale_color_gradientn(colors = jet(100)) +
  #theme_minimal() +  # Minimal theme
  theme(legend.position = "top")  # Move legend to the top

ggsave("GBR_local.png", 
       dpi=300, 
       units="in",
       height=6,
       width=8)
