###################################################################################################
#'  CALCULATE COMMUNITY WEIGHTED MEAN AND SUM OF AESTHETIC VALUE
#'  AND COMPARE VS. STANDARD AESTHETIC VALUE CALCULATION
#'
#'        
#' @date
##################################################################################################

rm(list = ls())
library(FD)
library(tidyverse)

# IMPORT NECESSARY DATA 

aesthe_species <- read.csv2(here::here("outputs", "aesthe_species.csv"))
survey_aesth <- read.csv(here::here("outputs", "survey_aesth.csv"))
abund_matrix <- readRDS(here::here("outputs", "sp_abund_matrix.rds"))

# COMMUNITY WEIDGHTED MEAN

aesthe_species <- aesthe_species %>%
  mutate(sp_name = gsub(" ", "_", aesthe_species$sp_name)) %>% 
  filter(sp_name %in% colnames(abund_matrix)) %>%
  column_to_rownames("sp_name") %>%
  select(aesthe_score)

abund_matrix <- abund_matrix %>%
  column_to_rownames("SurveyID")

log_abundance <- log1p(abund_matrix)

cwm_aes <- functcomp(as.matrix(aesthe_species), as.matrix(log_abundance))

# COMMUNITY WEIGHTED SUM 

# transpose abundance table

t_log_abundance <- as.data.frame(t(log_abundance))

identical(rownames(t_log_abundance), rownames(aesthe_species))

com_sum_aes <- t_log_abundance * aesthe_species$aesthe_score 
com_sum_aes <- as.data.frame(colSums(com_sum_aes)) 
colnames(com_sum_aes) <- "aesthe_score"

# SAVE OUTPUTS
com_sum_aes <- com_sum_aes %>%
  rownames_to_column("SurveyID")
saveRDS(com_sum_aes, "outputs/com_sum_aes.rds")

cwm_aes <- cwm_aes %>%
  rownames_to_column("SurveyID")
saveRDS(cwm_aes, "outputs/cwm_aes.rds")

# COMPARE THE 3 METHODS

par(mfrow=c(1,2))

plot(survey_aesth$aesthe_survey_abund, cwm_aes$aesthe_score,
     xlab="Abundance-Weighted Aesthetic Value",
     ylab="Community Weighted Mean Aesthetic Value")
plot(survey_aesth$aesthe_survey_abund, com_sum_aes$aesthe_score,
     xlab="Abundance-Weighted Aesthetic Value",
     ylab="Community Weighted Sum Aesthetic Value")

## GGPLOT FOR SUPP FIGURE




