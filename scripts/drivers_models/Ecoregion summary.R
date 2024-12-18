

library(dplyr)

survey_esth <- read.csv("outputs/survey_aesth.csv")

survey_esth$residual_abund <- survey_esth$aesthe_survey_abund - survey_esth$aesthe_SR_survey

survey_esth$residual_pres <- survey_esth$aesthe_survey_pres - survey_esth$aesthe_SR_survey

model_data <- readRDS("data/model_data.rds")

model_data <- model_data %>%
  select(SurveyID, SiteCode, Ecoregion)

survey_esth <- merge(survey_esth, model_data, by="SurveyID")

# MEAN BY SITE, THEN ECOREGION

eco_mean <- survey_esth %>%
  group_by(SiteCode, Ecoregion) %>%
  summarise_all(.funs=mean) %>%
  ungroup() %>%
  select(-SiteCode) %>%
  group_by(Ecoregion) %>%
  summarise_all(.funs=mean)

eco_mean

plot(eco_mean$aesthe_survey_abund ~ eco_mean$nb_species)
points(eco_mean$nb_species, eco_mean$aesthe_SR_survey,col=2)

cor.test(eco_mean$aesthe_survey_abund, eco_mean$residual_abund)


