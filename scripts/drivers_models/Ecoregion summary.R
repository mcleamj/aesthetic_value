

library(dplyr)

survey_esth <- read.csv("outputs/survey_aesth.csv")

survey_esth$residual <- survey_esth$aesthe_survey - survey_esth$aesthe_SR_survey

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

cor.test(eco_mean$aesthe_survey, eco_mean$residual)


