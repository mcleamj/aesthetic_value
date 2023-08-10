
# test trophic and mpa

test <- model_data %>%
  select(SurveyID, aesthe_survey)

test <- merge(model_data, trophic_structure, by="SurveyID")

boxplot(test$Diet_Grazer ~ test$MPA)

# DO TROPHIC GROUPS DIFFER BETWEEN MPA AND FISHED
# DIRICHLET MODEL?
# MODEL PC1 AND PC2

test_lm <- lme4::lmer(PC1_trophic ~ MPA +

                        as.factor(Temperature_Zone) +
                        (1 | Country/SiteCode),
                      data=standardized_data)

# DOES AESTHETIC VALUE CHANGE WITH TROPHIC LEVEL?

# POSITIVE RELATIONSHIP WITH PC1 AND PC2 - COOL

test_lm <- lme4::lmer(log(aesthe_survey) ~ PC1_trophic +
                        PC1_imputed +
                        PC2_imputed +
                        Biomass +
                        fun_entropy +
                        phylo_entropy + 
                        taxo_entropy +
                        
       as.factor(Temperature_Zone) +
       (1 | Country/SiteCode),
       data=standardized_data)

summary(test_lm)
visreg::visreg(test_lm, "PC1_trophic")
