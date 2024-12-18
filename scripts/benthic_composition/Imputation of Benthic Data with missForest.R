
#######################################################################################
#'  SCRIPT TO PEFORM RANDOM-FOREST BASED IMPUTATION OF BENTHIC DATA 
#'  FROM REEF LIFE SURVEY
#'
#' @author Matthew McLean, \email {mcleamj@gmail.com},
#'         Eva Maire
#'         
#' @date updated May 2024
#' 
########################################################################################

######################
#' LIBRARY PACKAGES ##
######################

library(tidyr)
library(vegan)
library(ape)
library(funrar)
library(FactoMineR)
library(factoextra)
library(corrgram)
library(corrplot)
library(RColorBrewer)
library(patchwork)
library(missForest)
library(dplyr)

########################################
#' LOAD SURVEY INFO AND BENTHIC DATA  ##
########################################

RLS_survey_info <- read.table("data/RLS_survey_info.txt")

RLS <- read.table("data/benthic_survey.txt", header=T)

###################
## DATA CLEANING ## 
###################

RLS$No.corals.kelp.scored <- NULL #remove empty column 

colnames(RLS) <- gsub("_Bleached", "", colnames(RLS)) # remove details on bleaching

category_table <- read.csv("data/Benthic_Categories.csv")

category_table$RLS.cat <- gsub("\\(", "\\.", category_table$RLS.cat)
category_table$RLS.cat <- gsub("\\)", "\\.", category_table$RLS.cat)
category_table$RLS.cat <- gsub(" ", "\\.", category_table$RLS.cat)
category_table$RLS.cat <- gsub("<", "\\.", category_table$RLS.cat)
category_table$RLS.cat <- gsub("/", "\\.", category_table$RLS.cat)
category_table$RLS.cat <- gsub("-", "\\.", category_table$RLS.cat)
category_table$RLS.cat <- gsub(",", "\\.", category_table$RLS.cat)

# RE-CLASSIFY RLS CATEGORIES TO 9 SIMPLIFIED CATEGORIES 
category_table <- category_table %>% mutate(Broad_Groups_EM = recode(Broad_Groups, 
                                                                     "ahermatypic coral"="coral",
                                                                     "Corymbose" = "coral",
                                                                     "Tabular"= "coral" ,
                                                                     "Branching" = "coral",
                                                                     "Encrusting"  = "coral",
                                                                     "Laminar" = "coral",
                                                                     "Hemispherical"= "coral",
                                                                     "coralline algae" = "coralline algae",
                                                                     "encrusting algae" = "microalgal_mats",
                                                                     "fleshy algae" = "algae",
                                                                     "Halimeda" = "algae",
                                                                     "canopy forming macroalgae"= "algae",
                                                                     "understory macroalgae"= "algae",
                                                                     "slime" = "microalgal_mats",
                                                                     "turf algae" = "algae",
                                                                     "seagrass" = "seagrass"))


# REMOVE OTHER CATEGORIES
colnames(RLS)[!colnames(RLS) %in% category_table$RLS.cat]

# AGGREGATE CATEGORIES
RLS_broad <- RLS
colnames(RLS_broad) <- category_table$Broad_Groups_EM[match(colnames(RLS_broad), category_table$RLS.cat)]
RLS_broad_names <- colnames(RLS_broad)
RLS_broad <- t(RLS_broad)
RLS_broad <- aggregate(RLS_broad, by=list(RLS_broad_names), FUN=sum)
rownames(RLS_broad) <- RLS_broad$Group.1; RLS_broad$Group.1 <- NULL
RLS_broad <- as.data.frame(t(RLS_broad))
dim(RLS_broad) #6515 transects

# SELECT TROPICAL TRANSECTS ONLY
#import list from N. LOISEAU - subset based on the temperature (17 degrees)
load("data/selected_surveys.RData")
selected_survey_info <- selected_surveys %>%
  select(SurveyID, SiteLongitude, SiteLatitude, SiteCode)
tropical <- selected_surveys
dim(tropical) #3784

RLStropical <- RLS_broad[which((rownames(RLS_broad) %in% tropical$SurveyID)==T),] 
dim(RLStropical) #2438 tropical transects scored 

# More cleaning: remove empty rows
length(which(rowSums(RLStropical)==0)) 
RLStropical <- RLStropical[-which(rowSums(RLStropical)==0),] # ONE TRANSECT WITH NO DATA 
RLStropical[is.na(RLStropical)] <- 0 #replace NAs
RLStropical <- make_relative(as.matrix(RLStropical)) #make % coverage
RLStropical <- as.data.frame(RLStropical)

###########################################
## APPLY A PRINCIPAL COMPONENET ANALYSIS ##
###########################################

## FIRST PERFORM ARC-SIN TRANSFORMATION FOR RIGHT-SKEWED PROPORTIONAL DATA

#transformation => ARCSIN
RLStropical_arc <- asin(sqrt(RLStropical))
RLStropical_PCA <- prcomp(RLStropical_arc,scale = F) # Run PCA with no scaling 
fviz_pca_var(RLStropical_PCA, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)


############################################################
## ROTATE THE FIRST PC AXIS 
## FURTHER ANALYSIS SHOWS A "NEGATIVE" CORRELATION BETWEEN 
## PC1 AND AESTHETIC VALUE 
## THIS RELATIONSHIP IS NOT ACTUALLY NEGATIVE
## (ONLY DIRECTIONAL FOR ORDINATION AXES)
## ROTATING THE PCA WILL MAKE THIS RELATIONSHIP POSITIVE 
## WITHOUT CHANGING ANYTHING, AND WILL MAKE 
## OTHER ANALYSES MUCH EASIER
############################################################

# MULTIPLY PC1 BY -1 - IS THERE A BETTER WAY TO DO THIS?
RLStropical_PCA$x[,1] <- RLStropical_PCA$x[,1] * -1
RLStropical_PCA$rotation[,1] <- RLStropical_PCA$rotation[,1] * -1
jpeg("figures_tables/benthic_PCA.jpeg")
fviz_pca_var(RLStropical_PCA, axes = c(1,2), col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE) 
dev.off()

#############################################################################################
#' PERFORM IMPUTATION WITH MISS FOREST PACKAGE
#' IMPUTATION WITH MISS FOREST VARIES FROM ONE ITERATION TO THE NEXT
#' THEREFORE, PEFORM N ITERATIONS AND TAKE THE AVERAGE VALUES
#' HERE, 10 ITERATIONS ARE PEFORMED
#' TO IMPROVE PEFORMANCE OF THE MISS FOREST, LATITUDE AND LONGTIDUE ARE INCLUDED AS COVARIATES
#' SENSTIVITY TESTS SHOW THAT INCLUDING ADDITIONAL COVARIATES DOES NOT IMPROVE PERFORMANCE
#############################################################################################

###########################################
## FIRST - IMPUTE THE MISSING PCA SCORES ##
###########################################

RLS_PC_scores <- data.frame(RLStropical_PCA$x[,1:3]) # AXES 1-3
RLS_PC_scores$SurveyID <- rownames(RLS_PC_scores)
RLS_PC_scores <- merge(selected_survey_info, RLS_PC_scores, all=T)
#write.table(RLS_PC_scores, "RLS_PCA_with_NA.txt")

# WHAT AMOUNT OF MISSING DATA?
sapply(RLS_PC_scores, function(x) paste(round(sum(is.na(x))/length(x),2)*100,"%",sep=""))
#36% FOR PC SCORES

PC1_imputation <- NULL
for(i in 1:10){
  PC1_data  <- RLS_PC_scores %>%
    select("SiteLongitude" , "SiteLatitude", "PC1")
  miss_PC1 <- missForest(PC1_data)$ximp
  PC1_imputation <- cbind(PC1_imputation, miss_PC1$PC1)
}
# SAVE ALL 10 RUNS FOR MULTIPLE IMPUTATION ANALYSIS
PC1_multiple <- data.frame(RLS_PC_scores[,1:4], PC1_imputation)
saveRDS(PC1_multiple, "outputs/PC1_multiple.rds")
PC1_imputation <- rowMeans(PC1_imputation)

PC2_imputation <- NULL
for(i in 1:10){
  PC2_data  <- RLS_PC_scores %>%
    select("SiteLongitude" , "SiteLatitude", "PC2")
  miss_PC2 <- missForest(PC2_data)$ximp
  PC2_imputation <- cbind(PC2_imputation, miss_PC2$PC2)
}
# SAVE ALL 10 RUNS FOR MULTIPLE IMPUTATION ANALYSIS
PC2_multiple <- data.frame(RLS_PC_scores[,1:4], PC2_imputation)
saveRDS(PC2_multiple, "outputs/PC2_multiple.rds")
PC2_imputation <- rowMeans(PC2_imputation)

PC3_imputation <- NULL
for(i in 1:10){
  PC3_data  <- RLS_PC_scores %>%
    select("SiteLongitude" , "SiteLatitude", "PC3")
  miss_PC3 <- missForest(PC3_data)$ximp
  PC3_imputation <- cbind(PC3_imputation, miss_PC3$PC3)
}
# SAVE ALL 10 RUNS FOR MULTIPLE IMPUTATION ANALYSIS
PC3_multiple <- data.frame(RLS_PC_scores[,1:4], PC3_imputation)
saveRDS(PC3_multiple, "outputs/PC3_multiple.rds")
PC3_imputation <- rowMeans(PC3_imputation)

RLS_PC_imputed <- data.frame(selected_survey_info, PC1_imputation, PC2_imputation, PC3_imputation)


# TEST ACCURACY WITH 35% OF DATA DELETED ##
test_data <- RLS_PC_scores %>%
  select(SurveyID, SiteLongitude, SiteLatitude, PC3) %>%
  na.omit
rownames(test_data) <- test_data$SurveyID
test_data$SurveyID <- NULL 

test_data_prime <- test_data
deleted_rows <- sample(1:nrow(test_data_prime), round(nrow(test_data_prime) * 0.35))
test_data_prime[deleted_rows,"PC3"] <- NA #DELETE 35% OF THE DATA AT RANDOM

test_data_prime <- missForest(test_data_prime)$ximp # RUN MISS FOREST

test_data <- test_data[deleted_rows,]
test_data_prime <- test_data_prime[deleted_rows,]

# CHECK CORRELATIONS BETWEEN TRUE AND IMPUTED DATA VAUES
cor.test(test_data$PC3, test_data_prime$PC3)

############################################################
## NEXT PERFORM IMPUTATION DIRECTLY ON BENTHIC CATEGORIES ##
############################################################

RLStropical$SurveyID <- rownames(RLStropical)

RLStropical_with_NA <- merge(selected_survey_info, RLStropical, by="SurveyID", all=T)
#write.table(RLStropical_with_NA, "RLS_benthic_with_NA.txt")

RLStropical_with_NA$SiteCode <- NULL
rownames(RLStropical_with_NA) <- RLStropical_with_NA$SurveyID
RLStropical_with_NA$SurveyID <- NULL

###############################################################
## FIRST, TRY IMPUTING ALL BENTHIC CATEGORIES SIMULTAENOUSLY ##
###############################################################

# TEST ACCURACY WITH 35% OF DATA DELETED ##

test_data <- RLStropical
test_data <- merge(RLStropical, selected_survey_info, by="SurveyID")
rownames(test_data) <- test_data$SurveyID
test_data$SurveyID <- NULL; test_data$SiteCode <- NULL 

test_data_prime <- test_data
deleted_rows <- sample(1:nrow(test_data_prime), round(nrow(test_data_prime) * 0.35))
test_data_prime[deleted_rows,1:(ncol(test_data)-2)] <- NA #DELETE 35% OF THE DATA AT RANDOM

test_data_prime <- missForest(test_data_prime)$ximp # RUN MISS FOREST

test_data <- test_data[deleted_rows,]
test_data_prime <- test_data_prime[deleted_rows,]

# CHECK CORRELATIONS BETWEEN TRUE AND IMPUTED DATA VAUES
for(i in 1:9){
  cor_val <- cor.test(test_data[,i], test_data_prime[,i])$estimate
  print(paste(colnames(test_data)[i], round(cor_val,2)))
}

## PERFORMANCE IS VERY POOR: 0.27 - 0.59 PEARSON CORRELATION

#RLStropical_imputed <- missForest(RLStropical_with_NA)
#RLStropical_imputed <- RLStropical_imputed$ximp

#########################################################
## NEXT, TRY IMPUTING EACH BENTHIC CATEGORY ONE BY ONE ##
#########################################################

benthic_categories <- RLStropical %>%
  select(-c(SurveyID)) %>%
  colnames()

cor_table <- NULL
for(i in unique(benthic_categories)){
  
  for(j in 1:10){
    
    test_data <- RLStropical %>%
      select(i , "SurveyID") # CHOOSE VARIABLE TO TEST
    test_data <- merge(test_data, selected_survey_info, by="SurveyID")
    rownames(test_data) <- test_data$SurveyID
    test_data$SurveyID <- NULL; test_data$SiteCode <- NULL 
    test_data_prime <- test_data
    deleted_rows <- sample(1:nrow(test_data_prime), round(nrow(test_data_prime) * 0.35))
    test_data_prime[deleted_rows,1:(ncol(test_data)-2)] <- NA #DELETE 35% OF THE DATA AT RANDOM
    
    test_data_prime <- missForest(test_data_prime)$ximp # RUN MISS FOREST
    
    test_data <- test_data[deleted_rows,]
    test_data_prime <- test_data_prime[deleted_rows,]
    
    cor_val <- cor.test(test_data[,1], test_data_prime[,1])$estimate
    cor_table <- rbind(cor_table, data.frame(variable=i,correlation=cor_val))
    
    remove(test_data, test_data_prime)
    
  }
  
}


# MEAN CORRELATION ACROSS 10 ITERATIONS
cor_table %>%
  group_by(variable) %>%
  summarise_all(.funs=mean) %>%
  arrange(correlation)

# MAX CORRELATION ACROSS 10 ITERATIONS
cor_table %>%
  group_by(variable) %>%
  summarise_all(.funs=max) %>%
  arrange(correlation)

# MIN CORRELATION ACROSS 10 ITERATIONS
cor_table %>%
  group_by(variable) %>%
  summarise_all(.funs=min) %>%
  arrange(correlation)

ggplot(cor_table) + geom_density(aes(x = correlation, fill = variable), alpha = 0.2)

ggplot(cor_table) + geom_density(aes(x = correlation)) + facet_wrap(~ variable)

# ****PERFORMANCE IS MUCH BETTER WHEN CATEGORIES ARE IMPUTED INDIVIDUALLY***

#' BECAUSE ACCURACY VARIES, RUN IMPUTATION N TIMES AND TAKE AVERAGE OF IMPUTATIONS
#' HERE 10 ITERATIONS ARE PEFORMED

coral_imputation <- NULL
for(i in 1:10){
  coral_data  <- RLStropical_with_NA %>%
    select("SiteLongitude" , "SiteLatitude", "coral")
  miss_coral <- missForest(coral_data)$ximp
  coral_imputation <- cbind(coral_imputation, miss_coral$coral)
}
coral_imputation <- rowMeans(coral_imputation)
min(coral_imputation)
coral_imputation[coral_imputation<0] <- 0 # SOME IMPUTATIONS CAN HAVE SLIGHTLY NEGATIVE NUMBERS - CONVERT TO 0

algae_imputation <- NULL
for(i in 1:10){
  algae_data  <- RLStropical_with_NA %>%
    select("SiteLongitude" , "SiteLatitude", "algae")
  miss_algae <- missForest(algae_data)$ximp
  algae_imputation <- cbind(algae_imputation, miss_algae$algae)
}
algae_imputation <- rowMeans(algae_imputation)
min(algae_imputation)
algae_imputation[algae_imputation<0] <- 0 # SOME IMPUTATIONS CAN HAVE SLIGHTLY NEGATIVE NUMBERS - CONVERT TO 0

rubble_imputation <- NULL
for(i in 1:10){
  rubble_data  <- RLStropical_with_NA %>%
    select("SiteLongitude" , "SiteLatitude", "coral rubble")
  miss_rubble <- missForest(rubble_data)$ximp
  rubble_imputation <- cbind(rubble_imputation, miss_rubble$`coral rubble`)
}
rubble_imputation <- rowMeans(rubble_imputation)
min(rubble_imputation)
rubble_imputation[rubble_imputation<0] <- 0 # SOME IMPUTATIONS CAN HAVE SLIGHTLY NEGATIVE NUMBERS - CONVERT TO 0

CCA_imputation <- NULL
for(i in 1:10){
  CCA_data  <- RLStropical_with_NA %>%
    select("SiteLongitude" , "SiteLatitude", "coralline algae")
  miss_CCA <- missForest(CCA_data)$ximp
  CCA_imputation <- cbind(CCA_imputation, miss_CCA$`coralline algae`)
}
CCA_imputation <- rowMeans(CCA_imputation)
min(CCA_imputation)
CCA_imputation[CCA_imputation<0] <- 0 # SOME IMPUTATIONS CAN HAVE SLIGHTLY NEGATIVE NUMBERS - CONVERT TO 0

micro_mats_imputation <- NULL
for(i in 1:10){
  micro_mats_data  <- RLStropical_with_NA %>%
    select("SiteLongitude" , "SiteLatitude", "microalgal_mats")
  miss_micro_mats <- missForest(micro_mats_data)$ximp
  micro_mats_imputation <- cbind(micro_mats_imputation, miss_micro_mats$microalgal_mats)
}
micro_mats_imputation <- rowMeans(micro_mats_imputation)
min(micro_mats_imputation)
micro_mats_imputation[micro_mats_imputation<0] <- 0 # SOME IMPUTATIONS CAN HAVE SLIGHTLY NEGATIVE NUMBERS - CONVERT TO 0

other_inverts_imputation <- NULL
for(i in 1:10){
  other_inverts_data  <- RLStropical_with_NA %>%
    select("SiteLongitude" , "SiteLatitude", "other sessile invert")
  miss_other_inverts <- missForest(other_inverts_data)$ximp
  other_inverts_imputation <- cbind(other_inverts_imputation, miss_other_inverts$`other sessile invert`)
}
other_inverts_imputation <- rowMeans(other_inverts_imputation)
min(other_inverts_imputation)
other_inverts_imputation[other_inverts_imputation<0] <- 0 # SOME IMPUTATIONS CAN HAVE SLIGHTLY NEGATIVE NUMBERS - CONVERT TO 0

rock_imputation <- NULL
for(i in 1:10){
  rock_data  <- RLStropical_with_NA %>%
    select("SiteLongitude" , "SiteLatitude", "Rock")
  miss_rock <- missForest(rock_data)$ximp
  rock_imputation <- cbind(rock_imputation, miss_rock$Rock)
}
rock_imputation <- rowMeans(rock_imputation)
min(rock_imputation)
rock_imputation[rock_imputation<0] <- 0 # SOME IMPUTATIONS CAN HAVE SLIGHTLY NEGATIVE NUMBERS - CONVERT TO 0

sand_imputation <- NULL
for(i in 1:10){
  sand_data  <- RLStropical_with_NA %>%
    select("SiteLongitude" , "SiteLatitude", "Sand")
  miss_sand <- missForest(sand_data)$ximp
  sand_imputation <- cbind(sand_imputation, miss_sand$Sand)
}
sand_imputation <- rowMeans(sand_imputation)
min(sand_imputation)
sand_imputation[sand_imputation<0] <- 0 # SOME IMPUTATIONS CAN HAVE SLIGHTLY NEGATIVE NUMBERS - CONVERT TO 0

seagrass_imputation <- NULL
for(i in 1:10){
  seagrass_data  <- RLStropical_with_NA %>%
    select("SiteLongitude" , "SiteLatitude", "seagrass")
  miss_seagrass <- missForest(seagrass_data)$ximp
  seagrass_imputation <- cbind(seagrass_imputation, miss_seagrass$seagrass)
}
seagrass_imputation <- rowMeans(seagrass_imputation)
min(seagrass_imputation)
seagrass_imputation[seagrass_imputation<0] <- 0 # SOME IMPUTATIONS CAN HAVE SLIGHTLY NEGATIVE NUMBERS - CONVERT TO 0

## AGGREGATE ALL BENTHIC CATEGORIES TOGETHER ##
imputed_benthic_data <- data.frame(selected_survey_info, coral_imputation, algae_imputation, CCA_imputation,
                                   rubble_imputation, micro_mats_imputation, rock_imputation, sand_imputation, 
                                  seagrass_imputation,  other_inverts_imputation)

#' AS A TEST TRY RUNNING A PCA ON THE IMPUTED BENTHIC CATEGORIES
#' THEN COMPARE THE CORRELATION TO THE IMPUTATION DONE DIRECTLY ON THE PC SCORES

imputed_benthic_transformed <-  imputed_benthic_data %>%
  select(coral_imputation:other_inverts_imputation)
range(imputed_benthic_transformed)
imputed_benthic_transformed <- asin(sqrt(imputed_benthic_transformed)) # ARC SIN TRANSFORMATION PRIOR TO PCA

RLS_PCA_using_imputed <- prcomp(imputed_benthic_transformed, scale. = FALSE)
# ROTATE FOR POSITIVE CORAL VALUES
RLS_PCA_using_imputed$rotation[,1] <- RLS_PCA_using_imputed$rotation[,1] * -1
fviz_pca_var(RLS_PCA_using_imputed, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)

cor.test(RLS_PC_imputed$PC1_imputation, RLS_PCA_using_imputed$x[,1]) # ALMOST IDENTICAL
cor.test(RLS_PC_imputed$PC2_imputation, RLS_PCA_using_imputed$x[,2])
cor.test(RLS_PC_imputed$PC3_imputation, RLS_PCA_using_imputed$x[,3])

## SAVE THE OUTPUTS ##

write.table(imputed_benthic_data, "RLS_benthic_data_imputed.txt")
write.table(RLS_PC_imputed, "RLS_benthic_PCA_imputed.txt")

