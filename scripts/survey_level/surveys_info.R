###################################################################################################
#'  Information at the survey scale
#'
#' @author Matthew McLean, \email {mcleamj@@gmail.com},
#'         Juliette Langlois, \email{juliette.a.langlois@@gmail.com},
#'         Nicolas Mouquet , \email{nicolas.mouquet@@cnrs.fr},
#'         
#'         
#'
#' @date 2021/08/11
##################################################################################################

rm(list = ls())

# See Mattia's script to uniformize

survey_spcompo <- readRDS(here::here("data", "RLS_fishdata.rds"))

# Loading the data sets ----

# the composition of the communities

# survey_spcompo <- survey_spcompo[-grep(pattern = "spp", x = survey_spcompo$TAXONOMIC_NAME),]
# survey_spcompo <- survey_spcompo[-grep(pattern = "\\.", x = survey_spcompo$TAXONOMIC_NAME),]
# survey_spcompo$TAXONOMIC_NAME <- gsub(" \\[green form\\]", "", survey_spcompo$TAXONOMIC_NAME)    

# check names of survey_spcompo against the worms database
#   name_ini <- unique(as.character(survey_spcompo$TAXONOMIC_NAME))
#       
# # create empty data frame for the new names
#   name_fin <- rep(x = NA, length(name_ini))
#   checked  <- data.frame(name_ini = name_ini, name_fin = name_fin)
# 
# # fill the dataframe
#   for (i in 1:length(name_ini)){
#     if (name_ini[i] == "Chrysophrys auratus"){checked$name_fin[i] <- "Pagrus auratus"
#     } else if (name_ini[i] == "Epinephelides guttatus"){checked$name_fin[i] <- "Epinephelus guttatus"
#     }else{
#       id_table <- taxize::get_wormsid(sci_com = name_ini[i], searchtype = "scientific")
#       u        <- worms::wormsbyid(as.numeric(id_table[1]))
#       if(u$status =="accepted"){
#         checked$name_fin[i] <- u$scientificname[which(u$status == "accepted")]
#       }else{
#         checked$name_fin[i] <- u$valid_name[which(u$status != "accepted")]
#       } # eo else
#     } # eo else
#   } # eo for i
#   
#   write.csv(checked, here::here("results", "00_species_names.csv"), row.names = FALSE)

# checked_names     <- read.csv(here::here("results", "00_species_names.csv"))
# survey_spcompo$sp_name <- survey_spcompo$TAXONOMIC_NAME
# 
# for(i in 1:nrow(checked_names)){
#   survey_spcompo$sp_name[
#     which(survey_spcompo$TAXONOMIC_NAME == as.character(checked_names$name_ini[i]))] <- 
#     as.character(checked_names$name_fin[i])
# } # eo for i
# 
# survey_spcompo <- survey_spcompo[which(survey_spcompo$sp_name %in% species_table$sp_name),
#                        c("SurveyID", "sp_name")]
# 
# write.csv(survey_spcompo, here::here("results", "00_survey_composition.csv"), row.names = FALSE)
# 
# rm(checked_names, i, survey_spcompo)
# ----

survey_compo         <- read.csv(here::here("results", "00_survey_composition.csv"))
survey_compo$sp_name <- as.character(survey_compo$sp_name)

# species_table <- species_table[which(species_table$sp_name %in% survey_compo$sp_name),]


# Get the number of species per survey ----
library(dplyr)
dataframe_survey <- survey_compo %>%
  group_by(SurveyID) %>%
  dplyr::summarize(nb_species = length(unique(sp_name)))

colnames(dataframe_survey)[1] <- "survey_id"
write.csv(dataframe_survey, here::here("results", "00_survey_table.csv"))

# ----

# Abundance matrix ----
ab_mat <- as.data.frame(matrix(nrow = length(unique(survey_compo$SurveyID)),
                               ncol = length(unique(survey_compo$sp_name))))
rownames(ab_mat) <- unique(survey_compo$SurveyID)
colnames(ab_mat) <- unique(survey_compo$sp_name)

parallel::mclapply(1:nrow(ab_mat),function(i){
  for(j in 1:ncol(ab_mat)){
    ab_mat[i,j] <<- length(which(survey_compo$SurveyID == rownames(ab_mat)[i] &
                                   survey_compo$sp_name == colnames(ab_mat)[j]))
  }
},mc.cores=7)



for(i in 1:nrow(ab_mat)){
  for(j in 1:ncol(ab_mat)){
    ab_mat[i,j] <- length(which(survey_compo$SurveyID == rownames(ab_mat)[i] &
                                  survey_compo$sp_name == colnames(ab_mat)[j]))
  }
}


# ----

# Relative abundance matrix ----

relab_mat <- funrar::make_relative(abund_matrix = ab_mat)
save(relab_mat, here::here("results", "00_relab.RData"))

# ----

# Presence absence matrix ----

presabs            <- relab_mat
presabs[presabs>0] <- 1

# # Matrix with presence and absence of each species per survey ----
# # List of species
# sp            <- as.character(unique(species_table$sp_name))
# surveyid_vect <- as.character(dataframe_survey$survey_id)
# 
# # Create empty matrix
# data_presab <- matrix(0,
#                       nrow = length(dataframe_survey$survey_id),
#                       ncol = length(sp))
# rownames(data_presab) <- surveyid_vect
# colnames(data_presab) <- sp

# # Fill the matrix
# start.time <- Sys.time()
# # length(surveyid_vect)
# for (i in 1:length(surveyid_vect)) {
#   # list the species in the survey
#   survey    <- surveyid_vect[i]
#   sp_survey <- unique(survey_compo$sp_name[survey_compo$SurveyID == survey])
#   sp_survey <- sp_survey[sp_survey %in% sp]
#   for (j in 1:length(sp_survey)) {
#     # identify the cell of this species in this survey and fill it with 1 if required
#     temp_sp <-  sp_survey[j]
#     k       <- which(colnames(data_presab) == temp_sp)
#     data_presab[i, k] <- 1
#     
#   } # eo for j
#   cat(paste0("survey ", i, " \n"))
# } # eo for i

# end.time   <- Sys.time()
# time.taken <- end.time - start.time
# time.taken

# # Format as data frame
# data_presab           <- as.data.frame(data_presab)

# Save
save(presabs, file = here::here("results", "00_presabs.RData"))

# Load
load(here::here("results", "00_presabs.RData"))


# ----