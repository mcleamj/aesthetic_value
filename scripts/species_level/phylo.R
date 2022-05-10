###################################################################################################
#' Compute the mean phylogenetic age of the communities
#'
#' @author Matthew McLean, \email {mcleamj@@gmail.com},
#'         Juliette Langlois, \email{juliette.a.langlois@@gmail.com},
#'         Nicolas Mouquet , \email{nicolas.mouquet@@cnrs.fr},
#'         
#' @date 2022/04/20
##################################################################################################

rm(list = ls())

# Load ----

require(dplyr) # for the pipe
source(here::here("R", "functions_age.R"))  

species_table <- read.csv(here::here("outputs", "esthe_base.csv"))

# ----

# ----

# Generate Phylogenetic tree ----

phylo_table          <- species_table[, c("sp_name", "family", "order")]
phylo_table$family   <- factor(phylo_table$family)
phylo_table$tip_name <- paste(phylo_table$family, phylo_table$sp_name, sep = "_")
phylo_table$sp_name <- gsub(" ","_",phylo_table$sp_name)

# Download the phylogenetic tree for all fishes to compute the age
set100_all  <- fishtree::fishtree_complete_phylogeny(species = phylo_table$sp_name)
# 34 species not found

# Dropping names not in list
set100 <- parallel::mclapply(set100_all,function(x){
  ape::drop.tip(phy = x, tip = x$tip.label[!is.element(x$tip.label, phylo_table$sp_name)])
}, mc.cores = 7)

# Save
save(set100, file =  here::here("outputs", "set100.RData"))

rm(set100_all)

# ----

# Species ages ----

load(here::here("outputs", "set100.RData"))

temp        <- set100[[1]]$tip.label
missing_all <- setdiff(phylo_table$sp_name, temp) # 34 species missing in the tree

# Remove the species not found by fishtree
phylo_table        <- phylo_table[phylo_table$sp_name %in% set100[[1]]$tip.label,]
phylo_table$family <- droplevels(phylo_table$family) 

# Compute the age of all species found by fishtree (for each of the 100 trees)
ages      <- do.call(cbind, lapply(set100, get_ages)) 
ages_mean <- apply(ages, 1, mean, na.rm = T) # mean across the trees
ages_mean <- ages_mean[match(phylo_table$sp_name, names(ages_mean))]
phylo_table$ages_mean <- ages_mean

# Add age info to species table
missing_all         <- data.frame(sp_name = missing_all, ages_mean = NA)
ages_mean           <- data.frame(sp_name = names(ages_mean), ages_mean = ages_mean)
rownames(ages_mean) <- NULL
ages_all            <- rbind(ages_mean, missing_all)
ages_all$sp_name <- gsub("_", " ", ages_all$sp_name)

write.csv(ages_all, here::here("outputs", "phylo_base.csv"), row.names = FALSE)


# ----



# # Mean, median and range of ages for surveys ----
# 
# species_table$sp_name <- gsub("_", " ", species_table$sp_name)
# survey_compo$sp_age   <- vector(mode = "numeric", length = nrow(survey_compo))
# 
# for(i in 1:nrow(species_table)){
#   x <- which(survey_compo$sp_name == species_table$sp_name[i])
#   survey_compo$sp_age[x] <- species_table$ages_mean[i]
# }
# 
# survey_table$age_mean <- vector(mode = "numeric", length = nrow(survey_table))
# survey_table$age_med  <- vector(mode = "numeric", length = nrow(survey_table))
# survey_table$age_sd   <- vector(mode = "numeric", length = nrow(survey_table))
# survey_table$age_min  <- vector(mode = "numeric", length = nrow(survey_table))
# survey_table$age_max  <- vector(mode = "numeric", length = nrow(survey_table))
# 
# for(i in 1:nrow(survey_table)){
#   sid                     <- survey_table$SurveyID[i]
#   
#   survey_table$age_mean[i] <- mean(survey_compo$sp_age[which(survey_compo$SurveyID == sid)],
#                                   na.rm = TRUE)
#   survey_table$age_med[i]  <- median(survey_compo$sp_age[which(survey_compo$SurveyID == sid)],
#                                     na.rm = TRUE)
#   survey_table$age_sd [i]  <- sd(survey_compo$sp_age[which(survey_compo$SurveyID == sid)],
#                                 na.rm = TRUE)
#   
#   if(is.na(survey_table$age_sd [i])){survey_table$age_sd [i] <- 0}
#   
#   survey_table$age_min[i]  <- min(survey_compo$sp_age[which(survey_compo$SurveyID == sid)],
#                                  na.rm = TRUE)
#   survey_table$age_max[i]  <- max(survey_compo$sp_age[which(survey_compo$SurveyID == sid)],
#                                  na.rm = TRUE)
# }
# 
# length(which(is.na(survey_table$age_max)))
# 
# survey_compo$SurveyID <- as.factor(survey_compo$SurveyID)
# 
# write.csv(survey_table, here::here("results", "02_survey_table.csv"), row.names = FALSE)
# 
# # ----


