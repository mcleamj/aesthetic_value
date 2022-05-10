###################################################################################################
#' Functional diversity
#'
#'
#'
#'@author Juliette Langlois, \email{juliette.a.langlois@@gmail.com},
#'        Nicolas Mouquet, \email{nicolas.mouquet@@cnrs.fr},

#' @date 2021/06/29  first created
##################################################################################################

rm(list = ls())

# Load data ----

require(dplyr) # for the pipe
species_table <- read.csv(here::here("outputs", "esthe_base.csv"))
funct_table   <- read.csv(here::here("data", "fonctio_table.csv"), sep = ";")

funct_table <- funct_table[funct_table$sp_name %in% species_table$sp_name,]

# ----

# Select the traits  ----
funct_table$habitat       <- tolower(funct_table$habitat)
funct_table$diel_activity <- tolower(funct_table$diel_activity)
funct_table$water_column  <- tolower(funct_table$water_column)
funct_table$trophic_group <- tolower(funct_table$trophic_group)

rownames(funct_table) <- funct_table$sp_name
funct_table         <- funct_table[, c('max_length', 'trophic_level','thermal_mp_5min_95max',
                                       'thermal_95thmax', 'trophic_group','water_column',
                                       'diel_activity','habitat')]

# Remove the species for which we have none of the traits
funct_table <- funct_table[rowSums(is.na(funct_table)) != ncol(funct_table), ]

# ---- 

# Predict missing value for traits ----

  # Which percentage of NAs?
  nbna <- length(which(is.na(funct_table)))/(nrow(funct_table)*ncol(funct_table))*100
  if(nbna <=5){
    cat(round(nbna, digits = 2),"% of NAs \n", ">>> You can use missForest!\n")
  } # eo if
  
  # Lines whith more than 50% NAs?
  temp <- funct_table[rowSums(is.na(funct_table)) <= ncol(funct_table)*0.5, ]
  ifelse(nrow(temp) == nrow(funct_table),
         cat("\n", "No line with more than 50% NAs\n"),
         cat("\n", "Lines to remove\n"))
  
  # Test efficiency of missForest
  # create artificial NAs and see how good it is to replace them
  # No need to run every time, once is ok
  
  # subset of no na
  nona <- funct_table %>% na.omit()
  
  # create artificial NAs
  artna           <- nona
  
  # the number of random values to replace
  artna_mis               <- missForest::prodNA(x = artna, noNA = 0.05)
  artna_mis$trophic_group <- as.factor(artna_mis$trophic_group)
  artna_mis$water_column  <- as.factor(artna_mis$water_column)
  artna_mis$diel_activity <- as.factor(artna_mis$diel_activity)
  artna_mis$habitat       <- as.factor(artna_mis$habitat)
  artna_recons <- missForest::missForest(xmis = artna_mis, verbose = TRUE)
  
  # Check if the predictions are good with correlation
  perf <- diag(cor(x = artna[, c("max_length", "trophic_level", "thermal_mp_5min_95max",
                                 "thermal_95thmax")],
                   y = artna_recons$ximp[, c("max_length", "trophic_level",
                                             "thermal_mp_5min_95max", "thermal_95thmax")]))
  trophic_group <- (length(which(as.factor(artna$trophic_group) ==
                                   artna_recons$ximp$trophic_group)))/nrow(artna)
  water_column  <- (length(which(as.factor(artna$water_column) ==
                                   artna_recons$ximp$water_column)))/nrow(artna)
  diel_activity <- (length(which(as.factor(artna$diel_activity) ==
                                   artna_recons$ximp$diel_activity)))/nrow(artna)
  habitat       <- (length(which(as.factor(artna$habitat) ==
                                   artna_recons$ximp$habitat)))/nrow(artna)
  perf <- c(perf, trophic_group, water_column, diel_activity, habitat)
  # all between 0.98 and 0.999 ==> ok method accepted
  
  rm(artna, nona, artna_recons, temp, nbna, artna_mis)
  
  if(length(which(is.na(funct_table))) !=0){
    funct_table               <- as.data.frame(funct_table)
    funct_table$trophic_group <- as.factor(funct_table$trophic_group)
    funct_table$water_column  <- as.factor(funct_table$water_column)
    funct_table$diel_activity <- as.factor(funct_table$diel_activity)
    funct_table$habitat       <- as.factor(funct_table$habitat)
    funct_table               <-  missForest::missForest(funct_table, verbose = TRUE)
    funct_table               <- funct_table$ximp
  } # eo if
  
  length(which(is.na(funct_table)))
  funct_table$sp_name <- rownames(funct_table)
  
  write.csv(funct_table, 
            file = here::here("outputs", "funk_base.csv"),
            row.names = FALSE)
  
# ----

# Log and normalize numeric variables ----

funct_table <- read.csv(here::here("outputs", "funk_base.csv"))


# Correlogram
GGally::ggpairs(funct_table[,c('max_length', 'trophic_level','thermal_mp_5min_95max',
                               'thermal_95thmax')], title = "correlogram with ggpairs()") 

# Log if skewed
log_var <- c('max_length')
for (id in log_var) funct_table[,id] <- log(funct_table[,id])

# Normalize all
normalize <- function(x){
  (x-min(x, na.rm = T))/(max(x, na.rm = T) - min(x, na.rm = T))
}

funct_table$max_length            <- normalize(funct_table$max_length)
funct_table$thermal_mp_5min_95max <- normalize(funct_table$thermal_mp_5min_95max)
funct_table$thermal_95thmax      <- normalize(funct_table$thermal_95thmax)
funct_table$trophic_level        <- normalize(funct_table$trophic_level)

rm(id, log_var)

# ----

# Compute the distance matrix ----
# Group the variables according to their type
# Quantitative variables
quant_var           <- cbind.data.frame(funct_table$max_length, funct_table$thermal_mp_5min_95max,
                                        funct_table$thermal_95thmax, funct_table$trophic_level)
rownames(quant_var) <- funct_table$sp_name
# Nominal variables
nom_var            <- cbind.data.frame(funct_table$trophic_group, funct_table$water_column,
                                       funct_table$diel_activity, funct_table$habitat)
rownames(nom_var)  <- funct_table$sp_name

# Compute distance matrix
disTraits <- ade4::dist.ktab(ade4::ktab.list.df(list(quant_var, nom_var)), c("Q","N"),
                             scan = FALSE) %>% as.matrix()
# Save
save(disTraits, file = here::here("outputs","disTraits.RData"))

# ----    

# Compute the functional distinctiveness and uniqueness ----

load(here::here("outputs", "disTraits.RData"))

# Build hypothetical community where all species are presents. 
Sim_commu           <- matrix(1, 1, ncol(disTraits))
colnames(Sim_commu) <- colnames(disTraits)

# Compute Ui & Di for each species in the hypothetical community
Di           <- t(funrar::distinctiveness(Sim_commu, disTraits))
colnames(Di) <- "Di"

Ui           <- funrar::uniqueness(Sim_commu, disTraits)
rownames(Ui) <- Ui$species

FR           <- merge(Di, Ui, by = "row.names", all.x = TRUE)

FR <- FR[,-1]

colnames(FR) <- c("Di","sp_name","Ui")

write.csv(FR, 
          file = here::here("outputs", "DiUi_base.csv"),
          row.names = FALSE)
# ----    
