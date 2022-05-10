###################################################################################################
#'  Combine all data at the species level
#'
#' @author Matthew McLean, \email {mcleamj@@gmail.com},
#'         Juliette Langlois, \email{juliette.a.langlois@@gmail.com},
#'         Nicolas Mouquet , \email{nicolas.mouquet@@cnrs.fr},

#' @date 2022/04/20
##################################################################################################

rm(list = ls())


# Load ----
  init_table  <- read.csv(here::here("outputs", "esthe_base.csv"))
  phylo_table  <- read.csv(here::here("outputs", "phylo_base.csv"))
  fishchips_table  <- read.csv(here::here("outputs", "fishchips_base.csv"))
  funk_table <- read.csv(here::here("outputs", "funk_base.csv"))
  DiUi_table <- read.csv(here::here("outputs", "DiUi_base.csv"))
  
# ----
  
# Merge and save ----
  
  all_species <- merge(init_table, phylo_table, by = "sp_name",all = TRUE) 
  all_species <- merge(all_species, fishchips_table, by = "sp_name",all = TRUE) 
  all_species <- merge(all_species, funk_table, by = "sp_name",all = TRUE) 
  all_species <- merge(all_species, DiUi_table, by = "sp_name",all = TRUE) 
  
  write.csv(all_species, here::here("outputs", "all_species.csv"), row.names = FALSE)
 
# ---- 
  
  all_species <- read.csv(here::here("results", "01_all_species.csv"))
  
  plot(all_species$esthe_score,all_species$Di)
  plot(all_species$esthe_score,log(all_species$Ui))
  plot(all_species$esthe_score,all_species$Fish_chips)
  
  plot(all_species$Fish_chips,all_species$Di)
  
  library(ggplot2)
  library(ggtern)
  ggtern(data=all_species, aes(x=Fish_chips,y=esthe_score, z=Di)) + geom_point()
  
  
  normalize <- function(x){
    (x-min(x, na.rm = T))/(max(x, na.rm = T) - min(x, na.rm = T))
  }
  
  test <- data.frame(esthe_score=normalize(all_species$esthe_score),Di=normalize(all_species$Di),Fish_Chips=normalize(all_species$Fish_chips))
  test <- test[complete.cases(test),]    
    
  library('Ternary')
  par(mar = rep(0.2, 4))
  TernaryPlot(axis.labels = seq(0, 1, by = 0.1),alab = 'esthe_score', blab = 'Di', clab = 'Fish_Chips')
  
  
  coordinates <- cbind(test$esthe_score,
                       test$Di,
                       test$Fish_Chips)
  
  #ColourTernary(TernaryDensity(coordinates, resolution = 10L))
  TernaryPoints(coordinates, col = 'red', pch = 1)
  TernaryDensityContour(coordinates, resolution = 30L)
  
  
  
  
  
  
  