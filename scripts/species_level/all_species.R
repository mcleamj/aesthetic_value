###################################################################################################
#'  Set up the main datafile with esthe scores and classification from taxsize
#'
#' @author Matthew McLean, \email {mcleamj@@gmail.com},
#'         Juliette Langlois, \email{juliette.a.langlois@@gmail.com},
#'         Nicolas Mouquet , \email{nicolas.mouquet@@cnrs.fr},
#'
#'  Produce the file all_species.csv with all species aesthe scores and classification   
#'         
#' @date 2022/04/20
##################################################################################################


#Get species list with esthe scores 

  esthe_table <- read.csv(here::here("data", "esthe_table.csv"))
  esthe_table$sp_name <- as.character(gsub("_"," ",esthe_table$sp_name))
  
  esthe_table <- esthe_table[!duplicated(esthe_table$sp_name), ] #remove duplicates if any
  
#Get the classification from taxsize

  get_classif <- function(names) {
    
    A <- data.frame(name = names, kingdom = NA, phylum = NA, class = NA, 
                    order = NA, family = NA, genus = NA, species = NA)
    
    for (i in 1:length(names)) {
      
      classi <- taxize::classification(names[i], db = "worms")
      classi <- classi[[1]]
      
      if (length(classi) == 1) {
        
        A[i, "kingdom"] <- NA
        A[i, "phylum"]  <- NA
        A[i, "class"]   <- NA
        A[i, "order"]   <- NA
        A[i, "family"]  <- NA
        A[i, "genus"]   <- NA
        A[i, "species"] <- NA
        
      } else {
        
        kingdom    <- classi$name[classi$rank == "Kingdom"]
        phylum     <- classi$name[classi$rank == "Phylum"]
        class      <- classi$name[classi$rank == "Class"]
        order      <- classi$name[classi$rank == "Order"]
        family     <- classi$name[classi$rank == "Family"]
        genus      <- classi$name[classi$rank == "Genus"]
        species    <- classi$name[classi$rank == "Species"]
        
        A[i, "kingdom"] <- kingdom
        A[i, "phylum"]  <- phylum    
        A[i, "class"]   <- class
        A[i, "order"]   <- order     
        A[i, "family"]  <- family    
        A[i, "genus"]   <- genus     
        A[i, "species"] <- species   
      }
    }
    
    A
  }
  
  classif       <- get_classif(esthe_table$sp_name)
  colnames(classif) <- c("sp_name", colnames(classif[-1]))
  classif <- classif[,-dim(classif)[2]]
  
  ## Scaridae is a subfamily of the Labridae (F. Leprieur personal comment)
    classif$family[classif$family == "Scaridae"] <- "Labridae"

#merge and save 
  
  species_table     <- merge(esthe_table, classif, by = "sp_name")  
  
  write.csv(x = species_table, file = here::here("output", "all_species.csv"),
            row.names = FALSE)









