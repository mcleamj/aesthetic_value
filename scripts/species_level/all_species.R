###################################################################################################
#'  Set up the main datafile with esthe scores and classification from taxsize
#'
#' @author Matthew McLean, \email {mcleamj@@gmail.com},
#'         Juliette Langlois, \email{juliette.a.langlois@@gmail.com},
#'         Nicolas Mouquet , \email{nicolas.mouquet@@cnrs.fr},
#'
#'  Produce the file all_species.csv with all species aesthe scores and classification   
#'  Produce FIG_1.png
#'         
#' @date 2022/04/20
##################################################################################################


#Get species list with esthe scores 

  aesthe_table <- read.csv(here::here("data", "aesthe_table.csv"))
  aesthe_table$sp_name <- as.character(gsub("_"," ",aesthe_table$sp_name))
  
  aesthe_table <- aesthe_table[!duplicated(aesthe_table$sp_name), ] #remove duplicates if any
  
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
  
  classif       <- get_classif(aesthe_table$sp_name)
  colnames(classif) <- c("sp_name", colnames(classif[-1]))
  classif <- classif[,-dim(classif)[2]]
  
  ## Scaridae is a subfamily of the Labridae (F. Leprieur personal comment)
    classif$family[classif$family == "Scaridae"] <- "Labridae"


# Compute the aesthe contribution of each species
  # with parameters from Tribot, A.S, Deter, J., Claverie, T., Guillhaumon, F., Villeger, S., & Mouquet, N. (2019). Species diversity and composition drive the aesthetic value of coral reef fish assemblages. Biology letters, 15, 20190703, doi:10.1098/rsbl.2019.0703
  # positive and negative effect are relative to the espected effect computed with the species 
    
  ## Computing the aesthe_effect
    aesthe_table$aesthe_effect <- (log(aesthe_table$esthe_score) - 7.3468679)/7.937672
  
#merge & save 
  
  all_species <- merge(aesthe_table, classif, by = "sp_name")  
  write.csv(x = all_species, file = here::here("outputs", "all_species.csv"),
            row.names = FALSE)
  
#FIG 1
  
  library(ggplot2)
  
  all_species$rank <- rank(all_species$aesthe_effect)
  
  nneg <- sum(all_species$aesthe_effect < 0) # 1805 species have a negative effect
  npos <- sum(all_species$aesthe_effect > 0) # 656 species have a positive effect
  
  Fig_1 <- ggplot(all_species, ggplot2::aes(y=aesthe_effect,x = rank)) +
    geom_point(col="royalblue1",alpha=0.5)+
    theme_bw()+
    xlab("Ranks") + ylab("Net species aesthetic effects")+
    ylim(-0.04,0.04)+
    geom_hline(yintercept=0, linetype="dashed", 
               color = "gray", size=1)+
    geom_vline(xintercept=nneg, linetype="dashed", 
               color = "gray", size=1)+
    geom_text(x=100, y=-0.003, label=paste("n=",nneg),size=4,col="gray")+
    geom_text(x=2400, y=+0.003, label=paste("n=",npos),size=4,col="gray")
  
  ggsave(here::here("figures_tables","fig_1.png"), plot = Fig_1,
         width = 10, height = 5, dpi = 300, units = "in", device='png')
  









