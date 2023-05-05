###################################################################################################
#'  Set up the main datafile with esthe scores and classification from taxsize
#'
#' @author Nicolas Mouquet , \email{nicolas.mouquet@@cnrs.fr}
#' @author Matthew McLean, \email {mcleamj@@gmail.com},
#' @author Juliette Langlois, \email{juliette.a.langlois@@gmail.com}
#'
#'  Produce the file all_species.csv with all species aesthe scores and classification   
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
  
 
  ex <- c("Chaetodontoplus septentrionalis","Calloplesiops altivelis","Paracanthurus hepatus",
          "Chelmonops curiosus","Acanthurus blochii","Scarus rivulatus","Pseudolabrus guentheri",
          "Scolopsis ghanam","Sarpa salpa","Salmo salar","Arripis truttacea","Sardina pilchardus")
  
  colsp <- "#40403F"
  size_cr <- 6
  font_cr <- 3
  
  ex_fl <- c("Chaetodontoplus_septentrionalis_A_1.png","Calloplesiops_altivelis_A_2.png","Paracanthurus_hepatus_A_3.png",
             "Chelmonops_curiosus_A_2.png","Acanthurus_blochii_A_1.png","Scarus_rivulatus_M_2.png",
             "Pseudolabrus_guentheri_A_2.png","Scolopsis_ghanam_A_3.png","Sarpa_salpa_A_2.png",
             "Salmo_salar_A_2.png","Arripis_truttacea_A_1.png","Sardina_pilchardus_A_1.png")

  Fig_1 <- ggplot(all_species, ggplot2::aes(y=aesthe_effect,x = rank)) +
    geom_point(col="royalblue1",alpha=0.5)+
    theme_bw()+
    xlab("Ranks") + ylab("Net species aesthetic effects")+
    theme(axis.text.x = element_text(size = 10),
          axis.title.x = element_text(size = 14),
          axis.text.y = element_text(size = 10),
          axis.title.y = element_text(size = 14))+
    ylim(-0.04,0.04)+
    geom_hline(yintercept=0, linetype="dashed", 
               color = "gray", size=0.8)+
    geom_vline(xintercept=nneg, linetype="dashed", 
               color = "gray", size=0.8)+
    geom_text(x=100, y=-0.003, label=paste("n=",nneg),size=4,col="#7D7D7C")+
    geom_text(x=2400, y=+0.003, label=paste("n=",npos),size=4,col="#7D7D7C")+
    geom_point(aes(x=rank[sp_name%in%ex[1]],
                   y=aesthe_effect[sp_name%in%ex[1]]),
                   shape = 21, colour = colsp, fill = "white",size=size_cr)+
    geom_text(x=all_species$rank[all_species$sp_name%in%ex[1]], 
              y=all_species$aesthe_effect[all_species$sp_name%in%ex[1]], 
              label="1",size=font_cr,col=colsp)+
    geom_point(aes(x=rank[sp_name%in%ex[2]],
                   y=aesthe_effect[sp_name%in%ex[2]]),
               shape = 21, colour = colsp, fill = "white",size=size_cr)+
    geom_text(x=all_species$rank[all_species$sp_name%in%ex[2]], 
              y=all_species$aesthe_effect[all_species$sp_name%in%ex[2]], 
              label="2",size=font_cr,col=colsp)+
    geom_point(aes(x=rank[sp_name%in%ex[3]],
                   y=aesthe_effect[sp_name%in%ex[3]]),
               shape = 21, colour = colsp, fill = "white",size=size_cr)+
    geom_text(x=all_species$rank[all_species$sp_name%in%ex[3]], 
              y=all_species$aesthe_effect[all_species$sp_name%in%ex[3]], 
              label="3",size=font_cr,col=colsp)+
    geom_point(aes(x=rank[sp_name%in%ex[4]],
                   y=aesthe_effect[sp_name%in%ex[4]]),
               shape = 21, colour = colsp, fill = "white",size=size_cr)+
    geom_text(x=all_species$rank[all_species$sp_name%in%ex[4]], 
              y=all_species$aesthe_effect[all_species$sp_name%in%ex[4]], 
              label="4",size=font_cr,col=colsp)+
    geom_point(aes(x=rank[sp_name%in%ex[5]],
                   y=aesthe_effect[sp_name%in%ex[5]]),
               shape = 21, colour = colsp, fill = "white",size=size_cr)+
    geom_text(x=all_species$rank[all_species$sp_name%in%ex[5]], 
              y=all_species$aesthe_effect[all_species$sp_name%in%ex[5]], 
              label="5",size=font_cr,col=colsp)+
    geom_point(aes(x=rank[sp_name%in%ex[6]],
                   y=aesthe_effect[sp_name%in%ex[6]]),
               shape = 21, colour = colsp, fill = "white",size=size_cr)+
    geom_text(x=all_species$rank[all_species$sp_name%in%ex[6]], 
              y=all_species$aesthe_effect[all_species$sp_name%in%ex[6]], 
              label="6",size=font_cr,col=colsp)+
    geom_point(aes(x=rank[sp_name%in%ex[7]],
                   y=aesthe_effect[sp_name%in%ex[7]]),
               shape = 21, colour = colsp, fill = "white",size=size_cr)+
    geom_text(x=all_species$rank[all_species$sp_name%in%ex[7]], 
              y=all_species$aesthe_effect[all_species$sp_name%in%ex[7]], 
              label="7",size=font_cr,col=colsp)+
    geom_point(aes(x=rank[sp_name%in%ex[8]],
                   y=aesthe_effect[sp_name%in%ex[8]]),
               shape = 21, colour = colsp, fill = "white",size=size_cr)+
    geom_text(x=all_species$rank[all_species$sp_name%in%ex[8]], 
              y=all_species$aesthe_effect[all_species$sp_name%in%ex[8]], 
              label="8",size=font_cr,col=colsp)+
    geom_point(aes(x=rank[sp_name%in%ex[9]],
                   y=aesthe_effect[sp_name%in%ex[9]]),
               shape = 21, colour = colsp, fill = "white",size=size_cr)+
    geom_text(x=all_species$rank[all_species$sp_name%in%ex[9]], 
              y=all_species$aesthe_effect[all_species$sp_name%in%ex[9]], 
              label="9",size=font_cr,col=colsp)+
    geom_point(aes(x=rank[sp_name%in%ex[10]],
                   y=aesthe_effect[sp_name%in%ex[10]]),
               shape = 21, colour = colsp, fill = "white",size=size_cr)+
    geom_text(x=all_species$rank[all_species$sp_name%in%ex[10]], 
              y=all_species$aesthe_effect[all_species$sp_name%in%ex[10]], 
              label="10",size=font_cr,col=colsp)+
    geom_point(aes(x=rank[sp_name%in%ex[11]],
                   y=aesthe_effect[sp_name%in%ex[11]]),
               shape = 21, colour = colsp, fill = "white",size=size_cr)+
    geom_text(x=all_species$rank[all_species$sp_name%in%ex[11]], 
              y=all_species$aesthe_effect[all_species$sp_name%in%ex[11]], 
              label="11",size=font_cr,col=colsp)+
    geom_point(aes(x=rank[sp_name%in%ex[12]],
                   y=aesthe_effect[sp_name%in%ex[12]]),
               shape = 21, colour = colsp, fill = "white",size=size_cr)+
    geom_text(x=all_species$rank[all_species$sp_name%in%ex[12]], 
              y=all_species$aesthe_effect[all_species$sp_name%in%ex[12]], 
              label="12",size=font_cr,col=colsp)
  
  
  ggsave(here::here("figures_tables","fig_1.png"), plot = Fig_1,
         width = 10, height = 5, dpi = 300, units = "in", device='png')
  









