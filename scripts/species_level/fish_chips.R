###################################################################################################
#' Scrap FISHBASE to collect intel on Fishery/Subsistence Importance (Fish & Chips)
#'
#' @author Nicolas Mouquet , \email{nicolas.mouquet@@cnrs.fr},
#'         Matthew McLean, \email {mcleamj@@gmail.com},
#'         
#' @date 2022/04/06
##################################################################################################

rm(list = ls())

#function norm01 

norm01 <- function(dat){
  
  min_dat <- min(dat,na.rm=TRUE)
  max_dat <- max(dat,na.rm=TRUE)
  dat <- dat+abs(min_dat)
  min_dat <- min(dat,na.rm=TRUE)
  max_dat <- max(dat,na.rm=TRUE)
  
  if (min_dat!=max_dat) dat <- (dat-min_dat)/(max_dat-min_dat) else { 
    if (min_dat==0) dat <- dat else dat[!is.na(dat)]=1}

  return(dat)
}



#get the species list
  esthe_table  <- read.csv(here::here("outputs", "esthe_base.csv"))

#import the intel (run only once)
  ct <- 0
  sp_list <- esthe_table$sp_name
  data_fishbase  <- do.call(rbind, lapply(sp_list, function(id){
    
    ct <<- ct+1
    ER <- FALSE
    dataFB <- tryCatch({
      rfishbase::species(id,fields=c('Importance','PriceCateg','UsedforAquaculture','UsedasBait','GameFish'))},
      error = function(e){ER <- TRUE},
      warning = function(w){ER <- TRUE})
    
    dataFB <- as.data.frame(dataFB)
    dataFB$sp_names_fishbase <- id
    if (ER==FALSE) dataFB
    cat("n=",ct," / ",length(sp_list), "\n")
    dataFB
    
  }))
  write.csv(x = data_fishbase, file = here::here("outputs/fishbase_base.csv"), row.names = FALSE)

#load the data and create the Fish_chips table 

  data_fishbase  <- read.csv(here::here("outputs", "fishbase_base.csv"))
  dat_fishchips <-data.frame(acc_sci_name=data_fishbase$sp_names_fishbase)

#fishery Importance 

  dat_fishchips$Importance[data_fishbase$Importance %in% "highly commercial"] <- 4
  dat_fishchips$Importance[data_fishbase$Importance %in% "commercial"] <- 3
  dat_fishchips$Importance[data_fishbase$Importance %in% "minor commercial"] <- 2
  dat_fishchips$Importance[data_fishbase$Importance %in% "subsistence fisheries"] <- 1
  dat_fishchips$Importance[data_fishbase$Importance %in% "of potential interest"] <- 1
  dat_fishchips$Importance[data_fishbase$Importance %in% "of no interest"] <- 0
  dat_fishchips$Importance[data_fishbase$Importance %in% c("unknown") ] <- NA
  dat_fishchips$Importance[data_fishbase$Importance %in% c("NA", NA) ] <- NA
  
  Importance <- norm01(as.numeric(dat_fishchips$Importance))

#PriceCateg

  dat_fishchips$PriceCateg[data_fishbase$PriceCateg %in% "very high"] <- 3
  dat_fishchips$PriceCateg[data_fishbase$PriceCateg %in% "high"] <- 2
  dat_fishchips$PriceCateg[data_fishbase$PriceCateg %in% "medium"] <- 1
  dat_fishchips$PriceCateg[data_fishbase$PriceCateg %in% "low"] <- 0
  dat_fishchips$PriceCateg[data_fishbase$PriceCateg %in% c("unknown","NA", NA) ] <- NA
  
  PriceCateg <- norm01(as.numeric(dat_fishchips$PriceCateg))

#UsedforAquaculture

  dat_fishchips$UsedforAquaculture[data_fishbase$UsedforAquaculture %in% "commercial"] <- 2
  dat_fishchips$UsedforAquaculture[data_fishbase$UsedforAquaculture %in% c("experimental","likely future use")] <- 1
  dat_fishchips$UsedforAquaculture[data_fishbase$UsedforAquaculture %in% "never/rarely"] <- 0
  dat_fishchips$UsedforAquaculture[data_fishbase$UsedforAquaculture %in% c("NA", NA) ] <- NA
  
  UsedforAquaculture <- norm01(as.numeric(dat_fishchips$UsedforAquaculture))
  

#UsedasBait

  dat_fishchips$UsedasBait[data_fishbase$UsedasBait %in% "usually"] <- 2
  dat_fishchips$UsedasBait[data_fishbase$UsedasBait %in% "occasionally"] <- 1
  dat_fishchips$UsedasBait[data_fishbase$UsedasBait %in% "never/rarely"] <- 0
  dat_fishchips$UsedasBait[data_fishbase$UsedasBait %in% c("NA", NA) ] <- NA

  UsedasBait <- norm01(as.numeric(dat_fishchips$UsedasBait))
  
  
#GameFish
  
  dat_fishchips$GameFish[data_fishbase$GameFish %in% 0] <- 0
  dat_fishchips$GameFish[data_fishbase$GameFish %in% -1] <- 1
  dat_fishchips$GameFish[data_fishbase$GameFish %in% c("NA", NA) ] <- NA
  
  GameFish <- norm01(as.numeric(dat_fishchips$GameFish))


#Combine allinto FishChips 
  
  dat_fishchips$Fish_chips <- norm01(apply(data.frame(Importance,PriceCateg,UsedforAquaculture,UsedasBait,GameFish),1,mean,na.rm=TRUE))
  dat_fishchips$Fish_chips[is.nan(dat_fishchips$Fish_chips)] <- NA
  
  plot(density(dat_fishchips$Fish_chips,na.rm=TRUE), lwd = 2, col = "red",
       main = "Density")
  rug(jitter(dat_fishchips$Fish_chips))
  
  colnames(dat_fishchips) <- c("sp_name", colnames(dat_fishchips[-1]))
  
  write.csv(x = dat_fishchips, file = here::here("outputs", "fishchips_base.csv"),
            row.names = FALSE)

  

