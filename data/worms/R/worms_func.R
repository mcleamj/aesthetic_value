#' @title GET AphiaRecordsByNames 
#' @importFrom httr  GET content
#' @importFrom plyr rbind.fill
#' 
#' @description takes character vector with taxon names and retrives AphiaRecords from WoRMS 
#'
#' @param taxon_names character vector with names of taxa to look up.
#' @param chunksize there is a limit to the number of taxa that can be looked up at once, so request are split up into chunks. This limit seems to be variable. 50 is very safe.  
#' @param verbose be verbose
#' @param ids add column "id" and "name" with running id and search names
#' @param match taxon_names that could not retrieved will be retried with \code{\link{wormsbymatchnames}}. Implies "id=TRUE"
#' @param like Add a "\%"-sign after the ScientificName (SQL LIKE function). Default=true
#' @param marine_only Limit to marine taxa. Default=true
#' @param sleep_btw_chunks_in_sec pause between requests 
#' 
#' @return a data frame.
#' @details This function will take a character vector with taxon names, 
#' retrive AphiaRecords (CC-BY) from www.marinespecies.org using the 
#' GET /AphiaRecordsByName/{ScientificName} Method described at
#' http://www.marinespecies.org/rest/.
#' Results will be output to a data.frame with each row being a record.
#' For each name given, only the one AphiaRecord will be retrived. AphiaRecord with "accepted" status are preferred.
#' If  not present last entry will be taken which seems to result in best results.
#' 
#' @examples
#' taxon_names <- c( "Westwodilla caecula" ,  "Abra alba", "Chaetozone cf. setosa",  "Algae" )
#' w <- wormsbynames(taxon_names)
#' ## print unrecognized returns
#' failed_species <- rownames(w[is.na(w[,1]),])
#' 
#' ## try again with fuzzy matching turned on
#' \donttest{w <- wormsbynames(taxon_names, match=TRUE)}
#' 
#' ## this is how to load taxon_names from file
#' write.csv(taxon_names , file = "tax.csv", 
#'         row.names = FALSE,
#'         na = "")
#' ## check it out, then load it
#' read.csv(file = "tax.csv",
#'         na = "", 
#'         stringsAsFactors = FALSE,
#'         col.names = FALSE)
#' ## save results to file to inspect with, e.g. spreadsheet software
#' write.csv(w,file = "aphiainfo.csv", 
#'         na = "", 
#'         col.names = TRUE,
#'         row.names = TRUE)
#' 
#'         
#' @export
wormsbynames <- function(taxon_names,ids=FALSE,match=FALSE,verbose=TRUE,chunksize=50,like="false", marine_only="true",sleep_btw_chunks_in_sec=0.1){
  #library(httr)
  #library(plyr)
  
  stopifnot(inherits(taxon_names,"character"))
  if(match){ids<-TRUE}
  search_options<-paste0("like=",like,"&marine_only=",marine_only)
  my_worms<-list()
  request<-"http://www.marinespecies.org/rest/AphiaRecordsByNames"
  #li<-setNames(    as.list(c(a[[1]],"false","true")) , c( rep("scientificnames[]",length(a[[1]])),"like","marine_only" ) )     
  # r<-GET("http://www.marinespecies.org/rest/AphiaRecordsByName", query = list("scientificnames[]" = "Abietinaria%20abietina","scientificnames[]" = "Acanthocardia%20echinata"))
  
  wrapname<-gsub(" ", "%20", taxon_names)
  chunk<-split(wrapname, ceiling(seq_along(taxon_names)/chunksize))
  chunkid<-split(1:length(taxon_names), ceiling(seq_along(taxon_names)/chunksize))
  cat("REQUESTING ",length(taxon_names)," ITEMS BY NAME from World Register of Marine Species (www.marinespecies.org), ",format(Sys.time(), "%d/%m/%Y %X")," (CC-BY)\n",sep = "")
  for (round in 1:length(chunk)){
    
    if(verbose){
      cat(sprintf("%62s", paste0("chunk ",round,"/",length(chunk))),"\n")
    }
    m<-paste(paste(paste0("scientificnames[]=",chunk[[round]]),collapse="&"),search_options,sep="&")
    r<-GET(paste(request,m,sep="?"))
    Sys.sleep(sleep_btw_chunks_in_sec)
    stopifnot(r$status_code==200)
    
    # gather lists in master list
    r_parsed<-content(r,as="parsed")
    for (i in 1:length(r_parsed)){
      w_index<-unlist(chunkid[round])[i]
      if(length(r_parsed[[i]])==0){
        my_worms[[w_index]]<-NA
        cat(sprintf("%-46s       %-40s", taxon_names[w_index] , "no match"),"\n")
      } else if(length(r_parsed[[i]])==1){
        my_worms[[w_index]]<-r_parsed[[i]][[1]]
      } else {
        l<-length(r_parsed[[i]])
        my_worms[[w_index]]<-r_parsed[[i]][[l]]
        for(j in 1:l){
          if(r_parsed[[i]][[j]]$status == "accepted"){
            my_worms[[w_index]]<-r_parsed[[i]][[j]]
          }
        }
        # alternate representation
      }
    }
  }
  # pull dataframe out of master list
  non.null.list <- lapply(my_worms, lapply, function(x)ifelse(is.null(x), NA, x))
  worms<-rbind.fill(lapply(non.null.list, as.data.frame,stringsAsFactors = F))
  worms$NA.<-NULL
  if(ids){
    worms<-cbind(data.frame(id=1:nrow(worms) , name=taxon_names,stringsAsFactors = F),worms)
  }
  if (verbose) {cat("by names ........................................... DONE\n")}
  if(match){
    nonefound<-is.na(worms[,"AphiaID"])
    failed_species<-taxon_names[nonefound]
    if(length(failed_species)>0){
      failed_worms<-wormsbymatchnames(failed_species,verbose=verbose,ids=FALSE)
      worms[nonefound,c(F,F,rep(T,ncol(failed_worms)))]<-failed_worms
    } else {
      cat("  Nothing to match.\n")
    }
    
  } 
  
  
  
  return(worms)
}


















#' @title GET AphiaRecordsByMatchNames
#' @importFrom httr  GET content
#' @importFrom plyr rbind.fill
#' 
#' @description takes character vector with taxon names and retrives AphiaRecords from WoRMS 
#'
#' @param taxon_names character vector with names of taxa to look up.
#' @param chunksize only 50 taxa can be looked up per request, so request are split up into chunks (should be 50 or less)
#' @param verbose be verbose
#' @param ids add column "id" and "name" with running id and search names
#' @param marine_only Limit to marine taxa. Default=true
#' @param sleep_btw_chunks_in_sec pause between requests 
#' 
#' @return a data frame.
#' @details This function will take a character vector with taxon names, 
#' retrive AphiaRecords (CC-BY) from www.marinespecies.org using the 
#' GET /AphiaRecordsByName/{ScientificName} Method described at
#' http://www.marinespecies.org/rest/.
#' Results will be output to a data.frame with each row being a record.
#' For each name given, only the one AphiaRecord will be retrived. AphiaRecord with "accepted" status are preferred.
#' If  not present last entry will be taken which seems to result in best results.
#' 
#' @export
wormsbymatchnames <- function(taxon_names,verbose=TRUE,ids=FALSE,chunksize=50, marine_only="true",sleep_btw_chunks_in_sec=0.1){
  #library(httr)
  #library(plyr)
  stopifnot(inherits(taxon_names,"character"))
  
  
  search_options<-paste0("marine_only=",marine_only)
  my_worms<-list()
  request<-"http://www.marinespecies.org/rest/AphiaRecordsByMatchNames"
  #li<-setNames(    as.list(c(a[[1]],"false","true")) , c( rep("scientificnames[]",length(a[[1]])),"like","marine_only" ) )     
  # r<-GET("http://www.marinespecies.org/rest/AphiaRecordsByName", query = list("scientificnames[]" = "Abietinaria%20abietina","scientificnames[]" = "Acanthocardia%20echinata"))
  wrapname<-gsub(" ", "%20", taxon_names)
  chunk<-split(wrapname, ceiling(seq_along(taxon_names)/chunksize))
  chunkid<-split(1:length(taxon_names), ceiling(seq_along(taxon_names)/chunksize))
  if(verbose){
    cat("REQUESTING ",length(taxon_names)," ITEMS USING FUZZY from World Register of Marine Species  (www.marinespecies.org), ",format(Sys.time(), "%d/%m/%Y %X")," (CC-BY)\n",sep = "")
  }
  for (round in 1:length(chunk)){
    if(verbose){
      cat(sprintf("%62s", paste0("chunk ",round,"/",length(chunk))),"\n")
    }
    m<-paste(paste(paste0("scientificnames[]=",chunk[[round]]),collapse="&"),search_options,sep="&")
    r<-GET(paste(request,m,sep="?"))
    Sys.sleep(sleep_btw_chunks_in_sec)
    stopifnot(r$status_code==200)
    
    # gather lists in master list
    count<-0
    r_parsed<-content(r,as="parsed")
    for (i in 1:length(r_parsed)){
      w_index<-unlist(chunkid[round])[i]
      if(is.null(r_parsed[[i]][[1]])){
        my_worms[[w_index]]<-NA
        count<-count+1
        if(verbose){
          cat(sprintf("%-46s       %-40s", taxon_names[w_index] , "no match"),"\n")
        }
      } else {
        
        if(length(r_parsed[[i]])==1){
          my_worms[[w_index]]<-r_parsed[[i]][[1]]
        } else {
          my_worms[[w_index]]<-r_parsed[[i]][[1]]
          
          for(j in 1:length(r_parsed[[i]])){
            if(r_parsed[[i]][[j]]$status == "unaccepted"){
              my_worms[[w_index]]<-r_parsed[[i]][[j]]
            }
          }
          for(j in 1:length(r_parsed[[i]])){
            if(r_parsed[[i]][[j]]$status == "accepted"){
              my_worms[[w_index]]<-r_parsed[[i]][[j]]
            }
          }
        }
        if(verbose){
          cat(sprintf("%-49s -> %-40s %7s", taxon_names[w_index],my_worms[[w_index]]$scientificname,my_worms[[w_index]]$match_type),"\n")
        }
      }
    }
  }
  if(count<length(r_parsed)){
    # pull dataframe out of master list
    non.null.list <- lapply(my_worms, lapply, function(x)ifelse(is.null(x), NA, x))
    worms<-rbind.fill(lapply(non.null.list, as.data.frame,stringsAsFactors = F))
    worms$NA.<-NULL
    if(ids){
      worms<-cbind(data.frame(id=1:nrow(worms) , name=taxon_names,stringsAsFactors = F),worms)
    }
  } else {
    worms<-my_worms
  }
  if (verbose) {
    cat("matching names ..................................... DONE \n")
    cat("WARNING, fuzzy matching may give unexpected results - check results!\n\n")
  }
  return(worms)
}









#' @title GET AphiaRecordByAphiaID
#' @importFrom httr  GET content
#' @importFrom plyr rbind.fill
#' 
#' @description takes more than one AphiaID and retrives AphiaRecords from WoRMS 
#'
#' @param x AphiaIDs
#' @param verbose be verbose
#' @param ids add column "id" and "name" with running id and search names
#' @param sleep_btw_chunks_in_sec pause between requests 
#' 
#' @return a data frame.
#' @details This function will take a integer vector with AphiaIDs, 
#' retrive AphiaRecords from www.marinespecies.org using the 
#' GET /AphiaRecordByAphiaID Method described at
#' http://www.marinespecies.org/rest/.
#' Results will be output to a data.frame with each row being a record.
#' 
#' For examples, see  \code{\link{wormsaccepted}}
#' 
#' @export
wormsbyid <- function(x,verbose=TRUE,ids=FALSE,sleep_btw_chunks_in_sec=0.01){
  #library(httr)
  #library(plyr)
  stopifnot(inherits(x,c("numeric","integer")))
  my_worms<-list()
  if(verbose){      
    cat("REQUESTING ",length(x)," ITEMS BY ID from World Register of Marine Species  (www.marinespecies.org), ",format(Sys.time(), "%d/%m/%Y %X")," (CC-BY)\n",sep = "")  
  }
  for (round in 1:length(x)){
    if(verbose){      cat(",",x[round],sep = "")    }
    r<-GET(paste("http://www.marinespecies.org/rest/AphiaRecordByAphiaID",x[round],sep="/"))
    Sys.sleep(sleep_btw_chunks_in_sec)
    if (r$status_code!=200) {
      if(verbose){        cat("(",r$status_code,")",sep = "")      }
      r_parsed<-NA
    } else {
      #if(verbose){        cat("success (Status ",r$status_code,") \n",sep = "")      }
      r_parsed<-content(r,as="parsed")
    }
    if(verbose & round%%10==0){        cat("\n",sep = "")      }
    my_worms[[round]]<-r_parsed
  }
  
  
  if (verbose) {cat("\n")}
  # pull dataframe out of master list
  non.null.list <- lapply(my_worms, lapply, function(x)ifelse(is.null(x), NA, x))
  worms<-rbind.fill(lapply(non.null.list, as.data.frame,stringsAsFactors = F))
  if(ids){
    worms<-cbind(data.frame(id=1:nrow(worms) , name=x,stringsAsFactors = F),worms)
  }
  if (verbose) {cat("by id .............................................. DONE\n")}
  return(worms)
}


#' @title Recursivly retrieves respective "accepted" AphiaRecords for all synonyms if not already there
#' 
#' @description takes data.frame as output by \code{\link{wormsbynames}} , 
#' \code{\link{wormsbymatchnames}}, or \code{\link{wormsbyid}} and retrieves  additional
#' Aphia records (CC-BY) for not-"accepted" records in order to ultimately have "accepted" synonyms for all 
#' records in the dataset.
#'
#' @param x data.frame
#' @param verbose be verbose
#' @param sleep_btw_chunks_in_sec pause between requests 
#' @param once only one retrival iteration. No concatination of output with result. (For debugging)
#' 
#' @return a data frame.
#' @details This function will take a integer vector with AphiaIDs, 
#' retrive AphiaRecords from www.marinespecies.org using the 
#' GET /AphiaRecordByAphiaID Method described at
#' http://www.marinespecies.org/rest/.
#' Results will be outbut to a data.frame with each row being a record.
#' 
#' For examples, see  \code{\link{wormsaccepted}}
#'
#' @export
wormsconsolidate <- function(x,verbose=TRUE,sleep_btw_chunks_in_sec=0.01,once=FALSE){
  if(FALSE){
    x<-w
    verbose=TRUE
    sleep_btw_chunks_in_sec=0.01
  }
  count<-0
  while(TRUE){
    count<-count+1
    ids<-ifelse(names(x)[1]=="id",TRUE,FALSE)
    #cat(ids,"\n")
    
    
    unexplained<-x[x$valid_AphiaID>0 & (!is.na(x$valid_AphiaID)) & (x$status!="accepted") & (!x$valid_AphiaID%in%x$AphiaID),]
    anz_unexplaind<-nrow(unexplained)
    if(   anz_unexplaind!=0   ){
      if (verbose) {    cat("\nstill ",anz_unexplaind,"unhappy worms\n")   }
      happyworms<-wormsbyid(unexplained$valid_AphiaID,verbose=verbose) 
      if(ids){
        happyworms<-cbind(data.frame(id=(max(x$id)+1):(max(x$id)+nrow(happyworms)) , name=happyworms$scientificname), happyworms)
      }
      if(once){
        cat("rbind result with input and run maybe again!\n")
        return(happyworms)
      } else {
        x<-rbind(x,happyworms)
      }
    } else {
      if (verbose) {cat("consolidating ...................................... DONE\n")}
      break
    }
    
  }
  return(x)
}




#' @title Constructs "accepted_id" column which contains the "AphiaID" of the respective "accepted" taxon
#' 
#' @description takes data.frame as output by \code{\link{wormsbynames}} , 
#' \code{\link{wormsbymatchnames}}, or \code{\link{wormsbyid}} 
#' and add field "accepted_id" wich contains the "AphiaID" 
#' of the respective "accepted" taxon
#'
#' @param x data.frame
#' @param verbose be verbose
#' @param n_iter maximum search depth. Usually 3 is sufficient. Safety feature for breaking the \code{while} loop
#' 
#' @return a data frame.
#' @details This function helps updating you taxon information and eliminates ambiguity 
#' because the valid AphiaID is nor neccessary the AphiaID of an accepted taxon. You should run 
#' \code{\link{wormsconsolidate}} bevorhand to enshure all "accepted" taxons are present.
#' 
#' @examples 
#' ## start with IDs that are no longer up to date 
#' # get the Aphia information
#' u<-wormsbyid(c(424548,340537))
#' 
#' #recursively retrive information on the taxa they refer to
#' v<-wormsconsolidate(u)
#' 
#' # what are the currently correct "accepted" taxa? Answer: "accepted_id".
#' w<-wormsaccepted(v)
#' w[,c("scientificname","AphiaID","status","valid_AphiaID","valid_name","accepted_id")]
#' 
#' 
#' @export
wormsaccepted <- function(x,verbose=TRUE,n_iter=10){
  # fixme make shure valid_AphiaIDs are unique
  x$accepted_id<-NA 
  noacceptedid<-noacceptedrow<-NULL
  na_entries<-noval_entries<-0
  for(i in 1:nrow(x)){
    count<-0
    if(is.na(x$valid_AphiaID[i])) {
      na_entries<-na_entries+1
      if(!is.na(x$scientificname[i])){
        noval_entries<-noval_entries+1
        cat(x$scientificname[i] , " has no valid_AphiaID \n")
      }
      next
    }
    x$accepted_id[i]<-x$valid_AphiaID[i]
    if(is.na(match(x$accepted_id[i],x$AphiaID))){
      cat("ERROR! AphiaID",x$accepted_id[i],"is missing ... rerun wormsconsolidate and join with input")
      return()
    }
    while(x$status[match(x$accepted_id[i],x$AphiaID)] != "accepted"){
      count<-count+1
      x$accepted_id[i]<-x$valid_AphiaID[match(x$accepted_id[i],x$AphiaID)]
      if(count>n_iter){
        cat("no accepted AphiaID for ",x$scientificname[i] ,x$AphiaID[i], " after", n_iter," interations\n")
        noacceptedid <- cbind(noacceptedid,x$AphiaID[i])
        noacceptedrow <- cbind(noacceptedrow,i)
        break
      }
    }
  }
  if (verbose) {
    if(length(noacceptedid)>0){
      cat("AphiaID without acceptance ........................... ",length(noacceptedid),"\n")
      cat(paste(noacceptedid,","))
      cat("\n")
    }
    if(length(noacceptedrow)>0){
      cat("AphiaID without acceptance (rowsnumbers)............... \n")
      cat(paste(noacceptedrow,","))
      cat("\n")
    }
    if(na_entries>0){
    cat("other entries with no valid ID (NAs)............... ",na_entries-noval_entries,"\n")
    }
    cat("DONE ............................................... accepting\n")
    }
  return(x)
}










