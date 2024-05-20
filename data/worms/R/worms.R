#' @title worms: Use WoRMS RESTful webservice to scrape Aphia information from World Register of Marine Species
#' 
#' @description The worms package provides two kinds of functions: 
#' 
#' a) retrieving taxonomic information using WoRMS' RESTful Webservice 
#' by using taxon name search, fuzzy matching, or  Aphia ID search
#' implementing methods documented at \url{http://www.marinespecies.org/rest/}
#' 
#' b) functions that parse the data for synonyms in order to complete the dataset so that
#' for every taxon in the dataset the respective taxon with status 'accepted' exists as well.
#' Constructed references to the respective taxon with status 'accepted' help aggregating biodiversity data 
#' without the use of synonyms, alternative representations, and common misspellings leading to errors.
#' 
#' Check out \url{https://github.com/janhoo/worms/} for the developement version.
#' 
#' @references 
#' This package is not connected or endorsed by WoRMS.
#' According to \href{http://www.marinespecies.org}{WoRMS},
#' information from World Register of Marine Species is free to use under the condition that 
#' they are cited (CC-BY). 
#' While no license model is specified for the webservice employed, 
#' we strongly recommended to give reference to WoRMS, e.g., www.marinespecies.org, 18/06/17 (CC-BY).
#' The citation for the full database is:  
#' 
#' WoRMS Editorial Board (2017). World Register of Marine Species. 
#' Available from http://www.marinespecies.org at VLIZ. 
#' Accessed <today>. 
#' doi:10.14284/170  
#' 
#' For single taxa, references are given in the citation column
#' Please give proper reference to them.
#' 
#' @keywords worms 
#' @docType package
#' @name worms
NULL