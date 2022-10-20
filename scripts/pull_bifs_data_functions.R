# This file stores function to obtain BIP data from the Flowcentral

#Load required packages
require(tidyverse)
require(plyr)
require(plotly)

## Get the file path of the flowcell when provided with the flowcell ID
getFCPath <- function(fcid){
  ivd_path <- "/ghds/ivd/flowcentral/"
  if (is.na(fcid)){ return(NULL)}
  ivd <- list.files(path = ivd_path, pattern = fcid)[1]
  if(!is.na(ivd)) {
    return(paste(ivd_path, ivd, "/",sep = ""))
  }
  else {
    print(paste(fcid, "flowcell path not found!"))
    return(NULL)
  }
}

## Pull the autoqc data from BIP output with provided flowcell ID
pull_autoqc_info <- function(fcid_list){
  autoqc_final <- list()
  for(x in 1:length(fcid_list)){
    fcid <- getFCPath(as.character(fcid_list[[x]]))
    autoqc_file <- list.files(fcid, pattern = glob2rx("autoqc_sample_qc.hdr.tsv"), full.names = TRUE)
    autoqc_list <- sapply(autoqc_file, read.delim, simplify = FALSE, USE.NAMES = TRUE)
    autoqc_list <- lapply(autoqc_list, transform, run_sample_id = as.character(run_sample_id) )
    autoqc <- rbind.fill(autoqc_list,  .id = NULL)
    autoqc_final <- rbind(autoqc_final, autoqc)
  }
  return(autoqc_final)
}

## Pull the ghcnv data from BIP output with provided flowcell ID
pull_ghcnv_info <- function(fcid_list){
  ghcnv_final <- list()
  for(x in 1:length(fcid_list)){
    fcid <- getFCPath(as.character(fcid_list[[x]]))
    ghcnv_file <- list.files(fcid, pattern = glob2rx("*.ghcnv_qc.hdr.tsv"), full.names = TRUE)
    ghcnv_list <- sapply(ghcnv_file, read.delim, simplify = FALSE, USE.NAMES = TRUE)
    ghcnv_list <- lapply(ghcnv_list, transform, run_sample_id = as.character(run_sample_id))
    ghcnv <- bind_rows(ghcnv_list, .id = NULL)
    ghcnv_final  <- rbind(ghcnv_final, ghcnv)
  }
  return(ghcnv_final)
}

# Pull the coverage informaiton from BIP output with provided flowcell ID
pull_coverage_info <- function(fcid_list){
  coverage_final <- list()
  for(x in 1:length(fcid_list)){
    fcid <- getFCPath(as.character(fcid_list[[x]]))
    coverage_file <- list.files(fcid, pattern = glob2rx("*.coverage.hdr.tsv"), full.names = TRUE)
    coverage_list <- sapply(coverage_file, read.delim, simplify = FALSE, USE.NAMES = TRUE)
    coverage_list <- lapply(coverage_list, transform, run_sample_id = as.character(run_sample_id) )
    coverage <- rbind.fill(coverage_list,  .id = NULL)
    coverage_final <- rbind(coverage_final, coverage)
  }
  return(coverage_final)
}

## Pull ontarget db data from BIP output
pull_otdb_info <- function(fcid_list){
  otdb_final <- list()
  for(x in 1:length(fcid_list)){
    fcid <- getFCPath(as.character(fcid_list[[x]]))
    otdb_file <- list.files(fcid, pattern = glob2rx("*.on_target_db.hdr.tsv"), full.names = TRUE)
    otdb_list <- sapply(otdb_file, read.delim, simplify = FALSE, USE.NAMES = TRUE)
    otdb_list <- lapply(otdb_list, transform, run_sample_id = as.character(run_sample_id))
    otdb <- bind_rows(otdb_list, .id = NULL)
    otdb_final <- rbind(otdb_final, otdb)
  }
  return(otdb_final)
}


## ggplot theme for nice looking/publication ready Plots 
# https://rpubs.com/Koundy/71792
theme_Publication <- function(base_size=14, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
}




