##Eleen Shum
##2019-03-25
##edited 2019-06-15 to add fragment isize
##edited v2 2019-06-25 to add mapd

##function to pull data of variant control with definied fcid based on whether the pipeline id's it as control
require(tidyverse)
require(jsonlite)
require(plyr)

## Get the path of the flowcell when providing the flowcell ID
getFCPath <- function(fcid){
  flowcentral.path <- "/ghds/flowcentral/"
  ivd.path <- "/ghds/ivd/flowcentral/"
  #return null if input is NA due to logging errors
  if (is.na(fcid)){ return(NULL)}
  flow <- list.files(path = flowcentral.path, pattern = fcid)[1]
  ivd <- list.files(path = ivd.path, pattern = fcid)[1]
  if(is.na(flow) & is.na(ivd)) {
    print(paste(fcid, "flowcell path not found!"))
    return(NULL)
  }
  if(!is.na(ivd)){
    return(paste(ivd.path, ivd, "/",sep = ""))
  }
  if(!is.na(flow)){
    return(paste(flowcentral.path, flow, "/", sep = ""))
  }
} 

getSampleMetrics <- function(fcid){
  if (is.null(getFCPath(fcid))){return(NULL)}
  fcPath <- paste(getFCPath(fcid), '/', sep = '')
  print(fcPath)
  outTable <- data.frame()
  
  if(!file.exists(paste(fcPath, 'autoqc_report.json', sep = ''))){return(NULL)}
  autoqc <- jsonlite::fromJSON(paste(fcPath, 'autoqc_report.json', sep = ''))
  
  for (i in 1:length(autoqc$samples)){
    if(is.null(autoqc$sample[[i]]$is_control)){return(NULL)}
    
    #check if bip recognizes sample as control
    if(autoqc$samples[[i]]$is_control == F){
      ghcnv <-  read.table(paste(fcPath, names(autoqc$sample)[i], '.ghcnv_qc.hdr.tsv', sep = ''), sep = '\t',header = T)
      metrics <-  autoqc$samples[[i]]$metrics[,c(1,4)]
      metrics <- rbind(metrics,
                       data.frame(metric = 'mapd',
                                  value = ghcnv$mapd))
    }
    if(autoqc$samples[[i]]$is_control == T){
      nsc <- read.table(paste(fcPath, names(autoqc$sample)[i], '.coverage.hdr.tsv', sep = ''), sep = '\t',header = T)
      onTarg <- read.table(paste(fcPath, names(autoqc$sample)[i], '.on_target_db.hdr.tsv', sep = ''), sep = '\t',header = T)
      ghcnv <-  read.table(paste(fcPath, names(autoqc$sample)[i], '.ghcnv_qc.hdr.tsv', sep = ''), sep = '\t',header = T)
      tmb <- read.table(paste(fcPath, names(autoqc$sample)[i], '.tmb_call.hdr.tsv', sep = ''), sep = '\t',header = T)
      umi <- jsonlite::fromJSON(paste(fcPath, names(autoqc$sample)[i], '.trimming_stats.json', sep = ''))
      umi.count <- c()
      for (j in 1:length(umi$tags)){
        umi.count <- c(umi.count, umi$tags[[j]][['count']])
      }
      metrics <-  autoqc$samples[[i]]$metrics[c(1:9),c(1,4)]
      metrics <- rbind(metrics,
                       data.frame(metric = 'tmb_score',
                                  value = tmb$tmb_score),
                       data.frame(metric = 'umi_cv',
                                  value = round(sd(umi.count)/mean(umi.count)*100,2)),
                       data.frame(metric = 'sample_non_singleton_families',
                                  value = round(median(nsc$median_coverage),2)),
                       data.frame(metric = 'sample_on_target_rate',
                                  value = round(onTarg$on_target_rate*100,2)),
                       data.frame(metric = 'sample_gc_bias',
                                  value = round(ghcnv$umol_gc_iqr,2)),
                       data.frame(metric = 'mapd',
                                  value = ghcnv$mapd))
    }
    metrics <- rbind(metrics, fragHist(names(autoqc$sample)[i], fcPath))
    x <- data.frame('fcid' = fcid,
                    'sampleName' = names(autoqc$sample)[i],
                    metrics)
    outTable <- rbind(outTable, x)
  }
  outTable <- spread(outTable, key = metric, value = value)
  outTable$sampleName <- as.character(outTable$sampleName)
  outTable <- outTable[order(outTable$sampleName),]
  return(outTable)
}
