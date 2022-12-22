# text mining pdf's 

install.packages("pdftools")
library(pdftools)
library(stringr)

setwd("/Volumes/GH Quality/Scanned Records/Materials Management/QC Forms Completed")

PN <- 'GHM0122'
lot_list <- c('2014', '658045', '658057', '658287')
# Get list of PN's within path
path <- "."
pdf_list <- list.files(path = path, pattern = PN)

pdfs <- list()
for(lot in lot_list){
  temp <- pdf_list[str_detect(pdf_list, lot)]
  pdfs <- append(pdfs, temp)
}
pdfs <- unlist(pdfs)
