# Insert LotID and Tapestation Run ID
# !!! File must be local with no duplicates !!!

TSRunID = "EIO-5291-P1_kk10" # change to your Run ID. Can be partial, as long as it's unique.

### 
list.of.packages <- c("tidyverse")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library(tidyverse)
###

# Find File
file_match <- list.files(path="/Users", pattern= TSRunID, recursive=TRUE, all.files=FALSE, full.names = TRUE, ignore.case=TRUE) %>% 
  stringr::str_subset(., "Library", negate = TRUE)
filepath <- grep(pattern="_compactPeakTable.csv", file_match, value=TRUE, ignore.case = TRUE)

fileName <- read.csv(filepath, header=FALSE)
names(fileName) <- fileName[1,]
fileName <- fileName[-1,]

# Select data
data <- fileName %>% select(2,4,5,7,8)
names(data) <- c("Well", "BP","Concentration", "Peak_Molarity", "Integrated")
data <- type.convert(data, as.is=TRUE)

# Remove any pattern EL* in data$Well
data <- data[-grep("EL", data$Well),]

data <- data[!(is.na(data$Integrated)), ] # Remove Upper/Lower Markers
data <- data[1:4]

# Heatmap Dataframe
heatmap_data <- data %>% select(1,2,4) %>% group_by(Well) %>% filter(BP >100 & BP <200) %>% filter(Peak_Molarity== max(Peak_Molarity))
write.csv(heatmap_data, paste0("~/Desktop/", TSRunID, "_filtered_data.csv"), row.names = FALSE)

# Print data frame
heatmap_data 
