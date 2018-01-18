######################
# Extract .tar files #
######################

setwd("T:/FS/RD/RMRS/Science/WTE/Research/Landscape_decision_tool/wildfire_data/tars")
library(utils)
library(stringr)

files <- list.files()
files <- str_subset(files, ".tar.gz")

for(i in 1:length(files)) untar(files[i])
