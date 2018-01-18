#########################################################################################
# Removes values from cmort.tif layers where RdNBR values are missing due to scan lines #
#########################################################################################

library(raster)
library(rgdal)
library(stringr)
library(dplyr)

setwd("T:/FS/RD/RMRS/Science/WTE/Research/Landscape_decision_tool/wildfire_data/")

years <- c(2004, 2007, 2008, 2009, 2012, 2013, 2014)

for(i in 1:length(years)) {
  files <- list.files(paste0(years[i], "/RdNBR/lines"))
  files <- files[which(str_detect(files, ".tif"))]
  ID <- files %>% str_sub(1, 39) %>% unique()
  for(j in 1:length(ID)) {
    cmort <- raster(paste0(years[i], "/RdNBR/lines/", ID[j], "_cmort.tif"))
    rdnbr <- raster(paste0(years[i], "/RdNBR/", ID[j], "_rdnbr.tif"))
    scan.lines <- getValues(rdnbr) == -32768
    cmort[which(scan.lines)] <- NA
    writeRaster(cmort, paste0(years[i], "/RdNBR/lines/", ID[j], "_cmort.tif"), format = "GTiff", overwrite = T)
  }
}
