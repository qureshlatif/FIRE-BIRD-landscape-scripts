####################################################################################
# Generates canopy cover layers for missing years using canopy mortality from fire #
# Imputes up to 10 years prior to earliest available layer from LANDFIRE           #
####################################################################################

library(raster)
library(rgdal)

setwd("T:/FS/RD/RMRS/Science/WTE/Research/Landscape_decision_tool/")

ccov.yrs <- c( seq(2002, 2007), 2009, 2011, 2013)
                  # years for which canopy cover needs to be generated
ccov.init <- c(seq(2001, 2006), 2008, 2010, 2012)
                  # can cov layer by which to multiply or divide cmort layer
cmort.yrs <- c(seq(2003, 2008), 2010, 2012, 2014)
                  # years for canopy mortality layers

for(i in 1:length(ccov.yrs)) {
  cmort <- raster(paste0("wildfire_data/", cmort.yrs[i], "/ccmort_mask.tif"))
  ccov.st <- raster(paste0("cancov_data/", ccov.init[i], "/ccov.tif"))
  ccov.new <- raster(ccov.st)
  pth <- paste0(getwd(), "/cancov_data/", ccov.yrs[i])
  if(!dir.exists(pth)) dir.create(pth)
  ccov_write <- writeStart(ccov.new, filename=
                             paste0("cancov_data/", ccov.yrs[i], "/ccov.tif"),
                         format = "GTiff", dataType = "INT1U", overwrite = T, update = T)
  nr<-round(nrow(ccov.new)/7)
  for (z in 0:6) {
    st <- z*nr+1
    if(z == 6) nr <- nrow(ccov.new) - (st + 1)
    cmvec <- getValues(cmort, row=st, nrows=nr)
    cmvec[which(cmvec == 128)] <- 0
    ccvec <- getValues(ccov.st, row=st, nrows=nr)
    vals <- round(ccvec * (1 - cmvec/100))
    ccov.new <- writeValues(ccov_write, v=vals, start=st)
  }
  ccov_write <- writeStop(ccov.new)
}
