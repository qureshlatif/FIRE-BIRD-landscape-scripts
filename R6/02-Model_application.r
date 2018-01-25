#########################################################################################################
# This script applies ensemble HSI model predictions for BBWO (Latif et al. 2013) to data for Region 6. #
# The result is HSI layers for each fire in each year.                                                  #
# A subsequent Python script mosaics data for all fires in a given year together.                       #
# Note: Initially tried doing this in Python and ran into memory issues, so switched to R.              #
# ***DO NOT DO THIS PART IN PYTHON***                                                                   #
#########################################################################################################

library(R.utils)
library(raster)

setwd("T:/FS/RD/RMRS/Science/WTE/Research/Landscape_decision_tool/R6")    #sets the working directory which should contain the file "Mahalanobis_source.r"
nest.samples <- loadObject("Nest_samples_mahalanobis") # Nest samples for Mahalanobis models

###############################
# Model application functions #
###############################

expit <- function(x) exp(x)/(1+exp(x))

Mahal <- function (Cal,Lnd) {
  HSI <- numeric(length=nrow(Lnd))
  
  mns <- apply(Cal,2,mean)
  sds <- apply(Cal,2,sd)
  cal <- t(apply(Cal,1,function(x) (x - mns)/sds))
  lnd <- t(apply(Lnd,1,function(x) (x - mns)/sds))
  
  cnt <- apply(cal,2,mean)
  cv <- cov(cal)
  dgf <- ncol(cal)          
  D.lnd <- mahalanobis(lnd,cnt,cv)
  p.lnd <- 1-pchisq(D.lnd,dgf)
  HSI <- p.lnd
  return(HSI)
}

Mxnt3v_scr <- function(cosasp,dnbr,loccc40) {
  #Define features and parameters
  l.cosasp <- c(0.9381002574647688, -1, 1)
  
  l.dnbr <- c(10.786746325527512, -590.47619629, 1024.12268066)
  dnbr[which(dnbr>l.dnbr[3])]<-l.dnbr[3]
  dnbr[which(dnbr<l.dnbr[2])]<-l.dnbr[2]
  
  l.loccc40 <- c(1.1312580836869854, 0.0, 1.0)
  
  l.cosasp2 <- c(-0.3980392534351936, 0.0, 1.0)
  cosasp2 <- cosasp^2
  
  l.dnbr2 <- c(-6.415155660077871, 0.0, 1048827.2650422244)
  dnbr2 <- dnbr^2
  dnbr2[which(dnbr2>l.dnbr2[3])]<-l.dnbr2[3]
  dnbr2[which(dnbr2<l.dnbr2[2])]<-l.dnbr2[2]
  
  l.loccc402 <- c(-1.9494267010037454, 0.0, 1.0)
  loccc402 <- loccc40^2
  
  l.ca_x_lcc40 <- c(-0.08235355467277128, -1.0, 1.0)
  ca_x_lcc40 <- cosasp*loccc40
  ca_x_lcc40[which(ca_x_lcc40>l.ca_x_lcc40[3])]<-l.ca_x_lcc40[3]
  ca_x_lcc40[which(ca_x_lcc40<l.ca_x_lcc40[2])]<-l.ca_x_lcc40[2]
  
  l.dnbr_x_lcc40 <- c(3.6550086289032, -213.7978722982045, 1024.0771640968703)
  dnbr_x_lcc40 <- dnbr*loccc40
  dnbr_x_lcc40[which(dnbr_x_lcc40>l.dnbr_x_lcc40[3])]<-l.dnbr_x_lcc40[3]
  dnbr_x_lcc40[which(dnbr_x_lcc40<l.dnbr_x_lcc40[2])]<-l.dnbr_x_lcc40[2]
  
  linPN <- 8.12466828053413
  densNorm <- 2846.1535672273485
  entropy <- 8.807545807433492
  
  exponent <- (l.cosasp[1]*((cosasp-l.cosasp[2])/(l.cosasp[3]-l.cosasp[2])) +
                 l.dnbr[1]*((dnbr-l.dnbr[2])/(l.dnbr[3]-l.dnbr[2])) +
                 l.loccc40[1]*((loccc40-l.loccc40[2])/(l.loccc40[3]-l.loccc40[2])) +
                 l.cosasp2[1]*((cosasp2-l.cosasp2[2])/(l.cosasp2[3]-l.cosasp2[2])) +
                 l.dnbr2[1]*((dnbr2-l.dnbr2[2])/(l.dnbr2[3]-l.dnbr2[2])) +
                 l.loccc402[1]*((loccc402-l.loccc402[2])/(l.loccc402[3]-l.loccc402[2])) +
                 l.ca_x_lcc40[1]*((ca_x_lcc40-l.ca_x_lcc40[2])/(l.ca_x_lcc40[3]-l.ca_x_lcc40[2])) +
                 l.dnbr_x_lcc40[1]*((dnbr_x_lcc40-l.dnbr_x_lcc40[2])/(l.dnbr_x_lcc40[3]-l.dnbr_x_lcc40[2]))
  ) - linPN
  mx.raw  <- exp(exponent)/densNorm
  HSI <- (mx.raw*exp(entropy))/(1+mx.raw*exp(entropy))
  return(HSI)
}

Mxntbrn_scr <- function(dnbr) {
  l.dnbr <- c(5.393373689029553, -590.47619629, 1024.12268066)
  dnbr[which(dnbr>l.dnbr[3])] <- l.dnbr[3]
  dnbr[which(dnbr<l.dnbr[2])] <- l.dnbr[2]
  linPN <- 5.393373689029553
  densNorm <- 1369.5883369712144
  entropy <- 8.899831912829168
  exponent <- l.dnbr[1]*((dnbr-l.dnbr[2])/(l.dnbr[3]-l.dnbr[2])) - linPN
  mx.raw  <- exp(exponent)/densNorm
  HSI <- (mx.raw*exp(entropy))/(1+mx.raw*exp(entropy))
  return(HSI)
}

#Function for Maxent and WLR HSIs
HSI_calc <- function(cosasp,dnbr,loccc40,lndcc40) {
  SG_wlr <- expit(-3.130234119 + 0.003704316*dnbr + 3.154616050*lndcc40)
  TP_wlr <- expit(-2.485159177 + 0.006151492*dnbr)
  TB_wlr <- expit(-0.743163 + 0.001546*dnbr)
  Maxent_3v <- Mxnt3v_scr(cosasp,dnbr,loccc40)
  Maxent_brn <- Mxntbrn_scr(dnbr)
  HSIs<-cbind(SG_wlr,TP_wlr,TB_wlr,Maxent_3v,Maxent_brn)
  return(HSIs)
}

###################################
# Iterate by year to apply models #
###################################
years <- c(1985:1993, 1995:2015)

for(yr in 1:length(years)) {
  yr_dir <- paste0(years[yr], "/")
  fire_dirs <- dir(yr_dir)
  for(f in 1:length(fire_dirs)){
    ## Load input raster layers ##
    CASP.r <- raster(paste0(yr_dir, fire_dirs[f], "/INPUTS/cosasp.tif"))
    DNBR.r <- raster(paste0(yr_dir, fire_dirs[f], "/INPUTS/locdnbr.tif"))
    LOCCC.r <- raster(paste0(yr_dir, fire_dirs[f], "/INPUTS/loccc.tif"))
    LANDCC.r <- raster(paste0(yr_dir, fire_dirs[f], "/INPUTS/landcc.tif"))
    
    HSI.r <- raster(CASP.r) # Copy layer to receive output values
    
    ## Get input values for calculating HSIs
    ind.vals <- which(!is.na(getValues(CASP.r)))
    COSASP <- getValues(CASP.r)[ind.vals]
    DNBR <- getValues(DNBR.r)[ind.vals]
    LOCCC40 <- getValues(LOCCC.r)[ind.vals]
    LANDCC40 <- getValues(LANDCC.r)[ind.vals]
    
    lndscp <- cbind(COSASP, DNBR, LOCCC40, LANDCC40)
    
    BC1.HSI <- TBC1.HSI <- TB1.HSI <- numeric(length=nrow(lndscp)) 
    for (y in 1:length(nest.samples)) {
      Recs <- nest.samples[[y]]
      #BC1 model ################################################################################
      Lnd <- as.matrix(lndscp[,c("DNBR","LOCCC40","LANDCC40")])          # Insert variable names here.
      Loc <- as.matrix(Recs[,c("DNBR","LOCCC40","LANDCC40")])     # Insert the same variable names as above.
      M <- Mahal(Loc,Lnd)
      BC1.HSI <- BC1.HSI + round(M,digits=4)
      #TBC1 model ################################################################################
      Lnd<-as.matrix(lndscp[,c("COSASP","DNBR","LOCCC40","LANDCC40")])          # Insert variable names here.
      Loc<-as.matrix(Recs[,c("COSASP","DNBR","LOCCC40","LANDCC40")])     # Insert the same variable names as above.
      M<-Mahal(Loc,Lnd)
      TBC1.HSI <- TBC1.HSI + round(M,digits=4)
      #TB1 model ################################################################################
      Lnd<-as.matrix(lndscp[,c("COSASP","DNBR")])          # Insert variable names here.
      Loc<-as.matrix(Recs[,c("COSASP","DNBR")])     # Insert the same variable names as above.
      M<-Mahal(Loc,Lnd)
      TB1.HSI <- TB1.HSI + round(M,digits=4)
    }
    #Calculate HSI means and s.e.'s across model replicates
    BC1.HSI <- BC1.HSI/100
    TBC1.HSI <- TBC1.HSI/100
    TB1.HSI <- TB1.HSI/100
    Dat <- cbind(BC1.HSI,TBC1.HSI,TB1.HSI)   #Bind matrix with HSI scores.
    dimnames(Dat)[[2]] <- c("BC1","TBC1","TB1")
    
    
    HSIs <- HSI_calc(COSASP, DNBR, LOCCC40, LANDCC40)
    Dat <- cbind(Dat, HSIs)
    rm(HSIs)
    
    ## Classify cells using the max-gain thresholds and compile ensemble predictions ##
    
    mods <- dimnames(Dat)[[2]][c(1:8)]
    thrhds <- c(0.17, 0.17, 0.32, 0.43, 0.43, 0.45, 0.41, 0.37)
    names(thrhds) <- mods
    
    for (j in 1:length(mods))
      Dat[,j] <- as.numeric(Dat[, mods[j]] >= thrhds[mods[j]])
    HSI <- apply(Dat, 1, sum)
    HSI.r[ind.vals] <- HSI
    writeRaster(HSI.r, paste0(yr_dir, fire_dirs[f], "/hsi.tif"), format = "GTiff", overwrite = T)
  }
}
