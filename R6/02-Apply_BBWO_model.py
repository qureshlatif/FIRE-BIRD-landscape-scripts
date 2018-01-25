# Ran into consistent memory error in 2003 at or4244112390420020713. Oddly, the BBWO model application tool ran fine with the data for this fire.
# Switched to R version of this script to see if memory errors can be avoided.

import arcpy
from arcpy import env
from arcpy.sa import *
arcpy.CheckOutExtension("Spatial")
import os, sys
import glob
import scipy
from string import *
from scipy import spatial
from scipy import stats
from numpy import *
from scipy import special
print 'IMPORTED'
arcpy.env.overwriteOutput = True

from BBWO_HSI import *
mHSI = Mahalonobis_HSI()

#_________________Inputs__________________#
# Base directory #
base = "T:/FS/RD/RMRS/Science/WTE/Research/Landscape_decision_tool/R6/"
years = range(1985, 1994) + range(1995, 2016)
#_________________________________________#

# function(s)
def checkPath(path):
    try: os.makedirs(path)
    except: pass

#Set workspace and scratch workspace for temp outputs/calculations
arcpy.env.workspace = base + "/_workspace"
checkPath(arcpy.env.workspace)
arcpy.env.scratchWorkspace = base + "/_scratch"
checkPath(arcpy.env.scratchWorkspace)

thrhds = array([0.17,0.17,0.32,0.43,0.43,0.45,0.41,0.37]) # Thresholds for classifying HSIs
# Note: Unfortunately, this loop will likely run into memory errors a couple times before completing. Restarting Python appears to
    # clear the memory and allow resumption of the loop. Unclear why these errors are arising and how to fix them.
for t in range(0, len(years)):
    yr_dir = base + str(years[t]) + '/'
    fire_dirs = os.listdir(yr_dir)
    HSI_rasters = []
    for f in range(0, len(fire_dirs)):
        INPUT_dir = yr_dir + fire_dirs[f] + '/INPUTS/'
        outputrast = yr_dir + fire_dirs[f] + '/hsi.tif'
        cosasp = r'' + INPUT_dir + 'cosasp'
        dnbr = r'' + INPUT_dir + 'locdnbr'
        loccc40 = r'' + INPUT_dir + 'loccc'
        landcc40 = r'' + INPUT_dir + 'landcc'
        ## Set environmental information
        BOTTOM  =arcpy.GetRasterProperties_management(cosasp, 'BOTTOM')
        BOTTOM  =BOTTOM.getOutput(0)
        LEFT = arcpy.GetRasterProperties_management(cosasp, 'LEFT')
        LEFT = LEFT.getOutput(0)
        lower_left_corner = arcpy.Point(LEFT, BOTTOM)
        COLUMNCOUNT= arcpy.GetRasterProperties_management(cosasp, 'COLUMNCOUNT')
        COLUMNCOUNT=int(COLUMNCOUNT.getOutput(0))
        ROWCOUNT= arcpy.GetRasterProperties_management(cosasp, 'ROWCOUNT')
        ROWCOUNT=int(ROWCOUNT.getOutput(0))
        sr = arcpy.Describe(cosasp).spatialReference
        # Convert Rasters to numpy arrays for use in statistical model...
        cosaspNP = arcpy.RasterToNumPyArray(cosasp, '', '', '', -9999.00)
        cosaspNP = reshape(cosaspNP, size(cosaspNP))
        dnbrNP = arcpy.RasterToNumPyArray(dnbr, '', '', '', -9999.00)
        dnbrNP = reshape(dnbrNP, size(dnbrNP))
        locccNP = arcpy.RasterToNumPyArray(loccc40, '', '', '', -9999.00)
        locccNP = reshape(locccNP, size(locccNP))
        landccNP = arcpy.RasterToNumPyArray(landcc40, '', '', '', -9999.00)
        landccNP = reshape(landccNP, size(landccNP))
        ### load nest_samples dimenstions: 0 = variable, 1 = record, 2 = sample
        #   VARIABLES: 0 = COSASP, 1 = DNBR, 2 = LOCCC40, 3 = LANDCC40
        dirnm = os.path.dirname(os.path.realpath(sys.argv[0]))
        nest_samples = load(dirnm+'\\nest_samples.npy') # Copy from FIRE-BIRD toolbox into this folder.
        ## Calculate BC1, TBC1 and TB1 HSI model values...
        BC1_HSI = zeros(len(cosaspNP))
        TBC1_HSI = zeros(len(cosaspNP))
        TB1_HSI = zeros(len(cosaspNP))
        for i in range(0,100):
            print 'Generating Sample ' + str(i+1) + ' out of 100...'
            arcpy.AddMessage('Generating Sample ' + str(i+1) +  ' out of 100...')
            Recs = nest_samples[:,:,i]
            ### BC1 Model
            Lnd = array([dnbrNP, locccNP, landccNP]).T
            Loc = copy(Recs[1:, :].T)
            M = mHSI.Mahal(Loc, Lnd)
            BC1_HSI += M.round(4)
            del(Lnd, Loc, M)
            ## TBC1 Model
            Lnd = array([cosaspNP, dnbrNP, locccNP, landccNP]).T
            Loc = copy(Recs[:, :].T)
            M = mHSI.Mahal(Loc, Lnd)
            TBC1_HSI += M.round(4)
            del(Lnd, Loc, M)
            ## TB1 model
            Lnd = array([cosaspNP, dnbrNP]).T
            Loc = copy(Recs[0:2, :].T)
            M = mHSI.Mahal(Loc, Lnd)
            TB1_HSI += M.round(4)
            del(Lnd, Loc, M)
        del(Recs)
        BC1_HSI /= 100.
        TBC1_HSI /= 100.
        TB1_HSI /= 100.
        HSIs = mHSI.HSI_calc(cosaspNP, dnbrNP, locccNP, landccNP)
        Model_data =  vstack([BC1_HSI,TBC1_HSI,TB1_HSI, HSIs]).T   #Bind cell IDs and spatial coordinates to matrix with HSI scores.
        tab = zeros([size(Model_data, 0),size(Model_data, 1)+1], dtype = 'int32')
        for j in range(0, size(Model_data, 1)):
            tab[:,j] = array(Model_data[:,j] >= thrhds[j], dtype = int).T
        tab[:,-1] = sum(tab[:,0:8], axis = 1)
        sumall = tab[:,-1]
        sumall[cosaspNP==-9999]=-9999
        HSI = reshape(sumall, (ROWCOUNT, COLUMNCOUNT))
        HSIrast = arcpy.NumPyArrayToRaster(HSI, lower_left_corner, 30.0,30.0, -9999.00)
        arcpy.DefineProjection_management(HSIrast, sr)
        HSIrast.save(outputrast)
        HSI_rasters = HSI_rasters + [outputrast]
    arcpy.MosaicToNewRaster_management(HSI_rasters, yr_dir, "HSI.tif", number_of_bands = 1)
    del(cosaspNP, dnbrNP, locccNP, landccNP, BC1_HSI, TBC1_HSI, TB1_HSI, Model_data, nest_samples, tab, HSIs)




















