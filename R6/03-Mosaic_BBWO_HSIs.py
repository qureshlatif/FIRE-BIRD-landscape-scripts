# Ran into consistent memory error in 2003 at or4244112390420020713. Oddly, the BBWO model application tool ran fine with the data for this fire.
# Switched to R version of this script to see if memory errors can be avoided.

import arcpy
from arcpy import env
from arcpy.sa import *
arcpy.CheckOutExtension("Spatial")
import os, sys
import glob
print 'IMPORTED'
arcpy.env.overwriteOutput = True

#_________________Inputs__________________#
# Base directory #
base = "T:/FS/RD/RMRS/Science/WTE/Research/Landscape_decision_tool/R6/"
years = range(1985, 1994) + range(1995, 2016)
#_________________________________________#

for t in range(0, len(years)):
    yr_dir = base + str(years[t]) + '/'
    fire_dirs = os.listdir(yr_dir)
    HSI_rasters = []
    for f in range(0, len(fire_dirs)):
        hsi = yr_dir + fire_dirs[f] + '/hsi.tif'
        HSI_rasters = HSI_rasters + [hsi]
    arcpy.MosaicToNewRaster_management(HSI_rasters, yr_dir, "HSI.tif", number_of_bands = 1)
    del(hsi)




















