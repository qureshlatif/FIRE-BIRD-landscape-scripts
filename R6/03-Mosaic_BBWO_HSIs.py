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
base = "T:/FS/RD/RMRS/Science/WTE/Research/Landscape_decision_tool/"
wildfire_data_dir = "wildfire_data/"
wildfire_years = range(1985, 1994) + range(1995, 2016)
HSI_dir = base + "R6/"
#_________________________________________#

for t in range(0, len(wildfire_years)):
    wf_dir = base + wildfire_data_dir + str(wildfire_years[t]) + '/'
    BB_polys = glob.glob(wf_dir + "brn_bounds/*.shp")
    BB_polys = [bb for bb in BB_polys if bb[89:91] == 'or' or bb[89:91] == 'wa']
    hsi_dir = HSI_dir + str(wildfire_years[t]) + '/'
    HSI_rasters = []
    for f in range(0, len(BB_polys)):
        fire = BB_polys[f][89:110]
        hsi = hsi_dir + fire + '/hsi.tif'
        HSI_rasters = HSI_rasters + [hsi]
    FP = arcpy.Merge_management(BB_polys, hsi_dir + "fire_perims.shp")
    HSI0 = arcpy.MosaicToNewRaster_management(HSI_rasters, hsi_dir, "HSI_unmask.tif", number_of_bands = 1)
    HSI = ExtractByMask(HSI0, FP)
    HSI.save(hsi_dir + "HSI.tif")
    del(hsi, HSI0, HSI)
    arcpy.Delete_management(hsi_dir + "HSI_unmask.tif")




















