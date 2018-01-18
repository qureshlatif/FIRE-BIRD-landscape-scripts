########################################################################
# Collects and processes cancov layers needed to impute missing values #
# in cmort.tif layers with scan lines.                                 #
########################################################################

# Note: Need to find and move (by hand) dNBR layers with scan lines to the "dNBR\\lines" folder before running this script.                 

import arcpy, os
from arcpy import env
import numpy
from numpy import *
from arcpy.sa import *

file_workspace = "T:\\FS\\RD\\RMRS\\Science\\WTE\\Research\\Landscape_decision_tool\\wildfire_data\\"

years = [2004, 2007, 2008, 2009, 2012, 2013, 2014]
years = map(str, years)

arcpy.CheckOutExtension("Spatial")
for i in range(0, len(years)):
    env.workspace = file_workspace + years[i] + "\\dNBR\\lines\\"
    arcpy.env.overwriteOutput = True # Set to overwrite pre-existing files
    DNBR = arcpy.ListRasters("*_dnbr.tif", "TIF")
    for layer in DNBR:
        env.workspace = file_workspace + years[i] + "\\dNBR\\lines\\"
        arcpy.env.overwriteOutput = True # Set to overwrite pre-existing files
        sr = arcpy.Describe(layer).spatialReference
        env.outputCoordinateSystem = sr
        ## Need to convert to numpy area in order to replace scan line with NoData values ##
        BOTTOM  =arcpy.GetRasterProperties_management(layer, 'BOTTOM')
        BOTTOM  =BOTTOM.getOutput(0)
        LEFT = arcpy.GetRasterProperties_management(layer, 'LEFT')
        LEFT = LEFT.getOutput(0)
        lower_left_corner = arcpy.Point(LEFT, BOTTOM)
        M = arcpy.RasterToNumPyArray(layer, '', '', '', -32768.00)
        R = arcpy.NumPyArrayToRaster(M, lower_left_corner, 30.0,30.0, -32768.00)
        R.save(layer[0:40] + "dnbr.tif")
        #####################################################################################
        dnbr = arcpy.Raster(r'' + layer)
        ID = str(layer)[0:39]
        bnd = file_workspace + years[i] + "\\brn_bounds\\" + ID + "_burn_bndy.shp"
        dnbr_msk = arcpy.gp.ExtractByMask_sa(dnbr, bnd)
        filled = arcpy.sa.Con(arcpy.sa.IsNull(dnbr_msk), arcpy.sa.FocalStatistics(dnbr_msk, arcpy.sa.NbrCircle(50, 'MAP'), 'MEAN'), dnbr_msk)
        filled = arcpy.sa.Con(arcpy.sa.IsNull(filled), arcpy.sa.FocalStatistics(dnbr_msk, arcpy.sa.NbrCircle(100, 'MAP'), 'MEAN'), filled)
        filled = arcpy.sa.Con(arcpy.sa.IsNull(filled), arcpy.sa.FocalStatistics(dnbr_msk, arcpy.sa.NbrCircle(150, 'MAP'), 'MEAN'), filled)
        filled = arcpy.sa.Con(arcpy.sa.IsNull(filled), arcpy.sa.FocalStatistics(dnbr_msk, arcpy.sa.NbrCircle(200, 'MAP'), 'MEAN'), filled)
        filled = arcpy.sa.Con(arcpy.sa.IsNull(filled), arcpy.sa.FocalStatistics(dnbr_msk, arcpy.sa.NbrCircle(250, 'MAP'), 'MEAN'), filled)
        filled = arcpy.sa.Con(arcpy.sa.IsNull(filled), arcpy.sa.FocalStatistics(dnbr_msk, arcpy.sa.NbrCircle(300, 'MAP'), 'MEAN'), filled)
        env.workspace = file_workspace + years[i] + "\\dNBR\\"
        arcpy.env.overwriteOutput = True # Set to overwrite pre-existing files
        filled = arcpy.gp.ExtractByMask_sa(filled, bnd, ID + "_dnbr.tif")
        del(M, R, dnbr, filled, ID)
    

    
