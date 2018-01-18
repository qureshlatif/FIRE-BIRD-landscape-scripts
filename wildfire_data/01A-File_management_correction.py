
#############################################################################################
## 1. Corrects originally implemented file organization so that June fires are counted as  ##
## contributing to habitat in the subsequent year (i.e., add 1 to year for all firest from ##
## June to December).                                                                      ##
## 2. Moves all June cmort.tif files (only run if cmort already generated).                ##
#############################################################################################

import arcpy
import os
from arcpy import env
import numpy as np

file_workspace = "T:\\FS\\RD\\RMRS\\Science\\WTE\\Research\\Landscape_decision_tool\\wildfire_data\\"

years = list(reversed(range(1985, 2015)))

for i in range(0, len(years)):
    env.workspace = file_workspace + str(years[i]) + "\\RdNBR\\"
    RdNBR = arcpy.ListRasters("*rdnbr.tif", "TIF")
    RdNBR_string = map(str, RdNBR)
    for j in range(0, len(RdNBR_string)):
        if int(RdNBR_string[j][17:19]) == 6:
            in_data = RdNBR[j]
            out_data = file_workspace + str(years[i] + 1) + "\\RdNBR\\" + RdNBR[j]
            arcpy.Copy_management(in_data, out_data, "TIF")
            arcpy.Delete_management(in_data)
    env.workspace = file_workspace + str(years[i]) + "\\brn_bounds\\"
    brn_bnds = arcpy.ListFeatureClasses("*burn_bndy.shp", "Polygon")
    bnd_string = map(str, brn_bnds)
    for j in range(0, len(bnd_string)):
        if int(bnd_string[j][17:19]) == 6:
            in_data = brn_bnds[j]
            out_data = file_workspace + str(years[i] + 1) + "\\brn_bounds\\" + brn_bnds[j]
            arcpy.Copy_management(in_data, out_data, "Polygon")
            arcpy.Delete_management(in_data)

for i in range(0, len(years)):
    env.workspace = file_workspace + str(years[i]) + "\\RdNBR\\"
    ccmort = arcpy.ListRasters("*cmort.tif", "TIF")
    ccmort_string = map(str, ccmort)
    for j in range(0, len(ccmort_string)):
        if int(ccmort_string[j][17:19]) == 6:
            in_data = ccmort[j]
            out_data = file_workspace + str(years[i] + 1) + "\\RdNBR\\" + ccmort[j]
            arcpy.Copy_management(in_data, out_data, "TIF")
            arcpy.Delete_management(in_data)
