
#############################################################################################
## Run this script after extracting MTBS tar files                                         ##
## This script places RdNBR rasters and burn boundary shapefiles into appropriate folders. ##
#############################################################################################

import arcpy
import os
from arcpy import env
import numpy as np
import shutil

file_workspace = "T:\\FS\\RD\\RMRS\\Science\\WTE\\Research\\Landscape_decision_tool\\wildfire_data\\"
env.workspace = file_workspace

years = list(range(1985, 2016))
years = map(str, years)

env.workspace = file_workspace + "tars\\"

#############################################
# Move remaining files into desired folders #
#############################################

## Fires missing dNBR (only moves .shp files; if running again, rescript to move all files)
dNBR_missing = arcpy.ListFeatureClasses("*burn_bndy.shp", "Polygon")
ID_missing = map(str, dNBR_missing)

sourcepath='T:\\FS\\RD\\RMRS\\Science\\WTE\\Research\\Landscape_decision_tool\\wildfire_data\\tars'
source = os.listdir(sourcepath)
destinationpath = 'T:\\FS\\RD\\RMRS\\Science\\WTE\\Research\\Landscape_decision_tool\\wildfire_data\\miss_dNBR'
for files in source:
    for i in range(0, len(ID_missing)):
        if files.startswith(ID_missing[i]):
            shutil.move(os.path.join(sourcepath, files), os.path.join(destinationpath, files))

## Remaining files ##
source = os.listdir(sourcepath)
destinationpath = 'T:\\FS\\RD\\RMRS\\Science\\WTE\\Research\\Landscape_decision_tool\\wildfire_data\\extra_files'
for files in source:
    if files.endswith('.aux'):
        shutil.move(os.path.join(sourcepath, files), os.path.join(destinationpath, files))
    if files.endswith('.rrd'):
        shutil.move(os.path.join(sourcepath, files), os.path.join(destinationpath, files))
    if files.endswith('.tif'):
        shutil.move(os.path.join(sourcepath, files), os.path.join(destinationpath, files))
    if files.endswith('.dbf'):
        shutil.move(os.path.join(sourcepath, files), os.path.join(destinationpath, files))
    if files.endswith('.prj'):
        shutil.move(os.path.join(sourcepath, files), os.path.join(destinationpath, files))
    if files.endswith('.sbn'):
        shutil.move(os.path.join(sourcepath, files), os.path.join(destinationpath, files))
    if files.endswith('.shp'):
        shutil.move(os.path.join(sourcepath, files), os.path.join(destinationpath, files))
    if files.endswith('.shx'):
        shutil.move(os.path.join(sourcepath, files), os.path.join(destinationpath, files))
    if files.endswith('.sbx'):
        shutil.move(os.path.join(sourcepath, files), os.path.join(destinationpath, files))
    if files.endswith('.txt'):
        shutil.move(os.path.join(sourcepath, files), os.path.join(destinationpath, files))
    if files.endswith('.pdf'):
        shutil.move(os.path.join(sourcepath, files), os.path.join(destinationpath, files))
    if files.endswith('.kmz'):
        shutil.move(os.path.join(sourcepath, files), os.path.join(destinationpath, files))
