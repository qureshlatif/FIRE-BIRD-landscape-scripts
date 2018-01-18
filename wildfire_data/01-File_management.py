
#############################################################################################
## Run this script after extracting MTBS tar files                                         ##
## This script places RdNBR rasters and burn boundary shapefiles into appropriate folders. ##
#############################################################################################

import arcpy
import os
from arcpy import env
import numpy as np

file_workspace = "T:\\FS\\RD\\RMRS\\Science\\WTE\\Research\\Landscape_decision_tool\\wildfire_data\\"
env.workspace = file_workspace

years = list(range(1985, 2016))
years = map(str, years)

# Create folder structure if it doesn't already exist.
for yr in years:
    if not os.path.exists(file_workspace + yr):
        os.mkdir(file_workspace + yr)
    if not os.path.exists(file_workspace + yr + "\\brn_bounds"):
        os.mkdir(file_workspace + yr + "\\brn_bounds")
    if not os.path.exists(file_workspace + yr + "\\dNBR"):
        os.mkdir(file_workspace + yr + "\\dNBR")

env.workspace = file_workspace + "tars\\"
dNBR = arcpy.ListRasters("*_dnbr.tif", "TIF")
brn_bnds = arcpy.ListFeatureClasses("*burn_bndy.shp", "Polygon")

# Retrieve IDs for RdNBR and burn boundary layers
dNBR_string = map(str, dNBR)
dNBR_ID = list()
for ID in dNBR_string:
    dNBR_ID.append(ID[0:39])

bnd_string = map(str, brn_bnds)
bnd_ID = list()
for ID in bnd_string:
    bnd_ID.append(ID[0:39])

# Retain only files with both a dNBR layer and boundary shapefile for tranfer
ID = [val for val in bnd_ID if val in dNBR_ID]
st = set(ID)
keep = [i for i, val in enumerate(bnd_ID) if val in st]
brn_bnds = [brn_bnds[i] for i in keep]
keep = [i for i, val in enumerate(dNBR_ID) if val in st]
dNBR = [dNBR[i] for i in keep]
dNBR_string = map(str, dNBR)
bnd_string = map(str, brn_bnds)

# Retrieve year and month for each fire
dNBR_year = list()
dNBR_mo = list()
bnd_year = list()
bnd_mo = list()
for i in range(0, len(ID)):
    dNBR_year.append(dNBR_string[i][13:17])
    dNBR_mo.append(dNBR_string[i][17:19])
    bnd_year.append(bnd_string[i][13:17])
    bnd_mo.append(bnd_string[i][17:19])

# Add one to year for fires that occurred June - December
dNBR_year = np.array(map(int, dNBR_year))
dNBR_mo = np.array(map(int, dNBR_mo))
dNBR_year[dNBR_mo > 5] = dNBR_year[dNBR_mo > 5] + 1

bnd_year = np.array(map(int, bnd_year))
bnd_mo = np.array(map(int, bnd_mo))
bnd_year[bnd_mo > 5] = bnd_year[bnd_mo > 5] + 1

# Move files to year-specific folders
for i in range(0, len(bnd_year)):
    in_data = dNBR[i]
    out_data = file_workspace + str(dNBR_year[i]) + "\\dNBR\\" + dNBR[i]
    arcpy.Copy_management(in_data, out_data, "TIF")
    arcpy.Delete_management(in_data)

for i in range(0, len(bnd_year)):
    in_data = brn_bnds[i]
    out_data = file_workspace + str(bnd_year[i]) + "\\brn_bounds\\" + brn_bnds[i]
    arcpy.Copy_management(in_data, out_data, "Polygon")
    arcpy.Delete_management(in_data)

#############################################
# Move remaining files into desired folders #
#############################################

import os, shutil

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
