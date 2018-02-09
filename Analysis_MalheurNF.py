#-------------------------------------------------------------------------------
# Name:        FIRE-BIRD Landscape analysis for Malheur NF
# Author:      qlatif
#-------------------------------------------------------------------------------

import os, arcpy
from arcpy.sa import *
import numpy as np

arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput = True

# function(s)
def checkPath(path):
    try: os.makedirs(path)
    except: pass

def checkExtentOverlap(rast, SA):
    check = True
    rast_TOP = float(arcpy.GetRasterProperties_management(rast, "TOP").getOutput(0))
    rast_BOTTOM = float(arcpy.GetRasterProperties_management(rast, "BOTTOM").getOutput(0))
    rast_LEFT = float(arcpy.GetRasterProperties_management(rast, "LEFT").getOutput(0))
    rast_RIGHT = float(arcpy.GetRasterProperties_management(rast, "RIGHT").getOutput(0))
    SA_TOP = arcpy.Describe(SA).extent.YMax
    SA_BOTTOM = arcpy.Describe(SA).extent.YMin
    SA_LEFT = arcpy.Describe(SA).extent.XMin
    SA_RIGHT = arcpy.Describe(SA).extent.XMax
    if (rast_TOP <= SA_BOTTOM) | (rast_BOTTOM >= SA_TOP) | (rast_LEFT >= SA_RIGHT) | (rast_RIGHT <= SA_LEFT):
        check = False
    return check

#_____________________Inputs____________________________#
base_dir = 'T:/FS/RD/RMRS/Science/WTE/Research/Landscape_decision_tool/'
NF_bound = u'' + base_dir + 'Malheur/NF_bound_NAD83alb.shp'
buffer_dist = 90000 # In meters, buffer distance to apply to NF bound
DCF_mask = r'T:/FS/RD/RMRS/Science/WTE/Research/HSI_applic_tool/TOOLBOX/masks/BBWO/dcfmask.img'
apply_mask = True # Set to False if you don't want DCF mask applied prior to tabulating habitat values by year.
HSI_dir = base_dir + 'R6/' # Data source - repository of BBWO HSI rasters
years = range(1985, 1994) + range(1995, 2016)

# Potential nest densities by suitability class (per acre)
D = {"low" :
     {"mean" : 1.28, "95lo" : 0.41, "95hi" : 2.26},
     "mod" :
     {"mean" : 5.92, "95lo" : 3.41, "95hi" : 8.84},
     "high" :
     {"mean" : 11.03, "95lo" : 8.38, "95hi" : 13.93}
     }
#**Note: currently set to values representing Toolbox, OR**

HSI_thresholds = {"low-mod" : 2.5, "mod-high" : 6.5} # Thresholds for classifying low, moderate, and high suitability habitat

#_______________________________________________________#

arcpy.env.outputCoordinateSystem = arcpy.Describe(DCF_mask).spatialReference

#Set workspace and scratch workspace for temp outputs/calculations
arcpy.env.workspace = base_dir + "_workspace"
checkPath(arcpy.env.workspace)
arcpy.env.scratchWorkspace = base_dir + "_scratch"
checkPath(arcpy.env.scratchWorkspace)

SA = arcpy.Buffer_analysis(NF_bound, buffer_distance_or_field = str(buffer_dist) + " meters", dissolve_option = "ALL")
dcf_msk = ExtractByMask(DCF_mask, SA)

YEAR = np.array(range(1985, 2016))
HAB_GEN_LOW = np.zeros(shape = len(YEAR))
HAB_GEN_MOD = np.zeros(shape = len(YEAR))
HAB_GEN_HIGH = np.zeros(shape = len(YEAR))

for t in range(0, len(years)):
    yr = years[t]
    hsiR6 = Raster(HSI_dir + str(yr) + '/HSI.tif')
    lo_count = 0; md_count = 0; hi_count = 0
    if checkExtentOverlap(hsiR6, SA):
        hsi = ExtractByMask(hsiR6, SA)
        if apply_mask:
            hsi = ExtractByMask(hsi, dcf_msk)
        AllNoDataTrue = int(arcpy.GetRasterProperties_management(hsi, "ALLNODATA").getOutput(0))
        if AllNoDataTrue == 0:
            rows = arcpy.SearchCursor(hsi, "", "", "VALUE; COUNT")
            for r in rows:
                v = r.getValue("VALUE")
                c = r.getValue("COUNT")
                if v == 0 or v == 1 or v == 2:
                    lo_count = lo_count + c
                if v == 3 or v == 4 or v == 5:
                    md_count = md_count + c
                if v == 6 or v == 7 or v == 8:
                    hi_count = hi_count + c
    HAB_GEN_LOW[YEAR == yr] = lo_count * (30 * 30) * 0.000247105
    HAB_GEN_MOD[YEAR == yr] = md_count * (30 * 30) * 0.000247105
    HAB_GEN_HIGH[YEAR == yr] = hi_count * (30 * 30) * 0.000247105

HAB_LOW = np.zeros(shape = len(YEAR))
HAB_MOD = np.zeros(shape = len(YEAR))
HAB_HIGH = np.zeros(shape = len(YEAR))
PotNest_mn = np.zeros(shape = len(YEAR))
PotNest_95lo = np.zeros(shape = len(YEAR))
PotNest_95hi = np.zeros(shape = len(YEAR))

for t in range(4, len(YEAR)):
    HAB_LOW[t] = sum(HAB_GEN_LOW[(t-4):t])
    HAB_MOD[t] = sum(HAB_GEN_MOD[(t-4):t])
    HAB_HIGH[t] = sum(HAB_GEN_HIGH[(t-4):t])
    PotNest_mn[t] = (HAB_LOW[t]/1000)*D['low']['mean'] + (HAB_MOD[t]/1000)*D['mod']['mean'] + (HAB_HIGH[t]/1000)*D['high']['mean']
    PotNest_95lo[t] = (HAB_LOW[t]/1000)*D['low']['95lo'] + (HAB_MOD[t]/1000)*D['mod']['95lo'] + (HAB_HIGH[t]/1000)*D['high']['95lo']
    PotNest_95hi[t] = (HAB_LOW[t]/1000)*D['low']['95hi'] + (HAB_MOD[t]/1000)*D['mod']['95hi'] + (HAB_HIGH[t]/1000)*D['high']['95hi']

#### Test ####
#x = np.array([1, 2])
#y = np.array([3, 4])
#output = np.vstack([np.array([x, y]).T])
#np.savetxt(r'' + base_dir + 'test.csv', output, delimiter = ',', header = 'X,Y', comments='', fmt = '%.2f,%d')
##############

output = np.hstack([np.array([YEAR, HAB_GEN_LOW, HAB_GEN_MOD, HAB_GEN_HIGH,
                              HAB_LOW, HAB_MOD, HAB_HIGH,
                              PotNest_mn, PotNest_95lo, PotNest_95hi]).T])
if apply_mask:
    out_file = r'' + base_dir + 'HabByYear_DCFonly.csv'
else:
    out_file = r'' + base_dir + 'HabByYear_nomask.csv'
np.savetxt(out_file, output, delimiter = ',',
           header = 'YEAR,HAB_GEN_LOW,HAB_GEN_MOD,HAB_GEN_HIGH,HAB_LOW,HAB_MOD,HAB_HIGH,PNest_mn,PNest_95lo,PNest_95hi',
           comments='', fmt = '%f')
