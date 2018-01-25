import arcpy
from arcpy import env
from arcpy.sa import *
arcpy.CheckOutExtension("Spatial")
import os, sys
import glob

#_________________Inputs__________________#
# Base directory #
base = "T:/FS/RD/RMRS/Science/WTE/Research/Landscape_decision_tool/"
wildfire_data_dir = "wildfire_data/"
wildfire_years = range(1985, 1994) + range(1995, 2016)
cc_data_dir = "cancov_data/GNN/"
cc_years = range(1984, 1993) + range(1994, 2013) + [2012]*2
COSASP = r'T:\FS\RD\RMRS\Science\WTE\Research\HSI_applic_tool\PA_RASTERS\cosasp'
HSI_dir = "R6/"
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

for t in range(0, len(wildfire_years)):
    wf_dir = base + wildfire_data_dir + str(wildfire_years[t]) + '/'
    bbnds = glob.glob(wf_dir + "brn_bounds/*.shp")
    bbnds = [bb for bb in bbnds if bb[89:91] == 'or' or bb[89:91] == 'wa']
    cc_raster = r'' + base + cc_data_dir + str(cc_years[t]) + '/cancov.tif'
    hsi_dir = base + HSI_dir + str(wildfire_years[t]) + '/'
    checkPath(hsi_dir)
    for f in range(0, len(bbnds)):
        fire = bbnds[f][89:110]
        dnbr = glob.glob(wf_dir + 'dNBR/' + fire + '*dnbr.tif')[0]
        sr = arcpy.Describe(dnbr).spatialReference
        env.outputCoordinateSystem = sr
        fire_dir = hsi_dir + fire + '/'
        checkPath(fire_dir)
        INPUTS_dir = fire_dir + 'INPUTS/'
        checkPath(INPUTS_dir)
        bb = bbnds[f]
        bb_buff = arcpy.Buffer_analysis(bb, INPUTS_dir + "fire_mask.shp", "1000", "", "", "ALL")
        LocdNBR = FocalStatistics(dnbr, NbrRectangle(3, 3, "CELL"), "MEAN", "DATA")
        LocdNBR_mask = ExtractByMask(LocdNBR, INPUTS_dir + "fire_mask.shp")
        LocdNBR_mask.save(INPUTS_dir + "locdnbr.tif")
        env.extent = LocdNBR_mask.extent
        env.snapRaster = LocdNBR_mask
        out = ExtractByMask(COSASP, LocdNBR_mask)
        out.save(INPUTS_dir + "cosasp.tif")
        cc = ExtractByMask(cc_raster, LocdNBR_mask)
        cc40 = Float(Con(cc >= 40, 1, 0))
        loccc = FocalStatistics(cc40, NbrRectangle(3, 3, "CELL"), "MEAN", "DATA")
        out =  ExtractByMask(loccc, LocdNBR_mask)
        out.save(INPUTS_dir + "loccc.tif")
        landcc = FocalStatistics(cc40, NbrCircle(1000, "MAP"), "MEAN", "DATA")
        out = ExtractByMask(landcc, LocdNBR_mask)
        out.save(INPUTS_dir + "landcc.tif")
        arcpy.ResetEnvironments()
        
