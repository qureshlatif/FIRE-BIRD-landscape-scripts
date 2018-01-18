##########################################################################
# Applies Lookup to extract cancov.tif layers from GNN time series data. #
##########################################################################

import arcpy, os
from arcpy import env

file_workspace = "T:\\FS\\RD\\RMRS\\Science\\WTE\\Research\\Landscape_decision_tool\\cancov_data\\GNN\\"
env.workspace = file_workspace
GNN_source = "T:/FS/NFS/Deschutes/Project/soopsDataMgmt2009/Silviculture/SO/Workspace/mlsimpson/GNN_TimeSeries/masked/"

years = list(range(1988, 2013))
years = map(str, years)

arcpy.env.overwriteOutput = True # Set to overwrite pre-existing files
arcpy.CheckOutExtension("Spatial")
for i in range(0, len(years)):
    env.workspace = file_workspace + years[i] + "\\"
    arcpy.env.overwriteOutput = True # Set to overwrite pre-existing files
    arcpy.gp.Lookup_sa(GNN_source + "mr200_" + years[i], "CANCOV", years[i] + "\\cancov.tif")
