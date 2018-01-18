################################################################################
# 1. Stitches together burn boundaries and RdNBR layers by year (1985 - 2015). #
# 2. Masks mosaic RdNBR layer by burn boundaries in each year.                 #
################################################################################

import arcpy, os
from arcpy import env

file_workspace = "T:\\FS\\RD\\RMRS\\Science\\WTE\\Research\\Landscape_decision_tool\\wildfire_data\\"
env.workspace = file_workspace

years = list(range(1985, 1991)) 
years = map(str, years)

## Mosaic RdNBR and boundary layers, then mask RdNBR by fire boundaries for each year.
arcpy.env.overwriteOutput = True # Set to overwrite pre-existing files
for yr in years:
    env.workspace = file_workspace + yr + "\\brn_bounds\\"
    brn_bnds = arcpy.ListFeatureClasses("*burn_bndy.shp", "Polygon")
    arcpy.Union_analysis(in_features = brn_bnds,
                         out_feature_class = file_workspace + yr + "\\brn_boundaries.shp",
                         join_attributes="ONLY_FID", cluster_tolerance="", gaps="GAPS")

    env.workspace = file_workspace + yr + "\\dNBR\\"
    dNBR = arcpy.ListRasters("*dnbr.tif", "TIF")
    arcpy.MosaicToNewRaster_management(input_rasters = dNBR,
                                       output_location = file_workspace + yr + "\\dNBR",
                                       raster_dataset_name_with_extension="dNBR_mosaic.tif",
                                       coordinate_system_for_the_raster="", pixel_type="16_BIT_SIGNED",
                                       cellsize="", number_of_bands="1", mosaic_method="LAST", mosaic_colormap_mode="FIRST")
arcpy.CheckOutExtension("spatial")
for yr in years:
    arcpy.gp.ExtractByMask_sa(file_workspace + yr + "\\dNBR\\dNBR_mosaic.tif",
                              file_workspace + yr + "\\brn_boundaries.shp",
                              file_workspace + yr + "\\dNBR.tif")
