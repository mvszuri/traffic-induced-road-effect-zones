# Import system modules
import arcpy, os, math, time, re
from arcpy import env
from arcpy.sa import *
from datetime import datetime
import pandas as pd
import numpy as np
import geopandas as gpd
import rasterio
from rasterio import mask
from rasterio.plot import show

# Set local variables
runNames = ["Run_2000_median_20250212","Run_2015_median_20250212"] #["Run_1975_median_20250212","Run_1990_median_20250212","Run_2000_median_20250212","Run_2015_median_20250212"]
baseDir = r"C:\PROCESSING\Kernel_density"
# Set the AADT thresholds
AADT_threshs = [10000] # [500,5000], [10000]
# Shapefile with all the countries in the world.
country_shp = "D:\\Publications\\Pub_2024_Global_traffic_flow\\Statistics\\geoBoundariesCGAZ_ADM0_MW.shp"
# Blueprint cell-size in meters
BP_cellsize = 5000 

# Arcpy setting
arcpy.SetProduct('ArcInfo')
arcpy.CheckOutExtension('Spatial')
arcpy.CheckOutExtension('GeoStats')
arcpy.env.overwriteOutput = True

# Create the global raster that will serve as blue-print for the following rasters
# Take the extent and coord_syst of the country_shp, which was used for creating the tiles.
desc = arcpy.Describe(country_shp)
xmin = desc.extent.XMin
ymin = desc.extent.YMin
xmax = desc.extent.XMax
ymax = desc.extent.YMax
# Set the output coordinate system to the same as country_shp
spatRef = desc.spatialReference
arcpy.env.outputCoordinateSystem = spatRef
# Create an empty raster with the specified coordinate system and resolution
emptyRast = CreateConstantRaster(0, "INTEGER", BP_cellsize, arcpy.Extent(xmin, ymin, xmax, ymax))
# Save the raster and set it as the snap_raster
emptyRast.save("C:\\PROCESSING\\Kernel_density\\Input_data.gdb\\EmptyRast")
arcpy.env.snapRaster = "C:\\PROCESSING\\Kernel_density\\Input_data.gdb\\EmptyRast"

# Loop over the runNames
for runName in runNames:
    print(runName)
    # Set the file geodatabase in baseDir
    gdbDir = baseDir+"\\"+runName+".gdb"
    
    # Create a new file geodatabase where the outputs will be saved
    if not os.path.exists(baseDir+"\\"+runName+"_globRast.gdb"):
        arcpy.CreateFileGDB_management(baseDir, runName+"_globRast.gdb")
    gdbOut = baseDir+"\\"+runName+"_globRast.gdb"
    #Set the workspace
    env.workspace = gdbDir
    
    # List all raster datasets
    tiles = arcpy.ListRasters()
    # Loop over the different thresholds
    for th in AADT_threshs:
        # Filter the list using the pattern
        pattern = re.compile(fr'{th}$')
        filtered_tiles = [s for s in tiles if re.search(pattern, s)]
        for tile in filtered_tiles:
            tileAggr = Aggregate(tile, int(BP_cellsize/20), "SUM", "EXPAND", "DATA")
            #Check that the aggregate raster is empty. If so, do not save.
            if ((tileAggr.maximum is not None) and (tileAggr.maximum > 0)): 
                tileAggr.save(gdbOut+"\\"+tile+"_aggr")
                print(tile)
        # List all raster datasets
        env.workspace = gdbOut
        tiles_aggr = arcpy.ListRasters()
        # Filter the list using the pattern
        pattern = re.compile(fr'{th}_aggr$')
        filtered_tiles_aggr = [s for s in tiles_aggr if re.search(pattern, s)]
        # Stitch all the new raster together and "SUM" overlapping rasters
        outRastName = runName+"_th_"+str(th)+"_glob"
        arcpy.management.MosaicToNewRaster(filtered_tiles_aggr, gdbOut, outRastName, spatRef , "32_BIT_UNSIGNED", "", 1, "SUM", "")
        PercRast = Raster(outRastName)/(BP_cellsize/20)**2
        # Save the output raster
        PercRast.save(outRastName+"_perc")
        # Reset the workspace back to the original workspace
        env.workspace = gdbDir
        




