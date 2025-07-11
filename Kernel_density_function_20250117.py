# Import system modules
import arcpy, os, math, time
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
runName = "Run_2000_median_20250212"
road_shp = r"C:\PROCESSING\Traffic_network\Traffic_analysis_2000\GRIP4_ExSet_2000_AADTpred_20240312.shp"
populationField = "median" #"median", "q__0_05", "q__0_95"
cellSize = 20
searchRadius = 300
baseDir = r"C:\PROCESSING\Kernel_density"
# Set the AADT thresholds
AADT_threshs = [500,5000,10000]
# Land use rasters from the ESA Climate Change Initiative
LU_rast = r"C:\PROCESSING\Kernel_density\input_data.gdb\ESACCI_2000_tif_MW" # "ESACCI_1992_tif_MW", "ESACCI_2000_tif_MW", "ESACCI_2015_tif_MW"
# The feature with the Key Biodiversity Areas (KBAs). The original file was dissolved and then intersected with the countries in country_shp.
KBA_shp = [r"C:\PROCESSING\Kernel_density\input_data.gdb","KBAsGlobal_2022_March_01_POL_MW_DissAll_pCountry"]
# Width and length of a Fishnet in meters
Fishnet_size = 500000 # Sensitivity analysis showed that 500000 was the fastest processing time.
# Shapefile with all the countries in the world.
country_shp = "D:\\Publications\\Pub_2024_Global_traffic_flow\\Statistics\\geoBoundariesCGAZ_ADM0_MW.shp"

# Create a new file geodatabase in the new baseDir
if not os.path.exists(baseDir+"\\"+runName+".gdb"):
    arcpy.CreateFileGDB_management(baseDir, runName+".gdb")
gdbDir = baseDir+"\\"+runName+".gdb"

# Arcpy setting
arcpy.SetProduct('ArcInfo')
arcpy.CheckOutExtension('Spatial')
arcpy.CheckOutExtension('GeoStats')
arcpy.env.overwriteOutput = True
env.workspace = gdbDir

# Generate a table with all the country names
ISOs = []
countries = []
countr_area = []
with arcpy.da.SearchCursor(country_shp, ['name','ISO_CODE','Area']) as cursor:
    for row in cursor:
        countries.append(row[0])
        ISOs.append(row[1])
        countr_area.append(row[2])
country_table = pd.DataFrame({'Country': countries, 'ISO_CODE': ISOs, 'CountrArea': countr_area})
# There is one country without a name in the list. Remove this one.
country_table = country_table[country_table['Country'] != '']

# Read the feature with the KBAs per country
KBA_pCountry = gpd.read_file(KBA_shp[0], layer=KBA_shp[1])
 
# Create several lists
ISO_list = [] # List with all the ISO codes of the countries
Country_list = [] # A list specifying all the countries
CntryHa_list = [] # A list with the surface area of each country (in ha)
KDE_rasts = [] # A list to put in the density rasters
Area_list = [] # A list with all the areas of the different KDE rasters
Area_KBA_list = [] # A list with all the areas of the different KDE rasters within Key Biodiversity Areas
AADT_th_list = [] # A list with all the AADT_thresholds
LUclass_counts = {} # A dictionairy of counts of classes in LU_rast
    
# Loop over the different countries
for index, row in country_table.iterrows():
    country = row['Country']
    ISO = row['ISO_CODE']
    Countr_Area = row['CountrArea']
    print(datetime.now().strftime("%d.%m %H:%M:%S"), ": PROCESSING: ", country, sep='')  
    # Clip the roads from a certain country.
    arcpy.MakeFeatureLayer_management(country_shp, "country_lyr", "ISO_CODE = '{}'".format(ISO)) # Make feature layer for input feature class
    arcpy.analysis.Clip(road_shp, "country_lyr", "in_memory\\road_country")
    
    # Select the relevant country from the KBA feature
    KBA_geo = KBA_pCountry[KBA_pCountry['ISO_CODE'] == ISO]
    
    # Get the extent of the input feature class
    desc = arcpy.Describe("in_memory\\road_country")
    xmin = desc.extent.XMin
    ymin = desc.extent.YMin
    xmax = desc.extent.XMax
    ymax = desc.extent.YMax
    spatial_ref = desc.spatialReference
    
    # Some countries will not have any roads. If this is the case the current iteration should be skipped.
    if math.isnan(xmin) or math.isnan(ymin) or math.isnan(ymax):
        print("Skipping due to no roads")
        continue
    
    # Create a fishnet feature class based on the extent of the input shapefile
    origin_coord = f"{xmin} {ymin}"
    y_axis_coord = f"{xmin} {ymax}"
    arcpy.CreateFishnet_management("Fishnet_temp", origin_coord, y_axis_coord, str(Fishnet_size), str(Fishnet_size), "", "", None, "NO_LABELS", "in_memory\\road_country", "POLYGON")
    # From the fishnet, select only those polygons that intersect with the inFeature
    arcpy.MakeFeatureLayer_management("Fishnet_temp", "Fishnet_lyr") # Make feature layer for input feature class
    arcpy.SelectLayerByLocation_management("Fishnet_lyr", "INTERSECT", "in_memory\\road_country") # Select features that intersect with intersect_feature
    arcpy.CopyFeatures_management("Fishnet_lyr", "Fishnet") # Copy selected features to output feature class
    arcpy.management.Delete("Fishnet_temp")
    
    # Buffer the fishnet with 300 m
    arcpy.analysis.PairwiseBuffer("Fishnet", "Fishnet_buf", "300 Meters", "", "", "", "")
    Fishnet_number = arcpy.management.GetCount("Fishnet_buf")[0]
    
    # Clip input features to each fishnet tile and run Kernel Density
    with arcpy.da.SearchCursor("Fishnet_buf", ["OID@", "SHAPE@"]) as cursor:
        for row in cursor:
            oid = row[0] # OID (OBJECTID) value, which is the same between the fishnet and buffered fishnet.
            print(datetime.now().strftime("%d.%m %H:%M:%S"), ": Tile: ", str(oid), "/", Fishnet_number, sep='')
            # Clip to the buffered fishnet
            clip = arcpy.Clip_analysis("in_memory\\road_country", row[1], "in_memory\\tile{}".format(oid))
            # Calculate the kernel density (KD).
            kd = KernelDensity(clip, populationField, cellSize, searchRadius, "HECTARES", "DENSITIES", "PLANAR")
            # Loop through the AADT thresholds
            for AADT_thresh in AADT_threshs:
                # Transform the AADT_thresh values to kd values.
                # This equation was derived from a sample of roads in South Africa
                # See D:\Publications\Prep_Road_effect_zone\Density_vs_AADT_sampling.xlsx
                KD_thresh = np.divide(AADT_thresh,2.9433201904)
                # Define the name for the raster
                kd_name = "KD_"+ISO+"_"+str(oid)+"_th"+str(AADT_thresh)
                # Select only those raster values from the KD raster that are above KD_thresh
                kd_bin = arcpy.sa.Con(kd > KD_thresh, 1)
                # Clip the KD with the original fishnet
                # Construct a SQL query to select the polygon based on the OID
                sql_query = "OBJECTID = {}".format(oid)
                arcpy.MakeFeatureLayer_management("Fishnet", "FN_sel", sql_query)
                arcpy.Clip_management(kd_bin, "#", kd_name, "FN_sel", "#", "ClippingGeometry")
                # Calculate the total area of cells above the KD_thresh
                # Read the raster attribute table where values are 1
                RAT_kd = arcpy.da.TableToNumPyArray(kd_name,"Count","Value = 1")
                # Set the two variables to 0. They will be updated under certain conditions in the following if-statements.
                area = 0  
                KBA_kd_area = 0
                # Check whether any cells are in the binary selection or not
                if RAT_kd.shape[0] > 0:
                    # Calculate the area of land (ha) that is covered by the road effect zone
                    area = (RAT_kd[0].astype('int64')*cellSize*cellSize)/10000 # Number of hectares
                    # Copy the raster to a tiff file
                    rast_path = r"C:\PROCESSING\Kernel_density\scratch2\temp_rast.tif"
                    arcpy.management.CopyRaster(kd_name, rast_path)
                    # Determine the total area of the road effect zone that fall within a Key Biodiversity Area (KBA)
                    # Calculate the area of the road effect zone that falls in a KBA
                    if KBA_geo.shape[0] > 0:
                        with rasterio.open(rast_path) as src:
                          KBA_array, out_transform = rasterio.mask.mask(src, KBA_geo.geometry, filled = True)      
                        # Reset the no-data values in the KBA raster to 0
                        KBA_array[KBA_array == 255] = 0
                        # Calculate the total ha of the road effect zone that fall within a KBA
                        KBA_kd_area = (np.sum(KBA_array)*cellSize*cellSize)/10000
                            
                    # Delete the temporary raster
                    arcpy.management.Delete(rast_path)
                        
                    # From the road-effect zone, calculate the frequency of land-use classes in LU_rast
                    out = ExtractByMask(LU_rast, kd_name)
                    out.save(kd_name+"_LU")
                    # Check whether the land-use raster is empty.
                    if ((out.maximum is not None) and (out.maximum > 0)):
                        # Add the counts to the dictionairy class_counts
                        LUclass_counts[kd_name] = {}
                        for key, val in arcpy.da.TableToNumPyArray(kd_name+"_LU",["Value","Count"]):
                            LUclass_counts[kd_name][str(key)] = val
                        
                
                # Fill in the lists
                AADT_th_list.append(AADT_thresh)
                Area_list.append(area)
                Area_KBA_list.append(KBA_kd_area)
                Country_list.append(country)
                ISO_list.append(ISO)
                CntryHa_list.append(Countr_Area)
                # Write the rasters to a list
                KDE_rasts.append(kd_name)


# Calculate the total area of a certain land-use that the road effect zone covers
# Determine the cell-size of the LU raster
LU_cellsize = float(arcpy.management.GetRasterProperties(LU_rast, "CELLSIZEY")[0])
# Transform the dict into a data frame and calculate the LU coverage per hectare
LU_tab = (pd.DataFrame.from_dict(LUclass_counts,'index')*LU_cellsize*LU_cellsize)/10000 # Number of hectares
LU_tab['Scenario'] = LU_tab.index.values
# Create a final table with the output
output_tab = pd.DataFrame(list(zip(Country_list,ISO_list,CntryHa_list,KDE_rasts,Area_list,Area_KBA_list,AADT_th_list)),columns =['Country','ISO_CODE','CountrArea_ha','Scenario','Area_ha','Area_KBA_ha','AADT_thresh'])
# Merge the output_tab with the LU_tab
output_tab = output_tab.merge(LU_tab, on='Scenario', how='left')
# Save the output tab
output_tab.to_excel(baseDir+"\\"+runName+".xlsx")
