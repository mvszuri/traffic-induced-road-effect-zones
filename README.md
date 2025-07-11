# traffic-induced-road-effect-zones
Python code to calculate traffic-induced Road Effect Zones (REZs)
The code was used to generate the results presented in the manuscript titled "Global expansion of the ecological impact of extra-urban road traffic" (version 11 July 2025).
The analysis makes used of the traffic volume predictions (Annual Average Daily Traffic [AADT]) documented in: https://doi.org/10.1038/s41597-024-03287-z
## Details of the individual files
#### Kernel_density_function_20250117.py
This code was used to calculate the REZs with the kernel density function from ArcGIS Pro 2.9.5 (Environmental Systems Research Institute, Redlands, USA).
At the beginning of the code, the following settings need to be specified:
runName: Choose a unique name for the run of the code
road_shp: Specify the location of the shapefile containing the traffic volume predictions. These can be downloaded from the ETH research collection: https://doi.org/10.3929/ethz-b-000666313
populationField: Specify the attribute in road_shp containing the traffic volume predictions ("median", "q__0_05" or "q__0_95")
cellSize: Specify the cell size in meters of the output kernel density maps. 
searchRadius: Specify the maximum search radius in the kernel density function.
baseDir: Specify the directory where the output should be stored.
AADT_threshs: Specify the different AADT thresholds with which the kernel density values should be truncated (e.g. 500, 5000 and 10000).
LU_rast: Specify the ESA Climate Change Initiative land cover rasters. These can be downloaded from: https://maps.elie.ucl.ac.be/CCI/viewer/download.php
KBA_shp: Specify the shapefile with the Key Biodiversity Areas. This file can be requested from: https://www.keybiodiversityareas.org/kba-data/request. The original file was dissolved and then intersected with the countries in country_shp.
Fishnet_size: To speed-up calculations the analyses are performed in different regions. Here the size of these regions in meters should be specified. Sensitivity analysis showed that 500000 m was the fastest processing time.
country_shp: Specify a shapefile with the countries in which the analysis should be performed. For all countries in the world, a shapefile can be downloaded from: https://www.geoboundaries.org/

#### Summarise_in_Rasters_20250211.py
This code was used to summarise the high-resolution REZs calculated with Kernel_density_function_20250117.py in a coarser global raster (e.g. 5000 x 5000 m raster).
At the beginning of the code, the following settings need to be specified:
runNames: Specify the runName from the code above.
baseDir: Specify the directory where the output should be stored.
AADT_threshs: Specify the different AADT thresholds that were used for this runName in Kernel_density_function_20250117.py
country_shp: Specify a shapefile with the countries in which the analysis should be performed. For all countries in the world, a shapefile can be downloaded from: https://www.geoboundaries.org/
BP_cellsize: Specify the raster cell size in meters for the aggregated global raster.
