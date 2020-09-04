# Schrimpf_Lynch_fetch
Extracted data and code for the Schrimpf and Lynch manuscript on fetch and pursuit-divers.  
Title: The role of wind fetch in structuring Antarctic seabird breeding occupancy   

Data on penguins came from MAPPPD (http://www.penguinmap.com/), but the code for turning those files into a summary file is included in the Extract_MAPPPD R script.

The data file 'pop.pursuit.csv' incorporates the data from MAPPPD as well as the Antarctic Shag population data from Schrimpf et al. 2018 (see https://github.com/lynch-lab/Antarctic_shags). Those sources were combined manually, resulting in the csv file included.

The rest of the analysis and figures for the paper are produced with the code found in the R script: Analysis_code.R.

Wind data files must be downloaded from the Copernicus website: https://cds.climate.copernicus.eu/about-c3s  
The link I used was this: https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels-monthly-means?tab=form.  
The selection used in the manuscript included the following options:  
Product type:  
--Monthly average reanalysis  
Variable (under “Popular”):  
--10m u-component of wind  
--10m v-component of wind  
Year:  
--1979-2019  
Month:  
--Jan, Feb, Nov, Dec  
Time:  
--00:00  
Format:  
--GRIB  

Other summary files produced by the R script are included in the project, and any shapefiles are stored in the Shapefile folder. Note, the polygon shapefile used in the fetch calculation is too large to put on GitHub, but is available online from the Antarctic Digital Database:  
https://www.add.scar.org/  
The coastline I used is the medium-resolution polygon coastline layer, and the most recent version is available here:  
https://data.bas.ac.uk/items/862f7159-9e0d-46e2-9684-df1bf924dabc/  
