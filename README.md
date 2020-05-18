# Schrimpf_Lynch_fetch
Extracted data and code for the Schrimpf and Lynch manuscript on fetch and pursuit-divers.
Title: The role of wind fetch in structuring Antarctic seabird breeding occupancy 

Data on penguins came from MAPPPD (http://www.penguinmap.com/), but the code for turning those files into a summary file is included in the Extract_MAPPPD R script.

The data file 'pop.pursuit.csv' incorporates the data from MAPPPD as well as the Antarctic Shag population data from Schrimpf et al. 2018 (see https://github.com/lynch-lab/Antarctic_shags). Those sources were combined manually, resulting in the csv file included.

The rest of the analysis and figures for the paper are produced with the code found in the R script: Analysis_code.R.

Wind data files extracted from https://cds.climate.copernicus.eu/about-c3s are included as .grib files.

Other summary files produced by the R script are included in the project, and any shapefiles are stored in the Shapefile folder.
