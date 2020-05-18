# Code to accompany the Schrimpf and Lynch manuscript: The role of wind fetch in structuring Antarctic seabird breeding occupancy 





# Libraries ---------------------------------------------------------------

# These are the packages required:

library(fetchR) # for calculating fetch

library(lessR) # for the to() function

library(MHTmult) # for the Sidak corrections

# Mapping and working in GIS:
library(rgdal) # reads in and manages shape files (requires GDAL)
library(sp) # more spatial analysis functions
library(raster) # more spatial analysis functions
library(rgeos) # for the gDifference() function

# Plotting:
library(ggplot2) # for actual plotting and for the fortify() function
library(scico) # Scientific color palettes



# NOT SURE IF THESE ARE NEEDED:
library(RColorBrewer) # for plotting colors
library(broom) # helps manage shapefiles





# Load Data ---------------------------------------------------------------

# Data on estimated population size at each site with pursuit-divers on the Peninsula:
pop.pursuit <- read.csv(file = "pop.pursuit.csv", header = T)
# This file was created by manually combining the Pygoscelis penguin data from MAPPPD (see the Extract_MAPPPD R script) and the shag sites from Schrimpf et al. 2018.
pop.pursuit <- pop.pursuit[order(pop.pursuit$Code),] # Sort by site code (just for consistency)


# Turn it into a shapefile
colonies <- SpatialPointsDataFrame(coords = pop.pursuit[,c("Lon", "Lat")],
                                   data = pop.pursuit)
proj4string(colonies) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") # applying WGS84 projection







# Custom CRS --------------------------------------------------------------

# Here I create a custom CRS, designed to have a projection that will accurately reflect distance between points as well as possible. I chose to use a Lambert Equal Area projection, centered as close as possible in the middle of the Antarctic Peninsula colonies.

# To center it I calculate the midpoint in both Lat and Lon:
mean(range(coordinates(colonies)[,"Lat"]))
mean(range(coordinates(colonies)[,"Lon"]))

# Then I use them to define the projection (calling it the 'MAPPPD' projection, to distinguish it from my other work on ASI presence/Absence data):
MAPPPD.laea.crs <- CRS('+proj=laea +lat_0=-64.411 +lon_0=-56.789 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs')


# Reproject my data using my custom CRS:
colonies.laea <- spTransform(colonies, MAPPPD.laea.crs)








# Summarize Wind ----------------------------------------------------------

# Note, it takes a while to load and manage the wind files, so you can skip to the end of this section and re-load the completed results in the next section.


### The first part of this section takes the grib files obtained online and recalculates the wind as a direction, for the area of interest and for each layer (i.e., each month x year in the file):.

# Load raw wind data:
wind_u <- brick("wind_u_JFND_1979-2019.grib") # u-component wind
wind_v <- brick("wind_v_JFND_1979-2019.grib") # v-component wind

# change longitude from 0:360 to standard -180:180
wind_u <- rotate(wind_u)
wind_v <- rotate(wind_v)

# Trim the extent:
ex <- extent(-72, -38, -71, -60) # Extent limits
wind_u_cp <- crop(wind_u, ex)
rm(wind_u)
wind_v_cp <- crop(wind_v, ex)
rm(wind_v)

# Optional - plot the first layer, just to see it
plot(wind_u_cp[[1]], main = "u") # east:west (towards east is +)
plot(wind_v_cp[[1]], main = "v") # north:south (towards north is +)

# Reproject into custom CRS
wind_u_cp_laea <- projectRaster(wind_u_cp, crs = MAPPPD.laea.crs)
wind_v_cp_laea <- projectRaster(wind_v_cp, crs = MAPPPD.laea.crs)

# Convert into direction the wind is coming from:
wind_dir <- atan2(-wind_u_cp_laea, -wind_v_cp_laea) / (2 * pi) + 1
wind_dir <- (wind_dir - floor(wind_dir)) * 360

# Optional - plot the first layer, just to see it
plot(wind_dir[[3]], col = rainbow(256))

# Save as an rds object (that can be loaded later)
saveRDS(wind_dir, file = "wind_dir_JFND_1979-2019.rds")


### The next part of this section extracts the wind direction at the location of each bird colony, and creates a histogram of the prevailing direction from all layers combined

# Create Matrix for extracted values:
wind_extr <- matrix(nrow = length(pop.pursuit$Code), ncol = dim(wind_dir)[3])

# Add dimnames:
yrs <- 1979:2019 # Year span downloaded
month <- c("N","D","J","F") # Note: only summer months, Nov-Feb. (change if different set of months downloaded)
cols <- vector(length = length(yrs) * length(month))
for (i in 1:length(yrs)) {
  for (j in 1:length(month)) {
    cols[(i-1)*(length(month)) + j] <- paste(month[j], as.character(yrs[i]), sep = "_")
  }
}
dimnames(wind_extr) <- list(pop.pursuit$Code, cols)
rm(yrs, month, cols) # Clean up

# Extract month x year average wind direction for each site, from each layer:
for (i in 1:nrow(wind_extr)) {
  wind_extr[which(rownames(wind_extr) == colonies.laea[i,]$Code),] <- 
    as.numeric(extract(wind_dir, colonies.laea[i,]))
}

# Create histogram object:
wind_hist <- matrix(nrow = length(pop.pursuit$Code), ncol = 36)
dimnames(wind_hist) <- list(pop.pursuit$Code, to(prefix = "d", until = 35, from = 0))

# Summarize histogram for each site (by 10 deg bins)
for (i in 1:nrow(wind_hist)) {
  temp_hist <- hist(wind_extr[i,], breaks = c(seq(0, 360, length.out = 37)), plot = F)
  wind_hist[i,] <- temp_hist$counts / sum(temp_hist$counts) # provides proportion of total
}

# Export data
write.csv(wind_hist, file = "wind_bysite.csv", row.names = T)
write.csv(wind_extr, file = "wind_extract_byyear.csv", row.names = T)









# Reload Wind -------------------------------------------------------------

wind_dir <- readRDS(file = "wind_dir_JFND_1979-2019.rds") # Brick object with a raster of wind direction for each month x year
wind_extr <- read.csv(file = "wind_extract_byyear.csv", header = T, row.names = 1) # Extracted wind directions for at each bird ste for each month x year
wind_extr <- as.matrix(wind_extr)
wind_hist <- read.csv(file = "wind_bysite.csv", header = T, row.names = 1) # Proportion of the total time series (i.e. layers) that had prevailing wind direction from each 10 degree bin at each bird site







# Plotting Wind -----------------------------------------------------------

library(spatstat) # for the rose function (masks the rotate() functions needed to extract the wind data)


par(mar=c(0.1, 0.1, 2.1, 0.1))
site <- "MIKK" # Enter any site code
temp_hist <- hist(wind_extr[site,], breaks = c(seq(0, 360, length.out = 37)), plot = F)
rose(temp_hist, col = 'gray', main = site, clockwise = T, start = 270)
par(mar=c(5.1, 4.1, 4.1, 2.1)) # Return to normal margins








# Subset sites ------------------------------------------------------------

# Select only the sites in one region (the region names are listed in the MAPPPD location file)

# Central West
colonies.laea.CW <- subset(colonies.laea, Reg=="Central-west Antarctic Peninsula")










# Load coastline data -----------------------------------------------------


# Load the full peninsula coastline layer from the ADD (Large file!) for the actual fetch calculation:
# To obtain the file, go to the Antarctic Digital Database: https://www.add.scar.org/ and download the medium resolution vector polygon of the Antarctic Coastline. I used GIS to trim it to only include the Peninsula region and then project it in WGS84. This trimmed version is still too large for GitHub, but I can provide upon request.
peninsula <- readOGR(dsn = "./Shapefiles", layer = "Peninsula_WGS84poly")
peninsula.laea <- spTransform(peninsula, MAPPPD.laea.crs) # Transform to custom CRS


# Load the more simple GADM (https://gadm.org/) coast layer to use for easier visualization:
peninsula_GADM <- readOGR(dsn = "./Shapefiles", layer = "GADM_peninsula")
# Turn the shapefile into an object that ggplot2 can use:
peninsula_GADM <- fortify(peninsula_GADM)







# Load coastline points ---------------------------------------------------

# This set of points is spaced evenly around the coastline (every 200 m), each 1m from the ADD coast layer.
# These were produced and subsetted in QGIS, and then some were labeled with a site_assoc attribute linking them to specific bird sites
# It's best to explore them in GIS

points.laea <- readOGR(dsn = "./Shapefiles", layer = "Peninsula_WGS84poly_laea_points")

# Subsetting to get all points in the CW
# (using the x/y units from the custom CRS that roughly correspond to the boundaries of the region)
# I needed to figure those out by hand in QGIS
points.laea.CW <- points.laea[which(coordinates(points.laea)[,"coords.x1"] < -170000 & 
                                      coordinates(points.laea)[,"coords.x1"] > -375000 &
                                      coordinates(points.laea)[,"coords.x2"] < 90000 &
                                      coordinates(points.laea)[,"coords.x2"] > -139000),]

# Points in the CW and labeled with a site:
points.laea.CW_sites <- points.laea.CW[which(!is.na(points.laea.CW$site_assoc)),]

# Subsets from the two example regions:
points.laea.Trinity <- readOGR(dsn = "./Shapefiles", layer = "Peninsula_WGS84poly_laea_points_Trinity")
points.laea.Paradise <- readOGR(dsn = "./Shapefiles", layer = "Peninsula_WGS84poly_laea_points_Paradise")









# Selecting random points -------------------------------------------------

# This section selects a random subset of coastline points that are not near bird sites.
# Note, after running this section, I manually chose which six points around each randomly selected point would be used to average together to use as the set of random "sites".
# So if recreating my analysis exactly, you can skip this section and reload that file in the next section.

# Create 1 km buffer around all points associated with bird sites
site_point_buff <- buffer(points.laea.CW_sites, width = 1000, dissolve = F)

# Select the points not in those buffers:
points.laea.CW_empty_full <- gDifference(points.laea.CW,site_point_buff)

# Save the full set of points away from sites as a shapefile (in the custom CRS):
points.laea.CW_empty_full <- as(points.laea.CW_empty_full,"SpatialPointsDataFrame")
points.laea.CW_empty_full@data <- data.frame(ID=1:length(points.laea.CW_empty_full))
writeOGR(dsn = "./Shapefiles", points.laea.CW_empty_full, layer = "points.laea.CW_empty_full", driver = "ESRI Shapefile", overwrite_layer = TRUE)

# Choose a random subset of those points (the same size as the number of sites):
set.seed(73)
rdm <- sample(1:length(points.laea.CW_empty_full), size = dim(colonies.laea.CW)[1], replace = F)
points.laea.CW_empty <- points.laea.CW_empty_full[rdm,]

# Save that random subset as a shapefile:
writeOGR(dsn = "./Shapefiles", points.laea.CW_empty, layer = "points.laea.CW_empty", driver = "ESRI Shapefile", overwrite_layer = TRUE)

# I then used QGIS to manually select which of the neighboring points in points.laea.CW_empty_full would be treated as each empty "site", indicated by the points in points.laea.CW_empty.







# Reload random points ----------------------------------------------------

points.laea.CW_empty_assoc <- readOGR(dsn = "./Shapefiles", layer = "points.laea.CW_empty_assoc") # Just points associated with the random points









# Run fetch ---------------------------------------------------------------

# Here I run the fetch function for each point in several point sets
# It takes a while (it would go faster if parallelized) but if desired you can skip to the next section to reload the csv files and shapefiles produced


#### Sites with birds

# Create Object to hold results:
fetch_results_CW_sites <- matrix(nrow = length(points.laea.CW_sites),
                                 ncol = 36+3)
dimnames(fetch_results_CW_sites) <- list(NULL, c("site_assoc", "coordX", "coordY", to(prefix = "d", until = 35, from = 0)))


# Run fetch() on each point:
for(i in 1:nrow(fetch_results_CW_sites)) {
  temp <- fetch(polygon_layer = peninsula.laea,
                site_layer = points.laea.CW_sites[i,],
                max_dist = 500)
  fetch_results_CW_sites[i,4:39] <- temp[[1]]$fetch
  fetch_results_CW_sites[i,2:3] <- coordinates(points.laea.CW_sites)[i,]
  fetch_results_CW_sites[i,1] <- as.character(points.laea.CW_sites$site_assoc[i])
}

fetch_results_CW_sites <- as.data.frame(fetch_results_CW_sites)


# Save raw fetch results for points:
write.csv(fetch_results_CW_sites, file = "fetch_results_CW_sites.csv", row.names = F)


# Turn into a layer:
fetch_results_CW_sites.laea <- SpatialPointsDataFrame(coords = fetch_results_CW_sites[,c("coordX", "coordY")],
                                                      data = fetch_results_CW_sites)
proj4string(fetch_results_CW_sites.laea) <- MAPPPD.laea.crs # Applying cuystom CRS
writeOGR(dsn = "./Shapefiles", fetch_results_CW_sites.laea, layer = "fetch_results_CW_sites.laea", driver = "ESRI Shapefile", overwrite_layer = TRUE)





#### Random sites without birds

# Create Object to hold results:
fetch_results_CW_empty_assoc <- matrix(nrow = length(points.laea.CW_empty_assoc),
                                       ncol = 36+3)
dimnames(fetch_results_CW_empty_assoc) <- list(NULL, c("assoc", "coordX", "coordY", to(prefix = "d", until = 35, from = 0)))

# Run fetch() on each point:
for(i in 1:nrow(fetch_results_CW_empty_assoc)) {
  temp <- fetch(polygon_layer = peninsula.laea,
                site_layer = points.laea.CW_empty_assoc[i,],
                max_dist = 500)
  fetch_results_CW_empty_assoc[i,4:39] <- temp[[1]]$fetch
  fetch_results_CW_empty_assoc[i,2:3] <- coordinates(points.laea.CW_empty_assoc)[i,]
  fetch_results_CW_empty_assoc[i,1] <- points.laea.CW_empty_assoc$rdm_assoc[i]
}

fetch_results_CW_empty_assoc <- as.data.frame(fetch_results_CW_empty_assoc)


# Save raw fetch results for points:
write.csv(fetch_results_CW_empty_assoc, file = "fetch_results_CW_empty_assoc.csv", row.names = F)


# Turn into a layer:
fetch_results_CW_empty_assoc.laea <- SpatialPointsDataFrame(coords = fetch_results_CW_empty_assoc[,c("coordX", "coordY")],
                                                            data = fetch_results_CW_empty_assoc)
proj4string(fetch_results_CW_empty_assoc.laea) <- MAPPPD.laea.crs # Applying cuystom CRS
writeOGR(dsn = "./Shapefiles", fetch_results_CW_empty_assoc.laea, layer = "fetch_results_CW_empty_assoc.laea", driver = "ESRI Shapefile", overwrite_layer = TRUE)





#### Points from specific regions

### Trinity Island

# Create Object to hold results:
fetch_results_Trinity <- matrix(nrow = length(points.laea.Trinity),
                                ncol = 36+3)
dimnames(fetch_results_Trinity) <- list(NULL, c("site_assoc", "coordX", "coordY", to(prefix = "d", until = 35, from = 0)))


# Run fetch() on each point:
for(i in 1:nrow(fetch_results_Trinity)) {
  temp <- fetch(polygon_layer = peninsula.laea,
                site_layer = points.laea.Trinity[i,],
                max_dist = 500)
  fetch_results_Trinity[i,4:39] <- temp[[1]]$fetch
  fetch_results_Trinity[i,2:3] <- coordinates(points.laea.Trinity)[i,]
  fetch_results_Trinity[i,1] <- as.character(points.laea.Trinity$site_assoc[i])
}

fetch_results_Trinity <- as.data.frame(fetch_results_Trinity)


# Save raw fetch results for points:
write.csv(fetch_results_Trinity, file = "fetch_results_Trinity.csv", row.names = F)
# Not sure why, but when it gets reloaded the numerical columns get treated better when it's reloaded
fetch_results_Trinity <- read.csv(file = "fetch_results_Trinity.csv", header = T)


# Turn into a layer:
fetch_results_Trinity.laea <- SpatialPointsDataFrame(coords = fetch_results_Trinity[,c("coordX", "coordY")],
                                                     data = fetch_results_Trinity)
proj4string(fetch_results_Trinity.laea) <- MAPPPD.laea.crs # Applying cuystom CRS
writeOGR(dsn = "./Shapefiles", fetch_results_Trinity.laea, layer = "fetch_results_Trinity.laea", driver = "ESRI Shapefile", overwrite_layer = TRUE)



### Paradise Harbour

# Create Object to hold results:
fetch_results_Paradise <- matrix(nrow = length(points.laea.Paradise),
                                 ncol = 36+3)
dimnames(fetch_results_Paradise) <- list(NULL, c("site_assoc", "coordX", "coordY", to(prefix = "d", until = 35, from = 0)))


# Run fetch() on each point:
for(i in 1:nrow(fetch_results_Paradise)) {
  temp <- fetch(polygon_layer = peninsula.laea,
                site_layer = points.laea.Paradise[i,],
                max_dist = 500)
  fetch_results_Paradise[i,4:39] <- temp[[1]]$fetch
  fetch_results_Paradise[i,2:3] <- coordinates(points.laea.Paradise)[i,]
  fetch_results_Paradise[i,1] <- as.character(points.laea.Paradise$site_assoc[i])
}

fetch_results_Paradise <- as.data.frame(fetch_results_Paradise)


# Save raw fetch results for points:
write.csv(fetch_results_Paradise, file = "fetch_results_Paradise.csv", row.names = F)
# Not sure why, but when it gets reloaded the numerical columns get treated better when it's reloaded
fetch_results_Paradise <- read.csv(file = "fetch_results_Paradise.csv", header = T)


# Turn into a layer:
fetch_results_Paradise.laea <- SpatialPointsDataFrame(coords = fetch_results_Paradise[,c("coordX", "coordY")],
                                                      data = fetch_results_Paradise)
proj4string(fetch_results_Paradise.laea) <- MAPPPD.laea.crs # Applying cuystom CRS
writeOGR(dsn = "./Shapefiles", fetch_results_Paradise.laea, layer = "fetch_results_Paradise.laea", driver = "ESRI Shapefile", overwrite_layer = TRUE)









# Reload fetch results ----------------------------------------------------

# Bird sites:
fetch_results_CW_sites <- read.csv(file = "fetch_results_CW_sites.csv", header = T) # results dataframe
fetch_results_CW_sites.laea <- readOGR(dsn = "./Shapefiles", layer = "fetch_results_CW_sites.laea") # Shapefile

# Random sites without birds
fetch_results_CW_empty_assoc <- read.csv(file = "fetch_results_CW_empty_assoc.csv", header = T) # results dataframe
fetch_results_CW_empty_assoc.laea <- readOGR(dsn = "./Shapefiles", layer = "fetch_results_CW_empty_assoc.laea") # Shapefile

# Trinity Island:
fetch_results_Trinity <- read.csv(file = "fetch_results_Trinity.csv", header = T) # results dataframe
fetch_results_Trinity.laea <- readOGR(dsn = "./Shapefiles", layer = "fetch_results_Trinity.laea") # Shapefile

# Paradise Harbour:
fetch_results_Paradise <- read.csv(file = "fetch_results_Paradise.csv", header = T) # results dataframe
fetch_results_Paradise.laea <- readOGR(dsn = "./Shapefiles", layer = "fetch_results_Paradise.laea") # Shapefile








# Wind Correction ---------------------------------------------------------


# Function to correct fetch for a set of points by a RasterBrick with layers of a directional environmental layer (i.e. wind)
# Note, the fetch points need to be in a shapefile, and it needs to have the same CRS as the environmental layer
# Also note, the environmental layer needs to have data split into 10 degree bins

WindCor <- function(pts, brick) {
  # Create Matrix for extracted values:
  extr <- matrix(nrow = nrow(pts), ncol = nlayers(brick))
  
  # Extract each layer's average wind direction for each point:
  for (i in 1:nrow(extr)) {
    extr[i,] <- as.numeric(extract(brick, pts[i,]))
  }
  
  # Create matrix for histogram results:
  histo <- matrix(nrow = nrow(extr), ncol = 36)
  dimnames(histo) <- list(NULL, to(prefix = "d", until = 35, from = 0))
  
  # Summarize histogram for each point (by 10 deg bins) - using proportions:
  for (i in 1:nrow(histo)) {
    temp_hist <- hist(extr[i,], breaks = c(seq(0, 360, length.out = 37)), plot = F)
    histo[i,] <- temp_hist$counts / sum(temp_hist$counts) # provides proportion of total
  }
  
  # Getting just the data from the fetch shapefile
  pts.mat <- as.matrix(pts[,to(prefix = "d", until = 35, from = 0)]@data)
  
  # Scaling fetch by the histogram:
  res <- pts.mat * histo
  res <- as.data.frame(res)
  
  return(res)
  
}


# Run function for the points associated with bird sites:
fetch_WindCor_CW_sites <- WindCor(pts = fetch_results_CW_sites.laea, brick = wind_dir)
# Adding site associations back in:
fetch_WindCor_CW_sites$site_assoc <- fetch_results_CW_sites$site_assoc


# Run function for the randomly selected points away from bird sites:
fetch_WindCor_CW_empty <- WindCor(pts = fetch_results_CW_empty_assoc.laea, brick = wind_dir)
# Adding point associations back in:
fetch_WindCor_CW_empty$rdm_assoc <- fetch_results_CW_empty_assoc$assoc


# Run function for the Trinity Island points:
fetch_WindCor_Trinity <- WindCor(pts = fetch_results_Trinity.laea, brick = wind_dir)
# Adding site associations back in:
fetch_WindCor_Trinity$site_assoc <- fetch_results_Trinity$site_assoc


# Run function for the Paradise Harbour points:
fetch_WindCor_Paradise <- WindCor(pts = fetch_results_Paradise.laea, brick = wind_dir)
# Adding site associations back in:
fetch_WindCor_Paradise$site_assoc <- fetch_results_Paradise$site_assoc







# Summarize by histograms -------------------------------------------------


#### Bird sites

# This function is for calculating the average corrected fetch at all points associated with named breeding sites:
SummarySiteMean <- function(cor_fetch) {
  
  sites <- unique(cor_fetch$site_assoc)
  res <- data.frame(sites)
  names(res) <- "Code"
  res$mean_fetch <- NA
  
  # Runs summary
  for (i in 1:length(sites)) {
    temp <- cor_fetch[which(cor_fetch$site_assoc == sites[i]),
                      to(prefix = "d", until = 35, from = 0)]
    res$mean_fetch[i] <- mean(apply(temp, MARGIN = 1, FUN = sum))
  }
  
  return(res)
}

# Run function:
res_site_mean <- SummarySiteMean(cor_fetch = fetch_WindCor_CW_sites)
# Combine results objects:
pop.pursuit.fetch.CW <- merge(pop.pursuit, res_site_mean, by = "Code") # Note: could use the all.x=T argument to make it a full list

# Record number of points per site:
n_points_site <- as.data.frame(table(fetch_WindCor_CW_sites$site_assoc))
colnames(n_points_site) <- c("Code", "n_points")
pop.pursuit.fetch.CW <- merge(pop.pursuit.fetch.CW, n_points_site)
mean(pop.pursuit.fetch.CW$n_points) # Average number of points per site
# This was how we made the decision to use 7 points when creating the random locations

# Creating a dataframe to plot
plot.colonies.CW <- data.frame(colonies.laea.CW)
colnames(plot.colonies.CW)[colnames(plot.colonies.CW) %in% c("Lon.1","Lat.1")] <- c("Lon_laea","Lat_laea")
plot.colonies.CW <- merge(plot.colonies.CW, pop.pursuit.fetch.CW[,c("Code","mean_fetch","n_points")], by = "Code")




### Empty sites:

# Function calculates the average for each set of associated points:
SummaryEmptyMean <- function(cor_fetch) {
  
  locs <- unique(cor_fetch$rdm_assoc)
  res <- data.frame(locs)
  names(res) <- "Loc_ID"
  res$mean_fetch <- NA
  
  # Runs summary
  for (i in 1:length(locs)) {
    temp <- cor_fetch[which(cor_fetch$rdm_assoc == locs[i]),
                      to(prefix = "d", until = 35, from = 0)]
    res$mean_fetch[i] <- mean(apply(temp, MARGIN = 1, FUN = sum))
  }
  
  return(res)
}

# Run function
fetch_WindCor_CW_empty_mean <- SummaryEmptyMean(cor_fetch = fetch_WindCor_CW_empty)




#### Sub-regions:

# Note, these are different, because rather than do quantitative analysis with them, they simply were used to create maps in QGIS. The last line of each section saves a shapefile, which are in the Shapefiles folder - no need to run this code again, it's just here for documentation.


## Trinity Island

# Calculate total for every point:
fetch_WindCor_Trinity$Sum <- apply(fetch_WindCor_Trinity[,to(prefix = "d", until = 35, from = 0)], MARGIN = 1, FUN = sum)
fetch_WindCor_Trinity$LogSum <- log(fetch_WindCor_Trinity$Sum) # Also calculate log-transformed value

# Turn into a layer:
fetch_WindCor_Trinity$coordX <- fetch_results_Trinity$coordX
fetch_WindCor_Trinity$coordY<- fetch_results_Trinity$coordY
fetch_WindCor_Trinity.laea <- SpatialPointsDataFrame(coords = fetch_WindCor_Trinity[,c("coordX", "coordY")],
                                                     data = fetch_WindCor_Trinity)
proj4string(fetch_WindCor_Trinity.laea) <- MAPPPD.laea.crs # Applying cuystom CRS
writeOGR(dsn = "./Shapefiles", fetch_WindCor_Trinity.laea, layer = "fetch_WindCor_Trinity.laea", driver = "ESRI Shapefile", overwrite_layer = TRUE)


## Paradise Harbour:

# Calculate total for every point:
fetch_WindCor_Paradise$Sum <- apply(fetch_WindCor_Paradise[,to(prefix = "d", until = 35, from = 0)], MARGIN = 1, FUN = sum)
fetch_WindCor_Paradise$LogSum <- log(fetch_WindCor_Paradise$Sum)  # Also calculate log-transformed value

# Turn into a layer:
fetch_WindCor_Paradise$coordX <- fetch_results_Paradise$coordX
fetch_WindCor_Paradise$coordY<- fetch_results_Paradise$coordY
fetch_WindCor_Paradise.laea <- SpatialPointsDataFrame(coords = fetch_WindCor_Paradise[,c("coordX", "coordY")],
                                                      data = fetch_WindCor_Paradise)
proj4string(fetch_WindCor_Paradise.laea) <- MAPPPD.laea.crs # Applying cuystom CRS
writeOGR(dsn = "./Data/Shapefiles", fetch_WindCor_Paradise.laea, layer = "fetch_WindCor_Paradise.laea", driver = "ESRI Shapefile", overwrite_layer = TRUE)








# Analysis ----------------------------------------------------------------

# Index for sites with each species
GEPE.ind <- which(pop.pursuit.fetch.CW$PopEstGEPE > 0) # with Gentoos
GEPEn.ind <- which(!pop.pursuit.fetch.CW$PopEstGEPE > 0) # without Gentoos
CHPE.ind <- which(pop.pursuit.fetch.CW$PopEstCHPE > 0) # with Chinstraps
CHPEn.ind <- which(!pop.pursuit.fetch.CW$PopEstCHPE > 0) # without Chinstraps
ADPE.ind <- which(pop.pursuit.fetch.CW$PopEstADPE > 0) # with Adelies
ADPEn.ind <- which(!pop.pursuit.fetch.CW$PopEstADPE > 0) # without Adelies
ANSH.ind <- which(pop.pursuit.fetch.CW$PopEstANSH > 0) # with Shags
ANSHn.ind <- which(!pop.pursuit.fetch.CW$PopEstANSH > 0) # without Shags


# Examining distributions of all bird sites:
hist(pop.pursuit.fetch.CW$mean_fetch, col = "cornflowerblue")
hist(log(pop.pursuit.fetch.CW$mean_fetch), col = "cornflowerblue")



#### Are the distributions actually log-normal?

# Function to run a test of log-normality:
LNormTest <- function(dat, species) {
  if(species == "GEPE") {
    dat.sub <- dat[GEPE.ind]
    dat.subn <- dat[GEPEn.ind]
  } else {
    if(species == "CHPE") {
      dat.sub <- dat[CHPE.ind]
      dat.subn <- dat[CHPEn.ind]
    } else {
      if(species == "ADPE") {
        dat.sub <- dat[ADPE.ind]
        dat.subn <- dat[ADPEn.ind]
      } else {
        if(species == "ANSH") {
          dat.sub <- dat[ANSH.ind]
          dat.subn <- dat[ANSHn.ind]
        } else {
          if(species == "All") {
            dat.sub <- dat
            dat.subn <- NULL
          } else {
            return(paste("Enter species"))
          }
        }
      }
    }
  }
  
  qqnorm(log(dat.sub), main = species)
  qqline(log(dat.sub), col = "gray", lwd = 2)
  
  sub <- ks.test(x = log(dat.sub),
                 y = pnorm,
                 mean = mean(log(dat.sub)),
                 sd = sd(log(dat.sub)))
  
  
  if(species != "All") {
    qqnorm(log(dat.subn), main = paste("no", species, sep = " "))
    qqline(log(dat.subn), col = "gray", lwd = 2)
    
    paste("Without:")
    subn <- ks.test(x = log(dat.subn),
                    y = pnorm,
                    mean = mean(log(dat.subn)),
                    sd = sd(log(dat.subn)))
    res <- list(sub, subn)
  } else {res <- sub}
  
  return(res)
  
}


# Run function on each species:
LNormTest(dat = pop.pursuit.fetch.CW$mean_fetch,
          species = "GEPE")
LNormTest(dat = pop.pursuit.fetch.CW$mean_fetch,
          species = "CHPE")
LNormTest(dat = pop.pursuit.fetch.CW$mean_fetch,
          species = "ADPE")
LNormTest(dat = pop.pursuit.fetch.CW$mean_fetch,
          species = "ANSH")
LNormTest(dat = pop.pursuit.fetch.CW$mean_fetch,
          species = "All")





#### Compare sites with and without each species:

# Function to create boxplots comparing sites with and without a species:
SpBoxFun <- function(dat, species) {
  
  if(species == "GEPE") {
    dat.sub <- dat[GEPE.ind]
    dat.subn <- dat[GEPEn.ind]
  } else {
    if(species == "CHPE") {
      dat.sub <- dat[CHPE.ind]
      dat.subn <- dat[CHPEn.ind]
    } else {
      if(species == "ADPE") {
        dat.sub <- dat[ADPE.ind]
        dat.subn <- dat[ADPEn.ind]
      } else {
        if(species == "ANSH") {
          dat.sub <- dat[ANSH.ind]
          dat.subn <- dat[ANSHn.ind]
        } else {
          return(paste("Enter species"))
        }
      }
    }
  }
  
  n <- c(length(dat.sub), length(dat.subn))
  
  
  boxplot(log(dat.sub), log(dat.subn),
          horizontal = F,
          col = "steelblue",
          names = c(paste("With\n(",n[1],")", sep = ""),
                    paste("Without\n(",n[2],")", sep = "")),
          las = 1,
          main = species,
          ylim = c(-2,6),
          yaxt = "n",
          cex = 0.75)
  
}


# Function to compare all sites with birds to the random sites wwithout birds:
FullBoxFun <- function(dat, rdm) {
  
  boxplot(log(dat), log(rdm),
          horizontal = F,
          col = "orchid3",
          names = c("Full","Random"),
          las = 1,
          main = "",
          ylim = c(-2,6),
          yaxt = "n",
          cex = 0.75)
}



# Production of figure (as svg):

svg(filename = "Graphs/barplots.svg", width = 6.5, height = 4)

par(mfrow = c(1, 5), oma=c(0,3,0,1), mgp = c(3, 1.5, 0))

par(mar=c(3.1, 1.1, 3.1, 0.1))
SpBoxFun(dat = pop.pursuit.fetch.CW$mean_fetch,
         species = "GEPE")
axis(side = 2, las = 2)
SpBoxFun(dat = pop.pursuit.fetch.CW$mean_fetch,
         species = "CHPE")
SpBoxFun(dat = pop.pursuit.fetch.CW$mean_fetch,
         species = "ADPE")
SpBoxFun(dat = pop.pursuit.fetch.CW$mean_fetch,
         species = "ANSH")
FullBoxFun(dat = pop.pursuit.fetch.CW$mean_fetch,
           rdm = fetch_WindCor_CW_empty_mean$mean_fetch)

dev.off()




## T tests

# Calculate:
t.GEPE <- t.test(log(pop.pursuit.fetch.CW$mean_fetch[GEPE.ind]),
                 log(pop.pursuit.fetch.CW$mean_fetch[GEPEn.ind]))
t.CHPE <- t.test(log(pop.pursuit.fetch.CW$mean_fetch[CHPE.ind]),
                 log(pop.pursuit.fetch.CW$mean_fetch[CHPEn.ind]))
t.ADPE <- t.test(log(pop.pursuit.fetch.CW$mean_fetch[ADPE.ind]),
                 log(pop.pursuit.fetch.CW$mean_fetch[ADPEn.ind]))
t.ANSH <- t.test(log(pop.pursuit.fetch.CW$mean_fetch[ANSH.ind]),
                 log(pop.pursuit.fetch.CW$mean_fetch[ANSHn.ind]))
t.All <- t.test(log(pop.pursuit.fetch.CW$mean_fetch),
                log(fetch_WindCor_CW_empty_mean$mean_fetch))

# Putting p-values together
t.res <- data.frame(c("GEPE","CHPE","ADPE","ANSH","All"))
colnames(t.res) <- "Species"
t.res$pval <- c(t.GEPE$p.value,
                t.CHPE$p.value,
                t.ADPE$p.value,
                t.ANSH$p.value,
                t.All$p.value)
# Function to calculate difference between means
DiffFun <- function(res) {
  return(as.numeric(res$estimate[1] - res$estimate[2]))
}
t.res$Diff <- c(DiffFun(t.GEPE),
                DiffFun(t.CHPE),
                DiffFun(t.ADPE),
                DiffFun(t.ANSH),
                DiffFun(t.All))
# t-statistics
t.res$t_stat <- c(t.GEPE$statistic,
                  t.CHPE$statistic,
                  t.ADPE$statistic,
                  t.ANSH$statistic,
                  t.All$statistic)
# Corrected with Bonferroni method
t.res$pval_Bon <- p.adjust(t.res$pval, method = "bonferroni")
# Corrected with Sidak method
t.res$pval_Sid <- gsidak.p.adjust(t.res$pval)









# Maps --------------------------------------------------------------------

# Transform points to WGS84 (to use Lat Long for plotting):
points.WGS84.CW <- spTransform(points.laea.CW, CRS("+init=epsg:4326"))
colonies.WGS84.CW <- spTransform(colonies.laea.CW, CRS("+init=epsg:4326"))

# Make points into dataframe for plotting:
plot.points.CW <- data.frame(points.WGS84.CW)

# Create background map:
basic <- ggplot() +
  geom_polygon(data = peninsula_GADM, aes(x = long, y = lat, group = group),
               size = 0.25, fill = "gray60", alpha = 1) +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "gray40", linetype = "dotted")) +
  labs(x = "", y = "") +
  coord_map(project="stereographic", orientation=c(-90,-62,0),
            ylim = c(-65.5, -63.5), xlim = c(-65, -60))


# Map with data (saved as an svg):
svg(filename = "Maps/Fetch_map.svg", width = 6.5, height = 6)

basic + 
  geom_point(data = plot.colonies.CW, aes(x = Lon, y = Lat, colour = log(mean_fetch)),
             size = 2.5) +
  scale_color_scico(palette = "imola") + # Which color scale to use
  labs(colour = "log(fetch)") # Legend title

dev.off()








# APPENDIX: Wind trends ---------------------------------------------------

# This is code used in the appendix to briefly examine wind trends

# First requires that the wind_extr object is loaded from the 'Reload Wind' section above.

# An index for the year for each column:
year_ind <- as.numeric(substr(colnames(wind_extr), start = 3, stop = 6))

# An index for the month for each column:
month_ind <- substr(colnames(wind_extr), start = 1, stop = 1)


library(circular) # Note, this masks sd() and var()

# Matrix for annual summer wind by site:
wind_summer <- matrix(nrow = nrow(wind_extr), ncol = length(unique(year_ind)),
                      dimnames = list(rownames(wind_extr), unique(year_ind)))


# Fill matrix with appropriate average:
for (i in 1:nrow(wind_summer)) {
  for (j in 1:ncol(wind_summer)) {
    temp <- wind_extr[i, which(year_ind == as.numeric(colnames(wind_summer)[j]))] # select  monthly values for the appropraite site x year
    temp_c <- circular(temp, units = "degrees") # converts selection to a circular object (keeping units as degrees)
    temp_d <- mean.circular(temp_c) # Calculates the circular mean for that site x year
    if(temp_d < 0) {temp_d <- temp_d + 360}
    wind_summer[i,j] <- temp_d
  }
}



# Many sites have the same underlying wind conditions, because the grid cell size for wind was fairly large.
# Only keep the unique rows:
wind_summer_uni <- unique(wind_summer)
nrow(wind_summer_uni) # The number of unique wind time series





### Circular regression for focal regions:


# Function to test for simple linear trend in wind direction over time series for specific site, and plot wind and trendline:
# (Note: it repeats the data twice, stacked on the y-axis, so that any trends that pass through North can be visualized the same as those that pass through any other portion of the compass rose)
Circ_LM <- function(site) {
  
  y <- circular(wind_summer[site,], units = "degrees") # Dependent (response) variable - circular
  x <- 1:length(y) # Independent (predictor) variable - linear (here, just year)
  
  mod <- lm.circular(type = "c-l", y = y, x = x, init = c(0)) # Model
  
  dat <- wind_summer[site,]
  dat_high <- wind_summer[site,] + 360 # repeat of the data, to plot above
  
  # Plotting lower data
  plot(dat, pch = 16,
       xaxt = 'n', yaxt = 'n',
       ann = F, ylim = c(0,720))
  
  # Plotting upper data
  points(x= 1:length(dat_high), y = dat_high, pch = 16)
  
  # Axes, etc.
  axis(side = 1)
  axis(side = 2, labels = c(seq(0,360,90), seq(90,360,90)),
       at = seq(0,720,90))
  title(xlab = paste("Years, site = ",site), ylab = "Mean Summer Wind Direction")
  abline(h = 360) # Zero line, separating upper and lower data
  
  # Trendline:
  abline(a = as.numeric(mod$mu), b = (((2*atan(mod$coefficients)) * 180) / pi), lty = 2 )
  abline(a = as.numeric(mod$mu)+360, b = (((2*atan(mod$coefficients)) * 180) / pi), lty = 2 )
  abline(a = as.numeric(mod$mu)+720, b = (((2*atan(mod$coefficients)) * 180) / pi), lty = 2 )
  
  return(mod)
  
}



## Paradise Harbor region:

# Wind data comparison (these are all of the sites in that area):
wind_summer[c("ONEI",
              "LAUT",
              "SANE",
              "DUPT",
              "WATE",
              "LEIC",
              "ALMI",
              "BRYE",
              "BRYS",
              "PABE"),1:3]
# So, only two distinct wind pixels: we'll use LAUT and ALMI for zones 1 and 2

circmod_LAUT <- Circ_LM("LAUT")
circmod_LAUT # Zone 1

circmod_ALMI <- Circ_LM("ALMI")
circmod_ALMI # Zone 2




## Trinity Island region:

# Wind data comparison:
wind_summer[c("MEGA",
              "TINW",
              "NWSP",
              "TIWE",
              "FARE",
              "SPER",
              "TISW",
              "SKOT",
              "MIKK",
              "TETR",
              "TIEA"),1:3]
# So, only three distinct wind pixels, we'll use TINW, FARE, and MIKK for zones 3-5.

circmod_TINW <- Circ_LM("TINW")
circmod_TINW # Zone 3

circmod_FARE <- Circ_LM("FARE")
circmod_FARE # Zone 4

circmod_MIKK <- Circ_LM("MIKK")
circmod_MIKK # Zone 5






# This is the same as the Circ_LM function above, but it does not return the model object, and adds the ability to include a custom range for the y-axis, so that the data do not need to be repeated.
# 'low' and 'high' should be the lower and upper y axis limits (high should be the angle + 360 if you need the second plot stacked above, and 'int' should be the interval for y axis tick marks

CircPlot_custom <- function(site, low, high, int) {
  
  y <- circular(wind_summer[site,], units = "degrees") # Dependent (response) variable - circular
  x <- 1:length(y) # Independent (predictor) variable - linear (here, just year)
  
  mod <- lm.circular(type = "c-l", y = y, x = x, init = c(0)) # Model
  
  dat <- wind_summer[site,]
  dat_high <- wind_summer[site,] + 360 # repeat of the data, to plot above
  
  # Plotting lower data
  plot(dat, pch = 16,
       xaxt = 'n', yaxt = 'n', yaxs = "i",
       ann = F, ylim = c(low,high))
  
  # Plotting upper data
  points(x= 1:length(dat_high), y = dat_high, pch = 16)
  
  # Axes, etc.
  axis(side = 1)
  axis(side = 2, labels = c(seq(0,360,int), seq(int,360,int)),
       at = seq(0,720,int))
  title(xlab = paste("Years, site = ",site), ylab = "Mean Summer Wind Direction")
  if (low < 360 & high > 360) {abline(h = 360, col = "gray")} # Zero line, separating upper and lower data
  
  # Trendline:
  abline(a = as.numeric(mod$mu), b = (((2*atan(mod$coefficients)) * 180) / pi), lty = 2 )
  abline(a = as.numeric(mod$mu)+360, b = (((2*atan(mod$coefficients)) * 180) / pi), lty = 2 )
  abline(a = as.numeric(mod$mu)+720, b = (((2*atan(mod$coefficients)) * 180) / pi), lty = 2 )
  
}



### Figure for the appendix

# Modified from above but without axis labels:
CircPlot_custom_noLab <- function(site, low, high, int) {
  
  y <- circular(wind_summer[site,], units = "degrees") # Dependent (response) variable - circular
  x <- 1:length(y) # Independent (predictor) variable - linear (here, just year)
  
  mod <- lm.circular(type = "c-l", y = y, x = x, init = c(0)) # Model
  
  dat <- wind_summer[site,]
  dat_high <- wind_summer[site,] + 360 # repeat of the data, to plot above
  
  # Plotting lower data
  plot(dat, pch = 16,
       xaxt = 'n', yaxt = 'n', yaxs = "i",
       ann = F, ylim = c(low,high))
  
  # Plotting upper data
  points(x= 1:length(dat_high), y = dat_high, pch = 16)
  
  # Axes, etc.
  axis(side = 1)
  axis(side = 2, labels = c(seq(0,360,int), seq(int,360,int)),
       at = seq(0,720,int), las = 1)
  if (low < 360 & high > 360) {abline(h = 360, col = "gray")} # Zero line, separating upper and lower data
  
  # Trendline:
  abline(a = as.numeric(mod$mu), b = (((2*atan(mod$coefficients)) * 180) / pi), lty = 2 )
  abline(a = as.numeric(mod$mu)+360, b = (((2*atan(mod$coefficients)) * 180) / pi), lty = 2 )
  abline(a = as.numeric(mod$mu)+720, b = (((2*atan(mod$coefficients)) * 180) / pi), lty = 2 )
  
}


svg(filename = "trends.svg", width = 6.5, height = 6.5)

par(mfrow = c(3, 2), oma=c(3,3,1,1))

plot(0, type='n', axes=FALSE, ann=FALSE) # I want the upper left to be blank

par(mar = c(3.1, 3.1, 1.1, 0.1))
# Zone 3 (using TINW)
CircPlot_custom_noLab(site = "TINW", low = 180, high = (360+180), int = 90)

# Zone 1 (using LAUT)
CircPlot_custom_noLab(site = "LAUT", low = 0, high = 360, int = 90)

# Zone 4 (using FARE)
CircPlot_custom_noLab(site = "FARE", low = 180, high = (360+180), int = 90)

# Zone 2 (using ALMI)
CircPlot_custom_noLab(site = "ALMI", low = 0, high = 360, int = 90)

# Zone 5 (using MIKK)
CircPlot_custom_noLab(site = "MIKK", low = 180, high = (360+180), int = 90)

dev.off()




### Compile results

library(MHTmult) # for the Sidak corrections


lm.res <- data.frame(c("Z1","Z2","Z3","Z4","Z5"))
colnames(lm.res) <- "Zone"

lm.res$slope_rad <- c(circmod_LAUT$coefficients,
                      circmod_ALMI$coefficients,
                      circmod_TINW$coefficients,
                      circmod_FARE$coefficients,
                      circmod_MIKK$coefficients) # Annual change in direction in radians

lm.res$slope_deg <- (((2*atan(lm.res$slope_rad)) * 180) / pi) # Annual change in direction in degrees

lm.res$t_stat <- c(circmod_LAUT$t.values,
                   circmod_ALMI$t.values,
                   circmod_TINW$t.values,
                   circmod_FARE$t.values,
                   circmod_MIKK$t.values)

lm.res$pval <- c(circmod_LAUT$p.values,
                 circmod_ALMI$p.values,
                 circmod_TINW$p.values,
                 circmod_FARE$p.values,
                 circmod_MIKK$p.values)

lm.res$pval_Sid <- gsidak.p.adjust(lm.res$pval)










