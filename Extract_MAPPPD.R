# Script to estimate Pygoscelis population size from Antarctic Peninsula MAPPPD sites



# Libraries ---------------------------------------------------------------

library(MCMCvis) # for the MCMCpstr function



# Load Data ---------------------------------------------------------------


# The first step here is to obtain data from the MAPPPD website: http://www.penguinmap.com/

# I then used three sources of information:

# Site metadata
site.meta <- read.csv("SiteLocations.csv")

# Raw count data
dat.ADPE.raw <- read.csv("ADPE_data_global_4January2016.csv", header = T)
dat.CHPE.raw <- read.csv("CHPE_data_global_4January2016.csv", header = T)
dat.GEPE.raw <- read.csv("GEPE_data_global_4January2016.csv", header = T)

# Model posteriors (code adapted from Lynch-Lab GitHub 'Tools' repository)
ExtractZ <- function(species){
  # First Download correct data
  if(species == "ADPE") {
    load("Adelie_MCMCzstate.RData")
    load("Adelie_SiteList.RData")
  } else if(species == "CHPE") {
    load("Chinstrap_MCMCzstate.RData")
    load("Chinstrap_SiteList.RData")
  } else if(species == "GEPE") {
    load("Gentoo_MCMCzstate.RData")
    load("Gentoo_SiteList.RData")
  } else {print("Not a correct species code.")}
  zstateList <- MCMCpstr(MCMCzstate, func = function(x) quantile(x, probs = c(.5))) # Extract the posteriors
  zstate <- zstateList[[1]] # same as a normal array
  rownames(zstate) <- SiteList$site_id # Add site names
  dimnames(zstate)[[2]] <- seq(1982,1982+dim(zstate)[2]-1) # Add season names
  return(zstate)
}

# Run extraction of model posteriors for each species
MAPPPD.ADPE.zstate <- ExtractZ("ADPE")
MAPPPD.CHPE.zstate <- ExtractZ("CHPE")
MAPPPD.GEPE.zstate <- ExtractZ("GEPE")





# Data Management ---------------------------------------------------------


# Function to extract population estimate for each site
PopEst <- function(raw, state) {
  
  sites <- unique(raw$CODE) # List of sites
  dat <- data.frame(rep(NA, times = length(sites))) # dataframe to build up
  names(dat) <- "Code"
  dat$Code <- sites
  
  # Pulls metadata from site.meta object (loaded seperately)
  for (i in 1:length(dat$Code)) {
    j <- which(raw$CODE == as.character(dat$Code[i]))[1]
    dat$Site[i] <- as.character(raw$SITE[j])
    k <- which(site.meta$site_id == as.character(dat$Code[i]))
    dat$Lat[i] <- site.meta$latitude[k]
    dat$Lon[i] <- site.meta$longitude[k]
    dat$Reg[i] <- as.character(site.meta$region[k])
  }
  dat <- dat[dat$Lat > -68.5 & dat$Lon > -70 & dat$Lon < -44.3,] # subsetting sites to Antarctic Peninsula area
  
  # Extract last count columns
  for (i in 1:nrow(dat)) {
    tmp.visits <- which(raw$CODE == paste(dat$Code[i]))
    tmp.last <- which.max(raw$YEAR[tmp.visits])
    dat$LastCt[i] <- as.numeric(raw$COUNT[tmp.visits[tmp.last]])
    dat$LastCtTp[i] <- paste(raw$WHAT[tmp.visits[tmp.last]])
    dat$LastCtYr[i] <- paste(raw$YEAR[tmp.visits[tmp.last]])
  }
  
  # Extract the posterior estimate (2016 median)
  dat$ModMed <- NA
  for (i in 1:length(dimnames(state)[[1]])) {
    tmp.index <- which(dat$Code == dimnames(state)[[1]][i])
    dat$ModMed[tmp.index] <- state[i,"2016"]
  }
  
  # Deciding what number to use for the population
  dat$PopEst <- NA
  for (i in 1:nrow(dat)) {
    
    if(dat$LastCtYr[i] >= 2006 &
       !is.na(dat$ModMed[i])) { # Recent counts with MAPPPD results:
      dat$PopEst[i] <- dat$ModMed[i] # use MAPPPD
      
    } else { # Older counts or those without MAPPPD results:
      
      if(dat$LastCtTp[i] == "N") { # Nests:
        dat$PopEst[i] <- dat$LastCt[i] # use last count
      }
      
      if(dat$LastCtTp[i] == "C" | dat$LastCtTp[i] == "A") { # Adults or Chicks:
        
        # Find earlier counts
        tmp.raw <- raw[which(raw$CODE == paste(dat$Code[i])),]
        if(any(tmp.raw$WHAT == "N")) { # if any ealier N counts exist, average all from most recent season:
          dat$PopEst[i] <- mean(tmp.raw$COUNT[which(
            tmp.raw$YEAR == max(tmp.raw$YEAR[which(tmp.raw$WHAT == "N")])
          )])
        } else { # If no earlier nest counts:
          dat$PopEst[i] <- dat$LastCt[i] / 1.5 # divide last count by 1.5
        }
        
      }
      
    }
    # Any remaining counts will simply remain NA
  }
  
  return(dat)
  
}


# Run function for each species
pop.ADPE <- PopEst(raw = dat.ADPE.raw, state = MAPPPD.ADPE.zstate)
pop.CHPE <- PopEst(raw = dat.CHPE.raw, state = MAPPPD.CHPE.zstate)
pop.GEPE <- PopEst(raw = dat.GEPE.raw, state = MAPPPD.GEPE.zstate)








# Combining all pygoscelid data -------------------------------------------
# Also don't need to do this - just load results (see below)


pop.ADPE <- read.csv(file = "pop.ADPE.csv", header = TRUE)
pop.CHPE <- read.csv(file = "pop.CHPE.csv", header = TRUE)
pop.GEPE <- read.csv(file = "pop.GEPE.csv", header = TRUE)
# Site metadata
site.meta <- read.csv("SiteLocations.csv")

# extracting needed info and recombining (this part is ugly, but it works):
tempADPE <- pop.ADPE[c(1,10)]
colnames(tempADPE) <- c("Code", "PopEstADPE")
tempCHPE <- pop.CHPE[c(1,10)]
colnames(tempCHPE) <- c("Code", "PopEstCHPE")
tempGEPE <- pop.GEPE[c(1,10)]
colnames(tempGEPE) <- c("Code", "PopEstGEPE")
tempSites <- site.meta[,c(1:3,5,6)]
colnames(tempSites) <- c("Code", "Site", "Reg", "Lat", "Lon")
tempPygo <- merge(merge(tempADPE, tempCHPE, by = "Code", all = TRUE),
                  tempGEPE, by = "Code", all = TRUE)
pop.pygo <- merge(tempSites, tempPygo, by = "Code")
rm(tempADPE, tempCHPE, tempGEPE, tempSites, tempPygo)

# Create sum column:
pop.pygo$PopEstPygo <- NA
for (i in 1:nrow(pop.pygo)) {
  pop.pygo$PopEstPygo[i] <- sum(pop.pygo$PopEstADPE[i],
                                pop.pygo$PopEstCHPE[i],
                                pop.pygo$PopEstGEPE[i],
                                na.rm = TRUE)
}

# Turn all population NAs to zeros:
pop.pygo$PopEstADPE[which(is.na(pop.pygo$PopEstADPE))] <- 0
pop.pygo$PopEstCHPE[which(is.na(pop.pygo$PopEstCHPE))] <- 0
pop.pygo$PopEstGEPE[which(is.na(pop.pygo$PopEstGEPE))] <- 0
pop.pygo$PopEstPygo[which(is.na(pop.pygo$PopEstPygo))] <- 0

# # Export as csv
# write.csv(pop.pygo, file = "pop.Pygo.csv", row.names = F)




# The pop.pygo object was then used to create the pop.pursuit.csv file in the main folder by combining with the data on Antarctic Shags manually.


