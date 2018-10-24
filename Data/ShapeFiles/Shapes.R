#---------------------------------------------------------------------
# 
# Shape File reading
# 
# Author: Samuel Johnson and Michelle Jones
# Date: 19/1/15
# 
# This script contains the functions for the simulation model project
# 
#---------------------------------------------------------------------

# setwd("/Users/michellejones/Dropbox/Sam_Inlets")
setwd("~/Dropbox/Work/Sam_Inlets")

# Load packages (these probably aren't all needed)
library(geosphere)
library(mapdata)
library(rgeos) 
library(maptools)
library(raster)
library (sp) #
library(ggplot2) 
library(ggmap) 
library(maps) 
library(rgdal) 

# Read in tagging data, subset to inlet releases
sable <- read.csv("recoveries.csv")
inlets <- subset(sable, release_set_type=="INLET STANDARDIZED")
# head(inlets)

# Read in shape file
management_areas <- readShapeSpatial("ShapeFiles/MajorMinorSQL_geo.shp")

# Extract major/minor area codes and the centroids of each
majCodes <- management_areas $ MAJOR_CODE
minCodes <- management_areas $ MINOR_CODE
centroids <- coordinates(management_areas)

# Create a data frame to hold them
MgmtAreaCentroids <- data.frame(  MAJOR = majCodes, MINOR = minCOdes, 
                                        centroids = centroids )

# Remove rows from inlets which have no recovery
rows <- which ( inlets $ recovery_Major != "" )
inlets <- inlets [ rows, ]
plot( inlets $ recovery_Major )

# okay, the NAs in the liberty years vector are screwing me up, so
rows <- which ( is.na ( inlets $ liberty_years ) )
inlets <- inlets [ -rows, ]

# YALfrequencies is a function which will create a table
# of frequencies for years at liberty (YAL) for recoveries in each major/minor
# area pair.
#     Args:     plusGroup = point at which to accumulate larger years
#               minYear = smallest value years at liberty
#               areas = a data frame containing maj/min area codes
#     Returns:  MgmtAreaFreqs = the areas df with the frequencies addended
YALfrequencies <- function (  plusGroup = 5, minYear = 0, 
                              areas = MgmtAreaCentroids, recoveries = inlets )
{
  # Get length of areas df
  nAreas <- nrow ( areas )
  bins <- length( minYear:plusGroup )

  # Define a matrix to hold the frequencies, and name the columns
  freqs <- matrix ( 0, nrow = nAreas, ncol = bins )
  pgName <- paste ( plusGroup, "+", sep = "" )
  freqCols <- c ( minYear:( plusGroup - 1 ), pgName )
  colnames ( freqs ) <- freqCols

  # Now to extract the frequencies for each major/minor pair
  areaCodes <- areas [ , 1:2 ]
  for ( mmArea in 1:nAreas ) {
    majCode <- areaCodes [ mmArea, 1 ]
    minCode <- areaCodes [ mmArea, 2 ]
    rows <- which ( inlets $ recovery_Major_cde == majCode & 
                    inlets $ recovery_Minor_cde == minCode )
    if ( length(rows) > 0 ) {
      areaRecs <- inlets [ rows, ]
      for ( bin in 1:(bins - 1)  ) {
        binRows <- which ( areaRecs $ liberty_years == (bin - 1) )
        freqs [ mmArea, bin ] <- nrow ( areaRecs [ binRows, ] )
      }
      plusRows <- which ( areaRecs $ liberty_years >= plusGroup )
      freqs [ mmArea, bins ] <- nrow ( areaRecs [ plusRows, ] )
    }
  }

  # Append frequencies areas
  areas <- cbind ( areas, freqs )

  # Return areas
  return ( areas )
}


releaseLocations <- data.frame  ( 
                                  long = -1 * inlets$release_longitude, 
                                  lat = inlets $ release_latitude 
                                )

# management_areas@data
plot(management_areas )
  points ( centroids )

