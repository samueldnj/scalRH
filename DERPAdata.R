# -----------------------------------------------------------------------------
#
# DERPAdata.R
#
# Reads in biological and survey CPUE data from the ./biology_with_age/
# and ./density/ folders inside the working directory. Used to plot CPUE
# time series, catch locations from the survey, and age compositions
# 
# -----------------------------------------------------------------------------

rm(list = ls())

library(dplyr)
library(ggmap)
library(maptools)
library(rgdal)
library( PBSmapping ) 
library( RColorBrewer )
library( RgoogleMaps)
library( sp )
library(raster)
library(MASS)
library(psyphy)
library(boot)
library(RCurl)

# Load shape file data
data(nepacLL)
load( file="bcAreas.RData" )

source("mseRtools.r")
source("DERPAfuns.R")

# Read in management areas shape file
readShapeSpatial("./ShapeFiles/MajorMinorSQL_geo.shp")
mgmtAreas <- readOGR(dsn = "./ShapeFiles/", layer = "MajorMinorSQL_geo")
proj4string(mgmtAreas)
mgmtAreas <- spTransform (mgmtAreas, CRS("+proj=longlat +datum=WGS84"))

# Read in density
doverDensity    <- read.csv("./density/dover_density.csv",header=T)
englishDensity  <- read.csv("./density/english_density.csv",header=T)
rockDensity     <- read.csv("./density/srock_density.csv",header=T)
petraleDensity  <- read.csv("./density/petrale_density.csv",header=T)
atoothDensity   <- read.csv("./density/atooth_density.csv",header=T)

# Read in strata areas
strata <- read.csv( "derpa_strata.csv", header=T )
stratArea <-  strata %>%
              dplyr::select( GROUPING_CODE, AREA_KM2 )
# Survey ids for plotting/legends
survIDs <-  c(  QCSyn = 1, HSAss = 2, HSSyn = 3, WCVISyn=4,
                QCShr = 6, WCVIShr = 7, WCHGSyn = 16 )

specNames <- c("Dover","English","Rock","Petrale","Arrowtooth")

dover <- dataPlots("dover", save = F, maxMatAge = 5 )
# english <- dataPlots("english", save = F, maxMatAge = 5 )
# srock <- dataPlots("srock", save = F, maxMatAge = 10)
# petrale <- dataPlots("petrale", save = F, maxMatAge = 15)
# atooth <- dataPlots("atooth", save = F, maxMatAge = 10)


# Combine density frames based on trip and trawl IDs
# First, rename the catch column in each df
doverDensity <- doverDensity %>% mutate(  doverCatch = CATCH_WEIGHT,
                                          doverDensity = DENSITY_KGPM2 )
englishDensity <- englishDensity %>% mutate(  englishCatch = CATCH_WEIGHT,
                                              englishDensity = DENSITY_KGPM2 )
rockDensity <- rockDensity %>% mutate(  rockCatch = CATCH_WEIGHT,
                                        rockDensity = DENSITY_KGPM2 )
petraleDensity <- petraleDensity %>% mutate(  petraleCatch = CATCH_WEIGHT,
                                              petraleDensity = DENSITY_KGPM2 )
atoothDensity <- atoothDensity %>% mutate(  atoothCatch = CATCH_WEIGHT,
                                            atoothDensity = DENSITY_KGPM2 )
# Now join and select the columns we want
surveyData <- doverDensity %>%
              left_join( englishDensity, by = "FISHING_EVENT_ID" ) %>%
              left_join( rockDensity, by = "FISHING_EVENT_ID" ) %>%
              left_join( petraleDensity, by = "FISHING_EVENT_ID" ) %>%
              left_join( atoothDensity, by = "FISHING_EVENT_ID" ) %>%
              left_join( stratArea, by = "GROUPING_CODE") %>%
              dplyr::select(  year = YEAR,
                              tripID = TRIP_ID,
                              eventID = FISHING_EVENT_ID,
                              lat = LATITUDE,
                              lon = LONGITUDE,
                              stratum = GROUPING_CODE,
                              majorArea = MAJOR_STAT_AREA_CODE,
                              minorArea = MINOR_STAT_AREA_CODE,
                              survey = SURVEY_DESC,
                              surveyID = SURVEY_ID,
                              survSeriesID = SURVEY_SERIES_ID,
                              stratArea = AREA_KM2,
                              doverDensity, doverCatch,
                              englishDensity, englishCatch,
                              rockDensity, rockCatch,
                              petraleDensity, petraleCatch,
                              atoothDensity, atoothCatch  ) %>%
              mutate( DERPAcatch =  doverCatch + 
                                    englishCatch + 
                                    rockCatch + 
                                    petraleCatch + 
                                    atoothCatch )

relativeBio <-  surveyData %>%
                group_by( survSeriesID, surveyID, year, stratum ) %>%
                dplyr::summarise( doverDensity = mean(doverDensity),
                                  englishDensity = mean(englishDensity),
                                  rockDensity = mean(rockDensity),
                                  petraleDensity = mean(petraleDensity),
                                  atoothDensity = mean(atoothDensity),
                                  area = mean(stratArea) ) %>%
                dplyr::summarize( relBioDover_Kt = sum( area * doverDensity),
                                  relBioEnglish_Kt = sum( area * englishDensity),
                                  relBioRock_Kt = sum( area * rockDensity),
                                  relBioPetrale_Kt = sum( area * petraleDensity),
                                  relBioAtooth_Kt = sum( area * atoothDensity) )



AEAproj <- CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m")
LLproj <- CRS("+proj=longlat +datum=WGS84")
UTMproj <- CRS("+proj=utm +zone=9 +datum=WGS84")

# Now create a raster brick object with 6 layers:
# one for each species, and one for the total
catchRaster   <- catchToRaster( saveFile = "surveyCatchRaster.RData" )
relBioRaster  <- densityToRaster( saveFile = "relBioRaster.RData", yrRange = c( 1984, 2016 ), dummRes = 50 )
names(catchRaster) <- c(specNames,"DERPA Total")
names(relBioRaster) <- c(specNames,"DERPA Total")



# Plot rasters with the map and management areas
plotCatchMap( catchRaster, saveFile = "surveyCatchHeatmap.png",
              units = "t", scale = 1000, quant = "Catch", leg.cex = 1 )
plotCatchMap( relBioRaster, saveFile = "relBioHeatmap.png",
              units = "kg", scale = 1, quant = "Tr. Biomass",
              width = 10, height = 5, axisTicks = F, leg.cex = 1 )

plotCatchMap( relBioRaster[[1]], label = FALSE, saveFile = "relBioHeatmapDover.png",
              units = "kg", scale = 1, quant = "B_trawl", multiPanel = F,
              width = 7, height = 7, axisTicks = F, leg.cex = 1.5, labLL = T )


# Read in commercial catch data
commCatch <- read.csv( "catch_by_maj.csv", header = T, stringsAsFactors = F)

# # Mung the commercial catch data
# catchAreaFleetSpecies <-  commCatch %>%
#                           group_by( YEAR, MAJOR_STAT_AREA_CODE, 
#                                     GEAR, SPECIES_CODE ) %>%
#                           summarise(  landedWt = sum( LANDED_KG ) / 1e6,
#                                       discardWt = sum( DISCARDED_KG ) / 1e6,
#                                       discardPc = sum( DISCARDED_PCS ),
#                                       specName = unique( SPECIES_COMMON_NAME ) ) %>%
#                           ungroup() %>%
#                           dplyr::select(  year = YEAR, 
#                                           majorStatArea = MAJOR_STAT_AREA_CODE,
#                                           gearType = GEAR,
#                                           specName,
#                                           specCode = SPECIES_CODE,
#                                           landedWt,
#                                           discardWt,
#                                           discardPc )

# catchSpecies <- catchAreaFleetSpecies %>%
#                 group_by(specName, year) %>%
#                 summarise( katch = sum(landedWt + discardWt) )

# plotCatch( df = catchAreaFleetSpecies )

