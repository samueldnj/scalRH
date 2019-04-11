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

library( dplyr )
library( ggmap )
library( maptools )
library( rgdal )
library( PBSmapping ) 
library( RColorBrewer )
library( RgoogleMaps)
library( sp )
library( raster )
library( MASS )
library( psyphy )
library( boot )
library( RCurl )
library( scales )

# Load shape file data
data(nepacLL)

source("mseRtools.r")
source("DERPAfuns.R")

# Need to pick a first year:
# Do we do 1956 (first year of comm CPUE)
# or 1954 (first year of catch)?

# Make most complete time series we can, so start
# at earliest year, which is from catch.

fYear <- 1954
lYear <- 2018

# Read in management areas shape file
readShapeSpatial("./Data/ShapeFiles/MajorMinorSQL_geo.shp")
mgmtAreas <- readOGR(dsn = "./Data/ShapeFiles/", layer = "MajorMinorSQL_geo")
proj4string(mgmtAreas)
mgmtAreas <- spTransform (mgmtAreas, CRS("+proj=longlat +datum=WGS84"))

# Read in strata areas
stratData <- read.csv( "Data/derpa_strata.csv", header=T )
stratData <-  stratData %>%
              dplyr::select( GROUPING_CODE, AREA_KM2 )
# Survey ids for plotting/legends
surveyIDs <-  c(  QCSyn = 1, 
                  HSAss = 2, 
                  HSSyn = 3, 
                  WCVISyn = 4,
                  WCHGSyn = 16 )

# Species codes
specCodes <- list(  "dover" = 626,
                    "english" = 628,
                    "rock" = 621,
                    "petrale" = 607,
                    "atooth" = 602 )

# Stock IDs for grouping data
stocksSurvey <- list( HSHG = c(2,3,16),
                      QCS = c(1),
                      WCVI = c(4) )
stocksCommBio <- list(  HSHG = c(7,8,9),
                        QCS = c(5,6),
                        WCVI = c(3,4) )
stocksCommCPUE  <- list(  HSHG = "5CDE",
                          QCS = "5AB",
                          WCVI = "3CD" )

# Species names for reading data
survSpecNames <- c( Dover = "dover",
                    English = "english",
                    Rock = "srock",
                    Petrale = "petrale",
                    Arrowtooth = "atooth" )

commSpecNames <- c( Dover = "dover-sole",
                    English = "english-sole",
                    Rock = "southern-rock-sole",
                    Petrale = "petrale-sole",
                    Arrowtooth = "arrowtooth-flounder" )



# Plots that we want to make - and may not
# need to reinvent the code for...
# 1. Catch and indices - how to plot comm CPUE?
relBioList_Survey <- lapply(  X = survSpecNames,
                              FUN = makeRelBioStocks,
                              years = c(fYear,lYear),
                              stocks = stocksSurvey,
                              survIDs = surveyIDs,
                              stratArea = stratData )
names(relBioList_Survey) <- names(survSpecNames)
save(relBioList_Survey, file = "./Data/Proc/surveyBio.RData")

commCPUEList <- lapply( X = commSpecNames,
                        FUN = readCommCPUE,
                        stocks = stocksCommCPUE )
names(commCPUEList) <- names(commSpecNames)
save(commCPUEList, file = "./Data/Proc/commCPUE.RData")

# Plot stock indices
plotIndices(save = TRUE)

# Create a version that has the catch bars in it - maybe scale
# catch by geometric mean also

# Read in bio data, join survey density by tripID
# to add data for year and location
bioData <- lapply(  X = survSpecNames,
                    FUN = readBioData )
names(bioData) <- names(commSpecNames)

# 2. Length at age plots - stock and sex - spit out age-length freq array
lenAge <- lapply( X = bioData, FUN = makeLenAge, stocks = names(stocksCommBio) )
names(lenAge) <- names(commSpecNames)
# Save data out
save(lenAge, file = "./Data/Proc/lenAge.RData")
plotLenAge(save = TRUE)

# 3. Length/wt plots - stock and sex
wtLen <- lapply( X = bioData, FUN = makeWtLen, stocks = names(stocksCommBio))
names(wtLen) <- names(commSpecNames)
# Save data out
save(wtLen, file = "./Data/Proc/wtLen.RData")
plotWtLen(save = TRUE)

# 4. Catch and discards - Species and area
catchData <- read.csv(  "./Data/Raw/catch_by_maj.csv", header = TRUE,
                        stringsAsFactors = FALSE )
plotCatch(save = TRUE)

# 5a. Age compositions by fleet, stock, and species - spit out comp data array
ageComps <- lapply( X = bioData, FUN = makeAgeComps )
save(ageComps, file = "./Data/Proc/ageComps.RData")
# plotAgeComps()

# 5b. length compositions by fleet, stock, and species - spit out comp data array


# 7. Maturity at age by stock and species



# dover <- dataPlots("dover", save = F, maxMatAge = 5 )
# english <- dataPlots("english", save = F, maxMatAge = 5 )
# srock <- dataPlots("srock", save = F, maxMatAge = 10)
# petrale <- dataPlots("petrale", save = F, maxMatAge = 15)
# atooth <- dataPlots("atooth", save = F, maxMatAge = 10)


# Now create a raster brick object with 6 layers:
# one for each species, and one for the total
# catchRaster   <- catchToRaster( saveFile = "surveyCatchRaster.RData" )
# relBioRaster  <- densityToRaster( saveFile = "relBioRaster.RData", yrRange = c( 1984, 2016 ), dummRes = 50 )
# names(catchRaster) <- c(specNames,"DERPA Total")
# names(relBioRaster) <- c(specNames,"DERPA Total")



# # Plot rasters with the map and management areas
# plotCatchMap( catchRaster, saveFile = "surveyCatchHeatmap.png",
#               units = "t", scale = 1000, quant = "Catch", leg.cex = 1 )
# plotCatchMap( relBioRaster, saveFile = "relBioHeatmap.png",
#               units = "kg", scale = 1, quant = "Tr. Biomass",
#               width = 10, height = 5, axisTicks = F, leg.cex = 1 )

# plotCatchMap( relBioRaster[[1]], label = FALSE, saveFile = "relBioHeatmapDover.png",
#               units = "kg", scale = 1, quant = "B_trawl", multiPanel = F,
#               width = 7, height = 7, axisTicks = F, leg.cex = 1.5, labLL = T )


# Read in commercial catch data
# commCatch <- read.csv( "catch_by_maj.csv", header = T, stringsAsFactors = F)

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

