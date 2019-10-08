# -----------------------------------------------------------------------------
#
# plotDERPAdata.R
#
# Reads in processed biological and index data from the ./Data/
# and folder inside the working directory, and plots standard
# data plots for use inside supplemental materials
# 
# -----------------------------------------------------------------------------

rm(list = ls())

library( dplyr )
# library( ggmap )
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

# # Read in management areas shape file
# readShapeSpatial("./Data/ShapeFiles/MajorMinorSQL_geo.shp")
# mgmtAreas <- readOGR(dsn = "./Data/ShapeFiles/", layer = "MajorMinorSQL_geo")
# proj4string(mgmtAreas)
# mgmtAreas <- spTransform (mgmtAreas, CRS("+proj=longlat +datum=WGS84"))

# # Read in strata areas
# stratData <- read.csv( "Data/derpa_strata.csv", header=T )
# stratData <-  stratData %>%
#               dplyr::select( GROUPING_CODE, AREA_KM2 )
# Survey ids for plotting/legends
loadStockSpecNameLists()


# Plots that we want to make - and may not
# need to reinvent the code for...
# 1. Indices
load("./Data/Proc/surveyBio.RData")
load("./Data/Proc/commCPUE.RData")

# Plot stock indices
# plotIndices(save = TRUE)

# Create a version that has the catch bars in it - maybe scale
# catch by geometric mean also

# 2. Length at age plots - stock and sex - spit out age-length freq array
# load("./Data/Proc/lenAge.RData")
# plotLenAge(save = TRUE)

# 3. Length/wt plots - stock and sex
load("./Data/Proc/wtLen.RData")
plotWtLen(save = TRUE)

# 4. Catch and discards - Species and area
catchData <- read.csv(  "./Data/Proc/catch_by_maj.csv", header = TRUE,
                        stringsAsFactors = FALSE )
plotCatch(save = TRUE)

# 5a. Age compositions by fleet, stock, and species - spit out comp data array
load("./Data/Proc/ageComps.RData")
plotComps(  comps = ageComps, save = TRUE,
            prefix = "age", saveDir = "Outputs/ageComps" )

# 5b. length compositions by fleet, stock, and species - spit out comp data array
load("./Data/Proc/lenComps.RData")
plotComps(  comps = lenComps, save = TRUE,
            prefix = "len", saveDir = "Outputs/lenComps" )

# 7. Maturity at age and length by stock and species
load("./Data/Proc/matOgives.RData" )
plotMatOgives( type = "age", save = TRUE )
plotMatOgives( type = "length", save = TRUE )


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

