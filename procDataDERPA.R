# -----------------------------------------------------------------------------
#
# procDataDERPA.R
#
# Reads in biological and index data from the ./Data/ folders inside the working 
# directory. Processes raw data into data.frame and array objects for
# easier plotting and passing to a TMB model. Will also plot data, but due
# to time and memory requirements this code is intended to be run only when
# updates to how the data is structured are required.
# 
# -----------------------------------------------------------------------------

rm(list = ls())

library( dplyr )
# library( ggmap )
library( maptools )
library( rgdal )
library( PBSmapping ) 
library( RColorBrewer )
# library( RgoogleMaps)
library( sp )
library( raster )
library( MASS )
# library( psyphy )
library( boot )
library( RCurl )
library( scales )
# library( rgeos )



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

# Read in strata areas
stratData <- read.csv( "Data/derpa_strata.csv", header=T )
stratData <-  stratData %>%
              dplyr::select( GROUPING_CODE, AREA_KM2, MIN_DEPTH, MAX_DEPTH )

# Survey ids for plotting/legends
surveyIDs <-  c(  QCSSyn = 1, 
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

# loadCRS()

# surveys     <- c("HS", "QCS", "WCHG", "WCVI")
# shapeFiles  <- paste(surveys,"_Synoptic_Survey_Active_Blocks",sep = "")
# shapePath   <- file.path("./Data/ShapeFiles/SynSurveyBlocks")

# Load grids in a list
# grids <- lapply(X = shapeFiles, FUN = openShapeFile, path = shapePath)
# names(grids) <- surveys

# These have been coerced to UTM so the grids are nice
# and rectangular.

# Make shrunken biomass indices by removing small fish
# source("removeSmallFish.R")

# Plots that we want to make - and may not
# need to reinvent the code for...
# 1. Catch and indices - how to plot comm CPUE?
# relBioList_Survey <- lapply(  X = survSpecNames,
#                               FUN = makeRelBioStocks,
#                               years = c(fYear,lYear),
#                               stocks = stocksSurvey,
#                               survIDs = surveyIDs,
#                               stratArea = stratData,
#                               grids = grids )
# names(relBioList_Survey) <- names(survSpecNames)
# save(relBioList_Survey, file = "./Data/surveyBio.RData")

commCPUEList <- lapply( X = commSpecNames,
                        FUN = readCommCPUE,
                        stocks = stocksCommCPUE )
names(commCPUEList) <- names(commSpecNames)
save(commCPUEList, file = "./Data/Proc/commCPUE.RData")


# Read in bio data, join survey density by tripID
# to add data for year and location
bioData <- lapply(  X = survSpecNames,
                    FUN = readProcBioData )
names(bioData) <- names(commSpecNames)
save(bioData, file = "./Data/Proc/bioData.RData")

# # 2. Length at age plots - stock and sex - spit out age-length freq array
# ALfreq <- lapply( X = bioData, FUN = makeALFreq_FleetYear )
# names(ALfreq) <- names(commSpecNames)
# # Save data out
# save(ALfreq, file = "./Data/Proc/ALfreq.RData")

# # 3. Length/wt plots - stock and sex
# wtLen <- lapply( X = bioData, FUN = makeWtLen, stocks = names(stocksCommBio))
# names(wtLen) <- names(commSpecNames)
# # Save data out
# save(wtLen, file = "./Data/Proc/wtLen.RData")
# wtLen.df <- makeWtLenDF(wtLen)
# write.csv(wtLen.df, file = "./Data/Proc/wtLen.csv")
# # plotWtLen(save = TRUE)

# # 4. Catch and discards - Species and area
# catchData <- read.csv(  "./Data/Raw/catch_by_maj.csv", header = TRUE,
#                         stringsAsFactors = FALSE )
# plotCatch(save = TRUE)
# # Get survey removals
# surveyCatch <- lapply( X = survSpecNames, FUN = makeSurveyCatchStocks )
# names(surveyCatch) <- names(survSpecNames)
# save( surveyCatch, file = "./Data/Proc/surveyCatch.RData" )

# 5a. Age compositions by fleet, stock, and species - spit out comp data array
ageComps <- lapply( X = bioData, FUN = makeAgeComps )
save(ageComps, file = "./Data/Proc/ageComps.RData")

# 5b. length compositions by fleet, stock, and species - spit out comp data array
lenComps <- lapply( X = bioData, FUN = makeLenComps )
save(lenComps, file = "./Data/Proc/lenComps.RData")



# 6. Maturity at age by stock and species
matOgives <- lapply( X = bioData, FUN = makeSpecMat )
save(matOgives, file = "./Data/Proc/matOgives.RData" )
# matAgeLen.df <- makeMatDF(matOgives)
# write.csv(matAgeLen.df, file = "./Data/Proc/matAgeLen.csv")


