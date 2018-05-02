# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# ageStructuredControl.R
# 
# Control script for development of the age-structured multistock
# OM estimation model. 
# 
# Model features:
#   - Multistock and multispecies (DERPA)
#       - Will have interaction/spawning mtx for individual
#         stocks within a species
#       - emigration/immigration? Sounds like OM world to me
#       - Species made up of a single stock have only speciew level 
#         priors
#   - Multi-level RH priors on:
#       - Growth (vonB or F-W undecided)
#       - Fishing mortality (correlation in REs if estimated)
#       - Natural mortality (M0 prior and correlated deviations)
#       - Selectivity
#       - Catchability
#       - S-R Steepness
#   - Length based observation model
#   - Integrated growth model to predict length dist
#   - Discarding
# 
# Data:
# 
# Data from the 2017 DERPA data request are provided to the assessment
# model. These include:
#   - yearly catch and discards by species, fleet and management 
#     area (1954 - 2016)
#   - survey data by year, species and management area for
#       - Synoptic survey (all 4 legs) (2003 - 2016)
#       - Hecate Strait Assemblage (1984 - 2003)
#       - Fine mesh surveys - not really used
#   - Biological observations from above surveys, including
#       - maturity
#       - length
#       - age
#       - weight
# 
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

# Load packages
library( "coda" )
library( "dplyr" )
library( "TMB" )
library( "raster" )
library( "grid" )
library( "RColorBrewer" )
library( "HapEstXXR" )
library( "parallel" )
library( "stringr" )

# compile and load msProd objective function.
compile ("scalRH.cpp")
dyn.load(dynlib("scalRH"))

# 1. Read control file that determines
#   a. Level of detail - multistock or multispecies?
#   b. What priors will we use?

fYear <- 1984
lYear <- 2016

# 2. Read in and groom data
#   a. Organise into the management units and species groups, 
#       lots of arrays here 
#   b. externally estimate any parameters that need it

stocksSurvDover <- list(  HG = c(2,3,16),
                          QCS = c(1),
                          WCVI = c(4) )

stocksSurvAtooth <- list( CW = c(1,2,3,4,16))
stocksCommAtooth <- list( CW = 3:9 )

stocksCommDover <- list(  HG = c(7,8,9),
                          QCS = c(5,6),
                          WCVI = c(3,4) )

# Number of gear types (fleets + each synoptic survey leg) - might restrict to 6
# initially to remove LL and trap - setting it up this way allows us to use 
# standardised fishery dep. indices later
nG <- 8
nP <- length(stocksSurvDover) + length(stocksSurvAtooth)
nS <- 2
s_p <- c(1,1,1,2)
type_f <- c(0,0,0,0,0,1,1,1)
A_s <- c(60,25)
nL_s <- c(5,5)
swRinit_p <- rep(1,nP)
lenD_s <- c(0,0)


# Go through bio data to get length bin midpoints and breaks

relBioDover  <- makeRelBioStocks( years = c(fYear,lYear), spec = "dover", collapseSyn = FALSE, 
                                  stocks = stocksSurvDover )
catchDover   <- makeStockCatch( years = c(fYear,lYear), spec = "dover", 
                                stocks = stocksCommDover )

relBioAtooth  <- makeRelBioStocks( years = c(iYear,2016), spec = "atooth", collapseSyn = FALSE, 
                                  stocks = stocksSurvAtooth )
catchAtooth   <- makeStockCatch( years = c(iYear,2016), spec = "atooth", 
                                stocks = stocksCommAtooth )





# 3. Create data, paramater and map lists for AM
#   a. Data is pretty self explanatory, collect 2a here
#   b. Parameter list contains initial values for all input pars, 
#       whether estimated or not - 2b stuff will go here
#   c. The map list collects parameters that are treated as identical
#       by the estimator, and also allows us to turn off parameters
#       (negative phase in ADMB)       
