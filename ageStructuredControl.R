# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# ageStructuredControl.R
# 
# Control script for development of the age-structured multistock
# OM estimation model. 
# 
# Model features:
#   - Multistock and multispecies (DERPA)
#   - Single sex? Or split sex?
#   - Multi-level RH priors on:
#       - Growth (vonB )
#       - Fishing mortality (correlation in REs if estimated)
#       - Natural mortality (M0 prior and correlated deviations)
#       - Selectivity
#       - Catchability for survey biomass - comm in kg/hr
#       - S-R Steepness
#   - Length composition observation model (used in years where ages unavailable)
#   - Age composition observation model (used when ages available)
#   - Integrated growth model to predict length dist
#   - Discarding - takes a grading length
# 
# Data:
# 
# Data from the 2017 DERPA data request are provided to the assessment
# model. These include:
#   - yearly catch and discards by species and management 
#     area (1954 - 2016)
#   - survey data by year, species and management area for
#       - Synoptic survey (all 4 legs) (2003 - 2016)
#       - Hecate Strait Assemblage (1984 - 2003)
#       - Fine mesh surveys - not really used
#       - Commercial CPUE: Modeled with tv q in 3 blocks:
#           - 1. Pre ASOP (1954 - 1996)
#           - 2. ASOP (1997-2005)
#           - 3. Integration (2006+)
#   - Biological observations from surveys on length and age, for 
#     integrated growth model
#   - Biological data from commercial fleets - length/age/sex
#   - Fixed maturity and weight-at-length values, estimated
#     separately from data
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

source("DERPAfuns.R")

# compile and load scalRH objective function.
compile ("scalRH.cpp")
dyn.load(dynlib("scalRH"))


# 1. Read control file that determines
#   a. Level of detail - multistock or multispecies?
#   b. What priors will we use?
fYear <- 1954
lYear <- 2018

# 2. Read in and groom data
#   a. Organise into the management units and species groups, 
#       lots of arrays here 
#   b. externally estimate any parameters that need it

stocksSurvDover <- list(  HG = c(2,3,16),
                          QCS = c(1),
                          WCVI = c(4) )

stocksCommDover <- list(  HG = c(7,8,9),
                          QCS = c(5,6),
                          WCVI = c(3,4) )

#stocksSurvAtooth <- list( CW = c(1,2,3,4,16))
#stocksCommAtooth <- list( CW = 3:9 )

stocksSurvAtooth <- stocksSurvDover
stocksCommAtooth <- stocksCommDover

relBioDover  <- makeRelBioStocks( years = c(fYear,lYear), spec = "dover", collapseSyn = FALSE, 
                                  stocks = stocksSurvDover )
catchDover   <- makeStockCatch( years = c(fYear,lYear), spec = "dover", 
                                stocks = stocksCommDover )

relBioAtooth  <- makeRelBioStocks( years = c(fYear,lYear), spec = "atooth", collapseSyn = FALSE, 
                                  stocks = stocksSurvAtooth )
catchAtooth   <- makeStockCatch( years = c(fYear,lYear), spec = "atooth", 
                                stocks = stocksCommAtooth )

bioDataAtooth <- makeBioDataStocks( years = c(fYear,lYear),
                                    spec = "atooth",
                                    surveyStocks = stocksSurvAtooth,
                                    commStocks = stocksCommAtooth )



# Number of gear types (fleets + each synoptic survey leg) - might restrict to 6
# initially to remove LL and trap - setting it up this way allows us to use 
# standardised fishery dep. indices later
nG <- 6
nP <- length(stocksSurvDover) + length(stocksSurvDover)
nS <- 2
s_p <- c(1,1,1,2,2,2)
type_f <- c(0,0,0,0,0,1)
A_s <- c(60,25)
nL_s <- c(5,5)
swRinit_p <- rep(1,nP)
lenD_s <- c(30,35)
nT <- lYear - fYear + 1

I_pft <- array( -1, dim = c(nP,nG,nT))
for( g in 1:5 )
{
  I_pft[1:3,g,] <- relBioDover$relBio.arr[g,,]
  I_pft[4:6,g,] <- relBioAtooth$relBio.arr[g,,]  
}


C_pft <- array( 0, dim = c(nP,nG,nT) )
C_pft[1:3,6,] <- catchDover$catch.arr
C_pft[4:6,6,] <- catchAtooth$catch.arr

# Go through bio data to get length bin midpoints and breaks





# 3. Create data, paramater and map lists for AM
#   a. Data is pretty self explanatory, collect 2a here
#   b. Parameter list contains initial values for all input pars, 
#       whether estimated or not - 2b stuff will go here
#   c. The map list collects parameters that are treated as identical
#       by the estimator, and also allows us to turn off parameters
#       (negative phase in ADMB)       
