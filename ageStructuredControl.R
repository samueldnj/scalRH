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
compile ("hierSCAL.cpp")
dyn.load(dynlib("hierSCAL"))


# 1. Read control file that determines
#   a. Level of detail - multistock or multispecies?
#   b. What priors will we use?
#   c. sex structured or not?
fYear <- 1954
lYear <- 2018

# 2. Read in and groom data - grooming is taken care
# of by procDataDERPA.R, and groomed objects are saved in 
# RData files in ./Data/ 
load("./Data/surveyBio.RData")
load("./Data/commCPUE.RData")
load("./Data/lenAge.RData")
load("./Data/wtLen.RData")
load("./Data/ageComps.RData")
load("./Data/lenComps.RData")
load("./Data/matOgives.RData" )

# Go through bio data to get length bin midpoints and breaks
# Need to decide on
# 1. plus groups for the age structure for each species
# Dover: 30
# English: 15
# Rock: 20
# Petrale: 20
# Arrowtooth: 20
# 2. maxLen for the length/age probability matrices

# 3. sex structured?
# To start with, no.



# 3. Create data, parameter and map lists for AM
#   a. Data is pretty self explanatory, collect 2a here
#   b. Parameter list contains initial values for all input pars, 
#       whether estimated or not - 2b stuff will go here
#   c. The map list collects parameters that are treated as identical
#       by the estimator, and also allows us to turn off parameters
#       (negative phase in ADMB)       
