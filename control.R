# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# ageStructuredControl.R
# 
# Control script for development of the age-structured multistock
# OM estimation model. 
# 
# Model features:
#   - Multistock and multispecies (DERPA)
#   - Sex structured
#   - Multi-level RH priors on:
#       - Growth (vonB)
#       - Fishing mortality (correlation in REs if estimated)
#       - Natural mortality (M0 prior and correlated deviations)
#       - Selectivity
#       - Catchability for survey biomass - comm in kg/hr
#       - S-R Steepness
#   - Length composition observation model
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
#       - Hecate Strait Assemblage (1984 - 2002)
#       - Commercial CPUE: Modeled with tv q in 2 blocks:
#           - 1. Pre ASOP (1954 - 1995)
#           - 2. ASOP (1996 +)
#   - Biological observations from surveys on length and age, for 
#     integrated growth model
#   - Biological data from commercial fleets - length/age/sex
#   - Fixed maturity and weight-at-length values, estimated
#     separately from data
# 
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>


source("loadPackages.R")
source("DERPAfuns.R")
source("SCALfuns.R")
source("SCALtools.R")
source("mseRtools.r")
source("batchTools.R")
source("plots.R")
source("refPts.R")
# Survey ids for plotting/legends/array dims
loadStockSpecNameLists()

# Make Outputs directory if it doesn't exists
if(!dir.exists("Outputs"))
  dir.create("Outputs")

# compile and load scalRH objective function.
compile ("hierSCAL.cpp")


