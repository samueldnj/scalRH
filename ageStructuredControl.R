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
# Survey ids for plotting/legends/array dims
loadStockSpecNameLists()

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

# Go through bio data to get length bin midpoints and breaks
# Need to decide on
# Plus groups for the age structure for each species
#   Dover: 30
#   English: 15
#   Rock: 20
#   Petrale: 20
#   Arrowtooth: 20
nS <- 5
nP <- 3
nF <- 7
nT <- length(fYear:lYear)
# 2. maxLen for the length/age probability matrices

load("./Data/surveyBio.RData")
load("./Data/commCPUE.RData")

I_spft <- makeIndexArray()

fleetIDs <- dimnames(I_spft)$fleets

catchData <- read.csv(  "./Data/catch_by_maj.csv", header = TRUE,
                        stringsAsFactors = FALSE )

catchDiscArrays <- makeCatchDiscArrays(fleetIDs = fleetIDs)

load("./Data/lenAge.RData")

ALK_spalx <- makeALFreq()

load("./Data/wtLen.RData")
load("./Data/ageComps.RData")

plusA_s <- c( Dover = 25, 
              English = 15, 
              Rock = 20, 
              Petrale = 20, 
              Arrowtooth = 20 )

plusL_s <- c( Dover = 65, 
              English = 55, 
              Rock = 55, 
              Petrale = 60, 
              Arrowtooth = 75 )

age_aspft <- makeCompsArray(  compList = ageComps,
                              plusGroups = plusA_s,
                              collapseComm = FALSE,
                              fleetIDs = fleetIDs,
                              combineSex = TRUE,
                              years = fYear:lYear,
                              xName = "ages" )

load("./Data/lenComps.RData")

len_lspft <- makeCompsArray(  compList = lenComps,
                              plusGroups = plusL_s,
                              collapseComm = FALSE,
                              fleetIDs = fleetIDs,
                              combineSex = TRUE,
                              years = fYear:lYear,
                              xName = "ages" )


load("./Data/matOgives.RData" )


# 3. sex structured?
# To start with, no.


# 3. Create data, parameter and map lists for AM
#   a. Data is pretty self explanatory, collect 2a here
dat <- list(  I_spft = I_spft,
              C_spft = catchDiscArrays$C_spft,
              D_spft = catchDiscArrays$D_spft,
              # ALK_spal = ALK_spalx,
              age_aspft = age_aspft,
              len_lspft = len_lspft,
              type_f = as.integer(rep(0,nF)),
              A_s = as.integer(plusA_s),
              L_s = as.integer(plusL_s),
              lenD_s = rep(0,nS),
              swRinit_sp = matrix(0,nrow = nS, ncol = nP) )
#   b. Parameter list contains initial values for all input pars, 
#       whether estimated or not - 2b stuff will go here
par <- list(  lnB0_sp = array(5,dim =c(nS,nP)),
              logitSteep_sp = array(0,dim =c(nS,nP)) ,
              lnM_sp = array(-2,dim =c(nS,nP)),
              lnRinit_sp = array(1,dim =c(nS,nP)),
              lnLinf_sp = array(4,dim =c(nS,nP)),
              lnvonK_sp = array(-2,dim =c(nS,nP)),
              lnL1_sp = array(2,dim =c(nS,nP)),
              LWa_s = c(6.1e-6,6.1e-6,6.1e-6,6.1e-6,6.1e-6),
              LWb_s = c(3.19,3.19,3.19,3.19,3.19),
              xMat50_s = c(26.2,26.2,26.2,26.2,26.2),
              xMat95_s = c(48,48,48,48,48),
              lnq_spf = array(-1,dim =c(nS,nP,nF)),
              lntau_spf = array(-1,dim =c(nS,nP,nF)),
              lnlenSel50_sf = array(3.4,dim =c(nS,nF)),
              lnlenSel95_sf = array(4.5,dim =c(nS,nF)),
              lnF_spft  = array(-4,dim =c(nS,nP,nF,nT)),
              lntauC_f  = rep(-3,nF),
              lntauD_f  = rep(-3,nF),
              muLinf_s  = c(70,70,70,70,70),
              sigmaLinf_s = rep(.2,nS),
              muvonK_s  = rep(.2,nS),
              sigmavonK_s = rep(.2,nS),
              muL1_s  = rep(25,nS),
              sigmaL1_s = rep(.2,nS),
              lnsigmaL_s = rep(2.5,nS),
              mulenSel50_f  = rep(3.4,nF),
              mulenSel95_f  = rep(4.2,nF),
              sigmalenSel50_f = rep(.2,nF),
              sigmalenSel95_f = rep(.2,nF),
              lnqbar_fs = array(-2,dim =c(nF,nS)),
              lntauq_fs = array(-2,dim =c(nF,nS)),
              lnqbar_f  = rep(-2,nF),
              lntauq_f  = rep(-2,nF),
              logitSteep_s  = rep(0,nS),
              lnsigmaSteep_s  = rep(-2,nS),
              logitSteep  = 0,
              lnsigmaSteep  = -2,
              lnM_s = rep(-2,nS),
              lnsigmaM_s  = rep( -2, nS ),
              ln_muM  = -2,
              lnsigmaM  = -2,
              omegaR_spt = array(0, dim = c(nS,nP,nT-1)) ,
              omegaRinit_spa = array(0, dim = c(nS,nP,max(plusA_s))) ,
              lnsigmaR_sp = array(0, dim = c(nS,nP))  ,
              logitRCorr_chol = rep(0, nS * nP),
              logitRgamma_sp  = array(0, dim = c(nS,nP))  )

#   c. The map list collects parameters that are treated as identical
#       by the estimator, and also allows us to turn off parameters
#       (negative phase in ADMB)
map <- list(  # lnB0_sp = factor(array(NA,dim =c(nS,nP))),
              logitSteep_sp = factor(array(NA,dim =c(nS,nP)) ),
              lnM_sp = factor(array(NA,dim =c(nS,nP))),
              lnRinit_sp = factor(array(NA,dim =c(nS,nP))),
              lnLinf_sp = factor(array(NA,dim =c(nS,nP))),
              lnvonK_sp = factor(array(NA,dim =c(nS,nP))),
              lnL1_sp = factor(array(NA,dim =c(nS,nP))),
              LWa_s = factor(rep(NA,nS)),
              LWb_s = factor(rep(NA,nS)),
              xMat50_s = factor(rep(NA,nS)),
              xMat95_s = factor(rep(NA,nS)),
              lnq_spf = factor(array(NA,dim =c(nS,nP,nF))),
              lntau_spf = factor(array(NA,dim =c(nS,nP,nF))),
              lnlenSel50_sf = factor(array(NA,dim =c(nS,nF))),
              lnlenSel95_sf = factor(array(NA,dim =c(nS,nF))),
              lnF_spft  = factor(array(NA,dim =c(nS,nP,nF,nT))),
              lntauC_f  = factor(rep(NA,nF)),
              lntauD_f  = factor(rep(NA,nF)),
              muLinf_s  = factor(rep(NA,nS)),
              sigmaLinf_s = factor(rep(NA,nS)),
              muvonK_s  = factor(rep(NA,nS)),
              sigmavonK_s = factor(rep(NA,nS)),
              muL1_s  = factor(rep(NA,nS)),
              sigmaL1_s = factor(rep(NA,nS)),
              lnsigmaL_s = factor(rep(NA,nS)),
              mulenSel50_f  = factor(rep(NA,nF)),
              mulenSel95_f  = factor(rep(NA,nF),),
              sigmalenSel50_f = factor(rep(NA,nF),),
              sigmalenSel95_f = factor(rep(NA,nF),),
              lnqbar_fs = factor(array(NA,dim =c(nF,nS))),
              lntauq_fs = factor(array(NA,dim =c(nF,nS))),
              lnqbar_f  = factor(rep(NA,nF)),
              lntauq_f  = factor(rep(NA,nF)),
              logitSteep_s  = factor(rep(NA,nS)),
              lnsigmaSteep_s  = factor(rep(-NA,nS)),
              logitSteep  = factor(NA),
              lnsigmaSteep  = factor(NA),
              lnM_s = factor(rep(NA,nS)),
              lnsigmaM_s  = factor(rep(NA, nS ) ),
              ln_muM  = factor(NA),
              lnsigmaM  = factor(NA),
              omegaR_spt = factor(array(NA, dim = c(nS,nP,nT-1))),
              omegaRinit_spa = factor(array(NA, dim = c(nS,nP,max(plusA_s)))),
              lnsigmaR_sp = factor(array(NA, dim = c(nS,nP)) ),
              logitRCorr_chol = factor(rep(NA, nS * nP)),
              logitRgamma_sp  = factor(array(NA, dim = c(nS,nP))  )) 


# Make the AD function
obj <- MakeADFun (  dat = dat, parameters = par, map = map,
                    silent = FALSE )

# Set max no of function/gradient evaluations
ctrl = list ( eval.max = 1000, iter.max = 1000 )

# # optimise the model
# fit <- try( nlminb (  start = obj$par,
#                       objective = obj$fn,
#                       gradient = obj$gr,
#                       control = ctrl,
#                       lower = -Inf,
#                       upper= Inf ) )
