# Functions for assessCA TMB model

# fitAssessCA()
# Wrapper function to fit the assessCA model under a data scenario 
# and model hypothesis, both defined in the ctlFile. If fit is succesful
# then model outputs are saved to a folder in the ./Outputs/fits/ directory
# inputs:   ctlFile=character with the name/path of the control file
#           folder=optional character name of output folder 
#                   ie saves to ./Outputs/fits/<folder>
# ouputs:   NULL
# usage:    from the console to run the procedure
fitHierSCAL <- function ( ctlFile = "fitCtlFile.txt", folder=NULL, quiet=TRUE )
{ 
  # read in control file
  controlList <- .readParFile ( ctlFile )
  controlList <- .createList  ( controlList )

  # Run simEst Procedure
  reports <- .runHierSCAL( obj = controlList )
  # save output to project folder
  # First, if a folder name isn't nominated, create a default sim folder
  if ( is.null(folder) )
  {
    stamp <- paste( format(Sys.time(),format="%d%m%Y%H%M%S" ),sep="" )
    folder <- paste ( "fit_",stamp, sep = "" )
  } else folder <- paste ("fit_", folder, sep = "")
  # Now paste together the path to the folder and create it
  path <- file.path (getwd(),"Outputs","fits",folder)
  dir.create ( path )
  cat( "\nMSG (saveSim) Created assessment fit folder ",folder,"in ./Outputs/fits/.\n" )
  # Save blob
  save(reports,file = file.path(path,paste(folder,".RData",sep="")))
  # Save plots

  ## Fill in plotting code later ##

  # Copy control file to sim folder for posterity
  file.copy(from=ctlFile,to=file.path(path,"fitCtlFile.txt"))
  # Done
}

# .runHierSCAL()
# Procedure to create data, parameter and map 
# lists, and fit the assessCA TMB model. There are
# some path specific features here that relate to YE
# assessment stuff, so take care if copying to other
# locations
.runHierSCAL <- function( obj = controlList )
{
  # Get data scenario and model hypothesis control lists
  dataCtl <- obj$data
  hypoCtl <- obj$hypo
  ctrlObj <- obj$ctrl

  # Get model dimensions
  nA      <- dataCtl$nA
  fYear   <- dataCtl$fYearData
  lYear   <- dataCtl$lYearData
  nT      <- lYear - fYear + 1
  years   <- fYear:lYear
  yrChar  <- as.character(years)

  # Now load the fleetIDs
  loadStockSpecNameLists()

  # Track used fleets
  allFleets     <- c( names(fleetsCommCatch), names(surveyIDs) )
  useFleets     <- c( dataCtl$commNames_g, dataCtl$survNames_g)
  useFleetsIdx  <- which( allFleets %in% useFleets )
  nF            <- length(useFleets)
  totF          <- length(allFleets)
  fTimeRange    <- dataCtl$fleetTimingRange

  # What stock are we assessing?
  stockID       <- dataCtl$stock

  # Load data from ./Data folder
  catchDataPath <- file.path("./Data",dataCtl$commCatchData )
  catchData <- read.csv(  catchDataPath, header = TRUE,
                          stringsAsFactors = FALSE )

  # Generate catch/disc arrays from using function from YEfuns.R
  catchDiscArrays <- makeCatchDiscArrays( data = catchData,
                                          stocks = stocksCommBio,
                                          fleets = fleetsCommCatch,
                                          years = fYear:lYear,
                                          surveys = surveyIDs )


  # Use historical reconstruction data
  C_pft <- catchDiscArrays$C_pftx[stockID,useFleets,,"recon"]

  # Load index data
  idxDataFiles <- dataCtl$idxFiles
  idxDataFilePaths <- file.path("./Data",idxDataFiles )

  # Load index data
  for( i in 1:length(idxDataFilePaths))
    load( idxDataFilePaths[i] )

  # Generate index data arrays using function from YEfuns.R
  I_ptf <- makeIndexArray(  trawlSurv = relBio_Trawl,
                            phma = phmaIdx,
                            iphc = iphcIdx_ysit,
                            years = fYear:lYear,
                            iphcThresh = dataCtl$iphcThreshold,
                            commIDs = names(fleetsCommCatch) )
  # Subset down to used fleets
  I_ptf <- I_ptf[stockID,,useFleets]

  # Now make age compositions
  ageDataPath <- file.path("./Data",dataCtl$ageData)
  load(ageDataPath)
  A_aptf <- makeCompsArray (  compList = ageComps,
                              plusGroups = c(nA),
                              fleets = fleetsCommCatch,
                              surveys = surveyIDs,
                              years = fYear:lYear,
                              combineSex = TRUE,
                              xName = "ages" )

  # Subset down to useFleets
  A_aptf <- A_aptf[,stockID,,useFleets]



  # Generate the model switch settings we're using
  fleetSwitches <- setFleetSwitches()
  
  calcIndex <- integer(length = nF)
  for( fIdx in 1:nF )
    if( any(I_ptf[,fIdx] > 0) )
      calcIndex[fIdx] <- 1

  # Initialised in a fished state (0 == unfished, 1 == fished)
  initUnfished    <- hypoCtl$initUnfished[stockID]
  yFirstRecDev  <- hypoCtl$yFirstRecDev[stockID]

  # Generate tFirstRecDev from yFirstRecDev and fYear
  tFirstRecDev <- yFirstRecDev - fYear
  if( tFirstRecDev < 1 ) tFirstRecDev <- 1

  # Generate the data list
  data <- list( I_tg = I_ptf,
                C_tg = t(C_pft),
                A_atg = A_aptf,
                survType_g = fleetSwitches$surveyType[useFleets],
                indexType_g = fleetSwitches$idxType[useFleets],
                calcIndex_g = calcIndex,
                selType_g = hypoCtl$selFun[1:nF],
                # Update the following line so age/len based sel can be changed
                # from control file
                selLen_g = rep(0,nF),
                fleetTiming = seq(fTimeRange[1], fTimeRange[2], length = nF),
                initCode = c(initUnfished),
                posPenFactor = c(dataCtl$posPenFactor),
                firstRecDev = tFirstRecDev )

  # Generate parameter list
  pars <- list( lnB0 = log(sum(C_pft)),
                logit_ySteepness = 0,
                lnM = log(hypoCtl$initMprior[1]),
                log_initN_mult = rep(0,nA),
                lnSelAlpha_g = rep(2.7,nF),
                lnSelBeta_g = rep(1,nF),
                lntauAge_g = rep(-1.6,nF),
                effSampleSize_g = rep(100,nF),
                recDevs_t = rep(0,nT-tFirstRecDev + 1 ),
                lnsigmaR = log(hypoCtl$sigmaR),
                omegaM_t = rep(0,nT-1),
                lnsigmaM = log(hypoCtl$sigmaM),
                obstau2IGa = rep(hypoCtl$tau2ObsIGa[1],nF),
                obstau2IGb = (hypoCtl$tau2ObsIGa[1]+1)*hypoCtl$tau2ObsPriorMode[useFleetsIdx],
                sig2RPrior = c(1,2),
                sig2MPrior = c(1,0.04),
                rSteepBetaPrior = c(11,9),
                initMPrior = hypoCtl$initMprior,
                mq = hypoCtl$mq[useFleetsIdx],
                sdq = hypoCtl$sdq[useFleetsIdx],
                # revise maturity ages
                aMat = c(15,25),
                Linf = c(67),
                L1 = c(27.63),
                vonK = c(0.05),
                lenWt = c(1.52e-5,3.05),
                mSelMode_g = hypoCtl$selModePriorMean[useFleetsIdx],
                selModeCV = hypoCtl$selModePriorCV )

  # Generate map list - this will have to be updated
  # so that it's sensitive to useFleetsIdx
  selAlphaMap <- c(1,NA,2,NA,rep(NA,5),3,4,5,NA)
  selBetaMap  <- c(11,NA,12,NA,rep(NA,5),13,14,15,NA)
  map <- list(  # lnB0 =factor(NA),
                # logit_ySteepness = factor(NA),
                lnM = factor(NA),
                # log_initN_mult = factor(rep(NA,nA)),
                lnSelAlpha_g = factor(selAlphaMap[useFleetsIdx]),
                lnSelBeta_g = factor(selBetaMap[useFleetsIdx]),
                recDevs_t = factor(rep(NA,nT-tFirstRecDev+1)),
                effSampleSize_g = factor(rep(NA,nF)),
                lnsigmaR = factor(NA),
                omegaM_t = factor(rep(NA,nT-1)),
                lnsigmaM = factor(NA),
                obstau2IGa = factor(rep(NA,nF)),
                obstau2IGb = factor(rep(NA,nF)),
                sig2RPrior = factor(c(NA,NA)),
                sig2MPrior = factor(c(NA,NA)),
                rSteepBetaPrior = factor(c(NA,NA)),
                initMPrior = factor(c(NA,NA)),
                aMat = factor(c(NA,NA)),
                Linf = factor(c(NA)),
                L1 = factor(c(NA)),
                vonK = factor(c(NA)),
                lenWt = factor(c(NA,NA)),
                mSelMode_g = factor(rep(NA,nF)),
                selModeCV  = factor(NA),
                mq = factor(rep(NA,nF)),
                sdq = factor(rep(NA,nF)) )  


  # Create a control list for the assessment model
  tmbCtrl <- list(  eval.max = ctrlObj$maxFunEval, 
                    iter.max = ctrlObj$maxIterations  )

  # Now create the AD fun object, assuming all fixed effects
  objFE <- MakeADFun( data = data,
                      parameters = pars,
                      map = map, 
                      random = NULL, 
                      silent = ctrlObj$quiet )

  # Now try fitting the model
  fitFE <- try( nlminb (  start     = objFE$par,
                          objective = objFE$fn,
                          gradient  = objFE$gr,
                          control   = tmbCtrl ) )

  # May hang here if there are problems fitting the model
  repFE <- objFE$report()
  sdrepFE <- sdreport(objFE)
  sdrepFE.full <- summary(sdrepFE)


  report <- list( repFE = repFE,
                  sdrepFE = sdrepFE,
                  sdrepFE.full = sdrepFE.full,
                  fYear = fYear, lYear = lYear,
                  gearLabs = useFleets )

  return(report)
} # END .runAssessCA()