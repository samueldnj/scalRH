# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
#
# SCALfuns.R
#
# Functions for hierSCAL TMB model.
# 
# To Do:
#   6. Add in joint priors - be clever about it, we need to avoid
#       joint q priors on 
#
#
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

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
  cat( "\nMSG (saveFit) Created assessment fit folder ",folder,"in ./Outputs/fits/.\n" )
  
  # Save path so we can access it later
  reports$path <- path

  # Save reports object
  save(reports,file = file.path(path,paste(folder,".RData",sep="")))

  # Write out fitReport
  fitRepPath <- file.path(path,"fitReport.csv")
  write.csv( reports$phaseList$fitReport, file = fitRepPath )

  # Save plots?
  if(controlList$ctrl$plots)
    savePlots(  fitObj = reports,
                useRep = reports$plotRep,
                saveDir = path  )

  # Copy control file to sim folder for posterity
  file.copy(from=ctlFile,to=file.path(path,"fitCtlFile.txt"))
  file.copy(from="hierSCAL.cpp",to=file.path(path,"hierSCAL.cpp"))
  # Done
} # END fitHierSCAL()

# loadFit()
# Loads the nominated fit reports object into memory, 
# so that plot functions can be called
# inputs:   fit=ordinal indicator of sim in project folder
# ouputs:   NULL
# usage:    Prior to plotting simulation outputs
.loadFit <- function( fit = 1 )
{
  # List directories in project folder, remove "." from list
  dirList <- list.dirs (path="./Outputs/fits",full.names = FALSE,
                        recursive=FALSE)
  # Restrict to sim_ folders, pick the nominated simulation
  fitList <- dirList[grep(pattern="fit",x=dirList)]
  folder <- fitList[fit]

  # Load the nominated blob
  reportsFileName <- paste(folder,".RData",sep="")
  reportsPath <- file.path(getwd(),"Outputs/fits",folder,reportsFileName)
  load ( file = reportsPath )

  assign( "reports",reports,pos=1 )
  cat("MSG (loadFit) Reports in ", folder, " loaded from ./Outputs/fits/\n", sep="" )

}


# rerunPlots()
# Function to redo all plots from a given
# fit report object
rerunPlots <- function( fitID = 1, rep = "FE" )
{
  # Load the fit object
  .loadFit(fitID)

  if(!is.null(reports$repFE))
    reports$repFE <- renameReportArrays( reports$repFE, reports$data )

  if(!is.null(reports$repInit))
    reports$repInit <- renameReportArrays( reports$repInit, reports$data )

  if(!is.null(reports$repRE))
    reports$repRE <- renameReportArrays( reports$repRE, reports$data )

  if(!is.null(reports$repOpt))
    reports$repOpt <- renameReportArrays( reports$repOpt, reports$data )

  if(is.null(reports$refPoints))
    reports$repOpt <- calcRefPts(reports$repOpt)

  # rerun savePlots
  savePlots(  fitObj = reports,
              useRep = reports$plotRep,
              saveDir = reports$path )

  cat("MSG (rerunPlots) Plots redone in ", reports$path, "\n", sep = "")
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
  dataObj <- obj$data
  hypoObj <- obj$hypo
  ctrlObj <- obj$ctrl
  phases  <- obj$phases

  # Get model dimensions
  nA_s      <- dataObj$nA_s
  minA_s    <- dataObj$minA_s
  A1_s      <- dataObj$A1_s
  A2_s      <- dataObj$A2_s
  nL_s      <- dataObj$nL_s
  minL_s    <- dataObj$minL_s
  fYear     <- dataObj$fYearData
  lYear     <- dataObj$lYearData
  fYearIdx  <- dataObj$fYearIdx
  lYearIdx  <- dataObj$lYearIdx
  nT        <- lYear - fYear + 1
  years     <- fYear:lYear
  yrChar    <- as.character(years)

  # Now load the fleetIDs
  loadStockSpecNameLists()

  # Track used fleets
  allFleets     <- c( names(commFleetYrRange), names(surveyIDs) )
  useFleets     <- c( dataObj$commNames_g, dataObj$survNames_g)
  useFleetsIdx  <- which( allFleets %in% useFleets )
  nF            <- length(useFleets)
  totF          <- length(allFleets)

  # Track included species
  allSpecies    <- names(survSpecNames)
  useSpecies    <- dataObj$species
  useSpecIdx    <- which( allSpecies %in% useSpecies )
  nS            <- length(useSpecies)
  # Track included stocks
  allStocks    <- names(stocksSurvey)
  useStocks    <- dataObj$stocks
  useStockIdx  <- which( allStocks %in% useStocks )
  nP           <- length(useStocks)  

  # Track fleets for which we'll use the index data
  # (this way we can turn off commercial indices)
  idxFleets    <- dataObj$idxFleets

  # Load data
  cat("Loading data for assessment\n")

  # Load index data
  idxDataFiles <- dataObj$idxFiles
  idxDataFilePaths <- file.path("./Data",idxDataFiles )
  for( i in 1:length(idxDataFilePaths))
    load( idxDataFilePaths[i] )

  # Update makeIndexArray to use the commFleetYrRange vector
  I_spft <- makeIndexArray( relBio = relBioList_Survey,
                            commCPUE = commCPUEList,
                            nP = length(allStocks), years = fYear:lYear,
                            collapseComm = FALSE,
                            scaleComm = dataObj$commCPUEscalar )
  I_spft <- I_spft[useSpecies,useStocks,useFleets,, drop = FALSE]

  # Cut out indices from outside of idxFleets
  # and from outside fYearIdx:lYearIdx
  outYears <- years[! years %in% fYearIdx:lYearIdx ]
  I_spft[,,!useFleets %in% idxFleets,] <- -1
  I_spft[,,,as.character(outYears)] <- -1

  minTimeIdx_spf <- array(nT-1, dim = c(nS,nP,nF) )
  for( sIdx in 1:nS )
    for( pIdx in 1:nP )
      for( fIdx in 1:nF )
        if( any(I_spft[sIdx,pIdx,fIdx,]>0))
          minTimeIdx_spf[sIdx,pIdx,fIdx] <- min(which(I_spft[sIdx,pIdx,fIdx,]>0)) - 1

  # Load catch data
  commCatch <- read.csv(  file.path("./Data/",dataObj$catchData["comm"]), 
                          header = TRUE,
                          stringsAsFactors = FALSE )

  survCatchPath  <- file.path("./Data",dataObj$catchData["survey"])
  load(  file = survCatchPath )

  # Update makeCatchDiscArrays to use the commFleetYrRange vector
  catchDiscArrays <- makeCatchDiscArrays( commData = commCatch,
                                          survData = surveyCatch,
                                          stocks = stocksCommBio,
                                          speciesCodes = specCodes,
                                          years = fYear:lYear,
                                          modernYear = 1996,
                                          nF = nF, 
                                          collapseComm = FALSE,
                                          fleetIDs = useFleets )
  # Extract catch and discards
  C_spft <- catchDiscArrays$C_spft[useSpecies,useStocks,,, drop = FALSE]
  D_spft <- catchDiscArrays$D_spft[useSpecies,useStocks,,, drop = FALSE]

  if( !dataObj$modelDisc )
    C_spft <- C_spft + D_spft

  # Sum catch for initial B0 estimate
  sumCat_sp <- apply( X = C_spft, FUN = sum, MARGIN = c(1,2) )

  # Now use C_spft to weed out fleets that shouldn't be included
  # for now
  noFleetCatch <- c()
  for( fleetID in useFleets )
    if( all(C_spft[,,fleetID,] <= 0) )
    {
      noFleetCatch <- c(noFleetCatch,fleetID)
    }

  # Now weed
  useFleets     <- useFleets[!useFleets %in% noFleetCatch]
  useFleetsIdx  <- which(allFleets %in% useFleets)
  nF            <- length(useFleets)

  # Remove no catch fleets from C, D and I
  C_spft <- C_spft[,,useFleets,,drop = FALSE]
  D_spft <- D_spft[,,useFleets,,drop = FALSE]
  I_spft <- I_spft[,,useFleets,,drop = FALSE]
  

  # Load growth data - use this
  # to inform the vonB model in the
  # assessment
  growthFiles <- dataObj$growthFiles
  growthFilePaths <- file.path("./Data",growthFiles )
  for( i in 1:length(growthFiles))
    load(growthFilePaths[i])

  # Either use modeled uncertainty in growth or fixed CV
  if( hypoObj$sigmaL == "Model" )
  {
    sigmaLa_s <- vonBFits$sigA_s[useSpecIdx]
    sigmaLb_s <- vonBFits$sigB_s[useSpecIdx]
  } else {
    sigmaLa_s = rep(hypoObj$sigmaL,nS)
    sigmaLb_s = rep(0,nS)
  }

  # Make the age-length key - this might be useful,
  # or we can include it so we can integrate a growth
  # model
  # Aggregate ages into plus groups later.
  ALFreq_spalftx <- makeALFreq( ALFreqList = ALfreq,
                                years = years,
                                gears = useFleets,
                                maxA = max(nA_s), maxL = max(nL_s) )

  # Pare down to specific fleets for growth
  # CAAL data - Lee et al 2019+ (shared privately)
  ALfleetNames    <- dimnames(ALFreq_spalftx)[[5]]
  growthFleetIdx  <- which(ALfleetNames %in% dataObj$growthFleets)
  
  # Set all observations outside those fleets to 0
  ALFreq_spalftx[,,,,-growthFleetIdx,,] <- 0


  # Load age and length compositions
  load(file.path("./Data",dataObj$ageData))
  load(file.path("./Data",dataObj$lenData))

  # Create compositional arrays
  age_aspftx <- makeCompsArray(  compList = ageComps,
                                plusGroups = nA_s,
                                minX = minA_s,
                                collapseComm = FALSE,
                                fleetIDs = useFleets,
                                years = fYear:lYear,
                                xName = "ages",
                                minSampSize = dataObj$minAgeSampSize )
  age_aspftx <- age_aspftx[,useSpecies,useStocks,useFleets,,, drop = FALSE]


  len_lspftx <- makeCompsArray( compList = lenComps,
                                plusGroups = nL_s,
                                minX = minL_s,
                                collapseComm = FALSE,
                                fleetIDs = useFleets,
                                years = fYear:lYear,
                                xName = "length",
                                minSampSize = dataObj$minLenSampSize )
  len_lspftx <- len_lspftx[,useSpecies,useStocks,useFleets,,, drop = FALSE]

  nX <- 2

  # Combine sexes if not sex structured
  if(!dataObj$sexStructured)
  {
    ALFreq_spalftx <- ALFreq_spalftx[,,,,,,1,drop = FALSE] + ALFreq_spalftx[,,,,,,2,drop = FALSE]
    dimnames(ALFreq_spalftx)[[7]] <- "both"
    age_aspftx <- age_aspftx[,,,,,1,drop = FALSE] + age_aspftx[,,,,,2,drop = FALSE] 
    dimnames(age_aspftx)[[6]] <- "both"
    len_lspftx <- len_lspftx[,,,,,1,drop = FALSE] + len_lspftx[,,,,,2,drop = FALSE] 
    dimnames(len_lspftx)[[6]] <- "both"

    nX <- 1
  }



  # Calculate the number of selectivity deviations
  # from length and age distributions
  nSelDevs <- 0
  for( fIdx in which( useFleets %in% hypoObj$tvSelFleets ) )
    for( sIdx in 1:nS )
      for( pIdx in 1:nP )
      {
        # Start with empty vectors of year indices with
        # observations
        lenYrs <- c()
        ageYrs <- c()
        # Cut down to those with age observations
        if( any(age_aspftx[1,sIdx,pIdx,fIdx,yrChar,] >= 0) &  dataObj$ageLikeWt > 0 )
          ageYrs <- which(age_aspftx[1,sIdx,pIdx,fIdx,yrChar,,drop = FALSE] >= 0,arr.ind = TRUE)[,5]
        # And if length observations are included, add them
        if( any(len_lspftx[1,sIdx,pIdx,fIdx,yrChar,] >= 0) &  dataObj$lenLikeWt > 0 )
          lenYrs <- which(len_lspftx[1,sIdx,pIdx,fIdx,yrChar,,drop = FALSE] >= 0, arr.ind = TRUE)[,5]

        # Take union of age/len years, then add number of years to 
        # estimate total number of selectivity deviations
        nSelDevs <- nSelDevs + length( union( lenYrs, ageYrs ) )
      } 

  # Calculate number of q deviations
  nqDevs <- 0
  tvqFleetIdx <- which(idxFleets %in% hypoObj$tvqFleets )
  for( fIdx in tvqFleetIdx )
    for( sIdx in 1:nS )
      for( pIdx in 1:nP )
      {
        nObs <- (length( which( I_spft[sIdx,pIdx,idxFleets[fIdx],] > 0) ) - 1)
        if(nObs > 0)
          nqDevs <- nqDevs + nObs
      }




  # Now make a vector to switch time varying selectivity
  # on and off
  tvSelFleets <- rep(0,nF)
  names(tvSelFleets) <- useFleets
  tvSelFleets[ useFleets %in% hypoObj$tvSelFleets ] <- 1 

  # And the same for time-varying catchability
  tvqFleets <- rep(0,nF)
  names(tvqFleets) <- useFleets
  tvqFleets[ useFleets %in% hypoObj$tvqFleets ] <- 1 

  # And the same for F regularisation
  regFfleets <- rep(0,nF)
  names(regFfleets) <- useFleets
  regFfleets[ useFleets %in% hypoObj$regFfleets ] <- 1 

  # Load maturity ogives
  load(file.path("./Data",dataObj$matFiles))

  # Get growth and maturity info
  # Vectors to hold
  # Mat
  xMat50 <- numeric(length = nS)
  xMat95 <- numeric(length = nS)
  # Wt
  LWa_s <- numeric(length = nS)
  LWb_s <- numeric(length = nS)

  matOgives <- matOgives[useSpecies]

  for( sIdx in 1:nS )
  {
    specName <- useSpecies[sIdx]
    if( hypoObj$matX == "length")
    {
      xMat50[sIdx] <- matOgives[[specName]]$coastWide$all$length$x50
      xMat95[sIdx] <- matOgives[[specName]]$coastWide$all$length$x95
    }
    if(hypoObj$matX == "age")
    {
      xMat50[sIdx] <- matOgives[[specName]]$coastWide$all$age$x50
      xMat95[sIdx] <- matOgives[[specName]]$coastWide$all$age$x95
    }

    LWa_s[sIdx] <- exp(coef(wtLen[[specName]]$wtLen$coastWide$wtLenAll)[1])
    LWb_s[sIdx] <- coef(wtLen[[specName]]$wtLen$coastWide$wtLenAll)[2]
  }

  matOgives <- NULL
  wtLen <- NULL
  gc()
  
  # Count number of free F parameters
  nPosCatch <- length(which(C_spft > 0) )

  # Initialised in a fished state (0 == unfished, 1 == fished)
  initFished_s     <- hypoObj$initFished[useSpecies]
  yFirstRecDev_s   <- hypoObj$yFirstRecDev[useSpecies]
  yLastRecDev_s    <- hypoObj$yLastRecDev[useSpecies]

  if( all(initFished_s == 0) )
    phases$omegaRinit_vec <- -1

  # Non-eq initialisation deviations
  nInitDevs <- sum( nP * (nA_s[useSpecIdx] * initFished_s ) )

  # Generate tFirstRecDev from yFirstRecDev and fYear - these are 
  # created adjusting for TMB zero initialisation
  tFirstRecDev_s <- yFirstRecDev_s - fYear 
  tFirstRecDev_s[ tFirstRecDev_s < 1 ] <- 1
  # same for last rec deviation
  tLastRecDev_s <- yLastRecDev_s - fYear
  tLastRecDev_s[ tLastRecDev_s > nT ] <- nT

  nRecDevs <- nP * sum( (tLastRecDev_s - tFirstRecDev_s + 1))


  calcIndex_spf <- array(0, dim = c(nS,nP,nF) )
  for( s in 1:nS )
    for( p in 1:nP )
      for( f in 1:nF )
        if( any(I_spft[s,p,f,] > 0) )
          calcIndex_spf[s,p,f] <- 1

  if( hypoObj$selX == "length")
  {
    # length at Sel50_sf initial value
    xSel50_sf <- matrix(  c(    36, 36, 29, 31, 33, 33, 33,
                                35, 35, 23, 25, 25, 25, 25,
                                33, 30, 18, 18, 20, 20, 20,
                                38, 38, 30, 30, 30, 30, 30,
                                35, 35, 35, 35, 35, 35, 35 ), 
                                nrow = length(allSpecies),
                                ncol = length(allFleets), 
                                byrow = TRUE )

    xSel50_sf <- xSel50_sf[useSpecIdx,useFleetsIdx]

    # length at SelStep_sf initial value
    xSelStep_sf <- matrix(  c(    2, 2, 2, 2, 2, 2, 2,
                                  2, 2, 2, 2, 2, 2, 2,
                                  2, 4, 5, 5, 5, 5, 5,
                                  2, 2, 6, 6, 6, 6, 6,
                                  2, 2, 2, 2, 2, 2, 2 ), 
                                  nrow = length(allSpecies),
                                  ncol = length(allFleets), 
                                  byrow = TRUE )

    xSelStep_sf <- xSelStep_sf[useSpecIdx,useFleetsIdx,drop = FALSE]
  }

  if( hypoObj$selX == "age")
  {
    # xSel50_sf initial value
    xSel50_sf <- matrix(  c(    8, 7, 5, 6, 5, 5, 5,
                                6, 5, 5, 5, 5, 5, 5,
                                5, 5, 5, 5, 5, 5, 5,
                                5, 5, 5, 5, 5, 5, 5,
                                5, 5, 5, 5, 5, 5, 5 ), 
                                nrow = length(allSpecies),
                                ncol = length(allFleets), 
                                byrow = TRUE )

    xSel50_sf <- xSel50_sf[useSpecIdx,useFleetsIdx,drop = FALSE]

    # xSelStep_sf initial value
    xSelStep_sf <- matrix(  c(    3, 3, 2, 3, 2, 2, 2,
                                  2, 2, 2, 2, 2, 2, 2,
                                  2, 2, 2, 2, 2, 2, 2,
                                  2, 2, 2, 2, 2, 2, 2,
                                  2, 2, 2, 2, 2, 2, 2 ), 
                                  nrow = length(allSpecies),
                                  ncol = length(allFleets), 
                                  byrow = TRUE )

    xSelStep_sf <- xSelStep_sf[useSpecIdx,useFleetsIdx,drop = FALSE]
  }

  lnxSel50_sf <- array(log(xSel50_sf),dim = c(nS,nF))
  lnxSelStep_sf <- array(log(xSelStep_sf),dim = c(nS,nF))

  # Use useSpec, useStocks, useFleets and yrChar to
  # reduce the data arrays so that the fits are dynamic

  # Pull prior mean values from hypoObj to keep pars list tidier
  mh    <- hypoObj$hPriorMean
  sdh   <- hypoObj$hPriorSD
  mq    <- hypoObj$mq
  sdq   <- hypoObj$sdq
  mM    <- hypoObj$mM
  sdM   <- hypoObj$sdM

  # Observation errors
  tau2ObsIGa  <- hypoObj$tau2ObsIGa[useFleetsIdx]
  tau2ObsMode <- hypoObj$tau2ObsPriorMode[useFleetsIdx]
  tau2ObsIGb  <- (tau2ObsIGa + 1) * tau2ObsMode

  nComm       <- length(dataObj$commNames_g)
  nSurv       <- length(dataObj$survNames_g)
  group_f     <- c(rep(0,nComm),1,rep(2,nSurv-1))

  # Generate initial L1 and L2 values
  initL1_s    <- numeric( length = nS )
  initL2_s    <- numeric( length = nS )

  survNames <- dataObj$survNames_g
  survNames <- survNames[survNames %in% useFleets]

  calcMeanLenAge <- function( sIdx = 1, ALK = ALFreq_spalftx,
                              age = A1_s, fleets = survNames )
  {
    # Pull length at the given age
    lenAtAge <- ALK[sIdx,,age[sIdx],,fleets,,,drop = FALSE]
    lenAtAge <- apply(X = lenAtAge, FUN = sum, MARGIN = c(4))
    # Calculate number of observations
    nLenAtAge <- sum(lenAtAge,na.rm = T)
    # Multiply length bin by freq
    totLenAtAge <- sum((1:length(lenAtAge)) * lenAtAge )
    # Calc mean
    meanLenAge <- sum(totLenAtAge)/nLenAtAge

    if(!is.finite(meanLenAge))
      browser(cat("Non-finite mean length-at-age."))

    meanLenAge
  }

  # Pull the fleets used for growth models
  growthFleets <- dataObj$growthFleets

  initL1_s <- sapply( X = 1:nS, 
                      FUN = calcMeanLenAge, 
                      ALK = ALFreq_spalftx, 
                      age = A1_s[useSpecies],
                      fleets = growthFleets[growthFleets %in% useFleets] )
  initL2_s <- sapply( X = 1:nS, 
                      FUN = calcMeanLenAge, 
                      ALK = ALFreq_spalftx, 
                      age = A2_s[useSpecies],
                      fleets = growthFleets[growthFleets %in% useFleets]  )


  # Generate the data list
  data <- list( I_spft            = I_spft,
                C_spft            = C_spft,
                D_spft            = D_spft,
                ALK_spalftx       = ALFreq_spalftx,
                age_aspftx        = age_aspftx,
                len_lspftx        = len_lspftx,
                group_f           = as.integer(group_f),
                A_s               = as.integer(nA_s[useSpecies]),
                minA_s            = as.integer(minA_s[useSpecies]),
                L_s               = as.integer(nL_s[useSpecies]),
                minL_s            = as.integer(minL_s[useSpecies]),
                lenD_s            = rep(0,nS),
                swRinit_s         = as.integer(initFished_s[useSpecies]),
                parSwitch         = 0,
                calcIndex_spf     = calcIndex_spf,
                tvSelFleets       = tvSelFleets,
                tvqFleets         = tvqFleets,
                regFfleets        = regFfleets,
                idxLikeWt         = dataObj$idxLikeWt,
                ageLikeWt         = dataObj$ageLikeWt,
                lenLikeWt         = dataObj$lenLikeWt,
                growthLikeWt      = dataObj$growthLikeWt,
                tFirstRecDev_s    = as.integer(tFirstRecDev_s),
                tLastRecDev_s     = as.integer(tLastRecDev_s),
                minAgeProp        = dataObj$minAgeProp,
                minLenProp        = dataObj$minLenProp,
                minPAAL           = dataObj$minPAAL,
                matX              = hypoObj$matX,
                selX              = hypoObj$selX,
                lenComps          = hypoObj$lenCompMethod,
                lambdaB0          = dataObj$lambdaB0,
                minTimeIdx_spf    = minTimeIdx_spf,
                nBaranovIter      = ctrlObj$nBaranovIter,
                lambdaBaranovStep = ctrlObj$lambdaBaranovStep,
                A1_s              = A1_s[useSpecies],
                A2_s              = A2_s[useSpecies] )


  # Generate parameter list
  pars <- list( ## Leading biological pars ##
                lnB0_sp           = log(sumCat_sp),
                logitSteep        = log((mh - .2)/(1 - mh)),
                lnM               = log(mM),
                lnL2step_s        = log(initL2_s - initL1_s),
                lnvonK_s          = log(hypoObj$initVonK[useSpecies]),
                lnL1_s            = log(initL1_s),
                # Stock specific growth pars
                deltaL2_sp        = array(0,dim = c(nS,nP)),
                lnsigmaL2_s       = rep(log(hypoObj$sigmaL2),nS),
                deltaVonK_sp      = array(0,dim = c(nS,nP)),
                lnsigmavonK_s     = rep(log(hypoObj$sigmavonK),nS),
                deltaL2_sx        = array(0,dim = c(nS,nX)),
                deltaVonK_sx      = array(0,dim = c(nS,nX)),
                lnsigmaL2         = log(hypoObj$sigmaL2),
                lnsigmavonK       = log(hypoObj$sigmavonK),
                # process error in growth model
                lnsigmaLa_s       = log(sigmaLa_s),
                sigmaLb_s         = sigmaLb_s,
                # L-W conversion
                LWa_s             = LWa_s,
                LWb_s             = LWb_s,
                # Maturity
                xMat50_s          = xMat50,
                xMat95_s          = xMat95,
                ## Observation models ##
                # fleet catchability and obs idx SD
                lnq_spf           = array(0,dim =c(nS,nP,nF)),
                lntauObs_spf      = array(-1,dim = c(nS,nP,nF)),
                # Selectivity
                lnxSel50_sf       = lnxSel50_sf,
                lnxSelStep_sf     = lnxSelStep_sf,
                epsxSel50_spf     = array(0,dim = c(nS,nP,nF)),
                epsxSelStep_spf   = array(0,dim = c(nS,nP,nF)),
                # Fishing mortality
                # lnF_spft          = rep(-1,length = nPosCatch),
                # Catch and discards obs SD
                lntauC_f          = rep(log(0.01),nF),
                lntauD_f          = rep(log(0.01),nF),
                ## Multilevel priors ##
                # Selectivity
                muxSel50_sg       = array(3.4, dim = c(nS,3)),
                muxSel95_sg       = array(4.2, dim = c(nS,3)),
                sigmaxSel50_sg    = array(.2, dim = c(nS,3)),
                sigmaxSel95_sg    = array(.2, dim = c(nS,3)),
                # Catchability
                lnqbarSyn_s       = rep(log(hypoObj$mq),nS),
                lntauqSyn_s       = rep(log(hypoObj$tauqSyn),nS),
                lnqbarSyn         = log(hypoObj$mq),
                lntauqSyn         = log(hypoObj$tauqSyn),
                mqSurveys         = hypoObj$mq,
                sdqSurveys        = hypoObj$sdq,
                # Steepness
                lnsigmah_s        = rep(log(sdh),nS),
                logit_muSteep     = log((mh - .2)/(1 - mh)),
                lnsigmah          = log(sdh),
                # Mortality
                lnsigmaM_s        = rep( log(sdM), nS ),
                ln_muM            = log(mM),
                lnsigmaM          = log(sdM),
                # IG Prior on obs error SD
                IGatau_f          = tau2ObsIGa,
                IGbtau_f          = tau2ObsIGb,
                # Species effect on steepness
                epsSteep_s        = rep(0,nS),
                # Species effect on M
                epsM_s            = rep(0,nS),
                # Species/stock effect on steepness
                epsSteep_sp       = array(0, dim = c(nS,nP)),
                # Species/stock effect on M
                epsM_sp           = array(0, dim = c(nS,nP)),
                epsM_sx          = array(0, dim = c(nS,nX)),
                # Recruitment resids
                omegaR_vec        = rep( 0, nRecDevs),
                omegaRinit_vec    = rep( 0, nInitDevs ),
                lnsigmaR_sp       = array( log(hypoObj$sigmaR), dim = c(nS,nP)),
                # Correlation in recruitment resids
                logitRCorr_chol   = rep(0, nS * nP),
                logitRgamma_sp    = array( 0, dim = c(nS,nP)),
                # Time-varying selectivity
                epsxSel50_vec     = rep( 0, nSelDevs ),
                epsxSelStep_vec   = rep( 0, nSelDevs ),
                lnsigmaSel        = log( hypoObj$sigmaSelDevs ),
                # Time-varying catchability
                epslnq_vec        = rep( 0, nqDevs ),
                lnsigmaepslnq     = log( hypoObj$sigmaqdevs ),
                ## Single-level priors ##
                # Priors on selectivity
                pmlnxSel50_sf     = array(log(xSel50_sf),dim = c(nS,nF) ),
                pmlnxSelStep_sf   = array(log(xSelStep_sf),dim = c(nS,nF) ),
                cvxSel            = hypoObj$cvxSel,
                pmlnL2_s          = log(initL2_s),
                pmlnL1_s          = log(initL1_s),
                cvL2              = .1,
                cvL1              = .1,
                pmlnVonK          = log(.3),
                cvVonK            = .1 )



  # Generate special entries for the 
  # base map list that are sensitive to useFleetsIdx
  qmap_spf <- array(  1 + 1:(nS*nP*nF),
                      dim = c(nS,nP,nF) )
  qmap_spf[calcIndex_spf == 0] <- NA

  taumap_spf <- qmap_spf + 1e3

  # Map selectivity at length
  selMap_spf <- array(NA, dim = c(nS,nP,nF) )
  # Make unique initially
  for(s in 1:nS)
    for( f in 1:nF)
    {
      for(p in 1:nP)  
        if( any(age_aspftx[,s,p,f,,] > 0) | any(len_lspftx[,s,p,f,,] > 0) )
          selMap_spf[s,p,f] <- 2e3 + s * (nF - 1) * (nP - 1) + f * (nP - 1) + p
    }

  # Collapse to species/fleet
  selMap_sf <- apply( X = selMap_spf, 
                      FUN = sum, MARGIN =c(1,3),
                      na.rm = T )
  selMap_sf[ selMap_sf == 0 ] <- NA

  if( hypoObj$identSel )
  {
    fleetGps <- hypoObj$fleetGroups
    nGps     <- length(fleetGps)
    for( gIdx in 1:length(fleetGps) )
    {
      gpFleets <- fleetGps[[gIdx]]
      gpFleets <- gpFleets[gpFleets %in% useFleets]
      estIdx <- which(!is.na(selMap_spf[,gpFleets,drop = FALSE]),arr.ind =T)
      selMap_spf[estIdx[,1],estIdx[,2],gpFleets[estIdx[,3]]] <- gIdx + 3e3 + estIdx[,1] * nGps
    }
  }

  # generate base map for TMBphase()
  map <- list(  lnq_spf           = factor(qmap_spf),
                lntauObs_spf      = factor(taumap_spf),
                epsxSel50_spf     = factor(selMap_spf),
                epsxSelStep_spf   = factor(selMap_spf + 105),
                lnxSel50_sf       = factor(selMap_sf + 50),
                lnxSelStep_sf     = factor(selMap_sf + 100) )

  # Turn off tv sel deviations if not being used
  if( !hypoObj$tvSel | nSelDevs == 0 )
  {
    phases$epsxSel50_vec    <- -1
    phases$epsxSelStep_vec  <- -1
  }

  # Turn off tvq deviations if not used
  if( !hypoObj$tvq | nqDevs == 0 )
    phases$epslnq_vec       <- -1

  if( nP == 1)
  {
    phases$epsxSel50_spf    <- -1
    phases$epsxSelStep_spf  <- -1
    phases$epsM_sp          <- -1
    phases$epsSteep_sp      <- -1
    phases$deltaL2_sp       <- -1
    phases$deltaVonK_sp     <- -1
  }

  if( nS == 1 )
  {
    phases$epsSteep_s       <- -1
    phases$epsM_s           <- -1
  }

  checkDat <- lapply( X = data, FUN = .checkNA )
  checkPar <- lapply( X = pars, FUN = .checkNA )

  if( any(checkDat) | any( checkPar ) )
  {
    browser(cat("Data or Pars have NAs\n") )
  }

  cat("\nFitting hierSCAL for ", paste(useSpecies,collapse = ","), "\n", sep = "")
  gc()

  phaseList <- TMBphase(  data = data, 
                          parameters = pars, 
                          random = hypoObj$RE, 
                          phases = phases, 
                          base_map = map,
                          maxPhase = ctrlObj$maxPhase,
                          model_name = "hierSCAL",
                          optimizer = "nlminb",
                          silent = ctrlObj$quiet,
                          maxEval = ctrlObj$maxFunEval,
                          maxIter = ctrlObj$maxIterations,
                          calcSD = ctrlObj$calcSD,
                          parBen = ctrlObj$parBen,
                          intMethod = ctrlObj$intMethod,
                          mcChainLength = ctrlObj$mcChainLength,
                          mcChains = ctrlObj$mcChains,
                          savePhases = ctrlObj$savePhases ) 

  repOpt <- phaseList$repOpt


  maxSuccPhz <- phaseList$maxPhaseComplete
  if( maxSuccPhz > 0)
    plotRep <- "opt"

  # Calculate refPts for repOpt
  if( ctrlObj$calcRefPts )
  {
    repOpt <- calcRefPts( phaseList$repOpt )
  }
  

  # Update names on report objects
  outList <- list(  repOpt = renameReportArrays(repOpt,data),
                    sdrepOpt = phaseList$sdrep,
                    sdOptOutput = phaseList$optOutput,
                    fYear = fYear, 
                    lYear = lYear,
                    gearLabs = useFleets,
                    species = useSpecies,
                    stocks = useStocks,
                    map = phaseList$map,
                    data = data,
                    pars = pars,
                    phaseList = phaseList )



  return(outList)
} # END .runHierSCAL()

# Custom TMBphase() function for running hierSCAL in phases. 
# Modified from the version in Kasper Kristensen's 
# TMB_contrib_R github repository 
# https://github.com/kaskr/TMB_contrib_R
# Author:Gavin Fay email: gfay42@gmail.com
# 
# Main modification adds a base map list for array based parameters that
# have (a) identified parameters or (b) only some elements
# activated due to missing data (e.g. selectivity pars for
# fleets in specific areas). Doing it this way reduces number
# of loops in the model, speeding up fitting time.
TMBphase <- function( data, 
                      parameters, 
                      random = NULL, 
                      phases, 
                      base_map = list(),
                      maxPhase = NULL,
                      model_name = "hierSCAL",
                      optimizer = "nlminb",
                      silent = FALSE,
                      calcSD = FALSE,
                      maxEval = 1000,
                      maxIter = 1000,
                      regFMaxPhase = 3,
                      parBen = FALSE,
                      intMethod = "RE",
                      mcChainLength = 100,
                      mcChains = 1,
                      savePhases = TRUE ) 
{
  # function to fill list component with a factor
  # of NAs
  fill_vals <- function(x,vals){ factor( rep( vals, length(x) ) ) }


  # compile the model
  DLL_use <- model_name  
  
  #loop over phases
  if(!is.null(maxPhase))
    maxPhase <- min( maxPhase, max(unlist(phases) ) )
  else maxPhase <- max(unlist(phases))


  # Make a data.frame that will hold the phase info
  fitReport <- matrix(NA, nrow = maxPhase + 1, ncol = 7 )
  colnames(fitReport) <- c("phase","objFun","maxGrad","nPar","convCode","convMsg", "time")
  fitReport <- as.data.frame(fitReport)

  fitReport$phase <- c(1:maxPhase,"RE")

  phaseReports <- vector(mode = "list", length = maxPhase)

  # generate a list of outputs to return
  # to runHierSCAL, initialise 
  # a success flag at TRUE
  outList <- list( success = TRUE )

  for( phase_cur in 1:maxPhase ) 
  {
    gc()
    # Start timing
    tBegin <- proc.time()[3]

    # work out the map for this phase
    # if the phase for a parameter is greater than the current phase 
    # or a negative value, then map will contain a factor filled with NAs
    map_use <- base_map
    j <- length(map_use)
    for( i in 1:length(parameters) ) 
    {
      parName <- names(parameters)[i]

      if( parName %in% names(phases) )
      {
        if( (phases[[parName]] > phase_cur) | phases[[parName]] < 0 ) 
        { 
          # Check if parName is included in the base_map
          if(parName %in% names(map_use))
            map_use[[parName]] <- fill_vals(parameters[[i]],NA)
          else
          {
            j <- j + 1
            map_use[[j]] <- fill_vals(parameters[[i]],NA)
            names(map_use)[j] <- parName
          }

        }
      } else {
        j <- j + 1
        map_use[[j]] <- fill_vals(parameters[[i]],NA)
        names(map_use)[j] <- parName
      }

    }
    
    if( phase_cur > regFMaxPhase )
      data$regFfleets <- rep(0,length(data$regFfleets))
  
    #remove the random effects if they are not estimated
    random_use <- random[ !random %in% names(map_use) ]
  
    # initialize the parameters at values in previous phase
    params_use <- parameters
    if( phase_cur > 1 ) 
      params_use <- obj$env$parList( opt$par )

    # Fit the model
    obj <- TMB::MakeADFun(  data = data,
                            parameters = params_use,
                            random = NULL,
                            DLL= DLL_use,
                            map= map_use,
                            silent = silent )

    # Run benchmark for parallel accumulation
    if(parBen &  phase_cur == 1)
      outList$phase1Benchmark <- benchmark(obj, cores = 1:(detectCores()-1))

    # Create a control list for the assessment model
    tmbCtrl <- list(  eval.max = maxEval, 
                      iter.max = maxIter  )

    if( phase_cur == 1 )
    {
      repInit <- obj$report()

      checkInitNaN    <- lapply( X = repInit, FUN = .checkNaN )
      checkInitFinite <- lapply( X = repInit, FUN = .checkFinite )
      if(any(unlist(checkInitNaN)) | any(checkInitFinite))
        browser(beep(expr=cat("NaN or Inf items in repInit\n")))
    }

    cat("\nStarting optimisation for phase ", phase_cur, "\n\n")
    # Try the optimisation
    opt <- try( nlminb (  start     = obj$par,
                          objective = obj$fn,
                          gradient  = obj$gr,
                          control   = tmbCtrl ) )

    # break if there is an issue
    if( class(opt) == "try-error" )
    {

      cat("\nOptimisation halted due to error\n")

      outList$success                   <- FALSE
      outList$maxPhaseComplete          <- phase_cur - 1

      if( savePhases )
      {
        phaseReports[[phase_cur]]$opt     <- opt
        phaseReports[[phase_cur]]$success <- FALSE
      }

      break
    }
    # Save max phase complete
    outList$maxPhaseComplete          <- phase_cur

    # Save reports and optimisation
    # output
    if(savePhases)
    {
      phaseReports[[phase_cur]]$report  <- obj$report()
      phaseReports[[phase_cur]]$opt     <- opt
      phaseReports[[phase_cur]]$success <- TRUE
      phaseReports[[phase_cur]]$map     <- map_use
      phaseReports[[phase_cur]]$hess    <- obj$he()
    }

    # Update fitReport
    if(class(opt) != "try-error")
    {
      fitReport[phase_cur,]$objFun      <- obj$fn()
      fitReport[phase_cur,]$maxGrad     <- max(abs(obj$gr()))
      fitReport[phase_cur,]$nPar        <- length(opt$par)
      fitReport[phase_cur,]$convCode    <- opt$convergence
      fitReport[phase_cur,]$convMsg     <- opt$message
    }
    fitReport[phase_cur,]$time           <- (proc.time()[3] - tBegin)/60



    cat(  "\nPhase ", phase_cur, " completed with code ",
          opt$convergence, " and following message:\n", sep = "" )
    cat("\n", opt$message, "\n\n", sep = "" )
    
    # close phase loop
  }

  # Fit the model
  if( intMethod == "RE" & !is.null(random) &  class(opt) != "try-error" )
  { 
    randEffList <- list()

    tBegin      <- proc.time()[3]
    params_use  <- obj$env$parList( opt$par )

    obj <- TMB::MakeADFun(  data = data,
                            parameters = params_use,
                            random = c(random),
                            DLL= DLL_use,
                            map= map_use,
                            silent = silent )  

    tol10 <- 0.01
    
    TMB::newtonOption(obj, tol10 = tol10)
    
    if( parBen )
      randEffList$randEffBenchmark <- benchmark( obj, cores = 1:detectCores() )

    # Try the optimisation
    opt <- try( nlminb (  start     = obj$par,
                          objective = obj$fn,
                          gradient  = obj$gr,
                          control   = tmbCtrl ) )

    randEffList$spHess  <- obj$env$spHess(random = TRUE)
    randEffList$opt     <- opt

    # Update fitReports
    if(class(opt) != "try-error")
    {
      fitReport[maxPhase + 1,]$objFun      <- obj$fn()
      fitReport[maxPhase + 1,]$maxGrad     <- max(abs(obj$gr()))
      fitReport[maxPhase + 1,]$nPar        <- length(opt$par)
      fitReport[maxPhase + 1,]$convCode    <- opt$convergence
      fitReport[maxPhase + 1,]$convMsg     <- opt$message
    }

    fitReport[maxPhase + 1,]$time          <- (proc.time()[3] - tBegin)/60

    outList$randEffList <- randEffList
  }

  if( intMethod == "MCMC" & class(opt) != "try-error" )
  {
    tBegin      <- proc.time()[3]
    params_use  <- obj$env$parList( opt$par )

    obj <- TMB::MakeADFun(  data = data,
                            parameters = params_use,
                            random = NULL,
                            DLL= DLL_use,
                            map= map_use,
                            silent = silent )  
    mcmc <- tmbstan(  obj, 
                      init = "last.par.best", 
                      iter = mcChainLength,
                      chains = mcChains )

    outList$mcmc <- mcmc

  }
  
  if(outList$success & calcSD )
    outList$sdrep <- TMB::sdreport(obj)

  if( savePhases )
    outList$phaseReports      <- phaseReports

  outList$repOpt            <- obj$report()
  outList$objfun            <- obj$fn()
  outList$optOutput         <- opt
  outList$map               <- map_use
  outList$maxGrad           <- max(obj$gr())
  outList$fitReport         <- fitReport
  outList$totTime           <- sum(fitReport$time)
  
  return( outList )  

} # END TMBphase()

# renameReportArrays()
# Updates the dimension names of the arrays in the 
# report lists, as a way of making plot code more efficient later.
renameReportArrays <- function( repObj = repInit, datObj = data )
{
  # Just go down the list, but first do the objects with the same names

  repNames <- names(repObj)
  datNames <- names(datObj)

  bothNames <- repNames[ repNames %in% datNames ]

  for( itemName in bothNames )
  {
    dimnames(repObj[[itemName]]) <- dimnames(datObj[[itemName]])
  }

  # Recover names
  specNames   <- dimnames(datObj$I_spft)[[1]]
  stockNames  <- dimnames(datObj$I_spft)[[2]]
  gearNames   <- dimnames(datObj$I_spft)[[3]]
  yearNames   <- dimnames(datObj$I_spft)[[4]]
  lenNames    <- dimnames(datObj$len_lspftx)[[1]]
  ageNames    <- dimnames(datObj$age_aspftx)[[1]]
  sexNames    <- dimnames(datObj$age_aspftx)[[6]]


  # Ok, that's the data taken care of. There are still all the
  # new arrays that we created
  # Predicted data
  dimnames(repObj$I_spft_hat) <- dimnames(datObj$I_spft)
  dimnames(repObj$aDist_aspftx_hat) <- dimnames(datObj$age_aspftx)
  dimnames(repObj$lDist_lspftx_hat) <- dimnames(datObj$len_lspft)
  # State arrays
  dimnames(repObj$B_asptx) <- dimnames(datObj$age_aspftx)[c(1:3,5)]
  dimnames(repObj$N_asptx) <- dimnames(datObj$age_aspftx)[c(1:3,5)]
  dimnames(repObj$B_spt) <- dimnames(datObj$age_aspftx)[c(2:3,5)]
  dimnames(repObj$R_spt) <- dimnames(datObj$age_aspftx)[c(2:3,5)]
  dimnames(repObj$SB_spt) <- dimnames(datObj$age_aspftx)[c(2:3,5)]
  dimnames(repObj$Bv_spft) <- dimnames(datObj$age_aspftx)[c(2:5)]
  dimnames(repObj$predC_spft) <- dimnames(datObj$age_aspftx)[c(2:5)]
  dimnames(repObj$predCw_spft) <- dimnames(datObj$age_aspftx)[c(2:5)]
  dimnames(repObj$C_aspftx) <- dimnames(datObj$age_aspftx)[c(1:5)]
  dimnames(repObj$Cw_aspftx) <- dimnames(datObj$age_aspftx)[c(1:5)]
  dimnames(repObj$F_aspftx) <- dimnames(datObj$age_aspftx)[c(1:5)]
  dimnames(repObj$F_spft) <- dimnames(datObj$age_aspftx)[c(2:5)]
  dimnames(repObj$Z_aspxt) <- dimnames(datObj$age_aspftx)[c(1:3,6,5)]
  # Biological parameters
  dimnames(repObj$R0_sp) <- dimnames(datObj$age_aspftx)[c(2:3)]
  dimnames(repObj$B0_sp) <- dimnames(datObj$age_aspftx)[c(2:3)]
  dimnames(repObj$h_sp) <- dimnames(datObj$age_aspftx)[c(2:3)]
  dimnames(repObj$M_sp) <- dimnames(datObj$age_aspftx)[c(2:3)]
  dimnames(repObj$phi_sp) <- dimnames(datObj$age_aspftx)[c(2:3)]
  dimnames(repObj$reca_sp) <- dimnames(datObj$age_aspftx)[c(2:3)]
  dimnames(repObj$recb_sp) <- dimnames(datObj$age_aspftx)[c(2:3)]

  # Observation models
  dimnames(repObj$q_spf)        <- dimnames(datObj$age_aspftx)[c(2:4)]  
  dimnames(repObj$tau2Obs_spf)  <- dimnames(datObj$age_aspftx)[c(2:4)]  
  dimnames(repObj$tauObs_spf)   <- dimnames(datObj$age_aspftx)[c(2:4)]  
  dimnames(repObj$sel_lfspt)    <- list(  len = lenNames, 
                                          fleet = gearNames,
                                          species = specNames,
                                          stock = stockNames,
                                          year = yearNames )
  dimnames(repObj$sel_afsptx)     <- list( age = ageNames, 
                                          fleet = gearNames,
                                          species = specNames,
                                          stock = stockNames,
                                          year = yearNames,
                                          sex = sexNames )


  dimnames(repObj$ageRes_aspftx)   <- list(  age = ageNames, 
                                            species = specNames,
                                            stock = stockNames,
                                            fleet = gearNames,
                                            year = yearNames,
                                            sex = sexNames )  
  dimnames(repObj$tau2Age_spf)   <- list( species = specNames,
                                          stock = stockNames,
                                          fleet = gearNames )  
  dimnames(repObj$lenRes_lspftx)   <- list( length = lenNames, 
                                            species = specNames,
                                            stock = stockNames,
                                            fleet = gearNames,
                                            year = yearNames,
                                            sex = sexNames ) 
  dimnames(repObj$tau2Len_spf)   <- list( species = specNames,
                                          stock = stockNames,
                                          fleet = gearNames )  

  # Growth model quants
  dimnames(repObj$probLenAge_laspx) <- list( len = lenNames, 
                                              age = ageNames,
                                              species = specNames,
                                              stock = stockNames,
                                              sex = sexNames )

  dimnames(repObj$lenAge_aspx) <- list(     age = ageNames,
                                            species = specNames,
                                            stock = stockNames,
                                            sex = sexNames  )

  dimnames(repObj$probAgeLen_alspftx) <- list(  age = ageNames,
                                                len = lenNames,
                                                species = specNames,
                                                stock = stockNames,
                                                fleet = gearNames,
                                                year = yearNames,
                                                sex = sexNames  ) 

  dimnames(repObj$ageAtLenResids_alspftx) <- list(  age = ageNames,
                                                    len = lenNames,
                                                    species = specNames,
                                                    stock = stockNames,
                                                    fleet = gearNames,
                                                    year = yearNames,
                                                    sex = sexNames  )
  # Growth model parameters
  # vectors
  names(repObj$A1_s)      <- specNames
  names(repObj$A2_s)      <- specNames
  names(repObj$L1_s)      <- specNames
  names(repObj$L2_s)      <- specNames
  names(repObj$sigmaLa_s) <- specNames
  names(repObj$sigmaLb_s) <- specNames
  # arrays
  dimnames(repObj$L1_sp)    <- list(  species = specNames,
                                      stock = stockNames )
  dimnames(repObj$L2_sp)    <- list(  species = specNames,
                                      stock = stockNames )
  dimnames(repObj$vonK_sp)  <- list(  species = specNames,
                                      stock = stockNames )

  dimnames(repObj$L1_spx)    <- list( species = specNames,
                                      stock = stockNames,
                                      sex = sexNames )
  dimnames(repObj$L2_spx)    <- list( species = specNames,
                                      stock = stockNames,
                                      sex = sexNames )
  dimnames(repObj$vonK_spx)  <- list( species = specNames,
                                      stock = stockNames,
                                      sex = sexNames )




  return(repObj)
}

# savePlots()
savePlots <- function(  fitObj = reports,
                        useRep = "opt",
                        saveDir = "./Outputs/fits/" )
{
  # Create a plots sub directory in the saveDir
  saveDir <- file.path(saveDir,"plots")
  if(!dir.exists(saveDir))
    dir.create(saveDir)

  fYear <- fitObj$fYear
  lYear <- fitObj$lYear

  report   <- fitObj$repOpt

  specNames   <- fitObj$species
  stockNames  <- fitObj$stocks
  fleetNames  <- fitObj$gearLabs

  nS <- length(specNames)
  nP <- length(stockNames)
  nF <- length(fleetNames)
  
  graphics.off()

  # Plot recruitments
  png(  file.path(saveDir,"plotRspt.png"),
        width = 11, height = 8.5, units = "in", res = 300)
  plotRspt( repObj = report, initYear = fYear )
  dev.off()

  # Plot recruitment deviations
  png(  file.path(saveDir,"plotRsptDev.png"),
        width = 11, height = 8.5, units = "in", res = 300)
  plotRsptDev( repObj = report, initYear = fYear )
  dev.off()

  # Plot spawning biomass with catch 
  png(  file.path(saveDir,"plotSBspt.png"),
        width = 11, height = 8.5, units = "in", res = 300)
  plotSBspt( repObj = report, initYear = fYear )
  dev.off()

  # Plot probLenAge 
  png(  file.path(saveDir,"plotProbLenAge.png"),
        width = 11, height = 8.5, units = "in", res = 300)
  plotProbLenAge_sp( repObj = report)
  dev.off()

  # Plot indices
  png( file.path(saveDir,"plotStdzedIndices.png"),
        width = 11, height = 8.5, units = "in", res = 300)
  plotIspft(  repObj = report,
              fYear = fYear, lYear = lYear,
              sIdx = 1:nS, pIdx = 1:nP,
              fIdx = 1:nF )
  dev.off()

  png(  file.path(saveDir,"plotYieldCurves.png"),
        width = 11, height = 8.5, units = "in", res = 300 )
  plotYeqF( repObj = report,
            sIdx = 1:nS, pIdx = 1:nP )
  dev.off()

  # Plot catch fits
  png(  file.path(saveDir,"plotCatchFit.png"),
        width = 11, height = 8.5, units = "in", res = 300)
  plotCatchFit_spt( repObj = report, initYear = fYear )
  dev.off()

  # Plot Fs
  png(  file.path(saveDir,"plotFspft.png"),
        width = 11, height = 8.5, units = "in", res = 300)
  plotFspft( repObj = report, initYear = fYear )
  dev.off()

  # Plot Stock recruit
  png(  file.path(saveDir,"plotSRsp.png"),
        width = 11, height = 8.5, units = "in", res = 300)
  plotSRsp( repObj = report, initYear = fYear )
  dev.off()

  for( sIdx in 1:nS )
  {
    # Make a directory to hold species specific plots
    specDir <- specNames[sIdx]
    specPath <- file.path(saveDir,specDir)
    if(!dir.exists(specPath))
      dir.create(specPath)

    fileName <- paste("plotAgeLenResids_",specDir,".png",sep = "")
    png(  file.path(specPath,fileName),
        width = 11, height = 8.5, units = "in", res = 300)
    plotHeatmapAgeLenResids(  repObj = report, 
                              sIdx = sIdx, pIdx = 1:nP,
                              fIdx = 1:nF )
    dev.off()  

    fileName <- paste("plotAgeResids_",specDir,".png",sep = "")
    png(  file.path(specPath,fileName),
        width = 11, height = 8.5, units = "in", res = 300)
    plotCompResids( repObj = report, 
                    sIdx = sIdx, pIdx = 1:nP,
                    fIdx = 1:nF, comps = "age" )
    dev.off()  

    fileName <- paste("plotLenResids_",specDir,".png",sep = "")
    png(  file.path(specPath,fileName),
        width = 11, height = 8.5, units = "in", res = 300)
    plotCompResids( repObj = report, 
                    sIdx = sIdx, pIdx = 1:nP,
                    fIdx = 1:nF, comps = "length" )
    dev.off()  

    # Plot probLenAge 
    png(  file.path(specPath,"plotLenAtAgeDist.png"),
          width = 11, height = 8.5, units = "in", res = 300)
    plotHeatmapProbLenAge( repObj = report, 
                            sIdx = sIdx, pIdx = 1:nP)
    dev.off()

    for( pIdx in 1:nP)
    {
      # Make a directory to hold stock specific plots
      stockDir <- stockNames[pIdx]
      stockPath <- file.path(specPath,stockDir)
      if(!dir.exists(stockPath))
        dir.create(stockPath)

      # Plot Index fits to vuln bio
      png(  file.path(stockPath,"plotIdxFits.png"),
            width = 8.5, height = 11, units = "in", res = 300)
      plotIdxFits(  repObj = report, initYear = fYear,
                    sIdx = sIdx, pIdx = pIdx )
      dev.off()

      # Plot Index residuals
      png(  file.path(stockPath,"plotIdxResids.png"),
            width = 8.5, height = 11, units = "in", res = 300)
      plotIdxResids(  repObj = report, initYear = fYear,
                      sIdx = sIdx, pIdx = pIdx )
      dev.off()

      # Plot selectivity at age
      png(  file.path(stockPath,"plotSelAge.png"),
            width = 8.5, height = 11, units = "in", res = 300)
      plotSelAge( repObj = report,
                  sIdx = sIdx, pIdx = pIdx )
      dev.off()

      # Plot selectivity at length
      png(  file.path(stockPath,"plotSelLen.png"),
            width = 8.5, height = 11, units = "in", res = 300)
      plotSelLen( repObj = report,
                  sIdx = sIdx, pIdx = pIdx )
      dev.off()

      # # Plot maturity at length
      # png(  file.path(stockPath,"plotMatLen.png"),
      #       width = 11, height = 8.5, units = "in", res = 300)
      # plotMatLength( repObj = report,
      #                sIdx = sIdx, pIdx = pIdx )
      # dev.off()

      # Plot maturity at Age
      png(  file.path(stockPath,"plotMatAge.png"),
            width = 11, height = 8.5, units = "in", res = 300)
      plotMatAge( repObj = report,
                  sIdx = sIdx, pIdx = pIdx )
      dev.off()


      # Plot length-at-age
      png(  file.path(stockPath,"plotLenAge.png"),
        width = 11, height = 8.5, units = "in", res = 300)
      plotHeatmapProbLenAge( repObj = report, 
                              sIdx = sIdx, pIdx = pIdx)

      dev.off()

      # Plot weight at Age
      png(  file.path(stockPath,"plotWtAge.png"),
            width = 11, height = 8.5, units = "in", res = 300)
      plotWtAge( repObj = report,
                  sIdx = sIdx, pIdx = pIdx )
      dev.off()

      # Plot age comp fits
      plotCompFitYrs( repObj = report,
                      initYear = fYear,
                      sIdx = sIdx, pIdx = pIdx,
                      comps = "age", save = TRUE,
                      savePath = file.path(stockPath,"plotFitAgeYrs") )

      # Plot age comp fits (average)
      png(  file.path(stockPath,"plotFitAgeAvg.png"),
            width = 8.5, height = 11, units = "in", res = 300)
      plotCompFitAvg( repObj = report,
                      initYear = fYear,
                      sIdx = sIdx, pIdx = pIdx,
                      comps = "age" )
      dev.off()

      # Plot length comp fits
      plotCompFitYrs( repObj = report,
                      initYear = fYear,
                      sIdx = sIdx, pIdx = pIdx,
                      comps = "length", save = TRUE,
                      savePath = file.path(stockPath,"plotFitLenYrs") )

      # Plot length comp fits (average)
      png(  file.path(stockPath,"plotFitLenAvg.png"),
            width = 8.5, height = 11, units = "in", res = 300)
      plotCompFitAvg( repObj = report,
                      initYear = fYear,
                      sIdx = sIdx, pIdx = pIdx,
                      comps = "length" )
      dev.off()


      # Plot Catch fits by gear
      png(  file.path(stockPath,"plotCatchFit_ft.png"),
            width = 8.5, height = 11, units = "in", res = 300)
      plotCatchFit_ft(  repObj = report,
                        initYear = fYear,
                        sIdx = sIdx, pIdx = pIdx )
      dev.off()

      png(  file.path(stockPath,"plotStockRecruit.png"),
            width = 8.5, height = 8.5, units = "in", res = 300)
      plotSR( repObj = report,
              sIdx = sIdx, pIdx = pIdx )
      dev.off()

      # Plot spawning biomass with catch and indices
      png(  file.path(stockPath,"plotSBt.png"),
            width = 11, height = 8.5, units = "in", res = 300)
      plotSBt( repObj = report, initYear = fYear,
                  sIdx = sIdx, pIdx = pIdx )
      dev.off()

      for( fIdx in 1:nF )
      {
        if(!dir.exists(file.path(stockPath,"resids")) )
          dir.create(file.path(stockPath,"resids"))
        
        fleetID <- fleetNames[fIdx]
        fileName <- paste("plotAgeLenResids_",fleetID,".png",sep = "")
        png(  file.path(stockPath,"resids",fileName),
            width = 11, height = 8.5, units = "in", res = 300)
        plotHeatmapAgeLenResids(  repObj = report, 
                                  sIdx = 1:nS, pIdx = 1:nP,
                                  fIdx = fIdx )
        dev.off()  

        fileName <- paste("plotAgeResids_",fleetID,".png",sep = "")
        png(  file.path(stockPath,"resids",fileName),
            width = 11, height = 8.5, units = "in", res = 300)
        plotCompResids( repObj = report, 
                        sIdx = 1:nS, pIdx = 1:nP,
                        fIdx = fIdx, comps = "age" )
        dev.off()  

        fileName <- paste("plotLenResids_",fleetID,".png",sep = "")
        png(  file.path(stockPath,"resids",fileName),
            width = 11, height = 8.5, units = "in", res = 300)
        plotCompResids( repObj = report, 
                        sIdx = 1:nS, pIdx = 1:nP,
                        fIdx = fIdx, comps = "length" )
        dev.off()  
      }

    }
  }


} # END savePlots()

# Check NAs/NaNs in lists
.checkNA <- function( listEntry )
{
  x <- any(is.na(listEntry))
  x
}

# Check for Infs in lists
.checkFinite <- function( listEntry )
{
  x <- any(!is.finite(listEntry))
  x
}

.checkNaN <- function( listEntry )
{
  x <- any(is.nan(listEntry))
  x
}
