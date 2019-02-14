# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
#
# SCALfuns.R
#
# Functions for hierSCAL TMB model.
# 
# To Do:
#   4. Post-processing for report objects: add dimension names
#       for use in plotting code
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
  # Save plots
  # Choose whether we're using initial parameters
  # or the fixed effects (extend to RE model later)
  if( controlList$ctrl$opt )
    rep <- "FE"
  else rep <- "init"

  savePlots(  fitObj = reports,
              useRep = rep,
              saveDir = path  )

  # Copy control file to sim folder for posterity
  file.copy(from=ctlFile,to=file.path(path,"fitCtlFile.txt"))
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

  # rerun savePlots
  savePlots(  fitObj = reports,
              useRep = rep,
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

  # Get model dimensions
  nA_s      <- dataObj$nA_s
  minA_s    <- dataObj$minA_s
  A1_s      <- dataObj$A1_s
  A2_s      <- dataObj$A2_s
  nL_s      <- dataObj$nL_s
  minL_s    <- dataObj$minL_s
  fYear     <- dataObj$fYearData
  lYear     <- dataObj$lYearData
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
  ALFreq_spalftx <- makeALFreq( ALFreqList = ALfreq,
                                years = years,
                                gears = useFleets,
                                maxA = max(nA_s), maxL = max(nL_s) )

  # Combine sexes for now
  ALFreq_spalft <- ALFreq_spalftx[,,,,,,1] + ALFreq_spalftx[,,,,,,2]
  ALFreq_spalft <- ALFreq_spalft[useSpecies,useStocks,,,useFleets,yrChar,drop = FALSE]


  # Load age and length compositions
  load(file.path("./Data",dataObj$ageData))
  load(file.path("./Data",dataObj$lenData))

  # Create compositional arrays
  age_aspft <- makeCompsArray(  compList = ageComps,
                                plusGroups = nA_s,
                                minX = minA_s,
                                collapseComm = FALSE,
                                fleetIDs = useFleets,
                                combineSex = TRUE,
                                years = fYear:lYear,
                                xName = "ages",
                                minSampSize = dataObj$minAgeSampSize )
  age_aspft <- age_aspft[,useSpecies,useStocks,useFleets,, drop = FALSE]


  len_lspft <- makeCompsArray(  compList = lenComps,
                                plusGroups = nL_s,
                                minX = minL_s,
                                collapseComm = FALSE,
                                fleetIDs = useFleets,
                                combineSex = TRUE,
                                years = fYear:lYear,
                                xName = "length",
                                minSampSize = dataObj$minLenSampSize )
  len_lspft <- len_lspft[,useSpecies,useStocks,useFleets,, drop = FALSE]

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
        if( any(age_aspft[1,sIdx,pIdx,fIdx,yrChar] >= 0) &  dataObj$ageLikeWt > 0 )
          ageYrs <- which(age_aspft[1,sIdx,pIdx,fIdx,yrChar] >= 0)
        # And if length observations are included, add them
        if( any(len_lspft[1,sIdx,pIdx,fIdx,yrChar] >= 0) &  dataObj$lenLikeWt > 0 )
          lenYrs <- which(len_lspft[1,sIdx,pIdx,fIdx,yrChar] >= 0)

        # Take union of age/len years, then add number of years to 
        # estimate total number of selectivity deviations
        nSelDevs <- nSelDevs + length( union( lenYrs, ageYrs ) )
      } 

  # Now make a vector to switch time varying selectivity
  # on and off
  tvSelFleets <- rep(0,nF)
  names(tvSelFleets) <- useFleets
  tvSelFleets[ useFleets %in% hypoObj$tvSelFleets ] <- 1 


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
  

  # Count number of free F parameters
  nPosCatch <- length(which(catchDiscArrays$C_spft > 0) )

  # Initialised in a fished state (0 == unfished, 1 == fished)
  initFished_s     <- hypoObj$initFished[useSpecies]
  yFirstRecDev_s   <- hypoObj$yFirstRecDev[useSpecies]
  yLastRecDev_s    <- hypoObj$yLastRecDev[useSpecies]

  # Non-eq initialisation deviations
  nInitDevs <- sum( nP * (nA_s[useSpecIdx] * initFished_s ) )

  # Generate tFirstRecDev from yFirstRecDev and fYear
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
    xSel50_sf <- matrix(  c(    41, 39, 29, 31, 33, 33,
                                35, 35, 35, 35, 35, 35,
                                35, 35, 35, 35, 35, 35,
                                35, 35, 35, 35, 35, 35,
                                35, 35, 35, 35, 35, 35 ), 
                                nrow = length(allSpecies),
                                ncol = nF, byrow = TRUE )

    xSel50_sf <- xSel50_sf[useSpecIdx,]

    # length at SelStep_sf initial value
    xSelStep_sf <- matrix(  c(    2, 2, 2, 2, 2, 2,
                                  2, 2, 2, 2, 2, 2,
                                  2, 2, 2, 2, 2, 2,
                                  2, 2, 2, 2, 2, 2,
                                  2, 2, 2, 2, 2, 2 ), 
                                  nrow = length(allSpecies),
                                  ncol = nF, byrow = TRUE )

    xSelStep_sf <- xSelStep_sf[useSpecIdx,]
  }

  if( hypoObj$selX == "age")
  {
    # xSel50_sf initial value
    xSel50_sf <- matrix(  c(    8, 7, 5, 6, 5, 5,
                                5, 5, 5, 5, 5, 5,
                                5, 5, 5, 5, 5, 5,
                                5, 5, 5, 5, 5, 5,
                                5, 5, 5, 5, 5, 5 ), 
                                nrow = length(allSpecies),
                                ncol = nF, byrow = TRUE )

    xSel50_sf <- xSel50_sf[useSpecIdx,]

    # xSelStep_sf initial value
    xSelStep_sf <- matrix(  c(    3, 3, 2, 3, 2, 2,
                                  2, 2, 2, 2, 2, 2,
                                  2, 2, 2, 2, 2, 2,
                                  2, 2, 2, 2, 2, 2,
                                  2, 2, 2, 2, 2, 2 ), 
                                  nrow = length(allSpecies),
                                  ncol = nF, byrow = TRUE )

    xSelStep_sf <- xSelStep_sf[useSpecIdx,]
  }


  # Use useSpec, useStocks, useFleets and yrChar to
  # reduce the data arrays so that the fits are dynamic

  # Pull prior mean values from hypoObj to keep pars list tidier
  mh    <- hypoObj$hPriorMean
  sdh   <- hypoObj$hPriorSD
  mq    <- hypoObj$mq
  sdq   <- hypoObj$sdq
  mM    <- hypoObj$mM
  sdM   <- hypoObj$sdM

  pmxSel95_s <- hypoObj$pmxSel95_s[useSpecIdx]
  cvxSel95_f <- hypoObj$cvxSel95_f[useFleetsIdx]

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

  calcMeanLenAge <- function( sIdx = 1, ALK = ALFreq_spalft,
                              age = A1_s )
  {
    # Pull length at the given age
    lenAtAge <- ALK[sIdx,,age[sIdx],,,]
    lenAtAge <- apply(X = lenAtAge, FUN = sum, MARGIN = c(2))
    # Calculate number of observations
    nLenAtAge <- sum(lenAtAge,na.rm = T)
    # Multiply length bin by freq
    totLenAtAge <- sum((1:length(lenAtAge)) * lenAtAge )
    # Calc mean
    meanLenAge <- sum(totLenAtAge)/nLenAtAge

    meanLenAge
  }

  initL1_s <- sapply( X = 1:nS, 
                      FUN = calcMeanLenAge, 
                      ALK = ALFreq_spalft, 
                      age = A1_s[useSpecies] )
  initL2_s <- sapply( X = 1:nS, 
                      FUN = calcMeanLenAge, 
                      ALK = ALFreq_spalft, 
                      age = A2_s[useSpecies] )


  # Generate the data list
  data <- list( I_spft          = I_spft,
                C_spft          = C_spft,
                D_spft          = D_spft,
                ALK_spalft      = ALFreq_spalft,
                age_aspft       = age_aspft,
                len_lspft       = len_lspft,
                group_f         = as.integer(group_f),
                A_s             = as.integer(nA_s[useSpecies]),
                minA_s          = as.integer(minA_s[useSpecies]),
                L_s             = as.integer(nL_s[useSpecies]),
                minL_s          = as.integer(minL_s[useSpecies]),
                lenD_s          = rep(0,nS),
                swRinit_s       = as.integer(initFished_s[useSpecies]),
                parSwitch       = 0,
                calcIndex_spf   = calcIndex_spf,
                tvSelFleets     = tvSelFleets,
                ageLikeWt       = dataObj$ageLikeWt,
                lenLikeWt       = dataObj$lenLikeWt,
                tFirstRecDev_s  = as.integer(tFirstRecDev_s),
                tLastRecDev_s   = as.integer(tLastRecDev_s),
                minAgeProp      = dataObj$minAgeProp,
                minLenProp      = dataObj$minLenProp,
                minPAAL         = dataObj$minPAAL,
                matX            = hypoObj$matX,
                selX            = hypoObj$selX,
                A1_s            = A1_s[useSpecies],
                A2_s            = A2_s[useSpecies] )


  # Generate parameter list
  pars <- list( lnB0_sp           = log(sumCat_sp),
                logitSteep        = log((mh - .2)/(1 - mh)),
                lnM               = log(mM),
                lnL2step_s        = log(initL2_s - initL1_s),
                lnvonK_s          = log(rep(.2,nS)),
                lnL1_s            = log(initL1_s[useSpecIdx]),
                deltaL2_sp        = log(array(0,dim = c(nS,nP))),
                deltaVonK_sp      = log(array(0,dim = c(nS,nP))),
                sigmaLa_s         = sigmaLa_s,
                sigmaLb_s         = sigmaLb_s,
                LWa_s             = LWa_s,
                LWb_s             = LWb_s,
                xMat50_s          = xMat50,
                xMat95_s          = xMat95,
                lnq_spf           = array(0,dim =c(nS,nP,nF)),
                lntau_spf         = array(log(0.04),dim =c(nS,nP,nF)),
                lnxSel50_sf       = array(t(log(xSel50_sf)),dim = c(nS,nF) ),
                lnxSelStep_sf     = array(t(log(xSelStep_sf)),dim = c(nS,nF)),
                lnF_spft          = rep(-2,length = nPosCatch),
                lntauC_f          = rep(log(0.01),nF),
                lntauD_f          = rep(log(0.01),nF),
                muxSel50_sg       = array(3.4, dim = c(nS,3)),
                muxSel95_sg       = array(4.2, dim = c(nS,3)),
                sigmaxSel50_sg    = array(.2, dim = c(nS,3)),
                sigmaxSel95_sg    = array(.2, dim = c(nS,3)),
                lnqbarSyn_s       = rep(0,nS),
                lntauqSyn_s       = rep(log(.2),nS),
                lnqbarSyn         = log(1),
                lntauqSyn         = log(.2),
                mqSurveys         = 1,
                sdqSurveys        = 1,
                epsSteep_s        = rep(0,nS),
                epsM_s            = rep(0,nS),
                epsSteep_sp       = array(0, dim = c(nS,nP)),
                epsM_sp           = array(0, dim = c(nS,nP)),
                lnsigmah_s        = rep(log(sdh),nS),
                logit_muSteep     = log((mh - .2)/(1 - mh)),
                lnsigmah          = log(sdh),
                lnsigmaM_s        = rep( log(sdM), nS ),
                ln_muM            = log(mM),
                lnsigmaM          = log(sdM),
                IGatau_f          = tau2ObsIGa,
                IGbtau_f          = tau2ObsIGb,
                omegaR_vec        = rep(0, nRecDevs),
                omegaRinit_vec    = rep( 0, nInitDevs ),
                lnsigmaR_sp       = array(0, dim = c(nS,nP)),
                logitRCorr_chol   = rep(0, nS * nP),
                logitRgamma_sp    = array(0, dim = c(nS,nP)),
                epsxSel50_vec     = rep(0,nSelDevs),
                epsxSelStep_vec   = rep(0,nSelDevs),
                lnsigmaSel        = log(hypoObj$sigmaSelDevs),
                pmxSel95_sf       = matrix(pmxSel95_s,nrow = nS, ncol = nF, byrow = F),
                cvxSel95_f        = cvxSel95_f )

  # Generate map list - this will have to be updated
  # so that it's sensitive to useFleetsIdx
  qmap_spf <- array(  101 + 1:(nS*nP*nF),
                      dim = c(nS,nP,nF) )
  qmap_spf[calcIndex_spf == 0] <- NA
  # map out the observation error SD
  taumap_spf <- array(  1 + 1:(nS*nP*nF),
                        dim = c(nS,nP,nF) )
  taumap_spf[calcIndex_spf == 0] <- NA

  # Map Rinit
  RinitMap <- array(NA, dim = c(nS,nP))
  for( sIdx in 1:nS )
    if( initFished_s[sIdx] == 1 )
      RinitMap[sIdx,] <- (1:nP) + nP*(sIdx-1) + 200

  # Map selectivity at length
  selMap <- array(NA, dim = c(nS,nF))
  # Turn on selectivity estimation if asked for
  if( hypoObj$estSel )
    for(s in 1:nS)
      for(f in 1:nF)
        if( any(age_aspft[,s,,f,] > 0) | any(len_lspft[,s,,f,] > 0) )
          selMap[s,f] <- 300 + s * (nF - 1) + f

  map <- list(  
                lnB0_sp           = factor(array(NA,dim =c(nS,nP))),
                logitSteep        = factor(NA),
                lnM               = factor(NA),
                lnL2step_s        = factor(array(NA,dim =c(nS,1))),
                # lnvonK_s          = factor(array(NA,dim =c(nS,1))),
                lnL1_s            = factor(array(NA,dim =c(nS,1))),
                deltaL2_sp        = factor(matrix(NA,nrow = nS, ncol = nP)),
                deltaVonK_sp      = factor(matrix(NA,nrow = nS, ncol = nP)),
                LWa_s             = factor(rep(NA,nS)),
                LWb_s             = factor(rep(NA,nS)),
                xMat50_s          = factor(rep(NA,nS)),
                xMat95_s          = factor(rep(NA,nS)),
                lnq_spf           = factor(qmap_spf),
                lntau_spf         = factor(taumap_spf),
                lnxSel50_sf       = factor(selMap),
                lnxSelStep_sf     = factor(selMap + 100),
                # lnF_spft          = factor(rep(NA,nPosCatch)),
                lntauC_f          = factor(rep(NA,nF)),
                lntauD_f          = factor(rep(NA,nF)),
                sigmaLa_s         = factor(rep(NA,nS)),
                sigmaLb_s         = factor(rep(NA,nS)),
                muxSel50_sg       = factor(array(NA, dim = c(nS,3))),
                muxSel95_sg       = factor(array(NA, dim = c(nS,3))),
                sigmaxSel50_sg    = factor(array(NA, dim = c(nS,3))),
                sigmaxSel95_sg    = factor(array(NA, dim = c(nS,3))),
                lnqbarSyn_s       = factor(rep(NA,nS)),
                lntauqSyn_s       = factor(rep(NA,nS)),
                lnqbarSyn         = factor(NA),
                lntauqSyn         = factor(NA),
                mqSurveys         = factor(NA),
                sdqSurveys        = factor(NA),
                logitSteep_s      = factor(rep(NA,nS)),
                lnsigmah_s        = factor(rep(-NA,nS)),
                logitSteep        = factor(NA),
                lnsigmah          = factor(NA),
                lnM_s             = factor(rep(NA,nS)),
                lnsigmaM_s        = factor(rep(NA, nS ) ),
                ln_muM            = factor(NA),
                lnsigmaM          = factor(NA),
                IGatau_f          = factor(rep(NA,nF)),
                IGbtau_f          = factor(rep(NA,nF)),
                # omegaR_vec        = factor( rep( NA, nRecDevs ) ),
                omegaRinit_vec    = factor( rep( NA, nInitDevs ) ),
                lnsigmaR_sp       = factor(array(NA, dim = c(nS,nP)) ),
                logitRCorr_chol   = factor(rep(NA, nS * nP)),
                logitRgamma_sp    = factor(array(NA, dim = c(nS,nP))),
                lnsigmaSel        = factor(NA),
                pmxSel95_sf     = factor(array(NA,dim = c(nS,nF))),
                cvxSel95_f      = factor(rep(NA,nF)) ) 

  # Turn off tv sel deviations if not being used
  if(!hypoObj$tvSel)
  {
    map$epsxSel50_vec   <- factor(rep(NA,nSelDevs))
    map$epsxSelStep_vec <- factor(rep(NA,nSelDevs))
  }

  # Create a control list for the assessment model
  tmbCtrl <- list(  eval.max = ctrlObj$maxFunEval, 
                    iter.max = ctrlObj$maxIterations  )


  # Now create the AD fun object, assuming all fixed effects
  cat("Creating objective function object\n")
  objFE <- MakeADFun( data = data,
                      parameters = pars,
                      map = map, 
                      random = NULL, 
                      silent = ctrlObj$quiet )

  repInit <- objFE$report()
  repInit <- renameReportArrays( repObj = repInit, datObj = data )

  # Update names on report objects

  outList <- list(  repInit = repInit,
                    repFE = NULL,
                    sdrepFE = NULL,
                    sdrepFE.full = NULL,
                    repRE = NULL,
                    sdrepRE = NULL,
                    sdrepRE.full = NULL,
                    fYear = fYear, 
                    lYear = lYear,
                    gearLabs = useFleets,
                    species = useSpecies,
                    stocks = useStocks,
                    map = map,
                    data = data,
                    pars = pars )

  # Now try fitting the model
  if( ctrlObj$opt )
  {
    cat("Optimising with all fixed effects\n")
    fitFE <- try( nlminb (  start     = objFE$par,
                            objective = objFE$fn,
                            gradient  = objFE$gr,
                            control   = tmbCtrl ) )

    # May hang here if there are problems fitting the model
    if( class(fitFE) != "try-error")
    {
      outList$repFE         <- renameReportArrays( repObj = objFE$report(), datObj = data)
      outList$sdrepFE       <- sdreport(objFE)
      outList$heFE          <- objFE$he()
      outList$fitFE         <- fitFE
      outList$map           <- map

    }


  }

  return(outList)
} # END .runHierSCAL()

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
  lenNames    <- dimnames(datObj$len_lspft)[[1]]
  ageNames    <- dimnames(datObj$age_aspft)[[1]]


  # Ok, that's the data taken care of. There are still all the
  # new arrays that we created
  # Predicted data
  dimnames(repObj$I_spft_hat) <- dimnames(datObj$I_spft)
  dimnames(repObj$aDist_aspft_hat) <- dimnames(datObj$age_aspft)
  dimnames(repObj$lDist_lspft_hat) <- dimnames(datObj$len_lspft)
  # State arrays
  dimnames(repObj$B_aspt) <- dimnames(datObj$age_aspft)[c(1:3,5)]
  dimnames(repObj$N_aspt) <- dimnames(datObj$age_aspft)[c(1:3,5)]
  dimnames(repObj$B_spt) <- dimnames(datObj$age_aspft)[c(2:3,5)]
  dimnames(repObj$R_spt) <- dimnames(datObj$age_aspft)[c(2:3,5)]
  dimnames(repObj$SB_spt) <- dimnames(datObj$age_aspft)[c(2:3,5)]
  dimnames(repObj$Bv_spft) <- dimnames(datObj$age_aspft)[c(2:5)]
  dimnames(repObj$predC_spft) <- dimnames(datObj$age_aspft)[c(2:5)]
  dimnames(repObj$predCw_spft) <- dimnames(datObj$age_aspft)[c(2:5)]
  dimnames(repObj$C_aspft) <- dimnames(datObj$age_aspft)[c(1:5)]
  dimnames(repObj$Cw_aspft) <- dimnames(datObj$age_aspft)[c(1:5)]
  dimnames(repObj$F_aspft) <- dimnames(datObj$age_aspft)[c(1:5)]
  dimnames(repObj$F_spft) <- dimnames(datObj$age_aspft)[c(2:5)]
  dimnames(repObj$Z_aspt) <- dimnames(datObj$age_aspft)[c(1:3,5)]
  # Biological parameters
  dimnames(repObj$R0_sp) <- dimnames(datObj$age_aspft)[c(2:3)]
  dimnames(repObj$B0_sp) <- dimnames(datObj$age_aspft)[c(2:3)]
  dimnames(repObj$h_sp) <- dimnames(datObj$age_aspft)[c(2:3)]
  dimnames(repObj$M_sp) <- dimnames(datObj$age_aspft)[c(2:3)]
  dimnames(repObj$phi_sp) <- dimnames(datObj$age_aspft)[c(2:3)]
  dimnames(repObj$reca_sp) <- dimnames(datObj$age_aspft)[c(2:3)]
  dimnames(repObj$recb_sp) <- dimnames(datObj$age_aspft)[c(2:3)]


  # Observation models
  dimnames(repObj$q_spf)      <- dimnames(datObj$age_aspft)[c(2:4)]  
  dimnames(repObj$tau_spf)    <- dimnames(datObj$age_aspft)[c(2:4)]  
  dimnames(repObj$sel_lfsp)   <- list(  len = lenNames, 
                                        fleets = gearNames,
                                        species = specNames,
                                        stocks = stockNames )
  dimnames(repObj$sel_afsp)   <- list(  age = ageNames, 
                                        fleets = gearNames,
                                        species = specNames,
                                        stocks = stockNames )



  return(repObj)
}

# savePlots()
savePlots <- function(  fitObj = reports,
                        useRep = "FE",
                        saveDir = "./Outputs/fits/" )
{
  # Create a plots sub directory in the saveDir
  saveDir <- file.path(saveDir,"plots")
  if(!dir.exists(saveDir))
    dir.create(saveDir)

  fYear <- fitObj$fYear
  lYear <- fitObj$lYear

  if(useRep == "init")
    report <- fitObj$repInit

  if( useRep == "FE" )
    report   <- fitObj$repFE

  nS <- report$nS
  nP <- report$nP

  specNames <- dimnames(report$R_spt)[[1]]
  stockNames <- dimnames(report$R_spt)[[2]]
  
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

      # Plot maturity at length
      png(  file.path(stockPath,"plotMatLen.png"),
            width = 11, height = 8.5, units = "in", res = 300)
      plotMatLength( repObj = report,
                     sIdx = sIdx, pIdx = pIdx )
      dev.off()

      # Plot maturity at Age
      png(  file.path(stockPath,"plotMatAge.png"),
            width = 11, height = 8.5, units = "in", res = 300)
      plotMatAge( repObj = report,
                  sIdx = sIdx, pIdx = pIdx )
      dev.off()

      # Plot length at Age
      png(  file.path(stockPath,"plotLenAge.png"),
            width = 11, height = 8.5, units = "in", res = 300)
      plotLenAge( repObj = report,
                  sIdx = sIdx, pIdx = pIdx )
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

    }
  }


} # END savePlots()



