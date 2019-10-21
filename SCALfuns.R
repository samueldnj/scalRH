# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
#
# SCALfuns.R
#
# Functions for fitting hierSCAL TMB model to DERPA data.
#
#
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

# fitHierSCAL()
# Wrapper function to fit the hierSCAL model under a data scenario 
# and model hypothesis, both defined in the ctlFile. If fit is succesful
# then model outputs are saved to a folder in the ./Outputs/fits/ directory
# inputs:   ctlFile=character with the name/path of the control file
#           folder=optional character name of output folder 
#                   ie saves to ./Outputs/fits/<folder>
#           quiet=logical to determine if optimisation progress and 
#                 other messages are returned
#           cplx=optional character vector to change the species
#                 complex. Meant for batching, but can be used
#                 from the console as well.
# ouputs:   NULL
# usage:    from the console to run the procedure
fitHierSCAL <- function ( ctlFile = "fitCtlFile.txt", 
                          folder=NULL, 
                          quiet=TRUE,
                          cplx = NULL,
                          groupFolder = "." )
{ 
  # read in control file
  controlTable <- .readParFile ( ctlFile )

  # Replace complex if provided
  if( !is.null(cplx) )
  {
    specRow <- which(controlTable$parameter == "data$species")

    cplx <- paste(cplx,collapse = "','")
    cplx <- paste("c('",cplx,"')",sep = "")

    controlTable[specRow,"value"] <- cplx
  }
  # Create control list
  controlList <- .createList  ( controlTable )

  # Load hierSCAL object
  dyn.load(dynlib("hierSCAL"))

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
  path <- here::here("Outputs","fits",groupFolder,folder)
  dir.create ( path )
  cat( "\nMSG (saveFit) Created assessment fit folder ",folder,"in ./Outputs/fits/.\n" )
  
  # Save path so we can access it later
  reports$path <- path

  # Save reports object
  save(reports,file = file.path(path,paste(folder,".RData",sep="")))

  # Make html fit report
  makeFitReport( fitID = folder, groupFolder = groupFolder )

  # Create a quick to read info file for the fit folder
  .makeInfoFile(reports)

  # Copy control file to sim folder for posterity
  cat(  "# fitCtlFile.txt, written to ", folder, "on ", Sys.time(),"\n", sep = "", 
        file = file.path(path,"fitCtlFile.txt"))
  write.table(  controlTable, 
                file = file.path(path,"fitCtlFile.txt"),
                row.names = FALSE,
                quote = FALSE, qmethod = "double",
                append = TRUE )
  cat(  "# <End File>", sep = "", 
        file = file.path(path,"fitCtlFile.txt"),
        append = TRUE)
  # Copy model files, so we can recreate report objects later.
  file.copy(from="hierSCAL.cpp",to=file.path(path,"hierSCAL.cpp"))
  file.copy(from="hierSCAL.so",to=file.path(path,"hierSCAL.so"))
  file.copy(from="hierSCAL.o",to=file.path(path,"hierSCAL.o"))

  dyn.unload(dynlib("hierSCAL"))

  # Save plots?
  if(controlList$ctrl$plots)
    savePlots(  fitObj = reports,
                saveDir = path  )

  beepr::beep("complete")
  
  # Done
} # END fitHierSCAL()

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
  nL_s      <- dataObj$nL_s
  minL_s    <- dataObj$minL_s
  fYear     <- dataObj$fYearData
  lYear     <- dataObj$lYearData
  fYearIdx  <- dataObj$fYearIdx
  lYearIdx  <- dataObj$lYearIdx
  nT        <- lYear - fYear + 1
  years     <- fYear:lYear
  yrChar    <- as.character(years)


  commNames_g <- dataObj$commNames_g
  survNames_g <- dataObj$survNames_g 

  # Now load the fleetIDs
  loadStockSpecNameLists()

  collapseSyn  <- dataObj$collapseSyn

  # Track used fleets
  allFleets     <- c( names(commFleetYrRange), names(surveyIDs) )
  useFleets     <- c( commNames_g, survNames_g)
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

  if( collapseSyn )
  {
    notSynSurv  <- survNames_g[!grepl("Syn",survNames_g)]
    notSynAll   <- allFleets[!grepl("Syn",allFleets)] 

    newSurvNames  <- c(notSynSurv,"Syn")
    useFleets     <- c(commNames_g, newSurvNames)
    allFleets     <- c(notSynAll,"Syn")

    useFleetIdx   <- which( allFleets %in% useFleets )
    nF            <- length(useFleets)
    totF          <- length(allFleets)
  }

  # Load data
  message("Loading data for assessment\n")

  # Load index data
  idxDataFiles <- dataObj$idxFiles
  idxDataFilePaths <- file.path("./Data/Proc",idxDataFiles )
  for( i in 1:length(idxDataFilePaths))
    load( idxDataFilePaths[i] )

  # Update makeIndexArray to use the commFleetYrRange vector
  I_spft <- makeIndexArray( relBio = relBioList,
                            commCPUE = commCPUEList,
                            nP = length(allStocks), years = fYear:lYear,
                            collapseComm = FALSE,
                            scaleComm = dataObj$commCPUEscalar )

  if( collapseSyn )
  {
    # collapseSyn assumes that the synoptic 
    # surveys are just HS, QCS and WCVI, from
    # N to S, and that Synoptic surveys are 
    # at the end of the fleetidx
    newI_spft <- I_spft

    whichSyn  <- which(grepl("Syn",dimnames(I_spft)[[3]]))
    notSyn    <- which(!grepl("Syn",dimnames(I_spft)[[3]]))

    notSynNames <- dimnames(I_spft)[[3]][notSyn]

    newI_spft <- I_spft[,,c(notSynNames,"HSSyn"),]
    dimnames(newI_spft)[[3]] <- c(notSynNames,"Syn")

    newI_spft[,"QCS","Syn",] <- I_spft[,"QCS","QCSSyn",]
    newI_spft[,"WCVI","Syn",] <- I_spft[,"WCVI","WCVISyn",]

    I_spft <- newI_spft
  }

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
  commCatch <- read.csv(  file.path("./Data/Raw/",dataObj$catchData["comm"]), 
                          header = TRUE,
                          stringsAsFactors = FALSE )

  survCatchPath  <- file.path("./Data/Proc",dataObj$catchData["survey"])
  load(  file = survCatchPath )

  # Update makeCatchDiscArrays to use the commFleetYrRange vector
  catchDiscArrays <- makeCatchDiscArrays( commData = commCatch,
                                          survData = surveyCatch,
                                          stocks = stocksCommBio,
                                          speciesCodes = specCodes,
                                          years = fYear:lYear,
                                          modernYear = 1996,
                                          nF = length(c(commNames_g,survNames_g)), 
                                          collapseComm = FALSE,
                                          fleetIDs = c(commNames_g,survNames_g) )
  # Extract catch and discards
  C_spft <- catchDiscArrays$C_spft[useSpecies,useStocks,,, drop = FALSE]
  D_spft <- catchDiscArrays$D_spft[useSpecies,useStocks,,, drop = FALSE]

  if( !dataObj$modelDisc )
    C_spft <- C_spft + D_spft

  if( collapseSyn )
  {
    # collapseSyn assumes that the synoptic 
    # surveys are just HS, QCS and WCVI, from
    # N to S, and that Synoptic surveys are 
    # at the end of the fleetidx
    newC_spft <- C_spft

    whichSyn  <- which(grepl("Syn",dimnames(C_spft)[[3]]))
    notSyn    <- which(!grepl("Syn",dimnames(C_spft)[[3]]))

    notSynNames <- dimnames(C_spft)[[3]][notSyn]

    newC_spft <- C_spft[,,c(notSynNames,"HSSyn"),]
    dimnames(newC_spft)[[3]] <- c(notSynNames,"Syn")

    newC_spft[,"QCS","Syn",] <- C_spft[,"QCS","QCSSyn",]
    newC_spft[,"WCVI","Syn",] <- C_spft[,"WCVI","WCVISyn",]

    C_spft <- newC_spft

    newD_spft <- D_spft

    whichSyn  <- which(grepl("Syn",dimnames(D_spft)[[3]]))
    notSyn    <- which(!grepl("Syn",dimnames(D_spft)[[3]]))

    notSynNames <- dimnames(C_spft)[[3]][notSyn]

    newD_spft <- D_spft[,,c(notSynNames,"HSSyn"),]
    dimnames(newD_spft)[[3]] <- c(notSynNames,"Syn")

    newD_spft[,"QCS","Syn",] <- D_spft[,"QCS","QCSSyn",]
    newD_spft[,"WCVI","Syn",] <- D_spft[,"WCVI","WCVISyn",]

    D_spft <- newD_spft
  }



  # Sum catch for initial B0 estimate
  sumCat_sp <- apply( X = C_spft, FUN = sum, MARGIN = c(1,2) )


  # Create an array of mean "effort" in a stock
  # area by fleet and time

  E_spft  <- C_spft / I_spft
  E_spft[E_spft < 0] <- NA
  E_pft   <- apply( X = E_spft, FUN = mean, MARGIN = c(2,3,4),
                    na.rm = TRUE)
  E_pft[is.na(E_pft)] <- -1

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
  growthFilePaths <- file.path("./Data/Proc",growthFiles )
  for( i in 1:length(growthFiles))
    load(growthFilePaths[i])


  # Growth parameters come from vonBFits
  L1_s        <- vonBFits$L1_s
  L2_spx      <- vonBFits$L2_spx
  A1_s        <- vonBFits$A1_s
  A2_s        <- vonBFits$A2_s
  sigA_s      <- vonBFits$sigA_s 
  sigB_s      <- vonBFits$sigB_s 
  vonK_spx    <- vonBFits$VonK_spx

  # Calculate L2step_spx
  L2step_spx <- L2_spx

  for( sIdx in 1:length(L1_s) )
    L2step_spx[sIdx,,] <- L2_spx[sIdx,,] - L1_s[sIdx] 



  # Make the age-length key - this might be useful,
  # or we can include it so we can integrate a growth
  # model
  # Aggregate ages into plus groups later.
  ALFreq_spalftx <- makeALFreq( ALFreqList = ALfreq,
                                years = years,
                                gears = c(commNames_g,survNames_g),
                                maxA = max(nA_s), maxL = max(nL_s) )

  if( collapseSyn )
  {
    # collapseSyn assumes that the synoptic 
    # surveys are just HS, QCS and WCVI, from
    # N to S, and that Synoptic surveys are 
    # at the end of the fleetidx
    newALFreq_spalftx <- ALFreq_spalftx

    whichSyn  <- which(grepl("Syn",dimnames(ALFreq_spalftx)[[5]]))
    notSyn    <- which(!grepl("Syn",dimnames(ALFreq_spalftx)[[5]]))

    notSynNames <- dimnames(newALFreq_spalftx)[[5]][notSyn]

    newALFreq_spalftx <- ALFreq_spalftx[,,,,c(notSynNames,"HSSyn"),,]
    dimnames(newALFreq_spalftx)[[5]] <- c(notSynNames,"Syn")
    
    newALFreq_spalftx[,"QCS",,,"Syn",,] <- ALFreq_spalftx[,"QCS",,,"QCSSyn",,]
    newALFreq_spalftx[,"WCVI",,,"Syn",,] <- ALFreq_spalftx[,"WCVI",,,"WCVISyn",,]

    ALFreq_spalftx <- newALFreq_spalftx
  }

  # Pare down to specific fleets for growth
  # CAAL data - Lee et al 2019+ (shared privately)
  ALfleetNames    <- dimnames(ALFreq_spalftx)[[5]]
  growthFleetIdx  <- which(ALfleetNames %in% dataObj$growthFleets)

  growthYears     <- dataObj$growthYears[1]:dataObj$growthYears[2]
  growthYrIdx     <- growthYears - fYear + 1

  
  # Set all observations outside those fleets to 0
  ALFreq_spalftx[,,,,-growthFleetIdx,-growthYrIdx,] <- 0



  # Load age and length compositions
  load(file.path("./Data/Proc",dataObj$ageData))
  load(file.path("./Data/Proc",dataObj$lenData))


  # Create compositional arrays
  age_aspftx <- makeCompsArray( compList = ageComps,
                                plusGroups = nA_s,
                                minX = minA_s,
                                collapseComm = FALSE,
                                fleetIDs = c(commNames_g,survNames_g),
                                years = fYear:lYear,
                                xName = "ages",
                                minSampSize = dataObj$minAgeSampSize )


  len_lspftx <- makeCompsArray( compList = lenComps,
                                plusGroups = nL_s,
                                minX = minL_s,
                                binWidth = dataObj$lenBinWidth,
                                collapseComm = FALSE,
                                fleetIDs = c(commNames_g,survNames_g),
                                years = fYear:lYear,
                                xName = "length",
                                minSampSize = dataObj$minLenSampSize )

  if( collapseSyn )
  {
    # collapseSyn assumes that the synoptic 
    # surveys are just HS, QCS and WCVI, from
    # N to S, and that Synoptic surveys are 
    # at the end of the fleetidx
    newage_aspftx <- age_aspftx

    whichSyn  <- which(grepl("Syn",dimnames(age_aspftx)[[4]]))
    notSyn    <- which(!grepl("Syn",dimnames(age_aspftx)[[4]]))

    notSynNames <- dimnames(newage_aspftx)[[4]][notSyn]

    newage_aspftx <- age_aspftx[,,,c(notSynNames,"HSSyn"),,]
    dimnames(newage_aspftx)[[4]] <- c(notSynNames,"Syn")
    
    newage_aspftx[,,"QCS","Syn",,] <- age_aspftx[,,"QCS","QCSSyn",,]
    newage_aspftx[,,"WCVI","Syn",,] <- age_aspftx[,,"WCVI","WCVISyn",,]

    age_aspftx <- newage_aspftx

    newlen_lspftx <- len_lspftx

    whichSyn  <- which(grepl("Syn",dimnames(len_lspftx)[[4]]))
    notSyn    <- which(!grepl("Syn",dimnames(len_lspftx)[[4]]))

    notSynNames <- dimnames(newlen_lspftx)[[4]][notSyn]

    newlen_lspftx <- len_lspftx[,,,c(notSynNames,"HSSyn"),,]
    dimnames(newlen_lspftx)[[4]] <- c(notSynNames,"Syn")
    
    newlen_lspftx[,,"QCS","Syn",,] <- len_lspftx[,,"QCS","QCSSyn",,]
    newlen_lspftx[,,"WCVI","Syn",,] <- len_lspftx[,,"WCVI","WCVISyn",,]

    len_lspftx <- newlen_lspftx
  }

  age_aspftx <- age_aspftx[,useSpecies,useStocks,useFleets,,1:2, drop = FALSE]
  len_lspftx <- len_lspftx[,useSpecies,useStocks,useFleets,,, drop = FALSE]
  
  lenBinMids  <- as.numeric(dimnames(len_lspftx)[[1]])
  maxLenBin_s <- ceiling( nL_s[useSpecies] / dataObj$lenBinWidth )
  minLenBin_s <- ceiling( minL_s[useSpecies] / dataObj$lenBinWidth )


  nX <- 2

  # Combine sexes if not sex structured
  if( dataObj$sexStructure == "Combined")
  {
    ALFreq_spalftx <- ALFreq_spalftx[,,,,,,1,drop = FALSE] + ALFreq_spalftx[,,,,,,2,drop = FALSE]
    dimnames(ALFreq_spalftx)[[7]] <- "both"
    age_aspftx <- age_aspftx[,,,,,1,drop = FALSE] + age_aspftx[,,,,,2,drop = FALSE] + age_aspftx[,,,,,3,drop = FALSE] 
    dimnames(age_aspftx)[[6]] <- "both"
    len_lspftx <- len_lspftx[,,,,,1,drop = FALSE] + len_lspftx[,,,,,2,drop = FALSE] + len_lspftx[,,,,,3,drop = FALSE] 
    dimnames(len_lspftx)[[6]] <- "both"

    nX <- 1
  }

  # Pare down to female only
  if( dataObj$sexStructure == "Female")
  {
    ALFreq_spalftx <- ALFreq_spalftx[,,,,,,2,drop = FALSE]
    age_aspftx <- age_aspftx[,,,,,2,drop = FALSE] 
    len_lspftx <- len_lspftx[,,,,,2,drop = FALSE]
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

  # And the catchability learning rate
  logitqFleets <- rep(0,nF)
  names(logitqFleets) <- useFleets
  logitqFleets[ useFleets %in% names(hypoObj$tq50_vec) ] <- 1 

  # Expand tq50 and tq95 into multiples for stock and species
  tq50_vec <- c()
  tq95_vec <- c()
  # Loop and fill tq50 and tq95
  for( f in 1:sum(logitqFleets))
  {
    tq50_vec <- c(tq50_vec, rep(hypoObj$tq50_vec[f], nS * nP) )
    tq95_vec <- c(tq95_vec, rep(hypoObj$tq95_vec[f], nS * nP) )
  }

  # And the same for F regularisation
  regFfleets <- rep(0,nF)
  names(regFfleets) <- useFleets
  regFfleets[ useFleets %in% hypoObj$regFfleets ] <- 1 

  # Load maturity ogives
  matOgives <- read.csv(file.path("./Data/Proc",dataObj$matFiles))
  wtLen     <- read.csv(file.path("./Data/Proc",dataObj$wtLenFile))

  # Get growth and maturity info
  # Vectors to hold
  # Mat
  xMat50 <- numeric(length = nS)
  xMat95 <- numeric(length = nS)
  # Wt
  LWa_s <- numeric(length = nS)
  LWb_s <- numeric(length = nS)

  for( sIdx in 1:nS )
  {
    specName <- useSpecies[sIdx]
    if( hypoObj$matX == "length")
    {
      xMat50[sIdx] <- matOgives[matOgives$Species == specName,"l50"]
      xMat95[sIdx] <- matOgives[matOgives$Species == specName,"l95"]
    }
    if(hypoObj$matX == "age")
    {
      xMat50[sIdx] <- matOgives[matOgives$Species == specName,"a50"]
      xMat95[sIdx] <- matOgives[matOgives$Species == specName,"a95"]
    }

    LWa_s[sIdx] <- wtLen[wtLen$Species == specName, 'lw.a']
    LWb_s[sIdx] <- wtLen[wtLen$Species == specName, 'lw.b']
  }

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

  # Count commercial and fishery indep fleets
  nComm       <- length(dataObj$commNames_g)
  nSurv       <- length(dataObj$survNames_g)
  # Group fleets
  fleetGroups <- hypoObj$fleetGroups
  nGroups     <- length(fleetGroups)
  
  group_f     <- integer(length = nF)
  for( f in 1:nF )
    for( g in 1:nGroups)
      if( useFleets[f] %in% fleetGroups[[g]] )
      {
        group_f[f] <- g
        break
      }
  

  if( hypoObj$selX == "length")
  {

    # length at Sel50_sf initial value
    xSel50_sf <- matrix(  25, 
                          nrow = length(allSpecies),
                          ncol = length(allFleets), 
                          byrow = TRUE )


    xSel50_sf <- xSel50_sf[useSpecIdx,useFleetsIdx]

    # length at SelStep_sf initial value
    xSelStep_sf <- matrix(  3, 
                            nrow = length(allSpecies),
                            ncol = length(allFleets), 
                            byrow = TRUE )

    xSelStep_sf <- xSelStep_sf[useSpecIdx,useFleetsIdx,drop = FALSE]



    # Selectivity by group
    xSel50_sg <- matrix(  c( 40, 23, 28,
                             30, 28, 22,
                             33, 35, 29,
                             35, 29, 35,
                             41, 35, 35 ), 
                          nrow = length(allSpecies),
                          ncol = nGroups, 
                          byrow = TRUE )


    xSelStep_sg <- matrix(  c(6, 6, 6,
                              5, 4, 4,
                              5, 5, 5,
                              4, 5, 7,
                              7, 15, 15 ), 
                            nrow = length(allSpecies),
                            ncol = nGroups, 
                            byrow = TRUE )

    xSel50_sg   <- xSel50_sg[useSpecIdx,]
    xSelStep_sg <- xSelStep_sg[useSpecIdx,]

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

  # by fleet
  lnxSel50_sf <- array(log(xSel50_sf),dim = c(nS,nF))
  lnxSelStep_sf <- array(log(xSelStep_sf),dim = c(nS,nF))


  # by fleet group
  lnxSel50_sg <- array(log(xSel50_sg),dim = c(nS,nGroups))
  lnxSelStep_sg <- array(log(xSelStep_sg),dim = c(nS,nGroups))

  calcStockSelDevs_spf <- array(0, dim = c(nS,nP,nF)) 
  calcStockQDevs_spf <- array(0, dim = c(nS,nP,nF)) 

  # Map selectivity at length
  nStockSelDevs <- 0
  nStockqDevs <- 0
  # Make unique initially
  for(s in 1:nS)
    for( f in 1:nF)
    {
      # Count number of fleets in the group, skip if == 1
      gp <- group_f[f]
      nFleetsInGp <- length(fleetGroups[[gp]])
      if( allFleets[f] == "HSAss" ) next

      for(p in 1:nP)  
      {
        if( any(age_aspftx[,s,p,f,,] > 0) | any(len_lspftx[,s,p,f,,] > 0) )
        {
          nStockSelDevs <- nStockSelDevs + 1
          calcStockSelDevs_spf[s,p,f] <- 1
        }
        if( calcIndex_spf[s,p,f] == 1 )
        {
          nStockqDevs <- nStockqDevs + 1
          calcStockQDevs_spf[s,p,f] <- 1
        }
      }
    }


  # Use useSpec, useStocks, useFleets and yrChar to
  # reduce the data arrays so that the fits are repsonsive
  # to changing data scenarios and species complex structures

  # Prior mean values for hierarchical parameters
  # steepness, catchability, and mortality
  hBetaPrior  <- hypoObj$hBetaPrior
  mq_f        <- hypoObj$mq_f
  sdq_f       <- hypoObj$sdq_f
  mM          <- hypoObj$mM
  sdM         <- hypoObj$sdM

  # Calculate prior mean h
  mh <- hBetaPrior[1] / sum(hBetaPrior)

  # Observation errors
  tau2ObsIGa  <- hypoObj$tau2ObsIGa[useFleetsIdx]
  tau2ObsMode <- hypoObj$tau2ObsPriorMode[useFleetsIdx]
  tau2ObsIGb  <- (tau2ObsIGa + 1) * tau2ObsMode

  # Catchability group level SD IG prior
  tau2qMode_f         <- (hypoObj$tauqMode_f)^2
  IGbtauq_f           <- (hypoObj$IGatauq_f + 1) * tau2qMode_f

  # Make shrinkage prior SD for 50% and step sel pars
  lnsigmaxSel50_sg    <- array(0, dim = c(nS,nGroups))
  lnsigmaxSelStep_sg    <- array(0, dim = c(nS,nGroups))

  for( s in 1:nS )
  {
    lnsigmaxSel50_sg[s,]    <- log(hypoObj$pmsigmaSel_g)
    lnsigmaxSelStep_sg[s,]  <- log(hypoObj$pmsigmaSel_g)
  }
  # IG Prior on recruitment variance
  sigma2R_mode        <- (hypoObj$IGsigmaRmode)^2
  IGbsigmaR           <- (hypoObj$IGasigmaR + 1) * sigma2R_mode
  # data <- .loadData( ctlList = obj )

  # Recruitment correlations
  IWmode <- diag(1,(nS * nP))
  IWnu   <- nS * nP + hypoObj$IWnu

  if( hypoObj$IWscale == "stockCorr" )
  {
    # Create a striped correlation matrix
    # that correlates the recruitments within stock areas
  }

  if( hypoObj$IWscale == "specCorr" )
  {
    # Create a striped correlation matrix
    # that correlates the recruitments within
    # species
  }

  if( hypoObj$IWscale == "distCorr" )
  {
    # Create a striped correlation matrix
    # that correlates the recruitments 
    # as a function of distance (GMRF style)
  }  


  # Generate the data list
  data <- list( I_spft                = I_spft,
                C_spft                = C_spft,
                D_spft                = D_spft,
                E_pft                 = E_pft,
                ALK_spalftx           = ALFreq_spalftx,
                age_aspftx            = age_aspftx,
                len_lspftx            = len_lspftx,
                lenBinMids_l          = lenBinMids,
                lenBinWidth           = dataObj$lenBinWidth,
                group_f               = as.integer(group_f - 1),
                A_s                   = as.integer(nA_s[useSpecies]),
                minA_s                = as.integer(minA_s[useSpecies]),
                L_s                   = maxLenBin_s,
                minL_s                = minLenBin_s,
                lenD_s                = rep(0,nS),
                nX                    = nX,
                swRinit_s             = as.integer(initFished_s[useSpecies]),
                parSwitch             = 0,
                calcIndex_spf         = calcIndex_spf,
                tvSelFleets           = tvSelFleets,
                tvqFleets             = tvqFleets,
                logitqFleets          = logitqFleets,
                regFfleets            = regFfleets,
                idxLikeWt_g           = dataObj$idxLikeWt_g,
                ageLikeWt_g           = dataObj$ageLikeWt_g,
                lenLikeWt_g           = dataObj$lenLikeWt_g,
                tFirstRecDev_s        = as.integer(tFirstRecDev_s),
                tLastRecDev_s         = as.integer(tLastRecDev_s),
                minAgeProp            = dataObj$minAgeProp,
                minLenProp            = dataObj$minLenProp,
                matX                  = hypoObj$matX,
                selX                  = hypoObj$selX,
                lenComps              = hypoObj$lenCompMethod,
                lambdaB0              = dataObj$lambdaB0,
                lambdaPropF           = dataObj$lambdaPropF,
                minTimeIdx_spf        = minTimeIdx_spf,
                nBaranovIter          = ctrlObj$nBaranovIter,
                lambdaBaranovStep     = ctrlObj$lambdaBaranovStep,
                calcFmethod           = ctrlObj$calcFmethod,
                A1_s                  = A1_s[useSpecIdx],
                A2_s                  = A2_s[useSpecIdx],
                postFitSR             = hypoObj$postFitSR,
                calcStockSelDevs_spf  = calcStockSelDevs_spf,
                calcStockQDevs_spf    = calcStockQDevs_spf,
                boundRecDevs          = hypoObj$boundRecDevs,
                recruitVariance       = hypoObj$recModel,
                recOption             = hypoObj$recOption,
                debugMode             = ctrlObj$debugMode,
                condMLEq              = hypoObj$condMLEq  )

  # Generate parameter list
  pars <- list( ## Leading biological pars ##
                lnB0_sp             = log(sumCat_sp),
                lnRbar_sp           = array(10, dim = c(nS,nP)),
                logitSteep          = log((mh - .2)/(1 - mh)),
                lnM                 = log(mM),
                lnL2step_spx        = log(L2step_spx[useSpecIdx,useStockIdx,,drop = FALSE]),
                lnvonK_spx          = log(vonK_spx[useSpecIdx,useStockIdx,,drop=FALSE]),
                lnL1_s              = log(L1_s[useSpecIdx]),
                # process error in growth model
                lnsigmaLa_s         = log(sigA_s[useSpecIdx]),
                sigmaLb_s           = sigB_s[useSpecIdx],
                # L-W conversion
                LWa_s               = LWa_s,
                LWb_s               = LWb_s,
                # Maturity
                xMat50_s            = xMat50,
                xMat95_s            = xMat95,
                ## Observation models ##
                # fleet catchability and obs idx SD
                lnq_sf              = array(0,dim = c(nS,nF)),
                lntq50_vec          = log(tq50_vec),
                lntq95_vec          = log(tq95_vec),
                lntauObs_spf        = array(0,dim = c(nS,nP,nF)),
                # Selectivity top level means and deviations for fleets/stocks
                lnxSel50_sg         = lnxSel50_sg,
                lnxSelStep_sg       = lnxSelStep_sg,
                epsxSel50spf_vec    = rep(0,nStockSelDevs),
                epsxSelStepspf_vec  = rep(0,nStockSelDevs),
                # discards obs SD
                lntauD_f            = rep(log(0.01),nF),
                ## Multilevel priors ##
                # Selectivity SDs
                lnsigmaxSel50_sg    = lnsigmaxSel50_sg,
                lnsigmaxSelStep_sg  = lnsigmaxSelStep_sg,
                pmsigmaSel_g        = hypoObj$pmsigmaSel_g,
                IGalphaSel          = hypoObj$IGalphaSel,
                # Catchability deviations and SDs
                # deltaq_sf           = array(0,dim=c(nS,nF)),
                # lntauq_f            = rep(log(hypoObj$pmtauq_f)),
                deltaqspf_vec       = rep(0, nStockqDevs),
                lntauq_sf           = matrix(log(hypoObj$pmtauq_f),nrow = nS, ncol = nF ,byrow = TRUE),
                mq_f                = hypoObj$mq_f,
                sdq_f               = hypoObj$sdq_f,
                pmtauq_f            = hypoObj$pmtauq_f,
                IGalphaq            = hypoObj$IGalphaq,
                # Steepness
                hBetaPrior          = hypoObj$hBetaPrior,
                lnsigmah_s          = rep(log(hypoObj$pmsigmah),nS),
                lnsigmah            = log(hypoObj$pmsigmah),
                pmsigmah            = hypoObj$pmsigmah,
                IGalphah            = hypoObj$IGalphah,
                # Mortality
                lnsigmaM_s          = rep( log(hypoObj$pmsigmaM), nS ),
                lnsigmaM            = log(hypoObj$pmsigmaM),
                ln_muM              = log(mM),
                sdM                 = sdM,
                pmsigmaM            = hypoObj$pmsigmaM,
                IGalphaM            = hypoObj$IGalphaM,
                # IG Prior on obs error SD
                pmtauObs_g          = hypoObj$pmtauObs_g,
                IGalphaObs          = hypoObj$IGalphaObs,
                # Species effect on steepness
                epsSteep_s          = rep(0,nS),
                # Species effect on M
                epsM_s              = rep(0,nS),
                # Species/stock effect on steepness
                epsSteep_sp         = array(0, dim = c(nS,nP)),
                # Species/stock effect on M
                epsM_sp             = array(0, dim = c(nS,nP)),
                epsM_sx             = array(0, dim = c(nS,nX)),
                # Recruitment resids
                omegaR_vec          = rep( 0, nRecDevs),
                omegaRinit_vec      = rep( 0, nInitDevs ),
                lnsigmaR_sp         = array( log(hypoObj$sigmaR), dim = c(nS,nP)),
                pmsigmaR            = hypoObj$pmsigmaR,
                IGalphaR            = hypoObj$IGalphaR,
                # Correlation in recruitment resids
                logitRgamma_sp      = array( 0, dim = c(nS,nP)),
                # Time-varying selectivity
                epsxSel50_vec       = rep( 0, nSelDevs ),
                epsxSelStep_vec     = rep( 0, nSelDevs ),
                lnsigmaSel          = log( hypoObj$sigmaSelDevs ),
                # Time-varying catchability
                epslnq_vec          = rep( 0, nqDevs ),
                lnsigmaepslnq       = log( hypoObj$sigmaqdevs ),
                # Correlations in recruitment
                recCorr_vec         = rep(0,(nS*nP)*(nS*nP - 1)/2),
                IWmode              = IWmode,
                IWnu                = IWnu,
                ## Single-level priors ##
                # Priors on top level selectivity
                pmlnxSel50_sg       = lnxSel50_sg,
                pmlnxSelStep_sg     = lnxSelStep_sg,
                cvxSel              = hypoObj$cvxSel,
                mF                  = hypoObj$mF,
                sdF                 = hypoObj$sdF )



  # Generate special entries for the 
  # base map list that are sensitive to useFleetsIdx
  taumap_spf <- array(  1 + 1:(nS*nP*nF),
                      dim = c(nS,nP,nF) )
  taumap_spf[calcIndex_spf == 0] <- NA


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
  map <- list( )

  # Adjust phase of B0/Rbar based on recOption
  if( hypoObj$recOption == "BH" )
    phases$lnRbar_sp  <- -1
  if( hypoObj$recOption == "avgR" )
  {
    phases$lnB0_sp        <- -1
    phases$logitSteep     <- -1
    phases$epsSteep_s     <- -1
    phases$epsSteep_sp    <- -1
  }

  # Turn off tv sel deviations if not being used
  if( !hypoObj$tvSel | nSelDevs == 0 )
  {
    phases$epsxSel50_vec    <- -1
    phases$epsxSelStep_vec  <- -1
  }

  # Turn off tvq deviations if not used
  if( !hypoObj$tvq | nqDevs == 0 )
    phases$epslnq_vec       <- -1

  # Turn off tvq deviations if not used
  if( !hypoObj$logitq )
  {
    phases$tq50_vec       <- -1
    phases$tq95_vec       <- -1
    data$logitqFleets     <- rep(0, nF)
  }

  # Turn off stock specific devs if nP == 1
  if( nP == 1)
  {
    phases$epsSteep_sp        <- -1
    phases$epsM_sp            <- -1
    phases$lnsigmah_s         <- -1
    phases$lntauq_sf          <- -1
  }

  # Turn off species specific devs if nS == 1
  if( nS == 1 )
  {
    phases$epsSteep_s       <- -1
    phases$epsM_s           <- -1
    phases$deltaq_sf        <- -1     
    phases$lnxSel50_sf      <- -1
    phases$lnxSelStep_sf    <- -1
    phases$lntauq_f         <- -1

  }
  # Turn off sexual dimorphism if nX == 1
  if( nX == 1 )
  {
    phases$epsM_sx          <- -1
  }
  # Turn off recruitment correlation cholesky factor 
  if( hypoObj$recModel == "uncorr" )
  {
    phases$recCorr_vec      <- -1
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
                          regFPhases = hypoObj$regFPhases,
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

  phaseList$repOpt <- NULL
  

  # Update names on report objects
  outList <- list(  repOpt = renameReportArrays(repOpt,data),
                    sdrepOpt = phaseList$sdrep,
                    optOutput = phaseList$optOutput,
                    fYear = fYear, 
                    lYear = lYear,
                    gearLabs = useFleets,
                    species = useSpecies,
                    stocks = useStocks,
                    map = phaseList$map,
                    data = data,
                    optPars = phaseList$optPars,
                    initPars = pars,
                    phaseList = phaseList,
                    ctlList = obj )



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
                      regFPhases = 3,
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

  regFfleets <- data$regFfleets  
  
  #loop over phases
  if(!is.null(maxPhase))
    maxPhase <- min( maxPhase, max(unlist(phases) ) )
  else maxPhase <- max(unlist(phases))


  # Make a data.frame that will hold the phase info
  fitReport <- matrix(NA, nrow = maxPhase + 1, ncol = 8 )
  colnames(fitReport) <- c( "phase",
                            "objFun",
                            "maxGrad",
                            "nPar",
                            "convCode",
                            "convMsg",
                            "time",
                            "mcmcTime" )
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

    
    if( phase_cur %in% regFPhases )
      data$regFfleets <- regFfleets
    else data$regFfleets <- rep(0,length(data$regFfleets))
  
    #remove the random effects if they are not estimated
    random_use <- random[ !random %in% names(map_use) ]
  
    # initialize the parameters at values in previous phase
    params_use <- parameters
    if( phase_cur > 1 ) 
      params_use <- obj$env$parList( opt$par )


    mapNames <- names(map_use)
    parNames <- names(params_use)

    # Check names in map correspond to par names
    if( any( ! names(map_use) %in% names(params_use) ) )
    {
      badNames <- names(map_use)[ !names(map_use) %in% names(params_use)]
      cat( badNames )
      browser()

    }


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

    repInit <- obj$report()



    if( phase_cur == 1 )
    {


      checkInitNaN    <- lapply( X = repInit, FUN = .checkNaN )
      checkInitFinite <- lapply( X = repInit, FUN = .checkFinite )
      if(any(unlist(checkInitNaN)) | any(unlist(checkInitFinite)))
      {
        browser(beep(expr=cat("NaN or Inf items in initial rep\n")))

        whichNaN <- which(unlist(checkInitNaN))
        whichInf <- which(unlist(checkInitFinite))

      }
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
      phaseReports[[phase_cur]]$pars    <- obj$env$parList(opt$par)
      phaseReports[[phase_cur]]$opt     <- opt
      phaseReports[[phase_cur]]$success <- TRUE
      phaseReports[[phase_cur]]$map     <- map_use
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
    

    # Want to test MCMC performance, so let's MCMC every phase!
    if( intMethod == "phaseMCMC" & class(opt) != "try-error" )
    {
      tBegin      <- proc.time()[3]
      params_use  <- obj$env$parList( opt$par )

      obj <- TMB::MakeADFun(  data = data,
                              parameters = params_use,
                              random = NULL,
                              DLL = DLL_use,
                              map = map_use,
                              silent = silent ) 
      mcmc <- tmbstan(  obj, 
                        init = "last.par.best", 
                        iter = mcChainLength,
                        chains = mcChains )

      phaseReports[[phase_cur]]$mcmc <- mcmc

      mcmcTime <- (proc.time()[3] - tBegin)/60

      fitReport[phase_cur,]$mcmcTime <- mcmcTime

    } # END phaseMCMC

  } # close phase loop

  # integration of the posterior using Laplace Approximation
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

    repInit <- obj$report()

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
    randEffList$par     <- obj$eng$parList( opt$par )

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

  # HMC using tmbstan package
  if( intMethod == "MCMC" & class(opt) != "try-error" )
  {
    tBegin      <- proc.time()[3]
    params_use  <- obj$env$parList( opt$par )

    obj <- TMB::MakeADFun(  data = data,
                            parameters = params_use,
                            random = NULL,
                            DLL = DLL_use,
                            map = map_use,
                            silent = silent ) 


    mcmc <- tmbstan(  obj, 
                      init = "last.par.best", 
                      iter = mcChainLength,
                      chains = mcChains )

    outList$mcmc <- mcmc

  }
  
  # Save sdreport object
  if(outList$success & calcSD )
    outList$sdrep <- TMB::sdreport(obj)

  # Save phase reports
  if( savePhases )
    outList$phaseReports      <- phaseReports

  # Now save report object
  if(outList$success)
  {
    outList$repOpt            <- obj$report()
    outList$optPar            <- obj$env$parList( opt$par )
  } else {
    outList$optPar         <- params_use
    outList$repOpt         <- repInit
  }

  # And the remainder of the details
  outList$objfun            <- obj$fn()
  outList$optOutput         <- opt
  outList$map               <- map_use
  outList$maxGrad           <- max(obj$gr())
  outList$fitReport         <- fitReport
  outList$totTime           <- sum(fitReport$time,na.rm = TRUE)
  
  return( outList )  

} # END TMBphase()
