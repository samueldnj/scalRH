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
  nL_s      <- dataObj$nL_s
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
                            nP = nP, years = fYear:lYear,
                            collapseComm = FALSE )[useSpecies,useStocks,useFleets,, drop = FALSE]

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
  # Sum catch
  sumCat_sp <- apply( X = C_spft, FUN = sum, MARGIN = c(1,2) )

  # Load growth data - use this
  # to inform the vonB model in the
  # assessment
  growthFiles <- dataObj$growthFiles
  growthFilePaths <- file.path("./Data",growthFiles )
  for( i in 1:length(growthFiles))
    load(growthFilePaths[i])
  # # Make the age-length key - this might be useful,
  # # or we can include it so we can integrate a growth
  # # model
  # ALK_spalx <- makeALFreq(  lenAgeList = lenAge,
  #                           combineSex = TRUE )[useSpecies, useStocks,,,drop = FALSE]


  # Load age and length compositions
  load(file.path("./Data",dataObj$ageData))
  load(file.path("./Data",dataObj$lenData))

  # Create compositional arrays
  age_aspft <- makeCompsArray(  compList = ageComps,
                                plusGroups = nA_s,
                                collapseComm = FALSE,
                                fleetIDs = useFleets,
                                combineSex = TRUE,
                                years = fYear:lYear,
                                xName = "ages" )
  age_aspft <- age_aspft[,useSpecies,useStocks,useFleets,, drop = FALSE]


  len_lspft <- makeCompsArray(  compList = lenComps,
                                plusGroups = nL_s,
                                collapseComm = FALSE,
                                fleetIDs = useFleets,
                                combineSex = TRUE,
                                years = fYear:lYear,
                                xName = "length" )
  len_lspft <- len_lspft[,useSpecies,useStocks,useFleets,, drop = FALSE]


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

  # Use useSpec, useStocks, useFleets and yrChar to
  # reduce the data arrays so that the fits are dynamic

  # Pull prior mean values from hypoObj to keep pars list tidier
  mh    <- hypoObj$hPriorMean
  sdh   <- hypoObj$hPriorSD

  # Observation errors
  tau2ObsIGa  <- hypoObj$tau2ObsIGa[useFleetsIdx]
  tau2ObsMode <- hypoObj$tau2ObsPriorMode[useFleetsIdx]
  tau2ObsIGb  <- (tau2ObsIGa + 1) * tau2ObsMode

  # Generate the data list
  data <- list( I_spft = I_spft,
                C_spft = C_spft,
                D_spft = D_spft,
                # ALK_spal = ALK_spalx,
                age_aspft = age_aspft,
                len_lspft = len_lspft,
                type_f = as.integer(rep(0,nF)),
                A_s = as.integer(nA_s),
                L_s = as.integer(nL_s),
                lenD_s = rep(0,nS),
                swRinit_s = as.integer(initFished_s[useSpecies]),
                parSwitch = 0,
                calcIndex_spf = calcIndex_spf,
                ageLikeWt = 1,
                lenLikeWt = 1,
                tFirstRecDev_s = as.integer(tFirstRecDev_s),
                tLastRecDev_s = as.integer(tLastRecDev_s) )

  

  # Generate parameter list
  pars <- list( lnB0_sp = log(sumCat_sp),
                logitSteep_sp = array(log((mh - .2)/(1 - mh)),dim =c(nS,nP)) ,
                lnM_sp = array(log(.04),dim =c(nS,nP)),
                lnLinf_sp = log(vonBFits$Linf_sp[useSpecIdx,useStockIdx,drop=FALSE]),
                lnvonK_sp = log(vonBFits$VonK_sp[useSpecIdx,useStockIdx,drop=FALSE]),
                lnL1_sp   = log(matrix(vonBFits$L1_s[useSpecIdx],nrow = nS, ncol = nP, byrow = FALSE)),
                LWa_s = LWa_s,
                LWb_s = LWb_s,
                xMat50_s = xMat50,
                xMat95_s = xMat95,
                lnq_spf = array(0,dim =c(nS,nP,nF)),
                lntau_spf = array(log(0.1),dim =c(nS,nP,nF)),
                lnlenSel50_sf = array(3.5,dim =c(nS,nF)),
                lnlenSelStep_sf = array(1.5,dim =c(nS,nF)),
                lnF_spft  = rep(-2,length = nPosCatch),
                lntauC_f  = rep(log(0.1),nF),
                lntauD_f  = rep(log(0.1),nF),
                muLinf_s  = vonBFits$Linf_s[useSpecIdx],
                sigmaLinf_s = rep(.1,nS),
                muvonK_s  = rep(.2,nS),
                sigmavonK_s = rep(.2,nS),
                muL1_s  = rep(25,nS),
                sigmaL1_s = rep(.2,nS),
                lnsigmaL_s = rep(-2.3,nS),
                mulenSel50_f  = rep(3.4,nF),
                mulenSel95_f  = rep(4.2,nF),
                sigmalenSel50_f = rep(.2,nF),
                sigmalenSel95_f = rep(.2,nF),
                lnqbar_fs = array(log(0.5),dim =c(nF,nS)),
                lntauq_fs = array(log(.5),dim =c(nF,nS)),
                lnqbar_f  = rep(log(0.5),nF),
                lntauq_f  = rep(log(0.5),nF),
                logitSteep_s  = rep(log((mh - .2)/(1 - mh)),nS),
                lnsigmah_s  = rep(log(sdh),nS),
                logitSteep  = log((mh - .2)/(1 - mh)),
                lnsigmah  = log(sdh),
                lnM_s = rep(log(0.04),nS),
                lnsigmaM_s  = rep( -2, nS ),
                ln_muM  = log(0.04),
                lnsigmaM  = -2,
                IGatau_f = tau2ObsIGa,
                IGbtau_f = tau2ObsIGb,
                omegaR_vec = rep(0, nRecDevs),
                omegaRinit_vec = rep( 0, nInitDevs ),
                lnsigmaR_sp = array(0, dim = c(nS,nP)),
                logitRCorr_chol = rep(0, nS * nP),
                logitRgamma_sp  = array(0, dim = c(nS,nP)),
                deltaLenSel50_p = rep(0,nP),
                deltaLenSelStep_p = rep(0,nP),
                lnsigmaSel = log(0.05),
                pmLenSel95_sf = array(45,dim = c(nS,nF)),
                cvLenSel95 = 0.3  )

  # Generate map list - this will have to be updated
  # so that it's sensitive to useFleetsIdx
  qmap_spf <- array(  sample(1:100,1) + 1:(nS*nP*nF),
                      dim = c(nS,nP,nF) )
  for( fIdx in 1:nF)
    qmap_spf[,,fIdx] <- fIdx + sample(1:100,1)

  qmap_spf[calcIndex_spf == 0] <- NA
  # map out the observation error SD
  taumap_spf <- array(  sample(1:100,1) + 1:(nS*nP*nF),
                        dim = c(nS,nP,nF) )

  for( fIdx in 1:nF)
    taumap_spf[,,fIdx] <- fIdx + sample(101:200,1)

  taumap_spf[calcIndex_spf == 0] <- NA

  # 

  map <- list(  # lnB0_sp = factor(array(NA,dim =c(nS,nP))),
                # logitSteep_sp = factor(array(NA,dim =c(nS,nP)) ),
                lnM_sp = factor(array(NA,dim =c(nS,nP))),
                lnLinf_sp = factor(array(NA,dim =c(nS,nP))),
                lnvonK_sp = factor(array(NA,dim =c(nS,nP))),
                lnL1_sp = factor(array(NA,dim =c(nS,nP))),
                LWa_s = factor(rep(NA,nS)),
                LWb_s = factor(rep(NA,nS)),
                xMat50_s = factor(rep(NA,nS)),
                xMat95_s = factor(rep(NA,nS)),
                lnq_spf = factor(qmap_spf),
                lntau_spf = factor(taumap_spf),
                # lnlenSel50_sf = factor(array(NA,dim =c(nS,nF))),
                # lnlenSelStep_sf = factor(array(NA,dim =c(nS,nF))),
                deltaLenSel50_p = factor(rep(NA,nP)),
                deltaLenSelStep_p = factor(rep(NA,nP)),
                # lnF_spft  = factor(rep(NA,nPosCatch)),
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
                # logitSteep_s  = factor(rep(NA,nS)),
                lnsigmah_s  = factor(rep(-NA,nS)),
                logitSteep  = factor(NA),
                lnsigmah  = factor(NA),
                lnM_s = factor(rep(NA,nS)),
                lnsigmaM_s  = factor(rep(NA, nS ) ),
                ln_muM  = factor(NA),
                lnsigmaM  = factor(NA),
                IGatau_f = factor(rep(NA,nF)),
                IGbtau_f = factor(rep(NA,nF)),
                # omegaR_vec = factor( rep( NA, nRecDevs ) ),
                # omegaRinit_vec = factor( rep( NA, nInitDevs ) ),
                lnsigmaR_sp = factor(array(NA, dim = c(nS,nP)) ),
                logitRCorr_chol = factor(rep(NA, nS * nP)),
                logitRgamma_sp  = factor(array(NA, dim = c(nS,nP))),
                lnsigmaSel = factor(NA),
                pmLenSel95_sf = factor(array(NA,dim = c(nS,nF))),
                cvLenSel95 = factor(NA) ) 


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
                    stocks = useStocks )

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
      outList$repFE         <- objFE$report()
      outList$sdrepFE       <- summary(sdreport(objFE))
      outList$heFE          <- objFE$he()
      outList$fitFE         <- fitFE
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

  # Ok, that's the data taken care of. There are still all the
  # new arrays that we created
  dimnames()

  return(report)
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

  # Plot catch fits
  png(  file.path(saveDir,"plotSRsp.png"),
        width = 11, height = 8.5, units = "in", res = 300)
  plotSRsp( repObj = report, initYear = fYear )
  dev.off()

  for( sIdx in 1:nS )
  {
    # Make a directory to hold species specific plots
    specDir <- paste("species",sIdx,sep = "")
    specPath <- file.path(saveDir,specDir)
    if(!dir.exists(specPath))
      dir.create(specPath)

    for( pIdx in 1:nP)
    {
      # Make a directory to hold stock specific plots
      stockDir <- paste("stock",pIdx,sep = "")
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



