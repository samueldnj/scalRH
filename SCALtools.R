# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
#
# SCALtools.R
#
# Functions for loading, saving, and modifying hierSCAL fit objects.
#
#
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

# Updating the info file.
.updateInfoFile <- function( fitID = 1, groupFolder =".")
{
  .loadFit(fit = fitID, groupFolder = groupFolder)

  .makeInfoFile(reports)
}

# postMakeBatchInfoFiles
postMakeBatchInfoFiles <- function(groupFolder = ".")
{
  # Get list of fits in the groupFolder
  fitDir  <- here::here("Outputs","fits",groupFolder)
  fitList <- list.dirs(fitDir, recursive = FALSE )
  fitList <- fitList[grepl(x = fitList,pattern = "fit_")]
  nFits   <- length(fitList)

  # Now lapply and make info files
  x <- lapply(X = 1:nFits, FUN = .updateInfoFile, groupFolder = groupFolder )

  return()
}

# readBatchInfo()
# Loads a data.frame of info files for
# a given groupFolder
readBatchInfo <- function(batchDir = here("Outputs","fits") )
{
  batchFitDirs <- list.dirs(batchDir,recursive = FALSE)
  batchFitDirs <- batchFitDirs[grepl(x = batchFitDirs, pattern = "fit_")]
  nFits <- length(batchFitDirs)

  # Read in the info files
  infoFiles <- file.path(batchFitDirs,"infoFile.txt")
  info.df <- lapply(X = infoFiles, FUN = lisread )
  info.df <- as.data.frame(do.call(rbind, info.df))

  info.df
}

# loadFit()
# Loads the nominated fit reports object into memory, 
# so that plot functions can be called
# inputs:   fit=ordinal indicator of sim in project folder
# ouputs:   NULL
# usage:    Prior to plotting simulation outputs
.loadFit <- function( fit = 1, groupFolder = "." )
{
  fitFolder <- here::here("Outputs","fits",groupFolder)

  # List directories in project folder, remove "." from list
  dirList <- list.dirs (path=fitFolder,full.names = FALSE,
                        recursive=FALSE)
  # Restrict to fit_ folders, pick the nominated simulation
  fitList <- dirList[grep(pattern="fit",x=dirList)]
  folder <- fitList[fit]

  # Load the nominated blob
  reportsFileName <- paste(folder,".RData",sep="")
  reportsPath <- file.path(fitFolder,folder,reportsFileName)
  load ( file = reportsPath )

  reports$repOpt <- renameReportArrays(reports$repOpt,reports$data)

  # Update the path object - not sure if this is necessary
  reports$path <- file.path(fitFolder, folder) 

  # Assign to global environment
  assign( "reports",reports,pos=1 )

  message("(.loadFit) Reports in ", folder, " loaded\n", sep="" )

  return( file.path(fitFolder, folder) )
} # END .loadFit()

# Make a vector of complex prefixes - need to
# pick off first letter of each name
makeCplxPrefix <- function( cplxVec )
{
  initials <- substr(cplxVec,1,1)

  prefix <- paste(initials, collapse = "")
}

# makeInfoFile()
# Creates an info list for the
# fit object in the argument
.makeInfoFile <- function( obj )
{
  infoList <- list()

  # Now start populating
  infoList$dataScenario <- obj$ctlList$ctrl$dataScenarioName
  infoList$modelHyp     <- obj$ctlList$ctrl$modelHypName
  infoList$path         <- obj$path
  infoList$complex      <- makeCplxPrefix(obj$ctlList$data$species)


  outFile <- file.path(obj$path, "infoFile.txt")


  # Create an output file
  cat('## Info file for hierSCAL model fit\n', file = outFile, append = FALSE, sep = "")
  cat('## Written ', Sys.time(), "\n", file = outFile, append = TRUE, sep = "")
  cat("", file = outFile, append = TRUE, sep = "")

  # Now start writing info list
  for( lIdx in 1:length(infoList) )
  {
    cat("# ", names(infoList)[lIdx], "\n", file = outFile, append = TRUE, sep = "")
    cat( infoList[[lIdx]], "\n", file = outFile, append = TRUE, sep = "")
    cat("", file = outFile, append = TRUE, sep = "")
  }

  cat("## End File\n", file = outFile, append = TRUE, sep = "")

  message("Info file created at ", outFile, ".\n", sep = "")
} # END .makeInfoFile()


# rerunPlots()
# Function to redo all plots from a given
# fit report object
rerunPlots <- function( fitID = 1 )
{
  # Load the fit object
  .loadFit(fitID)

  if(!is.null(reports$repOpt))
    reports$repOpt <- renameReportArrays( reports$repOpt, reports$data )

  if(is.null(reports$refPoints))
    reports$repOpt <- calcRefPts(reports$repOpt)

  # rerun savePlots
  savePlots(  fitObj = reports,
              saveDir = reports$path )

  cat("MSG (rerunPlots) Plots redone in ", reports$path, "\n", sep = "")
} # END rerunPlots()

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
  fltGrpNames <- c("commTrawl","HSAss","Synoptic")

  nS <- length(specNames)
  nP <- length(stockNames)

  # Create species/stock names vector
  specStock <- character(length = nS*nP)
  for( s in 1:nS )
    for( p in 1:nP )
      specStock[(s-1) * nP + p] <- paste( specNames[s],"_", stockNames[p],sep = "")




  # Ok, that's the data taken care of. There are still all the
  # new arrays that we created
  # Predicted data
  dimnames(repObj$I_spft_hat)       <- dimnames(datObj$I_spft)
  dimnames(repObj$aDist_aspftx_hat) <- dimnames(datObj$age_aspftx)
  dimnames(repObj$lDist_lspftx_hat) <- dimnames(datObj$len_lspftx)
  # State arrays
  dimnames(repObj$B_asptx)          <- dimnames(datObj$age_aspftx)[c(1:3,5)]
  dimnames(repObj$N_axspt)          <- dimnames(datObj$age_aspftx)[c(1,6,2:3,5)]
  dimnames(repObj$B_spt)            <- dimnames(datObj$age_aspftx)[c(2:3,5)]
  dimnames(repObj$R_spt)            <- dimnames(datObj$age_aspftx)[c(2:3,5)]
  dimnames(repObj$SB_spt)           <- dimnames(datObj$age_aspftx)[c(2:3,5)]
  dimnames(repObj$vB_spft)          <- list(  species = specNames,
                                              stock = stockNames,
                                              fleet = gearNames,
                                              year = yearNames )
  dimnames(repObj$predC_spft)       <- dimnames(datObj$age_aspftx)[c(2:5)]
  dimnames(repObj$predCw_spft)      <- dimnames(datObj$age_aspftx)[c(2:5)]
  dimnames(repObj$C_axspft)         <- dimnames(datObj$age_aspftx)[c(1,6,2:5)]
  dimnames(repObj$Cw_axspft)        <- dimnames(datObj$age_aspftx)[c(1,6,2:5)]
  dimnames(repObj$F_spft)           <- dimnames(datObj$age_aspftx)[c(2:5)]
  dimnames(repObj$Z_aspxt)          <- dimnames(datObj$age_aspftx)[c(1:3,6,5)]
  # Biological parameters
  dimnames(repObj$R0_sp)          <- dimnames(datObj$age_aspftx)[c(2:3)]
  dimnames(repObj$B0_sp)          <- dimnames(datObj$age_aspftx)[c(2:3)]
  dimnames(repObj$h_sp)           <- dimnames(datObj$age_aspftx)[c(2:3)]
  dimnames(repObj$phi_sp)         <- dimnames(datObj$age_aspftx)[c(2:3)]
  dimnames(repObj$reca_sp)        <- dimnames(datObj$age_aspftx)[c(2:3)]
  dimnames(repObj$recb_sp)        <- dimnames(datObj$age_aspftx)[c(2:3)]
  dimnames(repObj$sigmaR_sp)      <- dimnames(datObj$age_aspftx)[c(2:3)]
  dimnames(repObj$omegaRmat_spt)  <- list(  specStock = specStock, 
                                            year = yearNames )
  dimnames(repObj$recCorrMat_sp)  <- list(  specStock = specStock,
                                            specStock = specStock )

  # Observation models
  dimnames(repObj$q_spf)        <- dimnames(datObj$age_aspftx)[c(2:4)]  
  dimnames(repObj$q_spft)       <- dimnames(datObj$age_aspftx)[c(2:5)]  
  dimnames(repObj$tau2Obs_spg)  <- list(  species = specNames,
                                          stock = stockNames,
                                          group = fltGrpNames )  
  dimnames(repObj$tauObs_spg)   <- list(  species = specNames,
                                          stock = stockNames,
                                          group = fltGrpNames )
  dimnames(repObj$tau2Obs_spf)  <- list(  species = specNames,
                                          stock = stockNames,
                                          fleet = gearNames )  
  dimnames(repObj$tauObs_spf)   <- list(  species = specNames,
                                          stock = stockNames,
                                          fleet = gearNames )
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

  dimnames(repObj$L1_spx)   <- list(  species = specNames,
                                      stock = stockNames,
                                      sex = sexNames )
  dimnames(repObj$L2_spx)   <- list(  species = specNames,
                                      stock = stockNames,
                                      sex = sexNames )
  dimnames(repObj$vonK_spx) <- list(  species = specNames,
                                      stock = stockNames,
                                      sex = sexNames )
  dimnames(repObj$M_spx)    <- list(  species = specNames,
                                      stock = stockNames,
                                      sex = sexNames )




  return(repObj)
} # END renameReportArrays()


# makeFitReport()
# Makes an html fit report from the report
# object, based on a fit report template
makeFitReport <- function( fitID = 1, groupFolder = "." )
{
  fitFolder <- here::here("Outputs","fits",groupFolder)

  # List directories in project folder, remove "." from list
  dirList <- list.dirs (path=fitFolder,full.names = FALSE,
                        recursive=FALSE)
  # Restrict to fit_ folders, pick the nominated simulation
  fitList <- dirList[grep(pattern="fit",x=dirList)]

  if( is.character(fitID) )
    folder <- fitList[grepl(x = fitList, pattern = fitID) ]
  else
    folder <- fitList[fitID]

  # Load the nominated blob
  reportsFileName <- paste(folder,".RData",sep="")
  reportsPath <- file.path(fitFolder,folder,reportsFileName)
  fitFolderPath <- here::here("Outputs","fits",groupFolder,folder)

  # Create parameter list for rendering the document
  params <- list( rootDir= fitFolderPath,
                  RdataFile = reportsFileName)
  # Make an output file name
  outFile <- paste( "fitReport.html", sep = "")


  # Render
  rmarkdown::render(  input = here::here("Documentation","fitReportTemplate.Rmd"), 
                      output_file = outFile,
                      output_dir = fitFolderPath,
                      params = params,
                      envir = new.env(),
                      output_format = "bookdown::html_document2" )

  # remove temporary files
  dataReportFiles <- "fitReport_files"
  unlink(file.path(fitFolderPath,dataReportFiles), recursive = TRUE)
} # END makeFitReport()

# savePlots()
savePlots <- function(  fitObj = reports,
                        saveDir = "./Outputs/fits/" )
{
  # Create a plots sub directory in the saveDir
  saveDir <- file.path(saveDir,"plots")
  if(!dir.exists(saveDir))
    dir.create(saveDir)

  fYear <- fitObj$fYear
  lYear <- fitObj$lYear

  report   <- c(fitObj$repOpt,fitObj$data)

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
