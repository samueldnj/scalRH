# -----------------------------------------------------------------------------
#
# DERPAfuns.R
#
# Functions for DERPAdata.R and ageStructuredControl.R
# 
# Last revised: Nov 4, 2018
# 
# -----------------------------------------------------------------------------

# OK, we need functions to summarise the biological data 
# for a given species, given a determined stock structure

loadStockSpecNameLists <- function()
{
  # Survey ids for plotting/legends
  surveyIDs <<-  c( QCSSyn = 1, 
                    HSAss = 2, 
                    HSSyn = 3, 
                    WCVISyn = 4,
                    WCHGSyn = 16 )

  # Species codes
  specCodes <<- list( "Dover" = 626,
                      "English" = 628,
                      "Rock" = 621,
                      "Petrale" = 607,
                      "Arrowtooth" = 602 )

  # Stock IDs for grouping data
  stocksSurvey <<- list(  HSHG = c(2,3,16),
                          QCS = c(1),
                          WCVI = c(4) )
  stocksCommBio <<- list( HSHG = c(7,8,9),
                          QCS = c(5,6),
                          WCVI = c(3,4) )
  stocksCommCPUE  <<- list( HSHG = "5CDE",
                            QCS = "5AB",
                            WCVI = "3CD" )

  # Species names for reading data
  survSpecNames <<- c(  Dover = "dover",
                        English = "english",
                        Rock = "srock",
                        Petrale = "petrale",
                        Arrowtooth = "atooth" )

  commSpecNames <<- c(  Dover = "dover-sole",
                        English = "english-sole",
                        Rock = "southern-rock-sole",
                        Petrale = "petrale-sole",
                        Arrowtooth = "arrowtooth-flounder" )

  commFleetYrRange <<- list(  comm.hist = 1954:1995,
                              comm.mod = 1996:2018 )

  invisible(NULL)  
}


# Load CRS codes
loadCRS <- function()
{
  AEAproj <<- CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m")
  LLproj  <<- CRS("+proj=longlat +datum=WGS84")
  UTMproj <<- CRS("+proj=utm +zone=9 +datum=WGS84")


  invisible(NULL)
}


# openShapeFile()
openShapeFile <- function(  layer = "HS_Synoptic_Survey_Active_Blocks",
                            path = "./Data/ShapeFiles/SynSurveyBlocks/",
                            outCRS = UTMproj, inCRS = NULL )
{

  # Read in grid
  grid <- readOGR(dsn = path, layer = layer)

  # Convert to outCRS
  if(!is.null(inCRS))
    proj4string(grid) <- inCRS

  # Project to outCRS
  if(!is.null(outCRS))
  {
    grid <- spTransform( x = grid, CRSobj = outCRS )
  }

  grid
}

# makeIndexArray()
# Takes groomed data for each species and combines them into
# a data array for feeding to hierSCAL.
makeIndexArray <- function( relBio = relBioList_Survey,
                            commCPUE = commCPUEList,
                            nP = 3, years = 1954:2018,
                            collapseComm = FALSE,
                            scaleComm = 1e3
                          )
{
  nS <- length(relBio)
  nT <- length(years)

  nSurv   <- dim(relBio$Dover$relBio.arr)[1]
  survIDs <- dimnames(relBio$Dover$relBio.arr)[[1]]

  stockIDs <- dimnames(relBio$Dover$relBio.arr)[[2]]

  nComm   <- dim(commCPUE$Dover$cpue.arr)[1]
  commIDs <- dimnames(commCPUE$Dover$cpue.arr)[[1]]

  nF <- nSurv + nComm

  if( collapseComm )
    nF <- nF - 1


  I_spft <- array( NA,  dim = c( nS, nP, nF, nT ),
                        dimnames = list(  species = names(relBio),
                                          stocks = stockIDs,
                                          fleets = c(commIDs,survIDs),
                                          years = years ) )

  # Now loop over species and fill
  for( specIdx in 1:nS )
  {
    specID <- names(relBio)[specIdx]

    subRelBio   <- relBio[[specID]]$relBio.arr[,,as.character(years),1]
    subCommCPUE <- commCPUE[[specID]]$cpue.arr[,,as.character(years),1,1]

    subCommCPUE[subCommCPUE < 0] <- NA
    subCommCPUE <- exp(subCommCPUE)

    if( collapseComm )
    {
      # Add the two fleets together (should be zero overlap)
      subCommCPUE <- apply( X = subCommCPUE, FUN = sum, na.rm = T, MARGIN = c(2,3) )
      subCommCPUE[is.na(subCommCPUE)] <- -1
      I_spft[specID,,1,as.character(years)] <- subCommCPUE/scaleComm
    } else {
      subCommCPUE[is.na(subCommCPUE)] <- -1
      I_spft[specID,,1,as.character(years)] <- subCommCPUE["comm.hist",,as.character(years)]/scaleComm
      I_spft[specID,,2,as.character(years)] <- subCommCPUE["comm.mod",,as.character(years)]/scaleComm
    }

    for( survIdx in 1:nSurv )
    {
      survID <- survIDs[survIdx]
      I_spft[ specID,,survID,as.character(years) ] <- subRelBio[survID,,as.character(years)]
    }
  }

  return(I_spft)
} # END makeIndexArray()

# makeSurveyCatchStocks()
# Sum the catch in the surveys, arrange by year
# etc
makeSurveyCatchStocks <- function(  spec = "dover",
                                    years = c(fYear, lYear), 
                                    stocks = stocksSurvey,
                                    survIDs = surveyIDs )
{
  # Read in density
  specDensityFile <- paste(spec,"density.csv",sep = "_")
  specDensityPath <- file.path(getwd(),"Data","density",specDensityFile)
  densityTab <- read.csv(specDensityPath, header = T)

  # Now join and select the columns we want
  surveyCatch <-  densityTab %>%
                  dplyr::select(  year = YEAR,
                                  tripID = TRIP_ID,
                                  eventID = FISHING_EVENT_ID,
                                  majorArea = MAJOR_STAT_AREA_CODE,
                                  minorArea = MINOR_STAT_AREA_CODE,
                                  survey = SURVEY_DESC,
                                  surveyID = SURVEY_ID,
                                  survSeriesID = SURVEY_SERIES_ID,
                                  catch = CATCH_WEIGHT ) %>%
                  filter( survSeriesID %in% survIDs ) %>%
                  mutate( catch = catch/1e6,
                          stockName = sapply( X = survSeriesID, 
                                              FUN = appendName, 
                                              stocks ),
                          surveyName = sapply(  X = survSeriesID,
                                                FUN = appendName,
                                                survIDs ) ) %>%
                  group_by( stockName, surveyName, year ) %>%
                  summarise( catch = sum(catch) ) %>%
                  ungroup()


  surveyCatch
} # makeSurveyCatchStocks()


filterSurveyBlocks <- function( blocks = grids$HS,
                                density = synTab,
                                plot = TRUE,
                                species = "dover",
                                survey = "HS",
                                stratAreas = stratArea )
{
  # First, get survey series ID from blocks
  survSeriesID <- unique(blocks$SURVEY_SERI)

  # Filter density to that survey
  density <- density %>%
              filter(SURVEY_SERIES_ID == survSeriesID )
  

  # Now we want to 
  # 1. match survey sets to blocks,
  # Trying to use the over function, we get a pretty good 
  # coverage for the data set, but there are some missing points
  
  # This needs to be done by stratum (grouping code)
  # so that tows outside blocks are assigned to the correct
  # stratum
  grpCodes    <- unique(density$GROUPING_CODE)
  nGrpCodes   <- length(grpCodes)

  stratAreas <- stratAreas %>% 
                rename( grCode = GROUPING_CODE ) %>%
                filter( grCode %in% grpCodes )

  densByGrpCode   <-  vector(mode = "list", length = nGrpCodes )
  blocksByGrpCode <-  vector(mode = "list", length = nGrpCodes )
  summGrpCode     <-  vector(mode = "list", length = nGrpCodes )

  for( cIdx in 1:nGrpCodes )
  {
    # Subset to sets in the correct stratum
    grpCode <- grpCodes[cIdx]
    subDensity <- density %>%
                  filter( GROUPING_CODE == grpCode )

    subBlocks <- blocks[blocks@data$"GROUPING_CO" == grpCode, ]

    # Convert LL density data to UTM so it's the same as the grids
    setCoords.LL          <- subDensity[,c("LONGITUDE","LATITUDE")]
    setCoords.sp          <- SpatialPoints(coords = setCoords.LL, proj4string = LLproj )
    setCoords.UTM         <- spTransform( x = setCoords.sp, CRSobj = UTMproj )
    setCoords.UTM.df      <- as.data.frame(setCoords.UTM)

    subDensity$x          <- setCoords.UTM.df[,1]
    subDensity$y          <- setCoords.UTM.df[,2]
    subDensity$block      <- NA
    subDensity$assGrCd    <- NA
    subDensity$blockDepth <- NA

    tryPoints             <- over( x = setCoords.UTM, y = subBlocks )
    naSets                <- which(is.na(tryPoints[,1]))

    subDensity$block      <- tryPoints[,"BLOCK_DESIG"]
    subDensity$assGrCd    <- tryPoints[,"GROUPING_CO"]
    subDensity$blockDepth <- tryPoints[,"DEPTH_M"]

    if(length(naSets) > 0)
    {
      # Now find nearest blocks to unassigned points
      naPoints  <- setCoords.UTM[naSets]
      # Compute distances
      distMtx   <- gDistance( naPoints, subBlocks, byid = TRUE )
      # Find which.min
      minBlock  <- apply(X = distMtx, FUN = which.min, MARGIN = 2)

      subDensity[naSets,"block"]      <- subBlocks[minBlock,][["BLOCK_DESIG"]]
      subDensity[naSets,"assGrCd"]    <- subBlocks[minBlock,][["GROUPING_CO"]]
      subDensity[naSets,"blockDepth"] <- subBlocks[minBlock,][["DEPTH_M"]]
    }

    densByGrpCode[[cIdx]]     <- subDensity
    blocksByGrpCode[[cIdx]]   <- subBlocks


  }

  densityBlocks <- do.call(rbind, densByGrpCode) %>%
                    mutate(densityKGPM2 = CATCH_WEIGHT / areaFished_m2 )

  # 2. identify blocks with positive tows >= minPosTow
  posTowBlocks    <-  densityBlocks %>%
                      group_by( block, assGrCd ) %>%
                      summarise(  nTows = n(),
                                  catchWt = sum(CATCH_WEIGHT) ) %>%
                      filter( catchWt > 0 )
  zeroTowBlocks   <-  densityBlocks %>%
                      group_by( block, assGrCd ) %>%
                      summarise(  nTows = n(),
                                  catchWt = sum(CATCH_WEIGHT) ) %>%
                      filter( catchWt == 0 )

  posTowBlockIDs    <- posTowBlocks$block
  zeroTowBlockIDs   <- zeroTowBlocks$block

  grpSumm <- matrix( NA, ncol = 10, nrow = 1)  

  colnames(grpSumm) <- c( "grCode", 
                          "nBlocks", 
                          "pBlocks", 
                          "zBlocks",
                          "minDepth",
                          "maxDepth",
                          "minDepthP",
                          "maxDepthP",
                          "minDepthZ",
                          "maxDepthZ")

  grpSumm <- as.data.frame(grpSumm)

  if(plot)
  {
    # Load land masses
    data(nepacLL)
    nePac <- PolySet2SpatialPolygons(nepacLL)
    nePacUTM <- spTransform(x = nePac, CRSobj = UTMproj)

    nCols <- ceiling(sqrt(nGrpCodes))
    nRows <- ceiling(nGrpCodes/nCols)

    mapExtent <- extent(blocks)

    savePlotRoot <- paste(species,survey,sep = "")
    plotFile <- paste( savePlotRoot,"blockDesign.png", sep = "" )
    plotDir <- file.path("./Outputs/speciesSurveyData",species)

    if(!dir.exists(plotDir))
      dir.create(plotDir)

    plotPath <- file.path(plotDir,plotFile)

    png(plotPath, width = 11, height = 11, units = "in",
          res = 400 )

    par(  mfrow = c(nRows,nCols), 
          mar = c(0.5,0.5,0.5,0.5), 
          oma = c(3,3,3,3) )
    # Loop and plot each grouping code
    for( cIdx in 1:nGrpCodes )
    {
      subBlocks   <- blocksByGrpCode[[cIdx]]
      subDensity  <- densByGrpCode[[cIdx]]

      grpSumm$grCode <- grpCodes[cIdx]

      # Get some info about the stratum
      grpSumm$nBlocks     <- length(subBlocks[[1]])
      # Subset to positive and zero observations 
      posSubBlocks  <- subBlocks[subBlocks@data$BLOCK_DESIG %in% posTowBlockIDs, ]
      zeroSubBlocks <- subBlocks[subBlocks@data$BLOCK_DESIG %in% zeroTowBlockIDs, ]

      # Count number of blocks with positive and zero obs
      grpSumm$pBlocks    <- length(posSubBlocks[[1]])
      grpSumm$zBlocks    <- length(zeroSubBlocks[[1]])

      grpSumm[,c("minDepth","maxDepth")]    <- range(subBlocks$DEPTH_M)
      grpSumm[,c("minDepthP","maxDepthP")]  <- range(posSubBlocks$DEPTH_M)
      grpSumm[,c("minDepthZ","maxDepthZ")]  <- range(zeroSubBlocks$DEPTH_M)

      
      title     <- paste("GC", grpCodes[cIdx], sep = "" )


      # Plot maps of blocks
      plot( mapExtent, type = "n", xlab = "", ylab = "",
            axes = FALSE )
      plot( subBlocks, add = TRUE )
      box()
      mtext( side = 3, text = title, font = 2 )

      # Plot presences and absences
      plot(posSubBlocks, add =TRUE, col = "blue")
      plot(zeroSubBlocks, add =TRUE, col = "red")

      # Plot landmasses
      plot(nePacUTM, add = TRUE, col = "grey40", border = NA )

      if(cIdx == 1)
        legend( x = "topright",
                fill = c("blue","red"),
                legend = c("Positive tow observed","No positive tows") )   

      summGrpCode[[cIdx]] <- grpSumm

    }

    dev.off()
    
  }

  groupSummary <- do.call(rbind,summGrpCode)

  summDensSplitPos <- densityBlocks %>%
                      rename( grCode = GROUPING_CODE) %>%
                      group_by( grCode, CATCH_WEIGHT > 0 ) %>%
                      summarise(  meanWt = mean(CATCH_WEIGHT),
                                  minWt = min(CATCH_WEIGHT),
                                  maxWt = max(CATCH_WEIGHT),
                                  meanDens = mean(densityKGPM2),
                                  nTows = n() ) %>%
                      ungroup()

  summDens <- densityBlocks %>%
              rename( grCode = GROUPING_CODE) %>%
              group_by( grCode ) %>%
              summarise(  meanWt = mean(CATCH_WEIGHT),
                          maxWt = max(CATCH_WEIGHT),
                          meanDens = mean(densityKGPM2) ) %>%
              ungroup()

  groupSummary <- groupSummary %>%
                  left_join( summDens ) %>%
                  left_join( stratAreas )

  tabFile <- paste(species,survey,"StratSummary.csv",sep = "")
  tabPath <- file.path(plotDir,tabFile)
  write.csv( groupSummary, file = tabPath )

  # Now plot the depth profile of the kg/m2 density
  plotFile <- paste( savePlotRoot,"densByDepth.png", sep = "" )
  plotPath <- file.path(plotDir,plotFile)
  png(  plotPath, width = 11, height = 6, units = "in",
        res = 400 )  
  depRan  <- -1 * range(densityBlocks$blockDepth)
  densRan <- range(densityBlocks$densityKGPM2)
  wtRan   <- range(densityBlocks$CATCH_WEIGHT)
  cols    <- brewer.pal(n = nGrpCodes, "Dark2" )
  cols    <- alpha( cols, .4 )
  names(cols) <- grpCodes
  par( mfrow = c(2,1), mar =c(1,2,1,0), oma = c(3,3,3,3) )
  plot( x = depRan, y = densRan, type = "n",
        xlab = "Depth (m)", ylab = "Density (kg/m2)",
        las = 1 )
    points( x = -1 * densityBlocks$blockDepth,
            y = densityBlocks$densityKGPM2,
            col = cols[as.character(densityBlocks$assGrCd)],
            pch = 4, cex = .8 )
    abline( v = stratAreas$MIN_DEPTH, lty = 2, lwd = .8, col = "grey60" )
    abline( v = stratAreas$MAX_DEPTH, lty = 2, lwd = .8, col = "grey60" )

  plot( x = depRan, y = wtRan, type = "n",
        xlab = "Depth (m)", ylab = "Catch weight (kg)",
        las = 1 )
    points( x = -1 * densityBlocks$blockDepth,
            y = densityBlocks$CATCH_WEIGHT,
            col = cols[as.character(densityBlocks$assGrCd)],
            pch = 4, cex = .8 )
    abline( v = stratAreas$MIN_DEPTH, lty = 2, lwd = .8, col = "grey60" )
    abline( v = stratAreas$MAX_DEPTH, lty = 2, lwd = .8, col = "grey60" )

  mtext( side = 1, text = "Depth (m)", outer = T, line  = 2 )

  dev.off()


  # 3. Calculate strata area from new subset of blocks
  # 3. filter out sets from < minPosTow blocks



}


# Calculate relative biomass by species and arrange in an array
# for feeding to TMB model
makeRelBioStocks <- function( spec = "dover",
                              years = c(1975, 2016), 
                              stocks = list(  HGHS = c(2,3,16),
                                              QCS = c(1),
                                              WCVI = c(4) ),
                              survIDs = surveyIDs,
                              stratArea = stratData,
                              grids = grids   )
{

  # Read in density
  specDensityFile <- paste(spec,"density.csv",sep = "_")
  specDensityPath <- file.path(getwd(),"Data","density",specDensityFile)
  densityTab <- read.csv(specDensityPath, header = T)

  # first, calculate the tow length from speed and distance, if
  # it doesn't exist
  calcLengths <-  densityTab %>% 
                  filter( is.na(TOW_LENGTH_M) ) %>%
                  mutate( TOW_LENGTH_M = DURATION_MIN * SPEED_MPM )

  # Now rejoin those back to original data frame
  densityTab <- densityTab %>% 
                filter( !is.na( TOW_LENGTH_M ) ) %>%
                rbind( calcLengths )

  # Now fill in missing door spreads
  aveSpread <- mean(densityTab$DOORSPREAD_M, na.rm = T)
  missingSpreads <- densityTab %>% filter( is.na( DOORSPREAD_M ) ) %>%
                    mutate( DOORSPREAD_M = aveSpread )

  densityTab <- densityTab %>%
                filter( !is.na(DOORSPREAD_M) ) %>%
                rbind( missingSpreads )

  # Combine density frames based on trip and trawl IDs
  # First, rename the catch column in each df
  densityTab <- densityTab %>% mutate(  catch = CATCH_WEIGHT,
                                        areaFished_m2 = DOORSPREAD_M * TOW_LENGTH_M,
                                        density = CATCH_WEIGHT / areaFished_m2 )

  includedSurveys <- unlist(stocks)

  # Save survey names for use in filtering code
  surveys     <- c("HS", "QCS", "WCHG", "WCVI")  
  synSurveys <- paste(surveys,"Syn",sep ="")

  # Split density tab into two parts
  # 1. HSAss
  # 2. Synoptic surveys
  HSAssTab <- densityTab %>%
              filter(SURVEY_SERIES_ID == surveyIDs["HSAss"])

  synTab   <- densityTab %>%
              filter( SURVEY_SERIES_ID %in% surveyIDs[synSurveys] )


  # Now from here, we want to take the data for each survey
  # and determine positive blocks, from which
  # we will make a new stratArea table, as well as identify
  # observations to remove
  for( gridIdx in 1:length(grids) )
  {
    # Run filterSurveyBlocks to produce plots of block
    # design and density ~ depth for all species/surveys
    filterSurveyBlocks( blocks = grids[[gridIdx]],
                        density = synTab,
                        plot = TRUE,
                        species = spec,
                        survey = surveys[gridIdx],
                        stratAreas = stratArea )  
  }
  



  # Now join and select the columns we want
  surveyData <- densityTab %>%
                left_join( stratArea, by = "GROUPING_CODE") %>%
                dplyr::select(  year = YEAR,
                                tripID = TRIP_ID,
                                eventID = FISHING_EVENT_ID,
                                lat = LATITUDE,
                                lon = LONGITUDE,
                                stratum = GROUPING_CODE,
                                majorArea = MAJOR_STAT_AREA_CODE,
                                minorArea = MINOR_STAT_AREA_CODE,
                                survey = SURVEY_DESC,
                                surveyID = SURVEY_ID,
                                survSeriesID = SURVEY_SERIES_ID,
                                stratArea = AREA_KM2,
                                density, catch,
                                fishedArea = areaFished_m2  ) %>%
                filter( survSeriesID %in% includedSurveys ) %>%
                mutate( stockName = sapply( X = survSeriesID, 
                                            FUN = appendName, 
                                            stocks ) )



  relativeBio <-  surveyData %>%
                  group_by( year, stockName, survSeriesID, stratum ) %>%
                  dplyr::summarise( density.var = var(density),
                                    density.mean = mean(density),
                                    area = mean(stratArea),
                                    fishedArea = sum(fishedArea),
                                    nBlocks = n() ) %>%
                  ungroup() %>%
                  filter( !is.na(density.mean) & !is.na(density.var) ) %>%
                  group_by( year, stockName, survSeriesID ) %>%
                  mutate( varMult = area * (area - nBlocks) / nBlocks,
                          varSummand = varMult * density.var ) %>%
                  summarise(  relBio_Kt = sum( area * density.mean ),
                              relBio_Kt.var = sum( varSummand ),
                              surveyArea = sum( area ) ) %>%
                  ungroup() %>%
                  filter( year >= years[1] & year <= years[2] ) %>%
                  mutate( surveyName = sapply(  X = survSeriesID, 
                                                FUN = appendName,
                                                survIDs ) )

  # relativeBio <-  surveyData %>%
  #                 group_by( year, stockName, survSeriesID, stratum ) %>%
  #                 dplyr::summarise( density.var = var(density),
  #                                   density.mean = mean(density),
  #                                   area = mean(stratArea),
  #                                   fishedArea = sum(fishedArea),
  #                                   nBlocks = n() ) %>%
  #                 dplyr::summarize( relBio_Kt = sum( area * density.mean),
  #                                   relBio.var_Kt = sum( area * (area - fishedArea) * density.var / fishedArea ),
  #                                   surveyArea = sum(area) ) %>%
  #                 ungroup() %>%
  #                 filter( year >= years[1] & year <= years[2] )

  # if( collapseSyn )
  # {
  #   # isolate synoptic IDs
  #   synIDs <-  c( QCSyn = 1, HSSyn = 3, WCVISyn=4,
  #                 WCHGSyn = 16 )
  #   synStratAreas <- surveyData %>%
  #                   filter( survSeriesID %in% synIDs ) %>%
  #                   group_by( stockName, survSeriesID, stratum ) %>%
  #                   summarise( stratArea = unique( stratArea ) ) %>%
  #                   summarise( surveyArea = sum( stratArea ) ) %>%
  #                   ungroup()

  #   synStockAreas <-  synStratAreas %>%
  #                     group_by( stockName ) %>%
  #                     summarise( stockArea = sum( surveyArea ) ) %>%
  #                     ungroup()


  #   synRelativeBio <- relativeBio %>% 
  #                     filter( survSeriesID %in% synIDs ) %>%
  #                     left_join( synStockAreas, by = "stockName" ) %>%
  #                     group_by( stockName, year ) %>%
  #                     dplyr::summarize( area = sum(surveyArea),
  #                                       relBio_Kt = sum( relBio_Kt ),
  #                                       relBio_Kt.var = sum( relBio_Kt.var ),
  #                                       stockArea = unique(stockArea) ) %>%
  #                     dplyr::mutate(  relBio_Kt = relBio_Kt * stockArea / area,
  #                                     survSeriesID = 1 ) %>%
  #                     dplyr::select(  year,
  #                                     stockName,
  #                                     survSeriesID,
  #                                     relBio_Kt,
  #                                     relBio_Kt.var,
  #                                     surveyArea = area,
  #                                     stockArea = stockArea
  #                                   )

  #   relativeBio <-  relativeBio %>%
  #                   left_join( synStockAreas, by = "stockName" ) %>%
  #                   filter( !(survSeriesID %in% synIDs) )

  #   relativeBio <- rbind(as.data.frame(relativeBio),as.data.frame(synRelativeBio))

  #   survIDs <- intersect(survIDs[c(1,2,5,6)],includedSurveys)

  # }

  yrs <- years[1]:years[2]
  stockNames <- names(stocks)
  relativeBioArray <- array(  -1, dim = c(length(survIDs),length(stocks),length(yrs),2),
                              dimnames = list(  names(survIDs), stockNames,
                                                yrs, c("Mean", "SE") ) )

  for( survIdx  in 1:length(survIDs) )
  {
    for( yIdx in 1:length(yrs) )
    {
      for( sIdx in 1:length(stocks) )
      {
        stockID <- names(stocks)[sIdx]
        subRelBio <-  relativeBio %>%
                      filter( year == yrs[yIdx], survSeriesID == survIDs[survIdx], stockName == stockID )
        if( nrow(subRelBio) == 0 ) next
        if( nrow(subRelBio) > 1 ) browser()
        relativeBioArray[survIdx, sIdx, yIdx, "Mean" ]  <- subRelBio$relBio_Kt
        relativeBioArray[survIdx, sIdx, yIdx, "SE" ]    <- sqrt(subRelBio$relBio_Kt.var)  
      }
      
    }
  }

  out <- list(  relBio.df = relativeBio,
                relBio.arr = relativeBioArray)
  out
}

appendName <- function( dataCode, codeKey )
{
  codeKeyVec <- unlist(codeKey)
  nameVec <- c()
  for( k in 1:length(codeKey))
    nameVec <- c(nameVec,rep(names(codeKey)[k],length(codeKey[[k]])))

  stockCodeNum <- which( dataCode == codeKeyVec )

  if( length(stockCodeNum) > 0 )
    return(nameVec[stockCodeNum])
  else
    return(NA)
  
}


# Read in commercial CPUE, arrange in a data
# frame and an array for feeding to TMB.
# Need to decide if the commercial CPUE
# will be a single fleet with blocks for
# modern/historic catch, or 2 fleets 
# (or more blocks/fleets).
# inputs:   specName = species name char vector (root of filename)
#           stocks = character list of stock IDs for converting
#                     columns in data
readCommCPUE <- function( specName = "dover-sole",
                          stocks = stocksCommCPUE,
                          years = c(1954,2018) )
{
  # Set data path
  datPath       <- file.path(getwd(),"Data","comm-cpue-flatfish")
  modernName    <- paste("cpue-predictions-",specName,"-modern.csv",sep = "")
  historicName  <- paste("cpue-predictions-",specName,"-historic.csv",sep = "")
  # Read modern and historic data

  modData   <-  read.csv( file.path(datPath, modernName), header = T, 
                          stringsAsFactors = FALSE ) %>%
                mutate( period = "modern" ) %>%
                dplyr::select(  version = "formula_version",
                                est_link, se_link, area, year,
                                period ) %>%
                mutate( stockName = sapply( X = area, 
                                            FUN = appendName, 
                                            codeKey = stocks ) )

  histData  <-  read.csv( file.path(datPath, historicName), header = T,
                          stringsAsFactors = FALSE ) %>%
                mutate( period = "historic" ) %>%
                dplyr::select(  version = "formula_version",
                                est_link, se_link, area, year,
                                period ) %>%
                mutate( stockName = sapply( X = area, 
                                            FUN = appendName, 
                                            codeKey = stocks ) )

  # Combine data for DF returning
  allData <- rbind( modData, histData )

  # Pull years range
  modYears  <- range(modData$year)
  histYears <- range(histData$year)

  # create a years vector
  yrs <- years[1]:years[2]
  histYrs <- histYears[1]:histYears[2]
  modYrs  <- modYears[1]:modYears[2] 

  dataArray <- array( -1, dim = c(  2,               # Fleets (modern/hist)
                                    length(stocks),  # Number of stocks
                                    length(yrs),     # Number of time steps
                                    2,               # stdized/unstdized 
                                    2 ),             # Mean/SE
                          dimnames = list(  c("comm.hist","comm.mod"),
                                            names(stocks),
                                            yrs,
                                            c("log.mean","log.sd"),
                                            c("stdized","unstdized") ) )  

  for( stockIdx in 1:length(stocks) )
  {
    stockName <- names(stocks)[stockIdx]
    stockID   <- stocks[stockIdx]
    
    # Filter data by stock
    subModStd <-  modData %>%
                  filter( area == stockID,
                          version == "Full standardization")

    subModUnStd <-  modData %>%
                    filter( area == stockID,
                            version == "Unstandardized")

    subHistStd <- histData %>%
                  filter( area == stockID,
                          version == "Full standardization")

    subHistUnStd <- histData %>%
                    filter( area == stockID,
                            version == "Unstandardized" )

    # Place in array
    # historic data, mean values
    dataArray["comm.hist",stockName,as.character(histYrs),"log.mean","stdized"] <- subHistStd$est_link
    dataArray["comm.hist",stockName,as.character(histYrs),"log.mean","unstdized"] <- subHistUnStd$est_link
    # comm.hist data, standard errors
    dataArray["comm.hist",stockName,as.character(histYrs),"log.sd","stdized"] <- subHistStd$se_link
    dataArray["comm.hist",stockName,as.character(histYrs),"log.sd","unstdized"] <- subHistUnStd$se_link
    # Modern data, mean values
    dataArray["comm.mod",stockName,as.character(modYrs),"log.mean","stdized"] <- subModStd$est_link
    dataArray["comm.mod",stockName,as.character(modYrs),"log.mean","unstdized"] <- subModUnStd$est_link
    # comm.Mod data, standard errors
    dataArray["comm.mod",stockName,as.character(modYrs),"log.sd","stdized"] <- subModStd$se_link
    dataArray["comm.mod",stockName,as.character(modYrs),"log.sd","unstdized"] <- subModUnStd$se_link
  }

  return( list( cpue.df = allData,
                cpue.arr = dataArray ) ) 
}

# Quick function to scale by mean values
transByMean <- function(x)
{
  x <- x - mean(x,na.rm = T)

  x
}


# Read in bio data
readBioData <- function(  specName = "dover",
                          stocksSurv = stocksSurvey,
                          stocksComm = stocksCommBio,
                          years = c(1954,2018) )
{
  # First, read in survey bio data
  surveyDataName <- paste(specName,"_bio.csv",sep = "")
  surveyDataPath <- file.path(getwd(),"Data","biology_with_mat",surveyDataName)
  surveyBio      <- read.csv(surveyDataPath, header = TRUE, stringsAsFactors = FALSE )
  # Read in commercial bio data
  commDataName <- paste(specName,"_comm_biodata.csv",sep = "")
  commDataPath <- file.path(getwd(),"Data","DERPA_commercial_biodata",commDataName)
  commBio      <- read.csv(commDataPath, header = TRUE, stringsAsFactors = FALSE )

  # Read in density
  specDensityFile <- paste(specName,"density.csv",sep = "_")
  specDensityPath <- file.path(getwd(),"Data","density",specDensityFile)
  trawlInfoTab <- read.csv(specDensityPath, header = T) %>%
                  dplyr::select(  TRIP_ID, FE_MAJOR_LEVEL_ID,
                                  MAJOR_STAT_AREA_CODE, SURVEY_SERIES_ID,
                                  DEPTH_M, YEAR )
  # Join trawl info to survey bio data and
  # append stock names - I think that's all we want to do
  # for now
  surveyBio <-  surveyBio %>%
                left_join( trawlInfoTab, by = c("TRIP_ID", "FE_MAJOR_LEVEL_ID" ) ) %>%
                mutate( stockName = sapply( X = SURVEY_SERIES_ID, 
                                            FUN = appendName,
                                            codeKey = stocksSurv),
                        stockName = unlist(stockName) )

  
  # Now do the same for the commercial data
  commBio <-  commBio %>%
              mutate( stockName = sapply( X = MAJ, 
                                          FUN = appendName,
                                          codeKey = stocksComm ) )


  outList <- list(  survey = surveyBio,
                    comm = commBio )

  return(outList)
} # END readBioData()

# Define a NLL for optim
vonB_nll <- function( theta,
                      data )
{
  # Recover leading pars
  Linf  <- exp(theta[1])
  vonK  <- exp(-exp(theta[2]))
  L1    <- exp(theta[3])
  cvL   <- 1 / (1 + exp(-theta[4]))

  data <- data %>%
          mutate( expLt = vonB(age,theta = theta),
                  res   = expLt - length,
                  sigL  = cvL * expLt,
                  lik   = .5*(log(cvL^2) + (res/sigL/cvL)^2 ) )
  L50       <- quantile( data$length, probs = c(0.5) )
  
  # Get data to generate prior mean for Linf
  LinfData  <-  data %>% 
                filter( length > L50 ) 
  # Get data L1 for prior mean
  L1data    <-  data %>%
                filter(age == 1)


  # Add likelihood
  nll <- sum( data$lik, na.rm = TRUE ) 
  # add prior for Linf
  if( nrow(LinfData) > 0)
  {
    meanLinf  <- mean(LinfData$length)
    nll       <- nll  + 0.5*((Linf - meanLinf)/meanLinf)^2
  }
  # add prior for L1
  if(nrow(L1data) > 0)
  {
    meanL1  <- mean(L1data$length)
    nll     <- nll  + 0.5*((L1 - meanL1)/meanL1)^2
  }

  return(nll)
}

# von bertalanffy growth function
vonB <- function( age, vonK = NULL, Linf = NULL, L1 = NULL, theta = NULL )
{
  if(!is.null(theta))
  {
    Linf <- exp(theta[1])
    vonK <- exp(-exp(theta[2]))
    L1 <- exp(theta[3])
  }

  Lt <- Linf + (L1 - Linf) * exp( -vonK * (age - 1))
  Lt
}


# makeLenAge()
# Conducts the length at age analysis for a given set of data,
# splitting over stocks and estimating vonB functions
# for each stock area, species, and sex. Spits out a length
# at age frequency array (species,stock,sex)
# inputs: data = output of readBioData()
# outputs:  vonB = list of vonB fits, 
#           data = input data  
#           ALfreq = array of age-length freq, indexed by 
#                     stock and species
makeLenAge <- function( data = bioData$Dover,
                        stocks = names(stocksCommBio) )
{
  # Get survey and commercial data
  survData <- data$survey
  commData <- data$comm

  # We want to filter down to cases where there are age/length
  # observations
  survData <- survData %>%
              dplyr::select(  age = AGE, length = LENGTH_MM,
                              stockName, sex = SEX ) %>%
              mutate( length = length/10 ) %>%
              filter( !is.na(age),
                      !is.na(length),
                      length > 0 )
  commData <- commData %>%
              dplyr::select(  age = AGE, length = LENGTH_MM,
                              stockName, sex = SEX ) %>%
              filter( !is.na(age),
                      !is.na(length),
                      length > 0 ) %>%
              mutate( length = length/10 )

  data$survey <- survData
  data$comm   <- commData 

  # combine total data for likelihood
  combData <- rbind(survData,commData)

  bioDataBoys   <- combData %>% filter( sex == 1 )
  bioDataGirls  <- combData %>% filter( sex == 2 )

  # Now optimise model to create vonB fits
  Linf <- max(combData$length)
  L1init <- min(combData$length)
  theta <- c(log(Linf),0,log(L1init),0)
  lenAgeAll   <- optim(par = theta, fn = vonB_nll, data = combData, method = "Nelder-Mead" )
  lenAgeAll   <- optim(par = lenAgeAll$par, fn = vonB_nll, data = combData, method = "BFGS" )
  lenAgeBoys  <- optim(par = lenAgeAll$par, fn = vonB_nll, data = bioDataBoys, method = "Nelder-Mead" )
  lenAgeBoys  <- optim(par = lenAgeBoys$par, fn = vonB_nll, data = bioDataBoys, method = "BFGS" )
  lenAgeGirls <- optim(par = lenAgeAll$par, fn = vonB_nll, data = bioDataGirls, method = "Nelder-Mead" )
  lenAgeGirls <- optim(par = lenAgeGirls$par, fn = vonB_nll, data = bioDataGirls, method = "BFGS" )

  # save the coastwide model
  coastWide <- list(  lenAgeAll = lenAgeAll,
                      lenAgeBoys = lenAgeBoys,
                      lenAgeGirls = lenAgeGirls )

  # initialise age-length freq
  # First, get array dimensions
  lengths <- unique(combData$length)
  ages    <- unique(combData$age)

  allAges     <- 1:max(ages)
  allLengths  <- 1:max(lengths)

  # initialise array
  ALfreq <- array( NA,  dim = c(  length(stocks),
                                  length(allAges),
                                  length(allLengths),
                                  2 ),
                        dimnames = list(  stocks, 
                                          allAges, 
                                          allLengths, 
                                          c("boys","girls") ) )

  stockFits <- vector( mode = "list", length = length(stocks) )
  names(stockFits) <- stocks

  # Now loop over stocks
  for( stockIdx in 1:length(stocks) )
  {
    stockID <- stocks[stockIdx]
    stockData <-  combData %>% 
                  filter( stockName == stockID )

    stockDataBoys <- stockData %>% filter( sex == 1 )
    stockDataGirls <- stockData %>% filter( sex == 2 )

    if( nrow(stockDataBoys) == 0 | nrow(stockDataGirls) == 0 )
      next

    lenAgeAll   <- optim( par = lenAgeAll$par, fn = vonB_nll, data = stockData, method = "BFGS" )
    lenAgeBoys  <- optim( par = lenAgeBoys$par, fn = vonB_nll, data = stockDataBoys, method = "BFGS" )
    lenAgeGirls <- optim( par = lenAgeGirls$par, fn = vonB_nll, data = stockDataGirls, method = "BFGS" )

    stockFits[[stockIdx]] <- list(  lenAgeAll = lenAgeAll,
                                    lenAgeBoys = lenAgeBoys,
                                    lenAgeGirls = lenAgeGirls )
    # Now populate the array - round length to nearest cm
    # might need to aggregate further
    frqData <-  stockData %>%
                mutate( length = round(length) ) %>%
                group_by( age, length, sex ) %>%
                filter( sex %in% c(1,2) ) %>%
                summarise( nObs = n() ) %>%
                ungroup()

    for( a in allAges )
      for( l in allLengths )
      {
        boyObs <- frqData %>% filter( sex == 1, age == a, length == l )
        girlObs <- frqData %>% filter( sex == 2, age == a, length == l )

        if(nrow(boyObs) > 0)
          ALfreq[ stockIdx, a, l, "boys" ] <- boyObs$nObs
        if(nrow(girlObs) > 0)
          ALfreq[ stockIdx, a, l, "girls" ] <- girlObs$nObs

      }
  }

  outList <- list(  data = data,
                    vonB = list( coastWide = coastWide, stocks = stockFits ),
                    ALfreq = ALfreq )

  outList
} # END makeLenAge()

# plotLenAge()
# Plots vonBertalanffy growth models for each species
# and stock from the output of makeLenAge().
plotLenAge <- function( obj = lenAge,
                        save = FALSE,
                        saveDir = "Outputs",
                        stocks = names(stocksCommBio) )
{
  # Set up plotting area
  nSpec   <- length(obj)
  nStocks <- length(stocks)

  if(save)
  {
    graphics.off()
    saveDir <- file.path(getwd(),saveDir)
    if(! dir.exists(saveDir) )
      dir.create(saveDir)
    fileName <- "lengthAtAge.png"

    savePath <- file.path(saveDir, fileName )
    png(  savePath, width = 11, height = 8.5, units = "in",
          res = 300)
  }

  par(  mfcol = c(nStocks, nSpec ), 
        mar = c(0,2,0,2),
        oma = c(7,3,2,2) )

  sexCols       <- brewer.pal( n = 3, "Dark2" )
  sexPch        <- c(3,1,2,3)

  # Loop over species, pull data and fits
  # for that species
  for( specIdx in 1:length(obj))
  {
    specName      <- names(obj)[specIdx]
    specData      <- obj[[specIdx]]$data
    vonB          <- obj[[specIdx]]$vonB

    specSurv      <- specData$survey
    specComm      <- specData$comm

    maxAge        <- max(specSurv$age, specComm$age,na.rm = T)
    maxLen        <- max(specSurv$length, specComm$length, na.rm =T)

    coastWide     <- vonB$coastWide

    vonBAll.cw    <- vonB( age = 1:maxAge, theta = coastWide$lenAgeAll$par )
    vonBBoys.cw   <- vonB( age = 1:maxAge, theta = coastWide$lenAgeBoys$par )
    vonBGirls.cw  <- vonB( age = 1:maxAge, theta = coastWide$lenAgeGirls$par )

    # Then loop over stocks
    for( stockIdx in 1:length(stocks))
    {
      stockID       <- stocks[stockIdx]
      subData.surv  <- specSurv %>% filter( stockName == stockID )
      subData.comm  <- specComm %>% filter( stockName == stockID )

      stockFit      <- vonB$stocks[[stockID]]

      plot( x = c(1,maxAge),
            y = c(0,maxLen),
            type = "n", axes = F, xlab = "", ylab = "")
        mfg <- par("mfg")
        axis(side = 2, las = 1)
        if( mfg[1] == mfg[3] )
          axis( side = 1 )
        if( mfg[1] == 1 )
          mtext( side = 3, text = specName, font = 2 )
        if(mfg[2] == mfg[4] )
          mtext( side = 4, text = stockID, line = 1.5 )
        box()

        if( nrow( subData.surv ) != 0 | nrow(subData.comm) != 0)
        {
          vonBparsAll.stock <- stockFit$lenAgeAll$par
          vonBparsBoys.stock <- stockFit$lenAgeBoys$par
          vonBparsGirls.stock <- stockFit$lenAgeGirls$par

          Linf.stock  <- round(exp(c(vonBparsAll.stock[1], vonBparsBoys.stock[1], vonBparsGirls.stock[1]) ),2)
          K.stock     <- round(exp(-1*exp(c(vonBparsAll.stock[2], vonBparsBoys.stock[2], vonBparsGirls.stock[2]) )),2)
          L1.stock    <- round(exp(c(vonBparsAll.stock[3], vonBparsBoys.stock[3], vonBparsGirls.stock[3] )),2)

          vonBAll       <- vonB( age = 1:maxAge, theta = stockFit$lenAgeAll$par )
          vonBBoys      <- vonB( age = 1:maxAge, theta = stockFit$lenAgeBoys$par )
          vonBGirls     <- vonB( age = 1:maxAge, theta = stockFit$lenAgeGirls$par )
          points( x = subData.surv$age, y = subData.surv$length,
                  pch = sexPch[subData.surv$sex + 1], col = alpha("grey70",alpha = .3),
                  cex = .5 )
          points( x = subData.comm$age, y = subData.comm$length,
                  pch = sexPch[subData.comm$sex + 1], col = alpha("grey70",alpha = .3),
                  cex = .5 )
          lines( x = 1:maxAge, y = vonBAll, lty = 1, lwd = 3, col = "black" )
          lines( x = 1:maxAge, y = vonBBoys, lty = 1, lwd = 3, col = sexCols[1] )
          lines( x = 1:maxAge, y = vonBGirls, lty = 1, lwd = 3, col = sexCols[2] )
          text(  x = c(0.4,0.55,0.7,0.85)*maxAge, y = 0.25*maxLen,
                 label = c(" ", "All", "Boys", "Girls" ), cex = .6 )
          text(  x = c(0.4,0.55,0.7,0.85)*maxAge, y = 0.2*maxLen,
                 label = c("Linf ", Linf.stock ), cex = .6 )
          text(  x = c(0.4,0.55,0.7,0.85)*maxAge, y = 0.15*maxLen,
                 label = c("vonK ", K.stock ), cex = .6 )
          text(  x = c(0.4,0.55,0.7,0.85)*maxAge, y = 0.1*maxLen,
                 label = c("L1 ", L1.stock ), cex = .6 )
        } else {
          panLab( x = 0.3, y = 0.9, 
                  txt = "No Data" )
        }
        # Plot coastwide models
        lines( x = 1:maxAge, y = vonBAll.cw, lty = 2, lwd = 1, col = "grey70" )
        lines( x = 1:maxAge, y = vonBBoys.cw, lty = 2, lwd = 1, col = sexCols[1] )
        lines( x = 1:maxAge, y = vonBGirls.cw, lty = 2, lwd = 1, col = sexCols[2] )

    }
  }
  mtext( side = 1, text = "Age", outer = T, line = 3)
  mtext( side = 2, text = "Length (cm)", outer = T, line = 1.5)
  # Add a legend
  par( mfcol = c(1,1), oma = c(0,3,2,2)  )
  legend( x = "bottom",
          horiz = TRUE,
          legend = c("Male", "Female", "Unsexed", "Coastwide"),
          pch = c( 1, 2, 3, NA ),
          lty = c( 1, 1, 1, 2 ),
          lwd = c( 2, 2, 2, 1 ),
          pt.lwd = c( 2, 2, 2, NA),
          col = c( sexCols[1:2], "black", "grey70" ),
          bty = "n" )
  if( save )
  {
    dev.off()
    cat( "Length-at-age plots saved to ", savePath, "\n", sep = "")
  }
}

# makeWtLen()
# Conducts the weight at length analysis for a given set of data,
# splitting over stocks and estimating allometric 
# growth parameters functions for each stock area, 
# species, and sex.
# inputs: data = output of readBioData()
# outputs:  wtLen = list of wtLen fits, 
#           data = input data  
makeWtLen <- function(  data = bioData$English,
                        stocks = names(stocksCommBio) )
{
  # Get survey and commercial data
  survData <- data$survey
  commData <- data$comm

  # We want to filter down to cases where there are age/length
  # observations
  survData <- survData %>%
              dplyr::select(  weight = WEIGHT_G, length = LENGTH_MM,
                              stockName, sex = SEX ) %>%
              mutate( length = length / 10,
                      weight = weight / 1000,
                      logL = log(length),
                      logW = log(weight) ) %>%
              filter( !is.na(weight),
                      !is.na(length),
                      length > 0,
                      weight > 0 )
  commData <- commData %>%
              dplyr::select(  weight = WEIGHT_G, length = LENGTH_MM,
                              stockName, sex = SEX ) %>%
              filter( !is.na(weight),
                      !is.na(length),
                      length > 0,
                      weight > 0 ) %>%
              mutate( length = length/10,
                      weight = weight/1000,
                      logL = log(length),
                      logW = log(weight) )

  data$survey <- survData
  data$comm   <- commData 

  # combine total data for likelihood
  combData <- rbind(survData,commData)

  bioDataBoys   <- combData %>% filter( sex == 1 )
  bioDataGirls  <- combData %>% filter( sex == 2 )

  # Now optimise model to create vonB fits
  wtLenAll    <- lm( logW ~ logL, data = combData )
  wtLenBoys   <- lm( logW ~ logL, data = bioDataBoys )
  wtLenGirls  <- lm( logW ~ logL, data = bioDataGirls )


  # save the coastwide model
  coastWide <- list(  wtLenAll = wtLenAll,
                      wtLenBoys = wtLenBoys,
                      wtLenGirls = wtLenGirls )

  stockFits <- vector( mode = "list", length = length(stocks) )
  names(stockFits) <- stocks

  # Now loop over stocks
  for( stockIdx in 1:length(stocks) )
  {
    stockID <- stocks[stockIdx]
    stockData <-  combData %>% 
                  filter( stockName == stockID )

    stockDataBoys <- stockData %>% filter( sex == 1 )
    stockDataGirls <- stockData %>% filter( sex == 2 )

    if( nrow(stockDataBoys) == 0 & nrow(stockDataGirls) == 0 )
      next
    

    wtLenAll    <- lm( logW ~ logL, data = stockData )
    if( nrow(stockDataBoys) > 0 )
      wtLenBoys   <- lm( logW ~ logL, data = stockDataBoys )
    if( nrow(stockDataGirls) > 0)
      wtLenGirls  <- lm( logW ~ logL, data = stockDataGirls )

    stockFits[[stockIdx]] <- list(  wtLenAll = wtLenAll,
                                    wtLenBoys = wtLenBoys,
                                    wtLenGirls = wtLenGirls )
  }

  outList <- list(  data = data,
                    wtLen = list( coastWide = coastWide, stocks = stockFits ) )

  outList
} # END makeLenAge()

# plotLenAge()
# Plots vonBertalanffy growth models for each species
# and stock from the output of makeLenAge().
plotWtLen <- function(  obj = wtLen,
                        save = FALSE,
                        saveDir = "Outputs",
                        stocks = names(stocksCommBio) )
{
  # Set up plotting area
  nSpec   <- length(obj)
  nStocks <- length(stocks)

  if(save)
  {
    graphics.off()
    saveDir <- file.path(getwd(),saveDir)
    if(! dir.exists(saveDir) )
      dir.create(saveDir)
    fileName <- "wtLength.png"

    savePath <- file.path(saveDir, fileName )
    png(  savePath, width = 11, height = 8.5, units = "in",
          res = 300)
  }

  par(  mfcol = c(nStocks, nSpec ), 
        mar = c(0,2,0,2),
        oma = c(7,3,2,2)  )

  sexCols       <- brewer.pal( n = 3, "Dark2" )
  sexPch        <- c(3,1,2,3)
  # sexCols.alpha <- alpha( sexCols, alpha = .2 )

  # Loop over species, pull data and fits
  # for that species
  for( specIdx in 1:length(obj))
  {
    # First make coastwide models
    specName      <- names(obj)[specIdx]
    specData      <- obj[[specIdx]]$data
    wtLen         <- obj[[specIdx]]$wtLen

    specSurv      <- specData$survey
    specComm      <- specData$comm

    maxWt         <- max(specSurv$weight, specComm$weight,na.rm = T)
    maxLen        <- max(specSurv$length, specComm$length, na.rm =T)

    coastWide     <- wtLen$coastWide

    lenVec        <- seq(1,maxLen,length = 200)

    wtLenAll.cw   <- exp(coef(coastWide$wtLenAll)[1]) * lenVec^(coef(coastWide$wtLenAll)[2])
    wtLenBoys.cw  <- exp(coef(coastWide$wtLenBoys)[1]) * lenVec^(coef(coastWide$wtLenBoys)[2])
    wtLenGirls.cw <- exp(coef(coastWide$wtLenGirls)[1]) * lenVec^(coef(coastWide$wtLenGirls)[2])


    # Then loop over stocks
    for( stockIdx in 1:length(stocks))
    {
      stockID       <- stocks[stockIdx]
      
      subData.surv  <- specSurv %>% filter( stockName == stockID )
      subData.comm  <- specComm %>% filter( stockName == stockID )

      stockFit      <- wtLen$stocks[[stockID]]

      plot( x = c(1,maxLen),
            y = c(0,maxWt),
            type = "n", axes = F, xlab = "", ylab = "")
        mfg <- par("mfg")
        axis(side = 2, las = 1)
        if( mfg[1] == 1)
          mtext( side = 3, text = specName, font = 2, line = .5 )
        if( mfg[1] == mfg[3] )
          axis( side = 1 )
        if( mfg[2] == mfg[4] )
          mtext( side = 4, text = stockID, line = 1.5 )
        box()

        if( nrow( subData.surv ) != 0 | nrow(subData.comm) != 0)
        {
          wtLenAll.stock   <- coef(stockFit$wtLenAll)
          wtLenBoys.stock  <- coef(stockFit$wtLenBoys)
          wtLenGirls.stock <- coef(stockFit$wtLenGirls)

          a.stock <- round(exp(c(wtLenAll.stock[1],wtLenGirls.stock[1],wtLenBoys.stock[1])),2)
          b.stock <- round(c(wtLenAll.stock[2],wtLenGirls.stock[2],wtLenBoys.stock[2]),2)

          wtLenAll       <- exp(wtLenAll.stock[1]) * lenVec^wtLenAll.stock[2]
          wtLenBoys      <- exp(wtLenBoys.stock[1]) * lenVec^wtLenBoys.stock[2]
          wtLenGirls     <- exp(wtLenGirls.stock[1]) * lenVec^wtLenGirls.stock[2]
          points( x = subData.surv$length, y = subData.surv$weight,
                  pch = sexPch[subData.surv$sex + 1], col = alpha("grey70", alpha = .5) )
          points( x = subData.comm$length, y = subData.comm$weight,
                  pch = sexPch[subData.comm$sex + 1], col = alpha("grey70", alpha = .5) )
          lines( x = lenVec, y = wtLenAll, lty = 1, lwd = 3, col = "black" )
          lines( x = lenVec, y = wtLenBoys, lty = 1, lwd = 3, col = sexCols[1] )
          lines( x = lenVec, y = wtLenGirls, lty = 1, lwd = 3, col = sexCols[2] )
          text(  x = c(0.1,0.25,0.4,0.55)*maxLen, y = 0.95*maxWt,
                 label = c(" ", "All", "Boys", "Girls" ), cex = .8 )
          text(  x = c(0.1,0.25,0.4,0.55)*maxLen, y = 0.9*maxWt,
                 label = c("a ", a.stock ), cex = .8 )
          text(  x = c(0.1,0.25,0.4,0.55)*maxLen, y = 0.85*maxWt,
                 label = c("b ", b.stock ), cex = .8 )
        } else {
          panLab( x = 0.3, y = 0.9, 
                  txt = "No Data" )
        }
        # Plot coastwide models
        lines( x = lenVec, y = wtLenAll.cw, lty = 2, lwd = 1, col = "grey70" )
        lines( x = lenVec, y = wtLenBoys.cw, lty = 2, lwd = 1, col = sexCols[1] )
        lines( x = lenVec, y = wtLenGirls.cw, lty = 2, lwd = 1, col = sexCols[2] )

    }
  }
  mtext( side = 1, text = "Length (cm)", outer = T, line = 3)
  mtext( side = 2, text = "Weight (kg)", outer = T, line = 1.5)
  # Add a legend
  par( mfcol = c(1,1), oma = c(0,3,2,2)  )
  legend( x = "bottom",
          horiz = TRUE,
          legend = c("Male", "Female", "Unsexed", "Coastwide"),
          pch = c( 1, 2, 3, NA ),
          lty = c( 1, 1, 1, 2 ),
          lwd = c( 2, 2, 2, 1 ),
          pt.lwd = c( 2, 2, 2, NA),
          col = c( sexCols[1:2], "black", "grey70" ),
          bty = "n" )

  if( save )
  {
    dev.off()
    cat( "Weight-at-length plots saved to ", savePath, "\n", sep = "")
  }
} # END plotWtLen()


# plotIndices()
# Plots the commercial and survey biomass indices
# on a multipanel plot
plotIndices <- function(  survey = relBioList_Survey,
                          comm = commCPUEList,
                          stocksSurv = stocksSurvey,
                          stocksComm = stocksCommData,
                          save = FALSE,
                          saveDir = "Outputs" )
{
  # Get number of species
  nSpec <- length(survey)
  nStocks <- length( stocksSurv )

  # Get years from data
  survDimNames <- dimnames(survey[[1]]$relBio.arr)
  commDimNames <- dimnames(comm[[1]]$cpue.arr)
  yrs <- as.integer(survDimNames[[3]])

  vertLines <- seq(1960,2020, by = 10)

  # Survey and comm fleet names
  surveys <- survDimNames[[1]]
  commPeriods <- commDimNames[[1]]

  # Get colours for plotting
  cols <- brewer.pal( n = length(surveys) + length(commPeriods), "Dark2" )

  if( save )
  {
    graphics.off()
    saveDir <- file.path(getwd(),saveDir)
    if(! dir.exists(saveDir) )
      dir.create(saveDir)
    fileName <- "stockIndices.png"

    savePath <- file.path(saveDir, fileName )
    png(  savePath, width = 11, height = 8.5, units = "in",
          res = 300)
  }

  # Set up plotting window
  par(  mfcol = c(nStocks, nSpec), 
        oma = c(4,4,3,4), 
        mar = c(0,0,0,0) )

  # Now loop over species and stocks, plot each one's
  # indices
  for( specIdx in 1:nSpec )
  {
    # Pull survey data
    relBio    <-  survey[[specIdx]]$relBio.df %>%
                  group_by( stockName, surveyName ) %>%
                  mutate( logRelBio_Kt = log(relBio_Kt),
                          logRelBio_Kt.scaled = logRelBio_Kt - mean(logRelBio_Kt),
                          relBio_Kt.scaled = exp(logRelBio_Kt.scaled),
                          relBio_Kt.se = sqrt(relBio_Kt.var),
                          relBio_Kt.lwr = relBio_Kt - relBio_Kt.se,
                          relBio_Kt.upr = relBio_Kt + relBio_Kt.se,
                          relBio_Kt.lwr.scaled = exp(log(relBio_Kt.lwr) - mean(logRelBio_Kt)),
                          relBio_Kt.upr.scaled = exp(log(relBio_Kt.upr) - mean(logRelBio_Kt)) ) %>%
                  ungroup()

                  

    # Pull commercial CPUE
    commCPUE  <-  comm[[specIdx]]$cpue.df %>%
                  group_by( version, stockName, period) %>%
                  mutate( est_link.scaled = est_link - mean(est_link),
                          est.scaled = exp(est_link.scaled),
                          est.scaled.upr = exp(est_link.scaled + se_link),
                          est.scaled.lwr = exp(est_link.scaled - se_link) ) %>%
                  ungroup()

    for( stockIdx in 1:nStocks )
    {
      stockID   <- names(stocksSurv)[stockIdx]
      specName  <- names(survey)[specIdx]

      # Filter down to stock specific data
      relBio.stock  <-  relBio %>%
                        filter( stockName == stockID )
      commCPUE.stock <- commCPUE %>%
                        filter( stockName == stockID )


      # Start plot
      plot( x = range(yrs), y = range(0,6),
            type = "n", axes = F, xlab = "", ylab = "" )
        abline(v = vertLines, lty = 3, lwd = 1, col = "grey75")
        mfg <- par("mfg")
        if( mfg[1] == mfg[3])
          axis( side = 1 )
        if( mfg[2] == 1)
          axis( side = 2, las = 1 )
        if(mfg[1] == 1)
          mtext(side = 3, text = specName, line = 1.5, font = 2 )
        if(mfg[2] == mfg[4])
          mtext( side = 4, text = stockID, line = 1.5)
        box()
        # First plot stdized and unstdized commercial CPUE
        for( commIdx in 1:length(commPeriods))
        {
          commPeriod      <-  commPeriods[commIdx]
          subCommData.std <-  commCPUE.stock %>%
                              filter( period == commPeriod,
                                      version == "Full standardization" )
          subCommData.ust <-  commCPUE.stock %>%
                              filter( period == commPeriod,
                                      version == "Unstandardized" )
          # Plot standardized
          
          segments( x0 = subCommData.std$year - 0.1, x1 = subCommData.std$year - 0.1,
                    y0 = subCommData.std$est.scaled.lwr,
                    y1 = subCommData.std$est.scaled.upr, col = cols[commIdx], lwd = 2 )
          points( x = subCommData.std$year - 0.1,
                  y = subCommData.std$est.scaled,
                  col = cols[commIdx], pch = 21 )
          # Plot unstandardized 
          segments( x0 = subCommData.ust$year + 0.1, x1 = subCommData.ust$year + 0.1,
                    y0 = subCommData.ust$est.scaled.lwr,
                    y1 = subCommData.ust$est.scaled.upr, col = cols[commIdx], 
                    lwd = 2, lty = 2 )
          points( x = subCommData.ust$year + 0.1,
                  y = subCommData.ust$est.scaled,
                  col = cols[commIdx], pch = 22 )

        }

        # Now plot each survey
        for( survIdx in 1:length(surveys) )
        {
          colIdx <- length(commPeriods) + survIdx 

          survName <- surveys[survIdx]
          subRelBio.stock <-  relBio.stock %>%
                              filter( surveyName == survName )
          if(nrow(subRelBio.stock) == 0 )
            next

          segments( x0 = subRelBio.stock$year,
                    x1 = subRelBio.stock$year,
                    y0 = subRelBio.stock$relBio_Kt.lwr.scaled,
                    y1 = subRelBio.stock$relBio_Kt.upr.scaled,
                    lwd = 2, col = cols[colIdx])
          points( x = subRelBio.stock$year,
                  y = subRelBio.stock$relBio_Kt.scaled,
                  pch = 15 + survIdx, col = cols[colIdx] )
        }

    }

  }
  mtext(  side = 1, text = "Year", outer = TRUE,
          line = 2 )
  mtext(  side = 2, text = "Indices scaled by geometric mean", outer = TRUE,
          line = 2 )
  if(save)
  {
    dev.off()
    cat("Stock indices plot saved to ", savePath,"\n",sep = "")
  }
}

# makeCatchDiscArrays()
# Creates catch and dicard data arrays for hierSCAL model
makeCatchDiscArrays <- function(  commData = catchData,
                                  survData = surveyCatch,
                                  stocks = stocksCommBio,
                                  speciesCodes = specCodes,
                                  years = fYear:lYear,
                                  modernYear = 1996,
                                  nF = 7, collapseComm = FALSE,
                                  fleetIDs = NULL )
{
  commData <- commData %>%
              filter( FISHERY_SECTOR == "GROUNDFISH TRAWL" ) %>%
              mutate( stockName = sapply( X = MAJOR_STAT_AREA_CODE,
                                          FUN = appendName,
                                          stocks ),
                      species = sapply( X = SPECIES_CODE,
                                        FUN = appendName,
                                        speciesCodes ) ) %>%
              filter( stockName %in% names(stocks) ) %>%
              dplyr::select(  year          = YEAR,
                              catch         = LANDED_KG,
                              discardWt     = DISCARDED_KG,
                              discardNum  = DISCARDED_PCS,
                              species, stockName ) %>%
              group_by( species, stockName, year ) %>%
              summarise(  catch   = sum(catch)/1e6,
                          discWt  = sum(discardWt)/1e6,
                          discNum = sum(discardNum) ) %>%
              ungroup()

  # Array dimensions
  nS <- length(speciesCodes)
  nP <- length(stocks)
  nT <- length(years)

  specIDs <- names(speciesCodes)
  stockIDs <- names(stocks)

  if( collapseComm )
    nF <- nF - 1

  # initialise arrays
  C_spft <- array( 0, dim = c(nS, nP, nF, nT),
                      dimnames = list(  species = specIDs,
                                        stocks = stockIDs,
                                        fleets = fleetIDs,
                                        years = years ) )

  D_spft <- C_spft


  # Now loop over species
  for( specIdx in 1:nS )
  {
    specID <- specIDs[specIdx]
    specSurvData <- survData[[specID]] %>%
                    filter( surveyName %in% fleetIDs )
    for( stockIdx in 1:nP )
    {
      stockID <- stockIDs[stockIdx]
      subComm <-  commData %>%
                  filter( species == specID, stockName == stockID )

      obsYrs    <- subComm$year
      obsYrIdx  <- which(obsYrs %in% years)
      datYrs    <- obsYrs[obsYrs %in% years]
      if( collapseComm )
      {
        C_spft[specID, stockID, 1, as.character(datYrs) ] <- subComm$catch[obsYrIdx]
        D_spft[specID, stockID, 1, as.character(datYrs) ] <- subComm$discWt[obsYrIdx]
      } else {
        histObsIdx  <- intersect( which( obsYrs < modernYear ), obsYrIdx)
        histYrs     <- datYrs[ datYrs < modernYear ]
        modObsIdx   <- intersect( which( obsYrs >= modernYear ), obsYrIdx)
        modYrs      <- datYrs[ datYrs >= modernYear ]
        # historic fleet
        C_spft[specID, stockID, 1, as.character(histYrs) ] <- subComm$catch[histObsIdx ]
        D_spft[specID, stockID, 1, as.character(histYrs) ] <- subComm$discWt[histObsIdx ]

        # modern fleet
        C_spft[specID, stockID, 2, as.character(modYrs) ] <- subComm$catch[modObsIdx ]
        D_spft[specID, stockID, 2, as.character(modYrs) ] <- subComm$discWt[modObsIdx ]        
      }

      # Now do the survey catch
      specStockSurvData <- specSurvData %>%
                           filter( stockName == stockID )

      surveyGears <- unique(specStockSurvData$surveyName)

      for( gIdx in 1:length(surveyGears))
      {
        surveyID        <- surveyGears[gIdx]
        subSurvData_spg <-  specStockSurvData %>%
                            filter( surveyName == surveyID )

        # Get years
        obsYrs    <- subSurvData_spg$year
        obsYrIdx  <- which(obsYrs %in% years)
        datYrs    <- obsYrs[ obsYrs %in% years ]

        C_spft[ specID, stockID, surveyID, as.character(datYrs) ] <- subSurvData_spg$catch[obsYrIdx]
      }

    }
  } 

  # return arrays
  outList <- list(  C_spft = C_spft,
                    D_spft = D_spft)

  outList
} # END makeCatchDiscArrays()

# makeALFreq()
# Makes the time-averaged Age-length frequency
# array for the hierSCAL integrated vonB growth model
makeALFreq <- function( ALFreqList = ALfreq,
                        years = 1954:2018,
                        gears = gearNames,
                        maxA = 35, maxL = 80 )
{
  nS <- length(ALFreqList)
  specIDs <- names( ALFreqList )

  ALfreqList <- vector(mode = "list", length = nS)
  
  nGears <- length(gears)
  nYears <- length(years)

  for( sIdx in 1:nS)
  {
    ALfreqList[[sIdx]] <- ALFreqList[[sIdx]]$ALfreq
    maxA <- max(maxA, dim(ALfreqList[[sIdx]])[2])
    maxL <- max(maxL, dim(ALfreqList[[sIdx]])[3])
  }


  nP <- dim(ALfreqList[[1]])[1]
  stockIDs  <- dimnames(ALfreqList[[1]])[[1]]
  sexIDs    <-dimnames(ALfreqList[[1]])[[6]]

  years <- as.character(years)

  ALK_spalftx <- array( 0, dim = c( nS, nP, maxA, maxL, nGears, nYears, 2),
                            dimnames = list(  species = specIDs,
                                              stocks = stockIDs,
                                              ages = 1:maxA,
                                              lengths = 1:maxL,
                                              fleet = gears,
                                              years = years,
                                              sex = sexIDs ) )

  for( sIdx in 1:nS )
  {
    specID <- specIDs[sIdx]
    specAL <- ALfreqList[[sIdx]]
    specAL[specAL < 0] <- NA
    specA <- dim(specAL)[2]
    specL <- dim(specAL)[3]
    ALK_spalftx[specID, stockIDs, 1:specA, 1:specL, gears, years, ] <- specAL[stockIDs,1:specA,1:specL, gears, years, ]
  }

  return(ALK_spalftx)
} # END makeALFreq()

# makeCompsArray()
# Function to transform age or length comps into 
# multidimensional array for feeding to hierSCAL
makeCompsArray <- function( compList = ageComps,
                            plusGroups = plusA_s,
                            minX  = minA_s,
                            collapseComm = FALSE,
                            fleetIDs = fleetIDs,
                            combineSex = TRUE,
                            years = fYear:lYear,
                            xName = "ages",
                            minSampSize = 100 )
{
  # Get species names
  specIDs   <- names(compList)
  # Array dimensions
  nX      <- max(plusGroups)
  nS      <- length(compList)
  nP      <- dim(compList[[1]])[1]
  nF      <- length(fleetIDs)
  nT      <- length(years)

  stockIDs <- dimnames(compList[[1]])[[1]]
  # Create dimension names
  dimNames <- list( 1:nX, specIDs, stockIDs, fleetIDs, years )
  names(dimNames) <- c(xName, "species", "stock", "fleet", "years")

  # Initialise array
  comps_xspft <- array( 0, dim = c(nX, nS, nP, nF, nT ),
                            dimnames = dimNames )

  for( specIdx in 1:nS )
  {
    # Get species specific info
    specID <- specIDs[specIdx]
    specX <- plusGroups[specID]
    specMin <- minX[specID]

    specComps <- compList[[specID]]
    obsX      <- dim(specComps)[4]


    # Now loop over fleetIDs
    for( fleetIdx in 1:nF )
    { 
      fleetID <- fleetIDs[fleetIdx]
      # Need to aggregate plus group
      for( stockIdx in 1:nP )
      {
        stockID <- stockIDs[stockIdx]
        for( tIdx in 1:nT)
        {
          yearLab <- as.character(years[tIdx])
          sumComps <- sum(specComps[stockID,fleetID,yearLab,1:(specX-1),1],na.rm = T)
          if(sumComps <= minSampSize )
            comps_xspft[ , specID, stockID ,fleetID, yearLab ] <- -1
          else {
            comps_xspft[specMin:(specX-1), specID, stockID ,fleetID,yearLab] <- specComps[stockID,fleetID,yearLab,specMin:(specX-1),1]
            comps_xspft[specX, specID, stockID ,fleetID,yearLab] <- sum(specComps[stockID,fleetID,yearLab,specX:obsX,1],na.rm = T)
            if( specX < nX )
              comps_xspft[(specX+1):nX, specID, stockID ,fleetID,yearLab] <- 0
          }
          
        }
      }
    }
  }
  comps_xspft[is.na(comps_xspft)] <- -1
  return(comps_xspft)
} # END makeCompsArray()


# Plot catch bars with discarding stacked on top
# Need to do some work to improve these - maybe
# we just do lines instead
plotCatch <- function(  data = catchData,
                        stocks = stocksCommBio,
                        speciesCodes = specCodes,
                        save = FALSE,
                        saveDir = "Outputs",
                        years = fYear:lYear )
{
  # Append stock name
  data <- data %>%
          filter( FISHERY_SECTOR == "GROUNDFISH TRAWL" ) %>%
          mutate( stockName = sapply( X = MAJOR_STAT_AREA_CODE,
                                      FUN = appendName,
                                      stocks ),
                  species = sapply( X = SPECIES_CODE,
                                    FUN = appendName,
                                    speciesCodes ) ) %>%
          filter( stockName %in% names(stocks) ) %>%
          dplyr::select(  year          = YEAR,
                          catch         = LANDED_KG,
                          discardWt     = DISCARDED_KG,
                          discardNum  = DISCARDED_PCS,
                          species, stockName ) %>%
          group_by( species, stockName, year ) %>%
          summarise(  catch   = sum(catch)/1e6,
                      discWt  = sum(discardWt)/1e6,
                      discNum = sum(discardNum) ) %>%
          ungroup()

  # Count species and stocks
  nSpec   <- length( unique(data$species) )
  nStocks <- length( stocks )

  vertlines <- seq(1960,2020, by = 10)

  # Set up save location
  if(save)
  {
    graphics.off()
    saveDir <- file.path(getwd(),saveDir)
    if(! dir.exists(saveDir) )
      dir.create(saveDir)
    fileName <- "catchAndDiscards.png"

    savePath <- file.path(saveDir, fileName )
    png(  savePath, width = 11, height = 8.5, units = "in",
          res = 300)
  }

  # Now set up plotting area,
  par(  mfcol = c(nStocks, nSpec ), 
        mar = c(0,2,0,1),
        oma = c(4,4,3,3) )

  # loop over species, then stocks
  # Loop over species, pull data and fits
  # for that species
  for( specIdx in 1:length(speciesCodes))
  {
    specName <- names(speciesCodes)[specIdx]
    specData <- data %>%
                filter( species == specName )

    
    # Then loop over stocks
    for( stockIdx in 1:length(stocks))
    {
 
      stockID       <-  names(stocks)[stockIdx]
      subData       <-  specData %>%
                        filter( stockName == stockID )

      maxCatch      <- 1.05 * max( subData$catch + subData$discWt )

      # Plot window
      plot( x = range(years),
            y = c(0,maxCatch),
            type = "n", axes = F, xlab = "", ylab = "",
            yaxs = "i")
        mfg <- par("mfg")
        axis(side = 2, las = 1)
        if( mfg[1] == mfg[3] )
          axis( side = 1 )
        if( mfg[1] == 1 )
          mtext( side = 3, text = specName, line = .5, font = 2 )
        if( mfg[2] == mfg[4] )
          mtext( side = 4, text = stockID, line = 2 )
        box()
        # stacked bar plots for catch and discards
        rect( xleft = subData$year - 0.5,
              xright = subData$year + 0.5,
              ybottom = 0,
              ytop = subData$catch,
              col = "grey60", lwd = .8 )
        rect( xleft = subData$year - 0.5,
              xright = subData$year + 0.5,
              ybottom = subData$catch,
              ytop = subData$catch + subData$discWt,
              col = "red", lwd = .8 )
        abline( v = vertlines, lwd = .8, lty = 3,
                col = "grey80")
    }
  }
  # labels and legend
  mtext(  side = 1, text = "Year", outer = TRUE, line = 2.5 )
  mtext(  side = 2, text = "Catch and Discards (kt)", 
          outer = TRUE, line = 2 )
  

  if( save )
  {
    dev.off()
    cat( "Catch and discard plots saved to ", savePath, "\n", sep = "")
  } 
}

# appendCommFleetName
# Helper function to append historical/modern
# commercial fleet ID
appendCommFleetName <- function(year)
{
  if( year < 1996 )
    fleet <- "comm.hist"
  else fleet <- "comm.mod"

  fleet
}

# makeALFreq_FleetYear()
# Conducts the length at age analysis for a given set of data,
# splitting over stocks and estimating vonB functions
# for each stock area, species, and sex. Spits out a length
# at age frequency array (species,stock,sex)
# inputs: data = output of readBioData()
# outputs:  vonB = list of vonB fits, 
#           data = input data  
#           ALfreq = array of age-length freq, indexed by 
#                     stock and species
makeALFreq_FleetYear <- function( data = bioData$Dover,
                                  stocksComm = stocksCommBio,
                                  stocksSurv = stocksSurvey,
                                  survIDs = surveyIDs,
                                  years = 1954:2018  )
{
  # Get survey and commercial data
  survData <- data$survey %>%
              mutate(fleetID = sapply( X = SURVEY_SERIES_ID,
                                      FUN = appendName,
                                      survIDs) )
  commData <- data$comm %>%
              mutate( fleetID = sapply( X = YEAR, FUN = appendCommFleetName) )

  # We want to make an array to hold age comps,
  # so we need the largest observed age
  maxAge <- max(survData$AGE, commData$AGE, na.rm = T)
  maxLen <- round(max(survData$LENGTH_MM, commData$LENGTH_MM, na.rm = T)/10)

  # Make a vector of gear names - split
  # commercial data into modern and historic
  gearNames <- c(names(survIDs),"comm.hist","comm.mod")
  stockNames <- names(stocksComm)

  # Count dimensions
  nGears  <- length(gearNames)
  nYears  <- length(years)
  nStocks <- length(stocksSurv)


  # We want to filter down to cases where there are age/length
  # observations
  survData <- survData %>%
              filter( YEAR %in% years) %>%
              dplyr::select(  age = AGE, length = LENGTH_MM,
                              stockName, sex = SEX, fleetID,
                              year = YEAR ) %>%
              mutate( length = length/10 ) %>%
              filter( !is.na(age),
                      !is.na(length),
                      length > 0 )

  commData <- commData %>%
              filter( YEAR %in% years) %>%
              dplyr::select(  age = AGE, length = LENGTH_MM,
                              stockName, sex = SEX, fleetID,
                              year = YEAR ) %>%
              filter( !is.na(age),
                      !is.na(length),
                      length > 0 ) %>%
              mutate( length = length/10 )

  data$survey <- survData
  data$comm   <- commData 

  # combine total data for likelihood
  combData <- rbind(survData,commData)

  bioDataBoys   <- combData %>% filter( sex == 1 )
  bioDataGirls  <- combData %>% filter( sex == 2 )

  # initialise age-length freq
  # First, get array dimensions
  lengths <- unique(combData$length)
  ages    <- unique(combData$age)

  allAges     <- 1:max(ages)
  allLengths  <- 1:max(lengths)

  # initialise array
  ALfreq <- array( 0,  dim = c(  nStocks,
                                  length(allAges),
                                  length(allLengths),
                                  nGears,
                                  nYears,
                                  2 ),
                        dimnames = list(  stock = stockNames, 
                                          age = allAges, 
                                          length = allLengths, 
                                          fleet = gearNames,
                                          year = years,
                                          sex = c("boys","girls") ) )


  # Now populate the array - round length to nearest cm
  # might need to aggregate further
  frqData <-  combData %>%
              mutate( length = round(length) ) %>%
              group_by( stockName, year, age, length, sex, fleetID ) %>%
              filter( sex %in% c(1,2), !is.na(fleetID) ) %>%
              summarise( nObs = n() ) %>%
              ungroup()


  # Now loop over stocks
  for( stockIdx in 1:length(stockNames) )
  {
    for( gearIdx in 1:length(gearNames))
    {
      gearName <- gearNames[gearIdx]
      gearData <- frqData %>% 
                  filter( stockName == stockNames[stockIdx],
                          fleetID == gearNames[gearIdx] )

      gearAges <- unique( gearData$age )
      gearLengths <- unique( gearData$length )

      for( a in gearAges )
        for( l in gearLengths )
        {

          boyObs <- gearData %>% 
                    filter( sex == 1, 
                            age == a, 
                            length == l )
          girlObs <-  gearData %>% 
                      filter( sex == 2, 
                              age == a, 
                              length == l )

          if(nrow(boyObs) > 0)
          {
            ALfreq[ stockNames[stockIdx], a, l, gearName, as.character(boyObs$year), "boys" ] <- boyObs$nObs
          }
          if(nrow(girlObs) > 0)
          {
            ALfreq[ stockNames[stockIdx], a, l, gearName, as.character(girlObs$year), "girls" ] <- girlObs$nObs
          }

        }
    }
  }

  outList <- list(  data = data,
                    ALfreq = ALfreq )

  outList
} # END makeALFreq_FleetYear()




# makeAgeComps()
# Takes output of readBioData() for a species
# and generates an array of age observation
# frequencies.
makeAgeComps <- function( data = bioData$Dover,
                          stocksComm = stocksCommBio,
                          stocksSurv = stocksSurvey,
                          survIDs = surveyIDs,
                          years = 1954:2018  )
{
  # Pull survey and commercial data
  survData <- data$survey %>%
              mutate(survID = sapply( X = SURVEY_SERIES_ID,
                                      FUN = appendName,
                                      survIDs) )
  commData <- data$comm

  # We want to make an array to hold age comps,
  # so we need the largest observed age
  maxAge <- max(survData$AGE, commData$AGE, na.rm = T)

  # Make a vector of gear names - split
  # commercial data into modern and historic
  gearNames <- c(names(survIDs),"comm.hist","comm.mod")
  stockNames <- names(stocksComm)

  # Count dimensions
  nGears  <- length(gearNames)
  nYears  <- length(years)
  nStocks <- length(stocksSurv)

  # initialise array
  compFreq <- array(0,  dim = c(nStocks,nGears,nYears,maxAge,3),
                        dimnames = list(  stock = stockNames,
                                          fleets = gearNames,
                                          years = as.character(years),
                                          age = as.character(1:maxAge),
                                          sex = c("all", "boys", "girls" ) ) ) 
  for( stockIdx in 1:nStocks )
  {
    stockID <- stockNames[stockIdx]
    for( gearIdx in 1:nGears )
    {
      gearID <- gearNames[gearIdx]
      if( gearID %in% c("comm.hist","comm.mod") )
      {
        # Filter down to this gear 
        gearData <- commData %>%
                    filter( !is.na(AGE), 
                            stockName == stockID ) %>%
                    dplyr::select(  age = AGE,
                                    year = YEAR,
                                    sex = SEX ) 
        # Skip if zero observations
        if( nrow(gearData) == 0 )
          next

        # Summarise all data
        allData <-  gearData %>%
                    group_by( year, age ) %>%
                    summarise(  nObs = n(),
                                meanAge = mean(age),
                                sdAge = sd(age) ) %>%
                    ungroup()
        # Boys
        boyData <-  gearData %>%
                    filter( sex == 1 ) %>%
                    group_by( year, age ) %>%
                    summarise(  nObs = n(),
                                meanAge = mean(age),
                                sdAge = sd(age) ) %>%
                    ungroup()
        # Girls
        girlData <- gearData %>%
                    filter( sex == 2 ) %>%
                    group_by( year, age ) %>%
                    summarise(  nObs = n(),
                                meanAge = mean(age),
                                sdAge = sd(age) ) %>%
                    ungroup()

        # Pull those years with observations
        obsYrs <- unique(gearData$year)

        # 
        if( gearID == "comm.hist")
          obsYrs <- obsYrs[ obsYrs %in% 1954:1995 ]
        if( gearID == "comm.mod")
          obsYrs <- obsYrs[ obsYrs %in% 1996:2018 ]

        # Loop over years with observations
        for( yr in obsYrs )
        {
          yrChar <- as.character(yr)
          # Subset to that year's data
          yrData.all <- allData %>% filter( year == yr )
          # Fill all fish slice
          if( nrow(yrData.all) > 0)
          {
            # Replace NAs with 0s, as we have some observations
            compFreq[stockID, gearID, yrChar, ,"all" ] <- 0  
            # Get vector of ages with positive observations
            obsAges <- as.character(yrData.all$age)
            # Save into compFreq
            compFreq[stockID, gearID, yrChar, obsAges, "all" ] <- yrData.all$nObs            
          }
          
          # Subset to that year's data
          yrData.boy <- boyData %>% filter( year == yr )
          # Fill boy fish slice
          if( nrow(yrData.boy) > 0)
          {
            # Replace NAs with 0s, as we have some observations
            compFreq[stockID, gearID, yrChar, ,"boys" ] <- 0  
            # Get vector of ages with positive observations
            obsAges <- as.character(yrData.boy$age)
            # Save into compFreq
            compFreq[stockID, gearID, yrChar, obsAges, "boys" ]  <- yrData.boy$nObs
          }

          # Subset to that year's data
          yrData.girl <- girlData %>% filter( year == yr )
          # Fill girl fish slice
          if( nrow(yrData.girl) > 0)
          {
            # Replace NAs with 0s, as we have some observations
            compFreq[stockID, gearID, yrChar, ,"girls" ] <- 0  
            # Get vector of ages with positive observations
            obsAges <- as.character(yrData.girl$age)
            # Save into compFreq
            compFreq[stockID, gearID, yrChar, obsAges, "girls" ]  <- yrData.girl$nObs
          }   

        }
        # Go to next gearID
        next
      }

      gearData <- survData %>%
                  filter( !is.na(AGE), 
                          stockName == stockID ) %>%
                  dplyr::select(  age = AGE,
                                  year = YEAR,
                                  survID,
                                  sex = SEX ) %>%
                  filter( survID == gearID ) 
      # Skip if zero observations
      if( nrow(gearData) == 0 )
          next

      # Take mean and SD of ages for this gear
      gearMeanSD <- gearData %>%
                    summarise(  meanAge = mean(age),
                                sdAge = sd(age) )


      # Filter mean and sd by sex
      gearMeanSD.sex <- gearData %>%
                        filter( sex %in% c(1,2)) %>%
                        group_by( sex ) %>%
                        summarise(  meanAge = mean(age),
                                    sdAge = sd(age) )
      

      # Summarise combined sexes
      allData <-  gearData %>%
                  group_by( year, age ) %>%
                  summarise(  nObs = n(),
                              meanAge = mean(age),
                              sdAge = sd(age) ) %>%
                  ungroup()

      # Boys
      boyData <-  gearData %>%
                  filter( sex == 1 ) %>%
                  group_by( year, age ) %>%
                  summarise(  nObs = n(),
                              meanAge = mean(age),
                              sdAge = sd(age) ) %>%
                  ungroup()
      # Girls
      girlData <- gearData %>%
                  filter( sex == 2 ) %>%
                  group_by( year, age ) %>%
                  summarise(  nObs = n(),
                              meanAge = mean(age),
                              sdAge = sd(age) ) %>%
                  ungroup()

      # Pull those years with observations
      obsYrs <- unique(gearData$year)

      # Loop over years with observations
      for( yr in obsYrs )
      {
        yrChar <- as.character(yr)
        # Subset to that year's data
        yrData.all <- allData %>% filter( year == yr )
        # Fill all fish slice
        if( nrow(yrData.all) > 0)
        {
          # Replace NAs with 0s, as we have some observations
          compFreq[stockID, gearID, yrChar, ,"all" ] <- 0  
          # Get vector of ages with positive observations
          obsAges <- as.character(yrData.all$age)
          # Save into compFreq
          compFreq[stockID, gearID, yrChar, obsAges, "all"] <- yrData.all$nObs
        }
        
        # Subset to that year's data
        yrData.boy <- boyData %>% filter( year == yr )
        # Fill boy fish slice
        if( nrow(yrData.boy) > 0)
        {
          # Replace NAs with 0s, as we have some observations
          compFreq[stockID, gearID, yrChar, ,"boys" ] <- 0  
          # Get vector of ages with positive observations
          obsAges <- as.character(yrData.boy$age)
          # Save into compFreq
          compFreq[stockID, gearID, yrChar, obsAges, "boys" ]<- yrData.boy$nObs
        }

        # Subset to that year's data
        yrData.girl <- girlData %>% filter( year == yr )
        # Fill girl fish slice
        if( nrow(yrData.girl) > 0)
        {
          # Replace NAs with 0s, as we have some observations
          compFreq[stockID, gearID, yrChar, ,"girls" ] <- 0  
          # Get vector of ages with positive observations
          obsAges <- as.character(yrData.girl$age)
          # Save into compFreq
          compFreq[stockID, gearID, yrChar, obsAges, "girls" ]<- yrData.girl$nObs
        }   

      }

    }
  }

  # Return age frequencies
  return(compFreq)
} # END makeAgeComps()


# plotAgeComps()
# Takes the output from makeAgeComps() and plots
# age frequency distributions for each stock and gear,
# and all years that gear has observations of that stock.
# Will skip gears that don't fish a given stock.
plotComps <- function(  comps = ageComps,
                        save = FALSE,
                        saveDir = "Outputs",
                        prefix = "age" )
{
  # OK, we want to loop over species, then stocks
  # Set up save location
  if(save)
  {
    graphics.off()
    saveDir <- file.path(getwd(),saveDir)
    if(! dir.exists(saveDir) )
      dir.create(saveDir)
  }
  # Get species names
  nSpec     <- length(comps)
  specNames <- names(comps)

  # Loop over species
  for( specIdx in 1:nSpec )
  {
    # Get species name
    specID <- specNames[specIdx]

    if( save )
    {
      # Create save directory if it doesn't exist
      specDir <- file.path(saveDir,specID)
      if( !dir.exists( specDir ) )
        dir.create( specDir )
    }

    # pull age freq array
    freq    <- comps[[specIdx]]$freq
    meanSD  <- comps[[specIdx]]$meanSD

    # Count stocks, gears, age classes
    nStocks <- dim(ageFrq)[1]
    nGears  <- dim(ageFrq)[2]
    years   <- as.integer(dimnames(ageFrq)[[3]])
    nA      <- dim(ageFrq)[4]

    # Get stock and gear IDs
    stockNames  <- dimnames(ageFrq)[[1]]
    gearNames   <- dimnames(ageFrq)[[2]]

    # Loop over stocks
    for( stockIdx in 1:nStocks )
    {
      stockID <- stockNames[stockIdx]
      if( save )
      {
        # Create save directory if it doesn't exist
        stockDir <- file.path(specDir,stockID)
        if( !dir.exists( stockDir ) )
          dir.create( stockDir )
      }

      # Here I want to refactor and create 2 functions:
      # 1. Bubble plots
      if( save )
      {
        fileName <- paste( prefix,"BubblesByGear.png", sep = "" )
        png(  file = file.path(stockDir,fileName),
              width = 8.5, height = 11, res = 300,
              units = "in" )
      }


      # Set up multipanel region
      par(mfrow = c(nGears,3), oma = c(4,4,4,2), mar = c(0,0,0,0) )

      # Now loop over gears
      for( gearIdx in 1:nGears )
      {
        # Get gear name
        gearName <- gearNames[gearIdx]

        # Pull subset
        subFrq <- freq[ stockIdx, gearIdx,,,,]

        
        # Plot boys
        .plotBubbles( z = subFrq[,,"boys","freq"])
        mfg <- par("mfg")
        if(mfg[1] == 1)
          mtext( side = 3, text = "Male", line = 1 )
        # Plot girls
        .plotBubbles( z = subFrq[,,"girls","freq"])
        mfg <- par("mfg")
        if(mfg[1] == 1)
          mtext( side = 3, text = "Female", line = 1 )
        # Plot combined
        .plotBubbles( z = subFrq[,,"all","freq"])
        mfg <- par("mfg")
        if(mfg[1] == 1)
          mtext(  side = 3, text = "Both + Unsexed", line = 1 )

        mtext( side = 4, text = gearName, font = 2, line = 1 )
        # Add title
        mtext(  side = 3, text = paste(specID, stockID, sep = " - " ),
                outer = T, line = 3, cex = 1, font = 2 )

      }
      # Turn save output off
      if( save )
        dev.off()

      # 2. age freqency plots
      # Now loop over gears
      for( gearIdx in 1:nGears )
      {
        # Get gear name
        gearName <- gearNames[gearIdx]

        # Pull subset
        subFrq <- freq[ stockIdx, gearIdx,,,]

        if( sum(subFrq, na.rm = TRUE) == 0 )
          next

        # Now we want to plot the age frequency plots
        # but we can't plot the three groups (boys, girls, unsexed)
        # in one plot, so we need to plot each in a separate device
        # Plot males
        if(save)
        {
          fileName <- paste( prefix, "Freq_", gearName,"_Males.png", sep = "")
          png(  file = file.path(stockDir,fileName),
                width = 8.5, height = 11, res = 300,
                units = "in" )
        }
        
        .plotFreqBars( z = subFrq[,,"boys","freq"] )
        mtext(  side = 3, outer = TRUE, font = 2,
                text = paste(specID, " Males - ", stockID ," - ", gearName, sep = "") )

        if( save )
          dev.off()

        # Plot females
        if(save)
        {
          fileName <- paste( prefix, "Freq_", gearName,"_Females.png", sep = "")
          png(  file = file.path(stockDir,fileName),
                width = 8.5, height = 11, res = 300,
                units = "in" )
        }
        
        .plotFreqBars( z = subFrq[,,"girls","freq"] )
        mtext(  side = 3, outer = TRUE, font = 2,
                text = paste(specID, " Females - ", stockID ," - ", gearName, sep = "") )

        if( save )
          dev.off()

        # Plot combined
        if(save)
        {
          fileName <- paste( prefix, "Freq_", gearName,"_Combined.png", sep = "")
          png(  file = file.path(stockDir,fileName),
                width = 8.5, height = 11, res = 300,
                units = "in" )
        }
        
        .plotFreqBars( z = subFrq[,,"all","freq"] )
        mtext(  side = 3, outer = TRUE, font = 2,
                text = paste(specID, " Combined - ", stockID ," - ", gearName, sep = "") )

        if( save )
          dev.off()

        
      }
    }
  } 

} # END plotAgeComps()

# Hidden function to plot age bubbles - shamelessly stolen
# from sableOpMod.R (SPC and ARK), but simplified for 
# a specific application and our arrays with named dimensions.
# inputs:   z = years x ages array of age observation frequencies
.plotBubbles <- function( z,
                          years = NULL,
                          ages = NULL,
                          minAge = 1,
                          initYear = 1954,
                          hide0 = TRUE,
                          lwd = 1,
                          pwr = .5,
                          size = .1,
                          forceAxes = FALSE )
{
  # Get dimensions of z
  dz <- dim( z )

  # Get number of years
  nYears  <- dz[ 1 ]
  nAges   <- dz[ 2 ]

  # make colours
  clrs <- c( "grey40", "salmon", "steelblue" )

  # Create age and years vectors
  if( is.null(years) )
    years <- initYear + 1:nYears - 1

  if( is.null(ages) )
    ages <- minAge + 1:nAges - 1

  # Create vectors of x and y arguments
  # for a symbols call
  xArg <- rep( years, nAges )
  yArg <- rep( ages, each = nYears )

  # Sweep out sum of age observations
  # to convert to proportions
  zSum  <- apply( X = z, FUN = sum, MARGIN = 1, na.rm = TRUE )
  zz    <- sweep( x = z, MARGIN = 1, STATS = zSum, FUN = "/")


  # Separate positive, negative and zero values
  # for different plotting colours
  zNA <- is.na(zz) | is.nan(zz) | is.infinite(zz)
  zz[zNA] <- 0
  z0 <- sign(zz) * abs(zz)^abs(pwr)
  z1 <- z3 <- z0
  z1[z0 <= 0] <- NA
  z3[z0 < 0 | z0 > 0] <- NA
  z2 <- -z0
  z2[z0 >= 0] <- NA
  za <- max(z0, na.rm = TRUE)
  zb <- min(z0, na.rm = TRUE)
  zM <- max(abs(z0))
  sz1 <- max(za * size/zM, 0.001)
  sz2 <- max(-zb * size/zM, 0.001)

  plot( x = range(years), y = range(ages,nAges+1),
        type = "n", xlab = "", ylab = "", axes = F )
    mfg <- par("mfg")
    if( mfg[1] == mfg[3] | forceAxes )
      axis( side = 1 )
    if( mfg[2] == 1 | forceAxes )
      axis( side = 2, las = 1 )
    box()

    if( !hide0 && !all(is.na(z3)) ) 
    {
        PBSmodelling::evalCall(symbols, argu = list(x = xArg, y = yArg, circles = as.vector(z3), 
            inches = 0.001, fg = clrs[3], lwd = lwd, add = TRUE), 
            checkpar = TRUE)
    }
    if( !all(is.na(z2)) ) 
    {
        PBSmodelling::evalCall(symbols, argu = list(x = xArg, y = yArg, circles = as.vector(z2), 
            inches = sz2, fg = clrs[2], lwd = lwd, add = TRUE), 
            checkpar = TRUE)
    }
    if( !all(is.na(z1)) ) 
    {
        PBSmodelling::evalCall(symbols, argu = list(x = xArg, y = yArg, circles = as.vector(z1), 
            inches = sz1, fg = clrs[1], lwd = lwd, add = TRUE), 
            checkpar = TRUE)
    }
    if( any(zSum > 0) )
      text( x = years[zSum > 0], y = nAges + 1, labels = zSum[zSum > 0],
            srt = 45, cex = .5 )

} # END .plotBubbles()

# Hidden function to plot age frequencies. Shamelessly
# copied from ARK and SPC functions in SableOpMod.R,
# with simplifications for our purpose.
.plotFreqBars <- function(  z,
                            years = NULL,
                            ages = NULL,
                            minAge = 1,
                            initYear = 1954,
                            avg = FALSE,
                            delta = .4 )
{
  # Get dimensions of z
  dz <- dim( z )

  # Get number of years
  nYears  <- dz[ 1 ]
  nAges   <- dz[ 2 ]

  # Create age and years vectors
  if( is.null(years) )
    years <- initYear + 1:nYears - 1

  if( is.null(ages) )
    ages <- minAge + 1:nAges - 1

  # Take sum of rows to get the total
  zSum  <- apply( X = z, FUN = sum, MARGIN = 1, na.rm = TRUE )
  # Weed out negative sums for years of missing data
  posYrIdx  <- which(zSum > 0)
  posYrs    <- years[ posYrIdx ]
  zz        <- z[posYrIdx,] 
  
  # Sweep out sum to make proportions
  zz        <- sweep( x = zz, FUN = "/", STATS = zSum[posYrIdx], MARGIN = 1 )

  # Count years with positive observations
  nPanels  <- length(posYrIdx)

  myMar <- c( 2, 2, 1, 1 )
  myOma <- c( 4, 3, 3, 1 )

  # If not averaging, set up multi-panel
  if( !avg )
  {
    if( nPanels == 1 )
      par( oma = myOma, mar = myMar, mfcol = c( 1, 1 ) )
    if( nPanels <= 9 )
      par( oma=myOma, mar=myMar, mfcol=c(3,3) )
    else if ( nPanels > 9 & nPanels <= 12 )
      par( oma=myOma, mar=myMar, mfcol=c(4,3) )
    else if ( nPanels > 12 & nPanels <= 16 )
      par( oma=myOma, mar=myMar, mfcol=c(4,4) )
    else if ( nPanels > 16 & nPanels <= 20 )
      par( oma=myOma, mar=myMar, mfcol=c(5,4) )
    else if ( nPanels > 20 & nPanels <=24 )
      par( oma=myOma, mar=myMar, mfcol=c(6,4) )
    else
      par( oma=myOma, mar=myMar, mfcol=c(6,5) )
  }

  if( avg )
  {
    zz <- matrix( apply( X = zz, FUN = mean, 
                         MARGIN = 2, na.rm = T), 
                  nrow = 1 )
  }

  # Get xLim and yLim
  xLim <- range(ages)


  for( i in 1:nrow(zz) )
  {
    # Make labels for each plot
    if( avg ) 
    {
      yearLab <- "Averaged"
      numObs  <- sum( zSum )
    } else {
      yearLab <- posYrs[ i ]
      numObs  <- zSum[ posYrIdx[i] ]
    }
    # Plot away
    yLim <- c(0, 1.2*max(zz[i,], na.rm = T) )
    plot( xLim, yLim, type = "n",
          axes = F, xlab = "", ylab = "", yaxs = "i" )
      axis( side = 1 )
      axis( side = 2, las = 1 )
      box()
      text( x = 0.8*nAges, y = 0.95*yLim[2], label = yearLab, cex = .8 )
      text( x = 0.8*nAges, y = 0.9*yLim[2], label = paste("N =", numObs ), cex = .8 )
      rect( xleft = ages - delta, xright = ages + delta,
            ybottom = 0, ytop = zz[i,],
            col = "grey80" )

  }
  mtext( side = 1, line = 1, outer = TRUE, 
          text = "Age Class" )
  mtext( side = 2, line = 2, outer = TRUE, 
          text = "Proportion-at-age" )
} # END .plotFreqBars()


# makeLenComps()
# Takes output of readBioData() for a species
# and generates an array of length observation
# frequencies.
makeLenComps <- function( data = bioData$Dover,
                          stocksComm = stocksCommBio,
                          stocksSurv = stocksSurvey,
                          survIDs = surveyIDs,
                          years = 1954:2018  )
{
  # Pull survey and commercial data
  survData <- data$survey %>%
              filter(is.na(AGE)) %>%
              mutate( survID = sapply(  X = SURVEY_SERIES_ID,
                                        FUN = appendName,
                                        survIDs),
                      length = round(LENGTH_MM/10) )

  commData <- data$comm %>%
              filter(is.na(AGE)) %>%
              mutate( length = round(LENGTH_MM/10) )

  commTrips <-  commData %>%
                group_by(stockName,YEAR,SAMPLE_ID) %>%
                summarise(nSamples = n()) %>%
                ungroup() %>%
                filter( nSamples >= 4 ) 



  # We want to make an array to hold age comps,
  # so we need the largest observed age
  maxLen <- max(survData$length, commData$length, na.rm = T)

  # Make a vector of gear names - split
  # commercial data into modern and historic
  gearNames <- c(names(survIDs),"comm.hist","comm.mod")
  stockNames <- names(stocksComm)

  # Count dimensions
  nGears  <- length(gearNames)
  nYears  <- length(years)
  nStocks <- length(stocksSurv)

  # initialise array
  compFreq <- array(NA, dim = c(nStocks,nGears,nYears,maxLen,3),
                        dimnames = list(  stock = stockNames,
                                          fleets = gearNames,
                                          years = as.character(years),
                                          length = as.character(1:maxLen),
                                          sex = c("all", "boys", "girls" ) ) ) 
  
  


  for( stockIdx in 1:nStocks )
  {
    stockID <- stockNames[stockIdx]
    for( gearIdx in 1:nGears )
    {
      gearID <- gearNames[gearIdx]
      if( gearID %in% c("comm.hist","comm.mod") )
      {
        stockTrips <- commTrips %>%
                      filter(stockName == stockID)

        # Filter down to this gear 
        gearData <- commData %>%
                    filter( !is.na(length), 
                            stockName == stockID,
                            YEAR %in% stockTrips$YEAR ) %>%
                    dplyr::select(  length,
                                    year = YEAR,
                                    sex = SEX ) 

        if( nrow(gearData) == 0 )
          next

        # Summarise all data
        allData <-  gearData %>%
                    group_by( year, length ) %>%
                    summarise(  nObs = n(),
                                meanLen = mean(length),
                                sdLen = sd(length) ) %>%
                    ungroup()
        # Boys
        boyData <-  gearData %>%
                    filter( sex == 1 ) %>%
                    group_by( year, length ) %>%
                    summarise(  nObs = n(),
                                meanLen = mean(length),
                                sdLen = sd(length)  ) %>%
                    ungroup()
        # Girls
        girlData <- gearData %>%
                    filter( sex == 2 ) %>%
                    group_by( year, length ) %>%
                    summarise(  nObs = n(),
                                meanLen = mean(length),
                                sdLen = sd(length)  ) %>%
                    ungroup()

        # Pull those years with observations
        obsYrs <- unique(gearData$year)

        # 
        if( gearID == "comm.hist")
          obsYrs <- obsYrs[ obsYrs %in% 1954:1995 ]
        if( gearID == "comm.mod")
          obsYrs <- obsYrs[ obsYrs %in% 1996:2018 ]

        # Loop over years with observations
        for( yr in obsYrs )
        {
          yrChar <- as.character(yr)
          # Subset to that year's data
          yrData.all <- allData %>% filter( year == yr )
          # Fill all fish slice
          if( nrow(yrData.all) > 0)
          {
            # Replace NAs with 0s, as we have some observations
            compFreq[stockID, gearID, yrChar, ,"all" ] <- 0  
            # Get vector of ages with positive observations
            obsLens <- as.character(yrData.all$length)
            # Save into compFreq
            compFreq[stockID, gearID, yrChar, obsLens, "all" ]<- yrData.all$nObs
          }
          
          # Subset to that year's data
          yrData.boy <- boyData %>% filter( year == yr )
          # Fill boy fish slice
          if( nrow(yrData.boy) > 0)
          {
            # Replace NAs with 0s, as we have some observations
            compFreq[stockID, gearID, yrChar, ,"boys" ] <- 0  
            # Get vector of ages with positive observations
            obsLens <- as.character(yrData.boy$length)
            # Save into compFreq
            compFreq[stockID, gearID, yrChar, obsLens, "boys" ]<- yrData.boy$nObs
          }

          # Subset to that year's data
          yrData.girl <- girlData %>% filter( year == yr )
          # Fill girl fish slice
          if( nrow(yrData.girl) > 0)
          {
            # Replace NAs with 0s, as we have some observations
            compFreq[stockID, gearID, yrChar, ,"girls" ] <- 0  
            # Get vector of ages with positive observations
            obsLens <- as.character(yrData.girl$length)
            # Save into compFreq
            compFreq[stockID, gearID, yrChar, obsLens, "girls" ]<- yrData.girl$nObs
          }   

        }
        # Go to next gearID
        next
      }

      gearData <- survData %>%
                  filter( stockName == stockID,
                          !is.na(length) ) %>%
                  dplyr::select(  length,
                                  year = YEAR,
                                  survID,
                                  sex = SEX ) %>%
                  filter( survID == gearID ) 

      if( nrow(gearData) == 0 )
          next

      # Summarise combined sexes
      allData <-  gearData %>%
                  group_by( year, length ) %>%
                  summarise(  nObs = n(),
                              meanLen = mean(length),
                              sdLen = sd(length)  ) %>%
                  ungroup()

      # Boys
      boyData <-  gearData %>%
                  filter( sex == 1 ) %>%
                  group_by( year, length ) %>%
                  summarise(  nObs = n(),
                              meanLen = mean(length),
                              sdLen = sd(length)  ) %>%
                  ungroup()
      # Girls
      girlData <- gearData %>%
                  filter( sex == 2 ) %>%
                  group_by( year, length ) %>%
                  summarise(  nObs = n(),
                              meanLen = mean(length),
                              sdLen = sd(length)  ) %>%
                  ungroup()

      # Pull those years with observations
      obsYrs <- unique(gearData$year)

      # Loop over years with observations
      for( yr in obsYrs )
      {
        yrChar <- as.character(yr)
        # Subset to that year's data
        yrData.all <- allData %>% filter( year == yr )
        # Fill all fish slice
        if( nrow(yrData.all) > 0)
        {
          # Replace NAs with 0s, as we have some observations
          compFreq[stockID, gearID, yrChar, ,"all" ] <- 0  
          # Get vector of ages with positive observations
          obsLens <- as.character(yrData.all$length)
          # Save into compFreq
          compFreq[stockID, gearID, yrChar, obsLens, "all" ]<- yrData.all$nObs
        }
        
        # Subset to that year's data
        yrData.boy <- boyData %>% filter( year == yr )
        # Fill boy fish slice
        if( nrow(yrData.boy) > 0)
        {
          # Replace NAs with 0s, as we have some observations
          compFreq[stockID, gearID, yrChar, ,"boys" ] <- 0  
          # Get vector of ages with positive observations
          obsLens <- as.character(yrData.boy$length)
          # Save into compFreq
          compFreq[stockID, gearID, yrChar, obsLens, "boys" ]<- yrData.boy$nObs
          
        }

        # Subset to that year's data
        yrData.girl <- girlData %>% filter( year == yr )
        # Fill girl fish slice
        if( nrow(yrData.girl) > 0)
        {
          # Replace NAs with 0s, as we have some observations
          compFreq[stockID, gearID, yrChar, ,"girls" ] <- 0  
          # Get vector of ages with positive observations
          obsLens <- as.character(yrData.girl$length)
          # Save into compFreq
          compFreq[stockID, gearID, yrChar, obsLens, "girls" ]<- yrData.girl$nObs
          
        }   

      }

    }
  }

  # Return age frequencies
  return(compFreq)
} # END makeLenComps()

# Small function to switch maturity
# codes to boolean indicators of maturity
isMature <- function( mat, matThresh )
{
  mat[mat < matThresh]    <- 0
  mat[mat >= matThresh ]  <- 1

    
  mat
} # END isMature()


# makeMatOgives()
# Creates mat-age and mat-length ogives
# for a given set of data.
# Used as a part of makeSpecMat()
makeMatOgives <- function(  matData, 
                            maxPropMat = .99 )
{

  # First, let's make the age ogive
  ageProps  <-  matData %>%
                group_by( age ) %>%
                summarise(  propMat = mean(mat),
                            nObs = n() )
  
  # We want to filter out the ages that have all 100%
  # mature, as these bias the model

  # Get maximum age to be used in the logistic model
  subAgeProps <- ageProps %>% filter( propMat >= maxPropMat )
  maxAge      <- quantile(subAgeProps$age,prob = .5)
  # Now reduce data from the full set
  ageModelData <- matData %>%
                  filter( age <= maxAge )

  # Then fit the logistic model
  ageModel <- glm( mat ~ age, family = binomial(logit), data = ageModelData )
  a50 <- dose.p(ageModel,p=c(0.5))
  a95 <- dose.p(ageModel,p=c(0.95))

  # Collect in a list
  ageList <- list(  model   = ageModel,
                    props   = ageProps,
                    data    = ageModelData,
                    maxAge  = maxAge,
                    x50     = a50,
                    x95     = a95 )


  # Now do the same for length
  lenProps <- matData %>%
              group_by( length ) %>%
              summarise(  propMat = mean(mat),
                          nObs = n() )
  # We want to filter out the lengths that have all 100%
  # mature, as these bias the model
  # Get maximum length to be used in the logistic model
  subLenProps   <- lenProps %>% filter( propMat >= maxPropMat )
  maxLen        <- quantile(subLenProps$length,prob = .5)
  # reduce length model data
  lenModelData  <- matData %>% filter( length <= maxLen )

  # Then fit the logistic model
  lenModel <- glm( mat ~ length, family = binomial(logit), data = lenModelData )
  l50 <- dose.p(lenModel,p=c(0.5))
  l95 <- dose.p(lenModel,p=c(0.95))

  # Collect length model in a list
  lenList <- list(  model   = lenModel,
                    props   = lenProps,
                    data    = lenModelData,
                    maxLen  = maxLen,
                    x50     = l50,
                    x95     = l95 )

  # Return age and length models
  outList <- list( age = ageList, length = lenList )

  outList
} # END makeMatOgives()


# makeSpecMat()
# Takes output from readBioData() and creates maturity-at-length
# and maturity-at-age ogives for each stock and species.
makeSpecMat <- function(  data = bioData$Dover,
                          stocksComm = stocksCommBio,
                          stocksSurv = stocksSurvey,
                          survIDs = surveyIDs,
                          maxPropMat = .99 )
{
  # Pull survey and commercial data
  survData <- data$survey %>%
              dplyr::select(  age = AGE, mat = MAT,
                              length = LENGTH_MM,
                              stockName, sex = SEX ) %>%
              filter( !is.na(age), !is.na(mat),
                      !is.na(length) ) %>%
              mutate( length = round(length/10),
                      mat = isMature( mat, 3)  )

  commData <- data$comm %>%
              dplyr::select(  age = AGE, mat = MAT,
                              length = LENGTH_MM,
                              stockName, sex = SEX ) %>%
              filter( !is.na(age), !is.na(mat),
                      !is.na(length) ) %>%
              mutate( length = round(length/10),
                      mat = isMature( mat, 3) )

  # Combine data
  combData  <- rbind( survData, commData )

  stocks    <- names(stocksComm)
  nStocks   <- length(stocks)

  stockList <- vector(mode = "list", length = nStocks )
  names( stockList ) <- stocks

  # Make a coastwide model
  coastWideAll <- makeMatOgives(  matData = combData, 
                                  maxPropMat = maxPropMat )
  # Subset by sex
  cwBoys  <- combData %>% filter( sex == 1 )
  cwGirls <- combData %>% filter( sex == 2 )
  
  if( nrow( cwBoys ) > 0)
    coastWideBoys <- makeMatOgives( matData = cwBoys, 
                                    maxPropMat = maxPropMat )

  if( nrow( cwGirls ) > 0)
    coastWideGirls <- makeMatOgives( matData = cwGirls, 
                                    maxPropMat = maxPropMat )    

  cwList <- list( all   = coastWideAll,
                  boys  = coastWideBoys,
                  girls = coastWideGirls  )
  # Now loop over stocks
  for( stockIdx in 1:nStocks )
  {
    stockData <-  combData %>%
                  filter( stockName == stocks[stockIdx] )
    # Skip if no stock data
    if( nrow( stockData ) == 0 ) 
      next

    # Subset by sex
    stockDataBoys   <- stockData %>% filter( sex == 1 )
    stockDataGirls  <- stockData %>% filter( sex == 2 )

    stockList[[stockIdx]]$all <- makeMatOgives( matData = stockData, 
                                                maxPropMat = maxPropMat )

    if( nrow( stockDataBoys ) > 0 )
      stockList[[stockIdx]]$boys <- makeMatOgives(  matData = stockDataBoys,
                                                    maxPropMat = maxPropMat )

    if( nrow( stockDataGirls ) > 0 )
      stockList[[stockIdx]]$girls <- makeMatOgives( matData = stockDataGirls,
                                                    maxPropMat = maxPropMat )
  }

  outList <- list(  coastWide = cwList,
                    stocks = stockList )

  return(outList)
} # END makeSpecMat

# plotMatOgives()
# Takes output from makeSpecMat() and plots the ogives
# by species and stock, showing the data as proportions
# with confidence intervals, and showing the coastwide models
# overlaid on the stock specific models.
plotMatOgives <- function(  matList = matOgives,
                            save = FALSE,
                            saveDir = "Outputs",
                            type = "age" )
{
  # Get species names and stock names
  specNames   <- names(matList)
  stockNames  <- names( matList[[1]]$stocks )

  # Create save directory
  if(save)
  {
    graphics.off()
    saveDir <- file.path(getwd(),saveDir)
    if(! dir.exists(saveDir) )
      dir.create(saveDir)

    savePath <- file.path(saveDir,paste("mat",type,".png",sep = ""))

    png(  savePath, width = 11, height = 8.5, 
          units = "in", res = 300 )
  }

  # first, loop over species
  nSpec   <- length(specNames)
  nStocks <- length(stockNames)

  cols <- brewer.pal(n = 3, "Dark2")



  par( mfcol = c(nStocks, nSpec), mar = c(0,0,0,0), oma = c(4,4,3.5,3.5) )

  for( specIdx in 1:nSpec )
  {
    specID  <- specNames[specIdx]
    specMat <- matList[[specIdx]]

    coastWide <- specMat$coastWide
    stockFits <- specMat$stocks

    cwAll     <- coastWide$all[[type]]
    cwBoys    <- coastWide$boys[[type]]
    cwGirls   <- coastWide$girls[[type]]
    
    # Get max x value (age or length)
    maxX      <- max( cwAll$props[,type], na.rm = T )

    # Create coastwide ogives
    xSeq      <- seq(0,maxX, length = 100 )
    cwAll.y   <- logisticAtX(xSeq,cwAll$x50[1], cwAll$x95[1])
    cwBoys.y  <- logisticAtX(xSeq,cwBoys$x50[1], cwBoys$x95[1])
    cwGirls.y <- logisticAtX(xSeq,cwGirls$x50[1], cwGirls$x95[1])

    # Now loop over stocks
    for( stockIdx in 1:nStocks )
    {
      stockID <- stockNames[stockIdx]
      plot( x = c(0,maxX), y = c(0,1),
          type = "n", xlab = "", ylab = "",
          axes = FALSE )
        mfg <- par("mfg")
        # Add axes
        if( mfg[1] == 1 )
          mtext( side = 3, text = specID, font = 2, line = 2 )
        if( mfg[2] == mfg[4] )
          mtext( side = 4, text = stockID, font = 2, line = 2.5 )
        if( mfg[1] == mfg[3] )
          axis( side = 1 )
        if( mfg[2] == mfg[4] & mfg[1] %% 2 == 0 )
          axis( side = 4 )
        if( mfg[2] == 1 & mfg[1] %% 2 == 1 )
          axis( side = 2 )  
        box()
        # Plot coastwide models
        lines( x = xSeq, y = cwAll.y, lty = 2, lwd = .8,
                col = "grey50" )
        lines( x = xSeq, y = cwBoys.y, lty = 2, lwd = .8,
                col = cols[1] )
        lines( x = xSeq, y = cwGirls.y, lty = 2, lwd = .8,
                col = cols[2] )

      if( is.null(stockFits[[stockIdx]] ) )
      {
        panLab(x = 0.5, y = 0.5, txt = "No Data", font = 2 )
        next
        
      }

      # Pull stock specific models
      stockAll    <- stockFits[[stockIdx]]$all[[type]]
      stockBoys   <- stockFits[[stockIdx]]$boys[[type]]
      stockGirls  <- stockFits[[stockIdx]]$girls[[type]]

      boyProps    <-  stockBoys$props %>%
                      mutate( SE = sqrt(propMat*(1-propMat)/nObs) )
      girlProps   <- stockGirls$props %>%
                      mutate( SE = sqrt(propMat*(1-propMat)/nObs) )

      # Make stock specific lines
      stockAll.y    <- logisticAtX( xSeq, stockAll$x50[1], stockAll$x95[1]) 
      stockBoys.y   <- logisticAtX( xSeq, stockBoys$x50[1], stockBoys$x95[1]) 
      stockGirls.y  <- logisticAtX( xSeq, stockGirls$x50[1], stockGirls$x95[1]) 

      # Back to plotting
        # Plot data as distributions of prop mature at age
        points( x = (boyProps[[type]]) - .3, y = boyProps$propMat,
                pch = 16, cex = .8, col = cols[1]  )
        segments( x0 = boyProps[[type]] - .3, 
                  y0 = boyProps$propMat - boyProps$SE,
                  y1 = boyProps$propMat + boyProps$SE,
                  lwd = 1, col = cols[1]  )

        points( x = girlProps[[type]] + .3, y = girlProps$propMat,
                pch = 16, cex = .8, col = cols[2]  )
        segments( x0 = girlProps[[type]] + .3, 
                  y0 = girlProps$propMat - girlProps$SE,
                  y1 = girlProps$propMat + girlProps$SE,
                  lwd = 1, col = cols[2]  )

        # Plot stock specific models
        lines( x = xSeq, y = stockAll.y, lty = 1, lwd = 2,
                col = "grey50" )
        lines( x = xSeq, y = stockBoys.y, lty = 1, lwd = 2,
                col = cols[1] )
        lines( x = xSeq, y = stockGirls.y, lty = 1, lwd = 2,
                col = cols[2] )

      # Need to add data from props df - next
    }
  }
  if( type == "age" )
  {
    xLab <- "Age"
    yLab <- "Proportion Mature-at-age"
  }
  if( type == "length" )
  {
    xLab = "Length (cm)"
    yLab <- "Proportion Mature-at-length"
  }
  mtext( side = 1, outer = T, text = xLab, line = 3, font = 2, cex = 1.2 )
  mtext( side = 2, outer = T, text = yLab, line = 2, font = 2, cex = 1.2 )

  if(save)
    dev.off()
}

# logisticAtX()
# Takes a sequence of x values, and x at 50% and 95%
# and creates the logistic ogive across the x sequence
# logistic function with a50 and a95 as inputs
logisticAtX <- function( x, x50=5, x95=7 )
{
  # solve for step (scale), a50 is location
  s <- (x95 - x50) / log(19)
  # now compute logistic
  p <- 1 / (1 + exp((-x + x50)/s))
 
  p
}

