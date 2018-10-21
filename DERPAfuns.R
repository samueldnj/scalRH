# -----------------------------------------------------------------------------
#
# DERPAfuns.R
#
# Functions for DERPAdata.R and ageStructuredControl.R
# 
# -----------------------------------------------------------------------------

# OK, we need functions to summarise the biological data 
# for a given species, given a determined stock structure



# Reconstructing the time series from the fit to ensure that 
# everything is working ok
reconstruct <- function( s = 1, repObj )
{
  Bmsy <- repObj$Bmsy[s]
  Umsy <- repObj$Umsy[s]
  eps_t <- repObj$eps_t
  zeta_t <- repObj$zeta_st[s,]

  B0  <- 2*Bmsy
  r   <- 2*Umsy

  nT <- repObj$nT


  Bt <- numeric(length = nT)
  Ct <- numeric(length = nT)

  Ut <- repObj$Ut[s,]

  Bt[1] <- B0
  for( t in 2:nT )
  {
    Ct[t-1] <- Bt[t-1] * Ut[t-1]
    Bt[t] <- Bt[t-1] + r*Bt[t-1]*(1 - Bt[t-1]/B0) - Ct[t-1]
    Bt[t] <- Bt[t] * exp(eps_t[t-1] + zeta_t[t-1])
  }

  out <- cbind(1:nT,Bt,Ct)
  out
}

# Plotting function for plotting the assessment
# of DERPA data
plotAssess <- function( sdReport = sdrepObj,
                        report = repObj,
                        indices = relBio,
                        katch = Katch,
                        initYear = 1954,
                        depletion = TRUE )
{
  nSurv <- dim(indices)[1]
  nSpec <- dim(indices)[2]
  nT    <- dim(indices)[3]

  yrs <- seq(from = initYear, by = 1, length = nT )
  specNames <- dimnames(indices)[[2]]


  sdrep <- summary(sdReport)
  se.df <- data.frame(  par = dimnames(sdrep)[[1]],
                        est = sdrep[,"Estimate"],
                        se  = sdrep[,"Std. Error"] )
  se.df <-  se.df %>%
            mutate( CIlow = est - se,
                    CIhi  = est + se )

  # Now start arranging things in easy-to-plot-from
  # arrays
  # Get biomass
  Bt          <- array( 0, dim = c(nSpec,nT,3) )
  Bt[ ,,2 ]   <- matrix( se.df[which(se.df$par == "Bt"),"est"],
                         nrow = nSpec, ncol = nT, byrow = F )
  Bt[ ,,1 ]   <- matrix( se.df[which(se.df$par == "Bt"),"CIlow"],
                         nrow = nSpec, ncol = nT, byrow = F )
  Bt[ ,,3 ]   <- matrix( se.df[which(se.df$par == "Bt"),"CIhi"],
                         nrow = nSpec, ncol = nT, byrow = F )
  # Get Bmsy values
  Bmsy  <- report$Bmsy
  MSY   <- report$msy
  
  # Compute depletion
  Dt <- Bt[,,2]
  for( s in 1:nSpec ) Dt[s,] <- Dt[s,] / Bmsy[s]

  # Get process error values
  # Species specific
  zeta_st         <- array( 0, dim = c(nSpec,nT,3) )
  zeta_st[ , 2:nT, 2 ]  <- matrix(  se.df[which(se.df$par == "zeta_st"),"est"],
                              nrow = nSpec, ncol = nT-1, byrow = F )
  zeta_st[ , 2:nT, 1 ]  <- matrix( se.df[which(se.df$par == "zeta_st"),"CIlow"],
                             nrow = nSpec, ncol = nT-1, byrow = F )
  zeta_st[ , 2:nT, 3 ]   <- matrix( se.df[which(se.df$par == "zeta_st"),"CIhi"],
                              nrow = nSpec, ncol = nT-1, byrow = F )
  # shared mean term
  eps_t <- matrix( NA, nrow = 3, ncol = nT )
  if( any(grepl(pattern = "eps_t", x = se.df$par)) )
  {
    eps_t[2, 2:nT ] <- se.df[which(se.df$par == "eps_t"), "est"]
    eps_t[1, 2:nT ] <- se.df[which(se.df$par == "eps_t"), "CIlow"]
    eps_t[3, 2:nT ] <- se.df[which(se.df$par == "eps_t"), "CIhi"]  
  }
  

  # Now compute observation errors
  # replace missing indices with NAs and compute observation errors
  indices[ indices <= 0 ] <- NA
  q_os      <- report$q_os
  Itbar     <- array( 0, dim = c(nSurv, nSpec, nT ) )
  Itscaled  <- Itbar
  for( oIdx in 1:nSurv )
    for( t in 1:nT )
      {
        Itbar[oIdx,,t]    <- Bt[,t,2] * q_os[oIdx,]
        Itscaled[oIdx,,t] <- indices[oIdx,,t] / q_os[oIdx,]
      }

  delta_ost <- log(indices) - log(Itbar)

  # Make plotting area
  par( mfrow = c( 4, nSpec ), oma = c(3,3,3,1), mar = c(1,1,1,1) )
  for( sIdx in 1:nSpec )
  { 
    if( depletion ) range <- c(0,2)
    else  range <- c(0,max(Bt[sIdx,,]))
    DepScalar <- .5
    if( depletion ) DepScalar <- Bmsy[sIdx]
    plot( x = range(yrs) + c(-.4,.4), y = range, type = "n", axes = F )
      mtext( side = 3, text = specNames[sIdx], line = 2 )
      axis( side = 1 )
      axis( side = 2, las = 1 )
      polygon(  x = c(yrs,rev(yrs)), y = c(Bt[sIdx,,1],rev(Bt[sIdx,,3]))/DepScalar/2,
                border = NA, col = alpha("grey80",.5) )
      lines( x = yrs, y = Bt[sIdx,,2]/DepScalar/2, lwd = 2 )
      for( oIdx in 1:nSurv )
        points( x = yrs, y = Itscaled[oIdx,sIdx,]/DepScalar/2, pch = oIdx )
      if( depletion ) titleText <- "Relative Depletion" else titleText <- "Biomass (Kt)"
      if( sIdx == 1 ) mtext(  side = 2, text = titleText, 
                              line = 2 )
  }
  for( sIdx in 1:nSpec )
  {
    depScalar <- MSY[sIdx]
    titleText <- "Catch/MSY"
    range <- c(0,3)
    if( !depletion ) 
    {
      depScalar <- 1
      titleText <- "Catch (Kt)"
      range <- c(0,max(katch[sIdx,],na.rm=T))
    }
    plot( x = range(yrs) + c(-.4,.4), y = range, type = "n", axes = F )
      axis( side = 1 )
      axis( side = 2, las = 1 )
      rect( xleft = yrs - .4, xright = yrs + .4,
            ybottom = rep(0,length(yrs)), ytop = katch[ sIdx, ]/depScalar, border = NA, 
            col = "grey30" )
      if( sIdx == 1 ) mtext(  side = 2, text = titleText, 
                              line = 2 )
  }
  
  for( sIdx in 1:nSpec )
  {
    plot( x = range(yrs) + c(-.4,.4), y = c( -2,2 ), type = "n", axes = F )
      axis( side = 1 )
      axis( side = 2, las = 1 )
      for( oIdx in 1:nSurv )
      {
        points( x = yrs, y = delta_ost[oIdx,sIdx,], pch = oIdx, col = oIdx )
        surveyOn <- which( !is.na( delta_ost[ oIdx, sIdx, ] ) )
        lines( lowess(delta_ost[oIdx,sIdx,surveyOn] ~ yrs[surveyOn] ), col = oIdx, lwd = 2 )
      }
      if( sIdx == 1 ) mtext(  side = 2, text = "Observation Residual",
                              line = 2 )
  }
  for( sIdx in 1:nSpec )
  {
    plot( x = range(yrs) + c(-.4,.4), y = c( 0.8,1.2 ), type = "n", axes = F )
      axis( side = 1 )
      axis( side = 2, las = 1 )
      points( x = yrs, y = exp(zeta_st[sIdx,,2]), pch = 17, col = 2 )
      segments( x0 = yrs, y0 = exp(zeta_st[sIdx,,1]), y1 = exp(zeta_st[sIdx,,3]),
                lwd = 1, col = 2 )
      points( x = yrs, y = exp(eps_t[2,] ), pch = 16 )
      segments( x0 = yrs, y0 = exp(eps_t[1,]), y1 = exp(eps_t[3,]),
                lwd = 1 )
        
      if( sIdx == 1 ) mtext(  side = 2, text = "Proc. Err. Deviations",
                              line = 2 )
  }
  mtext( side = 1, text = "Years", outer = T, line = 2 )
}

# Makes data and parameter lists for fitting the 
# msProd model to DERPA data
makeDatPar <- function( relB = relBio, Cst = Katch )
{

  # Recover model dimensions
  nSurv   <- dim( relB )[ 1 ]
  nStocks <- dim( relB )[ 2 ]
  nT      <- dim( relB )[ 3 ]

  # Make data list
  data    <- list(  It              = relB, 
                    Ct              = Cst,
                    SigmaPriorCode  = 0,
                    sigUPriorCode   = 0,
                    tauqPriorCode   = 0,
                    lnqPriorCode    = 1,
                    lnUPriorCode    = 0,
                    BmsyPriorCode   = 0,
                    initT           = c(0,0,0,0,0),
                    initBioCode     = c(1,1,1,1,1) )

  
  
  # Sum catch for initial Bmsy guess
  sumCat <- apply( X = Cst, FUN = sum, MARGIN = 1 )

  # populate parameter list
  para <- list( lnBmsy            = log( 2*sumCat ),
                lnUmsy            = rep(  -2, nStocks ),
                lntau2_o          = rep(  0, nSurv ),
                lnq_os            = matrix( 0, nrow = nSurv, ncol = nStocks ),
                lnBinit           = log( sumCat ),        # Start every stock fished
                lnqbar_o          = rep( 0, nSurv ),
                lntauq_o          = rep( 0, nSurv ),
                mq                = .6,
                sq                = 0.2,
                lnUmsybar         = -2,
                lnsigUmsy         = -1,
                mUmsy             = 0.1,
                sUmsy             = 0.1,
                mBmsy             = sumCat,
                sBmsy             = sumCat/2,
                tau2IGa           = rep( 2, nSurv ),
                tau2IGb           = rep( 1, nSurv ),
                tauq2Prior        = c( 4, 0.5 ),
                sigU2Prior        = c( 0.5, 0.15 ),
                kappa2IG          = c( 1, 0.05 ),
                Sigma2IG          = c( 1, 0.05 ),
                wishScale         = diag( 1, nStocks ),
                nu                = nStocks,
                eps_t             = rep( 0, nT - 1 ),
                lnkappa2          = -10,
                zeta_st           = matrix( 0, nrow = nStocks, ncol = (nT - 1) ),
                lnSigmaDiag       = -3,
                SigmaDiagMult     = rep( 1, nStocks ),
                logitSigmaOffDiag = rep( 0, length = nStocks*(nStocks-1)/2 ),
                logit_gammaYr     = 0 )

  map <- list ( mBmsy             = factor( rep( NA, nStocks ) ),
                sBmsy             = factor( rep( NA, nStocks ) ),
                # lnBinit           = factor( c(NA,NA,13,13,NA) ),
                # lnq_os            = factor( matrix(NA, nSurv, nStocks) ),
                # lnqbar_o          = factor( rep(NA,nSurv)),
                # lntauq_o          = factor( rep(NA,nSurv)),
                lnUmsybar         = factor( NA ),
                lnsigUmsy         = factor( NA ),
                mUmsy             = factor( NA ),
                sUmsy             = factor( NA ),
                mq                = factor( NA ),
                sq                = factor( NA ),
                SigmaDiagMult     = factor( rep( NA, nStocks ) ),
                tau2IGa           = factor( rep( NA, nSurv ) ),
                tau2IGb           = factor( rep( NA, nSurv ) ),
                sigU2Prior        = factor( rep( NA, 2 ) ),
                tauq2Prior        = factor( rep( NA, 2 ) ),
                kappa2IG          = factor( rep( NA, 2 ) ),
                Sigma2IG          = factor( rep( NA, 2 ) ),
                wishScale         = factor( rep( NA, nStocks*nStocks ) ),
                nu                = factor( NA ),
                logit_gammaYr     = factor( NA ),
                logitSigmaOffDiag = factor( rep( NA, length = nStocks*(nStocks-1)/2 ) ), 
                lnkappa2          = factor( NA ),
                eps_t             = factor( rep(NA, nT - 1) )
              )

  out <- list( dat = data, par = para, map = map )
  out
}

# Split catch history by species and arrange in an array
# for feeding to TMB model
makeStockCatch <- function( years = c(1975,2016),
                            stocks = list(  HG    = c(7,8,9),
                                            QCS   = c(5,6),
                                            WCVI  = c(3,4) ),
                            spec = "dover"
                         )
{

  # function to replace codes with names in df
  dfCodeToName <- function( code )
  {
    # Create the same for the species
    specNames <- c("dover","english","rock","petrale","atooth")
    specCodes <- c( 626, 628, 621, 607, 602 )
    names( specCodes ) <- specNames

    cIdx <- which( specCodes == code[1] )
    name <- names( specCodes )[ cIdx ]
    name
  }

  # Read in commercial catch data
  commCatch <- read.csv( "catch_by_maj.csv", header = T, stringsAsFactors = F)

  appendStockName <- function( majorAreaCode, stockList )
  {
    majCodes <- unlist(stockList)
    stockVec <- rep(NA,max(majCodes))
    for( i in 1:length(stockList) )
    {
      stockVec[ stockList[[i]] ] <- names(stockList)[i]
    }

    stockVec[majorAreaCode]
  }


  # Mung the commercial catch data
  catchAreaFleetSpecies <-  commCatch %>%
                            group_by( YEAR, MAJOR_STAT_AREA_CODE, 
                                      GEAR, SPECIES_CODE ) %>%
                            summarise(  landedWt = sum( LANDED_KG ) / 1e6,
                                        discardWt = sum( DISCARDED_KG ) / 1e6,
                                        discardPc = sum( DISCARDED_PCS ),
                                        specName = dfCodeToName( SPECIES_CODE ) ) %>%
                            ungroup() %>%
                            dplyr::select(  year = YEAR, 
                                            majorStatArea = MAJOR_STAT_AREA_CODE,
                                            gearType = GEAR,
                                            specCode = SPECIES_CODE,
                                            specName,
                                            landedWt,
                                            discardWt,
                                            discardPc ) %>%
                            filter( majorStatArea %in% unlist(stocks) & specName == spec ) %>%
                            mutate( stockName = appendStockName( majorStatArea, stocks) ) %>%
                            filter( year >= years[1] & year <= years[2] )

  catchSpecies <- catchAreaFleetSpecies %>%
                  group_by(stockName, year) %>%
                  summarise( katch = sum(landedWt + discardWt) ) 


  # Array dimensions and labels
  stockIDs <- paste(spec,names(stocks),sep = "")
  yrs <- years[1]:years[2]

  # arrange in an array
  catchSpecArray <- array( 0, dim = c(length(stocks),length(yrs)),
                              dimnames = list( stockIDs, yrs ) )
  
  for( sIdx in 1:length(stockIDs) )
  {
    subCatch <- catchSpecies %>%
                filter(stockName == names(stocks)[sIdx] )

    catchSpecArray[stockIDs[sIdx],] <-  subCatch$katch
  }

  out <- list(  catch.df  = catchSpecies,
                catch.arr = catchSpecArray )

  out
}

# Calculate relative biomass by species and arrange in an array
# for feeding to TMB model
makeRelBioStocks <- function( years = c(1975, 2016), 
                              collapseSyn = TRUE,
                              stocks = list(  HG = c(2,3,16),
                                              QCS = c(1),
                                              WCVI = c(4) ),
                              spec = "dover"
                            )
{
  # Read in strata areas
  strata <- read.csv( "derpa_strata.csv", header=T )
  stratArea <-  strata %>%
                dplyr::select( GROUPING_CODE, AREA_KM2 )
  # Survey ids for plotting/legends
  survIDs <-  c(  QCSyn = 1, HSAss = 2, HSSyn = 3, WCVISyn=4,
                  QCShr = 6, WCVIShr = 7, WCHGSyn = 16 )

  # Read in density
  specDensityFile <- paste(spec,"density.csv",sep = "_")
  specDensityPath <- file.path(getwd(),"density",specDensityFile)
  densityTab <- read.csv(specDensityPath, header = T)

  # Combine density frames based on trip and trawl IDs
  # First, rename the catch column in each df
  densityTab <- densityTab %>% mutate(  catch = CATCH_WEIGHT,
                                        density = DENSITY_KGPM2 )

  appendStockName <- function( surveyID, stockList )
  {
    majCodes <- unlist(stockList)
    stockVec <- rep(NA,max(majCodes))
    for( i in 1:length(stockList) )
    {
      stockVec[ stockList[[i]] ] <- names(stockList)[i]
    }
    stockVec[surveyID]
  }

  includedSurveys <- unlist(stocks)

  survIDs <- survIDs[survIDs %in% includedSurveys]

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
                                density, catch  ) %>%
                filter( survSeriesID %in% includedSurveys ) %>%
                mutate( stockName = appendStockName( survSeriesID, stocks ) )

  relativeBio <-  surveyData %>%
                  filter( survSeriesID %in% survIDs) %>%
                  group_by( year, stockName, survSeriesID, stratum ) %>%
                  dplyr::summarise( density = mean(density),
                                    area = mean(stratArea) ) %>%
                  dplyr::summarize( relBio_Kt = sum( area * density),
                                    surveyArea = sum(area) ) %>%
                  ungroup() %>%
                  filter( year >= years[1] & year <= years[2] )

  if( collapseSyn )
  {
    # isolate synoptic IDs
    synIDs <-  c( QCSyn = 1, HSSyn = 3, WCVISyn=4,
                  WCHGSyn = 16 )
    synStratAreas <- surveyData %>%
                    filter( survSeriesID %in% synIDs ) %>%
                    group_by( stockName, survSeriesID, stratum ) %>%
                    summarise( stratArea = unique( stratArea ) ) %>%
                    summarise( surveyArea = sum( stratArea ) ) %>%
                    ungroup()

    synStockAreas <-  synStratAreas %>%
                      group_by( stockName ) %>%
                      summarise( stockArea = sum( surveyArea ) ) %>%
                      ungroup()


    synRelativeBio <- relativeBio %>% 
                      filter( survSeriesID %in% synIDs ) %>%
                      left_join( synStockAreas, by = "stockName" ) %>%
                      group_by( stockName, year ) %>%
                      dplyr::summarize( area = sum(surveyArea),
                                        relBio_Kt = sum( relBio_Kt ),
                                        stockArea = unique(stockArea) ) %>%
                      dplyr::mutate(  relBio_Kt = relBio_Kt * stockArea / area,
                                      survSeriesID = 1 ) %>%
                      dplyr::select(  year,
                                      stockName,
                                      survSeriesID,
                                      relBio_Kt,
                                      surveyArea = area,
                                      stockArea = stockArea
                                    )

    relativeBio <-  relativeBio %>%
                    left_join( synStockAreas, by = "stockName" ) %>%
                    filter( !(survSeriesID %in% synIDs) )

    relativeBio <- rbind(as.data.frame(relativeBio),as.data.frame(synRelativeBio))

    survIDs <- intersect(survIDs[c(1,2,5,6)],includedSurveys)

  }

  yrs <- years[1]:years[2]
  stockNames <- paste(spec,names(stocks),sep = "")
  relativeBioArray <- array(  -1, dim = c(length(survIDs),length(stocks),length(yrs)),
                              dimnames = list(  names(survIDs), stockNames,
                                                yrs ) )

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
        relativeBioArray[survIdx, sIdx, yIdx] <- subRelBio$relBio_Kt
        
      }
      
    }
  }

  out <- list(  relBio.df = relativeBio,
                relBio.arr = relativeBioArray)
  out
}

# Calculate relative biomass by species and arrange in an array
# for feeding to TMB model
makeRelativeBio <- function( years = c(1975, 2016), collapseSyn = TRUE )
{
  specNames <- c("dover","english","rock","petrale","atooth")
  
  # Read in strata areas
  strata <- read.csv( "derpa_strata.csv", header=T )
  stratArea <-  strata %>%
                dplyr::select( GROUPING_CODE, AREA_KM2 )
  # Survey ids for plotting/legends
  survIDs <-  c(  QCSyn = 1, HSAss = 2, HSSyn = 3, WCVISyn=4,
                  QCShr = 6, WCVIShr = 7, WCHGSyn = 16 )

  # Read in density
  doverDensity    <- read.csv("./density/dover_density.csv",header=T)
  englishDensity  <- read.csv("./density/english_density.csv",header=T)
  rockDensity     <- read.csv("./density/srock_density.csv",header=T)
  petraleDensity  <- read.csv("./density/petrale_density.csv",header=T)
  atoothDensity   <- read.csv("./density/atooth_density.csv",header=T)

  # Combine density frames based on trip and trawl IDs
  # First, rename the catch column in each df
  doverDensity <- doverDensity %>% mutate(  doverCatch = CATCH_WEIGHT,
                                            doverDensity = DENSITY_KGPM2 )
  englishDensity <- englishDensity %>% mutate(  englishCatch = CATCH_WEIGHT,
                                                englishDensity = DENSITY_KGPM2 )
  rockDensity <- rockDensity %>% mutate(  rockCatch = CATCH_WEIGHT,
                                          rockDensity = DENSITY_KGPM2 )
  petraleDensity <- petraleDensity %>% mutate(  petraleCatch = CATCH_WEIGHT,
                                                petraleDensity = DENSITY_KGPM2 )
  atoothDensity <- atoothDensity %>% mutate(  atoothCatch = CATCH_WEIGHT,
                                              atoothDensity = DENSITY_KGPM2 )
  # Now join and select the columns we want
  surveyData <- doverDensity %>%
                left_join( englishDensity, by = "FISHING_EVENT_ID" ) %>%
                left_join( rockDensity, by = "FISHING_EVENT_ID" ) %>%
                left_join( petraleDensity, by = "FISHING_EVENT_ID" ) %>%
                left_join( atoothDensity, by = "FISHING_EVENT_ID" ) %>%
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
                                doverDensity, doverCatch,
                                englishDensity, englishCatch,
                                rockDensity, rockCatch,
                                petraleDensity, petraleCatch,
                                atoothDensity, atoothCatch  ) %>%
                mutate( DERPAcatch =  doverCatch + 
                                      englishCatch + 
                                      rockCatch + 
                                      petraleCatch + 
                                      atoothCatch )

  relativeBio <-  surveyData %>%
                  group_by( year, survSeriesID, stratum ) %>%
                  dplyr::summarise( doverDensity = mean(doverDensity),
                                    englishDensity = mean(englishDensity),
                                    rockDensity = mean(rockDensity),
                                    petraleDensity = mean(petraleDensity),
                                    atoothDensity = mean(atoothDensity),
                                    area = mean(stratArea) ) %>%
                  dplyr::summarize( relBioDover_Kt = sum( area * doverDensity),
                                    relBioEnglish_Kt = sum( area * englishDensity),
                                    relBioRock_Kt = sum( area * rockDensity),
                                    relBioPetrale_Kt = sum( area * petraleDensity),
                                    relBioAtooth_Kt = sum( area * atoothDensity),
                                    surveyArea = sum(area) ) %>%
                  ungroup() %>%
                  filter( year >= years[1] & year <= years[2] )

  if( collapseSyn )
  {
    # isolate synoptic IDs
    synIDs <-  c( QCSyn = 1, HSSyn = 3, WCVISyn=4,
                  WCHGSyn = 16 )
    synStratAreas <- surveyData %>%
                    filter( survSeriesID %in% synIDs ) %>%
                    group_by( survSeriesID, stratum ) %>%
                    summarise( stratArea = unique( stratArea ) ) %>%
                    summarise( surveyArea = sum( stratArea ) ) %>%
                    ungroup()

    totalSynArea <- sum(synStratAreas $ surveyArea )

    synRelativeBio <- relativeBio %>% 
                      filter( survSeriesID %in% synIDs ) %>%
                      group_by( year ) %>%
                      dplyr::summarize( area = sum(surveyArea),
                                        relBioDover_Kt = sum( relBioDover_Kt ),
                                        relBioEnglish_Kt = sum( relBioEnglish_Kt ),
                                        relBioRock_Kt = sum( relBioRock_Kt ),
                                        relBioPetrale_Kt = sum( relBioPetrale_Kt ),
                                        relBioAtooth_Kt = sum( relBioAtooth_Kt ) ) %>%
                      dplyr::mutate(  relBioDover_Kt = relBioDover_Kt * totalSynArea / area,
                                      relBioEnglish_Kt = relBioEnglish_Kt * totalSynArea / area,
                                      relBioRock_Kt = relBioRock_Kt * totalSynArea / area,
                                      relBioPetrale_Kt = relBioPetrale_Kt * totalSynArea / area,
                                      relBioAtooth_Kt = relBioAtooth_Kt * totalSynArea / area,
                                      survSeriesID = 1 ) %>%
                      dplyr::select(  year,
                                      survSeriesID,
                                      relBioDover_Kt,
                                      relBioEnglish_Kt,
                                      relBioRock_Kt,
                                      relBioPetrale_Kt,
                                      relBioAtooth_Kt,
                                      surveyArea = area
                                    )

    relativeBio <-  relativeBio %>%
                    filter( !(survSeriesID %in% synIDs) ) %>%
                    rbind( synRelativeBio )


    survIDs <- survIDs[c(1,2,5,6)]

  }

  yrs <- years[1]:years[2]
  relativeBioArray <- array(  -1, dim = c(length(survIDs),5,length(yrs)),
                              dimnames = list(  names(survIDs), specNames,
                                                yrs ) )

  for( survIdx  in 1:length(survIDs) )
  {
    for( yIdx in 1:length(yrs) )
    {
      subRelBio <-  relativeBio %>%
                    filter( year == yrs[yIdx], survSeriesID == survIDs[survIdx] )
      if( nrow(subRelBio) == 0 ) next
      relativeBioArray[ survIdx, "dover", yIdx ] <- subRelBio$relBioDover_Kt
      relativeBioArray[ survIdx, "english", yIdx ] <- subRelBio$relBioEnglish_Kt
      relativeBioArray[ survIdx, "rock", yIdx ] <- subRelBio$relBioRock_Kt
      relativeBioArray[ survIdx, "petrale", yIdx ] <- subRelBio$relBioPetrale_Kt
      relativeBioArray[ survIdx, "atooth", yIdx ] <- subRelBio$relBioAtooth_Kt
    }
  }

  out <- list(  relBio.df = relativeBio,
                relBio.arr = relativeBioArray)
  out
}

# Start creating a function to plot the DERPA species extents as
# a raster grid
densityToRaster <- function(  catchDF = surveyData,
                              yrRange = c( 1975, 2016 ),
                              oldProj = CRS("+proj=longlat +datum=WGS84"),
                              newProj = UTMproj,
                              dummRes = 50,
                              plotExt = c(  xmin = 170000, 
                                            xmax = 820000,
                                            ymin = 5300000, 
                                            ymax = 6061000 ),
                              saveFile = "surveyCatchRaster.RData" )
{
  # First, convert the df to a spatial points data frame
  # Extract coordinates and data, filter by year
  catchDF <-  catchDF %>% 
              dplyr::filter( (year >= yrRange[1]) & (year <= yrRange[2]) )

  relBioData <- catchDF %>% dplyr::select(  doverDensity,
                                            englishDensity,
                                            rockDensity,
                                            petraleDensity,
                                            atoothDensity )
  relBioLoc  <- catchDF %>% dplyr::select(  x = lon, y = lat )

  spDF <- SpatialPointsDataFrame( coords = relBioLoc, data = relBioData,
                                  proj4string = oldProj )

  if( !is.null(newProj) )
    spDF <- spTransform ( spDF, newProj )

  if( is.null(newProj) ) newProj <- oldProj

  # Now aggregate into a raster, with the desired spatial grid.
  # First create a dummy raster
  dummy <- raster(  crs = newProj,
                    res = dummRes,
                    xmn = plotExt[1], xmx = plotExt[2],
                    ymn = plotExt[3], ymx = plotExt[4] )

  relBioRaster <- rasterize(  x = spDF, y = dummy,
                              fun = mean,
                              field = c(  "doverDensity",
                                          "englishDensity",
                                          "rockDensity",
                                          "petraleDensity",
                                          "atoothDensity" ) )

  relBioRaster <- relBioRaster * 2500

  relBioRaster <- raster::aggregate( x = relBioRaster, fun = sum, fact = 200 )

  relBioRaster$DERPADensity <- sum(relBioRaster)

  save( relBioRaster, file = saveFile )

  relBioRaster
  
}


# Start creating a function to plot the DERPA species extents as
# a raster grid
catchToRaster <- function(  catchDF = surveyData,
                              yrRange = c( 1984, 2016 ),
                              oldProj = CRS("+proj=longlat +datum=WGS84"),
                              newProj = UTMproj,
                              plotExt = c(  xmin = 170000, 
                                            xmax = 820000,
                                            ymin = 5300000, 
                                            ymax = 6061000 ),
                              saveFile = "surveyCatchRaster.RData" )
{
  # First, convert the df to a spatial points data frame
  # Extract coordinates and data, filter by year
  catchDF <-  catchDF %>% 
              dplyr::filter( (year >= yrRange[1]) & (year <= yrRange[2]) )

  catchData <- catchDF %>% dplyr::select( doverCatch,
                                          englishCatch,
                                          rockCatch,
                                          petraleCatch,
                                          atoothCatch,
                                          DERPAcatch )
  catchLoc  <- catchDF %>% dplyr::select(  x = lon, y = lat )

  spDF <- SpatialPointsDataFrame( coords = catchLoc, data = catchData,
                                  proj4string = oldProj )

  spDF <- spTransform ( spDF, newProj )

  # Now aggregate into a raster, with the desired spatial grid.
  # First create a dummy raster
  dummy <- raster(  crs = newProj,
                    res = 50,
                    xmn = plotExt[1], xmx = plotExt[2],
                    ymn = plotExt[3], ymx = plotExt[4] )

  catchRaster <- rasterize( x = spDF, y = dummy,
                            field = c(  "doverCatch",
                                        "englishCatch",
                                        "rockCatch",
                                        "petraleCatch",
                                        "atoothCatch",
                                        "DERPAcatch") )

  catchRaster <- raster::aggregate( x = catchRaster, fun = sum, fact = 200 )
  # Now replace zeroes with NAs
  catchRaster[catchRaster == 0] <- NA

  save( catchRaster, file = saveFile )

  catchRaster
  
}





plotCatchMap <- function( catchRaster = relBioRaster,
                          colBreaks = c(1,10,50,100),
                          base = nepacLL,
                          scale = 1,
                          newProj = UTMproj,
                          oldProj = LLproj,
                          plotExt = c(  xmin = 170000, 
                                        xmax = 820000,
                                        ymin = 5300000, 
                                        ymax = 6061000 ),
                          legend = TRUE, label = TRUE,
                          units = "kg",
                          quant = expression(B[trawl]),
                          mgmtAreas = NULL,
                          saveFile = NULL,
                          multiPanel = TRUE,
                          width = 4, height = 4,
                          axisTicks = FALSE,
                          leg.cex = 1, axes = T, labLL = T )
{
  # First, transform base to a spatialPolyDataFrame, project so
  # it's the same as the raster
  base.sp <- PolySet2SpatialPolygons(base)
  if( !is.null(newProj) ) 
    base.sp <- spTransform( base.sp, newProj )

  # Scale raster to get desired units
  catchRaster <- catchRaster / scale

  # Pad colBreaks
  colBreaks <- c(0,colBreaks,max(values(catchRaster),na.rm=T))
  cols <- brewer.pal( n = length(colBreaks)-1, name = "YlOrRd")
  nCols <- length(cols)
  legText <- vector(mode = "list", length = nCols)
  for( c in 1:(nCols-1) )
  {
    legText[[c]] <- substitute( paste( colB1, " kg < ", B[trawl], " < ", colB2, " kg", sep = "" ),
                              list( colB1 = colBreaks[c], colB2 = colBreaks[c+1] ) )
  }
  legText[[nCols]] <- substitute( paste( B[trawl],  " > ", colB, u, sep = "" ),
                                list(colB = colBreaks[nCols], u = units ) )

  if( !is.null(saveFile) )
  {
    png(  filename = saveFile,
          width = width, height = height, units = "in", res = 300 )
  }

  pullCol <- function( value, breaks, cols )
  {
    colIdx <- max( which( breaks <= value ) )
    return( cols[colIdx] )
  }
  
  # Now start plotting?
  if( multiPanel ) par( mfrow = c(2,3), oma = c(1,1,1,1), mar = c(.5,.5,.5,.5))
  else par( oma = c(3,3,3,3), mar = c(0,0,0,0))

  # Convert to LL labels
  if( labLL )
  {
    X   <- seq( -134, -124, by = 2 )
    Y   <- seq( 48, 54, by = 1 )

    Xcoords   <- data.frame(x = X, y = rep(48,length(X)) )
    spCoords  <- SpatialPoints(coords = Xcoords, proj4string = oldProj)
    if(!is.null(newProj))
      newCoords <- spTransform( spCoords, newProj) 
    else newCoords <- spCoords

    X.new <- as.data.frame(newCoords)[,1]

    Y.utm   <- Y * 111320

    Ylab <- paste( Y, "N")
    Xlab <- paste( -X, "W" )
  }


  for( k in 1:length(names(catchRaster)) )
  {
    catchPoly <- rasterToPolygons( x = catchRaster[[k]] )
    if(is.null(newProj)) catchPoly <- spTransform( catchPoly, oldProj )
    polyEntries <- catchPoly[[1]]
    colVec <- sapply( X = polyEntries, FUN = pullCol, breaks = colBreaks, cols = cols )

    plot( x = catchPoly, col = colVec, 
          las = 1, border = "grey30", frame.plot = T, axes = T,
          xaxt = "n", yaxt=  "n" )
      if( axes & multiPanel )
      {
        if( k == 1 | k == 4 )
        { 
          axis( side = 2, labels = Ylab, at = Y.utm, las = 1 )
          axis( side = 4, labels = Ylab, at = Y.utm, las = 1 )
        }
        if( k == 2 | k == 5 )
        {
          axis( side = 2, labels = Ylab, at = Y.utm, las = 1 )
          axis( side = 4, labels = Ylab, at = Y.utm, las = 1 )
        }
        if( k <= 3 ) 
        {
          axis( side = 3, labels = Xlab, at = X.new )
          axis( side = 1, labels = Xlab, at = X.new )
        }
        if( k == 3 | k == 6 ) 
        {
          axis( side = 4, labels = Ylab, at = Y.utm, las = 1 )
          axis( side = 2, labels = Ylab, at = Y.utm, las = 1 )
        }
        if( k > 3 )
        {
          axis( side = 1, labels = Xlab, at = X.new )
          axis( side = 3, labels = Xlab, at = X.new )
        }  
      }
      if( axes & !multiPanel )
      {
        axis( side = 1, labels = Xlab, at = X.new )
        axis( side = 2, labels = Ylab, at = Y.utm, las = 1 )
      }
      if(!is.null(mgmtAreas))
        plotMgmtAreas( mgmtAreas, newProj )
      plot( base.sp,
            col = "grey60", add = T, border = NA )
      if( legend & k == 1 )
        panLegend( "bottomleft", legTxt = as.expression(legText), fill = cols,
                   bty = "n", border = "grey30", cex = leg.cex )
      if( label ) panLab( x = .9, y = .9, txt = names(catchRaster)[k] )
  }

  if( !is.null(saveFile) ) dev.off()
}

plotMgmtAreas <- function( areaShp, proj )
{ 
  # Label positions in lat/lon
  x   = c( -128.5, -130.5, -131, -131.5, -134, -130.5, -131.5 )
  y   = c( 48.4,49.7,50.9,51.6,53, 54.25, 54.4)
  lab = c( "3C", "3D", "5A", "5B", "5E", "5C", "5D" )
  # lonlat projection
  LLproj <- CRS("+proj=longlat +datum=WGS84")

  if(class(areaShp)[1] == "PolySet")
    areaShp <- PolySet2SpatialPolygons(areaShp)

  areaShp <- spTransform( areaShp, proj )

  # Make a spatial points object
  areaLabels <- data.frame( x = x, y = y, label = lab )
  aLabels.xy <- areaLabels[,c("x","y")]
  aLabels.sp <- SpatialPoints(coords = aLabels.xy, proj4string = LLproj )
  # Convert to new projection
  aLabels.sp <- spTransform( aLabels.sp, proj )
  areaLabels$x <- aLabels.sp$x
  areaLabels$y <- aLabels.sp$y

  plot(areaShp, add = TRUE, lwd = 1 )

  # Plot labels
  for( aIdx in 1:nrow(aLabels.xy) )
  {
    text( x=aLabels.xy[aIdx,"x"], y=aLabels.xy[aIdx,"y"], aLabels.xy[aIdx,"label"], 
          cex=1.5, font=2, col="black" )
  }

  # Management areas
  text( x=-128.5, y=48.4, "3C", cex=1.5, font=2, col="black" )
  text( x=-130.5, y=49.7, "3D", cex=1.5, font=2, col="black" )
  text( x=-122.9, y=50.2, "4B", cex=1.5, font=2, col="black" )
  text( x=-131, y=50.9, "5A", cex=1.5, font=2, col="black" )
  text( x=-131.5, y=51.6, "5B", cex=1.5, font=2, col="black" )
  text( x=-134, y=53, "5E", cex=1.5, font=2, col="black" )

  # Rectangles under 5C and 5D
  rect( xleft=-130.85, ybottom=52.55, xright=-130.15, ytop=52.85, 
      col=rgb(1,1,1,0.9), border=NA )
  text( x=-130.5, y=52.7, "5C", cex=1.5, font=2, col="black" )
  rect( xleft=-131.85, ybottom=54.25, xright=-131.15, ytop=54.55, 
      col=rgb(1,1,1,0.9), border=NA )
  text( x=-131.5, y=54.4, "5D", cex=1.5, font=2, col="black" )

  # Rectangles under 5C and 5D
  # rect( xleft=-130.85, ybottom=52.55, xright=-130.15, ytop=52.85, 
  #     col=rgb(1,1,1,0.9), border=NA )
  # text( x=-130.5, y=52.7, "5C", cex=1.5, font=2, col="black" )
  # rect( xleft=-131.85, ybottom=54.25, xright=-131.15, ytop=54.55, 
  #     col=rgb(1,1,1,0.9), border=NA )
  # text( x=-131.5, y=54.4, "5D", cex=1.5, font=2, col="black" )



  # Vancouver Island
  # text( x=-125.45, y=49.54, "Vancouver    Isl.", srt=-43, cex=1.5, col="black" )

  # # North arrow
  # xAr = -119.9 ;  yAr = 50.65
  # points( x=xAr, y=yAr, cex=3, pch=24, col="black", bg="black" )
  # lines( x=c(xAr, xAr), y=c(yAr-0.5, yAr), lwd=3, col="black" )
  # text( x=xAr, y=yAr+0.5, col="black", "N", cex=1.5, font=2 )

  # # Point to 4B area
  # arrows( x0=-123.6, x1=-123.7, y0=49.35, y1=49.2, lwd=4, code=2, length=0.125,
  #     col=rgb(1,1,1,0.9) )  
  # arrows( x0=-123.2, x1=-123.7, y0=50, y1=49.2, lwd=2, code=2, length=0.125,
  #     col="black" )

}

plotCatch <- function(  df = catchAreaFleetSpecies,
                        saveRoot = "commCatch" )
{
  # Create a vector of areas and their codes in the data
  areas <- c( 'unspecified',
              '4B',
              '3C',
              '3D',
              '5A',
              '5B',
              '5C',
              '5D',
              '5E')
  areaCodes <- c(0,1,3:9)
  names(areaCodes) <- areas

  # Create the same for the species
  species <- c( "Dover Sole", "English Sole", "Rock Sole",
                "Petrale Sole", "Arrowtooth Flounder" )
  specCodes <- c( 626, 628, 621, 607, 602 )
  names( specCodes ) <- species

  # List fleets
  fleets <- unique( df$gearType )

  # Now plot a multipanel plot, with rows as species, columns as areas
  yrs       <- unique( df$year )

  # first, plot a grid of species/area catch
  par( mfrow = c(length(specCodes), length(areaCodes) + 1 ), 
       mar = c(1,1,1,1), oma = c(3,3,3,3) )
  for( sIdx in 1:length(specCodes) )
  {
    specCatch <-  df %>%
                  filter( specCode == specCodes[sIdx])
    specTotalCatch <- specCatch %>%
                      group_by(year) %>%
                      summarise(  landedWt = sum( landedWt ),
                                  discardWt = sum( discardWt ) ) %>%
                      ungroup()
    maxCatch <- max(specTotalCatch$landedWt)
    for( aIdx in 1:length(areaCodes) )
    {
      specAreaCatch <-  specCatch %>%
                        filter( majorStatArea == areaCodes[aIdx] ) %>%
                        group_by(year) %>%
                        summarise( landedWt = sum(landedWt),
                                   discardWt = sum(discardWt) ) %>%
                        ungroup()
      plot( x = range(yrs), y = c( 0, maxCatch ),
            type = "n", axes = F )
        axis( side = 1 )
        axis( side = 2, las = 1)
        if( sIdx == 1 ) mtext(  side = 3, text = paste( names(areaCodes)[aIdx], sep = "" ), 
                                cex = 0.8, line = 2 )
        if( aIdx == 1 ) mtext(  side = 2, las = 0, text = names(specCodes)[sIdx], 
                                cex = 0.8, line = 2 )
        if( nrow(specAreaCatch) == 0 ) next
        lines(  x = specAreaCatch$year,
                y = specAreaCatch$landedWt, lwd = 2 )
        lines( x = specAreaCatch$year, y = specAreaCatch$discardWt, lty = 2, lwd = 2 )
    }
    plot( x = range(yrs), y = c( 0, maxCatch ),
          type = "n", axes = F )
      axis( side = 1 )
      axis( side = 2, las = 1)
      if( sIdx == 1 ) mtext(  side = 3, text = "Total", 
                              cex = 0.8, line = 2 )
      lines(  x = specTotalCatch$year,
              y = specTotalCatch$landedWt, lwd = 2 )
      lines( x = specTotalCatch$year, y = specTotalCatch$discardWt, lty = 2, lwd = 2 )
  }

  dev.new()
  # now, plot a grid of species/fleet catch
  par( mfrow = c(length(specCodes), length(fleets) + 1 ), 
       mar = c(1,1,1,1), oma = c(3,3,3,3) )
  for( sIdx in 1:length(specCodes) )
  {
    specCatch <-  df %>%
                  filter( specCode == specCodes[sIdx])
    specTotalCatch <- specCatch %>%
                      group_by(year) %>%
                      summarise(  landedWt = sum( landedWt ),
                                  discardWt = sum( discardWt ) ) %>%
                      ungroup()
    maxCatch <- max(specTotalCatch$landedWt)
    for( fIdx in 1:length(fleets) )
    {
      specFleetCatch <- specCatch %>%
                        filter( gearType == fleets[fIdx] ) %>%
                        group_by(year) %>%
                        summarise( landedWt = sum(landedWt),
                                   discardWt = sum(discardWt) ) %>%
                        ungroup()
      plot( x = range(yrs), y = c( 0, maxCatch ),
            type = "n", axes = F )
        axis( side = 1 )
        axis( side = 2, las = 1)
        if( sIdx == 1 ) mtext(  side = 3, text = fleets[fIdx], 
                                cex = 0.8, line = 2 )
        if( fIdx == 1 ) mtext(  side = 2, las = 0, text = names(specCodes)[sIdx], 
                                cex = 0.8, line = 2 )
        if( nrow(specFleetCatch) == 0 ) next
        lines(  x = specFleetCatch$year,
                y = specFleetCatch$landedWt, lwd = 2 )
        lines( x = specFleetCatch$year, y = specFleetCatch$discardWt, lty = 2, lwd = 2 )
    }
    plot( x = range(yrs), y = c( 0, maxCatch ),
          type = "n", axes = F )
      axis( side = 1 )
      axis( side = 2, las = 1)
      if( sIdx == 1 ) mtext(  side = 3, text = "Total", 
                              cex = 0.8, line = 2 )
      lines(  x = specTotalCatch$year,
              y = specTotalCatch$landedWt, lwd = 2 )
      lines( x = specTotalCatch$year, y = specTotalCatch$discardWt, lty = 2, lwd = 2 )
  }
}


# function to sum all age observations for each year and return proportions
ageProps <- function(year, ageObs )
{
  for(y in unique(year) )
  {
    yIdx <- which(year == y)
    totObs <- sum(ageObs[yIdx])
    ageObs[yIdx] <- ageObs[yIdx]/totObs
  }
  ageObs
}


plotAgeBubbles <- function( bioData, yrs = allYrs )
{
  # ageBubbles()
  # Take biological data for one of the DERPA species and produces
  # an age bubble plot.
  # Todo: separate surveys and colour code.

  # Summarise dover sole ages into frequencies
  ageTbl <-   bioData[!is.na(bioData$AGE),] %>%
              group_by( YEAR, AGE ) %>%
              summarise( nAge = n() ) %>%
              mutate( pAge = nAge/sum(nAge),
                      tot = sum(nAge) )

  yrs <- yrs[1]:yrs[2]
  ages <- unique(ageTbl$AGE)
  ages <- range(ages)

  plot( x = range(yrs), y = range(c(-1,ages)), axes = F,
        xlab = "", ylab = "", type = "n")
    points( x = ageTbl$YEAR, y = ageTbl$AGE,
            cex = 3*sqrt(ageTbl$pAge), pch = 1 )
    text( x = unique(ageTbl$YEAR), y = -.5, 
          labels = unique(ageTbl$tot), cex = 0.7 )
    axis(side = 1, at = yrs )
    mtext( side = 2, text = "Age", line = 4, cex = 0.8 )
    axis(side = 2, las = 1 )
    
    
}

plotRelativeBio <- function(  densData, strataArea = stratArea, 
                              surveyIDs = survIDs, yrs )
{
  # plotRelativeBio()
  # Takes CPUE density estimates and produces weighted average relative biomass
  # estimates for each survey series. Averages are weighted by the survey strata
  # areas. 
  # inputs: densData: Density data from DFO data request
  #         strataArea: strata area with grouping codes for dens data
  #         surveyIDs:  survey ID codes with names for plot legend
  densData2 <-  densData %>%
                left_join( strataArea, by = "GROUPING_CODE") %>%
                group_by(SURVEY_SERIES_ID, SURVEY_ID, YEAR, GROUPING_CODE ) %>%
                dplyr::summarise( density = mean(DENSITY_KGPM2),
                                  area = mean(AREA_KM2) ) %>%
                dplyr::summarize( relBio_Kt = sum( area * density) )

  plot( x = yrs, y = range(densData2$relBio_Kt),
        axes = F, type = "n", xlab = "", ylab = "" )
    points( x = densData2$YEAR, y = densData2$relBio_Kt,
            pch = densData2$SURVEY_SERIES_ID )
    axis( side =1, at = seq(yrs[1],yrs[2], by = 5))
    axis( side = 2, las = 1 )
    mtext( side = 2, text = "Relative Biomass (Kt)", line = 4, cex = 0.8 )
    panLegend( x = 0, y = 1, legTxt = names(surveyIDs), pch = surveyIDs,
                bty = "n" )
}



lengthWt <- function( bioData, plot = T )
{
  # mung data
  # browser()
  bioData <-  bioData %>%
              dplyr::filter(!is.na(LENGTH_MM) & !is.na(WEIGHT_G) ) %>%
              mutate(logL = log(LENGTH_MM), logW = log(WEIGHT_G) ) %>%
              dplyr::select( LENGTH_MM, WEIGHT_G,logL, logW, SEX)
              

  bioDataBoys   <- bioData %>% filter( SEX == 1 )
  bioDataGirls  <- bioData %>% filter( SEX == 2 )
  bioDataBlank  <- bioData %>% filter( SEX == 0 | SEX == 3 )

  # Now fit lw relationship
  lwAll   <- lm( logW ~ logL, data = bioData )
  lwBoys  <- lm( logW ~ logL, data = bioDataBoys )
  lwGirls <- lm( logW ~ logL, data = bioDataGirls )
  lwBlank <- lm( logW ~ logL, data = bioDataBlank )

  # recover coefficients
  lwAll       <- coef(lwAll)
  lwAll[1]    <- exp(lwAll[1])
  lwBoys      <- coef(lwBoys)
  lwBoys[1]   <- exp(lwBoys[1])
  lwGirls     <- coef(lwGirls)
  lwGirls[1]  <- exp(lwGirls[1])
  lwBlank     <- coef(lwBlank)
  lwBlank[1]  <- exp(lwBlank[1])

  # Now plot
  if( plot )
  {
    plotCols <- brewer.pal(n = 4, name = "Dark2")
    lengthRan <- range(bioData$LENGTH_MM)
    lengthSeq <- seq(lengthRan[1],lengthRan[2], length = 100 )
    plot( x = range(bioData$LENGTH_MM),
          y = range(bioData$WEIGHT_G),
          axes = F, xlab = "", ylab = "", type = "n" )
      # DATA
      # Boys

      points( x = bioData$LENGTH_MM, y = bioData$WEIGHT_G,
              pch = bioData$SEX, col = "grey60" )
      # MODEL
      lines( x = lengthSeq, y = lwBoys[1] * lengthSeq^lwBoys[2],
              col = plotCols[1], lwd = 2 )
      lines(  x = lengthSeq, y = lwGirls[1] * lengthSeq^lwGirls[2],
              col = plotCols[2], lwd = 2 )
      lines(  x = lengthSeq, y = lwBlank[1] * lengthSeq^lwBlank[2],
              col = plotCols[3], lwd = 2 )
      lines(  x = lengthSeq, y = lwAll[1] * lengthSeq^lwAll[2],
              col = plotCols[4], lwd = 2 )
      # Legend
      panLegend(  x = 0, y = 0.9, legTxt = c("Boys", "Girls", "Blank", "All"),
                  col = plotCols, lwd = 2, pch = c(1:3,NA), bty = "n" )
      axis( side = 1 )
      axis( side = 2, las = 1)
      mtext( side = 1, text = "Length(mm)", line = 2, cex = 0.8)
      mtext( side = 2, text = "Weight (g)", line = 4, cex = 0.8)
  }

  outList <- list(  aggregate = lwAll, 
                    boys = lwBoys, 
                    girls = lwGirls, 
                    blank = lwBlank )
  outList 
}

# Define a NLL for optim
vonB_nll <- function( theta,
                      data )
{
  # Recover leading pars
  Linf  <- exp(theta[1])
  K     <- exp(theta[2])
  t0    <- theta[3]
  cvL   <- 1 / (1 + exp(-theta[4]))

  data <- data %>%
          mutate( expLt = vonB(AGE,K=K,Linf=Linf,t0=t0),
                  res   = expLt - LENGTH_MM )

  # browser()
  maxLen <- max(data$LENGTH_MM)
  lengths <- unique(data$LENGTH_MM)
  nll <- 0
  for( len in lengths )
  {
    subData <- data %>% filter( LENGTH_MM == len )
    sigL <- cvL * LENGTH_MM
    nll <- nll +  nrow(subData) * log(sigL*sigL) + 0.5*sum(subData$res^2)/sigL/sigL
  }
  
  nll <- nll + (Linf - maxLen)^2/(maxLen/2)

  return(nll)
}

# von bertalanffy growth function
vonB <- function( age, K, Linf, t0, theta = NULL )
{
  if(!is.null(theta))
  {
    Linf <- exp(theta[1])
    K <- exp(theta[2])
    t0 <- theta[3]
  }

  Lt <- Linf * ( 1 - exp( -K * (age - t0)))
  Lt
}

lengthAge <- function( bioData, plot = T, plotLeg = F )
{
  bioData <-  bioData %>%
              dplyr::filter(!is.na(LENGTH_MM) & !is.na(AGE) & LENGTH_MM > 0 ) %>%
              dplyr::select( LENGTH_MM, AGE, SEX)
              

  bioDataBoys   <- bioData %>% filter( SEX == 1 )
  bioDataGirls  <- bioData %>% filter( SEX == 2 )
  # bioDataBlank  <- bioData %>% filter( SEX == 0 | SEX == 3 )

  # Now optimise
  Linf <- max(bioData$LENGTH_MM)
  theta <- c(log(Linf),0,0,0)
  laAll <- optim(par = theta, fn = vonB_nll, data = bioData )
  laBoys <- optim(par = theta, fn = vonB_nll, data = bioDataBoys )
  laGirls <- optim(par = theta, fn = vonB_nll, data = bioDataGirls )
  # laBlank <- optim(par = c(6,0,0), fn = vonB_nll, data = bioDataBlank )

  # Now plot
  if( plot )
  {
    plotCols <- brewer.pal(n = 4, name = "Dark2")
    ageRan <- range(bioData$AGE)
    lengthRan <- range(bioData$LENGTH_MM)
    ageSeq <- seq(1,max(ageRan), length = 100 )
    plot( x = ageRan,
          y = lengthRan,
          axes = F, xlab = "", ylab = "", type = "n" )
      # DATA
      points( x = bioData$AGE, y = bioData$LENGTH_MM,
              pch = bioData$SEX, col = "grey60" )
      # points( x = bioData$LENGTH_MM, y = bioData$WEIGHT_G,
      #         pch = 4, col = "grey40" )
      # MODEL
      # browser()
      lines( x = ageSeq, y = vonB(ageSeq, theta=laBoys$par),
              col = plotCols[1], lwd = 2 )
      lines(  x = ageSeq, y = vonB(ageSeq, theta=laGirls$par),
              col = plotCols[2], lwd = 2 )
      # lines(  x = ageSeq, y = vonB(ageSeq, theta=laBlank$par),
      #         col = plotCols[3], lwd = 2 )
      lines(  x = ageSeq, y = vonB(ageSeq, theta=laAll$par),
              col = plotCols[4], lwd = 2 )
      # Legend
      if(plotLeg)
        panLegend(  x = 0, y = 0.9, legTxt = c("Boys", "Girls", "Blank", "All"),
                    col = plotCols, lwd = 2, pch = c(1:3,NA), bty = "n" )
      axis( side = 1 )
      axis( side = 2, las = 1)
      mtext( side = 1, text = "Age", line = 2, cex = 0.8)
      mtext( side = 2, text = "Length (mm)", line = 4, cex = 0.8)
  }

  outList <- list(  aggregate = laAll, 
                    boys = laBoys, 
                    girls = laGirls ) 
  outList 
}

makeWtMat <- function( bioData, plot = T, maxMatAge = 100 )
{
  # makeWtMat()
  # Creates weight and maturity at age arrays, for use in further
  # analyses, such as solving for maturity schedules and Ford-Walford
  # growth parameters. Need to think about how to use the sparse versions
  # of these that some species have - maybe solve for age and
  # weight from growth models, or use a synthetic cohort?

  # Mung data, 
  bioDataWt <-  bioData %>%
                filter( !is.na(AGE) & !is.na(MAT) & !is.na(WEIGHT_G)
                        & (SEX == 1 | SEX == 2) ) %>%
                group_by( YEAR, AGE, SEX ) %>%
                summarise(  weight = mean(WEIGHT_G), nFish = n() )
  
  pMatFun <- function( MAT, nMat, totAge )
  {
    if( MAT >= 3 )  pMatAge <- nMat / totAge
    if( MAT < 3 )   pMatAge <- 0
    pMatAge
  } 

  bioDatMat <-  bioData %>%
                filter( !is.na(AGE) & !is.na(MAT) & !is.na(WEIGHT_G) &
                        (SEX == 1 | SEX == 2 ) & ( AGE <= maxMatAge ) ) %>%
                group_by( YEAR, AGE, MAT ) %>%
                summarise( nMat = n() ) %>%
                mutate( totAge = sum(nMat),
                        pMatAge = mapply(FUN = pMatFun, MAT, nMat, totAge ) ) %>%
                summarise( pMatAge = sum(pMatAge) )


  yrs   <- range(bioData$YEAR)
  nYrs  <- yrs[2] - yrs[1] + 1
  ages  <- range(bioData$AGE,na.rm=T)
  nAges <- max(ages)

  # Separate weight by sex, to check for sexual dimorphism
  Wtas <- array(  NA, dim = c(nYrs, nAges, 2), 
                  dimnames = list(  yrs[1]:yrs[2], ages[1]:ages[2], 
                                    c("Boys","Girls" ) ) ) 
  for( t in yrs[1]:yrs[2] )
  {
    for( a in ages[1]:ages[2])
    {
      tbl <- bioDataWt %>% filter( YEAR == t, AGE == a)
      if(nrow(tbl) == 0) next
      Wtas[t-yrs[1]+1,a,] <- tbl$weight
    }
  }
  # browser()
  matModel <- glm( pMatAge ~ AGE, family = binomial(logit), data = bioDatMat )

  a50 <- dose.p(matModel,p=c(0.5))
  a95 <- dose.p(matModel,p=c(0.95))

  Mta <- array( NA, dim = c(nYrs, nAges ),
                dimnames = list( yrs[1]:yrs[2], ages[1]:ages[2]  )  ) 
  for( t in yrs[1]:yrs[2] )
  {
    for( a in ages[1]:ages[2])
    {
      tbl <- bioDatMat %>% filter( YEAR == t, AGE == a)
      if(nrow(tbl) == 0) next
      # browser()
      Mta[t-yrs[1]+1,a] <- tbl$pMatAge
    }
  }

  outList <- list ( Wtas = Wtas, Mta = Mta, 
                    matModel = list(  glm = matModel,
                                      a50 = a50, a95 = a95),
                    bioDatMat = bioDatMat )

  outList
}

plotMaturity <- function( matObj )
{
  a50 <- matObj$matModel$a50[1]
  a95 <- matObj$matModel$a95[1]
  bioDatMat <- matObj$bioDatMat

  ages <- c(1,ncol(matObj$Mta))

  # browser()
  logis <- logisticAtAge(ages[1]:ages[2], a50[1], a95[1] )
  plot( x = ages, y = c(0,1), type = "n", xlab = "",
        ylab = "", axes = F )
    text( x = bioDatMat$AGE, y = bioDatMat$pMatAge,
            labels = bioDatMat$YEAR, cex = 1, col = "grey70")
    lines(x = ages[1]:ages[2], y = logis, lwd = 2 )
    axis( side = 1, at = seq(ages[1],ages[2],by = 5) )
    axis( side = 2, las =1 )
    mtext(side = 1, text = "Age", cex = 0.8, line = 2)
    mtext(side = 2, text = "Proportion Mature", line = 4, cex = .8)
}

fitWalford <- function( mat, yrs )
{
  dim <- nrow(mat)
  yrs <- yrs[1]:yrs[2]
  # browser()
  coeffs <- matrix(NA,nrow=dim,ncol=3)
  coeffs[,1] <- yrs
  len    <- ncol(mat)
  for( t in 1:nrow(mat) )
  { 
    nObs <- length(which(!is.na(mat[t,])))
    if( nObs < 2 ) next
    la1 <- mat[t,1:(len-1)]
    la2 <- mat[t,2:len]
    if( t < nrow(mat) )
    {
      tmp <- lm(la2~la1)
      coeffs[t,2:3] <- tmp$coef[1:2]
    }
  }
  return( coeffs )
}

# logistic function with a50 and a95 as inputs
logisticAtAge <- function( a, a50=5, a95=7 )
{
  # solve for step (scale), a50 is location
  s <- (a95 - a50) / log(19)
  # now compute logistic
  p <- 1 / (1 + exp((-a + a50)/s))
  p
}

# Make proportions at age
makePAA <- function( bioData )
{
  # Mung data
  bioDataPAA <-   bioData %>%
                  filter( !is.na(AGE) ) %>%
                  group_by( YEAR, AGE ) %>%
                  summarise( nAge = n() ) %>%
                  mutate( pAge = nAge / sum(nAge) )

  totalPAA  <-    bioData %>%
                  filter( !is.na(AGE) ) %>%
                  group_by( AGE ) %>%
                  summarise( nAge = n() ) %>%
                  mutate( pAge = nAge / sum(nAge) )


  yrs <- range(bioDataPAA$YEAR)
  nYrs <- yrs[2] - yrs[1] + 1
  ages <- range(bioDataPAA$AGE)
  nAges <- max(ages)

  # arrange into an array
  Pta <- array( NA, dim = c(nYrs, nAges), 
                dimnames = list(  yrs[1]:yrs[2], ages[1]:ages[2] ) ) 
  for( t in yrs[1]:yrs[2] )
  {
    for( a in ages[1]:ages[2])
    {
      tbl <- bioDataPAA %>% filter( YEAR == t, AGE == a)
      if(nrow(tbl) == 0) next
      Pta[t-yrs[1]+1,a] <- tbl$pAge
    }
  }

  outList <- list( PAAmat = Pta, totPAA = totalPAA, PAAtbl = bioDataPAA )
}

plotPAA <- function( PAAobj )
{
  totalPAA <- PAAobj$totPAA
  bioDataPAA <- PAAobj$PAAtbl
  maxPAA <- max(bioDataPAA$pAge)
  ages <- range(bioDataPAA$AGE)
  plot( x = ages, y = c(-.05, 1.05*maxPAA + 0.05 ) ,
        axes = F, xlab = "", ylab = "", type = "n")
    yRan <- range(bioDataPAA$pAge)
    yxt <- seq(yRan[1],yRan[2],length=5)
    axis( side = 1, at = seq(ages[1],ages[2],by=5) )
    mtext(side = 1, text = "Age", line = 2)
    mtext(side = 2, text = "Proportion at age", line = 4, cex = .8)
    for(yr in bioDataPAA$YEAR)
    {
      tbl <-  bioDataPAA %>%
              filter(YEAR == yr)
      lines(x = tbl$AGE, y = tbl$pAge + 0.05, lwd = 0.1, col = "grey80")
    }
    
    lines( x = totalPAA$AGE, y = totalPAA$pAge + 0.05 )
    points( x = totalPAA$AGE, y = totalPAA$pAge + 0.05, pch = 16, cex = 0.5 )
    axis( side = 2, las = 1, at = seq(0,maxPAA,length = 3) + 0.05, labels = round(seq(0,maxPAA,length = 3),2) )
    text( x = totalPAA$AGE, y = -.02, labels = totalPAA$nAge, cex = 0.5,
          srt = 90)
}

# ZEstSeber() 
# uses the survivorship estimator defined on page 419 of
# Seber (1982), where 
# inputs:   a0 = age at full recruitment/selectivity 
#           J = max unrolled age class, plusgroup A = a0 + J + 1
#           nAges = vector of age proportions
# outputs:  S = estimate of survivorship
ZEstSeber <- function( nAges, a0, J = 0 )
{
  # roll up n(J) class
  # browser()
  if( J <= 0 ) J <- length( nAges ) - 1
  nJ <- sum(nAges[(J+1):length(nAges)], na.rm = T)
  n <- sum(nAges)
  
  # age class vector
  a <- seq(a0,J-1,1) - a0
  # geometric RV
  X <- sum(a * nAges[a0:(J-1)]) + nJ*J
  # MLE
  S <- X / (n - nJ + X)
  # Mortality
  M <- -1.*log(S)
  M
}

# Turns time-age matrices into cohort matrices for
# ford-walford parameters etc
# Source: SP Cox, plotWt.R
makeSizeAge <- function( Wta, lwRel )
{

  ages      <- 1:ncol(Wta)
  maxWt     <- max(Wta,na.rm=TRUE)
  maxAge    <- max(ages)
  xRange    <- c( 0,maxAge )
  yrs       <- range(as.numeric(rownames(Wta)))
  cohYrs    <- yrs
  cohYrs[1] <- yrs[1] - maxAge + 1

  dim <- cohYrs [ 2 ] - cohYrs [ 1 ] + 1

  c1 <- lwRel[1]
  c2 <- lwRel[2]

  matYr <- matrix(NA, nrow=dim, ncol=dim, byrow=T )
  matWt <- matrix(NA, nrow=dim, ncol=dim, byrow=T )
  matLt <- matrix(NA, nrow=dim, ncol=dim, byrow=T )
  for( cohRow in 1:dim )
  {
    maxCol <- min( cohRow + maxAge - 1, dim )
    a <- 0
    for( j in cohRow:maxCol )
    {
      a <- a + 1
      cohYr <- (cohYrs[1] + cohRow - 1) + a - 1
      matYr[ cohRow, j ] <- cohYr
      if( !( cohYr %in% as.numeric( rownames( Wta ) ) ) ) next
      matRow <- cohYr - yrs[1] + 1
      matWt[ cohRow, j ] <- Wta[ matRow, a ]
      logL               <- (log(Wta[matRow,a])-log(c1))/c2
      matLt[ cohRow, j ] <- exp(logL)
    }
  }
  result <- list()
  result$Yr <- matYr
  result$Wt <- matWt
  result$Lt <- matLt
  return( result )

}  

# makeCohortMat()
# Takes a matrix with time as rows and age as columns
# and produces a pair of cohort matrices with cohorts as rows
# and years as columns. Entries in one matrix are the years, 
# and the other showing the Mta entries
makeCohortMat <- function ( Mta )
{
  ages      <- 1:ncol(Mta)
  maxMt     <- max(Mta,na.rm=TRUE)
  maxAge    <- max(ages)
  xRange    <- c( 0,maxAge )
  yrs       <- range(as.numeric(rownames(Mta)))
  cohYrs    <- yrs
  cohYrs[1] <- yrs[1] - maxAge + 1

  dim <- cohYrs [ 2 ] - cohYrs [ 1 ] + 1

  matYr <- matrix(NA, nrow=dim, ncol=dim, byrow=T )
  matMt <- matrix(NA, nrow=dim, ncol=dim, byrow=T )

  for( cohRow in 1:dim )
  {
    maxCol <- min( cohRow + maxAge - 1, dim )
    a <- 0
    for( j in cohRow:maxCol )
    {
      a <- a + 1
      # browser()
      cohYr <- (cohYrs[1] + cohRow - 1) + a - 1
      matYr[ cohRow, j ] <- cohYr
      if( !( cohYr %in% as.numeric( rownames( Mta ) ) ) ) next
      # browser()
      matRow <- cohYr - yrs[1] + 1
      matMt[ cohRow, j ] <- Mta[ matRow, a ]
    }
  }

  result <- list()
  result$Yr <- matYr
  result$Mt <- matMt
  return( result )
}

# Plots cohort growth.
# Source: SP Cox, plotWt.R
plotMat <- function( matYr, mat, yrs )
{
  ages   <- 1:ncol(mat)
  maxX   <- max(mat,na.rm=TRUE)
  maxAge <- max(ages)
  xRange <- c( 0,maxAge )
   xLim  <- yrs
   yLim  <- c( 0,maxX )
   pCol  <- rep( c("white","white","white","white","black"), 8 )

  # Weight at age.
  plot( xLim,yLim,type="n",axes=FALSE,xlab="",ylab="" ) 
  
  for( t in 1:nrow(mat) )
  {
    lines( matYr[t,], mat[t,], lty=1, lwd=1 )
    points( matYr[t,], mat[t,], bg=pCol[t], col="black", cex=.8, pch=21 )
  }

  axis( side=1, cex.axis=1 )
  axis( side=2, cex.axis=1, las=1 )
} 

dataPlots <- function(  specLab = "dover", stratArea = strata,
                        survIDs = surveyIDs, minYr = 1995, save = F,
                        maxMatAge = 25 )
{
  # Read in density data
  densityFname  <- paste( specLab, "_density.csv", sep = "" )
  densityPath   <- file.path(getwd(),"density",densityFname)
  density       <- read.csv(densityPath, header = T)

  # Read in bio data
  bioFname  <- paste( specLab, "_bio.csv", sep = "" )
  bioPath   <- file.path(getwd(),"biology_with_mat",bioFname)
  bio       <- read.csv(bioPath, header = T)

  # Join desnsity and bio data by trip_ID so we can plot age bubbles etc
  survYearIDs <-  density %>%
                  dplyr::select(  SURVEY_SERIES_ID, 
                                  SURVEY_ID,
                                  YEAR,
                                  TRIP_ID,
                                  SURVEY_DESC )

  bio <-  bio %>% 
          inner_join( survYearIDs ) %>%
          distinct()

  # get length-at-age and length-wt relationships
  lw <- lengthWt( bio,plot = F )
  la <- lengthAge( bio, plot = F)
  

  # weight and maturity at age (from data, not growth model)
  # The following function estimates maturity ogives using a GLM,
  # but it's not very good - should really model by cohort, maybe
  # talk to Michelle
  # Also, check Seber 1982 for maturity

  # make arrays of average weight and proportion mature at time and age
  wtMatArrays <- makeWtMat( bio, plot = F, maxMatAge = maxMatAge )
  # Make proportions at age
  PAA <- makePAA(bio)
  PAAcohorts <- makeCohortMat( PAA$PAAmat)
  

  # make cohort matrices from Wta matrix
  sizeAge_Boys <- makeSizeAge(wtMatArrays$Wtas[,,1], lwRel = lw$boys )
  sizeAge_Girls <- makeSizeAge(wtMatArrays$Wtas[,,2], lwRel = lw$girls )

  # Ford-Walford
  # fw_Boys <- fitWalford(sizeAge_Boys$Lt, yrs = range(sizeAge_Boys$Yr, na.rm=T))
  # fw_Girls <-fitWalford(sizeAge_Girls$Lt, yrs = range(sizeAge_Girls$Yr, na.rm=T))
  
  Z <- ZEstSeber(nAges = PAA$totPAA$nAge, a0 = 9, J=20)

  allYrs <- range(density$YEAR)
  if(!is.null(minYr)) allYrs[1] <- minYr

  par(mar = c(2,2,2,2), oma = c(2,3,2,0))

  plotLayout <- matrix( c(1,1,5,
                          1,1,5,
                          2,2,6,
                          2,2,6,
                          3,3,7,
                          3,3,7,
                          4,4,8,
                          4,4,8), 
                        ncol = 3, byrow = T)

  if(save)
  {
    saveFile <- paste(specLab, ".pdf", sep = "" )
    savePath <- file.path(getwd(),"figs",saveFile)
    pdf(file = savePath )
  }
  
  layout(plotLayout)

  plotRelativeBio(density, yrs = allYrs)
  plotAgeBubbles(bio, yrs = allYrs)
  plotMat(PAAcohorts$Yr, PAAcohorts$Mt, yrs = allYrs)
  mtext(side =2, text = "Proportion at age", line = 4, cex = 0.8 )
  plotMat(sizeAge_Boys$Yr, sizeAge_Boys$Wt, yrs = allYrs)
  mtext(side =2, text = "Weight at age (g)", line = 4, cex = 0.8 )
  mtext(side = 1, text = "Year", outer = F, line = 3)
  lengthWt(bio, plot=T)
  lengthAge(bio,plot=T)
  plotPAA(PAA)
  plotMaturity( wtMatArrays )
  mtext(side = 3, outer = T, text = specLab )

  if(save) dev.off()

  result <- list()
  result$lw <- lw
  result$la <- la
  result$Z  <- Z
  result$PAA <- PAA
  result$sizeAge_Boys <- sizeAge_Boys
  result$sizeAge_Girls <- sizeAge_Girls
  result$bioData <- bio
  result$densityData <- density
  result$weightMaturity <- wtMatArrays

  return(result)
}
