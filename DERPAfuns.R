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

# Load CRS codes
loadCRS <- function()
{
  AEAproj <<- CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m")
  LLproj  <<- CRS("+proj=longlat +datum=WGS84")
  UTMproj <<- CRS("+proj=utm +zone=9 +datum=WGS84")

  invisible(NULL)
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



# Calculate relative biomass by species and arrange in an array
# for feeding to TMB model
makeRelBioStocks <- function( spec = "dover",
                              years = c(1975, 2016), 
                              stocks = list(  HGHS = c(2,3,16),
                                              QCS = c(1),
                                              WCVI = c(4) ),
                              survIDs = surveyIDs,
                              stratArea = stratData   )
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
                                        density = DENSITY_KGPM2,
                                        areaFished_km2 = DOORSPREAD_M * TOW_LENGTH_M )

  includedSurveys <- unlist(stocks)

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
                                fishedArea = areaFished_km2  ) %>%
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

appendName <- function( dataCode, stockCodeKey )
{
  stockCodeKeyVec <- unlist(stockCodeKey)
  nameVec <- c()
  for( k in 1:length(stockCodeKey))
    nameVec <- c(nameVec,rep(names(stockCodeKey)[k],length(stockCodeKey[[k]])))

  stockCodeNum <- which( dataCode == stockCodeKeyVec )

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
                                            stockCodeKey = stocks ) )

  histData  <-  read.csv( file.path(datPath, historicName), header = T,
                          stringsAsFactors = FALSE ) %>%
                mutate( period = "historic" ) %>%
                dplyr::select(  version = "formula_version",
                                est_link, se_link, area, year,
                                period ) %>%
                mutate( stockName = sapply( X = area, 
                                            FUN = appendName, 
                                            stockCodeKey = stocks ) )

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
                          dimnames = list(  c("historic","modern"),
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
    dataArray["historic",stockName,as.character(histYrs),"log.mean","stdized"] <- subHistStd$est_link
    dataArray["historic",stockName,as.character(histYrs),"log.mean","unstdized"] <- subHistUnStd$est_link
    # historic data, standard errors
    dataArray["historic",stockName,as.character(histYrs),"log.sd","stdized"] <- subHistStd$se_link
    dataArray["historic",stockName,as.character(histYrs),"log.sd","unstdized"] <- subHistUnStd$se_link
    # Modern data, mean values
    dataArray["modern",stockName,as.character(modYrs),"log.mean","stdized"] <- subModStd$est_link
    dataArray["modern",stockName,as.character(modYrs),"log.mean","unstdized"] <- subModUnStd$est_link
    # Modern data, standard errors
    dataArray["modern",stockName,as.character(modYrs),"log.sd","stdized"] <- subModStd$se_link
    dataArray["modern",stockName,as.character(modYrs),"log.sd","unstdized"] <- subModUnStd$se_link
  }

  return( list( cpue.df = allData,
                cpue.arr = dataArray ) ) 
}

# Quick function to scale by mean values
transByMean <- function(x)
{
  x <- x -mean(x,na.rm = T)

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
                                            stockCodeKey = stocksSurv),
                        stockName = unlist(stockName) )

  
  # Now do the same for the commercial data
  commBio <-  commBio %>%
              mutate( stockName = sapply( X = MAJ, 
                                          FUN = appendName,
                                          stockCodeKey = stocksComm ) )


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
        oma = c(4,4,1,3) )

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
                 label = c(" ", "All", "Boys", "Girls" ), cex = .8 )
          text(  x = c(0.4,0.55,0.7,0.85)*maxAge, y = 0.2*maxLen,
                 label = c("Linf ", Linf.stock ), cex = .8 )
          text(  x = c(0.4,0.55,0.7,0.85)*maxAge, y = 0.15*maxLen,
                 label = c("vonK ", K.stock ), cex = .8 )
          text(  x = c(0.4,0.55,0.7,0.85)*maxAge, y = 0.1*maxLen,
                 label = c("L1 ", L1.stock ), cex = .8 )
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
  # Add legend to bottom of plot
  par(mfcol = c(1,1) )
  legend( x = "bottom",
          horiz = TRUE,
          legend = c("Male", "Female", "Unsexed", "Coastwide"),
          pch = c( 1, 2, 3, NA ),
          lty = c( 1, 1, 1, 2 ),
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
        oma = c(4,4,1,3) )

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
        if( mfg[1] == mfg[3] )
          axis( side = 1 )
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
  # Add a legend
  par( mfcol = c(1,1) )
  legend( x = "bottom",
          horiz = TRUE,
          legend = c("Male", "Female", "Unsexed", "Coastwide"),
          pch = c( 1, 2, 3, NA ),
          lty = c( 1, 1, 1, 2 ),
          col = c( sexCols[1:2], "black", "grey70" ),
          bty = "n" )

  if( save )
  {
    dev.off()
    cat( "Weight-at-length plots saved to ", savePath, "\n", sep = "")
  }
}


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
  ageFreq <- array(NA,  dim = c(nStocks,nGears,nYears,maxAge,3),
                        dimnames = list(  stockNames,
                                          gearNames,
                                          as.character(years),
                                          as.character(1:maxAge),
                                          c("all", "boys", "girls" ) ) ) 

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

        if( nrow(gearData) == 0 )
          next

        # Summarise all data
        allData <-  gearData %>%
                    group_by( year, age ) %>%
                    summarise( nObs = n() ) %>%
                    ungroup()
        # Boys
        boyData <-  gearData %>%
                    filter( sex == 1 ) %>%
                    group_by( year, age ) %>%
                    summarise( nObs = n() ) %>%
                    ungroup()
        # Girls
        girlData <- gearData %>%
                    filter( sex == 2 ) %>%
                    group_by( year, age ) %>%
                    summarise( nObs = n() ) %>%
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
            ageFreq[stockID, gearID, yrChar, ,"all" ] <- 0  
            # Get vector of ages with positive observations
            obsAges <- as.character(yrData.all$age)
            # Save into ageFreq
            ageFreq[stockID, gearID, yrChar, obsAges, "all" ]<- yrData.all$nObs
          }
          
          # Subset to that year's data
          yrData.boy <- boyData %>% filter( year == yr )
          # Fill boy fish slice
          if( nrow(yrData.boy) > 0)
          {
            # Replace NAs with 0s, as we have some observations
            ageFreq[stockID, gearID, yrChar, ,"boys" ] <- 0  
            # Get vector of ages with positive observations
            obsAges <- as.character(yrData.boy$age)
            # Save into ageFreq
            ageFreq[stockID, gearID, yrChar, obsAges, "boys" ]<- yrData.boy$nObs
          }

          # Subset to that year's data
          yrData.girl <- girlData %>% filter( year == yr )
          # Fill girl fish slice
          if( nrow(yrData.girl) > 0)
          {
            # Replace NAs with 0s, as we have some observations
            ageFreq[stockID, gearID, yrChar, ,"girls" ] <- 0  
            # Get vector of ages with positive observations
            obsAges <- as.character(yrData.girl$age)
            # Save into ageFreq
            ageFreq[stockID, gearID, yrChar, obsAges, "girls" ]<- yrData.girl$nObs
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

      if( nrow(gearData) == 0 )
          next

      # Summarise combined sexes
      allData <-  gearData %>%
                  group_by( year, age ) %>%
                  summarise( nObs = n() ) %>%
                  ungroup()

      # Boys
      boyData <-  gearData %>%
                  filter( sex == 1 ) %>%
                  group_by( year, age ) %>%
                  summarise( nObs = n() ) %>%
                  ungroup()
      # Girls
      girlData <- gearData %>%
                  filter( sex == 2 ) %>%
                  group_by( year, age ) %>%
                  summarise( nObs = n() ) %>%
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
          ageFreq[stockID, gearID, yrChar, ,"all" ] <- 0  
          # Get vector of ages with positive observations
          obsAges <- as.character(yrData.all$age)
          # Save into ageFreq
          ageFreq[stockID, gearID, yrChar, obsAges, "all" ]<- yrData.all$nObs
        }
        
        # Subset to that year's data
        yrData.boy <- boyData %>% filter( year == yr )
        # Fill boy fish slice
        if( nrow(yrData.boy) > 0)
        {
          # Replace NAs with 0s, as we have some observations
          ageFreq[stockID, gearID, yrChar, ,"boys" ] <- 0  
          # Get vector of ages with positive observations
          obsAges <- as.character(yrData.boy$age)
          # Save into ageFreq
          ageFreq[stockID, gearID, yrChar, obsAges, "boys" ]<- yrData.boy$nObs
        }

        # Subset to that year's data
        yrData.girl <- girlData %>% filter( year == yr )
        # Fill girl fish slice
        if( nrow(yrData.girl) > 0)
        {
          # Replace NAs with 0s, as we have some observations
          ageFreq[stockID, gearID, yrChar, ,"girls" ] <- 0  
          # Get vector of ages with positive observations
          obsAges <- as.character(yrData.girl$age)
          # Save into ageFreq
          ageFreq[stockID, gearID, yrChar, obsAges, "girls" ]<- yrData.girl$nObs
        }   

      }

    }
  }

  # Return age frequencies
  return(ageFreq)
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
    ageFrq <- comps[[specIdx]]

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

      # Now loop over gears
      for( gearIdx in 1:nGears )
      {
        # Get gear name
        gearName <- gearNames[gearIdx]
        # Skip if there are no observations
        if( all(is.na(ageFrq[stockIdx,gearIdx,,,])) )
          next

        # Here I want to refactor and create 2 functions:
        # 1. Bubble plots
        if( save )
        {
          fileName <- paste( prefix,"Bubbles",gearName,".png", sep = "" )
          png(  file = file.path(stockDir,fileName),
                width = 8.5, height = 11, res = 300,
                units = "in" )
        }
        # Pull subset
        subFrq <- ageFrq[ stockIdx, gearIdx,,,]

        # Set up multipanel region
        par(mfrow = c(3,1), oma = c(4,4,3,2), mar = c(1,1,1,1) )
        # Plot boys
        .plotBubbles( z = subFrq[,,1])
        mtext( side = 4, text = "Male", line = 2 )
        # Plot girls
        .plotBubbles( z = subFrq[,,2])
        mtext( side = 4, text = "Female", line = 2 )
        # Plot combined
        .plotBubbles( z = subFrq[,,3])
        mtext(  side = 4, text = "Both + Unsexed", line = 2 )
        # Add title
        mtext(  side = 3, text = paste(specID, stockID, gearName, sep = " - " ),
                outer = T, line = 1, cex = 2, font = 2 )

        if( save )
          dev.off()

        if(!save)
          dev.new()

        # 2. age freqency plots

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
        
        .plotAgeFreq( z = subFrq[,,1] )
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
        
        .plotAgeFreq( z = subFrq[,,2] )
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
        
        .plotAgeFreq( z = subFrq[,,3] )
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
                          size = .1 )
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
    axis( side = 1 )
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
    text( x = years[zSum > 0], y = nAges + 1, labels = zSum[zSum > 0],
          srt = 45, cex = .5 )

} # END .plotBubbles()

# Hidden function to plot age frequencies. Shamelessly
# copied from ARK and SPC functions in SableOpMod.R,
# with simplifications for our purpose.
.plotAgeFreq <- function( z,
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
    if ( nPanels <= 9 )
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
  yLim <- c(0, 1.2*max(zz, na.rm = T) )

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
} # END .plotAgeFreq()


# makeAgeComps()
# Takes output of readBioData() for a species
# and generates an array of age observation
# frequencies.
makeLenComps <- function( data = bioData$Dover,
                          stocksComm = stocksCommBio,
                          stocksSurv = stocksSurvey,
                          survIDs = surveyIDs,
                          years = 1954:2018  )
{
  # Pull survey and commercial data
  survData <- data$survey %>%
              mutate( survID = sapply(  X = SURVEY_SERIES_ID,
                                        FUN = appendName,
                                        survIDs),
                      length = round(LENGTH_MM/10) )
  commData <- data$comm %>%
              mutate( length = round(LENGTH_MM/10) )

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
  lenFreq <- array(NA,  dim = c(nStocks,nGears,nYears,maxLen,3),
                        dimnames = list(  stockNames,
                                          gearNames,
                                          as.character(years),
                                          as.character(1:maxLen),
                                          c("all", "boys", "girls" ) ) ) 

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
                    filter( !is.na(length), 
                            stockName == stockID ) %>%
                    dplyr::select(  length,
                                    year = YEAR,
                                    sex = SEX ) 

        if( nrow(gearData) == 0 )
          next

        # Summarise all data
        allData <-  gearData %>%
                    group_by( year, length ) %>%
                    summarise( nObs = n() ) %>%
                    ungroup()
        # Boys
        boyData <-  gearData %>%
                    filter( sex == 1 ) %>%
                    group_by( year, length ) %>%
                    summarise( nObs = n() ) %>%
                    ungroup()
        # Girls
        girlData <- gearData %>%
                    filter( sex == 2 ) %>%
                    group_by( year, length ) %>%
                    summarise( nObs = n() ) %>%
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
            lenFreq[stockID, gearID, yrChar, ,"all" ] <- 0  
            # Get vector of ages with positive observations
            obsLens <- as.character(yrData.all$length)
            # Save into lenFreq
            lenFreq[stockID, gearID, yrChar, obsLens, "all" ]<- yrData.all$nObs
          }
          
          # Subset to that year's data
          yrData.boy <- boyData %>% filter( year == yr )
          # Fill boy fish slice
          if( nrow(yrData.boy) > 0)
          {
            # Replace NAs with 0s, as we have some observations
            lenFreq[stockID, gearID, yrChar, ,"boys" ] <- 0  
            # Get vector of ages with positive observations
            obsLens <- as.character(yrData.boy$length)
            # Save into lenFreq
            lenFreq[stockID, gearID, yrChar, obsLens, "boys" ]<- yrData.boy$nObs
          }

          # Subset to that year's data
          yrData.girl <- girlData %>% filter( year == yr )
          # Fill girl fish slice
          if( nrow(yrData.girl) > 0)
          {
            # Replace NAs with 0s, as we have some observations
            lenFreq[stockID, gearID, yrChar, ,"girls" ] <- 0  
            # Get vector of ages with positive observations
            obsLens <- as.character(yrData.girl$length)
            # Save into lenFreq
            lenFreq[stockID, gearID, yrChar, obsLens, "girls" ]<- yrData.girl$nObs
          }   

        }
        # Go to next gearID
        next
      }

      gearData <- survData %>%
                  filter( !is.na(AGE), 
                          stockName == stockID ) %>%
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
                  summarise( nObs = n() ) %>%
                  ungroup()

      # Boys
      boyData <-  gearData %>%
                  filter( sex == 1 ) %>%
                  group_by( year, length ) %>%
                  summarise( nObs = n() ) %>%
                  ungroup()
      # Girls
      girlData <- gearData %>%
                  filter( sex == 2 ) %>%
                  group_by( year, length ) %>%
                  summarise( nObs = n() ) %>%
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
          lenFreq[stockID, gearID, yrChar, ,"all" ] <- 0  
          # Get vector of ages with positive observations
          obsLens <- as.character(yrData.all$length)
          # Save into lenFreq
          lenFreq[stockID, gearID, yrChar, obsLens, "all" ]<- yrData.all$nObs
        }
        
        # Subset to that year's data
        yrData.boy <- boyData %>% filter( year == yr )
        # Fill boy fish slice
        if( nrow(yrData.boy) > 0)
        {
          # Replace NAs with 0s, as we have some observations
          lenFreq[stockID, gearID, yrChar, ,"boys" ] <- 0  
          # Get vector of ages with positive observations
          obsLens <- as.character(yrData.boy$length)
          # Save into lenFreq
          lenFreq[stockID, gearID, yrChar, obsLens, "boys" ]<- yrData.boy$nObs
        }

        # Subset to that year's data
        yrData.girl <- girlData %>% filter( year == yr )
        # Fill girl fish slice
        if( nrow(yrData.girl) > 0)
        {
          # Replace NAs with 0s, as we have some observations
          lenFreq[stockID, gearID, yrChar, ,"girls" ] <- 0  
          # Get vector of ages with positive observations
          obsLens <- as.character(yrData.girl$length)
          # Save into lenFreq
          lenFreq[stockID, gearID, yrChar, obsLens, "girls" ]<- yrData.girl$nObs
        }   

      }

    }
  }

  # Return age frequencies
  return(lenFreq)
} # END makeLenComps()

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

  appendName <- function( majorAreaCode, stockList )
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
                            mutate( stockName = appendName( majorStatArea, stocks) ) %>%
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
