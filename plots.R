# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
#
# plots.R
#
# Plotting functions for hierSCAL TMB model.
# 
# Plots to add:
#     2. Age/length residual bubbles - use plotBubble functions
#     3. Pred/obs catch
#
#
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

# plotSRsp()
# Plots stock-recruit curve for all stocks and species. 
# Wraps plotSR(), and overwrites that functions margins.
plotSRsp <- function( repObj = repInit,
                      initYear = fYear )
{
  # Pull model dimensions
  nS <- repObj$nS
  nP <- repObj$nP

  specNames <- dimnames(repObj$R0_sp)[[1]]
  stockNames <- dimnames(repObj$R0_sp)[[2]]
  
  # Set up plotting environment
  par(mfcol = c(nP, nS), mar = c(2,2,0,0), oma = c(4,4,3,3) )

  for( sIdx in 1:nS )
    for( pIdx in 1:nP )
    {
      plotSR( repObj = repObj,
              sIdx = sIdx, pIdx = pIdx,
              nopar = TRUE )
      # Detect where we are in the mfg
      mfg <- par("mfg")
      # x axis in bottom row
      axis( side = 1)
      # y axis on all plots
      axis( side = 2, las = 1)
      box()

      # Add species names
      if( mfg[1] == 1 )
        mtext( side = 3, text = specNames[sIdx] )
      if( mfg[2] == mfg[4] )
        mtext( side = 4, text = stockNames[pIdx] )
    } 
  mtext(  side = 1, text = "Spawning Stock Biomass (kt)", outer = TRUE, line = 2 )
  mtext(  side = 2, text = "Recruitment (1e6)", outer = TRUE,
          line = 2 )

} # END plotSBspt()



# plotSR()
# Plot stock-recruitment curve for a given species/stock
plotSR <- function( repObj = repInit,
                    sIdx = 1, pIdx = 1,
                    nopar = FALSE )
{
  # Pull stuff from repObj
  Rt      <- repObj$R_spt[sIdx,pIdx,]
  Bt      <- repObj$SB_spt[sIdx,pIdx,]
  R0      <- repObj$R0_sp[sIdx,pIdx]
  B0      <- repObj$B0_sp[sIdx,pIdx]
  rec.a   <- repObj$reca_sp[sIdx,pIdx]
  rec.b   <- repObj$recb_sp[sIdx,pIdx]
  h       <- round(repObj$h_sp[sIdx,pIdx],2)
  phi     <- repObj$phi_sp[sIdx,pIdx]
  
  # Get number of time steps
  nT      <- repObj$nT

  SB <- seq(0,1.2*B0,length = 1000 )
  R  <- rec.a * SB / (1 + rec.b*SB)

  B20 <- 0.2*B0
  R20 <- rec.a * B20 / (1 + rec.b*B20)

  plot( x = range(SB), y = range(R,Rt,na.rm = T ),
        type = "n", las = 1, xlab = "",
        ylab = "", axes = FALSE )
    if(!nopar)
    {
      mtext( side = 1, text = "Spawning Biomass (kt)", line = 2.5)
      mtext( side = 2, text = "Recruits (1e6)", line = 2.5)
      axis( side = 1 )
      axis( side = 2, las = 1)
      box()
    }
    lines( x = SB, y = R, lwd = 3 )
    points( x = Bt[1:nT], y = Rt[2:(nT+1)], pch = 16, col = "grey60" )
    # Plot B0,R0
    segments( x0 = B0, x1 = B0, y0 = 0, y1 = R0,
              lty = 2, lwd = 2 )
    segments( x0 = 0, x1 = B0, y0 = R0, y1 = R0,
              lty = 2, lwd = 2 )
    # Plot B20,R20
    segments( x0 = B20, x1 = B20, y0 = 0, y1 = R20,
              lty = 2, lwd = 2 )
    segments( x0 = 0, x1 = B20, y0 = R20, y1 = R20,
              lty = 2, lwd = 2 )
    # Label with steepness
    panLab( x = 0.8, y = 0.95, txt = paste("h = ", h, sep = "") )

} # END plotSR()

# plotFspft()
plotFspft <- function( repObj = repInit,
                        initYear = fYear )
{
  # get observed and precdicted catch
  F_spft      <- repObj$F_spft
  specNames   <- dimnames(F_spft)[[1]]
  stockNames  <- dimnames(F_spft)[[2]]
  gearNames   <- dimnames(F_spft)[[3]]

  # species/stocks
  nS <- repObj$nS
  nP <- repObj$nP
  nT <- repObj$nT
  nF <- repObj$nF

  gearCols    <- brewer.pal( nF, "Dark2")

  # Create years vector
  years     <- seq(from = initYear, length = nT + 1, by = 1)
  vertLines <- seq(from = initYear, to = max(years), by = 10)

  # Set up plotting environment
  par(mfcol = c(nP, nS), mar = c(0,1,0,1), oma = c(5,4,3,3) )

  for( sIdx in 1:nS )
    for( pIdx in 1:nP )
    {
      
      plot( x = range(years), y = c(0,max(F_spft[sIdx,pIdx,,])),
            type = "n", axes = FALSE )
      # Detect where we are in the mfg
      mfg <- par("mfg")
      # x axis in bottom row
      if( mfg[1] == mfg[3])
        axis( side = 1)
      # y axis on all plots
      axis( side = 2, las = 1)
      box()
      for( fIdx in 1:nF)
      {
        points( x = years[1:nT], y = F_spft[sIdx,pIdx,fIdx,], 
                pch = 16, col = gearCols[fIdx])
        lines(  x = years[1:nT], y = F_spft[sIdx,pIdx,fIdx,], 
                lwd = 2, lty = 1, col = gearCols[fIdx] )
      }

      # Add species names
      if( mfg[1] == 1 )
        mtext( side = 3, text = specNames[sIdx] )
      if( mfg[2] == mfg[4] )
        mtext( side = 4, text = stockNames[pIdx] )
    } 
  mtext(  side = 1, text = "Year", outer = TRUE, line = 2 )
  mtext(  side = 2, text = "Fishing Mortality (/yr)", outer = TRUE,
          line = 2 )
  par( oma =c(0,1,1,1) )
  legend( x = "bottom", horiz =TRUE, bty = "n",
          legend = gearNames,
          pch = c(16),
          col = gearCols,
          lty = c(1), lwd = c(2) )
} # END plotCatchFit_spt


# plotIspft()
plotIspft <- function(  repObj = repOpt,
                        fYear = fYear, lYear = lYear,
                        sIdx = 1, pIdx = 1:3,
                        fIdx = 1:6 )
{
  
  # get observed and precdicted catch
  I_spft      <- repObj$I_spft[sIdx,pIdx,fIdx,,drop = FALSE]
  I_spft[I_spft < 0 ] <- NA 
  specNames   <- dimnames(I_spft)[[1]]
  stockNames  <- dimnames(I_spft)[[2]]
  gearNames   <- dimnames(I_spft)[[3]]

  # species/stocks
  nS <- repObj$nS
  nP <- repObj$nP
  nT <- repObj$nT
  nF <- repObj$nF

  meltI <- melt(I_spft) %>%
            rename( index = value ) %>%
            group_by( species, stocks, fleets ) %>%
            mutate( index = index / mean(index,na.rm = T) ) %>%
            ungroup()

  tmpPlot <-  ggplot( data = meltI, aes(x = years, y = index, col = fleets ) ) +
              facet_grid( stocks ~ species, scale = "fixed" ) +
              geom_point() +
              geom_line() + 
              theme_sleek()

  print(tmpPlot)
  
} # END plotCatchFit_spt

# plotCatchFit()
plotCatchFit_spt <- function( repObj = repInit,
                              initYear = fYear )
{
  # get observed and precdicted catch
  C_spft      <- repObj$C_spft
  predC_spft  <- repObj$predCw_spft

  
  specNames   <- dimnames(C_spft)[[1]]
  stockNames  <- dimnames(C_spft)[[2]]


  # Sum catch over fleets
  C_spt     <- apply(X = C_spft, FUN = sum, MARGIN = c(1,2,4) )
  predC_spt <- apply(X = predC_spft, FUN = sum, MARGIN = c(1,2,4) )

  # species/stocks
  nS <- repObj$nS
  nP <- repObj$nP
  nT <- repObj$nT

  # Create years vector
  years     <- seq(from = initYear, length = nT + 1, by = 1)
  vertLines <- seq(from = initYear, to = max(years), by = 10)

  # Set up plotting environment
  par(mfcol = c(nP, nS), mar = c(0,1,0,1), oma = c(5,4,3,3) )

  for( sIdx in 1:nS )
    for( pIdx in 1:nP )
    {
      
      plot( x = range(years), y = c(0,max(C_spt[sIdx,pIdx,],predC_spt[sIdx,pIdx,])),
            type = "n", axes = FALSE )
      # Detect where we are in the mfg
      mfg <- par("mfg")
      # x axis in bottom row
      if( mfg[1] == mfg[3])
        axis( side = 1)
      # y axis on all plots
      axis( side = 2, las = 1)
      box()
      points( x = years[1:nT], y = C_spt[sIdx,pIdx,], pch = 21, cex = 1.2,
              bg = "white", lwd = 1.5, col = "grey20" )
      lines( x = years[1:nT], y = predC_spt[sIdx,pIdx,], lwd = 2, lty = 2, col = "steelblue" )
      points( x = years[1:nT], y = predC_spt[sIdx,pIdx,], pch = 16, col = "steelblue"  )

      # Add species names
      if( mfg[1] == 1 )
        mtext( side = 3, text = specNames[sIdx] )
      if( mfg[2] == mfg[4] )
        mtext( side = 4, text = stockNames[pIdx] )
    } 
  mtext(  side = 1, text = "Year", outer = TRUE, line = 2 )
  mtext(  side = 2, text = "Catch (kt)", outer = TRUE,
          line = 2 )
  par( oma =c(0,1,1,1) )
  legend( x = "bottom", horiz =TRUE, bty = "n",
          legend = c("Observed","Predicted"),
          pch = c(21,16),
          col = c("grey20","steelblue"),
          lty = c(NA,2), lwd = c(NA,2),
          bg = c("white",NA) )
} # END plotCatchFit_spt

# plotCatchFit_ft()
plotCatchFit_ft <- function(  repObj = repInit,
                              initYear = fYear,
                              sIdx = 1, pIdx = 1 )
{
  # get observed and precdicted catch
  C_ft      <- repObj$C_spft[sIdx,pIdx,,]
  predC_ft  <- repObj$predCw_spft[sIdx,pIdx,,]

  gearNames <- dimnames(C_ft)[[1]]

  # species/stocks
  nT <- repObj$nT
  nF <- repObj$nF

  # Create years vector
  years     <- seq(from = initYear, length = nT + 1, by = 1)
  vertLines <- seq(from = initYear, to = max(years), by = 10)

  # Now count gears with any positive catch (should just be comm right now)
  posCatchGears <- c()
  for( fIdx in 1:nF )
    if( any(C_ft[fIdx,] > 0) )
      posCatchGears <- c(posCatchGears,fIdx)
    
  nPos <- length(posCatchGears)

  gearCols  <- brewer.pal(nF, "Dark2")

  # Set up plotting environment
  par(mfcol = c(nPos, 1), mar = c(0,1,0,1), oma = c(5,4,3,3) )

  for( fIdx in posCatchGears)
  {
    plot( x = range(years), y = c(0,max(C_ft[fIdx,],predC_ft[fIdx,])),
            type = "n", axes = FALSE )
      # Detect where we are in the mfg
      mfg <- par("mfg")
      # x axis in bottom row
      if( mfg[1] == mfg[3])
        axis( side = 1)
      # y axis on all plots
      axis( side = 2, las = 1)
      box()
      points( x = years[1:nT], y = C_ft[fIdx,], pch = 21, cex = 1.2,
              bg = "white", lwd = 1.5, col = "grey20" )
      lines( x = years[1:nT], y = predC_ft[fIdx,], lwd = 2, lty = 2, col = gearCols[fIdx] )
      points( x = years[1:nT], y = predC_ft[fIdx,], pch = 16, col = gearCols[fIdx]  )

      mtext( side = 4, text = gearNames[fIdx], line = 2)
  }
  mtext(  side = 1, outer = T, text = "Year", line = 2 )
  mtext(  side = 2, outer = T, text = "Catch (kt)", 
          line = 2 )
  par( oma = c(1,1,1,1) )

} # END plotCatchFit_ft()

# Plot comp fits
plotCompFitYrs <- function( repObj = repInit,
                            initYear = fYear,
                            sIdx = 1, pIdx = 1,
                            sex = "female",
                            comps = "age",
                            save = FALSE,
                            savePath = "plotFitYrs" )
{
  # Pull predicted and observed ages
  if( comps == "age" )
  {
    max       <- repObj$A_s[sIdx]
    pred_xft <- repObj$aDist_aspftx_hat[1:max,sIdx,pIdx,,,sex]
    obs_xft  <- repObj$age_aspftx[1:max,sIdx,pIdx,,,sex]  
    xLab      <- "Age"
    minProp   <- repObj$minAgeProp
  }
  if( comps == "length" )
  {
    max       <- repObj$L_s[sIdx]  
    pred_xft  <- repObj$lDist_lspftx_hat[1:max,sIdx,pIdx,,,sex]
    obs_xft   <- repObj$len_lspftx[1:max,sIdx,pIdx,,,sex] 
    xLab      <- "Length"
    minProp   <- repObj$minLenProp
  }

  dimNames  <- dimnames(pred_xft)
  compNames <- dimNames[[1]]
  gearNames <- dimNames[[2]]
  yearNames <- dimNames[[3]]

  # Pull model dims
  nF      <- repObj$nF
  nT      <- repObj$nT  

  # Make colours vector
  cols    <- brewer.pal( n = nF, "Dark2" )

  # Make years vector
  years   <- seq(initYear, length = nT+1, by = 1)

  # Now, we want to loop over gear types now, 
  # and record the gIdxes for which there are
  # observations
  obsGears <- c()
  gearTimes <- vector(mode = "list", length = nF)
  for( fIdx in 1:nF )
  {
    if( any(obs_xft[1,fIdx,] >= 0) )
    {
      obsGears <- c(obsGears,fIdx)
      gearTimes[[fIdx]] <- which(obs_xft[1,fIdx,] >= 0)
    }
  }


  # ok, obsGears are the ones we want to plot,
  # and gearTimes is the list of time indices
  for( fIdx in obsGears )
  {
    if(!save)
      dev.new()

    if(save)
    {
      gearPath <- paste(savePath,gearNames[fIdx],".png",sep = "")
      png(  gearPath, 
            width = 11, height = 8.5,
            units = "in", res = 300 )
    }

    times <- gearTimes[[fIdx]]
    # Count the number of age observations
    # there are, and make the plotting window
    nObs <- length(times)
    nCols <- round(sqrt(nObs))
    nRows <- ceiling(nObs/nCols)

    par(  mfcol = c(nRows,nCols), 
          mar = c(1,1,1,1),
          oma = c(3,3,3,3) )

    for( tIdx in times )
    { 
      # get age obs and preds
      nObs           <- sum(obs_xft[,fIdx,tIdx])
      compObsProp_x  <- obs_xft[,fIdx,tIdx]/sum(obs_xft[,fIdx,tIdx])
      compPred_x     <- pred_xft[,fIdx,tIdx]

      plot( x = c(1,max), y = c(0,max(compObsProp_x,compPred_x,na.rm = T) ),
            xlab = "", ylab = "", type = "n", las = 1 )
        rect( xleft = 1:max - .3, xright = 1:max + .3,
              ybottom = 0, ytop = compObsProp_x,
              col = "grey40", border = NA )
        abline( h = minProp, lty = 3, lwd = .8, col = "grey30" )
        lines(  x = 1:max, y = compPred_x, lwd = 1,
                col = cols[fIdx] )
        points(  x = 1:max, y = compPred_x,
                col = cols[fIdx], pch = 21 )
        panLab( x=.5, y = .95, txt = years[tIdx] )
        panLab( x=.1, y = .95, txt = paste("N = ", nObs, sep = "") )

    }
    mtext( side = 3, outer = T, text = gearNames[fIdx], line = 2, font = 2)
    mtext(  side = 1, outer = T, text = xLab, line = 2 )
    mtext(  side = 2, outer = T, text = paste("Proportion-at-",xLab,sep=""), 
            line = 2 )
    if(save)
      dev.off()
  }

} # END plotCompFitYrs()

# Plot comp fits averaged over time (good for early diagnostics)
plotCompFitAvg <- function( repObj = repInit,
                            initYear = fYear,
                            sIdx = 1, pIdx = 1,
                            sex = "female",
                            comps = "age" )
{
  # Pull predicted and observed ages
  if( comps == "age" )
  {
    max       <- repObj$A_s[sIdx]
    pred_xft  <- repObj$aDist_aspftx_hat[1:max,sIdx,pIdx,,,sex]
    obs_xft   <- repObj$age_aspftx[1:max,sIdx,pIdx,,,sex]  
    xLab      <- "Age"
    minProp   <- repObj$minAgeProp
  }
  if( comps == "length" )
  {
    max       <- repObj$L_s[sIdx]  
    pred_xft  <- repObj$lDist_lspftx_hat[1:max,sIdx,pIdx,,,sex]
    obs_xft   <- repObj$len_lspftx[1:max,sIdx,pIdx,,,sex] 
    xLab      <- "Length"
    minProp   <- repObj$minLenProp
  }

  dimNames  <- dimnames(pred_xft)
  compNames <- dimNames[[1]]
  gearNames <- dimNames[[2]]
  yearNames <- dimNames[[3]]


  # Pull model dims
  nF      <- repObj$nF
  nT      <- repObj$nT  

  # Make colours vector
  cols    <- brewer.pal( n = nF, "Dark2" )

  # Make years vector
  years   <- seq(initYear, length = nT+1, by = 1)

  # Now, we want to loop over gear types now, 
  # and record the gIdxes for which there are
  # observations
  obsGears <- c()
  gearTimes <- vector(mode = "list", length = nF)
  for( fIdx in 1:nF )
  {
    if( any(obs_xft[1,fIdx,] >= 0) )
    {
      obsGears <- c(obsGears,fIdx)
      gearTimes[[fIdx]] <- which(obs_xft[1,fIdx,] >= 0)
    }
  }

  # Some stocks don't have age observations,
  # so skip it
  if(length(obsGears) == 0)
    return()

  par(  mfcol = c(length(obsGears),1), 
        mar = c(1,1,1,1),
        oma = c(3,3,3,3) )


  # ok, obsGears are the ones we want to plot,
  # and gearTimes is the list of time indices
  for( fIdx in obsGears )
  {
    times <- gearTimes[[fIdx]]
    
    # Average the observations and
    fleetObs_xt  <- obs_xft[,fIdx,times,drop = FALSE]
    fleetPred_xt <- pred_xft[,fIdx,times,drop = FALSE]
    
    # Average observations and predictions over time
    compObs_x     <- apply( X = fleetObs_xt, FUN = sum, MARGIN = 1 )
    compObsProp   <- compObs_x / sum(compObs_x)
    compPred_x    <- apply( X = fleetPred_xt, FUN = mean, MARGIN = 1 )

    plot( x = c(1,max), y = c(0,max(compObsProp,compPred_x,na.rm = T) ),
            xlab = "", ylab = "", type = "n", las = 1 )
      rect( xleft = 1:max - .3, xright = 1:max + .3,
            ybottom = 0, ytop = compObsProp,
            col = "grey40", border = NA )
      abline( h = minProp, lty = 3, lwd = .8, col = "grey30" )
      lines(  x = 1:max, y = compPred_x, lwd = 1,
              col = cols[fIdx] )
      points(  x = 1:max, y = compPred_x,
              col = cols[fIdx], pch = 21 )

      mtext( side = 4, text = gearNames[fIdx], line = 2)
  }
  mtext(  side = 1, outer = T, text = xLab, line = 2 )
  mtext(  side = 2, outer = T, text = paste("Proportion-at-",xLab,sep=""), 
          line = 2 )

} # END plotCompFitAvg()

# plotRtDev()
# Plot recruitment deviations over time for a given
# species/stock combo. Used individually and wrapped
# in plotRtDevspt()
plotRtDev <- function(  repObj = repInit,
                        initYear = fYear,
                        sIdx = 1, pIdx = 1,
                        nopar = FALSE )
{
  # Pull model dimensions
  nS <- repObj$nS
  nP <- repObj$nP
  nT <- repObj$nT

  # Pull spawning biomass
  omegaR_t    <- repObj$omegaR_spt[sIdx, pIdx, ]

  # Create years vector
  years     <- seq(from = initYear, length = nT + 1, by = 1)
  vertLines <- seq(from = initYear, to = max(years), by = 10)

  if( !nopar )
    par(mfcol = c(1, 1), mar = c(2,2,1,1), oma = c(1,1,1,1) )

  plot( x = range(years), y = range(omegaR_t, na.rm =T),
        type = "n", xlab = "", ylab = "",
        las = 1, axes = FALSE )
    # Plot axes if not wrapped
    if(!nopar)
    {
      axis(side = 2, las =1 )
      axis(side = 1)
      mtext(side = 2, text = "Recruitment Residuals", line = 3)
    }
    box()
    # Plot recruitment
    abline( v = vertLines, lwd = .8, lty = 3, col = "grey80")
    abline( h = 0, lty = 2, lwd = 1, col = "grey50")
    lines( x = years[1:nT], y = omegaR_t[1:nT], lwd = 2, col = "grey30" )
    points( x = years[1:nT], y = omegaR_t[1:nT], pch = 21, bg = "white" )
    
} # END plotRtDev()

# plotRDevspt()
# Plots spawning biomass for all species and stocks
# in a given rep file. Wraps plotRt(), and overwrites
# that functions margins.
plotRsptDev <- function(  repObj = repInit,
                          initYear = fYear )
{
  # Pull model dimensions
  nS <- repObj$nS
  nP <- repObj$nP
  nT <- repObj$nT

  specNames <- dimnames(repObj$R_spt)[[1]]
  stockNames <- dimnames(repObj$R_spt)[[2]]

  # Create years vector
  years     <- seq(from = initYear, length = nT + 1, by = 1)
  vertLines <- seq(from = initYear, to = max(years), by = 10)

  # Set up plotting environment
  par(mfcol = c(nP, nS), mar = c(0,1,0,1), oma = c(4,4,3,3) )

  for( sIdx in 1:nS )
    for( pIdx in 1:nP )
    {
      plotRtDev(  repObj = repObj,
                  initYear = initYear,
                  sIdx = sIdx, pIdx = pIdx,
                  nopar = TRUE )
      # Detect where we are in the mfg
      mfg <- par("mfg")
      # x axis in bottom row
      if( mfg[1] == mfg[3])
        axis( side = 1)
      # y axis on all plots
      axis( side = 2, las = 1)

      # Add species names
      if( mfg[1] == 1 )
        mtext( side = 3, text = specNames[sIdx] )
      if( mfg[2] == mfg[4] )
        mtext( side = 4, text = stockNames[pIdx] )
    } 
  mtext(  side = 1, text = "Year", outer = TRUE, line = 2 )
  mtext(  side = 2, text = "Recruitment Residuals", outer = TRUE,
          line = 2 )

} # END plotSBspt()



# plotRt()
# Plots the Recruitments over time for a given
# species/stock combo. Used individually and wrapped
# in plotSBspt()
plotRt <- function( repObj = repInit,
                    initYear = fYear,
                    sIdx = 1, pIdx = 1,
                    nopar = FALSE )
{
  # Pull model dimensions
  nS <- repObj$nS
  nP <- repObj$nP
  nT <- repObj$nT

  # Pull spawning biomass
  R_t   <- repObj$R_spt[sIdx, pIdx, ]
  R0    <- repObj$R0_sp[sIdx,pIdx]

  # Create years vector
  years     <- seq(from = initYear, length = nT + 1, by = 1)
  vertLines <- seq(from = initYear, to = max(years), by = 10)

  if( !nopar )
    par(mfcol = c(1, 1), mar = c(2,2,1,1), oma = c(1,1,1,1) )

  plot( x = range(years), y = c(0,max(R_t, na.rm =T) ),
        type = "n", xlab = "", ylab = "",
        las = 1, axes = FALSE )
    # Plot axes if not wrapped
    if(!nopar)
    {
      axis(side = 2, las =1 )
      axis(side = 1)
      mtext(side = 2, text = "Recruitment (1e6)", line = 3)
    }
    box()
    # Plot recruitment
    abline( v = vertLines, lwd = .8, lty = 3, col = "grey80")
    abline( h = R0, lty = 2, lwd = 1, col = "grey50")
    panLab( x = 0.8, y = 0.85, txt = paste("R0 = ", round(R0,2), sep = "") )
    lines( x = years[1:nT], y = R_t[1:nT], lwd = 2, col = "grey30" )
    points( x = years[1:nT], y = R_t[1:nT], pch = 21, bg = "white")
    
} # END plotRt()


# plotRspt()
# Plots recruitment for all species and stocks
# in a given rep file. Wraps plotRt(), and overwrites
# that functions margins.
plotRspt <- function( repObj = repInit,
                      initYear = fYear )
{
  # Pull model dimensions
  nS <- repObj$nS
  nP <- repObj$nP
  nT <- repObj$nT

  # Pull spawning biomass
  R_spt  <- repObj$SB_spt

  specNames <- dimnames(repObj$R_spt)[[1]]
  stockNames <- dimnames(repObj$R_spt)[[2]]
  
  # Create years vector
  years     <- seq(from = initYear, length = nT + 1, by = 1)
  vertLines <- seq(from = initYear, to = max(years), by = 10)

  # Set up plotting environment
  par(mfcol = c(nP, nS), mar = c(0,1,0,1), oma = c(4,4,3,3) )

  for( sIdx in 1:nS )
    for( pIdx in 1:nP )
    {
      plotRt( repObj = repObj,
              initYear = initYear,
              sIdx = sIdx, pIdx = pIdx,
              nopar = TRUE )
      # Detect where we are in the mfg
      mfg <- par("mfg")
      # x axis in bottom row
      if( mfg[1] == mfg[3])
        axis( side = 1)
      # y axis on all plots
      axis( side = 2, las = 1)

      # Add species names
      if( mfg[1] == 1 )
        mtext( side = 3, text = specNames[sIdx] )
      if( mfg[2] == mfg[4] )
        mtext( side = 4, text = stockNames[pIdx] )
    } 
  mtext(  side = 1, text = "Year", outer = TRUE, line = 2 )
  mtext(  side = 2, text = "Recruitment (1e6)", outer = TRUE,
          line = 2 )

} # END plotSBspt()


plotALFreq <- function( repObj = repOpt,
                        sIdx = 1:5, pIdx = 1:3, fIdx = 1:7 )
{
  ALFreq_spalft <- repObj$ALK_spalft[sIdx,pIdx,,,fIdx,,drop = FALSE]

  ALFreq.df <- melt(ALFreq_spalft) %>%
                group_by( species, stocks, ages, lengths ) %>%
                summarise( value = sum(value) ) %>%
                ungroup() %>%
                filter( value > 0 )

  # Get colours
  cols <- wes_palette("Zissou1", 50, type = "continuous")

  tmpPlot <- ggplot( data = ALFreq.df ) +
              geom_tile( aes( x = ages, y = lengths, fill = value ) ) +
              facet_grid( stocks ~ species,
                      scale = "fixed" ) +
              scale_fill_gradientn(colours = cols) + 
              theme_sleek()

  print(tmpPlot)

}


# plotSBt()
# Plots the spawning biomass over time for a given
# species/stock combo. Used individually and wrapped
# in plotSBspt()
plotSBt <- function(  repObj = repInit,
                      initYear = fYear,
                      sIdx = 1, pIdx = 1,
                      nopar = FALSE )
{
  # Pull model dimensions
  nS <- repObj$nS
  nP <- repObj$nP
  nT <- repObj$nT
  nF <- repObj$nF

  # Pull spawning biomass
  SB_t  <- repObj$SB_spt[sIdx, pIdx, ]
  B0    <- round(repObj$B0_sp[sIdx,pIdx],2)
  M_x   <- round(repObj$M_spx[sIdx,pIdx,],2)
  C_ft  <- repObj$C_spft[sIdx, pIdx, , ]
  Bv_ft <- repObj$Bv_spft[sIdx,pIdx,,]
  q_ft  <- repObj$q_spft[sIdx,pIdx,,]
  I_ft  <- repObj$I_spft[sIdx,pIdx,,]


  # replace -ve indices with NA
  I_ft[I_ft < 0] <- NA

  # Scale indices for each fleet
  scaledI_ft <- I_ft
  for( f in 1:nF )
  {
    scaledI_ft[f,] <- I_ft[f,] / q_ft[f,] *( SB_t / Bv_ft[f,])
  }

  gearCols <- brewer.pal(nF, "Dark2")
  gearCols <- alpha(gearCols, alpha = .6)

  # Sum catch over fleets
  C_t   <- apply( X = C_ft, FUN = sum, MARGIN = c(2), na.rm = T )

  # Let's rescale the catch so we can see it comparable to the 
  # biomass
  maxCt <- max(C_t)
  maxBt <- max(SB_t)

  # make a vector of unscaled catch axis ticks
  catchTicks <- seq(0, ceiling(maxCt), by = .5 )

  # Compute the floor of the scalar to get maxCt = minBt,
  # and subtract 1 so catch doesn't get too close
  catchScalar <- max(floor(maxBt * 0.2 / maxCt) - 1, 1)

  # Create years vector
  years     <- seq(from = initYear, length = nT + 1, by = 1)
  vertLines <- seq(from = initYear, to = max(years), by = 10)

  if( !nopar )
    par(mfcol = c(1, 1), mar = c(2,2,1,1), oma = c(1,1,1,1) )

  plot( x = range(years), y = c(0,max(SB_t, C_t, scaledI_ft, na.rm =T) ),
        type = "n", xlab = "", ylab = "",
        las = 1, axes = FALSE )
    # Plot axes if not wrapped
    if(!nopar)
    {
      axis(side = 2, las =1 )
      axis(side = 1)
      mtext(side = 2, text = "Biomass (kt)", line = 3)
    }

    box()
    # Plot recruitment
    abline( v = vertLines, lwd = .8, lty = 2, col = "grey80")
    abline( h = B0, lty = 2, lwd = .8, col = "red")
    panLab( x = 0.8, y = c(0.9,0.85), txt = paste("M = ", M_x, sep = "") )
    panLab( x = 0.8, y = 0.95, txt = paste("B0 = ", B0, sep = "") )
    lines( x = years[1:nT], y = SB_t[1:nT], lwd = 2, col = "red" )
    rect( xleft = years[1:nT] - .3,
          xright = years[1:nT] + .3,
          ybottom = 0, ytop = catchScalar * C_t,
          col = "grey40", border=NA )
    for( f in 1:nF)
    {
      points( x = years[1:nT], y = scaledI_ft[f,],
              col = gearCols[f], pch = 16, cex = .8 )
    }
    axis( side = 4, at = catchScalar * catchTicks,
          labels = catchTicks, las = 1 )
    
} # END plotSBt()

# plotSBspt()
# Plots spawning biomass for all species and stocks
# in a given rep file. Wraps plotSBt(), and overwrites
# that functions margins.
plotSBspt <- function(  repObj = repInit,
                        initYear = fYear )
{
  # Pull model dimensions
  nS <- repObj$nS
  nP <- repObj$nP
  nT <- repObj$nT

  # Pull spawning biomass
  SB_spt  <- repObj$SB_spt
  C_spft  <- repObj$C_spft

  specNames <- dimnames(repObj$SB_spt)[[1]]
  stockNames <- dimnames(repObj$SB_spt)[[2]]

  # Sum catch over fleets
  C_spt   <- apply( X = C_spft, FUN = sum, MARGIN = c(1,2,4) )

  # Create years vector
  years     <- seq(from = initYear, length = nT + 1, by = 1)
  vertLines <- seq(from = initYear, to = max(years), by = 10)

  # Set up plotting environment
  par(mfcol = c(nP, nS), mar = c(0,1,0,1), oma = c(4,4,3,3) )

  for( sIdx in 1:nS )
    for( pIdx in 1:nP )
    {
      plotSBt(  repObj = repObj,
                initYear = initYear,
                sIdx = sIdx, pIdx = pIdx,
                nopar = TRUE )
      # Detect where we are in the mfg
      mfg <- par("mfg")
      # x axis in bottom row
      if( mfg[1] == mfg[3])
        axis( side = 1)
      # y axis on all plots
      axis( side = 2, las = 1)

      # Add species names
      if( mfg[1] == 1 )
        mtext( side = 3, text = specNames[sIdx] )
      if( mfg[2] == mfg[4] )
        mtext( side = 4, text = stockNames[pIdx] )
    } 
  mtext(  side = 1, text = "Year", outer = TRUE, line = 2 )
  mtext(  side = 2, text = "Biomass and Catch (kt)", outer = TRUE,
          line = 2 )

} # END plotSBspt()

# plotProbLenAge()
plotProbLenAge_sp <- function( repObj = repInit )
{
  # Count species and stocks
  nS <- repObj$nS
  nP <- repObj$nP

  specNames <- dimnames(repObj$R_spt)[[1]]
  stockNames <- dimnames(repObj$R_spt)[[2]]

  # Set up plotting environment
  par(mfcol = c(nP, nS), mar = c(0,1,0,1), oma = c(4,4,3,3) )

  for( s in 1:nS )
    for( p in 1:nP )
    {
      plotProbLenAge( repObj = repObj, pIdx = p, sIdx = s,
                      nopar = TRUE )
    
      # Detect where we are in the mfg
      mfg <- par("mfg")
      # x axis in bottom row
      if( mfg[1] == mfg[3])
        axis( side = 1)
      # y axis on all plots
      axis( side = 2, las = 1)

      box()

      # Add species names
      if( mfg[1] == 1 )
        mtext( side = 3, text = specNames[s] )
      if( mfg[2] == mfg[4] )
        mtext( side = 4, text = stockNames[p] )
    }
  
} # END plotProbLenAge()

# plotProbLenAge()
plotProbLenAge <- function( repObj = repInit,
                            sIdx = 1, pIdx = 1,
                            nopar = FALSE )
{
  # Get probability matrix
  probLenAge_lax <- repObj$probLenAge_laspx[,,sIdx,pIdx,]

  # Get max length and ages classes
  L   <- repObj$L_s[sIdx]
  A   <- repObj$A_s[sIdx]
  nX  <- repObj$nX

  cols <- matrix(NA, nrow = nX, ncol = A)
  if( nX > 1)
    sexCols <- c("steelblue","salmon")
  else sexCols <- "grey30"
  for( x in 1:nX )
    cols[x,] <- alpha(sexCols[x], alpha = seq(from = .2, to = .8, length = A) )

  # Open plotting window
  plot( x = c(0,L), y = range(probLenAge_lax), type = "n", xlab = "",
        ylab = "", axes = FALSE )
    if( !nopar )
    {
      # x axis in bottom row
      if( mfg[1] == mfg[3])
        axis( side = 1)
      # y axis on all plots
      axis( side = 2, las = 1)
      box()
    }
    for(x in 1:nX)
      for( aIdx in 1:A )
        lines(  x = 1:L, y = probLenAge_lax[1:L,aIdx,x],
                col = cols[x,aIdx], lwd = 1 )
} # END plotProbLenAge()

# plotYeq
# Plot equilibrium yield reference curve as a function
# of F
plotYeqF <- function( repObj = repOpt,
                      sIdx = 1, pIdx = 1 )
{
  # Pull reference points object
  refPoints <- repObj$refPts
  refCurves <- refPoints$refCurves

  # Pull eqbm yields
  YeqF_spf     <- refCurves$Yeq_sp[sIdx,pIdx,,drop = FALSE]
  Fmsy_sp      <- refPoints$FmsyRefPts$Fmsy_sp[sIdx,pIdx,drop = FALSE]
  YeqFmsy_sp   <- refPoints$FmsyRefPts$YeqFmsy_sp[sIdx,pIdx,drop = FALSE]

  # Now melt
  YeqF.df     <- melt(YeqF_spf) %>%
                  rename(yield = value) %>%
                  filter( yield >= 0 )
  Fmsy.df     <- melt(Fmsy_sp) %>%
                  rename( Fmsy = value)
  YeqFmsy.df  <- melt(YeqFmsy_sp) %>%
                  rename( yield = value) %>%
                  left_join(Fmsy.df)

  # Figure out how to stack multiple YPR/SPR ref point
  # tables together, with a column for the ref point
  # label, for easier plotting of multiple ref points

  tmp <-  ggplot(data = YeqF.df, aes(x=F, y=yield)) + 
          geom_line() +
          facet_grid( stock ~ species, scale = "fixed") +
          theme_sleek() +
          geom_point( data = YeqFmsy.df, inherit.aes = FALSE,
                      mapping = aes( x = Fmsy, y = yield ),
                      col = "red", size = 1.5 )

  print(tmp)

  return(tmp)
} # END plotYeqF()

# plotProbLenAge()
plotHeatmapProbLenAge <- function(  repObj = repInit,
                                    sIdx = 1, pIdx = 1,
                                    facets = c("stock","species") )
{
  # Get probability matrix
  probLenAge_lax <- repObj$probLenAge_laspx[,,sIdx,pIdx, , drop = FALSE]

  # dimnames(probLenAge_la) <- list(  length = 1:nrow(probLenAge_la),
  #                                   age = 1:ncol(probLenAge_la) )

  # Get max length and ages classes
  l <- repObj$minL_s[sIdx]
  L <- repObj$L_s[sIdx]
  A <- repObj$A_s[sIdx]

  nS <- length(sIdx)
  nP <- length(pIdx)

  # Get vonB parameters
  A1_s      <- repObj$A1_s[sIdx,drop = FALSE]
  A2_s      <- repObj$A2_s[sIdx,drop = FALSE]
  vonK_spx  <- repObj$vonK_spx[sIdx,pIdx,, drop = FALSE]
  L1_spx    <- repObj$L1_spx[sIdx,pIdx,, drop = FALSE]
  L2_spx    <- repObj$L2_spx[sIdx,pIdx,, drop = FALSE]

  assignAgePar <- function( spec, ageVec )
  {
    A <- ageVec[spec]

    A
  }

  # Combine vonB parameters into a data.frame
  vonK.df <- melt(vonK_spx) %>%
              mutate( par = "vonK")
  L1.df   <- melt(L1_spx) %>%
              mutate( par = "L1",
                      zero = 0 ) %>%
              group_by( species ) %>%
              mutate( age = assignAgePar( species, A1_s ) ) %>%
              ungroup()
  L2.df   <- melt(L2_spx) %>%
              mutate( par = "L2",
                      zero = 0 ) %>%
              group_by( species ) %>%
              mutate( age = assignAgePar( species, A2_s ) ) %>%
              ungroup()



  cols <- wes_palette("Zissou1", 50, type = "continuous")


  melted_probMat <- melt(probLenAge_lax)
  melted_probMat <- melted_probMat %>% filter( value > 1e-4 )

  predLenAtAge_ax <- repObj$lenAge_aspx[,sIdx,pIdx,, drop = FALSE]
  predLenAtAge_ax <- melt(predLenAtAge_ax) 




  tmp <-  ggplot(data = melted_probMat, aes(x=age, y=len, fill=value)) + 
          geom_tile() +
          facet_grid( stock ~ species + sex,
                      scale = "fixed" ) +
          scale_fill_gradientn(colours = cols) +
          geom_line(  data = predLenAtAge_ax, 
                      mapping = aes(x = age, y = value),
                      inherit.aes = FALSE, size = 1) +
          geom_segment( data = L1.df, inherit.aes = FALSE,
                        mapping = aes(  x = age, y = zero,
                                        xend = age, yend = value), 
                        size = .8, alpha = .4,
                        linetype = 2 ) +
          geom_segment( data = L1.df, inherit.aes = FALSE,
                        mapping = aes(  x = zero, y = value,
                                        xend = age, yend = value), 
                        size = .8, alpha = .4,
                        linetype = 2 ) +
          geom_segment( data = L2.df, inherit.aes = FALSE,
                        mapping = aes(  x = age, y = zero,
                                        xend = age, yend = value), 
                        size = .8, alpha = .4,
                        linetype = 2 ) +
          geom_segment( data = L2.df, inherit.aes = FALSE,
                        mapping = aes(  x = zero, y = value,
                                        xend = age, yend = value), 
                        size = .8, alpha = .4,
                        linetype = 2 ) +
          theme_sleek()

  print(tmp)
  
} # END plotProbLenAge()

# plotHeatmapAgeLenResids()
plotHeatmapAgeLenResids <- function(  repObj = repInit,
                                      sIdx = 1, pIdx = 1,
                                      fIdx = 1 )
{
  # Get resids matrix
  resids <- repObj$ageAtLenResids_alspftx[,,sIdx,pIdx,fIdx,,,drop = FALSE]
  resids <- apply( X = resids, FUN = mean, MARGIN = c(1,2,3,4,5,7))

  resids[resids == 0] <- NA
  

  # Get max length and ages classes
  L <- repObj$maxL_s[sIdx]
  l <- repObj$minL_s[sIdx]
  A <- repObj$maxA_s[sIdx]

  nS <- length(sIdx)
  nP <- length(pIdx)
  nF <- length(fIdx)

  cols <- wes_palette("Zissou1", 50, type = "continuous")
  
  melted_residsMat <- melt( resids ) 
  
  nL <- repObj$nL
  nA <- repObj$nA

  predLenAtAge <- repObj$lenAge_aspx[,sIdx,pIdx,, drop = FALSE]
  predLenAtAge <- melt(predLenAtAge) 


  tmp <- ggplot(data = melted_residsMat, aes(x=age, y=len, fill=value)) + 
          geom_tile() +
          facet_grid( stock ~ species + fleet + sex, scale = "fixed" ) +
          scale_fill_gradientn(colours = cols, na.value = NA) +
          geom_line(  data = predLenAtAge, 
                      mapping = aes(x = age, y = value),
                      inherit.aes = FALSE,
                      size = 2 ) +
          theme_sleek()

  print(tmp)

} # END plotHeatmapAgeLenResids()

# plotCompResids()
# Generic function to plot standardised residuals
# from age or length compositions
plotCompResids <- function( repObj = repInit,
                            sIdx = 1, pIdx = 1,
                            fIdx = 1, sex = "male",
                            comps = "age" )
{
  # First, get residuals
  if( comps == "age" )
  {
    # Pull resids array
    resids_xspft <- repObj$ageRes_aspftx[,sIdx,pIdx,fIdx,,sex,drop = FALSE]
    tauRes_spf   <- sqrt(repObj$tau2Age_spf[sIdx,pIdx,fIdx,drop = FALSE])
  }

  if( comps == "length")
  {
    # Pull resids array
    resids_xspft <- repObj$lenRes_lspftx[,sIdx,pIdx,fIdx,,sex,drop = FALSE]
    tauRes_spf   <- sqrt(repObj$tau2Len_spf[sIdx,pIdx,fIdx,drop = FALSE])
  }
  # Get dimension
  nS <- length(sIdx)
  nP <- length(pIdx)
  nF <- length(fIdx)
  # Colours
  cols <- wes_palette("Zissou1", 3, type = "continuous")
  cols <- rev(as.vector(cols))
  cols[2] <- NA
  names(cols) <- c("-","Zero","+")



  # Melt residuals into a data.frame
  melted_residsArr <- melt( resids_xspft ) %>%
                      rename( resids = value )
  melted_tauArr    <- melt( tauRes_spf ) %>%
                      rename( sd = value) 

  resids.df <- left_join( melted_residsArr, melted_tauArr, 
                          by = c("species","stock","fleet") ) %>%
                mutate( stdRes = resids/sd,
                        size = abs(stdRes),
                        sign = ifelse(resids > 0, "+", ifelse( resids == 0, "Zero","-" ) ) )
  


  tmpPlot <-  ggplot( data = resids.df,
                      aes(  x = year,
                            y = !!ensym(comps) ) ) +
              facet_grid(stock ~ species + fleet, scale = "fixed" ) +
              geom_point( aes(  size = size,
                                colour = sign ),
                          shape = 21,
                          fill = NA ) +
              scale_colour_manual(values = cols) +
              theme_sleek()


  print(tmpPlot)

}


# plotMatAge()
plotMatAge <- function( repObj = repInit,
                        sIdx = 1, pIdx = 1,
                        sex = "male" )
{
  # Get probability matrix
  probLenAge_la <- repObj$probLenAge_laspx[,,sIdx,pIdx,sex]
  matLen_l      <- repObj$matLen_ls[,sIdx]
  matAge_a      <- repObj$matAge_asp[,sIdx,pIdx]

  # replace -1s with NA
  matLen_l[matLen_l < 0] <- NA

  # Get max length and ages classes
  L <- repObj$L_s[sIdx]
  A <- repObj$A_s[sIdx]


  # Open plotting window
  plot( x = c(0,A), y = c(0,1), type = "n", xlab = "Age",
        ylab = "Maturity-at-age", las = 1 )
    lines(  x = 1:A, y = matAge_a[1:A],
            col = "steelblue", lwd = 1 )
} # END plotMatAge()

# plotMatLength()
plotMatLength <- function( repObj = repInit,
                           sIdx = 1, pIdx = 1,
                           sex = "male" )
{
  # Get probability matrix
  probLenAge_la <- repObj$probLenAge_laspx[,,sIdx,pIdx,sex]
  matLen_l      <- repObj$matLen_ls[,sIdx]
  matAge_a      <- repObj$matAge_asp[,sIdx,pIdx]

  # replace -1s with NA
  matLen_l[matLen_l < 0] <- NA

  # Get max length and ages classes
  L <- repObj$L_s[sIdx]


  # Open plotting window
  plot( x = c(0,L), y = c(0,1), type = "n", xlab = "Length (cm)",
        ylab = "Maturity-at-length", las = 1 )
    lines(  x = 1:L, y = matLen_l[1:L],
            col = "steelblue", lwd = 1 )
} # END plotMatLength()

# plotWtAge()
plotWtAge <- function(  repObj = repInit,
                        sIdx = 1, pIdx = 1 )
{
  # Get probability matrix
  meanWtAge_ax <- repObj$meanWtAge_aspx[,sIdx,pIdx,]

  # Get max age classes
  A <- repObj$A_s[sIdx]

  nX <- repObj$nX

  sexCols <- c("steelblue","salmon")

  # Open plotting window
  plot( x = c(0,A), y = c(0,max(meanWtAge_ax,na.rm =T)), type = "n", 
        xlab = "Age",
        ylab = "Weight (kg)", las = 1 )
    for( x in 1:nX)
      lines(  x = 1:A, y = meanWtAge_ax[1:A,x],
              col = sexCols[x], lwd = 2 )
} # END plotWtAge()

# plotLenAge()
plotLenAge <- function( repObj = repInit,
                        sIdx = 1, pIdx = 1 )
{
  # Get probability matrix
  lenAge_ax <- repObj$lenAge_aspx[,sIdx,pIdx,]

  # Get max age classes
  A <- repObj$A_s[sIdx]

  nX <- repObj$nX
  sexCols <- c("steelblue","salmon")

  # Open plotting window
  plot( x = c(0,A), y = c(0,max(lenAge_ax,na.rm =T)), type = "n", 
        xlab = "Age",
        ylab = "Length (cm)", las = 1 )
    for( x in 1:nX)
      lines(  x = 1:A, y = lenAge_ax[1:A,x],
              col = sexCols[x], lwd = 2 )
} # END plotLenAge()

# plotIdxFits()
plotIdxFits <- function ( repObj = repInit,
                          initYear = fYear,
                          sIdx = 1, pIdx = 1 )
{
  # Pull vulnerable biomass and indices
  Bv_ft       <- repObj$Bv_spft[sIdx,pIdx,,]
  # vulnNtg     <- repObj$vulnN_tg
  I_ft        <- repObj$I_spft[sIdx,pIdx,,]
  q_f         <- repObj$q_spf[sIdx,pIdx,]
  q_ft        <- repObj$q_spft[sIdx,pIdx,,]
  tau_f       <- repObj$tau2Idx_spf[sIdx,pIdx,]

  gearNames   <- dimnames(Bv_ft)[[1]]

  # Model dimensions
  nF          <- repObj$nF
  nT          <- repObj$nT

  # replace missing indices with NAs
  I_ft[ I_ft < 0 ] <- NA
  # Now loop over gears, choose 
  surveyGears <- c()
  for( f in 1:nF )
    if(any(!is.na(I_ft[f,])))
      surveyGears <- c(surveyGears,f)

  par(  mfrow = c(length(surveyGears),1), 
        oma = c(3,4,1,1),
        mar = c(1,1,1,1) )

  years <- seq(initYear, by = 1, length = nT + 1 )
  vertLines <- seq(from = initYear, to = max(years), by = 10)

  gearCols <- brewer.pal(nF, "Dark2")

  for( fIdx in surveyGears )
  {
    scaledIndices <- I_ft[fIdx,]/q_ft[fIdx,]
    maxRange <- max(Bv_ft[fIdx,], scaledIndices,na.rm = T)
    plot( x = range(years), c(0,maxRange),
          xlab = "", ylab = "", type = "n",
          axes = F )
      axis(side = 1)
      axis(side = 2, las = 1)
      box()
      lines( x = years[1:nT], y = Bv_ft[fIdx,], lwd = 2, col = "grey30" )
      points( x = years[1:nT], y = scaledIndices, col = gearCols[fIdx], pch = 16 )
      panLab( x = 0.9, y = 0.8, txt = gearNames[fIdx], cex = 2 )
      panLab( x = 0.9, y = .7, txt = paste("q = ",  round(q_f[fIdx],2),sep = ""), cex = 1.5 )
      panLab( x = 0.9, y = .65, txt = paste("tau = ", round(tau_f[fIdx],2),sep = ""), cex = 1.5 )
  }
  mtext(  side = 1, text = "Year", outer = T )
  mtext(  side = 2, text = "Vulnerable Biomass", outer = T, 
          line = 2)
} # END plotIdxFits()

# plotIdxResids()
plotIdxResids <- function ( repObj = repInit,
                            initYear = fYear,
                            sIdx = 1, pIdx = 1 )
{
  # Pull vulnerable biomass and indices
  I_ft_hat    <- repObj$I_spft_hat[sIdx,pIdx,,]
  I_ft        <- repObj$I_spft[sIdx,pIdx,,]

  gearNames   <- dimnames(I_ft_hat)[[1]]

  # Model dimensions
  nF          <- repObj$nF
  nT          <- repObj$nT

  # replace missing indices with NAs
  I_ft[ I_ft < 0 ] <- NA
  I_ft_hat[ I_ft_hat < 0 ] <- NA
  # Now loop over gears, choose 
  surveyGears <- c()
  for( f in 1:nF )
    if(any(!is.na(I_ft[f,])))
      surveyGears <- c(surveyGears,f)

  par(  mfrow = c(length(surveyGears),1), 
        oma = c(3,4,1,1),
        mar = c(1,1,1,1) )

  years <- seq(initYear, by = 1, length = nT + 1 )
  vertLines <- seq(from = initYear, to = max(years), by = 10)


  for( g in surveyGears )
  {
    skipLo <- FALSE
    resids <- log(I_ft_hat[g,] / I_ft[g,])
    yLim <- range(resids,na.rm = T)
    if(any(!is.finite(range(resids,na.rm = T))))
    { 
      skipLo <- TRUE
      yLim <- c(-1,1)
    }
    nonNA <- which(!is.na(resids))

    plot( x = range(years), y = yLim,
          xlab = "", ylab = "", type = "n",
          axes = F )
      axis(side = 1)
      axis(side = 2, las = 1)
      box()
      abline(h = 0, lty = 2, lwd = .8 )
      points( x = years[1:nT], y = resids )
      if(!skipLo)
        lines(  loess.smooth(  x = years[nonNA], y = resids[nonNA]), 
                lwd = 2, col = "salmon" )
      
      panLab( x = 0.9, y = 0.8, txt = gearNames[g], cex = 2 )
  }
  mtext(  side = 1, text = "Year", outer = T )
  mtext(  side = 2, text = "Vulnerable Biomass", outer = T, 
          line = 2)
} # END plotIdxFits()

# plotSelLen()
# Plot selectivity-at-length for a given species/stock
# over all fleets
plotSelLen <- function( repObj = repInit,
                        sIdx = 1, pIdx = 1 )
{
  # Get selectivity functions
  sel_lft <- repObj$sel_lfspt[ , , sIdx, pIdx, ]

  gearNames <- dimnames(sel_lft)[[2]]

  # Dimensions
  nF <- repObj$nF
  L  <- repObj$L_s[sIdx]
  nT <- repObj$nT


  # Get species/stock/fleet selectivity parameter values
  xSel50_f    <- repObj$xSel50_spf[sIdx,pIdx,]
  xSel95_f    <- repObj$xSel95_spf[sIdx,pIdx,]
  xSelStep_f  <- repObj$xSelStep_spf[sIdx,pIdx,]

  # get species/fleet
  if(!is.null(repObj$xSel50_sf))
  {
    xSel50_sf   <- repObj$xSel50_sf[sIdx,]
    xSel95_sf   <- repObj$xSel95_sf[sIdx,]
    xSelStep_sf <- repObj$xSelStep_sf[sIdx,]

    sel_lfs  <- matrix(0, nrow = L, ncol = nF)
  }

  # Get prior mean pars
  pmxSel50_sf   <- exp(repObj$pmlnxSel50_sf[sIdx,])
  pmxSelStep_sf <- exp(repObj$pmlnxSelStep_sf[sIdx,])
  pmxSel95_sf   <- pmxSel50_sf + pmxSelStep_sf

  
  sel_lf    <- matrix(0, nrow = L, ncol = nF)
  pmSel_lf  <- matrix(0, nrow = L, ncol = nF)

  for( fIdx in 1:nF )
  {
    sel_lf[,fIdx] <- 1 / (1 + exp(-log(19) * (1:L - xSel50_f[fIdx])/xSelStep_f[fIdx]))
    pmSel_lf[,fIdx] <- 1 / (1 + exp(-log(19) * (1:L - pmxSel50_sf[fIdx])/pmxSelStep_sf[fIdx]))
    if(!is.null(sel_lfs))
      sel_lfs[,fIdx] <- 1 / (1 + exp(-log(19) * (1:L - xSel50_sf[fIdx])/xSelStep_sf[fIdx]))
  }

  par(  mfrow =c(nF,1), 
        oma = c(3,4,1,3),
        mar = c(0,1,0,1) )
  for( fIdx in 1:nF )
  {
    plot( x = c(1,L), y = c(0,1), type = "n",
          axes = FALSE, xlab = "", ylab = "" )
      # plot axes at the bottom and left sides
      mfg <- par("mfg")
      if(mfg[1] == mfg[3])
        axis( side = 1)
      axis( side = 2, las = 1)
      box()
      # Plot time-varying selectivity for this fleet on this stock
      for( t in 1:nT )
        lines(  x = 1:L, y = sel_lft[1:L,fIdx,t], col = "grey75",
                lwd = 1 )
      # Plot mean selectivity for this species/fleet
      lines( x = 1:L, y = sel_lf[1:L,fIdx], col = "salmon", lwd = 2)
      lines( x = 1:L, y = pmSel_lf[1:L,fIdx], col = "steelblue", lwd = 2, lty = 2 )
      if(!is.null(sel_lfs))
        lines( x = 1:L, y = sel_lfs[1:L,fIdx], col = "grey50", lwd = 2, lty = 3 )
      mtext( side = 4, text = gearNames[fIdx], line = 2, font = 2)
      if(fIdx == 1)
        legend(x = "bottomright",
                legend = c("Prior","Estimate"),
                col = c("salmon","steelblue"),
                lty = c(1,2),
                lwd = c(2,2), bty = "n" )
  }
  mtext( side = 1, text= "Length", outer = T, line = 2)
  mtext( side = 2, text= "Selectivity-at-length", outer = T, line = 2)

} # END plotSelLen()


# plotSelAge()
# Plot selectivity-at-age for a given species/stock
# over all fleets
plotSelAge <- function( repObj = repInit,
                        sIdx = 1, pIdx = 1 )
{
  # Get selectivity functions
  sel_aftx <- repObj$sel_afsptx[ , , sIdx, pIdx, ,  ]

  # Get gear labels
  gearNames <- dimnames(sel_aftx)[[2]]

  # Dimensions
  nF <- repObj$nF
  A  <- repObj$A_s[sIdx]
  nT <- repObj$nT
  nX <- repObj$nX

  sexCols <- c("steelblue","salmon")
  sexCols <- alpha( sexCols, alpha = .5 )

  par(  mfrow =c(nF,1), 
        oma = c(3,4,1,3),
        mar = c(0,1,0,1) )
  for( fIdx in 1:nF )
  {
    plot( x = c(1,A), y = c(0,1), type = "n",
          axes = FALSE, xlab = "", ylab = "" )
      # plot axes at the bottom and left sides
      mfg <- par("mfg")
      if(mfg[1] == mfg[3])
        axis( side = 1)
      axis( side = 2, las = 1)
      box()
      for( x in 1:nX )
        for( t in 1:nT)
          lines( x = 1:A, y = sel_aftx[1:A,fIdx,t,x], col = sexCols[x],
                  lwd = 2 )
      mtext( side = 4, text = gearNames[fIdx], line = 2, font = 2)
  }
  mtext( side = 1, text= "Age", outer = T, line = 2)
  mtext( side = 2, text= "Selectivity-at-age", outer = T, line = 2)

} # END plotSelAge()


