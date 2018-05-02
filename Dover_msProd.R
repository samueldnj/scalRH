# -----------------------------------------------------------------------------
#
# DERPAmsProd.R
#
# SDN Johnson
# May 23, 2017
# 
# Reads in survey CPUE and commercial catch data from the
# working directoty and ./density/ folders. Used to create CPUE
# CPUE and catch data arrays for fitting the msProd production model
# in TMB
# 
# -----------------------------------------------------------------------------

rm(list = ls())

library(dplyr)
library(ggmap)
library(maptools)
library(rgdal)
library( PBSmapping ) 
library( RColorBrewer )
library( RgoogleMaps)
library( sp )
library(raster)
library(MASS)
library(psyphy)
library(boot)
library(RCurl)
library(TMB)
# library(TMBphase)

source("mseRtools.r")
source("DERPAfuns.R")

surveys <- c(2,1)

iYear <- 1954


stocksSurv = list(  HG = c(2,3),
                    QCS = c(1),
                    WCVI = c(4) )

stocksComm = list(  HG = c(7,8,9),
                    QCS = c(5,6),
                    WCVI = c(3,4) )

# Make relative biomass and catch arrays for feeding to msProd simulation
# framework
relBio  <- makeRelBioStocks( years = c(iYear,2016), spec = "dover", collapseSyn = TRUE, stocks = stocksSurv )
relBio  <- relBio$relBio.arr[surveys,,]
Katch   <- makeStockCatch(years = c(iYear,2016), spec = "dover", stocks = stocksComm )

# save to an Rdata file
data <-  list( indices = relBio, katch = Katch$catch.arr )
save( data, file = "doverStocks.RData" )


relBio[relBio < 0] <- NA

nStocks <- dim(relBio)[2]

years <- c(iYear,2016)

par( mfrow = c(nStocks,1), oma = c(4,5,3,5), mar = c(1,1,0,1) )
for( sIdx in 1:nStocks )
{
  subKatch  <- Katch$catch.arr[sIdx,]
  subBio    <- relBio[,sIdx,]
  plot( x = years, y = c( 0, max( relBio, na.rm=T ) ), axes = F, type = "n",
        xlab = "", ylab = "" )
    if(sIdx == nStocks )axis( side = 1 )
    axis( side = 2, las = 1 )
    mtext( side = 1, outer = T, text = "Year", line = 2)
    mtext( side = 4, outer = T, text = "Minimum Trawlable Biomass (Kt)", line = 3, las = 0 )

    rect( xleft = years[1]:years[2] - .3, xright = years[1]:years[2] + .3,
              ybottom = 0, ytop = subKatch, lwd = 3, border =NA, col = "grey70")
    axis( side = 4, las = 1 )
    mtext( side = 2, outer = T, text = "Landings and Discards (Kt)", line = 3, las = 0 )
    for( i in 1:length(surveys) )
    {
      survIdx <- surveys[i]
      points( x = years[1]:years[2], y = subBio[survIdx,], pch = survIdx, cex = 2 )
    }
  panLab( x = 0.1, y=0.9, txt = names(stocksSurv)[sIdx], cex = 2 ) 
  if( sIdx == 1 ) panLegend(  x = 0.3, y = 1,
                              pch = 1:length(surveys), 
                              legTxt = c(paste("Survey",1:length(surveys))),
                              fill = c(rep(NA,length(surveys))),
                              border = c(NA),
                              col = c(rep(1,length(surveys))),
                              bty = "n",
                              cex = 1.5 )
}
# # relBio[ relBio == 0 ] <- -1

# lists <- makeDatPar()

# compile ("msProd.cpp")
# dyn.load(dynlib("msProd"))

# modelName <- "msProd"
# # RE <- c("lnq_os","lnUmsy","zeta_st","eps_t")

# obj <- MakeADFun (  dat = lists$dat, parameters = lists$par, map = lists$map,
#                     random = c("zeta_st","lnq_os"), silent = FALSE )

# ctrl = list ( eval.max = 100, iter.max = 100 )

# # optimise the model
# fit <- nlminb ( start = obj$par,
#                 objective = obj$fn,
#                 gradient = obj$gr,
#                 control = ctrl,
#                 lower = -Inf,
#                 upper= Inf )

# # get the SDs in leading and derived pars, and their
# # MLEs
# sdrepObj  <- sdreport(obj)
# repObj    <- obj$report()

# plotAssess(initYear = iYear, depletion = F)
