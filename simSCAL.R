# <><><><><><><><><><><><><><><><><><><><><><><><><><><>
# simSCAL.R
#
# Simulation model for generating data
# for testing a statistical catch at age 
# and length (SCAL) model.
#
# Author: SDN Johnson
# Date: July 25, 2019
#
# Last Update: July 25, 2019
#
# Ideas from Tom for hierSCAL:
# 1. Master index/effective effort for catch
# 2. ALK with fishing out the percentiles
#
# <><><><><><><><><><><><><><><><><><><><><><><><><><><>


simAgeLenPop <- function( B0        = 100,
                          h         = .78,
                          M_x       = c(.2,.2),
                          nG        = 2,
                          nX        = 2,
                          nT        = 65,
                          nA        = 35,
                          nL        = 1,
                          q_g       = c(.1,.3),
                          aMat50    = 4,
                          aMat95    = 7.5,
                          aSel50_g  = c(2,5),
                          aSel95_g  = c(5,9),
                          Linf_x    = c(34,34),
                          K_x       = c(.25,.25),
                          L1        = 14,
                          sigmaL    = .15 )
{
  # Set up arrays to hold simulated states
  B_alxt <- array(NA,  dim = c(nA,nL,nX,nT) )
  N_alxt <- array(NA,  dim = c(nA,nL,nX,nT) )

  # Catch
  C_alxgt   <- array(NA,  dim = c(nA,nL,nX,nG,nT) ) # Numbers
  Cw_alxgt  <- array(NA,  dim = c(nA,nL,nX,nG,nT) ) # Weight

  # probability of being in a growth group
  probGrowth_lx <- array(NA, dim = c(nL,nX) )

  # Calculate growth models - 
  # need one curve for each GG

  # Calculate ssbpr

  # Calculate R0, BH SR pars

  # Start simulating 


  




}