// ><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><><>><><>><><>><>
// scalRH.cpp
// 
// A multi-stock age structured stock assessment model with Robin Hood 
// hierarchical priors on biological and fishery parameters
// 
// Author: Samuel Johnson
// Date: 20 August, 2017
// Purpose: This model is intended to estimate leading parameters
//          for an age structured multi-stock operating model, which can
//          be used to test this assessment model in a sim-est 
//          procedure, and later to test multispecies management
//          systems in closed loop simulation
// 
// 
// To Do:
//  4. Consider correlation in the compositional likelihoods
// 
// 
// Intended features:
// - Multistock and multispecies (DERPA)
//     - Species made up of a single stock have only species level 
//       priors
//     - Different initialisation times for each stock - wishlist
// - Multi-sector
//     - Fleets are fisheries by gear type (including surveys)
//     - Selectivity is length based at fleet/species level,
//        but a random stock effect for each species/stock combo
//        is added - could correlate if time-varying
//     - Fishing mortality for commercial fleets is effort
//        based, with a single effort index for each 
//        fleet/stock-area combo being converted to species
//        specific F by comm fleet catchability
// - Multi-level RH priors on:
//     - Growth (vonB)
//     - Fishing mortality (correlation in REs if estimated)
//     - Natural mortality (M0 prior and correlated deviations)
//     - Selectivity
//     - Catchability
//     - S-R Steepness
// - Age or Length observations are used (Logistic Normal comp likelihood)
//     - Age observations are used where available
//     - length comps are included, but not from same fish that
//        provided age comps
//     - Currently, age-length distribution is stationary, could
//        augment to explicitly model growth-groups (like Atl Halibut)
// - Integrated growth model to predict length dist when age data missing
//     - RAL likelihood, to explicitly account for biased sampling
//        from different fleets
//     - This is too much, need to switch to an ALK
// - Discarding - need to talk to fisheries about discarding behaviour
//     - "Grading length" for fisheries with no size limit
// - Two catch equation solution options
//     - Pope's approximation for age structured populations and
//        multifleet systems implemented
//     - NR solver for Baranov equation: can be turned on in last
//        phase, will use Pope's approx as initial value.
// 
// ><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><><>><><>><><>><>

#include <TMB.hpp>                                // Links in the TMB libraries
#include <iostream>
// link personal library of TMB stock assessment functions
// #include "/Users/sdnjohnson/Work/code/tmbFuns/stockAssessmentFuns.hpp"  

// square()
// Shortcut wrapper for squaring quantities
// inputs:    x = quantity to be squared
// ouputs:    y = x^2
// Usage:     when squaring.
// Source:    Almost surely stolen from somewhere.
// square
template<class Type>
Type square(Type x)
{
  return x*x;
}

// VECTORIZE1_t(square)

// dhalfnorm()
// R-style half-normal probability distribution density
// function.
// inputs:    x = real number to calculate density of
//            mean = desired expected value of the half-normal dist
//            log = int determining whether function returns log (1) or 
//                  natural (0) scale density
// outputs:   dens = density on natural or log scale
// Usage:     usually for variance priors
// Source:    Some Gelman paper, function written by S. D. N. Johnson
template<class Type>
Type dhalfnorm( Type x,
                Type mean,
                int logscale )
{
  // First, scale mean to distribution SD
  Type pi = 3.1415926535;
  Type sd = sqrt(pi) * mean / sqrt(2);

  // Calculate density
  Type dens = sqrt(2) * exp( - 0.5 * square(x) / square(sd) ) / sqrt(pi) / sd;

  // log transform
  if( logscale == 1 )
    dens = log(dens);

  return(dens);
}
// VECTORIZE1_t(dhalfnorm)

// dinvgamma()
// R-style inverse gamma probability distribution density
// function.
// inputs:    x = real number to calculate density of
//            alpha = shape parameter
//            beta = scale parameter
//            logscale = int determining whether function returns log (1) or 
//                        natural (0) scale density
// outputs:   dens = density on natural or log scale
// Usage:     usually for variance priors
// Source:    Wikipedia,  by S. D. N. Johnson
template<class Type>
Type dinvgamma( Type x,
                Type alpha,
                Type beta,
                int logscale )
{

  Type dens = -1 * ( (alpha + 1) * log(x) + beta / x );

  if( logscale == 0 )
    dens = exp( dens );

  return(dens);
}

// posfun
template<class Type>
Type posfun(Type x, Type eps, Type &pen){
  pen += CppAD::CondExpLt(x, eps, Type(0.01) * pow(x-eps,2), Type(0));
  return CppAD::CondExpGe(x, eps, x, eps/(Type(2)-eps/x));
}


// calcPopesApprox()
// Takes numbers at age/species/pop/sex, and fleet/spec/pop
// catch and selectivity and produces pope's approximations 
// of fishing mortality rates.
// inputs:    N_axsp    = Numbers at age/spec/pop/sex
//            sel_axspf = selectivity at age/spec/pop/fleet/sex
//            C_spf     = catch in biomass for spec/pop/fleet
//            M_xsp     = Natural mortality by spec/pop/sex
//            A_s       = Number of age classes for each species
//            wt_axsp    = weight at age/species/pop/sex
template<class Type>
array<Type> calcPopesApprox(  array<Type> N_axsp,     // Numbers at age/spec/pop/sex
                              array<Type> sel_axspf,  // Selectivity-at-age (spec/pop/fleet/sex)
                              array<Type> C_spf,      // Catch in biomass units
                              array<Type> M_xsp,      // Natural mortality
                              vector<int> A_s,        // number of age classes by species
                              array<Type> wt_axsp,    // Mean weight at age
                              array<Type> vB_axspf,   // Vulnerable biomass at age/spec/pop/fleet/sex
                              array<Type> vB_spf,     // Total vulnerable biomass for each spec/pop/fleet
                              array<Type>& F_spf,     // Fishing mortality for each spec/pop/fleet
                              array<Type>& Z_axsp,    // Total mortality for each age/spec/pop/sex
                              array<Type>& Cw_axspf,  // Catch weight at sex/age/spec/pop/fleet
                              array<Type>& C_axspf )  // Catch in numbers at sex/age/spec/pop/fleet
{
  // Get dimensions
  int nA = N_axsp.dim(0);
  int nX = N_axsp.dim(1);
  int nS = N_axsp.dim(2);
  int nP = N_axsp.dim(3);
  int nF = vB_axspf.dim(4);

  array<Type> remN_axsp(nA,nX,nS,nP);
  array<Type> newN_axsp(nA,nX,nS,nP);
  remN_axsp.setZero();
  newN_axsp.setZero();
  array<Type> pvB_axspf(nA,nX,nS,nP,nF);
  pvB_axspf.setZero();
  pvB_axspf = vB_axspf;

  for( int sIdx = 0; sIdx < nS; sIdx ++ )
  {
    for( int pIdx = 0; pIdx < nP; pIdx++ )
    {
      for( int fIdx = 0; fIdx < nF; fIdx++ )
      {
        for( int xIdx = 0; xIdx < nX; xIdx++)
        {
          // Pull biomass at age, convert to proportions
          pvB_axspf.col(fIdx).col(pIdx).col(sIdx).col(xIdx) /= vB_spf(sIdx,pIdx,fIdx);

          // Calculate numbers to remove by converting catch to weight
          Cw_axspf.col(fIdx).col(pIdx).col(sIdx).col(xIdx)  = C_spf(sIdx, pIdx, fIdx) * pvB_axspf.col(fIdx).col(pIdx).col(sIdx).col(xIdx);
          C_axspf.col(fIdx).col(pIdx).col(sIdx).col(xIdx)   = Cw_axspf.col(fIdx).col(pIdx).col(sIdx).col(xIdx) / wt_axsp.col(pIdx).col(sIdx).col(xIdx);

          // Now add the catch at age in numbers to the removed fish
          remN_axsp.col(pIdx).col(sIdx).col(xIdx) += C_axspf.col(fIdx).col(pIdx).col(sIdx).col(xIdx);
        } // END x loop

      } // END f loop

      for( int xIdx = 0; xIdx < nX; xIdx++ )
      {
        // Apply half M
        newN_axsp.col(pIdx).col(sIdx).col(xIdx) = N_axsp.col(pIdx).col(sIdx).col(xIdx) * exp(-M_xsp(xIdx,sIdx,pIdx)/2);
        // Remove the catch
        newN_axsp.col(pIdx).col(sIdx).col(xIdx) -= remN_axsp.col(pIdx).col(sIdx).col(xIdx);
        // Apply remaining mortality
        newN_axsp.col(pIdx).col(sIdx).col(xIdx) *= exp(-M_xsp(xIdx,sIdx,pIdx)/2);
      } // END x loop
    } // END p loop
  } // END s loop

  Z_axsp = log(N_axsp / newN_axsp);


  for( int sIdx = 0; sIdx < nS; sIdx ++ )
  {
    for( int pIdx = 0; pIdx < nP; pIdx++ )
    {
      for( int fIdx = 0; fIdx < nF; fIdx++ )
      {
        F_spf(sIdx,pIdx,fIdx) = C_axspf(A_s(sIdx)-1,nX-1,sIdx,pIdx,fIdx) * Z_axsp(A_s(sIdx)-1,nX-1,sIdx,pIdx);
        F_spf(sIdx,pIdx,fIdx) /= N_axsp(A_s(sIdx)-1,nX-1,sIdx,pIdx) / sel_axspf(A_s(sIdx)-1,nX - 1,sIdx,pIdx,fIdx)/(1 - exp(-Z_axsp(A_s(sIdx)-1,nX-1,sIdx,pIdx)));
      }
    }
  }

  
  return(newN_axsp);
}


// calcLogistNormLikelihood()
// Calculates the logistic normal likelihood for compositional data.
// Automatically accumulates proportions below a given threshold. Takes
// cumulative sum of squared resids and number of classes for which
// residuals are calculated as pointers, so that these can be accumulated
// across years for conditional MLE of variance. 
// Will extend to correlated residuals and time-varying variance later.
// inputs:    yObs    = Type vector of observed compositions (samples or proportions)
//            pPred   = Type vector of predicted parameters (true class proportions)
//            minProp = minimum proportion threshold to accumulate classes above
//            etaSumSq= cumulative sum of squared 0-mean resids
//            nResids = cumulative sum of composition classes (post-accumulation)
// outputs:   resids = vector of resids (accumulated to match bins >= minProp)
// Usage:     For computing the likelihood of observed compositional data
// Source:    S. D. N. Johnson
// Reference: Schnute and Haigh, 2007; Francis, 2014
template<class Type>
vector<Type> calcLogistNormLikelihood(  vector<Type>& yObs, 
                                        vector<Type>& pPred,
                                        Type minProp,
                                        Type& etaSumSq,
                                        Type& nResids )

{
  // Get size of vector 
  int nX = yObs.size();

  // Normalise the observed samples in case they are numbers
  // and not proportions
  yObs /= yObs.sum();
  pPred /= pPred.sum();

  // Create vector of residuals to return
  vector<Type> resids(nX);
  vector<Type> aboveInd(nX);
  resids.setZero();
  aboveInd.setZero();

  // Count number of observations less
  // than minProp
  int nAbove = 0;
  for( int x = 0; x < nX; x++)
    if(yObs(x) > minProp)
    {
      nAbove++;
      aboveInd(x) = 1;
    }

  // Now loop and fill
  vector<Type> newObs(nAbove);
  vector<Type> newPred(nAbove);
  newObs.setZero();
  newPred.setZero();
  int k = 0;
  for( int x = 0; x < nX; x++)
  {
    // accumulate observations
    newObs(k) += yObs(x);
    newPred(k) += pPred(x);

    // Increment k if we reach a bin
    // with higher than minProp observations
    // and we aren't yet in the last bin
    if(yObs(x) >= minProp & k < nAbove - 1)
      k++;    
  }

  // Create a residual vector
  vector<Type> res(nAbove);
  res.setZero(); 
  
  for(int j = 0; j < nAbove; j ++)
    if( newObs(j) > 0 &  newPred(j) > 0)
      res(j) = log(newObs(j)) - log(newPred(j));

  // Calculate mean residual
  Type meanRes = res.sum()/nAbove;
  // centre residuals
  res -= meanRes;


  // Now loop again, and fill in
  // the residuals vector
  // for plotting later
  k = 0;
  for( int x =0; x < nX; x++)
    if( aboveInd(x) == 1)
    {
      resids(x) = res(k);
      k++;
    }

  
  // Now add squared resids to etaSumSq and nRes to nResids
  etaSumSq  += (res*res).sum();
  nResids   += nAbove;

  return(resids);
} // end calcLogistNormLikelihood()

// calcCorrLogistNormLikelihood()
// Calculates the logistic normal likelihood for compositional data.
// Automatically accumulates proportions below a given threshold. Takes
// cumulative sum of squared resids and number of classes for which
// residuals are calculated as pointers, so that these can be accumulated
// across years for conditional MLE of variance. 
// Will extend to correlated residuals and time-varying variance later.
// inputs:    yObs    = Type vector of observed compositions (samples or proportions)
//            pPred   = Type vector of predicted parameters (true class proportions)
//            minProp = minimum proportion threshold to accumulate classes above
//            etaSumSq= cumulative sum of squared resids (possibly correlated)
//            nResids = cumulative sum of composition classes (post-accumulation)
//            meanSampSize = mean size of annual composition samples
//            compLikeFun = switch for compositional likelihood from Francis 2014
//                            0 => no correlation (AR0, LN1)
//                            1 => AR1 model (LN2)
//                            2 => AR2 model (LN3)
//                            3 => ARMA model (LN3m)
//            phi1    = Correlation parameter, depends on compLikeFun
//            psi     = Correlation parameter, depends on compLikeFun
//            compObsNLL = compositional observation NLL contribution (external var)
// outputs:   resids = vector of standardised resids (accumulated to match bins >= minProp)
// Usage:     For computing the likelihood of observed compositional data
// Source:    S. D. N. Johnson
// Reference: Schnute and Haigh, 2007; Francis, 2014
template<class Type>
vector<Type> calcCorrLogistNormLikelihood(  vector<Type>   yObs, 
                                            vector<Type>   pPred,
                                            Type           minProp,
                                            Type&          etaSumSq,
                                            Type&          nResids,
                                            Type           meanSampSize,
                                            matrix<Type>   Corr_aa,
                                            Type&          compObsNLL)
{

  // NOTE: some of the notation departs from
  // that defined in the model header. This is
  // to keep it close to Francis 2014 while I work
  // out the model

  // Get size of vector 
  int nX = yObs.size();

  // Get size of this sample
  Type thisSampSize = yObs.sum();

  // Calculate this year's weighting
  Type Wy = sqrt(meanSampSize / thisSampSize);

  // Normalise the samples and predictions in case they are numbers
  // and not proportions
  yObs /= yObs.sum();
  pPred /= pPred.sum();

  // Create vector of residuals to return
  vector<Type> resids(nX);
  vector<Type> aboveInd(nX);
  resids.setZero();
  aboveInd.setZero();

  // Count number of observations less
  // than minProp
  int nAbove = 0;
  for( int x = 0; x < nX; x++)
    if(yObs(x) > minProp)
    {
      nAbove++;
      aboveInd(x) = 1;
    }

  vector<int> idxAbove(nAbove);
  int k = 0;
  for( int x =0; x < nX; x++ )
    if( aboveInd(x) == 1)
    {
      idxAbove(k) = x;
      k++;
    }

  // Now loop and fill
  vector<Type> newObs(nAbove);
  vector<Type> newPred(nAbove);
  newObs.setZero();
  newPred.setZero();

  // Tail compression 
  // go up from the left
  k = 0;
  for( int x = 0; x < nX; x++)
  {
    // accumulate observations
    newObs(k) += yObs(x);
    newPred(k) += pPred(x);

    // Increment k if we reach a bin
    // with higher than minProp observations
    // and we aren't yet in the last bin
    if( (yObs(x) > minProp) & (k < nAbove - 1) )
      k++;    
  }


  // OK, now we create a nAbove x nAbove correlation matrix Corr
  // from the above correlation vector, the same
  // size matrix V, and the multiplier matrix K
  matrix<Type> Corr(nAbove,nAbove);
  matrix<Type> V(nAbove,nAbove);
  matrix<Type> Vinv(nAbove,nAbove);
  matrix<Type> Gamma(nAbove,nAbove);
  matrix<Type> K(nAbove-1,nAbove);
  matrix<Type> F(nAbove-1,nAbove);
  matrix<Type> H(nAbove-1,nAbove-1);
  matrix<Type> Hinv(nAbove-1,nAbove-1);


  // Fill with zeroes
  Corr.setZero();
  V.setZero();
  Vinv.setZero();
  K.setZero();
  F.setZero();
  H.fill(1);

  // Get submatrix of Corr_aa
  for( int k = 0; k < nAbove; k++ )
  {
    Corr(k,k) += Corr_aa(idxAbove(k),idxAbove(k));
    for( int kk = k + 1; kk < nAbove; kk++ )
    {
      Corr(k,kk) += Corr_aa(idxAbove(k),idxAbove(kk));
      Corr(kk,k) += Corr_aa(idxAbove(kk),idxAbove(k));
    }
  }
  

  // Now we fill
  for( int rIdx = 0; rIdx < nAbove; rIdx ++ )
  {

    // Now fill in K
    if( rIdx < nAbove-1 )
    {
      K( rIdx, rIdx ) = 1;
      F( rIdx, rIdx ) = 1;
      H( rIdx, rIdx) += 1;
    }

    if( rIdx == nAbove-1 )
    {
      K.col(rIdx).fill(-1);
      F.col(rIdx).fill(1);
    }

    
  }

  // Generate V and its inverse
  V = K * Corr * K.transpose();
  Type logdetV = 0;
  Vinv = atomic::matinvpd( V, logdetV );
  Hinv = atomic::matinv( H );

  Gamma = F.transpose() * Hinv * V * Hinv * F;

  // Calculate w_b, independent resids
  matrix<Type> w_b(nAbove-1,1);
  w_b.setZero();
  for( int vIdx = 0; vIdx < nAbove-1; vIdx ++)
    if( newObs(vIdx) > 0 & newPred(vIdx) > 0)
      w_b(vIdx,0) = log( newObs(vIdx) / newObs(nAbove-1) ) - 
                  log( newPred(vIdx) / newPred(nAbove-1) );

  
  // Now add correlated squared resids to etaSumSq and nRes to nResids
  matrix<Type> tmpeta(1,1);
  tmpeta = (w_b.transpose() * Vinv * w_b) / square(Wy); 
  etaSumSq  += tmpeta(0,0);
  nResids   += nAbove - 1;

  compObsNLL += 0.5 * logdetV + (nAbove-1) * log(Wy); // + log( newObs.segment(0,nAbove -1) ).sum();

  // Now expand residuals back out
  // Now loop again, and fill in
  // the residuals vector
  // for plotting later
  vector<Type> res(nAbove);
  res.setZero();
  Type gmObs  = 1;
  Type gmPred = 1;
  gmObs  =  exp( log(newObs).sum()  / nAbove  );
  gmPred =  exp( log(newPred).sum() / nAbove  );
  // gmPred =  pow(newPred.prod(),1/nAbove);

  for( int k = 0; k < nAbove; k++ )
  {
    res    =  log( newObs / gmObs ) - 
              log( newPred / gmPred );  
  }
  
  k = 0;
  for( int x =0; x < nX; x++)
    if( aboveInd(x) == 1)
    {
      resids(x)   = res(k) / Wy / sqrt(Gamma(k,k)) ;
      k++;
    }

  return(resids);
}  // END calcCorrLogistNormLikelihood()



// solveBaranov_spfx()
// Newton-Rhapson solver for Baranov catch equation for a multi-stock
// multi-species, age and sex structured population fished by multiple
// fleets
// inputs:    nIter = number of NR iterations
//            Bstep = fraction of NR step (Jacobian) to take at each iteration
//            A_s = number of age classes by species
//            C_spf = Catch
//            M_xsp = natural mortality
//            B_aspx = Biomass at age/sex
//            B_spx = biomass for each sex
//            sel_axspf = selectivity for each age/sex in each fleet
//            & Z = total mortality (external variable)
//            & F = Fishing mortality (external variable)
// returns:   C_xaspf, NR solver estimate of catch at sex/age/spec/pop/fleet
// Side-effs: variables passed as Z, F overwritten with total, fishing mortality
// Author:    Modified by S. D. N. Johnson from S. Rossi and S. P. Cox
template<class Type>
array<Type> solveBaranov_spf( int   nIter,
                              Type  Bstep,
                              vector<int>  A_s,         // number of age classes by species
                              array<Type>  C_spf,       // Total observed catch
                              array<Type>  M_xsp,       // Mortality rate
                              array<Type>  B_axsp,      // Biomass at age/sex
                              array<Type>  vB_axspf,    // vuln biomass at age
                              array<Type>  vB_spf,      // vuln biomass in each fleet
                              array<Type>  sel_axspf,   // selectivity at age/sex
                              array<Type>& Z_axsp,      // total mortality at age/se
                              array<Type>& F_spf,       // fleet F
                              array<Type>  N_axsp )     // Numbers at age
{
  int nS = C_spf.dim(0);
  int nP = C_spf.dim(1);
  int nF = C_spf.dim(2);
  int nA = B_axsp.dim(0);
  int nX = B_axsp.dim(1);

  array<Type> f_spf(nS,nP,nF);              // Function value
  array<Type> J_spf(nS,nP,nF);              // Jacobian
  array<Type> newZ_axsp(nA,nX,nS,nP);       // Updated Z
  array<Type> tmpZ_axsp(nA,nX,nS,nP);       // Updated Z
  array<Type> tmp_spf(nS,nP,nF);            // predicted catch given F
  array<Type> tmp_axspf(nA,nX,nS,nP,nF);    // predicted catch-at-age/sex given F

  newZ_axsp.setZero();
  tmpZ_axsp.setZero();
  f_spf.setZero();
  J_spf.setZero();
  tmp_spf.setZero();
  tmp_axspf.setZero();

  array<Type> newN_axsp(nA,nX,nS,nP);       // End of time step N
  newN_axsp.setZero();
  array<Type> C_axspf(nA,nX,nS,nP,nF);      // Array to hold catch at age in numbers
  C_axspf.setZero();

  // // Initialise F using Pope's approximation
  // newN_axsp = calcPopesApprox(  N_axsp,       // Numbers at age/spec/pop/sex
  //                               sel_axspf,    // Selectivity-at-age (spec/pop/fleet/sex)
  //                               C_spf,        // Catch in biomass units
  //                               M_xsp,        // Natural mortality
  //                               A_s,          // number of age classes by species
  //                               wt_axsp,      // Mean weight at age
  //                               vB_axspf,     // Vulnerable biomass at age/spec/pop/fleet
  //                               vB_spf,       // Vulnerable biomass at spec/pop/fleet
  //                               F_spf,        // Fishing mortality for each spec/pop/fleet  
  //                               Z_axsp,       // Total mortality at age/spec/pop/sex
  //                               tmp_axspf,    // Estimated catch in biomass for sex/age/spec/pop/fleet
  //                               C_axspf );    // Catch in numbers 

  F_spf = C_spf / vB_spf;

  // Initial approximation of F
  for( int s = 0; s < nS; s++ )
    for( int p = 0; p < nP; p++)
    {
      for( int x = 0; x < nX; x++)
        newZ_axsp.col(p).col(s).col(x).fill(M_xsp(x,s,p));

      for( int f = 0; f < nF; f++ )
      {

        for( int x = 0; x < nX; x++)
          newZ_axsp.col(p).col(s).col(x) +=  F_spf(s,p,f) * sel_axspf.col(f).col(p).col(s).col(x);
      }
    }
  
  // Only run this calculation if nIter > 0
  if( nIter > 0)
  {
    // Refine F
    for( int i=0; i<nIter; i++ )
    {
      // Now reset newZ
      tmp_spf.setZero();
      tmp_axspf.setZero();
      J_spf.setZero();
      f_spf.setZero();
      tmpZ_axsp.setZero();


      // Reset objective function and Z
      f_spf += C_spf;


      // Calculate predicted catch and Jacobian
      for( int s = 0; s < nS; s++ )
        for( int p = 0; p < nP; p++ )
        {
          for( int x = 0; x < nX; x++ )
          {
            for( int a = 0; a < A_s(s); a++ )
            {
              for(int f = 0; f < nF; f++ )  
              {  
                tmp_axspf(a,x,s,p,f) += B_axsp(a,x,s,p) * (1 - exp(-newZ_axsp(a,x,s,p)) ) * sel_axspf(a,x,s,p,f) *  F_spf(s,p,f) / newZ_axsp(a,x,s,p);
                tmp_spf(s,p,f) += tmp_axspf(a,x,s,p,f);

                // Calculate jacobian
                Type tmpJ1  = vB_axspf(a,x,s,p,f) / newZ_axsp(a,x,s,p);
                Type tmpJ2  = F_spf(s,p,f) * sel_axspf(a,x,s,p,f) * exp( - newZ_axsp(a,x,s,p));
                Type tmpJ3  = (1 - exp( -newZ_axsp(a,x,s,p)));
                Type tmpJ4  = (newZ_axsp(a,x,s,p) - F_spf(s,p,f) * sel_axspf(a,x,s,p,f))/newZ_axsp(a,x,s,p);

                J_spf(s,p,f) -= tmpJ1 * ( tmpJ2 + tmpJ3 * tmpJ4);
              }
            }
          }
        }

      // Subtract predicted catch
      f_spf -= tmp_spf;

      // Updated fishing mortality
      F_spf -= Bstep * f_spf / J_spf;

      newZ_axsp.setZero();

      // Updated total mortality
      for( int s = 0; s < nS; s++ )
        for( int p = 0; p < nP; p++ )
          for( int x = 0; x < nX; x++ )
          {
            newZ_axsp.col(p).col(s).col(x).fill(M_xsp(x,s,p));
            for( int f= 0 ; f < nF; f ++)
              newZ_axsp.col(p).col(s).col(x) += sel_axspf.col(f).col(p).col(s).col(x) * F_spf(s,p,f) ;  
          }

    }  // end i

    // Now save the tmpZ_axsp array out to the Z_axsp input
    Z_axsp = newZ_axsp;
  }


  return( tmp_axspf );

}  // end solveBaranov_spfx()

// solveBaranov_spfx()
// Newton-Rhapson solver for Baranov catch equation for a multi-stock
// multi-species, age and sex structured population fished by multiple
// fleets
// inputs:    nIter = number of NR iterations
//            Bstep = fraction of NR step (Jacobian) to take at each iteration
//            A_s = number of age classes by species
//            C_spf = Catch
//            M_xsp = natural mortality
//            B_aspx = Biomass at age/sex
//            B_spx = biomass for each sex
//            sel_axspf = selectivity for each age/sex in each fleet
//            & Z = total mortality (external variable)
//            & F = Fishing mortality (external variable)
// returns:   C_xaspf, NR solver estimate of catch at sex/age/spec/pop/fleet
// Side-effs: variables passed as Z, F overwritten with total, fishing mortality
// Author:    Modified by S. D. N. Johnson from S. Rossi and S. P. Cox
template<class Type>
array<Type> solveBaranovEff_spf(  int   nIter,
                                  Type  Bstep,
                                  vector<int>  A_s,         // number of age classes by species
                                  array<Type>  C_spf,       // Total observed catch
                                  array<Type>  M_xsp,       // Mortality rate
                                  array<Type>  B_axsp,      // Biomass at age/sex
                                  array<Type>  vB_axspf,    // vuln biomass at age
                                  array<Type>  vB_spf,      // vuln biomass in each fleet
                                  array<Type>  sel_axspf,   // selectivity at age/sex
                                  array<Type>& Z_axsp,      // total mortality at age/se
                                  array<Type>& F_spf,       // fleet F
                                  array<Type>& E_pf,        // Fleet effort on stock p
                                  array<Type>& q_spf,       // Fleet catchability for species/stock
                                  array<Type>  N_axsp )     // Numbers at age
{
  int nS = C_spf.dim(0);
  int nP = C_spf.dim(1);
  int nF = C_spf.dim(2);
  int nA = B_axsp.dim(0);
  int nX = B_axsp.dim(1);

  array<Type> f_spf(nS,nP,nF);              // Function value
  array<Type> J_spf(nS,nP,nF);              // Jacobian
  array<Type> newZ_axsp(nA,nX,nS,nP);       // Updated Z
  array<Type> tmpZ_axsp(nA,nX,nS,nP);       // Updated Z
  array<Type> tmp_spf(nS,nP,nF);            // predicted catch given F
  array<Type> tmp_axspf(nA,nX,nS,nP,nF);    // predicted catch-at-age/sex given F

  newZ_axsp.setZero();
  tmpZ_axsp.setZero();
  f_spf.setZero();
  J_spf.setZero();
  tmp_spf.setZero();
  tmp_axspf.setZero();

  array<Type> newN_axsp(nA,nX,nS,nP);       // End of time step N
  newN_axsp.setZero();
  array<Type> C_axspf(nA,nX,nS,nP,nF);      // Array to hold catch at age in numbers
  C_axspf.setZero();

  // // Initialise F using Pope's approximation
  // newN_axsp = calcPopesApprox(  N_axsp,       // Numbers at age/spec/pop/sex
  //                               sel_axspf,    // Selectivity-at-age (spec/pop/fleet/sex)
  //                               C_spf,        // Catch in biomass units
  //                               M_xsp,        // Natural mortality
  //                               A_s,          // number of age classes by species
  //                               wt_axsp,      // Mean weight at age
  //                               vB_axspf,     // Vulnerable biomass at age/spec/pop/fleet
  //                               vB_spf,       // Vulnerable biomass at spec/pop/fleet
  //                               F_spf,        // Fishing mortality for each spec/pop/fleet  
  //                               Z_axsp,       // Total mortality at age/spec/pop/sex
  //                               tmp_axspf,    // Estimated catch in biomass for sex/age/spec/pop/fleet
  //                               C_axspf );    // Catch in numbers 

  F_spf = C_spf / vB_spf;


  // Initial approximation of F
  for( int s = 0; s < nS; s++ )
  {
    q_spf.col(s) = F_spf.col(s) / E_pf;
    for( int p = 0; p < nP; p++)
    {
      for( int x = 0; x < nX; x++)
        newZ_axsp.col(p).col(s).col(x).fill(M_xsp(x,s,p));

      for( int f = 0; f < nF; f++ )
      {

        for( int x = 0; x < nX; x++)
          newZ_axsp.col(p).col(s).col(x) +=  F_spf(s,p,f) * sel_axspf.col(f).col(p).col(s).col(x);
      }
    }
  }
  
  // Only run this calculation if nIter > 0
  if( nIter > 0)
  {
    // Refine F
    for( int i=0; i<nIter; i++ )
    {
      // Now reset newZ
      tmp_spf.setZero();
      tmp_axspf.setZero();
      J_spf.setZero();
      f_spf.setZero();
      tmpZ_axsp.setZero();


      // Reset objective function and Z
      f_spf += C_spf;


      // Calculate predicted catch and Jacobian
      for( int s = 0; s < nS; s++ )
        for( int p = 0; p < nP; p++ )
        {
          for( int x = 0; x < nX; x++ )
          {
            for( int a = 0; a < A_s(s); a++ )
            {
              for(int f = 0; f < nF; f++ )  
              {  
                tmp_axspf(a,x,s,p,f) += B_axsp(a,x,s,p) * (1 - exp(-newZ_axsp(a,x,s,p)) ) * sel_axspf(a,x,s,p,f) *  F_spf(s,p,f) / newZ_axsp(a,x,s,p);
                tmp_spf(s,p,f) += tmp_axspf(a,x,s,p,f);

                // Calculate jacobian
                Type tmpJ1  = vB_axspf(a,x,s,p,f) / newZ_axsp(a,x,s,p);
                Type tmpJ2  = F_spf(s,p,f) * sel_axspf(a,x,s,p,f) * exp( - newZ_axsp(a,x,s,p));
                Type tmpJ3  = (1 - exp( -newZ_axsp(a,x,s,p)));
                Type tmpJ4  = (newZ_axsp(a,x,s,p) - F_spf(s,p,f) * sel_axspf(a,x,s,p,f))/newZ_axsp(a,x,s,p);

                J_spf(s,p,f) -= tmpJ1 * ( tmpJ2 + tmpJ3 * tmpJ4);
              }
            }
          }
        }

      // Subtract predicted catch
      f_spf -= tmp_spf;

      // Updated fishing mortality
      F_spf -= Bstep * f_spf / J_spf;

      newZ_axsp.setZero();

      // Updated total mortality
      for( int s = 0; s < nS; s++ )
      {
        q_spf.col(s) = F_spf.col(s) / E_pf;
        for( int p = 0; p < nP; p++ )
          for( int x = 0; x < nX; x++ )
          {
            newZ_axsp.col(p).col(s).col(x).fill(M_xsp(x,s,p));
            for( int f= 0 ; f < nF; f ++)
              newZ_axsp.col(p).col(s).col(x) += sel_axspf.col(f).col(p).col(s).col(x) * F_spf(s,p,f) ;  
          }
      }

    }  // end i

    // Now save the tmpZ_axsp array out to the Z_axsp input
    Z_axsp = newZ_axsp;
  }


  return( tmp_axspf );

}  // end solveBaranov_spfx()


// objective function
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Call namespaces //
  using namespace density;

  /* ----------------- Data Section ------------------------------- */
  // Data Structures
  DATA_ARRAY(I_spft);                   // CPUE data
  DATA_ARRAY(C_spft);                   // Catch data (biomass)
  DATA_ARRAY(D_spft);                   // Discard data (biomass)
  DATA_ARRAY(E_pft);                    // Mean effective effort by stock/fleet/time
  DATA_ARRAY(age_axspft);               // Age observations for each population and fleet over time (-1 missing)
  DATA_ARRAY(len_lxspft);               // Length observations for each population and fleet over time (-1 missing)
  // DATA_ARRAY(ALFreq_laxspft);           // Length observations for each population and fleet over time (-1 missing)
  DATA_VECTOR(lenBinMids_l);            // Midpoints of length bins
  DATA_SCALAR(lenBinWidth);             // Width of length bins (UB - LB)
  DATA_IVECTOR(group_f);                // Fleet group (0=comm, 1=HSAss,2=Synoptic)
  DATA_IVECTOR(A_s);                    // +Group age by species
  DATA_IVECTOR(L_s);                    // Max length bin by species
  DATA_IVECTOR(minA_s);                 // Min age class by species
  DATA_IVECTOR(minL_s);                 // Min length bin by species
  DATA_IVECTOR(lenD_s);                 // Discard length by species
  DATA_INTEGER(nX);                     // Number of sex classes
  DATA_IARRAY(age_table);               // age data in table format
  DATA_IARRAY(len_table);               // length data in table format
  DATA_IARRAY(aal_table);               // age-at-length data in table format


  // Model dimensions
  int nS = I_spft.dim(0);               // No. of species
  int nP = I_spft.dim(1);               // No. of stocks in each species
  int nF = I_spft.dim(2);               // No. of fleets (surveys + trawl)  
  int nT = I_spft.dim(3);               // No of time steps
  int nL = len_lxspft.dim(0);           // max no of length bins (for creating state arrays)
  int nA = age_axspft.dim(0);           // max no of age bins (for creating state arrays)
  int nAgeObs = age_table.dim(0);       // Number of unique age comp observations
  int nLenObs = len_table.dim(0);       // Number of unique length comp observations
  int nAgeLenObs = aal_table.dim(0);    // Number of age-at-length observations


  // Model switches - fill these in when we know what we want
  DATA_IVECTOR(swRinit_s);              // species fished initialisation switch (0 == unfished, 1 == fished)
  DATA_INTEGER(parSwitch);              // parallel accumulator switch for faster run time
  DATA_ARRAY(calcIndex_spf);            // Switches to calculate indices
  DATA_IVECTOR(tvSelFleets);            // Switch for fleet to have time varying selectivity (0 = off, 1 = on)
  DATA_IVECTOR(tvqFleets);              // Switch for fleet to have time varying catchability (0 = off, 1 = on)
  DATA_IVECTOR(logitqFleets);           // Switch for fleet to have time logistic increase in catchability (0 = off, 1 = on)
  DATA_IVECTOR(regFfleets);             // Switch for each fleet to have Fs regularised (Fbar penalised)
  DATA_IVECTOR(commSurv_f);             // Switch for fleet type (commercial == 1 or survey == 0)
  DATA_IVECTOR(solveQ_f);               // switch for each fleet to have q_spft solved from effort and fishing mort.
  DATA_VECTOR(idxLikeWt_g);             // Scalar modifying the weight of the age likelihood by fleet group
  DATA_VECTOR(ageLikeWt_g);             // Scalar modifying the weight of the age likelihood by fleet group
  DATA_VECTOR(lenLikeWt_g);             // Scalar modifying the weight of the length likelihood by fleet group
  DATA_SCALAR(growthLikeWt);            // Scalar modifying the weight of the RAL growth model likelihood
  DATA_IVECTOR(tFirstRecDev_s);         // Initial rec devs for each species
  DATA_IVECTOR(tLastRecDev_s);          // Last recruitment devs for each species
  DATA_SCALAR(minAgeProp);              // Minimum observed proportion in age comps
  DATA_SCALAR(minLenProp);              // Minimum observed proportion in length comps
  DATA_SCALAR(minPropAAL);              // Minimum observed proportion in length comps
  DATA_STRING(matX);                    // Age/Len switch for maturity model
  DATA_STRING(selX);                    // Age/Len switch for selectivity model
  DATA_STRING(lenComps);                // Switch for computation of expected length comps
  DATA_SCALAR(lambdaB0);                // lambda scalar for exponential prior on B0
  DATA_SCALAR(lambdaPropF);             // lambda scalar for prior on propF
  DATA_ARRAY(minTimeIdx_spf);           // array of earliest times for each index
  DATA_INTEGER(nBaranovIter);           // number of baranov steps
  DATA_SCALAR(lambdaBaranovStep);       // fractional step size for Newton-Rhapson Baranov iteration
  DATA_STRING(calcFmethod);             // Indicator for F calculation method ("popesApprox" or "solveBaranov")
  DATA_INTEGER(postFitSR);              // Switch for post-fitting SR model (0 == NO, 1 == YES)
  DATA_ARRAY(calcStockSelDevs_spf);     // Switch for calculating stock-specific selectivity parameters
  DATA_ARRAY(calcStockQDevs_spf);       // Switch for calculating stock-specific catchability parameters
  DATA_INTEGER(boundRecDevs);           // Switch for bounding recruitment deviations between +/- 5 sds
  DATA_STRING(recruitVariance);         // Character-switch for recruitment variance model
  DATA_STRING(recOption);               // Character-switch for recruitment variance model
  DATA_INTEGER(debugMode);              // 1 = debug mode on, return all arrays, 0 = fit mode, return only relevant arrays
  DATA_IVECTOR(condMLEq_f);             // calculate conditional MLE of q for fleet f; 1 => YES, 0 => NO (estimate freely or fix)
  DATA_SCALAR(fixVonB);                 // Fix vonB pars? 0 == No/estimate, 1 == use inputs below
  DATA_VECTOR(inputL1_s);               // input L1 values
  DATA_ARRAY(inputL2_xsp);              // input L2 values
  DATA_ARRAY(inputVonK_xsp);            // input VonK values

  // Fixed values
  DATA_IVECTOR(A1_s);                   // Age at which L1_s is estimated
  DATA_IVECTOR(A2_s);                   // Age at which L2_s is estimated

  // Compositional Likelihood values
  // Compositional likelihood inputs
  DATA_INTEGER(compLikeFun);            // Compositional data likelihood function (correlation)
  DATA_ARRAY(meanAgeSampSize_spf);      // Mean sample sizes for sex-species-stock-specific compositional data
  DATA_ARRAY(meanLenSampSize_spf);      // Mean sample sizes for sex-species-stock-specific compositional data



  /* -------------------------- Parameter Section ------------------------ */
  // Leading Parameters
  // Biological
  PARAMETER_ARRAY(lnB0_sp);             // Biomass at unfished for ind. stocks
  PARAMETER_ARRAY(lnRbar_sp);           // Average recruitment - for post-fitting SR model
  PARAMETER(logitSteep);                // pop specific steepness
  PARAMETER(lnM);                       // Complex M 

  // Growth
  // PARAMETER_VECTOR(thetaL2_s);          // Logit-scale L2 parameter for species s
  // PARAMETER_VECTOR(thetaL1_s);          // Logit-scale L1 parameter for species s
  // PARAMETER_VECTOR(thetaVonK_s);        // Logit-scale vonK parameter for species s
  PARAMETER_VECTOR(lnL2step_s);         // log-scale L2-step parameter for species s
  PARAMETER_VECTOR(lnL1_s);             // log-scale L1 parameter for species s
  PARAMETER_VECTOR(lnVonK_s);           // log-scale vonK parameter for species s
  PARAMETER_ARRAY(deltaL2_sp);          // Stock-deviations for L2
  PARAMETER_ARRAY(deltaVonK_sp);        // Stock deviations for VonK
  PARAMETER_ARRAY(deltaL2_xs);          // sex deviations for L2
  PARAMETER_ARRAY(deltaVonK_xs);        // sex deviations for vonK
  PARAMETER_VECTOR(lnsigmaLa_s);        // species-specific individual growth SD intercept
  PARAMETER_VECTOR(sigmaLb_s);          // species-specific individual growth SD slope
  PARAMETER_VECTOR(LWa_s);              // species L-W conversion a par (units conv: cm -> kt)
  PARAMETER_VECTOR(LWb_s);              // species L-W conversion b par (length to volume)

  // Maturity
  PARAMETER_VECTOR(xMat50_s);           // x (age/len) at 50% maturity
  PARAMETER_VECTOR(xMat95_s);           // x (age/len) at 95% maturity

  // Observation model
  // survey //
  PARAMETER_ARRAY(lnqComm_spf);         // Species-stock-fleet specific catchability
  PARAMETER_ARRAY(lnqSurv_sf);          // fleet-species specific mean catchability
  PARAMETER_VECTOR(lntq50_vec);         // fleet specific time at 50% catchability (improving fishing)
  PARAMETER_VECTOR(lntq95_vec);         // fleet specific time at 95% catchability (improving fishing)
  PARAMETER_VECTOR(lntauObs_vec);       // fleet-species-stock specific observation error variance (vector for avoiding missing observation fleets)
  // Fishery model
  // selectivity by fleet/species/
  PARAMETER_ARRAY(thetaSelAlpha_sg);    // Selectivity Alpha parameter by species/fleetgroup
  PARAMETER_ARRAY(thetaSelBeta_sg);     // Selectivity Beta parameter by species/fleetgroup
  // Now use deviations to get the fleet/stock specific values
  PARAMETER_VECTOR(epsSelAlphaspf_vec); // Selectivivity devs for stocks
  PARAMETER_VECTOR(epsSelBetaspf_vec);  // Selectivivity devs for stocks

  // Count fleet groups
  int nGroups = thetaSelAlpha_sg.dim(1);
  int nComm = lnqComm_spf.dim(2);
  int nSurv = lnqSurv_sf.dim(1);

  // Fishing mortality
  PARAMETER_VECTOR(lntauD_f);           // Discards observation SD by fleet
  PARAMETER_ARRAY(lnqF_spf);            // Mean scaling coefficent for transformin Effort to Fishing mortality
  PARAMETER(lntauqFdev);                // SD on deltalnqF devs
  PARAMETER_VECTOR(deltalnqFspft_vec);  // Annual deviation in scaling coefficent for transforming Effort to Fishing mortality

  // Priors
  // selectivity //
  PARAMETER_ARRAY(lnsigmaxSel50_sg);      // SD in x-at-50% selectivity by fleet group
  PARAMETER_ARRAY(lnsigmaxSelStep_sg);    // SD in step from x-at-50% to x-at-95% selectivity by fleet group
  PARAMETER_VECTOR(pmsigmaSel_g);         // IG prior mode sel dev SD
  PARAMETER(IGalphaSel);                  // IG alpha parameter for sel dev var
  PARAMETER_ARRAY(LB_xSelAlpha_sg);       // Lower bound on alpha par for species s and fleetgp g
  PARAMETER_ARRAY(UB_xSelAlpha_sg);       // Upper bound on alpha par for species s and fleetgp g
  PARAMETER(cvSel);                       // Prior SD on ThetaSel values

  // Growth
  // PARAMETER_VECTOR(UBL2_s);             // Upper bound on L2_s
  // PARAMETER_VECTOR(LBL2_s);             // Upper bound on L2_s
  // PARAMETER_VECTOR(UBL1_s);             // Upper bound on L1_s
  // PARAMETER_VECTOR(LBL1_s);             // Upper bound on L1_s
  // PARAMETER_VECTOR(UBVonK_s);           // Upper bound on VonK_s
  // PARAMETER_VECTOR(LBVonK_s);           // Upper bound on VonK_s
  PARAMETER(lnsigmaL2);                 // log SD in L2 deviations (sex and stock) 
  PARAMETER(lnsigmaVonK);               // log SD in VonK deviations (sex and stock) 
  PARAMETER(cvL2);                      // Prior CV on L2
  PARAMETER(cvL1);                      // Prior CV on L1
  PARAMETER(cvVonK);                    // Prior CV on VonK


  // catchability //
  PARAMETER_VECTOR(deltaqspf_vec);      // Stock-specific fleet index catchability
  PARAMETER_ARRAY(lntauqSurv_sf);       // SD of Species-specific fleet group catchability dist
  PARAMETER_VECTOR(pmtauqSurv_f);       // IG prior mode q dev SD
  PARAMETER(IGalphaq);                  // IG alpha parameter for q dev variance
  PARAMETER_VECTOR(mq_f);               // Normal prior mean for catchability
  PARAMETER_VECTOR(sdq_f);              // Normal prior sd for catchability


  // steepness //
  PARAMETER_VECTOR(hBetaPrior);         // steepness beta hyperprior parameters
  PARAMETER_VECTOR(lnsigmah_s);         // species level steepness SD
  PARAMETER(lnsigmah);                  // complex level steepness SD
  PARAMETER(pmsigmah);                  // IG prior mode h dev SD
  PARAMETER(IGalphah);                  // IG alpha parameter for h dev variance
  // Natural Mortality //
  PARAMETER_VECTOR(lnsigmaM_s);         // species M SD
  PARAMETER(lnsigmaM);                  // Assemblage M SD
  PARAMETER(ln_muM);                    // Multispecies assemblage prior mean M
  PARAMETER(sdM);                       // Multispecies assemblage prior M sd
  PARAMETER(pmsigmaM);                  // IG prior mode M dev SD
  PARAMETER(IGalphaM);                  // IG alpha parameter for M dev variance
  // Index observation errors SD
  PARAMETER_VECTOR(pmtauObs_g);         // IG prior mode obs error SD, by group 
  PARAMETER(IGalphaObs);                // IG alpha parameter for observation error variance

  // Random Effects
  // Steepness and mortality //
  PARAMETER_VECTOR(epsSteep_s);         // Species level steepness effect
  PARAMETER_VECTOR(epsM_s);             // Species level natural mortality effect
  PARAMETER_ARRAY(epsSteep_sp);         // stock level steepness effect
  PARAMETER_ARRAY(epsM_sp);             // stock level natural mortality effect
  PARAMETER_ARRAY(epsM_sx);             // sex specific natural mortality effect
  // Recruitment //
  PARAMETER_VECTOR(omegaR_vec);         // species-stock specific recruitment errors 2:nT
  PARAMETER_VECTOR(omegaRinit_vec);     // stock-age specific recruitment initialisation errors
  PARAMETER_ARRAY(lnsigmaR_sp);         // stock recruitment errors sd (sqrt cov matrix diag)
  PARAMETER(pmsigmaR);                  // IG prior mode for sigma R
  PARAMETER(IGalphaR);                  // IG shape parameter for proc error variance
  PARAMETER_ARRAY(logitRgamma_sp);      // species-stock-specific AR1 auto-corr on year effect (omega)
  PARAMETER(sd_omegaRbar);              // omegaRbar regularisation SD
  // Time varying selectivity and catchability //
  PARAMETER_VECTOR(epsSelAlpha_vec);    // random effect on length at 50% selectivity
  PARAMETER_VECTOR(epsSelBeta_vec);     // random effect on length at 95% selectivity
  PARAMETER(lnsigmaSel);                // SD on time varying random effect for selectivity
  PARAMETER_VECTOR(epslnq_vec);         // random effect for time varying catchability parameter
  PARAMETER( lnsigmaepslnq );           // SD on time varying RE for lnq
  
  // Recruitment correlations
  PARAMETER_VECTOR(recCorr_vec);        // vector for cholesky decomposition of recruitment correlation mtx
  PARAMETER_MATRIX(IWmode);             // IW scale matrix for Sigma prior
  PARAMETER(IWnu);                      // IW degrees of freedom for Sigma prior    

  // Priors fishing mortality
  PARAMETER(mF);                        // Prior mean value for mean F
  PARAMETER(sdF);                       // Prior sd for meanF
  
  // Compositional data likelihood correlation matrix parameters
  PARAMETER_ARRAY(logitphi1Age_sf);     // AR1 correlation coefficient, or AR2 parameter
  PARAMETER_ARRAY(logitpsiAge_sf);      // AR2 parameter
  PARAMETER_ARRAY(logitphi1Len_sf);     // AR1 correlation coefficient, or AR2 parameter
  PARAMETER_ARRAY(logitpsiLen_sf);      // AR2 parameter


  /* ------------------------------ Procedure Section --------------------- */
  // Exponentiate leading parameters
  // Biological //
  array<Type>  B0_sp(nS,nP);
  array<Type>  Rbar_sp(nS,nP);
  array<Type>  h_sp(nS,nP);
  array<Type>  M_sp(nS,nP);
  array<Type>  sigmaR_sp(nS,nP);
  // Sex/stock specific
  array<Type>  M_xsp(nX,nS,nP);
  array<Type>  L2_xsp(nX,nS,nP);
  array<Type>  L1_xsp(nX,nS,nP);
  array<Type>  vonK_xsp(nX,nS,nP);
  // Sex specific
  array<Type>  L1_xs(nX,nS);
  array<Type>  L2_xs(nX,nS);
  array<Type>  vonK_xs(nX,nS);
  // Stock specific
  array<Type>  L1_sp(nS,nP);
  array<Type>  L2_sp(nS,nP);
  array<Type>  vonK_sp(nS,nP);

  L2_xsp.setZero();
  L1_xsp.setZero();
  vonK_xsp.setZero();
  L2_xs.setZero();
  L1_xs.setZero();
  vonK_xs.setZero();
  L2_sp.setZero();
  L1_sp.setZero();
  vonK_sp.setZero();

  // Growth model vectors
  vector<Type>  sigmaLa_s(nS);
  vector<Type>  L1_s(nS);
  vector<Type>  L2_s(nS);   
  vector<Type>  vonK_s(nS);   
  // Average growth model devs
  vector<Type>  deltaVonKbar_p(nP);   
  vector<Type>  deltaL2bar_p(nP);   
  deltaVonKbar_p.setZero();
  deltaL2bar_p.setZero();

  // derived variables
  // Stock recruitment //
  array<Type>   R0_sp(nS,nP);            // eqbm recruitment
  array<Type>   phi_sp(nS,nP);           // eqbm SSB per recruit
  array<Type>   reca_sp(nS,nP);          // BH a parameter for pops
  array<Type>   recb_sp(nS,nP);          // BH b parameter for pops
  array<Type>   gammaR_sp(nS,nP);        // AR-1 autocorrelation factor for SR devs

  R0_sp.setZero();
  phi_sp.setZero();
  reca_sp.setZero();
  recb_sp.setZero();
  gammaR_sp.setZero();

  // Growth //
  array<Type>   Wlen_ls(nL,nS);                   // weight-at-length by species
  array<Type>   lenAge_axsp(nA,nX,nS,nP);         // Mean length-at-age by population
  array<Type>   probLenAge_laxsp(nL,nA,nX,nS,nP); // Probability of length-at-age
  array<Type>   meanWtAge_axsp(nA,nX,nS,nP);      // Mean weight-at-age
  array<Type>   probAgeLen_laxspft(nL,nA,nX,nS,nP,nF,nT);   // Probability of age-at-length
  probAgeLen_laxspft.setZero();


  // Maturity //
  array<Type>   matAge_as(nA,nS);       // Proportion mature at age by species
  array<Type>   matAge_asp(nA,nS,nP);   // Proportion mature at age by population
  array<Type>   matLen_ls(nL,nS);       // Proportion mature at Length by species


  // Catchability and index observation errors
  // vector<Type>  qSurv_sf(nF);
  array<Type>   tq50_spf(nS,nP,nF);
  array<Type>   tq95_spf(nS,nP,nF);
  // vector<Type>  tauq_f(nF);
  // vector<Type>  tau2q_f(nF);
  array<Type>   qSurv_sf(nS,nSurv);
  array<Type>   qComm_spf(nS,nP,nComm);
  array<Type>   tauqSurv_sf(nS,nSurv);
  array<Type>   tau2qSurv_sf(nS,nSurv);
  array<Type>   q_spf(nS,nP,nF);
  array<Type>   qF_spf(nS,nP,nF);
  array<Type>   deltaq_spf(nS,nP,nF);
  array<Type>   lntauObs_spf(nS,nP,nF);
  array<Type>   tauObs_spf(nS,nP,nF);
  array<Type>   tau2Obs_spf(nS,nP,nF);
  
  tq50_spf.setZero();
  tq95_spf.setZero();
  qSurv_sf.setZero();
  q_spf.setZero();
  qF_spf.setZero();
  deltaq_spf.setZero();
  tau2qSurv_sf.setZero();
  tauqSurv_sf.setZero();
  // tauq_f.setZero();
  // tau2q_f.setZero();
  lntauObs_spf.setZero();
  tauObs_spf.setZero();
  tau2Obs_spf.setZero();


  // Time-varying catchability
  array<Type>   q_spft(nS,nP,nF,nT);
  // Fishing mortality
  array<Type>   qF_spft(nS,nP,nF,nT);
  array<Type>   F_spft(nS,nP,nF,nT);
  array<Type>   paF_spft(nS,nP,nF,nT);
  array<Type>   totCatch_sp(nS,nP);
  array<Type>   Fbar_spf(nS,nP,nF);
  array<Type>   calcSel_spf(nS,nP,nF);

  // selectivity //
  // Fleet/species
  array<Type>   xSel50_sf(nS,nF);
  array<Type>   xSelStep_sf(nS,nF);
  array<Type>   xSel95_sf(nS,nF);
  // species/fleet group
  array<Type>   xSel50_sg(nS,nGroups);
  array<Type>   xSelStep_sg(nS,nGroups);
  array<Type>   xSel95_sg(nS,nGroups);
  array<Type>   sigmaxSel50_sg(nS,nGroups);
  array<Type>   sigmaxSelStep_sg(nS,nGroups);
  // Fleet/species/Stock
  array<Type>   xSel50_spf(nS,nP,nF);
  array<Type>   xSelStep_spf(nS,nP,nF);
  array<Type>   xSel95_spf(nS,nP,nF);
  array<Type>   epsxSel50_spf(nS,nP,nF);
  array<Type>   epsxSelStep_spf(nS,nP,nF);
  // Now time-varying by stock
  array<Type>   xSel50_spft(nS,nP,nF,nT);
  array<Type>   xSelStep_spft(nS,nP,nF,nT);
  array<Type>   xSel95_spft(nS,nP,nF,nT);
  array<Type>   epsxSel50_spft(nS,nP,nF,nT);
  array<Type>   epsxSelStep_spft(nS,nP,nF,nT);
  array<Type>   epslnq_spft(nS,nP,nF,nT);
  // Zero-initialise
  xSel50_sg.setZero();
  xSel95_sg.setZero();
  xSelStep_sg.setZero();
  xSel50_spf.setZero();
  xSel95_spf.setZero();
  xSelStep_spf.setZero();
  xSel50_spft.setZero();
  xSel95_spft.setZero();
  xSelStep_spft.setZero();
  epsxSel50_spft.setZero();
  epsxSelStep_spft.setZero();
  epslnq_spft.setZero();
  F_spft.setZero();
  paF_spft.setZero();
  Fbar_spf.setZero();
  calcSel_spf.setZero();
  q_spft.fill(1.);
  qF_spft.fill(0);

  // Compositional likelihood correlation matrix
  // parameters
  array<Type> phi1Age_sf(nS,nF);
  array<Type> psiAge_sf(nS,nF);
  array<Type> phi2Age_sf(nS,nF);
  array<Type> phi1Len_sf(nS,nF);
  array<Type> psiLen_sf(nS,nF);
  array<Type> phi2Len_sf(nS,nF);

  phi2Age_sf.fill(0);
  phi2Len_sf.fill(0);
  phi1Age_sf.fill(0);
  phi1Len_sf.fill(0);
  psiAge_sf.fill(.5);
  psiLen_sf.fill(.5);

  array<Type> corrAge_sfa(nS,nF,nA);
  array<Type> corrLen_sfl(nS,nF,nL);
  corrAge_sfa.fill(0);
  corrAge_sfa.col(0) += 1;
  corrLen_sfl.fill(0);
  corrLen_sfl.col(0) += 1;

  // Now fill in the rest: 
  if( compLikeFun > 0)
  {
    // The LN2 function
    if( compLikeFun == 1 )
    {
      for( int f = 0; f < nF; f++ )
      {
        for( int s = 0; s < nS; s++ )
        {
          phi1Age_sf(s,f) = 2 /(1 + exp(-logitphi1Age_sf(s,f))) - 1;
          phi1Len_sf(s,f) = 2 /(1 + exp(-logitphi1Len_sf(s,f))) - 1;
        }
        for( int a = 1; a < nA; a++ )
          corrAge_sfa.col(a).col(f) = corrAge_sfa.col(a-1).col(f) * phi1Age_sf.col(f);
        for( int l = 1; l < nL; l++ )
          corrLen_sfl.col(l).col(f) = corrLen_sfl.col(l-1).col(f) * phi1Len_sf.col(f);
      }
    }  

    // The LN3 function
    if( compLikeFun == 2 )
    { 
      for( int f = 0; f < nF; f++ )
      {
        for( int s = 0; s < nS; s++ )
        {
          phi1Age_sf(s,f) = 4 * 1/(1 + exp(-logitphi1Age_sf(s,f))) - 2;
          psiAge_sf(s,f)  = 1/(1 + exp(-logitpsiAge_sf(s,f)));
          phi2Age_sf(s,f) = -1 + (2 - sqrt(square(phi1Age_sf(s,f))) ) * psiAge_sf(s,f);

          phi1Len_sf(s,f) = 4 * 1/(1 + exp(-logitphi1Len_sf(s,f))) - 2;
          psiLen_sf(s,f)  = 1/(1 + exp(-logitpsiLen_sf(s,f)));
          phi2Len_sf(s,f) = -1 + (2 - sqrt(square(phi1Len_sf(s,f))) ) * psiLen_sf(s,f);
        }
      
        corrAge_sfa.col(1).col(f) = phi1Age_sf.col(f) / phi2Age_sf.col(f);
        corrLen_sfl.col(1).col(f) = phi1Len_sf.col(f) / phi2Len_sf.col(f);
        for( int a = 2; a < nA; a++ )
          corrAge_sfa.col(a).col(f) = phi1Age_sf.col(f) * corrAge_sfa.col(a-1).col(f) + phi2Age_sf.col(f) * corrAge_sfa.col(a-2).col(f);
        for( int l = 2; l < nL; l++)
          corrLen_sfl.col(l).col(f) = phi1Len_sf.col(f) * corrLen_sfl.col(l-1).col(f) + phi2Len_sf.col(f) * corrLen_sfl.col(l-2).col(f);
      }
      
    }

    // LN3m, ARMA model
    if( compLikeFun == 3 )
    {
      for( int f = 0; f < nF; f++ )
      {
        for( int s = 0; s < nS; s++ )
        {
          phi1Age_sf(s,f)    = 2 /(1 + exp(-logitphi1Age_sf(s,f))) - 1;
          psiAge_sf(s,f)    += logitpsiAge_sf(s,f);
          phi1Len_sf(s,f)    = 2 /(1 + exp(-logitphi1Len_sf(s,f))) - 1;
          psiLen_sf(s,f)    += logitpsiLen_sf(s,f);
        }

        vector<Type> sumPhiPsiAge = phi1Age_sf.col(f) + psiAge_sf.col(f);
        vector<Type> sumPhiPsiLen = phi1Len_sf.col(f) + psiLen_sf.col(f);
      
        corrAge_sfa.col(1).col(f) = phi1Age_sf.col(f) + psiAge_sf.col(f) * (1 + (sumPhiPsiAge * sumPhiPsiAge )/(1 - phi1Age_sf.col(f) * phi1Age_sf.col(f) ) );
        corrLen_sfl.col(1).col(f) = phi1Len_sf.col(f) + psiLen_sf.col(f) * (1 + (sumPhiPsiLen * sumPhiPsiLen )/(1 - phi1Len_sf.col(f) * phi1Age_sf.col(f) ) );
      
        for( int a = 2; a < nA; a++ )
          corrAge_sfa.col(a).col(f) = phi1Age_sf.col(f) * corrAge_sfa.col(a-2).col(f);
        for( int l = 2; l < nL; l++ )
          corrLen_sfl.col(l).col(f) = phi1Len_sf.col(f) * corrLen_sfl.col(l-2).col(f);
      }
    }
  }

  // Now make correlation matrices
  array<Type> CorrAge_sfaa(nS,nF,nA,nA);
  array<Type> CorrLen_sfll(nS,nF,nL,nL);
  CorrAge_sfaa.setZero();
  CorrLen_sfll.setZero();
  

  for( int s = 0; s < nS; s++)
    for( int f = 0; f < nF; f++ )
    {
      for( int rIdx = 0; rIdx < nA; rIdx ++ )
      {
        CorrAge_sfaa(s,f,rIdx,rIdx) = 1;
        for( int cIdx = rIdx + 1; cIdx < nA; cIdx ++)
        {
          int idxDiff = cIdx - rIdx;

          CorrAge_sfaa(s,f,rIdx,cIdx) = corrAge_sfa(s,f,idxDiff);   
          CorrAge_sfaa(s,f,cIdx,rIdx) = corrAge_sfa(s,f,idxDiff);

        }
      }

      for( int rIdx = 0; rIdx < nL; rIdx ++ )
      {
        CorrLen_sfll(s,f,rIdx,rIdx) = 1;
        for( int cIdx = rIdx + 1; cIdx < nL; cIdx ++)
        {
          int idxDiff = cIdx - rIdx;

          CorrLen_sfll(s,f,rIdx,cIdx) = corrLen_sfl(s,f,idxDiff);   
          CorrLen_sfll(s,f,cIdx,rIdx) = corrLen_sfl(s,f,idxDiff);

        }
      }    
    }

  // Make a vector of ages/lengths for use in 
  // age and length based quantities later
  vector<Type> age(nA);

  for( int a = 0; a < nA; a++)
    age(a) = a+1;

  // parallel_accumulator<Type> f(this);
  Type f = 0;
  Type joint_nlp = 0.0;
  

  // Recruitment deviations
  // correlation matrix
  matrix<Type> recCorrMat_sp(nS*nP,nS*nP);
  recCorrMat_sp.setZero();
  // AR1 factor
  gammaR_sp = -1 + 1.96 /( 1 + exp (-1 * (logitRgamma_sp) ) );


  // Prior Hyperparameters //
  // Steepness
  vector<Type>  sigmah_s      = exp(lnsigmah_s);
  Type          sigmah        = exp(lnsigmah);
  vector<Type>  h_s           = .21 + Type(0.78) / (1 + exp( -1 * ( logitSteep + sigmah * epsSteep_s ) ));
  Type          h             = .21 + .78 * invlogit(logitSteep);           

  // Natural mortality
  Type          sigmaM        = exp(lnsigmaM);
  vector<Type>  sigmaM_s      = exp(lnsigmaM_s);
  vector<Type>  M_s           = exp(lnM + sigmaM * epsM_s);
  Type          muM           = exp(ln_muM);
  
  // Time varying RE SDs
  Type          sigmaSel      = exp(lnsigmaSel);
  Type          sigmalnq      = exp(lnsigmaepslnq);

  // qFdevs
  Type          tauqFdev      = exp(lntauqFdev);


  // Growth model 
  // deviation SDs
  Type sigmaL2    = exp(lnsigmaL2);
  Type sigmaVonK  = exp(lnsigmaVonK);

  // species level values
  L1_s.setZero();
  L2_s.setZero();
  if( fixVonB == 0)
  {
    L1_s      = exp(lnL1_s);
    L2_s      = exp(lnL1_s) + exp(lnL2step_s);
    // L1_s      = LBL1_s + (UBL1_s - LBL1_s) * ( 1/20 + 17/20/(1 + exp(-thetaL1_s) ) ); 
    // L2_s      = LBL2_s + (UBL2_s - LBL2_s) * ( 1/20 + 17/20/(1 + exp(-thetaL2_s) ) ); 
  }
  if( fixVonB == 1 )
  {
    L1_s = inputL1_s;
    L2_xsp = inputL2_xsp;
    vonK_xsp = inputVonK_xsp;
  }
  sigmaLa_s = exp(lnsigmaLa_s);

  // Fill arrays that hold stock-specific parameters
  B0_sp      = exp(lnB0_sp);
  Rbar_sp    = exp(lnRbar_sp);
  sigmaR_sp  = exp(lnsigmaR_sp);


  int stockGrowthVecIdx   = 0;
  int stockSpecSelVecIdx  = 0;
  int stockSpecQVecIdx    = 0;
  int tauObsIdx           = 0;

  // Exponentiate and build catchability parameters
  // q_f       = exp(lnq_f);
  // tauq_f    = exp(lntauq_f);
  // tau2q_f   = exp(2*lntauq_f);
  qComm_spf     = exp(lnqComm_spf);
  tauqSurv_sf   = exp(lntauqSurv_sf);
  tau2qSurv_sf  = exp(2*lntauqSurv_sf);

  // Create group mean catchabilities for each
  // species
  qSurv_sf = exp(lnqSurv_sf);
  // for( int f =0; f < nF; f++ )
  //   qSurv_sf.col(f) = exp(lnq_f(f) + tauq_f(f) * deltaqSurv_sf.col(f));

  // Calculate propQ_ft over fleets and times
  array<Type> propQ_spft(nS,nP,nF,nT);
  propQ_spft.fill(1);
  int vIdx =0;


  for( int s = 0; s < nS; s++ )
  {
    for( int p = 0; p < nP; p++ )
    {
      
      M_sp(s,p)   = exp(lnM) * exp( sigmaM * epsM_s(s) + sigmaM_s(s) * epsM_sp(s,p) );
      h_sp(s,p)   = .21 + 0.78 / (1 + exp(-1 * (logitSteep + sigmah * epsSteep_s(s) + sigmah_s(s) * epsSteep_sp(s,p) ) ));

      // Stock specific (not species/group) sel pars
      for( int f = 0; f < nF; f++ )
      {
        if(calcStockSelDevs_spf(s,p,f) == 1)
        {
          epsxSel50_spf(s,p,f)    = epsSelAlphaspf_vec(stockSpecSelVecIdx);
          epsxSelStep_spf(s,p,f)  = epsSelBetaspf_vec(stockSpecSelVecIdx);
          stockSpecSelVecIdx++;
        }

        if( calcStockQDevs_spf(s,p,f) == 1)
        {
          deltaq_spf(s,p,f)       = deltaqspf_vec(stockSpecQVecIdx);
          stockSpecQVecIdx++;
        }

        if( C_spft.transpose().col(s).col(p).col(f).sum() > 0 )
        {
          calcSel_spf(s,p,f) = 1.;
          if( calcFmethod == "effortModel" )
            qF_spf(s,p,f) = exp(lnqF_spf(s,p,f));
        }

        if( logitqFleets(f) == 1 )
        {
          for( int s = 0; s < nS; s++ )
            for( int p = 0; p < nP; p++ )
            {
              tq50_spf(s,p,f) = exp(lntq50_vec(vIdx));
              tq95_spf(s,p,f) = exp(lntq95_vec(vIdx));    

              Type tqStep = tq95_spf(s,p,f) - tq50_spf(s,p,f);
              vIdx ++;
              for( int t = 0; t < nT; t++ )
                propQ_spft(s,p,f,t) = .2 + .8 / (1 + exp( -1 * log(19) * (t + 1 - tq50_spf(s,p,f) )/tqStep) );
            }
          
          
        }

        if( calcIndex_spf(s,p,f) == 1 )
        {
          lntauObs_spf(s,p,f) = lntauObs_vec(tauObsIdx);
          tauObsIdx++;
        }
      }

      // And sex specific devs
      for( int x = 0; x < nX; x++)
      {
        M_xsp(x,s,p) = M_sp(s,p) * exp( sigmaM_s(s) * epsM_sx(s,x) );

        // Growth model pars

        L1_xsp(x,s,p)     += L1_s(s);                    

        if( fixVonB == 0 )
        {
          L2_xsp(x,s,p)     += L1_s(s) + 
                              exp( lnL2step_s(s) + sigmaL2 * ( deltaL2_xs(x,s) + deltaL2_sp(s,p) ) );

          vonK_xsp(x,s,p)   += exp( lnVonK_s(s) + sigmaVonK * ( deltaVonK_xs(x,s) + deltaVonK_sp(s,p) ) );
        }

        // Make species/sex versions
        if( p == 0  )
        {
          L1_xs(x,s)    += L1_s(s);                    

          if( fixVonB == 0)
          {
            L2_xs(x,s)    += L1_s(s) + 
                            exp( lnL2step_s(s) + sigmaL2 * ( deltaL2_xs(x,s) ) );

            vonK_xs(x,s)  += exp(lnVonK_s(s) + sigmaVonK * ( deltaVonK_xs(x,s) ) ); 
          }
        }
      }

      // Make species/stock growth pars w/out sex devs
      // Growth model pars
      if( fixVonB == 0)
      {
        L1_sp(s,p)    = L1_s(s);                    

        if( fixVonB == 0 )
        {
          L2_sp(s,p)    = L1_s(s) + exp( lnL2step_s(s) + sigmaL2 * (  deltaL2_sp(s,p) ) );

          vonK_sp(s,p)  = exp(lnVonK_s(s) + sigmaVonK * (  deltaVonK_sp(s,p) ) );
        }
      }
    }
  }


  // Exponentiate species/fleetGroup sel
  xSel50_sg = LB_xSelAlpha_sg + (UB_xSelAlpha_sg - LB_xSelAlpha_sg) * (1/20 + 18 * exp(thetaSelAlpha_sg) / 20 / (1 + exp(thetaSelAlpha_sg)));
  xSelStep_sg = log(19) * (UB_xSelAlpha_sg - LB_xSelAlpha_sg) * (0.005 + 0.5 * exp(thetaSelBeta_sg) / (1 + exp(thetaSelBeta_sg)));;
  xSel95_sg = xSel50_sg + xSelStep_sg;  

  // SDs for species/fleetgroup sel
  sigmaxSel50_sg    = exp(lnsigmaxSel50_sg);
  sigmaxSelStep_sg  = exp(lnsigmaxSelStep_sg);
  tauObs_spf        = exp(lntauObs_spf);
  tau2Obs_spf       = exp(2 * lntauObs_spf);

  // Loop over gears, species, stocks and time steps, create a matrix
  // of deviations in selectivity pars and calculate the pars
  // for each time step
  int epsSelVecIdx  = 0;
  int epslnqVecIdx  = 0;
  int qFVecIdx      = 0;
  for(int f = 0; f < nF; f++ )
  {
    xSel50_sf.col(f) = xSel50_sg.col(group_f(f));
    xSelStep_sf.col(f) = xSelStep_sg.col(group_f(f));
    xSel95_sf.col(f) = xSel95_sg.col(group_f(f));
    for( int p = 0; p < nP; p++ )
    {
      vector<Type> tmpSelAlpha    = thetaSelAlpha_sg.col(group_f(f)) + sigmaxSel50_sg(group_f(f)) * epsxSel50_spf.col(f).col(p);
      vector<Type> tmpSelBeta     = thetaSelBeta_sg.col(group_f(f)) + sigmaxSel50_sg(group_f(f)) * epsxSelStep_spf.col(f).col(p);
      xSel50_spf.col(f).col(p)    = LB_xSelAlpha_sg.col(group_f(f)) + (UB_xSelAlpha_sg.col(group_f(f)) - LB_xSelAlpha_sg.col(group_f(f))) *(1/20 + 18 * exp(tmpSelAlpha) / 20 / (1 + exp(tmpSelAlpha)) );
      xSelStep_spf.col(f).col(p)  = log(19) * (UB_xSelAlpha_sg.col(group_f(f)) - LB_xSelAlpha_sg.col(group_f(f))) *(0.005 + 0.5 * exp(tmpSelBeta) /  (1 + exp(tmpSelBeta)) );
      xSel95_spf.col(f).col(p)    = xSel50_spf.col(f).col(p) + xSelStep_spf.col(f).col(p);
      for( int s = 0; s < nS; s++)
      {
        // Exponentiate q and tau for the index observations
        if( commSurv_f(f) == 0)
        {
          int fleetIdx = f - nComm;
          q_spf(s,p,f) = exp(lnqSurv_sf(s,fleetIdx) + tauqSurv_sf(s,fleetIdx) * deltaq_spf(s,p,f) );
        }

        if( commSurv_f(f) == 1 )
        {
          q_spf(s,p,f) = exp(lnqComm_spf(s,p,f));
        } 

        xSel50_spft(s,p,f,0)      = xSel50_spf(s,p,f);
        xSelStep_spft(s,p,f,0)    = xSelStep_spf(s,p,f);
        // initialise tv q
        q_spft(s,p,f,0)           = q_spf(s,p,f);

        for( int t = 0; t < nT; t++ )
        { 
          // Calculate fishin mortality from effort and
          // catchability scalars
          if( C_spft(s,p,f,t) > 0 )
          {
            if( calcFmethod == "effortModel" )
            {
              qF_spft(s,p,f,t) = exp(lnqF_spf(s,p,f) + tauqFdev * deltalnqFspft_vec(qFVecIdx));
              F_spft(s,p,f,t) += qF_spft(s,p,f,t) * E_pft(p,f,t);
              qFVecIdx++;
            }
          }
          // Bring previous time-step selecivity pars
          // forward
          if( t > 0 )
          {
            xSel50_spft(s,p,f,t)    = xSel50_spft(s,p,f,t-1);
            xSelStep_spft(s,p,f,t)  = xSelStep_spft(s,p,f,t-1);
            q_spft(s,p,f,t)         = q_spft(s,p,f,t-1);  
          }
          // Update the epsilon array          
          // Check if ages or lengths are observed this time step
          if( tvSelFleets(f) == 1)
          {
            for( int x = 0; x < nX; x++)
            {
              if( ((age_axspft(0,x,s,p,f,t) >= 0) & (ageLikeWt_g(group_f(f)) > 0)) |
                  ((len_lxspft(0,x,s,p,f,t) >= 0) & (lenLikeWt_g(group_f(f)) > 0)) )
              {
                epsxSel50_spft(s,p,f,t)   += sigmaSel * epsSelAlpha_vec(epsSelVecIdx);
                epsxSelStep_spft(s,p,f,t) += sigmaSel * epsSelBeta_vec(epsSelVecIdx);
                epsSelVecIdx++;
                break;
              }
            }

            // do the same for time varying q
            if( tvqFleets(f) == 1 )
              if( I_spft(s,p,f,t) > 0 & t > minTimeIdx_spf(s,p,f) )
              {
                epslnq_spft(s,p,f,t) += epslnq_vec(epslnqVecIdx);
                epslnqVecIdx++;
              }

            // Now update xSel50 and xSelStep - can happen
            // at first time step, since we're using the gear specific
            // selectivity curve when we don't have age observations
            // for a stock within a fleet/species
            xSel50_spft(s,p,f,t)    *= exp(epsxSel50_spft(s,p,f,t));
            xSelStep_spft(s,p,f,t)  *= exp(epsxSelStep_spft(s,p,f,t));

            // Now update q
            q_spft(s,p,f,t)         *= exp(epslnq_spft(s,p,f,t));
          } // END tvSelFleet condition
   
        } // END t loop
      } // END s loop
    } // END p loop
  } // END f loop

  // Add learning rate for catchability
  q_spft *= propQ_spft;

  // Compute xSel95 for the time-varying arrays
  xSel95_spft = xSel50_spft + xSelStep_spft;


  // recruitment deviations and initRecDevs - convert from vector
  // to array of deviations
  // Also want to include a matrix of recruitment deviations
  // for estimating correlation
  matrix<Type> omegaRmat_spt(nS*nP,nT);
  omegaRmat_spt.setZero();
  array<Type> omegaRbar_sp(nS,nP);
  omegaRbar_sp.setZero();

  array<Type> omegaRinit_asp(nA,nS,nP);
  omegaRinit_asp.setZero();
  array<Type> omegaR_spt(nS,nP,nT);
  omegaR_spt.setZero();
  int initVecIdx = 0;
  int devVecIdx = 0;
  for( int s = 0; s < nS; s++ )
    for( int p = 0; p < nP; p++ )
    {
      for( int t = tFirstRecDev_s(s); t <= tLastRecDev_s(s); t++ )
      {
        if(boundRecDevs == 0 )
          omegaR_spt(s,p,t) = omegaR_vec(devVecIdx);

        if( boundRecDevs == 1 )
          omegaR_spt(s,p,t) = -5. + 10. * invlogit(omegaR_vec(devVecIdx));
          
        // And fill the matrix of deviations
        omegaRmat_spt( s*nP + p, t ) = omegaR_spt(s,p,t);

        omegaRbar_sp(s,p) += omegaR_spt(s,p,t) / (tLastRecDev_s(s) - tFirstRecDev_s(s) + 1);

        devVecIdx++;
      }
      // And calculate initialisation devs if initialised
      // in a fished state.
      if(swRinit_s(s) == 1)
        for( int a = 0; a < A_s(s); a++ )
        {
          if( boundRecDevs == 0 )
            omegaRinit_asp(a,s,p) = omegaRinit_vec(initVecIdx) ;
          if( boundRecDevs == 1 )
            omegaRinit_asp(a,s,p) = -5. + 10. * invlogit( omegaRinit_vec(initVecIdx) );

          initVecIdx++;
        }
    }

  
  // Probably concentrate these out later....
  vector<Type>  tauD_f        = exp(lntauD_f);    
  
  // Set up model state arrays
  array<Type> B_axspt(nA,nX,nS,nP,nT);      // Biomass at age-pop-time
  array<Type> N_axspt(nA,nX,nS,nP,nT);      // Numbers at age-pop-time
  array<Type> R_spt(nS,nP,nT);              // Recruits by pop-time
  array<Type> bhR_spt(nS,nP,nT);            // Expected BH Recruits by pop-time
  array<Type> Z_axspt(nA,nX,nS,nP,nT);      // Total mortality by age, population and time
  array<Type> paZ_axspt(nA,nX,nS,nP,nT);    // Approximate total mortality by age, population and time
  array<Type> C_axspft(nA,nX,nS,nP,nF,nT);  // Predicted catch-at-age (in numbers), fleet, population and time
  array<Type> paC_axspft(nA,nX,nS,nP,nF,nT);// Predicted catch-at-age (in numbers), fleet, population and time by Pope's Approx
  array<Type> Cw_axspft(nA,nX,nS,nP,nF,nT); // Predicted catch-at-age (in weight), population and time
  array<Type> predCw_spft(nS,nP,nF,nT);     // Predicted catch (in weight) by population and time
  array<Type> predC_spft(nS,nP,nF,nT);      // Predicted catch (in numbers) by population and time
  array<Type> B_spt(nS,nP,nT);              // Total biomass by species, pop, time
  array<Type> SB_spt(nS,nP,nT);             // Spawning biomass by species, pop, time
  array<Type> vB_spft(nS,nP,nF,nT);         // Vulnerable biomass by species, pop, time
  array<Type> Nv_spft(nS,nP,nF,nT);         // Vulnerable numbers by species, pop, time
  array<Type> Nv_axspft(nA,nX,nS,nP,nF,nT); // Vulnerable numbers at age by species, pop, time
  array<Type> Surv_axsp(nA,nX,nS,nP);       // Eqbm survival at age
  array<Type> SSBpr_asp(nA,nS,nP);          // Spawning stock biomass per rec. at age



  // ---------------- Growth, Maturity, Stock-Recruitment ----------- //
  Wlen_ls.fill(-1.0);
  matLen_ls.fill(-1.0);
  for( int s = 0; s < nS; s++ )
  {
    // weight-at-length for each length bin midpoint
    Wlen_ls.col(s) = LWa_s(s) * pow( lenBinMids_l, LWb_s(s) );
    
    // proportion mature-at-length
    if(matX == "length")
      matLen_ls.col(s) = 1/(1 + exp( -1. * log(Type(19.0)) * ( lenBinMids_l - xMat50_s(s)) / (xMat95_s(s) - xMat50_s(s)) ) );  
  }

  // Now produce the probability of length-at-age matrix for each stock
  // (taken from LIME by Rudd and Thorson, 2017)
  // In the same loop, we're going to compute the ssbpr parameter phi
  // for each stock
  lenAge_axsp.setZero();
  probLenAge_laxsp.setZero();
  meanWtAge_axsp.setZero();
  phi_sp.setZero();
  Surv_axsp.setZero();
  SSBpr_asp.setZero();
  matAge_asp.setZero();

  
  for( int s = 0; s < nS; s++)
  {
    int A1 = A1_s(s);
    int A2 = A2_s(s);
    for( int p = 0; p < nP; p++ )
    {
      for( int x = 0; x < nX; x++)
      {

        lenAge_axsp.col(p).col(s).col(x) += L2_xsp(x,s,p) - L1_xsp(x,s,p);
        lenAge_axsp.col(p).col(s).col(x) *= (exp(-vonK_xsp(x,s,p) * A1 ) - exp(-vonK_xsp(x,s,p) * age) );
        lenAge_axsp.col(p).col(s).col(x) /= (exp(-vonK_xsp(x,s,p) * A1 ) - exp(-vonK_xsp(x,s,p) * A2) );
        lenAge_axsp.col(p).col(s).col(x) += L1_xsp(x,s,p);

        for( int a = 0; a < A_s(s); a ++ )
        {
          Type  sumProbs  = 0;
          Type sigmaL = sigmaLa_s(s) + sigmaLb_s(s) * lenAge_axsp(a,x,s,p);
          for( int l = 0; l < L_s(s); l++ )
          {
            // Get length bin ranges
            Type lenHi = lenBinMids_l(l) + .5*lenBinWidth;
            Type lenLo = lenBinMids_l(l) - .5*lenBinWidth;
            // Compute density in each bin
            // Compute LH tail of the distribution
            
            if( (l >= 0) & (l < L_s(s) - 1 ) )
            {
              probLenAge_laxsp(l,a,x,s,p) = pnorm( log(lenHi), log(lenAge_axsp(a,x,s,p)), sigmaL );
              // Subtract LH tail for interior bins
              if(l > 0)
                probLenAge_laxsp(l,a,x,s,p) -= pnorm( log(lenLo), log(lenAge_axsp(a,x,s,p)), sigmaL );
              // Accumulate probability
              sumProbs += probLenAge_laxsp(l,a,x,s,p);
            }
            
            // Now the RH tail
            if( l == L_s(s) - 1 ) 
              probLenAge_laxsp(l,a,x,s,p) = Type(1.0) - sumProbs;

            // Calculate mean weight at age from probLenAge
            meanWtAge_axsp(a,x,s,p) += probLenAge_laxsp(l,a,x,s,p) * Wlen_ls(l,s);
            // Calculate maturity at age from probLenAge if maturity
            // is length based
            if( matX == "length")
              if( (matLen_ls(l,s) > 0) )
                matAge_asp(a,s,p) += probLenAge_laxsp(l,a,x,s,p) * matLen_ls(l,s);


          } // END l loop

          // Calculate maturity at age if age based
          if( matX == "age" )
            matAge_asp(a,s,p) = 1 / (1 + exp( -1. * log(Type(19.0)) * ( age(a) - xMat50_s(s)) / (xMat95_s(s) - xMat50_s(s)) ) );

          // To compute ssbpr, we need to reduce by the fraction of spawners
          Surv_axsp(a,x,s,p) = exp(-a * M_xsp(x,s,p));        
          if( a == A_s(s) - 1 ) 
            Surv_axsp(a,x,s,p) /= (1. - exp( -M_xsp(x,s,p)) );
      

          
        } // END a loop
      } // END x loop
      // Compute ssbpr
      SSBpr_asp.col(p).col(s) = Surv_axsp.col(p).col(s).col(nX-1) * matAge_asp.col(p).col(s) * meanWtAge_axsp.col(p).col(s).col(nX-1);
    } // END p loop
  } // END s loop

  


  // --------- Selectivity -------- //
  array<Type> sel_lspft(nL,nS,nP,nF,nT);
  array<Type> sel_axspft(nA,nX,nS,nP,nF,nT);
  array<Type> maxSel_spf(nS,nP,nF);
  sel_lspft.setZero();
  sel_axspft.setZero();
  maxSel_spf.setZero();
  for( int t = 0; t < nT; t++ )
    for( int s = 0; s < nS; s++ )
      for( int p = 0; p < nP; p++ )
      {
        // Loop over fleets, create selectivities
        for( int f = 0; f < nF; f++ )
        {
          // selectivity-at-length
          if( selX == "length" )
          {
            sel_lspft.col(t).col(f).col(p).col(s).fill(1.);
            sel_lspft.col(t).col(f).col(p).col(s) /= ( 1 + exp( - log(Type(19.)) * ( lenBinMids_l - xSel50_spft(s,p,f,t) ) / xSelStep_spft(s,p,f,t)  ) );
            // // reduce selectivity below min L to 0
            sel_lspft.col(t).col(f).col(p).col(s).segment(0,minL_s(s)-1).fill(0);
            // convert selectivity-at-length to selectivity-at-age using probability matrix
            for( int x = 0; x < nX; x++)
            {
              for( int a = 0; a < A_s(s); a++ )
              {
                vector<Type> probLenAgea_l(nL) ;
                probLenAgea_l.setZero();
                probLenAgea_l = probLenAge_laxsp.col(p).col(s).col(x).col(a);
                sel_axspft(a,x,s,p,f,t) = (probLenAgea_l*sel_lspft.col(t).col(f).col(p).col(s)).sum();
              }
              // Reduce selectivity below minA to 0
              sel_axspft.col(t).col(f).col(p).col(s).col(x).segment(0,minA_s(s)-1).fill(0);
            } 
          }
          if( selX == "age" )
          {
            for( int x = 0; x < nX; x++)
            {
              sel_axspft.col(t).col(f).col(p).col(s).col(x) = 1/( 1 + exp( - log(Type(19.)) * ( age - xSel50_spft(s,p,f,t) ) / xSelStep_spft(s,p,f,t)  ) );
              // Reduce selectivity below minA to 0
              sel_axspft.col(t).col(f).col(p).col(s).col(x).segment(0,minA_s(s)-1).fill(0);
            }
          }

          // Rescale selectivity at age by the highest value in the plusgroup across modeled sexes
          if(nX == 2)
          {
            // Might be worth adding a the CppAD function CondExpGt here, to avoid
            // non-differentiability
            if( sel_axspft(A_s(s)-1,0,s,p,f,t) > sel_axspft(A_s(s)-1,1,s,p,f,t) )
            {
              sel_axspft.col(t).col(f).col(p).col(s).col(0) /= sel_axspft(A_s(s)-1,0,s,p,f,t);
              sel_axspft.col(t).col(f).col(p).col(s).col(1) /= sel_axspft(A_s(s)-1,0,s,p,f,t);

              maxSel_spf(s,p,f) = sel_axspft(A_s(s)-1,0,s,p,f,t);
            }
            if( sel_axspft(A_s(s)-1,0,s,p,f,t) <= sel_axspft(A_s(s)-1,1,s,p,f,t) )
            {
              sel_axspft.col(t).col(f).col(p).col(s).col(0) /= sel_axspft(A_s(s)-1,1,s,p,f,t);
              sel_axspft.col(t).col(f).col(p).col(s).col(1) /= sel_axspft(A_s(s)-1,1,s,p,f,t);
              maxSel_spf(s,p,f) = sel_axspft(A_s(s)-1,1,s,p,f,t);
            }
          }
          if( nX == 1 )
          {
            sel_axspft.col(t).col(f).col(p).col(s).col(0) /= sel_axspft(A_s(s)-1,0,s,p,f,t);
            maxSel_spf(s,p,f) = sel_axspft(A_s(s)-1,0,s,p,f,t);
          }
        } // END f loop
      } // END p loop

  // --------- Population Dynamics ----------- //
  // First, compute eqbm recruitment
  // Derived parameters
  for( int p = 0; p < nP; p++ )
    for( int s = 0; s < nS; s++)
      phi_sp(s,p) = SSBpr_asp.col(p).col(s).sum();

  // Calculate B0 from Rbar
  if( recOption == "avgR")
    B0_sp = Rbar_sp * phi_sp;

  R0_sp   = B0_sp / phi_sp;


  reca_sp = 4.*h_sp*R0_sp / ( B0_sp*(1.-h_sp) );
  recb_sp = (5.*h_sp-1.) / ( B0_sp*(1.-h_sp) );

  // Set all state arrays to zero
  B_axspt.setZero();
  R_spt.setZero();
  bhR_spt.setZero();
  Z_axspt.setZero();
  paZ_axspt.setZero();
  C_axspft.setZero();
  paC_axspft.setZero();
  Cw_axspft.setZero();
  B_spt.setZero();
  SB_spt.setZero();
  vB_spft.setZero();
  Nv_spft.setZero();
  Nv_axspft.setZero();
  N_axspt.setZero();
  predCw_spft.setZero();
  predC_spft.setZero();


  // Copies of state arrays for testing Baranov solver
  array<Type> B_xspt(nX,nS,nP,nT);
  array<Type> vB_axspft(nA,nX,nS,nP,nF,nT);
  array<Type> vB_xspft(nX,nS,nP,nF,nT);
  array<Type> U_spft(nS,nP,nF,nT);
  array<Type> endN_axspt(nA,nX,nS,nP,nT);

  endN_axspt.setZero();
  B_xspt.setZero();
  vB_axspft.setZero();
  vB_xspft.setZero();
  vB_spft.setZero();

  array<Type> nPosCatch_spf(nS,nP,nF);
  nPosCatch_spf.fill(1e-9);

  // Age/Length observations and CPUE
  array<Type> I_spft_hat(nS,nP,nF,nT);
  array<Type> aDist_axspft_hat(nA,nX,nS,nP,nF,nT);
  array<Type> lDist_lxspft_hat(nL,nX,nS,nP,nF,nT);
  array<Type> propF_lspft_hat(nL,nS,nP,nF,nT);
  array<Type> propF_lspft(nL,nS,nP,nF,nT);
  array<Type> CPUEnll_spf(nS,nP,nF);          // CPUE
  array<Type> residCPUE_spft(nS,nP,nF,nT);    // CPUE std resids
  residCPUE_spft.setZero();
  CPUEnll_spf.setZero();
  I_spft_hat.fill(-1);
  aDist_axspft_hat.fill(0);
  lDist_lxspft_hat.fill(0);
  propF_lspft_hat.fill(-1);
  propF_lspft.fill(-1);
  
  array<Type> simI_spft(nS,nP,nF,nT);
  simI_spft.fill(-1);

  // vonB RAL growth model likelihood //
  // Currently has lengths binned the same as the integration above, 
  // and uses a random-at-length likelihood to account for differences 
  // in fleet selectivity, and at different times
  array<Type> vonBnll_spf(nS,nP,nF);
  array<Type> ageAtLenResids_laxspft(nL,nA,nX,nS,nP,nF,nT);
  array<Type> nObsAgeAtLen_spf(nS,nP,nF);
  array<Type> nResidsAgeAtLen_spf(nS,nP,nF);
  array<Type> etaSumSqAgeAtLen_spf(nS,nP,nF);
  array<Type> tau2AgeAtLenObs_spf(nS,nP,nF);

  // Zero-init all arrays
  tau2AgeAtLenObs_spf.setZero();
  nObsAgeAtLen_spf.setZero();
  nResidsAgeAtLen_spf.setZero();
  etaSumSqAgeAtLen_spf.setZero();
  vonBnll_spf.setZero();
  ageAtLenResids_laxspft.setZero();

  // Loop over time, species
  for( int t = 0; t < nT; t++ )
  { 
    for( int s = 0; s < nS; s++ )
    {
      // Loop over stocks
      for( int p = 0; p < nP; p++ )
      {
        // Set plus group age
        int A    = A_s(s);
        int L    = L_s(s);
        for( int x = 0; x < nX; x++)
        {
          // Initial time-step
          if( t == 0 )
          {
            // Initialise population either at unfished
            // or using the estimated fished initialisation
            // Populate first year - split into nX groups
            if( recOption == "BH")
              N_axspt.col(t).col(p).col(s).col(x) += R0_sp(s,p) * Surv_axsp.col(p).col(s).col(x);

            if( recOption == "avgR")
              N_axspt.col(t).col(p).col(s).col(x) += Rbar_sp(s,p) * Surv_axsp.col(p).col(s).col(x);

            // // Non-eqbm initialisation
            if(swRinit_s(s) == 1)
              N_axspt.col(t).col(p).col(s).col(x) *= exp( sigmaR_sp(s,p) * omegaRinit_asp.col(p).col(s) );

          }


          // time series history
          if( t > 0 )
          {
            // Generate recruitment
            Type SBt = SB_spt(s,p,t-1);
            if( recOption == "avgR" )
              N_axspt(0,x,s,p,t) = Rbar_sp(s,p);
            if( recOption == "BH" )
              N_axspt(0,x,s,p,t) = reca_sp(s,p) * SBt / (1 + recb_sp(s,p) * SBt);
            
            N_axspt(0,x,s,p,t) *= exp( sigmaR_sp(s,p) * omegaR_spt(s,p,t) );

            // Now loop over ages and apply fishing mortality (no discarding yet)
            for( int a = 1; a < A; a ++ )
            {
              // Update numbers at age using Pope's Approximation
              if( calcFmethod == "calcPopesApprox" )
              {
                if( a > 0 )
                  N_axspt(a,x,s,p,t) = endN_axspt( a-1, x, s, p, t-1);
                if( a == A - 1 )
                  N_axspt(a,x,s,p,t) += endN_axspt( a, x, s, p, t-1);
              }
              // Update numbers at age using mortality
              if( calcFmethod == "solveBaranov" | calcFmethod == "effortModel")
              {
                if( a > 0 )     
                  N_axspt(a,x,s,p,t) = N_axspt( a-1, x, s, p, t-1) *  exp( - Z_axspt(a-1,x,s,p,t-1) );
                if( a == A - 1) 
                  N_axspt(a,x,s,p,t) += N_axspt(A-1, x, s, p, t-1) *  exp( - Z_axspt(A-1,x,s,p,t-1) );
              }

              
            }
          }


          // Save recruits in R_pt
          R_spt(s,p,t) += N_axspt(0,x,s,p,t) / 2;
          

          // Compute biomass at age and total biomass
          B_axspt.col(t).col(p).col(s).col(x) += N_axspt.col(t).col(p).col(s).col(x) * meanWtAge_axsp.col(p).col(s).col(x);

          B_xspt(x,s,p,t) += B_axspt.col(t).col(p).col(s).col(x).sum();


          // Loop over fleets and compute catch
          for( int f =0; f < nF; f++ )
          {
            // Calculate vulnerable numbers and biomass for this fleet
            Nv_axspft.col(t).col(f).col(p).col(s).col(x) = N_axspt.col(t).col(p).col(s).col(x) * sel_axspft.col(t).col(f).col(p).col(s).col(x);
            Nv_spft(s,p,f,t) += Nv_axspft.col(t).col(f).col(p).col(s).col(x).sum();

            vB_axspft.col(t).col(f).col(p).col(s).col(x) += B_axspt.col(t).col(p).col(s).col(x) * sel_axspft.col(t).col(f).col(p).col(s).col(x);
            vB_xspft(x,s,p,f,t) = vB_axspft.col(t).col(f).col(p).col(s).col(x).sum();
            vB_spft(s,p,f,t) += vB_xspft(x,s,p,f,t);
            

            if( x == nX - 1 )
            {
              // Calculate expected biomass indices
              if( calcIndex_spf(s,p,f) == 1)
              {
                I_spft_hat(s,p,f,t) = q_spft(s,p,f,t) * vB_spft(s,p,f,t);  
                
                // Simulate some index data
                SIMULATE{
                  simI_spft(s,p,f,t) = I_spft(s,p,f,t) * exp(rnorm(Type(0),tauObs_spf(s,p,f)));
                }

                // Now add contribution to CPUE likelihood
                if( I_spft(s,p,f,t) > 0)
                {
                  residCPUE_spft(s,p,f,t) = ( log(I_spft(s,p,f,t)) - log(I_spft_hat(s,p,f,t)) )/tauObs_spf(s,p,f);

                  CPUEnll_spf(s,p,f) -= idxLikeWt_g(group_f(f)) * dnorm( residCPUE_spft(s,p,f,t), Type(0), Type(1), true);
                }

              }
            } // END x == nX - 1 for calculating expected indices

            for( int l = minL_s(s)-1; l < L; l++ )
            {
              Type sumProbs = 0;

              for(int a = minA_s(s)-1; a < A; a++ )
              {
                probAgeLen_laxspft(l,a,x,s,p,f,t) = probLenAge_laxsp(l,a,x,s,p) * N_axspt(a,x,s,p,t) * sel_lspft(l,s,p,f,t);
                sumProbs += probAgeLen_laxspft(l,a,x,s,p,f,t);
              }
              for( int a = minA_s(s)-1; a < A; a++)
              {
                probAgeLen_laxspft(l,a,x,s,p,f,t) /= sumProbs;  
              }

              // // Get ALFreq obs and pred if observations are +ve
              // if( growthLikeWt > 0)
              //   if( ALFreq_laxspft.rotate(-1).col(l).col(t).col(f).col(p).col(s).col(x).sum() > 0 )
              //   {
              //     vector<Type> pred = ALFreq_laxspft.rotate(-1).col(l).col(t).col(f).col(p).col(s).col(x).segment(minA_s(s)-1,A - minA_s(s) + 1);
              //     vector<Type> obs  = probAgeLen_laxspft.rotate(-1).col(l).col(t).col(f).col(p).col(s).col(x).segment(minA_s(s)-1,A - minA_s(s) + 1);

                  
              //     ageAtLenResids_laxspft.rotate(-1).col(l).col(t).col(f).col(p).col(s).col(x).segment(minA_s(s)-1,A - minA_s(s) + 1) = 
              //       calcLogistNormLikelihood( obs, 
              //                                 pred,
              //                                 minPropAAL,
              //                                 etaSumSqAgeAtLen_spf(s,p,f),
              //                                 nResidsAgeAtLen_spf(s,p,f) );

              //     nObsAgeAtLen_spf(s,p,f) += 1;


              //   }
            }


          } // END f loop 

          // Compute total and spawning biomass
          B_spt(s,p,t) += B_axspt.col(t).col(p).col(s).col(x).sum();
          // The females are the only important
          // mature individuals
          if( SB_spt(s,p,t) == 0)        
            SB_spt(s,p,t) = (B_axspt.col(t).col(p).col(s).col(nX - 1) * matAge_asp.col(p).col(s)).sum();


        } // END x loop
        // // Need another f loop, can reorder loops later maybe?
        // for( int f = 0; f <  nF; f ++ )
        // {
        //   if( nX == 2 & lDist_lxspft_hat.col(t).col(f).col(p).col(s).col(nX).sum() > 0 )
        //     lDist_lxspft_hat.col(t).col(f).col(p).col(s).col(nX) /= lDist_lxspft_hat.col(t).col(f).col(p).col(s).col(nX).sum();
        // }

      } // END p loop
    } // END s loop

    array<Type> tmpZ_axsp(nA,nX,nS,nP);
    array<Type> tmp_paZ_axsp(nA,nX,nS,nP);
    tmpZ_axsp = Z_axspt.col(t);
    tmp_paZ_axsp.setZero();

    array<Type> tmpF_spf(nS,nP,nF);
    array<Type> tmp_paF_spf(nS,nP,nF);
    tmpF_spf = F_spft.col(t);
    tmp_paF_spf = paF_spft.col(t);

    // Catch in biomass
    array<Type> tmpCw_axspf(nA,nX,nS,nP,nF);
    tmpCw_axspf.setZero();
    // Catch in numbers
    array<Type> tmpC_axspf(nA,nX,nS,nP,nF);
    tmpC_axspf.setZero();

    // if( calcFmethod == "calcPopesApprox" )
    // {
    //   endN_axspt.col(t) = calcPopesApprox(  N_axspt.col(t),  
    //                                         sel_axspft.col(t),
    //                                         C_spft.col(t),   
    //                                         M_xsp,   
    //                                         A_s,     
    //                                         meanWtAge_axsp, 
    //                                         vB_axspft.col(t),
    //                                         vB_spft.col(t),
    //                                         tmp_paF_spf,  
    //                                         tmp_paZ_axsp,
    //                                         tmpCw_axspf,
    //                                         tmpC_axspf );


    //   paF_spft.col(t) = tmp_paF_spf;
    //   paZ_axspt.col(t) = tmp_paZ_axsp;
    //   paC_axspft.col(t) = tmpC_axspf;

    //   // Save outputs from F calculation
    //   Cw_axspft.col(t)  = tmpCw_axspf;
    //   Z_axspt.col(t)    = tmpZ_axsp;
    //   F_spft.col(t)     = tmpF_spf;

    // }
    

    if( calcFmethod == "solveBaranov")
    {
      // Baranov solver returns predicted catch in biomass units
      tmpCw_axspf += solveBaranov_spf(  nBaranovIter,
                                        lambdaBaranovStep,
                                        A_s,
                                        C_spft.col(t),
                                        M_xsp,
                                        B_axspt.col(t),
                                        vB_axspft.col(t),
                                        vB_spft.col(t),
                                        sel_axspft.col(t),
                                        tmpZ_axsp,
                                        tmpF_spf,
                                        N_axspt.col(t) );      // Numbers at age
      // Save outputs from F calculation
      Cw_axspft.col(t)  = tmpCw_axspf;
      Z_axspt.col(t)    = tmpZ_axsp;
      F_spft.col(t)     = tmpF_spf;
    }

    
    

    

    for( int s = 0; s < nS; s++)
      for( int p = 0; p < nP; p++ )
        for( int f = 0; f < nF; f++ )
        {
          // Set plus group age
          int A    = A_s(s);
          int L    = L_s(s);
          if( C_spft(s,p,f,t) > 0 )
          {
            if( calcFmethod == "effortModel" )
            {
              // Compute total mortality by age/sex for Baranov eqn
              for( int x = 0; x < nX; x++ )
                for( int a = 0; a < A_s(s); a++)
                {
                  Z_axspt(a,x,s,p,t) += M_xsp(x,s,p) + sel_axspft(a,x,s,p,f,t) * F_spft(s,p,f,t);
                }
            }

            for( int x = 0; x < nX; x++ )
            {
              for( int a = 0; a < A_s(s); a++)
              {
                if( calcFmethod == "effortModel" )
                {
                  // Now apply Baranov equation
                  Type tmpF_axspft        = sel_axspft(a,x,s,p,f,t) * F_spft(s,p,f,t);
                  Cw_axspft(a,x,s,p,f,t)  = vB_axspft(a,x,s,p,f,t) * (1 - exp(-Z_axspt(a,x,s,p,t))) * tmpF_axspft / Z_axspt(a,x,s,p,t);
                }


                // Calc catch-at-age
                C_axspft(a,x,s,p,f,t)  = Cw_axspft(a,x,s,p,f,t) / meanWtAge_axsp(a,x,s,p);  
              }

              // Calculate expected age and length comps
              // Create array of predicted age distributions - this should just be catch-at-age props
              // Calculate total catch
              Type totCatch = C_axspft.col(t).col(f).col(p).col(s).col(x).sum();
                
              // Convert catch-at-age to proportions-at-age
              if(totCatch > 0)
                aDist_axspft_hat.col(t).col(f).col(p).col(s).col(x) = C_axspft.col(t).col(f).col(p).col(s).col(x) / totCatch;
                  


              // Create array of predicted length distributions
              // Need to do some thinking here... follow LIME to start with...
              // Probability of being harvested at age
              vector<Type> probHarvAge(nA);
              probHarvAge.setZero();
              probHarvAge = Nv_axspft.col(t).col(f).col(p).col(s).col(x) / N_axspt.col(t).col(p).col(s).col(x).sum();
              
              // Probability of sampling a given length bin
              vector<Type> probHarvLen(L);
              probHarvLen.setZero();
              // loop over length and age, calculate probability of harvest at length
              for( int l = 0; l < L; l++ )
              {
                if( lenComps == "LIME" )
                {
                  for( int a = 0; a < A; a++ )
                    probHarvLen(l) += probHarvAge(a)*probLenAge_laxsp(l,a,x,s,p);
                }

                if( lenComps == "Francis" )
                {
                  for( int a = 0; a < A; a++ )
                    probHarvLen(l) += probLenAge_laxsp(l,a,x,s,p) * N_axspt(a,x,s,p,t) * sel_lspft(l,s,p,f,t); 
                }

                // Save to length dists
                // lDist_lxspft_hat(l,x,s,p,f,t) = probHarvLen(l);
                // // Add total length dist so we can fit
                // // to unsexed observations
                // if( nX == 2)
                // {
                //   // Add observations for unsexed comps
                //   lDist_lxspft_hat(l,nX,s,p,f,t) += lDist_lxspft_hat(l,x,s,p,f,t);

                //   // Calculate proportion female
                //   if( x == 1 )
                //     propF_lspft_hat(l,s,p,f,t) = lDist_lxspft_hat(l,x,s,p,f,t) /  lDist_lxspft_hat(l,nX,s,p,f,t);
                // }

              }
              // renormalise
              if(probHarvLen.sum() > 0)
                lDist_lxspft_hat.col(t).col(f).col(p).col(s).col(x).segment(0,L) = probHarvLen / probHarvLen.sum();


            }

            predCw_spft(s,p,f,t) += Cw_axspft.col(t).col(f).col(p).col(s).sum();
            predC_spft(s,p,f,t)  += C_axspft.col(t).col(f).col(p).col(s).sum();
          
            nPosCatch_spf(s,p,f) += 1;


          }
        }

    Fbar_spf += F_spft.col(t);
  } // END t loop for pop dynamics

  Fbar_spf /= nPosCatch_spf;
    
  // Compute exploitation rate
  U_spft = predCw_spft / (1e-9 + vB_spft);

  // --------- Statistical Model ----------- //

  // Post-fit recruitment model
  array<Type> deltaRec_spt(nS,nP,nT-1);
  array<Type> expRec_spt(nS,nP,nT-1);
  array<Type> SRnlp_sp(nS,nP);
  SRnlp_sp.setZero();

  if( (recOption == "avgR") & (postFitSR == 1) )
  {
    for( int s = 0; s < nS; s++ )
      for( int p = 0; p < nP; p++ )
        for( int t = 1; t < nT; t++ )
        {
          expRec_spt(s,p,t-1) = reca_sp(s,p) * SB_spt(s,p,t-1) / (1 + recb_sp(s,p) * SB_spt(s,p,t-1));
          deltaRec_spt(s,p,t-1) = log( R_spt(s,p,t) / expRec_spt(s,p,t-1) );  
          Type resid = deltaRec_spt(s,p,t-1) / sigmaR_sp(s,p);
          SRnlp_sp(s,p) -= dnorm( resid, Type(0.), Type(1), true);
        } // END t loop
  } // END postfit SR model


  // Process model prior //
  // Recruitment deviations (initial and ongoing) - add correlation later
  // both auto-correlation and between-pop correlation can be added.
  Type recnlp = 0;
  Type corrRecnlp = 0.;
  recnlp -= dnorm( omegaRinit_vec,Type(0), Type(1), true).sum();
  if( recruitVariance == "uncorr" )
  {
    recnlp -= dnorm( omegaR_vec, Type(0), Type(1), true).sum();
    for( int p = 0; p < nP; p++)
    {
      vector<Type> omegaRbar_s(nS);
      omegaRbar_s = omegaRbar_sp.col(p);
      recnlp -= dnorm( omegaRbar_s, Type(0), sd_omegaRbar, true).sum();
    }
    
  }

  if( recruitVariance == "corr" )
  {
    // Use the unstructured correlation matrix here
    UNSTRUCTURED_CORR_t<Type> recCorrnlp(recCorr_vec);

    for( int t = 0; t < nT; t++ )
      corrRecnlp += recCorrnlp( omegaRmat_spt.col(t) );

    recCorrMat_sp += recCorrnlp.cov();

    int nPops = nS * nP;

    matrix<Type> IWscale(nPops,nPops);
    IWscale = IWmode * (IWnu + nPops + 1);

    // Now apply an IW prior to the recruitment correlation matrix
    matrix<Type> traceMat = IWscale * recCorrMat_sp.inverse();
    Type trace = 0.0;
    for( int sp = 0; sp < nPops ; sp++ ) 
      trace += traceMat(sp,sp);

    corrRecnlp += Type(0.5) *( (IWnu + nPops + 1) * atomic::logdet(recCorrMat_sp) + trace);
    

    // And add an inverse wishart prior
    recnlp += corrRecnlp;
  }

  for( int s = 0; s < nS; s++ )
    for( int p = 0; p < nP; p++ )
    {
      recnlp -= dinvgamma(pow(sigmaR_sp(s,p),2),IGalphaR,square(pmsigmaR)*(IGalphaR + 1),1);
    }

  // Observation likelihoods
  array<Type> ageCompsnll_spf(nS,nP,nF);      // Age observations
  array<Type> lenCompsnll_spf(nS,nP,nF);      // length observations
  array<Type> Ctnll_spf(nS,nP,nF);            // Catch
  array<Type> Dtnll_sp(nS,nP);                // Discards
  // Conditional MLEs of age/length observation error variance
  array<Type> propFnll_spf(nS,nP,nF);         // Proportion female NLL
  array<Type> tau2PropF_spf(nS,nP,nF);        // Proportion female  conditional variance
  array<Type> nObsPropF_spf(nS,nP,nF);        // Proportion female number of observations
  array<Type> resPropF_spf(nS,nP,nF);         // Proportion female number of log-residue
  array<Type> tau2Age_spf(nS,nP,nF);
  array<Type> tau2Len_spf(nS,nP,nF);
  array<Type> etaSumSqLen_spf(nS,nP,nF);
  array<Type> etaSumSqAge_spf(nS,nP,nF);
  array<Type> nResidsAge_spf(nS,nP,nF);
  array<Type> nResidsLen_spf(nS,nP,nF);
  array<Type> nObsAge_spf(nS,nP,nF);
  array<Type> nObsLen_spf(nS,nP,nF);
  array<Type> intAgeLikeCpt_spf(nS,nP,nF);
  array<Type> intLenLikeCpt_spf(nS,nP,nF);
  array<Type> ageRes_axspft(nA,nX,nS,nP,nF,nT);
  array<Type> lenRes_lxspft(nL,nX,nS,nP,nF,nT);
  // Zero-init
  ageCompsnll_spf.setZero();
  lenCompsnll_spf.setZero();
  propFnll_spf.setZero();
  tau2PropF_spf.setZero();
  nObsPropF_spf.setZero();
  resPropF_spf.setZero();
  Ctnll_spf.setZero();
  Dtnll_sp.setZero();
  tau2Age_spf.setZero();
  tau2Len_spf.setZero();
  etaSumSqLen_spf.setZero();
  etaSumSqAge_spf.setZero();
  nResidsAge_spf.setZero();
  nResidsLen_spf.setZero();
  nObsAge_spf.setZero();
  nObsLen_spf.setZero();
  ageRes_axspft.setZero();
  lenRes_lxspft.setZero();


  // Loop over age table, calculate residuals and likelihood
  for( int aIdx = 0; aIdx < nAgeObs; aIdx ++  )
  {
    // Get array dims from table
    int t = age_table(aIdx,0)-1;
    int s = age_table(aIdx,1)-1;
    int p = age_table(aIdx,2)-1;
    int f = age_table(aIdx,3)-1;
    int x = age_table(aIdx,4)-1;

    // Get species ages
    int A      = A_s(s) - minA_s(s) + 1;    

    // tmp vectors to hold age and length observed and
    // predicted values
    vector<Type> fleetAgeObs(A);
    vector<int>  intFleetAgeObs(A);
    vector<Type> fleetAgePred(A);

    // Set tmp vectors to zero 
    fleetAgeObs.setZero();
    fleetAgePred.setZero();

    intFleetAgeObs   = age_table.transpose().col(aIdx).segment(5 + minA_s(s)-1,A);
    fleetAgeObs      = intFleetAgeObs.cast<Type>(); 
    fleetAgePred    += aDist_axspft_hat.col(t).col(f).col(p).col(s).col(x).segment(minA_s(s)-1,A);    

    // Now for compositional data. First, check if ages are 
    // being used (or exist), if so
    // compute logistic normal likelihood
    if( (fleetAgeObs.sum() > 0) & (fleetAgePred.sum() > 0) & (ageLikeWt_g(group_f(f)) > 0) )
    {
    
      matrix<Type> tmpCorr_aa = CorrAge_sfaa.transpose().col(s).col(f).matrix();

      ageRes_axspft.col(t).col(f).col(p).col(s).col(x).segment(minA_s(s) - 1,A) =
        calcCorrLogistNormLikelihood( fleetAgeObs,
                                      fleetAgePred,
                                      minAgeProp,
                                      etaSumSqAge_spf(s,p,f),
                                      nResidsAge_spf(s,p,f),
                                      meanAgeSampSize_spf(s,p,f),
                                      tmpCorr_aa,
                                      intAgeLikeCpt_spf(s,p,f) );                

      // Increment number of age observations
      nObsAge_spf(s,p,f) += 1;

    }     
  } // END aIdx loop
  // End age obs table likelihood loop

  // Loop over len table, calc residuals and likelihood
  for( int lIdx = 0; lIdx < nLenObs; lIdx ++  )
  {
    // Get array dims from table
    int t = len_table(lIdx,0)-1;
    int s = len_table(lIdx,1)-1;
    int p = len_table(lIdx,2)-1;
    int f = len_table(lIdx,3)-1;
    int x = len_table(lIdx,4)-1;

    // Get species ages
    int L      = L_s(s) - minL_s(s) + 1;

    vector<Type> fleetLenObs(L);
    vector<int>  intFleetLenObs(L);
    vector<Type> fleetLenPred(L);  

    // Set tmp vectors to zero 
    fleetLenObs.setZero();
    fleetLenPred.setZero();

    intFleetLenObs   = len_table.transpose().col(lIdx).segment(5 + minL_s(s)-1,L);
    fleetLenObs      = intFleetLenObs.cast<Type>();
    fleetLenPred    += lDist_lxspft_hat.col(t).col(f).col(p).col(s).col(x).segment(minL_s(s)-1,L);


    if( (fleetLenObs.sum() > 0) & (fleetLenPred.sum() > 0) &  (lenLikeWt_g(group_f(f)) > 0) )
    {
      matrix<Type> tmpCorr_ll = CorrLen_sfll.transpose().col(s).col(f).matrix();

      lenRes_lxspft.col(t).col(f).col(p).col(s).col(x).segment(minL_s(s)-1,L) = 
        calcCorrLogistNormLikelihood( fleetLenObs, 
                                      fleetLenPred,
                                      minLenProp,
                                      etaSumSqLen_spf(s,p,f),
                                      nResidsLen_spf(s,p,f),
                                      meanLenSampSize_spf(s,p,f),
                                      tmpCorr_ll,
                                      intAgeLikeCpt_spf(s,p,f) ); 

      nObsLen_spf(s,p,f) += 1;
    }


  } // END lIdx loop
  // END length obs likelihood table 


  // Final age and length comp calcs
  for( int s = 0; s < nS; s++ )
    for( int p = 0; p < nP; p++ )
    {
      int A      = A_s(s) - minA_s(s) + 1;
      int L      = L_s(s) - minL_s(s) + 1;
      // Loop over fleets
      for( int f = 0; f < nF; f++ )
      {   
        if( nResidsAge_spf(s,p,f) > 0)
        {
          tau2Age_spf(s,p,f)      += etaSumSqAge_spf(s,p,f) / nResidsAge_spf(s,p,f);
          ageCompsnll_spf(s,p,f)  += ageLikeWt_g(group_f(f)) *( 
                                        0.5 * nResidsAge_spf(s,p,f) * log(tau2Age_spf(s,p,f)) +
                                        intAgeLikeCpt_spf(s,p,f) +
                                        0.5 * etaSumSqAge_spf(s,p,f) / tau2Age_spf(s,p,f)) ;

        }
        if( nResidsLen_spf(s,p,f) > 0)
        {
          tau2Len_spf(s,p,f)      += etaSumSqLen_spf(s,p,f) / nResidsLen_spf(s,p,f);
          lenCompsnll_spf(s,p,f)  += lenLikeWt_g(group_f(f)) * 0.5 * (nResidsLen_spf(s,p,f) - nObsLen_spf(s,p,f))  * log(tau2Len_spf(s,p,f));
        }

      }
    } // END p loop
  // END s loop
  // End data likelihoods



  // vonB RAL growth model likelihood //
  // Currently has lengths binned the same as the integration above, 
  // and uses a random-at-length likelihood to account for differences 
  // in fleet selectivity, and at different times
  // array<Type> vonBnll_spf(nS,nP,nF);
  // array<Type> ageAtLenResids_laxspft(nL,nA,nX,nS,nP,nF,nT);
  // array<Type> nObsAgeAtLen_spf(nS,nP,nF);
  // array<Type> nResidsAgeAtLen_spf(nS,nP,nF);
  // array<Type> etaSumSqAgeAtLen_spf(nS,nP,nF);
  // array<Type> tau2AgeAtLenObs_spf(nS,nP,nF);


  // Zero-init all arrays
  // tau2AgeAtLenObs_spf.setZero();
  // nObsAgeAtLen_spf.setZero();
  // nResidsAgeAtLen_spf.setZero();
  // etaSumSqAgeAtLen_spf.setZero();
  // vonBnll_spf.setZero();
  // ageAtLenResids_laxspft.setZero();

  // Simulate the age-at-length data
  array<Type> simAgeLen_alspftx(nA,nL,nS,nP,nF,nT,nX);
  simAgeLen_alspftx.fill(0);

  // Now loop and compute
  if( growthLikeWt > 0 )
  {
    for( int rIdx = 0; rIdx < nAgeLenObs; rIdx++ )
    {
      // Pull indices
      int t = aal_table(rIdx,0)-1;
      int s = aal_table(rIdx,1)-1;
      int p = aal_table(rIdx,2)-1;
      int f = aal_table(rIdx,3)-1;
      int x = aal_table(rIdx,4)-1;
      int l = aal_table(rIdx,5)-1;

      // containers for a given year/fleet/length/species/pop - function
      // expects vector<Type>
      int minA =  minA_s(s);
      int maxA =  A_s(s);
      int A = maxA - minA + 1;
      vector<Type> obs(A);
      vector<int>  intObs(A);
      vector<Type> pred(A);
      vector<Type> resids(A);

      // Zero init containers
      intObs.setZero();
      obs.setZero();
      pred.setZero();
      resids.setZero();

      intObs      += aal_table.transpose().col(rIdx).segment(6 + minA - 1,A);
      obs         = intObs.cast<Type>();
      for( int a = 0; a < A; a++)
        pred(a)  = probAgeLen_laxspft(l,a + minA - 1,x,s,p,f,t);
      // pred      += probAgeLen_laxspft.rotate(-1).col(l).col(t).col(f).col(p).col(s).col(x).segment(minA-1,A);
      // Compute resids etc. if observations aren't missing
      resids = 
        calcLogistNormLikelihood( obs, 
                                  pred,
                                  minPropAAL,
                                  etaSumSqAgeAtLen_spf(s,p,f),
                                  nResidsAgeAtLen_spf(s,p,f) );

      for( int a = 0; a < A; a++ )
        ageAtLenResids_laxspft(l,a+minA-1,x,s,p,f,t) = resids(a);

      nObsAgeAtLen_spf(s,p,f) += 1;

    } // END rIdx loop
  
    // Calculatte RAL likelihood
    for( int s = 0; s < nS; s++)
      for( int p = 0; p < nP; p++ )
        for( int f = 0; f < nF; f++)
          if( nResidsAgeAtLen_spf(s,p,f) > 0)
          {
            tau2AgeAtLenObs_spf(s,p,f) += etaSumSqAgeAtLen_spf(s,p,f) / nResidsAgeAtLen_spf(s,p,f);
            vonBnll_spf(s,p,f) += 0.5 * (nResidsAgeAtLen_spf(s,p,f) - nObsAgeAtLen_spf(s,p,f));
            vonBnll_spf(s,p,f) *= log(tau2AgeAtLenObs_spf(s,p,f));
          }
  }
  // END RAL growth model likelihood 

  // Now a regularisation prior on F_spft and
  // an IG prior on observation error variance
  vector<Type> tauObsnlp_f(nF);
  tauObsnlp_f.setZero();
  Type Fnlp = 0;
  Type qFdev_nlp = 0;
  for( int s = 0; s < nS; s++ )
    for( int p = 0; p < nP; p++ )
    {
      for( int f = 0; f < nF; f++ )
      {
        tauObsnlp_f(f) -= dinvgamma( pow(tauObs_spf(s,p,f),2), IGalphaObs, square(pmtauObs_g(group_f(f))) * (IGalphaObs + 1), 1);
        if( regFfleets(f) > 0)
        {
          vector<Type> Fvec(nT);
          vector<Type> Uvec(nT);
          Fvec = F_spft.rotate(1).col(f).col(p).col(s);
          Uvec = U_spft.rotate(1).col(f).col(p).col(s);
          // Penalize deviations of F from exploitation rate.
          // Fnlp -= dnorm(Fvec, Uvec, sdF, true).sum();
          Fnlp -= dnorm(Fbar_spf(s,p,f), mF, sdF, true);
          // for( int t = 0; t < nT; t++)
          // {
          //   Type tmpF = F_spft(s,p,f,t);
          //   Fnlp += CppAD::CondExpGe( tmpF, Type(1), Type(1e1), 0 );
          // }
        }

        if( calcFmethod == "effortModel" )
        {
          for( int t = 0; t < nT; t++ )
          {
            if( C_spft(s,p,f,t) > 0 &  predCw_spft(s,p,f,t) > 0 )
            {
              Type catRes = log(C_spft(s,p,f,t)) - log(predCw_spft(s,p,f,t));
              Ctnll_spf(s,p,f) -= dnorm( catRes, Type(0), Type(0.02), true );
            }
          }

          qFdev_nlp -= dnorm( deltalnqFspft_vec, Type(0), 1.0, true ).sum();

        }
        
      }
        
    }

  
  // ------------------ Shared Hierarchical Prior Distributions ------------------ //
  // Steepness
  vector<Type>  steepnessnlp_s(nS);
  array<Type>   steepnessnlp_sp(nS,nP);
  Type steepnessnlp = 0.;

  // Mortality
  vector<Type>  Mnlp_s(nS);
  array<Type>   Mnlp_sp(nS,nP);
  array<Type>   Mnlp_sx(nS,nX);
  Type Mnlp = 0.;

    
  // 0-initialise
  steepnessnlp_s.setZero();
  steepnessnlp_sp.setZero();
  Mnlp_s.setZero();
  Mnlp_sp.setZero();
  Mnlp_sx.setZero();
  Type sel_nlp = 0.;
  // VonB deviations
  Type L1nlp = 0;
  Type L2nlp = 0;
  Type VonKnlp = 0;

  // vague priors on theta pars for VonB model
  // L1nlp   -= dnorm( thetaL1_s, Type(0), cvL1, true ).sum();
  // L2nlp   -= dnorm( thetaL2_s, Type(0), cvL2, true ).sum();
  // VonKnlp -= dnorm( thetaVonK_s, Type(0), cvVonK, true ).sum();

  // Biological parameters
  for( int s = 0; s < nS; s++ )
  {

    steepnessnlp_s(s)     -= dinvgamma( pow(sigmah_s(s),2), IGalphah, square(pmsigmah) * (IGalphah + 1), 1);
    Mnlp_s(s)             -= dinvgamma( pow(sigmaM_s(s),2), IGalphaM, square(pmsigmaM) * (IGalphaM + 1), 1);

    

    for( int p = 0; p < nP; p ++)
    {
      // Steepness
      steepnessnlp_sp(s,p)  -= dnorm( epsSteep_sp(s,p), Type(0), Type(1), true);
      // Mortality
      Mnlp_sp(s,p)          -= dnorm( epsM_sp(s,p), Type(0), Type(1), true);

      // vonB pars
      VonKnlp -= dnorm( deltaVonK_sp(s,p), Type(0), Type(1), true);
      L2nlp   -= dnorm( deltaL2_sp(s,p), Type(0), Type(1), true);
    }

    for( int x = 0; x < nX; x++ )
    {
      Mnlp_sx(s,x) -= dnorm(epsM_sx(s,x), Type(0), Type(1), true);

      VonKnlp -= dnorm( deltaVonK_xs(x,s), Type(0), Type(1), true);
      L2nlp   -= dnorm( deltaL2_xs(x,s), Type(0), Type(1), true);

    }

  }

  // add species level h prior
  steepnessnlp_s  -= dnorm( epsSteep_s, 0., Type(1), true );
  Mnlp_s          -= dnorm( epsM_s, 0., Type(1), true );
  
  // Now prior on complex mean M and steepness
  Mnlp            -= dnorm( lnM, log(muM), sdM, true);
  Mnlp            -= dinvgamma( pow(sigmaM,2), IGalphaM, square(pmsigmaM) * (IGalphaM + 1), true);
  steepnessnlp    -= dbeta( h, hBetaPrior(0), hBetaPrior(1), true);
  steepnessnlp    -= dinvgamma( pow(sigmah,2), IGalphah, square(pmsigmah) * (IGalphah + 1), 1);

  // catchability
  Type qnlp_tv = 0;       // time varying catchability
  Type qnlp_stock = 0;    // stock specific catchability devs
  Type qnlp_gps = 0;      // species specific group level catchability devs

  qnlp_tv     -= dnorm( epslnq_vec, Type(0), Type(1), true).sum();
  qnlp_stock  -= dnorm( deltaqspf_vec, Type(0), Type(1), true).sum();

  // Selectivity
  // Add time-varying selectivity deviations
  sel_nlp -= dnorm( epsSelAlpha_vec, Type(0), 1., true).sum();
  sel_nlp -= dnorm( epsSelBeta_vec, Type(0), 1., true).sum();
  // Stock specific sel deviations
  sel_nlp -= dnorm( epsSelAlphaspf_vec, Type(0), Type(1), true).sum();
  sel_nlp -= dnorm( epsSelBetaspf_vec, Type(0), Type(1), true).sum();  

  // Priors on fleetGroup leading pars
  for( int g = 0; g < nGroups; g++ )
  {
    vector<Type> value  = thetaSelAlpha_sg.col(g);
    sel_nlp -= dnorm( value, Type(0) ,cvSel, true).sum();

    value  = thetaSelAlpha_sg.col(g);
    sel_nlp -= dnorm( value, Type(0), cvSel, true).sum();

    // Add IG prior for tauq_g and tauq_sg
    for( int s = 0; s < nS; s++)
    {
      sel_nlp   -= dinvgamma( pow(sigmaxSel50_sg(s,g),2), IGalphaSel, square(pmsigmaSel_g(g)) * (IGalphaSel + 1),1);
      sel_nlp   -= dinvgamma( pow(sigmaxSelStep_sg(s,g),2), IGalphaSel, square(pmsigmaSel_g(g)) * (IGalphaSel + 1),1);  
    }
  }

  for( int f = 0; f < nF; f++ )
  {

    // Catchability
    // vector<Type> qDevVec = deltaqSurv_sf.col(f);
    // qnlp_gps -= dnorm( qDevVec, Type(0), Type(1), true).sum();
    // qnlp_gps -= dinvgamma( pow(tauq_f(f),2), IGalphaq, square(pmtauqSurv_f(f)) * (IGalphaq + 1), 1);
    for( int s = 0; s < nS; s++)
    {
      // Normal prior on catchability
      if( commSurv_f(f) == 0)
      {
        int fleetIdx = f - nComm;
        qnlp_gps  -= dnorm( qSurv_sf(s,fleetIdx), mq_f(fleetIdx), sdq_f(fleetIdx), true);
        qnlp_gps  -= dinvgamma( pow(tauqSurv_sf(s,fleetIdx),2), IGalphaq, square(pmtauqSurv_f(fleetIdx)) * (IGalphaq + 1), 1);
      }

      if( solveQ_f(f) == 1)
        for( int p = 0; p < nP; p++ )
        {
          vector<Type> q_t(nT); 
          q_t = q_spft.transpose().col(s).col(p).col(f);
          qnlp_gps -= dnorm(log(q_t), log(q_spf(s,p,f)), tauqSurv_sf(s,f), true).sum();
        }
    }
  }
  
  // Observations
  f += Ctnll_spf.sum();       // Catch
  f += qFdev_nlp;
  f += lnqF_spf.sum();        // Jeffries prior on qF
  f += lenCompsnll_spf.sum(); // Length compositions
  f += ageCompsnll_spf.sum(); // Age compositions
  f += CPUEnll_spf.sum();     // Survey CPUE
  f += sel_nlp;
  f += Fnlp;
  f += tauObsnlp_f.sum();
  f += lambdaPropF * propFnll_spf.sum();
  // Recruitment errors
  f += recnlp;      // recruitment process errors
  // q groups
  f += qnlp_stock + qnlp_gps + qnlp_tv;
  // Biological parameters
  f += steepnessnlp_s.sum() + steepnessnlp_sp.sum() + steepnessnlp;
  f += Mnlp_s.sum() + Mnlp_sp.sum() + Mnlp_sx.sum() + Mnlp;
  f += lambdaB0 * lnB0_sp.sum();
  f += SRnlp_sp.sum();

  // Growth model
  f += growthLikeWt * vonBnll_spf.sum();
  f += growthLikeWt * VonKnlp;
  f += growthLikeWt * L2nlp;
  f += growthLikeWt * L1nlp;


  if( compLikeFun == 1 )
  {
    f += (logitphi1Age_sf*logitphi1Age_sf).sum();
    f += (logitphi1Len_sf*logitphi1Len_sf).sum();
  }

  if( compLikeFun >= 2 )
  {
    f += (logitpsiAge_sf*logitpsiAge_sf).sum();
    f += (logitpsiLen_sf*logitpsiLen_sf).sum();
  }

  joint_nlp += f;

  
  // // Return quantities
  // Model dimensions
  REPORT( nS );
  REPORT( nP );
  REPORT( nF );
  REPORT( nT );
  REPORT( nL );
  REPORT( nA );
  REPORT( nX );

  // Natural scale leading parameters
  REPORT( B0_sp ); 
  REPORT( Rbar_sp );   
  REPORT( h_sp );     
  REPORT( M_xsp );     
  REPORT( L2_xsp );  
  REPORT( vonK_xsp );  
  REPORT( L1_xsp ); 
  REPORT( sigmaLa_s );
  REPORT( sigmaLb_s );
  REPORT( sigmaR_sp );

  // Fishery/Survey model pars
  REPORT( q_spf );          // Catchability
  REPORT( qSurv_sf );       // Species/Fleet catchability for surveys
  REPORT( qComm_spf );      // Species/Stock/Fleet catchability for Comm fleets
  REPORT( q_spft );         // Time-varying catchability
  REPORT( xSel50_spft );    // length at 50% sel by species/fleet/pop
  REPORT( xSel95_spft );    // length at 95% sel by species/fleet/pop
  REPORT( xSelStep_spft );  // length at 95% sel by species/fleet/pop
  REPORT( xSel50_spf );     // length at 50% sel by species/fleet
  REPORT( xSel95_spf );     // length at 95% sel by species/fleet
  REPORT( xSelStep_spf );   // length at 95% sel by species/fleet
  REPORT( xSel50_sf );      // length at 50% sel by species/fleet
  REPORT( xSel95_sf );      // length at 95% sel by species/fleet
  REPORT( xSelStep_sf );    // length at 95% sel by species/fleet
  REPORT( F_spft );         // Fishing mortality
  REPORT( paF_spft );         // Fishing mortality
  REPORT( U_spft );         // Harvest rate
  REPORT( sel_lspft );      // Selectivity at length
  REPORT( sel_axspft );     // Selectivity at age
  REPORT( Fbar_spf );       // Average fishing mortality
  REPORT( tauD_f );
  REPORT( sigmaxSel50_sg );
  REPORT( sigmaxSelStep_sg );
  REPORT( propQ_spft );       // Proportion of catchability (learning to fish)


  // Model states
  REPORT( B_axspt );
  REPORT( N_axspt );
  REPORT( R_spt );
  REPORT( bhR_spt );
  REPORT( C_axspft );
  REPORT( paC_axspft );
  REPORT( Cw_axspft );
  REPORT( predCw_spft );
  REPORT( predC_spft );
  REPORT( vB_spft );
  REPORT( Nv_spft );
  REPORT( Nv_axspft );
  REPORT( SB_spt );
  REPORT( B_spt );
  REPORT( Z_axspt );
  REPORT( paZ_axspt );
  REPORT( endN_axspt );


  // Stock recruit parameters
  REPORT( Surv_axsp );
  REPORT( SSBpr_asp );
  REPORT( R0_sp );
  REPORT( phi_sp );
  REPORT( reca_sp );
  REPORT( recb_sp );

  // Growth and maturity
  REPORT( probLenAge_laxsp );
  REPORT( lenAge_axsp );
  REPORT( Wlen_ls );
  REPORT( matLen_ls );
  REPORT( matAge_asp );
  REPORT( xMat50_s );
  REPORT( xMat95_s );
  REPORT( meanWtAge_axsp );
  REPORT( L1_s );
  REPORT( L2_xs );
  REPORT( L2_sp );
  REPORT( vonK_xs );
  REPORT( vonK_sp );
  REPORT( A1_s );
  REPORT( A2_s );

  // Species and complex level hyperparameters
  REPORT( h_s );
  REPORT( h );
  REPORT( sigmah_s );
  REPORT( sigmah );

  // Random effects
  REPORT( epsxSel50_spft );
  REPORT( epsxSelStep_spft );
  REPORT( epslnq_spft );
  REPORT( omegaR_spt );
  REPORT( omegaRinit_asp );
  REPORT( omegaRmat_spt );
  REPORT( omegaRbar_sp );

  // Species/stock effects
  REPORT( epsM_sx );
  REPORT( epsM_sp );
  REPORT( epsM_s );
  REPORT( epsSteep_sp );
  REPORT( epsSteep_s );
  REPORT( sigmaM_s );
  REPORT( sigmaM );
  REPORT( sigmah );
  REPORT( sigmah_s );
  REPORT( sigmalnq );


  // Likelihood values
  REPORT( joint_nlp );
  REPORT( ageRes_axspft );
  REPORT( SRnlp_sp );
  REPORT( tau2Age_spf );
  REPORT( tau2PropF_spf );
  REPORT( propFnll_spf );
  REPORT( lenRes_lxspft );
  REPORT( tau2Len_spf );
  REPORT( recCorrMat_sp );
  REPORT( tau2Obs_spf );
  REPORT( tauObs_spf );
  REPORT( L1nlp );
  REPORT( L2nlp );
  REPORT( VonKnlp );
  if( debugMode == 1 )
  {
    REPORT( recnlp );
    REPORT( vonBnll_spf );
    REPORT( CPUEnll_spf );
    REPORT( tauObsnlp_f );
    REPORT( residCPUE_spft );
    REPORT( ageCompsnll_spf );
    REPORT( etaSumSqAge_spf );
    REPORT( nResidsAge_spf );
    REPORT( nObsAge_spf );
    REPORT( lenCompsnll_spf );
    REPORT( etaSumSqLen_spf );
    REPORT( nResidsLen_spf );
    REPORT( nObsLen_spf );
    REPORT( Fnlp );
    REPORT( steepnessnlp_sp );
    REPORT( steepnessnlp_s );
    REPORT( steepnessnlp );
    REPORT( Mnlp_sp );
    REPORT( Mnlp_s );
    REPORT( Mnlp_sx );
    REPORT( Mnlp );
    REPORT( qnlp_tv );
    REPORT( qnlp_gps );
    REPORT( qnlp_stock );
    REPORT( minAgeProp );
    REPORT( minLenProp );
    REPORT( sel_nlp );
    REPORT( corrRecnlp );

    // REPORT( UBL1_s );
    // REPORT( UBL2_s );
    // REPORT( UBVonK_s );
    // REPORT( LBL1_s );
    // REPORT( LBL2_s );
    // REPORT( LBVonK_s );
  }

  // Switches and constants
  REPORT( group_f );
  REPORT( A_s );
  REPORT( minA_s );
  REPORT( L_s );
  REPORT( minL_s );
  REPORT( lenD_s );
  REPORT( calcIndex_spf );
  REPORT( tvSelFleets );
  REPORT( tvqFleets );
  REPORT( logitqFleets );
  REPORT( regFfleets );

  // Predicted observations
  REPORT( lDist_lxspft_hat );
  REPORT( aDist_axspft_hat );
  REPORT( I_spft_hat );
  REPORT( probAgeLen_laxspft );
  REPORT( ageAtLenResids_laxspft ); 

  // growth model
  REPORT( nObsAgeAtLen_spf );
  REPORT( etaSumSqAgeAtLen_spf );
  REPORT( tau2AgeAtLenObs_spf );
  REPORT( nAgeLenObs );


  // Compositional correlation matrices
  REPORT( CorrAge_sfaa );
  REPORT( CorrLen_sfll );

  // Vulnerable biomass at age
  REPORT( vB_axspft );
  REPORT( vB_xspft );

  // Effort model pars
  REPORT( qFVecIdx );
  REPORT( qF_spf );
  REPORT( qF_spft );
  REPORT( Ctnll_spf );
  REPORT( qFdev_nlp );


  // Echo switches
  if( debugMode ==1 )
  {
    REPORT( swRinit_s );
    REPORT( parSwitch );
    REPORT( calcIndex_spf );
    REPORT( tvSelFleets );
    REPORT( tvqFleets );
    REPORT( regFfleets );
    REPORT( idxLikeWt_g );
    REPORT( ageLikeWt_g );
    REPORT( lenLikeWt_g );
    REPORT( tFirstRecDev_s );
    REPORT( tLastRecDev_s );
    REPORT( minAgeProp );
    REPORT( minLenProp );
    REPORT( lambdaB0 );
    REPORT( minTimeIdx_spf );
    REPORT( maxSel_spf );
  }

  // Return simulated data
  SIMULATE{
    REPORT( simI_spft );
  }


  // sd report for standard errors
  ADREPORT(lnB0_sp);    
  // ADREPORT(logitSteep_sp);     
  // ADREPORT(lnM_xsp);     
  // ADREPORT(lnVonK_sp);  
  // ADREPORT(lnL1_sp); 
  // ADREPORT(lnsigmaL_s);
  // ADREPORT(lnsigmaR_sp);
  // ADREPORT(lnq_spf)          // Catchability
  // ADREPORT(lntauObs_spf)        // Observation error
  // ADREPORT(lnxSel50_sf)   // length at 50% sel
  // ADREPORT(lnxSelStep_sf)   // length at 95% sel


  return( f );

}

