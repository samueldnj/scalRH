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
//  5. Hierarchical selectivity:
//      a. Species/fleet for pooled compositions
//      b. Species/fleet/stock for stock specific but pooled over years?
//      c. Time-varying selectivity for species/fleet/stock/year
//  6. Refactor Mortality and Steepness parameters:
//      - Complex level means are leading pars
//      - Species and stock level values are modeled as random 
//        effects (deviations from higher level mean values)
//  7. Investigate whether integrating growth model will improve 
//      length composition fits.
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
// - Multi-level RH priors on:
//     - Growth (vonB or F-W undecided)
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
// - Discarding - need to talk to fisheries about discarding behaviour
//     - "Grading length" for fisheries with no size limit
// 
// ><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><><>><><>><><>><>

#include <TMB.hpp>                                // Links in the TMB libraries
#include <iostream>
// link personal library of TMB stock assessment functions
// #include "/Users/sdnjohnson/Work/code/tmbFuns/stockAssessmentFuns.hpp"  



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
  
  res = log(newObs) - log(newPred);

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




// objective function
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Call namespaces //
  using namespace density;

  /*data section*/
  // Data Structures
  DATA_ARRAY(I_spft);               // CPUE data
  DATA_ARRAY(C_spft);               // Catch data (biomass)
  DATA_ARRAY(D_spft);               // Discard data (biomass)
  DATA_ARRAY(ALK_spalft);           // Age-length observations (in freq) by pop, separated by fleet and year
  DATA_ARRAY(age_aspft);            // Age observations for each population and fleet over time (-1 missing)
  DATA_ARRAY(len_lspft);            // Length observations for each population and fleet over time (-1 missing)
  DATA_IVECTOR(group_f);            // Fleet group (0=comm, 1=HSAss,2=Synoptic)
  DATA_IVECTOR(A_s);                // +Group age by species
  DATA_IVECTOR(L_s);                // +Group length by species
  DATA_IVECTOR(minA_s);             // Min age class by species
  DATA_IVECTOR(minL_s);             // Min length bin (cm) by species
  DATA_IVECTOR(lenD_s);             // Discard length by species

  // Model dimensions
  int nS = I_spft.dim(0);           // No. of species
  int nP = I_spft.dim(1);           // No. of stocks in each species
  int nF = I_spft.dim(2);           // No. of fleets (surveys + trawl)  
  int nT = I_spft.dim(3);           // No of time steps
  int nL = len_lspft.dim(0);        // max no of length bins (for creating state arrays)
  int nA = age_aspft.dim(0);        // max no of age bins (for creating state arrays)

  // Model switches - fill these in when we know what we want
  DATA_IVECTOR(swRinit_s);          // species fished initialisation switch (0 == unfished, 1 == fished)
  DATA_INTEGER(parSwitch);          // parallel accumulator switch for faster run time
  DATA_ARRAY(calcIndex_spf);        // Switches to calculate indices
  DATA_IVECTOR(tvSelFleets);        // Switch for fleet to have time varying selectivity (0 = off, 1 = on)
  DATA_IVECTOR(tvqFleets);          // Switch for fleet to have time varying catchability (0 = off, 1 = on)
  DATA_SCALAR(idxLikeWt);           // Scalar modifying the weight of the age likelihood
  DATA_SCALAR(ageLikeWt);           // Scalar modifying the weight of the age likelihood
  DATA_SCALAR(lenLikeWt);           // Scalar modifying the weight of the length likelihood
  DATA_SCALAR(growthLikeWt);        // Scalar modifying the weight of the growth model likelihood
  DATA_IVECTOR(tFirstRecDev_s);     // Initial rec devs for each species
  DATA_IVECTOR(tLastRecDev_s);      // Last recruitment devs for each species
  DATA_SCALAR(minAgeProp);          // Minimum observed proportion in age comps
  DATA_SCALAR(minLenProp);          // Minimum observed proportion in length comps
  DATA_SCALAR(minPAAL);             // Minimum observed proportion age at length
  DATA_STRING(matX);                // Age/Len switch for maturity model
  DATA_STRING(selX);                // Age/Len switch for selectivity model
  DATA_STRING(lenComps);            // Switch for computation of expected length comps
  DATA_SCALAR(lambdaB0);            // lambda scalar for exponential prior on B0
  DATA_ARRAY(minTimeIdx_spf);       // array of earliest times for each index

  // Fixed values
  DATA_IVECTOR(A1_s);               // Age at which L1_s is estimated
  DATA_IVECTOR(A2_s);               // Age at which L2_s is estimated


  /*parameter section*/
  // Leading Parameters
  // Biological
  PARAMETER_ARRAY(lnB0_sp);         // Biomass at unfished for ind. stocks
  PARAMETER(logitSteep);            // pop specific steepness
  PARAMETER(lnM);                   // stock-specific M
  PARAMETER_VECTOR(lnL2step_s);     // species-spec. step in Schnute-vonB L2 parameter
  PARAMETER_VECTOR(lnvonK_s);       // species-specific vonB K parameter
  PARAMETER_VECTOR(lnL1_s);         // species-specific vonB Length-at-age 1
  PARAMETER_ARRAY(deltaL2_sp);      // species-pop specific deviation in L2
  PARAMETER_ARRAY(deltaVonK_sp);    // species-pop specific deviation in VonK parameter
  PARAMETER_VECTOR(sigmaLa_s);      // species-specific individual growth SD intercept
  PARAMETER_VECTOR(sigmaLb_s);      // species-specific individual growth SD slope
  PARAMETER_VECTOR(LWa_s);          // species L-W conversion a par (units conv: cm -> kt)
  PARAMETER_VECTOR(LWb_s);          // species L-W conversion b par (length to volume)
  PARAMETER_VECTOR(xMat50_s);       // x (age/len) at 50% maturity
  PARAMETER_VECTOR(xMat95_s);       // x (age/len) at 95% maturity
  // Observation model
  // survey //
  PARAMETER_ARRAY(lnq_spf);         // fleet-species-stock specific catchability
  PARAMETER_ARRAY(lntau_spf);       // fleet-species-stock specific obs error sd
  // Fishery model
  // selectivity by fleet //
  PARAMETER_ARRAY(lnxSel50_sf);   // Selectivity Alpha parameter by species/fleet
  PARAMETER_ARRAY(lnxSelStep_sf); // Selectivity Beta parameter by species/fleet

  // Fishing mortality
  PARAMETER_VECTOR(lnF_spft);       // Fishing mortality as a long vector by species/population/fleet/time
  PARAMETER_VECTOR(lntauC_f);       // Catch observation SD by fleet
  PARAMETER_VECTOR(lntauD_f);       // Discards observation SD by fleet

  // Priors
  // selectivity //
  PARAMETER_ARRAY(muxSel50_sg);    // mean length at 50% selectivity by fleet group
  PARAMETER_ARRAY(muxSel95_sg);    // mean selectivity-at-length step by fleet group 
  PARAMETER_ARRAY(sigmaxSel50_sg); // SD in length at 50% selectivity by fleet group
  PARAMETER_ARRAY(sigmaxSel95_sg); // SD in length at 95% selectivity by fleet group

  // Count fleet groups
  int nGroups = muxSel95_sg.dim(1);

  // catchability //
  PARAMETER_VECTOR(lnqbarSyn_s);      // Synoptic survey q by species (each leg surveys one stock)
  PARAMETER_VECTOR(lntauqSyn_s);      // Synoptic survey q SD by species (each leg survey one stock)
  PARAMETER(lnqbarSyn);               // Synoptic survey complex q
  PARAMETER(lntauqSyn);               // Synoptic survey complex q SD
  PARAMETER(mqSurveys);               // prior mean survey q (HSAss included)
  PARAMETER(sdqSurveys);              // Prior SD survey q
  // steepness //
  PARAMETER_VECTOR(lnsigmah_s);       // species level steepness SD
  PARAMETER(logit_muSteep);           // complex steepness
  PARAMETER(lnsigmah);                // complex level steepness SD
  // Natural Mortality //
  PARAMETER_VECTOR(lnsigmaM_s);       // species M SD
  PARAMETER(ln_muM);                  // Multispecies assemblage mean M
  PARAMETER(lnsigmaM);                // Assemblage M SD
  // Index observation errors SD
  PARAMETER_VECTOR(IGatau_f);         // IG a parameter for tau_spf observation errors
  PARAMETER_VECTOR(IGbtau_f);         // IG b parameter for tau_spf observation errors


  // Random Effects
  // recruitment //
  PARAMETER_VECTOR(epsSteep_s);         // Species level steepness effect
  PARAMETER_VECTOR(epsM_s);             // Species level natural mortality effect
  PARAMETER_ARRAY(epsSteep_sp);         // stock level steepness effect
  PARAMETER_ARRAY(epsM_sp);             // stock level natural mortality effect
  PARAMETER_VECTOR(omegaR_vec);         // species-stock specific recruitment errors 2:nT
  PARAMETER_VECTOR(omegaRinit_vec);     // stock-age specific recruitment initialisation errors
  PARAMETER_ARRAY(lnsigmaR_sp);         // stock recruitment errors sd (sqrt cov matrix diag)
  PARAMETER_VECTOR(logitRCorr_chol);    // stock recruitment errors corr chol factor off diag  
  PARAMETER_ARRAY(logitRgamma_sp);      // species-stock-specific AR1 auto-corr on year effect (omega)
  PARAMETER_VECTOR(epsxSel50_vec);      // random effect on length at 50% selectivity
  PARAMETER_VECTOR(epsxSelStep_vec);    // random effect on length at 95% selectivity
  PARAMETER(lnsigmaSel);                // SD on time varying random effect for selectivity
  PARAMETER_VECTOR(epslnq_vec);         // random effect for time varying catchability parameter
  PARAMETER( lnsigmaepslnq );           // SD on time varying RE for lnq


  // Priors on growth and selectivity
  PARAMETER_ARRAY(pmlnxSel50_sf);       // Prior mean x-at-50% sel
  PARAMETER_ARRAY(pmlnxSelStep_sf);     // Prior mean xSelStep from 50% to 95% selectivity
  PARAMETER(cvxSel);                    // CV on xSel50/xSelStep prior
  PARAMETER_VECTOR(pmlnL2_s);           // Prior mean L2 par
  PARAMETER_VECTOR(pmlnL1_s);           // Prior mean L1 par
  PARAMETER(cvL2);                      // Prior CV on L2
  PARAMETER(cvL1);                      // Prior CV on L1
  PARAMETER(pmlnVonK);                  // Prior complex vonK 
  PARAMETER(cvVonK);                    // Prior CV on VonK
  // mortality deviations //
  /*
  PARAMETER_ARRAY(epsilonM_spt);    // M deviations by population, 2:nT
  PARAMETER_VECTOR(lnsigmaM_sp);    // M deviation pop-specific SD
  PARAMETER_VECTOR(logitMgamma_sp); // AR1 auto-corr on M devs effect (omega)
  */

  /* Procedure Section */
  // Exponentiate leading parameters
  // Biological //
  array<Type>  B0_sp(nS,nP);
  array<Type>  h_sp(nS,nP);
  array<Type>  M_sp(nS,nP);
  array<Type>  L2_sp(nS,nP);
  array<Type>  vonK_sp(nS,nP);
  array<Type>  L1_sp(nS,nP);
  array<Type>  sigmaR_sp(nS,nP);

  // Growth model vectors
  vector<Type>  L1_s(nS);
  vector<Type>  L2_s(nS);
  vector<Type>  vonK_s(nS);


  // derived variables
  // Stock recruitment //
  array<Type>  R0_sp(nS,nP);            // eqbm recruitment
  array<Type>  phi_sp(nS,nP);           // eqbm SSB per recruit
  array<Type>  reca_sp(nS,nP);          // BH a parameter for pops
  array<Type>  recb_sp(nS,nP);          // BH b parameter for pops

  // Growth //
  array<Type>   Wlen_ls(nL,nS);         // weight-at-length by species
  array<Type>   lenAge_asp(nA,nS,nP);   // Mean length-at-age by population
  array<Type>   probLenAge_lasp(nL,nA,nS,nP);
  array<Type>   meanWtAge_asp(nA,nS,nP);
  array<Type>   probAgeLen_alspft(nA,nL,nS,nP,nF,nT);

  // Maturity //
  array<Type>   matAge_as(nA,nS);       // Proportion mature at age by species
  array<Type>   matAge_asp(nA,nS,nP);   // Proportion mature at age by population
  array<Type>   matLen_ls(nL,nS);       // Proportion mature at Length by species

  
  
  // Make a vector of ages/lengths for use in 
  // age and length based quantities later
  vector<Type> age(nA);
  vector<Type> len(nL);
  for( int l = 0; l < nL; l++)
    len(l) = l+1;

  for( int a = 0; a < nA; a++)
    age(a) = a+1;



  // Prior Hyperparameters //
  // Steepness
  vector<Type>  h_s           = .2 + Type(0.8) / ( Type(1.0) + exp( -1 * (logitSteep + epsSteep_s ) ) );
  vector<Type>  sigmah_s      = exp(lnsigmah_s);
  Type          mh            = .2 + Type(0.8) / ( Type(1.0) + exp( -1 * logit_muSteep) );
  Type          sigmah        = exp(lnsigmah);
  // Natural mortality
  vector<Type>  M_s           = exp(lnM + epsM_s);
  vector<Type>  sigmaM_s      = exp(lnsigmaM_s);
  Type          muM           = exp(ln_muM);
  Type          sigmaM        = exp(lnsigmaM);
  // RE SDs
  Type          sigmaSel      = exp(lnsigmaSel);
  Type          sigmalnq      = exp(lnsigmaepslnq);
  // Catchability for synoptic survey
  vector<Type>  qbarSyn_s     = exp(lnqbarSyn_s);
  vector<Type>  tauqSyn_s     = exp(lnqbarSyn_s);
  Type          qbarSyn       = exp(lnqbarSyn);
  Type          tauqSyn       = exp(lntauqSyn);

  L1_s    = exp(lnL1_s); 
  vonK_s  = exp(lnvonK_s);
  L2_s    = exp(lnL1_s) + exp(lnL2step_s);
  

  // Fill arrays that hold stock-specific parameters
  for( int pIdx = 0; pIdx < nP; pIdx++ )
  {
    B0_sp.col(pIdx)      = exp(lnB0_sp.col(pIdx));
    h_sp.col(pIdx)       = .2 + 0.8 / ( Type(1.0) + exp( -1 * (logitSteep + epsSteep_s + epsSteep_sp.col(pIdx) ) ) );
    M_sp.col(pIdx)       = exp(lnM) * exp(epsM_s) * exp(epsM_sp.col(pIdx) );
    L1_sp.col(pIdx)      = L1_s;
    sigmaR_sp.col(pIdx)  = exp(lnsigmaR_sp.col(pIdx));
    L2_sp.col(pIdx)      = L2_s * exp(deltaL2_sp.col(pIdx));
    vonK_sp.col(pIdx)    = vonK_s * exp(deltaVonK_sp.col(pIdx));
  }


  // Observation model
  array<Type>   q_spf(nS,nP,nF);
  array<Type>   tau_spf(nS,nP,nF);
  q_spf.setZero();
  tau_spf.setZero();

  // Time-varying catchability
  array<Type>   q_spft(nS,nP,nF,nT);

  // selectivity //
  // Fleet/species
  array<Type>   xSel50_sf(nS,nF);
  array<Type>   xSelStep_sf(nS,nF);
  array<Type>   xSel95_sf(nS,nF);
  // Now time-varying by stock
  array<Type>   xSel50_spft(nS,nP,nF,nT);
  array<Type>   xSelStep_spft(nS,nP,nF,nT);
  array<Type>   xSel95_spft(nS,nP,nF,nT);
  array<Type>   epsxSel50_spft(nS,nP,nF,nT);
  array<Type>   epsxSelStep_spft(nS,nP,nF,nT);
  array<Type>   epslnq_spft(nS,nP,nF,nT);
  // Zero-initialise
  xSel50_sf.setZero();
  xSel95_sf.setZero();
  xSelStep_sf.setZero();
  xSel50_spft.setZero();
  xSel95_spft.setZero();
  xSelStep_spft.setZero();
  epsxSel50_spft.setZero();
  epsxSelStep_spft.setZero();
  epslnq_spft.setZero();


  // Loop over gears, species, stocks and time steps, create a matrix
  // of deviations in selectivity pars and calculate the pars
  // for each time step
  int epsSelVecIdx = 0;
  int epslnqVecIdx = 0;
  for(int f = 0; f < nF; f++ )
  {
    xSel50_sf.col(f)    = exp(lnxSel50_sf.col(f));
    xSelStep_sf.col(f)  = exp(lnxSelStep_sf.col(f));
    xSel95_sf.col(f)    = xSel50_sf.col(f) + xSelStep_sf.col(f);
    for( int s = 0; s < nS; s++)
      for( int p = 0; p < nP; p++ )
      {
        // Exponentiate q and tau for the index observations
        q_spf(s,p,f) = exp(lnq_spf(s,p,f));
        tau_spf(s,p,f) = exp(lntau_spf(s,p,f));

        xSel50_spft(s,p,f,0)      = xSel50_sf(s,f);
        xSelStep_spft(s,p,f,0)    = xSelStep_sf(s,f);
        q_spft(s,p,f,0)           = q_spf(s,p,f);
        for( int t = 0; t < nT; t++ )
        { 
          // Bring previous time-step selecivity pars
          // forward
          if( t > 0 )
          {
            xSel50_spft(s,p,f,t)    = xSel50_spft(s,p,f,t-1);
            xSelStep_spft(s,p,f,t)  = xSelStep_spft(s,p,f,t-1);
            q_spft(s,p,f,t)         = q_spft(s,p,f,t-1);
          }
          // Update the epsilon array          
          if( tvSelFleets(f) == 1)
            if( ((age_aspft(0,s,p,f,t) >= 0) & (ageLikeWt > 0)) |
                ((len_lspft(0,s,p,f,t) >= 0) & (lenLikeWt > 0)) )
            {
              epsxSel50_spft(s,p,f,t)   += epsxSel50_vec(epsSelVecIdx);
              epsxSelStep_spft(s,p,f,t) += epsxSelStep_vec(epsSelVecIdx);
              epsSelVecIdx++;
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
          
        }
      }
  }
  // Compute xSel95 for the time-varying arrays
  xSel95_spft = xSel50_spft + xSelStep_spft;
  
  // Fishing mortality
  array<Type>   F_spft(nS,nP,nF,nT);
  F_spft.setZero();
  int vecIdx = 0;
  for( int s = 0; s < nS; s++ )
    for( int f = 0; f < nF; f++ )
      for( int p = 0; p < nP; p++ )
        for( int t = 0; t < nT; t++ )
          if( C_spft(s,p,f,t) > 0)
          {
            F_spft(s,p,f,t) = exp( lnF_spft(vecIdx) );
            vecIdx++;
          }

  // recruitment deviations and initRecDevs - convert from vector
  // to array of deviations
  array<Type> omegaRinit_spa(nS,nP,nA);
  omegaRinit_spa.setZero();
  array<Type> omegaR_spt(nS,nP,nT);
  omegaR_spt.setZero();
  int initVecIdx = 0;
  int devVecIdx = 0;
  for( int s = 0; s < nS; s++ )
    for( int p = 0; p < nP; p++ )
    {
      for( int t = tFirstRecDev_s(s); t < tLastRecDev_s(s); t++ )
      {
        omegaR_spt(s,p,t) = omegaR_vec(devVecIdx);
        devVecIdx++;
      }
      if(swRinit_s(s) == 1)
        for( int a = 0; a < A_s(s); a++ )
        {
          omegaRinit_spa(s,p,a) = omegaRinit_vec(initVecIdx);
          initVecIdx++;
        }
    }

  
  // Probably concentrate these out later....
  vector<Type>  tauC_f        = exp(lntauC_f);    
  vector<Type>  tauD_f        = exp(lntauD_f);    
  
  // Set up model state arrays
  array<Type> B_aspt(nA,nS,nP,nT);          // Biomass at age-pop-time
  array<Type> N_aspt(nA,nS,nP,nT);          // Numbers at age-pop-time
  array<Type> R_spt(nS,nP,nT);              // Recruits by pop-time
  array<Type> F_aspft(nA,nS,nP,nF,nT);      // Fishing mortality by age, fleet, population and time
  array<Type> Z_aspt(nA,nS,nP,nT);          // Total mortality by age, population and time
  array<Type> C_aspft(nA,nS,nP,nF,nT);      // Predicted catch-at-age (in numbers), fleet, population and time
  array<Type> Cw_aspft(nA,nS,nP,nF,nT);     // Predicted catch-at-age (in weight), population and time
  array<Type> predCw_spft(nS,nP,nF,nT);     // Predicted catch (in weight) by population and time
  array<Type> predC_spft(nS,nP,nF,nT);      // Predicted catch (in numbers) by population and time
  array<Type> B_spt(nS,nP,nT);              // Total biomass by species, pop, time
  array<Type> SB_spt(nS,nP,nT);             // Spawning biomass by species, pop, time
  array<Type> Bv_spft(nS,nP,nF,nT);         // Vulnerable biomass by species, pop, time
  array<Type> Nv_spft(nS,nP,nF,nT);         // Vulnerable numbers by species, pop, time
  array<Type> Nv_aspft(nA,nS,nP,nF,nT);     // Vulnerable numbers at age by species, pop, time
  array<Type> Surv_asp(nA,nS,nP);           // Eqbm survival at age
  array<Type> SSBpr_asp(nA,nS,nP);          // Spawning stock biomass per rec. at age



  // ---------------- Growth, Maturity, Stock-Recruitment ----------- //
  Wlen_ls.fill(-1.0);
  matLen_ls.fill(-1.0);
  for( int s = 0; s < nS; s++ )
  {
    for( int l = 0; l < L_s(s); l++ )
    {
      // weight-at-length for each length bin midpoint
      Wlen_ls(l,s) = LWa_s(s) * pow( l+1, LWb_s(s) );
      // proportion mature-at-length
      if(matX == "length")
        matLen_ls(l,s) = 1/(1 + exp( -1. * log(Type(19.0)) * (l+1 - xMat50_s(s)) / (xMat95_s(s) - xMat50_s(s)) ) );
      
    }
  }

  // Now produce the probability of length-at-age matrix for each stock
  // (taken from LIME by Rudd and Thorson, 2017)
  // In the same loop, we're going to compute the ssbpr parameter phi
  // for each stock
  lenAge_asp.setZero();
  probLenAge_lasp.setZero();
  meanWtAge_asp.setZero();
  phi_sp.setZero();
  Surv_asp.setZero();
  SSBpr_asp.setZero();
  matAge_asp.setZero();

  for( int s = 0; s < nS; s++)
  {
    int A1 = A1_s(s);
    int A2 = A2_s(s);
    for( int p = 0; p < nP; p++ )
    {
      Type tmpPhi = 0;
      for( int a = 0; a < A_s(s); a ++ )
      {
        Type  sumProbs  = 0;
        Type vonK_tmp = vonK_sp(s,p);
        lenAge_asp(a,s,p) += L2_sp(s,p) - L1_sp(s,p);
        lenAge_asp(a,s,p) *= (exp(-vonK_tmp * A1 ) - exp(-vonK_tmp * (a+1)) );
        lenAge_asp(a,s,p) /= (exp(-vonK_tmp * A1 ) - exp(-vonK_tmp * A2) );
        lenAge_asp(a,s,p) += L1_sp(s,p);
        Type sigmaL = sigmaLa_s(s) + sigmaLb_s(s) * lenAge_asp(a,s,p);
        for( int l = 0; l < L_s(s); l++ )
        {
          // Get length bin ranges
          Type len = l + 1;
          Type lenHi = len + .5;
          Type lenLo = len - .5;
          // Compute density in each bin
          // Compute LH tail of the distribution
          
          if( (l >= 0) & (l < L_s(s) - 1 ) )
          {
            probLenAge_lasp(l,a,s,p) = pnorm( log(lenHi), log(lenAge_asp(a,s,p)), sigmaL );
            // Subtract LH tail for interior bins
            if(l > 0)
              probLenAge_lasp(l,a,s,p) -= pnorm( log(lenLo), log(lenAge_asp(a,s,p)), sigmaL );
            // Accumulate probability
            sumProbs += probLenAge_lasp(l,a,s,p);
          }
          
          // Now the RH tail
          if( l == L_s(s) - 1 ) 
            probLenAge_lasp(l,a,s,p) = Type(1.0) - sumProbs;

          // Calculate mean weight at age from probLenAge
          meanWtAge_asp(a,s,p) += probLenAge_lasp(l,a,s,p) * Wlen_ls(l,s);
          // Calculate maturity at age from probLenAge if maturity
          // is length based
          if( matX == "length")
            if( (matLen_ls(l,s) > 0) )
              matAge_asp(a,s,p) += probLenAge_lasp(l,a,s,p) * matLen_ls(l,s);


        }

        // Calculate maturity at age if age based
        if( matX == "age" )
          matAge_asp(a,s,p) = 1 / (1 + exp( -1. * log(Type(19.0)) * ( a+1 - xMat50_s(s)) / (xMat95_s(s) - xMat50_s(s)) ) );

        // To compute ssbpr, we need to borrow N_aspt for a moment 
        Surv_asp(a,s,p) = exp(-a * M_sp(s,p));        
        if( a == A_s(s) - 1 ) 
          Surv_asp(a,s,p) /= (1. - exp( -M_sp(s,p)) );
        // Compute ssbpr
        
      }
    }
  }

  SSBpr_asp = Surv_asp * matAge_asp * meanWtAge_asp;

  // --------- Selectivity and Fishing Mortality -------- //
  array<Type> sel_lfspt(nL,nF,nS,nP,nT);
  array<Type> sel_afspt(nA,nF,nS,nP,nT);
  sel_lfspt.setZero();
  sel_afspt.setZero();
  for( int t = 0; t < nT; t++ )
    for( int s = 0; s < nS; s++ )
      for( int p = 0; p < nP; p++ )
      {
        // Loop over fleets, create selectivities
        for( int f = 0; f < nF; f++ )
        {
          Type maxSel = 0.0;
          // selectivity-at-length
          if( selX == "length" )
          {
            for( int l = 0; l < L_s(s); l++ )
            {
              sel_lfspt(l,f,s,p,t) = 1;
              sel_lfspt(l,f,s,p,t) /= ( 1 + exp( - log(Type(19.)) * ( l + 1 - xSel50_spft(s,p,f,t) ) / xSelStep_spft(s,p,f,t)  ) );
            }
            // convert selectivity-at-length to selectivity-at-age using probability matrix
            for( int a = 0; a < A_s(s); a++ )
            {
              vector<Type> probLenAgea_l(nL) ;
              probLenAgea_l.setZero();
              probLenAgea_l = probLenAge_lasp.col(p).col(s).col(a);
              sel_afspt(a,f,s,p,t) = (probLenAgea_l*sel_lfspt.col(t).col(p).col(s).col(f)).sum();

              if( sel_afspt(a,f,s,p,t) > maxSel )
                maxSel = sel_afspt(a,f,s,p,t);
            } 
          }
          if( selX == "age" )
          {
            for( int a = 0; a < A_s(s); a++)
            {
              sel_afspt(a,f,s,p,t) = 1/( 1 + exp( - log(Type(19.)) * ( a + 1 - xSel50_spft(s,p,f,t) ) / xSelStep_spft(s,p,f,t)  ) );
              if( sel_afspt(a,f,s,p,t) > maxSel )
                maxSel = sel_afspt(a,f,s,p,t);
            }
          }
          // if(maxSel > 0)
          //   sel_afspt.col(t).col(p).col(s).col(f) /= maxSel;
        }
      }

  // --------- Population Dynamics ----------- //
  // First, compute eqbm recruitment
  // Derived parameters
  
  for( int p = 0; p < nP; p++ )
  {
    for( int s = 0; s < nS; s++)
      phi_sp(s,p) = SSBpr_asp.col(p).col(s).sum();

    R0_sp.col(p) = B0_sp.col(p) / phi_sp.col(p);
    reca_sp.col(p) = 4.*h_sp.col(p)*R0_sp.col(p) / ( B0_sp.col(p)*(1.-h_sp.col(p)) );
    recb_sp.col(p) = (5.*h_sp.col(p)-1.) / ( B0_sp.col(p)*(1.-h_sp.col(p)) );
  }

  // Set all state arrays to zero
  B_aspt.setZero();
  R_spt.setZero();
  F_aspft.setZero();
  Z_aspt.setZero();
  C_aspft.setZero();
  Cw_aspft.setZero();
  B_spt.setZero();
  SB_spt.setZero();
  Bv_spft.setZero();
  Nv_spft.setZero();
  Nv_aspft.setZero();
  N_aspt.setZero();

  // Loop over species
  for( int s = 0; s < nS; s++ )
  {
    // Loop over stocks
    for( int p = 0; p < nP; p++ )
    {
      // Set plus group age
      int A    = A_s(s);
      for( int t = 0; t < nT; t++ )
      {
        // Initial time-step
        if( t == 0 )
        {
          // Initialise population either at unfished
          // or using the estimated fished initialisation
          // Populate first year
          N_aspt.col(t).col(p).col(s) += R0_sp(s,p) * Surv_asp.col(p).col(s);
          for( int a = 0; a < A; a++ )
          {            
            // Produce init eqbm age structure
            // N_aspt(a,s,p,t) += R0_sp(s,p) * Surv_asp(a,s,p);
            // // Non-eqbm initialisation
            if(swRinit_s(s) == 1)
              N_aspt(a,s,p,t) *= exp( omegaRinit_spa(s,p,a));
            
            // Compute total mortality at age (for catch later)
            Z_aspt(a,s,p,t) = M_sp(s,p);
            for( int f = 0; f < nF; f++ )
            {
              // Compute fishing mortality at age
              F_aspft(a,s,p,f,t) = sel_afspt(a,f,s,p,t) * F_spft(s,p,f,t);
              // Add to Z
              Z_aspt(a,s,p,t) += F_aspft(a,s,p,f,t);  
            }
          }            
        }

        // time series history
        if( t > 0 )
        {
          // Generate recruitment
          Type SBt = SB_spt(s,p,t-1);
          N_aspt(0,s,p,t) = reca_sp(s,p) * SBt / (1 + recb_sp(s,p) * SBt);
          N_aspt(0,s,p,t) *= exp( omegaR_spt(s,p,t) );

          // Now loop over ages and apply fishing mortality (no discarding yet)
          for( int a = 0; a < A; a ++ )
          {
            // Compute fishing mortality at age
            // Compute total mortality at age
            Z_aspt(a,s,p,t) = M_sp(s,p);
            for( int f = 0; f < nF; f++ ) 
            {
              F_aspft(a,s,p,f,t) = sel_afspt(a,f,s,p,t) * F_spft(s,p,f,t);
              Z_aspt(a,s,p,t) += F_aspft(a,s,p,f,t);
            }

            // Update numbers according to age
            if( a > 0 )     
              N_aspt(a,s,p,t) = N_aspt( a-1, s, p, t-1 ) *  exp( - Z_aspt(a,s,p,t-1) );
            if( a == A - 1) 
              N_aspt(a,s,p,t) += N_aspt(a,s, p, t-1 ) *  exp( - Z_aspt(a,s,p,t-1) );
          }
        }

        // Save recruits in R_pt
        R_spt(s,p,t) = N_aspt(0,s,p,t);

        // Loop over fleets and compute catch
        for( int f =0; f < nF; f++ )
        {
          predCw_spft(s,p,f,t) = 0.;
          predC_spft(s,p,f,t) = 0.;

          // Calculate vulnerable numbers and biomass for this fleet
          Nv_aspft.col(t).col(f).col(p).col(s) = N_aspt.col(t).col(p).col(s) * sel_afspt.col(t).col(p).col(s).col(f);
          Nv_spft(s,p,f,t) = Nv_aspft.col(t).col(f).col(p).col(s).sum();

          // Have to loop over ages here to avoid adding NaNs
          for(int a = 0; a < A; a++ )
          {
            C_aspft(a,s,p,f,t) = N_aspt(a,s,p,t);
            C_aspft(a,s,p,f,t) *= (Type(1) - exp( -Z_aspt(a,s,p,t) )) * F_aspft(a,s,p,f,t) / Z_aspt(a,s,p,t);  
            Cw_aspft(a,s,p,f,t) = C_aspft(a,s,p,f,t) * meanWtAge_asp(a,s,p);
            predCw_spft(s,p,f,t) += Cw_aspft(a,s,p,f,t);
            predC_spft(s,p,f,t) += C_aspft(a,s,p,f,t);
            Bv_spft(s,p,f,t) += Nv_aspft(a,s,p,f,t) * meanWtAge_asp(a,s,p);

            // Compute probAgeAtLen for each time step
            for( int l = minL_s(s) - 1; l < L_s(s); l++ )
              probAgeLen_alspft(a,l,s,p,f,t) += probLenAge_lasp(l,a,s,p) * N_aspt(a,s,p,t) * sel_lfspt(l,f,s,p,t);
              

          }
          // Renormalise probAgeAtLen
          for( int l = minL_s(s) - 1; l < L_s(s); l++ )
          {
            Type sumProbs = probAgeLen_alspft.col(t).col(f).col(p).col(s).col(l).sum();
            probAgeLen_alspft.col(t).col(f).col(p).col(s).col(l) /= sumProbs;
          }
          

        }

        // Compute total and spawning biomass
        B_aspt.col(t).col(p).col(s) = N_aspt.col(t).col(p).col(s) * meanWtAge_asp.col(p).col(s);
        B_spt(s,p,t) = B_aspt.col(t).col(p).col(s).sum();
        SB_spt(s,p,t) = (B_aspt.col(t).col(p).col(s) * matAge_asp.col(p).col(s)).sum();


      }

    }
  }

  
  // --------- Observation Model --------- //
  // Age/Length observations and CPUE
  array<Type> I_spft_hat(nS,nP,nF,nT);
  array<Type> aDist_aspft_hat(nA,nS,nP,nF,nT);
  array<Type> lDist_lspft_hat(nL,nS,nP,nF,nT);
  I_spft_hat.fill(-1);
  for(int s = 0; s < nS; s++ )
    for( int p = 0; p < nP; p++ )
    {
      int A      = A_s(s);
      int L      = L_s(s);
      for( int f = 0; f < nF; f++ )
      {
        for( int t = 0; t < nT; t++ )
        {
          // Predict CPUE for surveys
          if( calcIndex_spf(s,p,f) == 1)
          {
            I_spft_hat(s,p,f,t) = q_spft(s,p,f,t) * Bv_spft(s,p,f,t);  
          }

          // Create array of predicted age distributions - this should just be catch-at-age props
          // Calculate total catch
          Type totCatch = C_aspft.col(t).col(f).col(p).col(s).sum();
            
          // Convert catch-at-age to proportions-at-age
          if(totCatch > 0)
            aDist_aspft_hat.col(t).col(f).col(p).col(s) = C_aspft.col(t).col(f).col(p).col(s) / totCatch;
              

          // Create array of predicted length distributions
          // Need to do some thinking here... follow LIME to start with...
          // Probability of being harvested at age
          vector<Type> probHarvAge(A);
          probHarvAge.setZero();
          probHarvAge += Nv_aspft.col(t).col(f).col(p).col(s) / N_aspt.col(t).col(p).col(s).sum();
          
          // Probability of sampling a given length bin
          vector<Type> probHarvLen(L);
          probHarvLen.setZero();
          // loop over length and age, calculate probability of harvest at length
          for( int l = 0; l < L; l++ )
          {
            if( lenComps == "LIME" )
            {
              for( int a = 0; a < A; a++ )
                probHarvLen(l) += probHarvAge(a)*probLenAge_lasp(l,a,s,p);
            }

            if( lenComps == "Francis" )
            {
              for( int a = 0; a < A; a++ )
                probHarvLen(l) += probLenAge_lasp(l,a,s,p) * N_aspt(a,s,p,t) * sel_lfspt(l,f,s,p,t); 
            }

            // Save to length dists
            lDist_lspft_hat(l,s,p,f,t) = probHarvLen(l);
          }
          // renormalise
          if(probHarvLen.sum() > 0)
            lDist_lspft_hat.col(t).col(f).col(p).col(s) /= probHarvLen.sum();
        }
      }
    }

  // --------- Statistical Model ----------- //

  // Process model likelihood //
  // Recruitment deviations (initial and ongoing) - add correlation later
  // both auto-correlation and between-pop correlation can be added.
  array<Type> recnll_sp(nS,nP);
  recnll_sp.setZero();
  for( int s = 0; s < nS; s++ )
    for( int p = 0; p < nP; p++ )
    {
      // First initialisation deviations
      int A = A_s(s);
      for( int a = 0; a < A; a++)
        recnll_sp(s,p) -= dnorm( omegaRinit_spa(s,p,a),Type(0),sigmaR_sp(s,p),true);

      // then yearly recruitment deviations
      for( int t = 1; t <  nT; t++ )
        recnll_sp(s,p) -= dnorm( omegaR_spt(s,p,t-1), Type(0.), sigmaR_sp(s,p),true);
    }

  // vonB growth model likelihood //
  // Currently has lengths binned the same as the integration above, 
  // and uses a random-at-length likelihood to account for differences 
  // in fleet selectivity, and at different times
  array<Type> vonBnll_spf(nS,nP,nF);
  array<Type> ageAtLenResids_alspft(nA,nL,nS,nP,nF,nT);
  array<Type> nObsAgeAtLen_spf(nS,nP,nF);
  array<Type> nResidsAgeAtLen_spf(nS,nP,nF);
  array<Type> etaSumSqAgeAtLen_spf(nS,nP,nF);
  array<Type> tau2AgeAtLenObs_spf(nS,nP,nF);
  // containers for a given year/fleet/length/species/pop - function
  // expects vector<Type>
  vector<Type> obs(nA);
  vector<Type> pred(nA);
  vector<Type> resids(nA);

  // Zero-init all arrays
  tau2AgeAtLenObs_spf.setZero();
  nObsAgeAtLen_spf.setZero();
  nResidsAgeAtLen_spf.setZero();
  etaSumSqAgeAtLen_spf.setZero();
  vonBnll_spf.setZero();

  // Now loop and compute
  for( int s = 0; s < nS; s++ )
    for( int p = 0; p < nP; p++ )
      for( int f = 0; f < nF; f++)
      {
        for( int t = 0; t < nT; t++)
        {
          for( int l = minL_s(s) - 1; l < L_s(s); l++ )
          {
            // Zero init containers
            obs.setZero();
            pred.setZero();
            resids.setZero();

            // Fill containers
            for(int a = 0; a < nA; a++)
            {
              obs(a) += ALK_spalft(s,p,a,l,f,t);
              pred(a) += probAgeLen_alspft(a,l,s,p,f,t);
            }
            // Compute resids etc. if observations aren't missing
            if(obs.sum() > 0 &  pred.sum() > 0)
            {
              resids = calcLogistNormLikelihood(  obs, 
                                                  pred,
                                                  minPAAL,
                                                  etaSumSqAgeAtLen_spf(s,p,f),
                                                  nResidsAgeAtLen_spf(s,p,f) );
              nObsAgeAtLen_spf(s,p,f) += 1;

              ageAtLenResids_alspft.col(t).col(f).col(p).col(s).col(l) += resids;
            }
          } 
        }
        if( nResidsAgeAtLen_spf(s,p,f) > 0)
        {
          tau2AgeAtLenObs_spf(s,p,f) += etaSumSqAgeAtLen_spf(s,p,f) / nResidsAgeAtLen_spf(s,p,f);
          vonBnll_spf(s,p,f) += 0.5 * (nResidsAgeAtLen_spf(s,p,f) - nObsAgeAtLen_spf(s,p,f));
          vonBnll_spf(s,p,f) *= log(tau2AgeAtLenObs_spf(s,p,f));
        }
      }

  // Observation likelihoods
  array<Type> CPUEnll_spft(nS,nP,nF,nT);  // CPUE
  array<Type> ageCompsnll_spf(nS,nP,nF);  // Age observations
  array<Type> lenCompsnll_spf(nS,nP,nF);  // length observations
  array<Type> Ctnll_sp(nS,nP);            // Catch
  array<Type> Dtnll_sp(nS,nP);            // Discards
  // Conditional MLEs of age/length observation error variance
  array<Type> tau2Age_spf(nS,nP,nF);
  array<Type> tau2Len_spf(nS,nP,nF);
  array<Type> etaSumSqLen_spf(nS,nP,nF);
  array<Type> etaSumSqAge_spf(nS,nP,nF);
  array<Type> nResidsAge_spf(nS,nP,nF);
  array<Type> nResidsLen_spf(nS,nP,nF);
  array<Type> nObsAge_spf(nS,nP,nF);
  array<Type> nObsLen_spf(nS,nP,nF);
  array<Type> ageRes_aspft(nA,nS,nP,nF,nT);
  array<Type> lenRes_lspft(nL,nS,nP,nF,nT);
  // Zero-init
  CPUEnll_spft.setZero();
  ageCompsnll_spf.setZero();
  lenCompsnll_spf.setZero();
  Ctnll_sp.setZero();
  Dtnll_sp.setZero();
  tau2Age_spf.setZero();
  tau2Len_spf.setZero();
  etaSumSqLen_spf.setZero();
  etaSumSqAge_spf.setZero();
  nResidsAge_spf.setZero();
  nResidsLen_spf.setZero();
  nObsAge_spf.setZero();
  nObsLen_spf.setZero();
  ageRes_aspft.setZero();
  lenRes_lspft.setZero();

  for( int s = 0; s < nS; s++ )
    for( int p = 0; p < nP; p++ )
    {
      int A      = A_s(s) - minA_s(s) + 1;
      int L      = L_s(s) - minL_s(s) + 1;
      // Loop over fleets
      for( int f = 0; f < nF; f++ )
      {
        // tmp vectors to hold age and length observed and
        // predicted values
        vector<Type> fleetAgeObs(A);
        vector<Type> fleetAgePred(A);
        vector<Type> fleetAgeResids(A);

        vector<Type> fleetLenObs(L);
        vector<Type> fleetLenPred(L);  
        vector<Type> fleetLenResids(L);   


        // Now loop over time steps

        for( int t = 0; t < nT; t++ )
        {

          // Set tmp vectors to zero 
          fleetAgeObs.setZero();
          fleetAgePred.setZero();
          fleetAgeResids.setZero();
          fleetLenObs.setZero();
          fleetLenPred.setZero();
          fleetLenResids.setZero();

          // Only use a year if the data exists and that fleet is used as a survey
          if( (I_spft(s,p,f,t) > 0) & (calcIndex_spf(s,p,f) > 0) )
            CPUEnll_spft(s,p,f,t) += -dnorm( log(I_spft(s,p,f,t)), log(I_spft_hat(s,p,f,t)),tau_spf(s,p,f),true);

          // Now compute the catch and discards likelihoods
          if( C_spft(s,p,f,t) > 0 & predCw_spft(s,p,f,t) > 0 )
            Ctnll_sp(s,p) += - dnorm( log(C_spft(s,p,f,t)), log(predCw_spft(s,p,f,t)), tauC_f(f),true);

          // Loop and fill age obs and pred matrices
          for( int a = 0; a < A; a++ )
          {
            fleetAgeObs(a)  = age_aspft(minA_s(s) + a - 1,s,p,f,t);
            fleetAgePred(a) = aDist_aspft_hat(minA_s(s) + a - 1,s,p,f,t);
          }

          // Loop and fill length obs and pred matrices
          for( int l = 0; l < L; l++ )
          {
            fleetLenObs(l)  = len_lspft(minL_s(s) + l - 1,s,p,f,t);
            fleetLenPred(l) = lDist_lspft_hat(minL_s(s) + l,s,p,f,t);
          }

          // Now for compositional data. First, check if ages are 
          // being used (or exist), if so
          // compute logistic normal likelihood
          if( (fleetAgeObs.sum() > 0) & (fleetAgePred.sum() > 0) & (ageLikeWt > 0) )
          {
            fleetAgeResids = calcLogistNormLikelihood(  fleetAgeObs, 
                                                        fleetAgePred,
                                                        minAgeProp,
                                                        etaSumSqAge_spf(s,p,f),
                                                        nResidsAge_spf(s,p,f) );

            // Now place resids in the ageRes array
            for( int a = minA_s(s) - 1; a < A_s(s); a++)
              ageRes_aspft(a,s,p,f,t) = fleetAgeResids(a - minA_s(s) + 1);

            // Increment number of ages
            nObsAge_spf(s,p,f) += 1;

          }     

          // If ages aren't being used, but lengths exist, then compute 
          // logistic normal likelihood
          if( (fleetLenObs.sum() > 0) & (fleetLenPred.sum() > 0) &  (lenLikeWt > 0) )
          {
            fleetLenResids = calcLogistNormLikelihood(  fleetLenObs, 
                                                        fleetLenPred,
                                                        minLenProp,
                                                        etaSumSqLen_spf(s,p,f),
                                                        nResidsLen_spf(s,p,f) );
            
            for( int l = minL_s(s) - 1; l < L_s(s); l++)
              lenRes_lspft(l,s,p,f,t) = fleetLenResids(l - minL_s(s) + 1);
            
            nObsLen_spf(s,p,f) += 1;
          }     
          
        }
        if( nResidsAge_spf(s,p,f) > 0)
        {
          tau2Age_spf(s,p,f)      += etaSumSqAge_spf(s,p,f) / nResidsAge_spf(s,p,f);
          ageCompsnll_spf(s,p,f)  += 0.5 * (nResidsAge_spf(s,p,f) - nObsAge_spf(s,p,f)) * log(tau2Age_spf(s,p,f));
        }
        if( nResidsLen_spf(s,p,f) > 0)
        {
          tau2Len_spf(s,p,f)      += etaSumSqLen_spf(s,p,f) / nResidsLen_spf(s,p,f);
          lenCompsnll_spf(s,p,f)  += 0.5 * (nResidsLen_spf(s,p,f) - nObsLen_spf(s,p,f))  * log(tau2Len_spf(s,p,f));
        }
        // Set tmp vectors to zero 
        fleetAgeObs.setZero();
        fleetAgePred.setZero();
        fleetLenObs.setZero();
        fleetLenPred.setZero();

      }
    }


  // Shared Hierarchical Prior Distributions //
  // loop over species and stocks, add penalty for stock deviations from
  // species mean
  vector<Type> steepnessnlp_s(nS);
  array<Type> steepnessnlp_sp(nS,nP);

  // Synoptic catchability prior
  vector<Type> qnlpSyn_s(nS);
  Type qnlpSyn = 0.;
  qnlpSyn_s.setZero();
  Type qnlpSurv = 0.;

  // Now a prior on catchability and observation error SD
  // Surveys only really act on one stock each, 
  // so their catchabiliity prior will be across species
  vector<Type> tauObsnlp_f(nF);
  tauObsnlp_f.setZero();
  for( int s = 0; s < nS; s++ )
    for( int p = 0; p < nP; p++ )
      // Penalise fleet specific catchabilities if
      // they are a calculated index
      for( int f = 0; f < nF; f++ )
      {
        // Synoptic survey
        if( group_f(f) == 2)
          if( calcIndex_spf(s,p,f) > 0)
            qnlpSyn_s(s) -= dnorm( log(q_spf(s,p,f)), lnqbarSyn_s(s), tauqSyn_s(s), true);
    
        // HS Assemblage survey
        if( group_f(f) == 1)
          if( calcIndex_spf(s,p,f) > 0)
            qnlpSurv -= dnorm( log(q_spf(s,p,f)), log(mqSurveys), sdqSurveys,true );
        
        // Observation error SD
        tauObsnlp_f(f) += (IGatau_f(f)+Type(1))*lntau_spf(s,p,f)+IGbtau_f(f)/tau_spf(s,p,f);
      }
    
  // Synoptic species qs are penaliesd against a complex q
  // if there are more than 1 species
  if(nS > 1)
  {
    qnlpSyn -= dnorm( lnqbarSyn_s, lnqbarSyn, tauqSyn, true ).sum();
    // Synoptic complex qs have a prior
    qnlpSurv -= dnorm( lnqbarSyn, log(mqSurveys), sdqSurveys, true);
  }
  // Synoptic species qs are penalised against
  // the overal q prior if there is only 1 species
  if(nS == 1)
    qnlpSyn -= dnorm( lnqbarSyn_s, log(mqSurveys), sdqSurveys, true ).sum();
  
  // time-varying catchability
  Type qnlp_tv = 0;
  qnlp_tv -= dnorm( epslnq_vec, 0, sigmalnq, true).sum();



  // Mortality prior
  vector<Type>  Mnlp_s(nS);
  vector<Type>  Mnlp_p(nP);
  array<Type>   Mnlp_sp(nS,nP);
  // 0-initialise
  steepnessnlp_s.setZero();
  steepnessnlp_sp.setZero();
  Mnlp_s.setZero();
  Mnlp_p.setZero();
  Mnlp_sp.setZero();
  Type sel_nlp = 0.;


  for( int s = 0; s < nS; s++ )
    for(int p = 0; p < nP; p++ )
    {
      // Steepness
      steepnessnlp_sp(s,p) -= dnorm( log(h_sp(s,p)), log(h_s(s)), sigmah_s(s), true);
      Mnlp_sp(s,p) -= dnorm( log(M_sp(s,p)), log(M_s(s)), sigmaM_s(s), true);
    }
  
  // add species level h prior
  steepnessnlp_s -= dnorm( log(h_s), log(mh), sigmah, true );
  Mnlp_s -= dnorm( log(M_s), log(muM), sigmaM, true );

  // Add time-varying selectivity deviations
  sel_nlp -= dnorm( epsxSel50_vec, Type(0), sigmaSel, true).sum();
  sel_nlp -= dnorm( epsxSelStep_vec, Type(0), sigmaSel, true).sum();
  // Prior on xSel95_sf
  for(int s = 0; s < nS; s++)
    for( int f = 0; f < nF; f++)
    {
      sel_nlp -= dnorm( lnxSel50_sf(s,f), pmlnxSel50_sf(s,f), cvxSel, true);
      sel_nlp -= dnorm( lnxSelStep_sf(s,f), pmlnxSelStep_sf(s,f), cvxSel, true);
    }


  // VonB priors
  vector<Type>  L2nlp_s(nS);
  vector<Type>  L1nlp_s(nS);
  vector<Type>  vonKnlp_s(nS);
  // Zero-init
  L2nlp_s.setZero();
  L1nlp_s.setZero();
  vonKnlp_s.setZero();

  // Compute
  L1nlp_s   -= dnorm( lnL1_s, pmlnL1_s, cvL1, true);
  L2nlp_s   -= dnorm( log(L2_s), pmlnL2_s, cvL2, true);
  vonKnlp_s -= dnorm( lnvonK_s, pmlnVonK, cvVonK, true);
        

  // Loop over populations, penalise pop-specific deviations from spec-specific values
  // vector<Type> vonBnlp_p(nP);
  // vector<Type> qnlp_p(nP);
  // vector<Type> steepnessnlp_p(nP);
  // vector<Type> Mnlp_p(nP);
  // vonBnlp_p.setZero();
  // qnlp_p.setZero();
  // steepnessnlp_p.setZero();
  // Mnlp_p.setZero();
  // for( int p = 0; p < nP; p ++ )
  // {
  //   int specID = s_p(p);
  //   // Growth
  //   // asymptotic length
  //   vonBnlp_p(p) += Type(0.5) * ( log(sigmaL2_s(specID)) + square(L2_p(p) - muL2_s(specID) ) / square(sigmaL2_s (specID)) );
  //   // growth coefficient
  //   vonBnlp_p(p) += Type(0.5) * ( log(sigmavonK_s(specID)) + square(vonK_p(p) - muvonK_s(specID) ) / square(sigmavonK_s (specID) ));
  //   // time at length 0
  //   vonBnlp_p(p) += Type(0.5) * ( log(sigmat0_s(specID)) + square(vont0_p(p) - mut0_s(specID) ) / square(sigmat0_s (specID) ));
  //   // Steepness
  //   Type sigmaSteep_spec = exp(lnsigmaSteep_s(specID));
  //   steepnessnlp_p(p) += Type(0.5) * ( lnsigmaSteep_s(specID) + square(h_p(p) - h_s(specID)) / square(sigmaSteep_spec) );
  //   // Natural mortality?
  //   Mnlp_p(p) += Type(0.5) * ( lnsigmaM_s(specID) + square(M_p(p) - M_s(specID)) / square(sigmaM_s(specID)) );

  //   // For selectivity and catchability, loop over fleets
  //   for( int f = 0; f < nF; f++ )
  //   {      
  //     // Catchability
  //     if( group_f(f) == 0) 
  //       qnlp_p(p) += Type(0.5) * ( lntauq_fs(f,specID) + square( lnq_fp(f,p) - lnqbar_fs(f,specID) ) / exp(2*lntauq_fs(f,specID)) );
  //   }
  // }
  // Type pop_nlp = 0.0;
  // pop_nlp += vonBnlp_p.sum() + steepnessnlp_p.sum() + Mnlp_p.sum() + qnlp_p.sum();

  // // Now loop over species and penalise by "complex level" means - will have to think
  // // carefully about what parameters have a complex mean, and what don't
  // // vector<Type> vonBnlp_p(nP);
  // // vector<Type> selnlp_p(nP);
  // vector<Type> qnlp_s(nS);
  // vector<Type> selnlp_s(nS);
  // vector<Type> steepnessnlp_s(nS);
  // vector<Type> Mnlp_s(nS);
  // // vonBnlp_p.setZero();
  // selnlp_s.setZero();
  // qnlp_s.setZero();
  // steepnessnlp_s.setZero();
  // Mnlp_s.setZero();
  // for( int s = 0; s < nS; s++ )
  // {
  //   // Steepness
  //   Type sigmaSteep = exp(lnsigmaSteep);
  //   steepnessnlp_s(s) += Type(0.5) * ( lnsigmaSteep + square(h_s(s) - muSteep) / square(sigmaSteep) );

  //   // Natural mortality
  //   Mnlp_s(s) += Type(0.5) * ( lnsigmaM + square(M_s(s) - muM) / square(sigmaM) );

  //   // Catchability
  //   for( int f = 0; f < nF; f++ )
  //   {
  //     // Catchability
  //     if( group_f(f) == 0) 
  //       qnlp_s(s) += Type(0.5) * ( lntauq_f(f) + square( lnqbar_fs(f,s) - lnqbar_f(f) ) / exp(2*lntauq_f(f)) );

  //     // Selectivity
  //     selnlp_s(s) += Type(0.5) * ( log(sigmaxSel50_f(f) ) + square( muxSel50_f(f) - xSel50_sf(s,f) ) / square(sigmaxSel50_f(f)) );
  //     selnlp_s(s) += Type(0.5) * ( log(sigmaxSel95_f(f) ) + square( muxSel95_f(f) - xSel95_sf(s,f) ) / square(sigmaxSel95_f(f)) );
  //   }
  // }
  // Type spec_nlp = 0.0;  
  // spec_nlp += steepnessnlp_s.sum() + Mnlp_s.sum() + qnlp_s.sum();
  

  // Now take the sum of all the likelihoods, keep separate
  // in case we want to weight data later
  // if( parSwitch == 1 )
  //   parallel_accumulator<Type> joint_nlp(this);
  // else
  Type joint_nlp = 0.0;

  array<Type> B0nlp_sp(nS,nP);
  B0nlp_sp.setZero();
  for( int p = 0; p < nP; p++)
    B0nlp_sp.col(p) += lambdaB0 * (B0_sp.col(p));

  
  // Observations
  joint_nlp += Ctnll_sp.sum();       // Catch
  joint_nlp += ageLikeWt * lenCompsnll_spf.sum(); // Length compositions
  joint_nlp += lenLikeWt * ageCompsnll_spf.sum(); // Age compositions
  joint_nlp += idxLikeWt * CPUEnll_spft.sum();     // Survey CPUE
  joint_nlp += steepnessnlp_s.sum() + steepnessnlp_sp.sum();
  joint_nlp += sel_nlp;
  joint_nlp += tauObsnlp_f.sum();
  // Growth model
  joint_nlp += growthLikeWt * vonBnll_spf.sum();     // Growth model
  joint_nlp += growthLikeWt * (L1nlp_s.sum() + L2nlp_s.sum() + vonKnlp_s.sum());
  // Recruitment errors
  joint_nlp += recnll_sp.sum();      // recruitment process errors
  // Commercial q
  joint_nlp += qnlpSurv + qnlpSyn_s.sum() + qnlpSyn;
  joint_nlp += Mnlp_s.sum() + Mnlp_s.sum();
  joint_nlp += B0nlp_sp.sum();
  // joint_nlp += pop_nlp + spec_nlp;

  
  // // Return quantities
  // Model dimensions
  REPORT( nS );
  REPORT( nP );
  REPORT( nF );
  REPORT( nT );
  REPORT( nL );
  REPORT( nA );

  // Leading parameters
  REPORT( B0_sp );    
  REPORT( h_sp );     
  REPORT( M_sp );     
  REPORT( L2_sp );  
  REPORT( vonK_sp );  
  REPORT( L1_sp ); 
  REPORT( sigmaLa_s );
  REPORT( sigmaLb_s );
  REPORT( sigmaR_sp );
  REPORT( tauC_f );
  REPORT( tauD_f );
  // Fishery/Survey model pars
  REPORT( q_spf );           // Catchability
  REPORT( q_spft );
  REPORT( tau_spf );         // Observation error
  REPORT( xSel50_spft );   // length at 50% sel by species/fleet/pop
  REPORT( xSel95_spft );   // length at 95% sel by species/fleet/pop
  REPORT( xSelStep_spft ); // length at 95% sel by species/fleet/pop
  REPORT( xSel50_sf );     // length at 50% sel by species/fleet
  REPORT( xSel95_sf );     // length at 95% sel by species/fleet
  REPORT( xSelStep_sf );   // length at 95% sel by species/fleet
  REPORT( F_spft );          // Fishing mortality
  REPORT( sel_lfspt );
  REPORT( sel_afspt );
  // Model states
  REPORT( B_aspt );
  REPORT( N_aspt );
  REPORT( R_spt );
  REPORT( F_aspft );
  REPORT( Z_aspt );
  REPORT( C_aspft );
  REPORT( Cw_aspft );
  REPORT( predCw_spft );
  REPORT( predC_spft );
  REPORT( Bv_spft );
  REPORT( Nv_spft );
  REPORT( Nv_aspft );
  REPORT( SB_spt );
  REPORT( B_spt );
  // Stock recruit parameters
  REPORT( Surv_asp );
  REPORT( SSBpr_asp );
  REPORT( R0_sp );
  REPORT( phi_sp );
  REPORT( reca_sp );
  REPORT( recb_sp );

  // Growth and maturity
  REPORT( probLenAge_lasp );
  REPORT( probAgeLen_alspft );
  REPORT( lenAge_asp );
  REPORT( Wlen_ls );
  REPORT( matLen_ls );
  REPORT( matAge_asp );
  REPORT( xMat50_s );
  REPORT( xMat95_s );
  REPORT( meanWtAge_asp );
  REPORT( ageAtLenResids_alspft );
  REPORT( nObsAgeAtLen_spf );
  REPORT( etaSumSqAgeAtLen_spf );
  REPORT( nResidsAgeAtLen_spf );
  REPORT( L1_s );
  REPORT( L2_s );
  REPORT( vonK_s );
  REPORT( A1_s );
  REPORT( A2_s );
  REPORT( deltaVonK_sp );
  REPORT( deltaL2_sp );

  // Species and complex level hyperparameters
  REPORT( h_s );
  REPORT( sigmah_s );
  REPORT( mh );
  REPORT( sigmah );
  REPORT( qbarSyn );
  REPORT( tauqSyn );
  REPORT( qbarSyn_s );
  REPORT( tauqSyn_s );

  // Echo input priors
  REPORT( pmlnxSel50_sf );
  REPORT( pmlnxSelStep_sf );
  REPORT( pmlnL1_s );
  REPORT( pmlnL1_s );
  REPORT( pmlnVonK );
  REPORT( muxSel50_sg );
  REPORT( muxSel95_sg );
  REPORT( sigmaxSel50_sg );
  REPORT( sigmaxSel95_sg );

  // Random effects
  REPORT( epsxSel50_spft );
  REPORT( epsxSelStep_spft );
  REPORT( epslnq_spft );
  REPORT( omegaR_spt );
  REPORT( omegaRinit_spa );

  // Species/stock effects
  REPORT( epsM_sp );
  REPORT( epsM_s );
  REPORT( epsSteep_sp );
  REPORT( epsSteep_s );


  // Likelihood values
  REPORT( joint_nlp );
  REPORT( recnll_sp );
  REPORT( vonBnll_spf );
  REPORT( CPUEnll_spft );
  REPORT( ageCompsnll_spf );
  REPORT( etaSumSqAge_spf );
  REPORT( nResidsAge_spf );
  REPORT( nObsAge_spf );
  REPORT( tau2Age_spf );
  REPORT( ageRes_aspft );
  REPORT( lenCompsnll_spf );
  REPORT( etaSumSqLen_spf );
  REPORT( nResidsLen_spf );
  REPORT( nObsLen_spf );
  REPORT( tau2Len_spf );
  REPORT( lenRes_lspft );
  REPORT( steepnessnlp_sp );
  REPORT( steepnessnlp_s );
  REPORT( Ctnll_sp );
  REPORT( qnlpSurv );
  REPORT( qnlpSyn );
  REPORT( qnlpSyn_s );
  REPORT( tauObsnlp_f );
  REPORT( minAgeProp );
  REPORT( minLenProp );
  REPORT( L1nlp_s );
  REPORT( L2nlp_s );
  REPORT( vonKnlp_s );
  REPORT( sel_nlp );

  // Data
  REPORT( I_spft );
  REPORT( C_spft );
  REPORT( D_spft );
  REPORT( ALK_spalft );
  REPORT( age_aspft );
  REPORT( len_lspft );
  REPORT( group_f );
  REPORT( A_s );
  REPORT( minA_s );
  REPORT( L_s );
  REPORT( minL_s );
  REPORT( lenD_s );
  REPORT( calcIndex_spf );

  // Predicted observations
  REPORT( lDist_lspft_hat );
  REPORT( aDist_aspft_hat );
  REPORT( I_spft_hat );

  // REPORT( matX );
  // REPORT( selX );


  // sd report for standard errors
  ADREPORT(lnB0_sp);    
  // ADREPORT(logitSteep_sp);     
  // ADREPORT(lnM_sp);     
  // ADREPORT(lnRinit_sp); 
  // ADREPORT(lnL2_sp);  
  // ADREPORT(lnVonK_sp);  
  // ADREPORT(lnL1_sp); 
  // ADREPORT(lnsigmaL_s);
  // ADREPORT(lnsigmaR_sp);
  // ADREPORT(lnq_spf)          // Catchability
  // ADREPORT(lntau_spf)        // Observation error
  // ADREPORT(lnxSel50_sf)   // length at 50% sel
  // ADREPORT(lnxSelStep_sf)   // length at 95% sel
  // ADREPORT(lnF_spft)         // Fishing mortality


  return( joint_nlp );

}

