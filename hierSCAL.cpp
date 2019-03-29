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
//  6. Un-center hierarchical distributions
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

// square()
// Shortcut wrapper for squaring quantities
// inputs:    x = quantity to be squared
// ouputs:    y = x^2
// Usage:     when squaring.
// Source:    Almost surely stolen from somewhere.
template <class Type> 
Type square(Type x){return pow(x,2);}
// VECTORIZE1_t(square)

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

// solveBaranov_spfx()
// Newton-Rhapson solver for Baranov catch equation for a multi-stock
// multi-species, age and sex structured population fished by multiple
// fleets
// inputs:    nIter = number of NR iterations
//            Bstep = fraction of NR step (Jacobian) to take at each iteration
//            C_spf = Catch
//            M_spx = natural mortality
//            B_aspx = Biomass at age/sex
//            B_spx = biomass for each sex
//            sel_aspfx = selectivity for each age/sex in each fleet
//            & Z = total mortality (external variable)
//            & F = Fishing mortality (external variable)
// returns:   NA, void function
// Side-effs: variables passed as Z, F overwritten with total, fishing mortality
// Author:    Modified by S. D. N. Johnson from S. Rossi and S. P. Cox
template<class Type>
array<Type> solveBaranov_spf( int   nIter,
                              Type  Bstep,
                              vector<int>  A_s,         // number of age classes by species
                              array<Type>  C_spf,       // Total observed catch
                              array<Type>  M_spx,       // Mortality rate
                              array<Type>  B_aspx,      // Biomass at age/sex
                              array<Type>  vB_aspfx,    // vuln biomass at age
                              array<Type>  vB_spfx,     // vuln biomass for each sex
                              array<Type>  vB_spf,      // vuln biomass in each fleet
                              array<Type>  sel_aspfx,   // selectivity at age/sex
                              array<Type>& Z_aspx,      // total mortality at age/se
                              array<Type>& F_spf )       // fleet F
{
  int nS = C_spf.dim(0);
  int nP = C_spf.dim(1);
  int nF = C_spf.dim(2);
  int nA = B_aspx.dim(0);
  int nX = B_aspx.dim(3);

  array<Type> f_spf(nS,nP,nF);              // Function value
  array<Type> J_spf(nS,nP,nF);              // Jacobian
  array<Type> newZ_aspx(nA,nS,nP,nX);       // Updated Z
  array<Type> tmpZ_aspx(nA,nS,nP,nX);       // Updated Z
  array<Type> tmp_spf(nS,nP,nF);            // predicted catch given F
  array<Type> tmp_axspf(nA,nX,nS,nP,nF);    // predicted catch-at-age/sex given F
  array<Type> F_aspfx(nA,nS,nP,nF,nX);      // Fishing mortality at age/sex

  F_aspfx.setZero();
  newZ_aspx.setZero();
  tmpZ_aspx.setZero();
  f_spf.setZero();
  J_spf.setZero();
  tmp_spf.setZero();
  tmp_axspf.setZero();


  
  // Initial approximation of F
  for( int s = 0; s < nS; s++ )
    for( int p = 0; p < nP; p++)
    {
      for( int x = 0; x < nX; x++)
        newZ_aspx.col(x).col(p).col(s).fill(M_spx(s,p,x));

      for( int f = 0; f < nF; f++ )
      {
        // Use catch plus bio for years where bio is dangerously small
        F_spf(s,p,f) = C_spf(s,p,f) / (vB_spf(s,p,f));
        for( int x = 0; x < nX; x++)
          newZ_aspx.col(x).col(p).col(s) +=  F_spf(s,p,f) * sel_aspfx.col(x).col(f).col(p).col(s);
      }
    }
  

  // Refine F
  for( int i=0; i<nIter; i++ )
  {
    // Now reset newZ
    tmp_spf.setZero();
    J_spf.setZero();
    f_spf.setZero();
    tmpZ_aspx.setZero();


    // Reset objective function and Z
    f_spf += C_spf;


    // Calculate predicted catch and Jacobian
    for( int s = 0; s < nS; s++ )
      for( int p = 0; p < nP; p++ )
      {
        for( int x = 0; x < nX; x++ )
        {
          // Fill newZ up with M
          for( int a = 0; a < A_s(s); a++ )
          {
            // Overwrite newZ into Z array
            // tmpZ_aspx(a,s,p,x)     += newZ_aspx(a,s,p,x);

            for(int f = 0; f < nF; f++ )  
            {             
              tmp_spf(s,p,f) += B_aspx(a,s,p,x) * (1 - exp(-newZ_aspx(a,s,p,x)) ) * sel_aspfx(a,s,p,f,x) *  F_spf(s,p,f) / newZ_aspx(a,s,p,x);

              Type tmpJ1  = vB_aspfx(a,s,p,f,x) / newZ_aspx(a,s,p,x);
              Type tmpJ2  = F_spf(s,p,f) * sel_aspfx(a,s,p,f,x) * exp( - newZ_aspx(a,s,p,x));
              Type tmpJ3  = (1 - exp( -newZ_aspx(a,s,p,x)));
              Type tmpJ4  = (newZ_aspx(a,s,p,x) - F_spf(s,p,f) * sel_aspfx(a,s,p,f,x))/newZ_aspx(a,s,p,x);

              J_spf(s,p,f) -= tmpJ1 * ( tmpJ2 + tmpJ3 * tmpJ4);

            }
          }
        }
      }

    // Subtract predicted catch
    f_spf -= tmp_spf;

    // Updated fishing mortality
    F_spf -= Bstep * f_spf / J_spf;

    newZ_aspx.setZero();

    // Updated total mortality
    for( int s = 0; s < nS; s++ )
      for( int p = 0; p < nP; p++ )
        for( int x = 0; x < nX; x++ )
        {
          newZ_aspx.col(x).col(p).col(s).fill(M_spx(s,p,x));
          for( int f= 0 ; f < nF; f ++)
            newZ_aspx.col(x).col(p).col(s) += sel_aspfx.col(x).col(f).col(p).col(s) * F_spf(s,p,f) ;  
        }

  }  // end i

  // Now save the tmpZ_aspx array out to the Z_aspx input
  Z_aspx = newZ_aspx;

  return( J_spf );

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
  DATA_ARRAY(ALK_spalftx);              // Age-length observations (in freq) by pop, separated by fleet and year
  DATA_ARRAY(age_aspftx);               // Age observations for each population and fleet over time (-1 missing)
  DATA_ARRAY(len_lspftx);               // Length observations for each population and fleet over time (-1 missing)
  DATA_IVECTOR(group_f);                // Fleet group (0=comm, 1=HSAss,2=Synoptic)
  DATA_IVECTOR(A_s);                    // +Group age by species
  DATA_IVECTOR(L_s);                    // +Group length by species
  DATA_IVECTOR(minA_s);                 // Min age class by species
  DATA_IVECTOR(minL_s);                 // Min length bin (cm) by species
  DATA_IVECTOR(lenD_s);                 // Discard length by species

  // Model dimensions
  int nS = I_spft.dim(0);               // No. of species
  int nP = I_spft.dim(1);               // No. of stocks in each species
  int nF = I_spft.dim(2);               // No. of fleets (surveys + trawl)  
  int nT = I_spft.dim(3);               // No of time steps
  int nL = len_lspftx.dim(0);           // max no of length bins (for creating state arrays)
  int nA = age_aspftx.dim(0);           // max no of age bins (for creating state arrays)
  int nX = age_aspftx.dim(5);           // No of sex classes (0 == male, 1 == female)

  // Model switches - fill these in when we know what we want
  DATA_IVECTOR(swRinit_s);              // species fished initialisation switch (0 == unfished, 1 == fished)
  DATA_INTEGER(parSwitch);              // parallel accumulator switch for faster run time
  DATA_ARRAY(calcIndex_spf);            // Switches to calculate indices
  DATA_IVECTOR(tvSelFleets);            // Switch for fleet to have time varying selectivity (0 = off, 1 = on)
  DATA_IVECTOR(tvqFleets);              // Switch for fleet to have time varying catchability (0 = off, 1 = on)
  DATA_IVECTOR(regFfleets);             // Switch for each fleet to have Fs regularised (Fbar penalised)
  DATA_SCALAR(idxLikeWt);               // Scalar modifying the weight of the age likelihood
  DATA_SCALAR(ageLikeWt);               // Scalar modifying the weight of the age likelihood
  DATA_SCALAR(lenLikeWt);               // Scalar modifying the weight of the length likelihood
  DATA_SCALAR(growthLikeWt);            // Scalar modifying the weight of the growth model likelihood
  DATA_IVECTOR(tFirstRecDev_s);         // Initial rec devs for each species
  DATA_IVECTOR(tLastRecDev_s);          // Last recruitment devs for each species
  DATA_SCALAR(minAgeProp);              // Minimum observed proportion in age comps
  DATA_SCALAR(minLenProp);              // Minimum observed proportion in length comps
  DATA_SCALAR(minPAAL);                 // Minimum observed proportion age at length
  DATA_STRING(matX);                    // Age/Len switch for maturity model
  DATA_STRING(selX);                    // Age/Len switch for selectivity model
  DATA_STRING(lenComps);                // Switch for computation of expected length comps
  DATA_SCALAR(lambdaB0);                // lambda scalar for exponential prior on B0
  DATA_ARRAY(minTimeIdx_spf);           // array of earliest times for each index
  DATA_INTEGER(nBaranovIter);           // number of baranov steps
  DATA_SCALAR(lambdaBaranovStep);       // fractional step size for Newton-Rhapson Baranov iteration
  // DATA_STRING(recruitVariance);         // Character-switch for recruitment variance model

  // Fixed values
  DATA_IVECTOR(A1_s);                   // Age at which L1_s is estimated
  DATA_IVECTOR(A2_s);                   // Age at which L2_s is estimated


  /* -------------------------- Parameter Section ------------------------ */
  // Leading Parameters
  // Biological
  PARAMETER_ARRAY(lnB0_sp);             // Biomass at unfished for ind. stocks
  PARAMETER(logitSteep);                // pop specific steepness
  PARAMETER(lnM);                       // stock-specific M
  PARAMETER_VECTOR(lnL2step_s);         // species-spec. step in Schnute-vonB L2 parameter
  PARAMETER_VECTOR(lnvonK_s);           // species-specific vonB K parameter
  PARAMETER_VECTOR(lnL1_s);             // species-specific vonB Length-at-age 1
  PARAMETER_ARRAY(deltaL2_sp);          // species-pop specific deviation in L2
  PARAMETER_VECTOR(lnsigmaL2_s);        // species level L2 dev SD
  PARAMETER_ARRAY(deltaVonK_sp);        // species-pop specific deviation in VonK parameter
  PARAMETER_VECTOR(lnsigmavonK_s);      // species level vonK dev SD
  PARAMETER_ARRAY(deltaL2_sx);          // sex-species-pop specific deviation in L2
  PARAMETER_ARRAY(deltaVonK_sx);        // sex-species-pop specific deviation in VonK parameter
  PARAMETER(lnsigmaL2);                 // Complex level L2 dev SD (Is this necessary??)
  PARAMETER(lnsigmavonK);               // Complex level vonK dev SD (Is this necessary??)
  PARAMETER_VECTOR(lnsigmaLa_s);        // species-specific individual growth SD intercept
  PARAMETER_VECTOR(sigmaLb_s);          // species-specific individual growth SD slope
  PARAMETER_VECTOR(LWa_s);              // species L-W conversion a par (units conv: cm -> kt)
  PARAMETER_VECTOR(LWb_s);              // species L-W conversion b par (length to volume)
  PARAMETER_VECTOR(xMat50_s);           // x (age/len) at 50% maturity
  PARAMETER_VECTOR(xMat95_s);           // x (age/len) at 95% maturity
  // Observation model
  // survey //
  PARAMETER_ARRAY(lnq_spf);             // fleet-species-stock specific catchability
  PARAMETER_ARRAY(lntauObs_spf);        // fleet-species-stock specific observation error variance
  // Fishery model
  // selectivity by fleet/species //
  PARAMETER_ARRAY(lnxSel50_sf);         // Selectivity Alpha parameter by species/fleet
  PARAMETER_ARRAY(lnxSelStep_sf);       // Selectivity Beta parameter by species/fleet
  PARAMETER_ARRAY(epsxSel50_spf);       // Selectivivity devs for stocks
  PARAMETER_ARRAY(epsxSelStep_spf);     // Selectivivity devs for stocks

  // Fishing mortality
  PARAMETER_VECTOR(lntauC_f);           // Catch observation SD by fleet
  PARAMETER_VECTOR(lntauD_f);           // Discards observation SD by fleet

  // Priors
  // selectivity //
  PARAMETER_ARRAY(muxSel50_sg);         // mean length at 50% selectivity by fleet group
  PARAMETER_ARRAY(muxSel95_sg);         // mean selectivity-at-length step by fleet group 
  PARAMETER_ARRAY(sigmaxSel50_sg);      // SD in length at 50% selectivity by fleet group
  PARAMETER_ARRAY(sigmaxSel95_sg);      // SD in length at 95% selectivity by fleet group

  // Count fleet groups
  int nGroups = muxSel95_sg.dim(1);

  // catchability //
  PARAMETER_VECTOR(lnqbarSyn_s);        // Synoptic survey q by species (each leg surveys one stock)
  PARAMETER_VECTOR(lntauqSyn_s);        // Synoptic survey q SD by species (each leg survey one stock)
  PARAMETER(lnqbarSyn);                 // Synoptic survey complex q
  PARAMETER(lntauqSyn);                 // Synoptic survey complex q SD
  PARAMETER(mqSurveys);                 // prior mean survey q (HSAss included)
  PARAMETER(sdqSurveys);                // Prior SD survey q
  // steepness //
  PARAMETER_VECTOR(lnsigmah_s);         // species level steepness SD
  PARAMETER(logit_muSteep);             // complex steepness
  PARAMETER(lnsigmah);                  // complex level steepness SD
  // Natural Mortality //
  PARAMETER_VECTOR(lnsigmaM_s);         // species M SD
  PARAMETER(ln_muM);                    // Multispecies assemblage mean M
  PARAMETER(lnsigmaM);                  // Assemblage M SD
  // Index observation errors SD
  PARAMETER_VECTOR(IGatau_f);           // IG a parameter for tauObs_spf observation errors
  PARAMETER_VECTOR(IGbtau_f);           // IG b parameter for tauObs_spf observation errors


  // Random Effects
  // recruitment //
  PARAMETER_VECTOR(epsSteep_s);         // Species level steepness effect
  PARAMETER_VECTOR(epsM_s);             // Species level natural mortality effect
  PARAMETER_ARRAY(epsSteep_sp);         // stock level steepness effect
  PARAMETER_ARRAY(epsM_sp);             // stock level natural mortality effect
  PARAMETER_ARRAY(epsM_sx);             // sex specific natural mortality effect
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
  // PARAMETER_VECTOR(recCorr_vec);        // vector for cholesky decomposition of recruitment correlation mtx


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
  PARAMETER_ARRAY(epsilonM_spt);        // M deviations by population, 2:nT
  PARAMETER_VECTOR(lnsigmaM_sp);        // M deviation pop-specific SD
  PARAMETER_VECTOR(logitMgamma_sp);     // AR1 auto-corr on M devs effect (omega)
  */

  /* ------------------------------ Procedure Section --------------------- */
  // Exponentiate leading parameters
  // Biological //
  array<Type>  B0_sp(nS,nP);
  array<Type>  h_sp(nS,nP);
  array<Type>  M_sp(nS,nP);
  array<Type>  L2_sp(nS,nP);
  array<Type>  vonK_sp(nS,nP);
  array<Type>  L1_sp(nS,nP);
  array<Type>  sigmaR_sp(nS,nP);
  // Sex specific
  array<Type>  M_spx(nS,nP,nX);
  array<Type>  L2_spx(nS,nP,nX);
  array<Type>  L1_spx(nS,nP,nX);
  array<Type>  vonK_spx(nS,nP,nX);

  // Growth model vectors
  vector<Type>  sigmaLa_s(nS);
  vector<Type>  L1_s(nS);
  vector<Type>  L2_s(nS);
  vector<Type>  vonK_s(nS);

  // Growth model stock priors
  vector<Type>  deltaVonKbar_p(nP);
  vector<Type>  deltaL2bar_p(nP);


  // derived variables
  // Stock recruitment //
  array<Type>  R0_sp(nS,nP);            // eqbm recruitment
  array<Type>  phi_sp(nS,nP);           // eqbm SSB per recruit
  array<Type>  reca_sp(nS,nP);          // BH a parameter for pops
  array<Type>  recb_sp(nS,nP);          // BH b parameter for pops

  // Growth //
  array<Type>   Wlen_ls(nL,nS);                 // weight-at-length by species
  array<Type>   lenAge_aspx(nA,nS,nP,nX);           // Mean length-at-age by population
  array<Type>   probLenAge_laspx(nL,nA,nS,nP,nX);
  array<Type>   meanWtAge_aspx(nA,nS,nP,nX);
  array<Type>   probAgeLen_alspftx(nA,nL,nS,nP,nF,nT,nX);

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

  // parallel_accumulator<Type> f(this);
  Type f = 0;
  Type joint_nlp = 0.0;

  // Let's test the multi-dimensional containers thingy
  // I think it will speed things up later
  vector<int> d(3); d << 2,2,2;
  int a = d(0)*d(1)*d(2);

  matrix<Type> m(2,3); 
  m << 1,2,3,4,5,6;

  vector< matrix<Type> > V(a);
  array< matrix<Type> > A(V,d);

  for(int i=0; i<2; i++)
    for(int j=0; j<2; j++)
      for(int k=0; k<2; k++)
        A(i,j,k) = m*i*j*k;



  // Prior Hyperparameters //
  // Steepness
  vector<Type>  sigmah_s      = exp(lnsigmah_s);
  Type          sigmah        = exp(lnsigmah);
  vector<Type>  h_s           = .2 + Type(0.78) / ( Type(1.0) + exp( -1 * (logitSteep + sigmah * epsSteep_s ) ) );
  Type          mh            = .2 + Type(0.78) / ( Type(1.0) + exp( -1 * logit_muSteep) );
  // Need to add an sdh

  // Natural mortality
  Type          sigmaM        = exp(lnsigmaM);
  vector<Type>  sigmaM_s      = exp(lnsigmaM_s);
  vector<Type>  M_s           = exp(lnM + sigmaM * epsM_s);
  Type          muM           = exp(ln_muM);
  // need to add an sdM
  
  // RE SDs
  Type          sigmaSel      = exp(lnsigmaSel);
  Type          sigmalnq      = exp(lnsigmaepslnq);
  // Catchability for synoptic survey
  vector<Type>  qbarSyn_s     = exp(lnqbarSyn_s);
  vector<Type>  tauqSyn_s     = exp(lnqbarSyn_s);
  Type          qbarSyn       = exp(lnqbarSyn);
  Type          tauqSyn       = exp(lntauqSyn);


  // Growth model hierarchical priors
  L1_s      = exp(lnL1_s); 
  vonK_s    = exp(lnvonK_s);
  L2_s      = exp(lnL1_s) + exp(lnL2step_s);
  sigmaLa_s = exp(lnsigmaLa_s);

  // Fix sigmaVonK and sigmaL2 for now
  Type sigmavonK  = exp(lnsigmavonK);
  Type sigmaL2    = exp(lnsigmaL2);

  vector<Type>    sigmavonK_s = exp(lnsigmavonK_s);
  vector<Type>    sigmaL2_s   = exp(lnsigmaL2_s);

  // Fill arrays that hold stock-specific parameters
  for( int pIdx = 0; pIdx < nP; pIdx++ )
  {
    B0_sp.col(pIdx)      = exp(lnB0_sp.col(pIdx));
    h_sp.col(pIdx)       = .2 + 0.8 / ( Type(1.0) + exp( -1 * (logitSteep + sigmah * epsSteep_s + sigmah_s * epsSteep_sp.col(pIdx) ) ) );
    M_sp.col(pIdx)       = exp(lnM) * exp(sigmaM * epsM_s) * exp(sigmaM_s * epsM_sp.col(pIdx) );
    L1_sp.col(pIdx)      = L1_s;
    L2_sp.col(pIdx)      = L2_s * exp(sigmaL2_s * deltaL2_sp.col(pIdx));
    vonK_sp.col(pIdx)    = vonK_s * exp(sigmavonK_s * deltaVonK_sp.col(pIdx));
    sigmaR_sp.col(pIdx)  = exp(lnsigmaR_sp.col(pIdx));

    // And sex specific
    for( int x = 0; x < nX; x++)
    {
      L1_spx.col(x).col(pIdx)   = L1_s;                    
      M_spx.col(x).col(pIdx)    = M_sp.col(pIdx) * exp( sigmaM_s * epsM_sx.col(x));
      L2_spx.col(x).col(pIdx)   = L2_sp.col(pIdx) * exp( sigmaL2_s * deltaL2_sx.col(x));
      vonK_spx.col(x).col(pIdx) = vonK_sp.col(pIdx) * exp(sigmavonK_s * deltaVonK_sx.col(x));
    }
  }


  // Observation model
  array<Type>   q_spf(nS,nP,nF);
  array<Type>   tauObs_spf(nS,nP,nF);
  array<Type>   tau2Obs_spf(nS,nP,nF);
  q_spf.setZero();
  tauObs_spf.setZero();
  tau2Obs_spf.setZero();

  // Time-varying catchability
  array<Type>   q_spft(nS,nP,nF,nT);
  // Fishing mortality
  array<Type>   F_spft(nS,nP,nF,nT);
  array<Type>   Fbar_spf(nS,nP,nF);

  // selectivity //
  // Fleet/species
  array<Type>   xSel50_sf(nS,nF);
  array<Type>   xSelStep_sf(nS,nF);
  array<Type>   xSel95_sf(nS,nF);
  // Fleet/species/Stock
  array<Type>   xSel50_spf(nS,nP,nF);
  array<Type>   xSelStep_spf(nS,nP,nF);
  array<Type>   xSel95_spf(nS,nP,nF);
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
  Fbar_spf.setZero();

  // Exponential species/fleet sel
  xSel50_sf = exp(lnxSel50_sf);
  xSelStep_sf = exp(lnxSelStep_sf);
  xSel95_sf = xSel50_sf + xSelStep_sf;

  // Loop over gears, species, stocks and time steps, create a matrix
  // of deviations in selectivity pars and calculate the pars
  // for each time step
  int epsSelVecIdx = 0;
  int epslnqVecIdx = 0;
  int vecIdx = 0;
  for(int f = 0; f < nF; f++ )
  {
    for( int p = 0; p < nP; p++ )
    {
      xSel50_spf.col(f).col(p)    = exp(lnxSel50_sf.col(f) + epsxSel50_spf.col(f).col(p));
      xSelStep_spf.col(f).col(p)  = exp(lnxSelStep_sf.col(f) + epsxSelStep_spf.col(f).col(p));
      xSel95_spf.col(f).col(p)    = xSel50_spf.col(f).col(p) + xSelStep_spf.col(f).col(p);
      for( int s = 0; s < nS; s++)
      {
        int nPosCatch = 0;
        // Exponentiate q and tau for the index observations
        q_spf(s,p,f)        = exp(lnq_spf(s,p,f));
        tauObs_spf(s,p,f)   = exp(lntauObs_spf(s,p,f));
        // Compute obs err variance
        tau2Obs_spf(s,p,f)  = square(tauObs_spf(s,p,f));

        xSel50_spft(s,p,f,0)      = xSel50_spf(s,p,f);
        xSelStep_spft(s,p,f,0)    = xSelStep_spf(s,p,f);
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
          // Check if ages or lengths are observed this time step
          if( tvSelFleets(f) == 1)
            for( int x = 0; x < nX; x++)
            {
              if( ((age_aspftx(0,s,p,f,t,x) >= 0) & (ageLikeWt > 0)) |
                  ((len_lspftx(0,s,p,f,t,x) >= 0) & (lenLikeWt > 0)) )
              {
                epsxSel50_spft(s,p,f,t)   += epsxSel50_vec(epsSelVecIdx);
                epsxSelStep_spft(s,p,f,t) += epsxSelStep_vec(epsSelVecIdx);
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
          
        }
      }
    }
  }
  // Compute xSel95 for the time-varying arrays
  xSel95_spft = xSel50_spft + xSelStep_spft;


  // recruitment deviations and initRecDevs - convert from vector
  // to array of deviations
  // Also want to include a matrix of recruitment deviations
  // for estimating correlation
  matrix<Type> omegaRmat_spt(nS*nP,nT);
  omegaRmat_spt.setZero();

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
        omegaR_spt(s,p,t) = -5. + 10. / (1. + exp(-omegaR_vec(devVecIdx)));
        // And fill the matrix of deviations
        omegaRmat_spt( s*nP + p, t ) = -5. + 10. / (1. + exp(-omegaR_vec(devVecIdx)));
        devVecIdx++;
      }
      if(swRinit_s(s) == 1)
        for( int a = 0; a < A_s(s); a++ )
        {
          omegaRinit_asp(a,s,p) = -5. + 10./( 1 + exp(-omegaRinit_vec(initVecIdx) ));
          initVecIdx++;
        }
    }

  
  // Probably concentrate these out later....
  vector<Type>  tauC_f        = exp(lntauC_f);    
  vector<Type>  tauD_f        = exp(lntauD_f);    
  
  // Set up model state arrays
  array<Type> B_asptx(nA,nS,nP,nT,nX);      // Biomass at age-pop-time
  array<Type> N_asptx(nA,nS,nP,nT,nX);      // Numbers at age-pop-time
  array<Type> R_spt(nS,nP,nT);              // Recruits by pop-time
  array<Type> F_aspftx(nA,nS,nP,nF,nT,nX);  // Fishing mortality by age, fleet, population and time
  array<Type> Z_aspxt(nA,nS,nP,nX,nT);      // Total mortality by age, population and time
  array<Type> C_aspftx(nA,nS,nP,nF,nT,nX);  // Predicted catch-at-age (in numbers), fleet, population and time
  array<Type> Cw_aspftx(nA,nS,nP,nF,nT,nX); // Predicted catch-at-age (in weight), population and time
  array<Type> predCw_spft(nS,nP,nF,nT);     // Predicted catch (in weight) by population and time
  array<Type> predC_spft(nS,nP,nF,nT);      // Predicted catch (in numbers) by population and time
  array<Type> B_spt(nS,nP,nT);              // Total biomass by species, pop, time
  array<Type> SB_spt(nS,nP,nT);             // Spawning biomass by species, pop, time
  array<Type> Bv_spft(nS,nP,nF,nT);         // Vulnerable biomass by species, pop, time
  array<Type> Nv_spft(nS,nP,nF,nT);         // Vulnerable numbers by species, pop, time
  array<Type> Nv_aspftx(nA,nS,nP,nF,nT,nX); // Vulnerable numbers at age by species, pop, time
  array<Type> Surv_aspx(nA,nS,nP,nX);       // Eqbm survival at age
  array<Type> SSBpr_asp(nA,nS,nP);          // Spawning stock biomass per rec. at age



  // ---------------- Growth, Maturity, Stock-Recruitment ----------- //
  Wlen_ls.fill(-1.0);
  matLen_ls.fill(-1.0);
  for( int s = 0; s < nS; s++ )
  {
    // weight-at-length for each length bin midpoint
    Wlen_ls.col(s) = LWa_s(s) * pow( len, LWb_s(s) );
    
    // proportion mature-at-length
    if(matX == "length")
      matLen_ls.col(s) = 1/(1 + exp( -1. * log(Type(19.0)) * ( len - xMat50_s(s)) / (xMat95_s(s) - xMat50_s(s)) ) );  
  }

  // Now produce the probability of length-at-age matrix for each stock
  // (taken from LIME by Rudd and Thorson, 2017)
  // In the same loop, we're going to compute the ssbpr parameter phi
  // for each stock
  lenAge_aspx.setZero();
  probLenAge_laspx.setZero();
  meanWtAge_aspx.setZero();
  phi_sp.setZero();
  Surv_aspx.setZero();
  SSBpr_asp.setZero();
  matAge_asp.setZero();

  for( int x = 0; x < nX; x++)
    for( int s = 0; s < nS; s++)
    {
      int A1 = A1_s(s);
      int A2 = A2_s(s);
      for( int p = 0; p < nP; p++ )
      {
        Type tmpPhi = 0;

        lenAge_aspx.col(x).col(p).col(s) += L2_spx(s,p,x) - L1_sp(s,p);
        lenAge_aspx.col(x).col(p).col(s) *= (exp(-vonK_spx(s,p,x) * A1 ) - exp(-vonK_spx(s,p,x) * age) );
        lenAge_aspx.col(x).col(p).col(s) /= (exp(-vonK_spx(s,p,x) * A1 ) - exp(-vonK_spx(s,p,x) * A2) );
        lenAge_aspx.col(x).col(p).col(s) += L1_sp(s,p);

        for( int a = 0; a < A_s(s); a ++ )
        {
          Type  sumProbs  = 0;
          // Type vonK_tmp = vonK_spx(s,p,x);
          // lenAge_aspx(a,s,p,x) += L2_spx(s,p,x) - L1_sp(s,p);
          // lenAge_aspx(a,s,p,x) *= (exp(-vonK_tmp * A1 ) - exp(-vonK_tmp * (a+1)) );
          // lenAge_aspx(a,s,p,x) /= (exp(-vonK_tmp * A1 ) - exp(-vonK_tmp * A2) );
          // lenAge_aspx(a,s,p,x) += L1_sp(s,p);
          Type sigmaL = sigmaLa_s(s) + sigmaLb_s(s) * lenAge_aspx(a,s,p,x);
          for( int l = 0; l < L_s(s); l++ )
          {
            // Get length bin ranges
            Type lenHi = len(l) + .5;
            Type lenLo = len(l) - .5;
            // Compute density in each bin
            // Compute LH tail of the distribution
            
            if( (l >= 0) & (l < L_s(s) - 1 ) )
            {
              probLenAge_laspx(l,a,s,p,x) = pnorm( log(lenHi), log(lenAge_aspx(a,s,p,x)), sigmaL );
              // Subtract LH tail for interior bins
              if(l > 0)
                probLenAge_laspx(l,a,s,p,x) -= pnorm( log(lenLo), log(lenAge_aspx(a,s,p,x)), sigmaL );
              // Accumulate probability
              sumProbs += probLenAge_laspx(l,a,s,p,x);
            }
            
            // Now the RH tail
            if( l == L_s(s) - 1 ) 
              probLenAge_laspx(l,a,s,p,x) = Type(1.0) - sumProbs;

            // Calculate mean weight at age from probLenAge
            meanWtAge_aspx(a,s,p,x) += probLenAge_laspx(l,a,s,p,x) * Wlen_ls(l,s);
            // Calculate maturity at age from probLenAge if maturity
            // is length based
            if( matX == "length")
              if( (matLen_ls(l,s) > 0) )
                matAge_asp(a,s,p) += probLenAge_laspx(l,a,s,p,x) * matLen_ls(l,s);


          }

          // Calculate maturity at age if age based
          if( matX == "age" )
            matAge_asp(a,s,p) = 1 / (1 + exp( -1. * log(Type(19.0)) * ( a+1 - xMat50_s(s)) / (xMat95_s(s) - xMat50_s(s)) ) );

          // To compute ssbpr, we need to reduce by the fraction of spawners
          Surv_aspx(a,s,p,x) = exp(-a * M_spx(s,p,x));        
          if( a == A_s(s) - 1 ) 
            Surv_aspx(a,s,p,x) /= (1. - exp( -M_spx(s,p,x)) );
          // Compute ssbpr
          
        }
      }
    }

  SSBpr_asp = Surv_aspx.col(nX-1) * matAge_asp * meanWtAge_aspx.col(nX-1) / nX;

  // --------- Selectivity and Fishing Mortality -------- //
  array<Type> sel_lfspt(nL,nF,nS,nP,nT);
  array<Type> sel_afsptx(nA,nF,nS,nP,nT,nX);
  sel_lfspt.setZero();
  sel_afsptx.setZero();
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
            sel_lfspt.col(t).col(p).col(s).col(f).fill(1.);
            sel_lfspt.col(t).col(p).col(s).col(f) /= ( 1 + exp( - log(Type(19.)) * ( len - xSel50_spft(s,p,f,t) ) / xSelStep_spft(s,p,f,t)  ) );
            // convert selectivity-at-length to selectivity-at-age using probability matrix
            for( int x = 0; x < nX; x++)
              for( int a = 0; a < A_s(s); a++ )
              {
                vector<Type> probLenAgea_l(nL) ;
                probLenAgea_l.setZero();
                probLenAgea_l = probLenAge_laspx.col(x).col(p).col(s).col(a);
                sel_afsptx(a,f,s,p,t,x) = (probLenAgea_l*sel_lfspt.col(t).col(p).col(s).col(f)).sum();

              } 
          }
          if( selX == "age" )
            for( int x = 0; x < nX; x++)
              sel_afsptx.col(x).col(t).col(p).col(s).col(f) = 1/( 1 + exp( - log(Type(19.)) * ( age - xSel50_spft(s,p,f,t) ) / xSelStep_spft(s,p,f,t)  ) );


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
  B_asptx.setZero();
  N_asptx.setZero();
  R_spt.setZero();
  F_aspftx.setZero();
  Z_aspxt.setZero();
  C_aspftx.setZero();
  Cw_aspftx.setZero();
  B_spt.setZero();
  SB_spt.setZero();
  Bv_spft.setZero();
  Nv_spft.setZero();
  Nv_aspftx.setZero();
  N_asptx.setZero();
  predCw_spft.setZero();
  predC_spft.setZero();


  // Copies of state arrays for testing Baranov solver
  array<Type> B_aspxt(nA,nS,nP,nX,nT);
  array<Type> B_spxt(nS,nP,nX,nT);
  array<Type> vB_aspfxt(nA,nS,nP,nF,nX,nT);
  array<Type> vB_spfxt(nS,nP,nF,nX,nT);
  array<Type> vB_spft(nS,nP,nF,nT);
  array<Type> sel_aspfxt(nA,nS,nP,nF,nX,nT);
  array<Type> J_spft(nS,nP,nF,nT);

  B_aspxt.setZero();
  B_spxt.setZero();
  vB_aspfxt.setZero();
  vB_spfxt.setZero();
  vB_spft.setZero();
  sel_aspfxt.setZero();

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
        for( int x = 0; x < nX; x++)
        {
          // Initial time-step
          if( t == 0 )
          {
            // Initialise population either at unfished
            // or using the estimated fished initialisation
            // Populate first year - split into nX groups
            N_asptx.col(x).col(t).col(p).col(s) += R0_sp(s,p) * Surv_aspx.col(x).col(p).col(s) / nX;

            // // Non-eqbm initialisation
            if(swRinit_s(s) == 1)
              N_asptx.col(x).col(t).col(p).col(s) *= exp( omegaRinit_asp.col(p).col(s));

          }


          // time series history
          if( t > 0 )
          {
            // Generate recruitment
            Type SBt = SB_spt(s,p,t-1);
            N_asptx(0,s,p,t,x) = reca_sp(s,p) * SBt / (1 + recb_sp(s,p) * SBt)/nX;
            N_asptx(0,s,p,t,x) *= exp( omegaR_spt(s,p,t) );

            // Now loop over ages and apply fishing mortality (no discarding yet)
            for( int a = 1; a < A; a ++ )
            {
              // Update numbers according to age
              if( a > 0 )     
                N_asptx(a,s,p,t,x) = N_asptx( a-1, s, p, t-1 ,x) *  exp( - Z_aspxt(a,s,p,x,t-1) );
              if( a == A - 1) 
                N_asptx(a,s,p,t,x) += N_asptx(a,s, p, t-1 ,x) *  exp( - Z_aspxt(a,s,p,x,t-1) );
            }
          }

          // Save recruits in R_pt
          R_spt(s,p,t) += N_asptx(0,s,p,t,x);

          // Compute biomass at age and total biomass
          B_aspxt.col(t).col(x).col(p).col(s) += N_asptx.col(x).col(t).col(p).col(s) * meanWtAge_aspx.col(x).col(p).col(s);

          B_spxt(s,p,x,t) += B_aspxt.col(t).col(x).col(p).col(s).sum();


          // Loop over fleets and compute catch
          for( int f =0; f < nF; f++ )
          {
            // Calculate vulnerable numbers and biomass for this fleet
            Nv_aspftx.col(x).col(t).col(f).col(p).col(s) = N_asptx.col(x).col(t).col(p).col(s) * sel_afsptx.col(x).col(t).col(p).col(s).col(f);
            Nv_spft(s,p,f,t) += Nv_aspftx.col(x).col(t).col(f).col(p).col(s).sum();

            vB_aspfxt.col(t).col(x).col(f).col(p).col(s) += B_aspxt.col(t).col(x).col(p).col(s) * sel_afsptx.col(x).col(t).col(p).col(s).col(f);
            vB_spfxt(s,p,f,x,t) = vB_aspfxt.col(t).col(x).col(f).col(p).col(s).sum();
            vB_spft(s,p,f,t) = Bv_spft(s,p,f,t);
            
            
            // Calculate vulnerable biomass
            Bv_spft(s,p,f,t) += (Nv_aspftx.col(x).col(t).col(f).col(p).col(s) * meanWtAge_aspx.col(x).col(p).col(s)).sum();

            // Have to loop over ages here to avoid adding NaNs
            for(int a = 0; a < A; a++ )
            {
              vB_aspfxt(a,s,p,f,x,t) = sel_afsptx(a,f,s,p,t,x) * B_aspxt(a,s,p,x,t);
              sel_aspfxt(a,s,p,f,x,t) = sel_afsptx(a,f,s,p,t,x);
              // Compute probAgeAtLen for each time step
              for( int l = minL_s(s) - 1; l < L_s(s); l++ )
                probAgeLen_alspftx(a,l,s,p,f,t,x) += probLenAge_laspx(l,a,s,p,x) * N_asptx(a,s,p,t,x) * sel_lfspt(l,f,s,p,t);
            }

            

            // Renormalise probAgeAtLen
            for( int l = minL_s(s) - 1; l < L_s(s); l++ )
            {
              Type sumProbs = probAgeLen_alspftx.col(x).col(t).col(f).col(p).col(s).col(l).sum();
              if( sumProbs > 0 )
                probAgeLen_alspftx.col(x).col(t).col(f).col(p).col(s).col(l) /= sumProbs;
            }
            

          }

          // Compute total and spawning biomass
          B_asptx.col(x).col(t).col(p).col(s) = N_asptx.col(x).col(t).col(p).col(s) * meanWtAge_aspx.col(x).col(p).col(s);
          B_spt(s,p,t) += B_asptx.col(x).col(t).col(p).col(s).sum();
          // The females are the only important
          // mature individuals
          if( SB_spt(s,p,t) == 0)        
            SB_spt(s,p,t) = (B_asptx.col(nX - 1).col(t).col(p).col(s) * matAge_asp.col(p).col(s)).sum();


        }

      }
    }

    array<Type> tmpZ_aspx(nA,nS,nP,nX);
    tmpZ_aspx = Z_aspxt.col(t);

    array<Type> tmpF_spf(nS,nP,nF);
    tmpF_spf = F_spft.col(t);

    // Apply Baranov solver
    J_spft.col(t) = solveBaranov_spf( nBaranovIter,
                                      lambdaBaranovStep,
                                      A_s,
                                      C_spft.col(t),
                                      M_spx,
                                      B_aspxt.col(t),
                                      vB_aspfxt.col(t),
                                      vB_spfxt.col(t),
                                      vB_spft.col(t),
                                      sel_aspfxt.col(t),
                                      tmpZ_aspx,
                                      tmpF_spf);

    Z_aspxt.col(t) += tmpZ_aspx;
    F_spft.col(t)  += tmpF_spf;

    for( int s = 0; s < nS; s++)
      for( int p = 0; p < nP; p++ )
        for( int f = 0; f < nF; f++ )
          for( int x = 0; x < nX; x++ )
            for( int a = 0; a < A_s(s); a++)
            {
              F_aspftx(a,s,p,f,t,x)   = sel_afsptx(a,f,s,p,t,x) * F_spft(s,p,f,t);
              C_aspftx(a,s,p,f,t,x)   = N_asptx(a,s,p,t,x);
              C_aspftx(a,s,p,f,t,x)  *= (1.  - exp( -1. * Z_aspxt(a,s,p,x,t))); 
              C_aspftx(a,s,p,f,t,x)  *= F_aspftx(a,s,p,f,t,x); 
              C_aspftx(a,s,p,f,t,x)  /= Z_aspxt(a,s,p,x,t);  
              Cw_aspftx(a,s,p,f,t,x) = C_aspftx(a,s,p,f,t,x) * meanWtAge_aspx(a,s,p,x);
              // Generate predicted total catch in weight and numbers
              predCw_spft(s,p,f,t) += Cw_aspftx(a,s,p,f,t,x);
              predC_spft(s,p,f,t) += C_aspftx(a,s,p,f,t,x);

            }

  }

 

  
  // --------- Observation Model --------- //
  // Age/Length observations and CPUE
  array<Type> I_spft_hat(nS,nP,nF,nT);
  array<Type> aDist_aspftx_hat(nA,nS,nP,nF,nT,nX);
  array<Type> lDist_lspftx_hat(nL,nS,nP,nF,nT,nX);
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

          for( int x = 0; x < nX; x++)
            {
              // Create array of predicted age distributions - this should just be catch-at-age props
              // Calculate total catch
              Type totCatch = C_aspftx.col(x).col(t).col(f).col(p).col(s).sum();
                
              // Convert catch-at-age to proportions-at-age
              if(totCatch > 0)
                aDist_aspftx_hat.col(x).col(t).col(f).col(p).col(s) = C_aspftx.col(x).col(t).col(f).col(p).col(s) / totCatch;
                  

              // Create array of predicted length distributions
              // Need to do some thinking here... follow LIME to start with...
              // Probability of being harvested at age
              vector<Type> probHarvAge(nA);
              probHarvAge.setZero();
              probHarvAge = Nv_aspftx.col(x).col(t).col(f).col(p).col(s) / N_asptx.col(x).col(t).col(p).col(s).sum();
              
              // Probability of sampling a given length bin
              vector<Type> probHarvLen(L);
              probHarvLen.setZero();
              // loop over length and age, calculate probability of harvest at length
              for( int l = 0; l < L; l++ )
              {
                if( lenComps == "LIME" )
                {
                  for( int a = 0; a < A; a++ )
                    probHarvLen(l) += probHarvAge(a)*probLenAge_laspx(l,a,s,p,x);
                }

                if( lenComps == "Francis" )
                {
                  for( int a = 0; a < A; a++ )
                    probHarvLen(l) += probLenAge_laspx(l,a,s,p,x) * N_asptx(a,s,p,t,x) * sel_lfspt(l,f,s,p,t); 
                }

                // Save to length dists
                lDist_lspftx_hat(l,s,p,f,t,x) = probHarvLen(l);
              }
              // renormalise
              if(probHarvLen.sum() > 0)
                lDist_lspftx_hat.col(x).col(t).col(f).col(p).col(s) /= probHarvLen.sum();
            }
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
      vector<Type> initRecDevVec = omegaRinit_asp.col(p).col(s);
      recnll_sp(s,p) -= dnorm( initRecDevVec,Type(0), sigmaR_sp(s,p),true).sum();

       // then yearly recruitment deviations
      for( int t = 0; t < nT; t++)
        recnll_sp(s,p) -= dnorm( omegaR_spt(s,p,t), Type(0.), sigmaR_sp(s,p),true);
        
    }

  // vonB growth model likelihood //
  // Currently has lengths binned the same as the integration above, 
  // and uses a random-at-length likelihood to account for differences 
  // in fleet selectivity, and at different times
  array<Type> vonBnll_spf(nS,nP,nF);
  array<Type> ageAtLenResids_alspftx(nA,nL,nS,nP,nF,nT,nX);
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
  ageAtLenResids_alspftx.setZero();

  // Now loop and compute
  for( int s = 0; s < nS; s++ )
  {
    // containers for a given year/fleet/length/species/pop - function
    // expects vector<Type>
    vector<Type> obs(A_s(s));
    vector<Type> pred(A_s(s));
    vector<Type> resids(nA);
    for( int p = 0; p < nP; p++ )
      for( int f = 0; f < nF; f++)
      {
        for( int t = 0; t < nT; t++)
        {
          for( int l = minL_s(s) - 1; l < L_s(s); l++ )
          {
            for( int x = 0;  x < nX; x ++)
            {
              // Zero init containers
              obs.setZero();
              pred.setZero();
              resids.setZero();

              // Fill containers
              for(int a = 0; a < A_s(s); a++)
              {
                obs(a) += ALK_spalftx(s,p,a,l,f,t,x);
                pred(a) += probAgeLen_alspftx(a,l,s,p,f,t,x);
              }
              // Compute resids etc. if observations aren't missing
              if(obs.sum() > 0 &  pred.sum() > 0)
              {
                ageAtLenResids_alspftx.col(x).col(t).col(f).col(p).col(s).col(l).segment(0,A_s(s)) = calcLogistNormLikelihood(  obs, 
                                                                                                                                pred,
                                                                                                                                minPAAL,
                                                                                                                                etaSumSqAgeAtLen_spf(s,p,f),
                                                                                                                                nResidsAgeAtLen_spf(s,p,f) );

                nObsAgeAtLen_spf(s,p,f) += 1;
              }
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
  }

  // Observation likelihoods
  array<Type> CPUEnll_spf(nS,nP,nF);          // CPUE
  array<Type> ageCompsnll_spf(nS,nP,nF);      // Age observations
  array<Type> lenCompsnll_spf(nS,nP,nF);      // length observations
  // array<Type> Ctnll_sp(nS,nP);                // Catch
  array<Type> Dtnll_sp(nS,nP);                // Discards
  array<Type> residCPUE_spft(nS,nP,nF,nT);    // CPUE resids (concentrating variance parameter)
  array<Type> validIdxObs_spf(nS,nP,nF);      // valid observations for condMLE variance
  array<Type> ssrIdx_spf(nS,nP,nF);           // valid observations for condMLE variance
  // Conditional MLEs of age/length observation error variance
  array<Type> tau2Age_spf(nS,nP,nF);
  array<Type> tau2Len_spf(nS,nP,nF);
  array<Type> etaSumSqLen_spf(nS,nP,nF);
  array<Type> etaSumSqAge_spf(nS,nP,nF);
  array<Type> nResidsAge_spf(nS,nP,nF);
  array<Type> nResidsLen_spf(nS,nP,nF);
  array<Type> nObsAge_spf(nS,nP,nF);
  array<Type> nObsLen_spf(nS,nP,nF);
  array<Type> ageRes_aspftx(nA,nS,nP,nF,nT,nX);
  array<Type> lenRes_lspftx(nL,nS,nP,nF,nT,nX);
  // Zero-init
  CPUEnll_spf.setZero();
  ageCompsnll_spf.setZero();
  lenCompsnll_spf.setZero();
  // Ctnll_sp.setZero();
  Dtnll_sp.setZero();
  tau2Age_spf.setZero();
  tau2Len_spf.setZero();
  etaSumSqLen_spf.setZero();
  etaSumSqAge_spf.setZero();
  nResidsAge_spf.setZero();
  nResidsLen_spf.setZero();
  nObsAge_spf.setZero();
  nObsLen_spf.setZero();
  ageRes_aspftx.setZero();
  lenRes_lspftx.setZero();
  residCPUE_spft.setZero();
  validIdxObs_spf.setZero();
  ssrIdx_spf.setZero();

  for( int s = 0; s < nS; s++ )
    for( int p = 0; p < nP; p++ )
    {
      int A      = A_s(s) - minA_s(s) + 1;
      int L      = L_s(s) - minL_s(s) + 1;
      // Loop over fleets
      for( int f = 0; f < nF; f++ )
      {   
        // Now loop over time steps
        for( int t = 0; t < nT; t++ )
        {

          // Only use a year if the data exists and that fleet is used as a survey
          if( (I_spft(s,p,f,t) > 0) & (calcIndex_spf(s,p,f) > 0) )
          {
            // validIdxObs_spf(s,p,f) += 1;
            residCPUE_spft(s,p,f,t) = ( log(I_spft(s,p,f,t)) - log(I_spft_hat(s,p,f,t)) );
            // ssrIdx_spf(s,p,f) += square(residCPUE_spft(s,p,f,t));

            CPUEnll_spf(s,p,f) -= dnorm( log(I_spft(s,p,f,t)), log(I_spft_hat(s,p,f,t)), tauObs_spf(s,p,f),true);
          }

          for( int x = 0; x < nX; x++)
          {
            // tmp vectors to hold age and length observed and
            // predicted values
            vector<Type> fleetAgeObs(A);
            vector<Type> fleetAgePred(A);
            vector<Type> fleetAgeResids(A);

            vector<Type> fleetLenObs(L);
            vector<Type> fleetLenPred(L);  
            vector<Type> fleetLenResids(L);

            // Set tmp vectors to zero 
            fleetAgeObs.setZero();
            fleetAgePred.setZero();
            fleetAgeResids.setZero();
            fleetLenObs.setZero();
            fleetLenPred.setZero();
            fleetLenResids.setZero();

            // Get observations and predicted age comps
            fleetAgeObs   = age_aspftx.col(x).col(t).col(f).col(p).col(s).segment(minA_s(s)-1,A);
            fleetAgePred  = aDist_aspftx_hat.col(x).col(t).col(f).col(p).col(s).segment(minA_s(s)-1,A);
            // And length comps
            fleetLenObs   = len_lspftx.col(x).col(t).col(f).col(p).col(s).segment(minL_s(s)-1,L);
            fleetLenPred  = lDist_lspftx_hat.col(x).col(t).col(f).col(p).col(s).segment(minL_s(s)-1,L);

      
            // Now for compositional data. First, check if ages are 
            // being used (or exist), if so
            // compute logistic normal likelihood
            if( (fleetAgeObs.sum() > 0) & (fleetAgePred.sum() > 0) & (ageLikeWt > 0) )
            {
              ageRes_aspftx.col(x).col(t).col(f).col(p).col(s).segment(minA_s(s)-1,A)  = calcLogistNormLikelihood(  fleetAgeObs, 
                                                                                                                    fleetAgePred,
                                                                                                                    minAgeProp,
                                                                                                                    etaSumSqAge_spf(s,p,f),
                                                                                                                    nResidsAge_spf(s,p,f) );

              // Increment number of ages
              nObsAge_spf(s,p,f) += 1;

            }     

            // If ages aren't being used, but lengths exist, then compute 
            // logistic normal likelihood
            if( (fleetLenObs.sum() > 0) & (fleetLenPred.sum() > 0) &  (lenLikeWt > 0) )
            {
              lenRes_lspftx.col(x).col(t).col(f).col(p).col(s).segment(minL_s(s)-1,L) = calcLogistNormLikelihood( fleetLenObs, 
                                                                                                                  fleetLenPred,
                                                                                                                  minLenProp,
                                                                                                                  etaSumSqLen_spf(s,p,f),
                                                                                                                  nResidsLen_spf(s,p,f) );
              
              nObsLen_spf(s,p,f) += 1;
            }
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

      }
    }


  // Shared Hierarchical Prior Distributions //
  // loop over species and stocks, add penalty for stock deviations from
  // species mean

  // Synoptic catchability prior
  vector<Type> qnlpSyn_s(nS);
  Type qnlpSyn = 0.;
  qnlpSyn_s.setZero();
  Type qnlpSurv = 0.;

  Type Fnlp = 0;



  // Now a prior on catchability and observation error SD
  // Surveys only really act on one stock each, 
  // so their catchabiliity prior will be across species
  vector<Type> tauObsnlp_f(nF);
  tauObsnlp_f.setZero();
  for( int s = 0; s < nS; s++ )
    for( int p = 0; p < nP; p++ )
      for( int f = 0; f < nF; f++ )
      {

        // if( Fbar_spf(s,p,f) > 0)
        //   Fnlp += log(Fbar_spf(s,p,f));
        // Penalise fleet specific catchabilities if
        // they are a calculated index
        // Synoptic survey
        if( group_f(f) == 2)
          if( calcIndex_spf(s,p,f) > 0)
            qnlpSyn_s(s) -= dnorm( log(q_spf(s,p,f)), lnqbarSyn_s(s), tauqSyn_s(s), true);
    
        // HS Assemblage survey
        if( group_f(f) == 1)
          if( calcIndex_spf(s,p,f) > 0)
            qnlpSurv -= dnorm( log(q_spf(s,p,f)), log(mqSurveys), sdqSurveys,true );
        
        // Observation error SD
        if( calcIndex_spf(s,p,f) > 0 & tauObs_spf(s,p,f) > 0 )
          tauObsnlp_f(f) += (IGatau_f(f)+Type(1))*2*lntauObs_spf(s,p,f)+IGbtau_f(f)/(tau2Obs_spf(s,p,f));
        
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


  // Steepness
  vector<Type>  steepnessnlp_s(nS);
  array<Type>   steepnessnlp_sp(nS,nP);
  Type steepnessnlp = 0.;

  // Mortality
  vector<Type>  Mnlp_s(nS);
  array<Type>   Mnlp_sp(nS,nP);
  array<Type>   Mnlp_spx(nS,nP,nX);
  Type Mnlp = 0.;

  // VonB priors
  // Growth rate
  vector<Type>  vonKnlp_s(nS);
  vector<Type>  vonKnlp_p(nP);
  Type          vonKnlp = 0.;
  // L2
  vector<Type>  L2nlp_s(nS);
  vector<Type>  L2nlp_p(nP);
  Type          L2nlp = 0.;

  // Zero initialise
  vonKnlp_s.setZero();
  vonKnlp_p.setZero();
  L2nlp_s.setZero();
  L2nlp_p.setZero();


  deltaL2bar_p.setZero();
  deltaVonKbar_p.setZero();
  // Calculate stock mean delta values for
  // growth model, to extend growth pars to 
  // species/stock combinations without data
  if(nS > 1)
    for( int p = 0; p < nP; p++ )
    {
      deltaVonKbar_p(p) = deltaVonK_sp.col(p).sum()/nP;
      deltaL2bar_p(p) = deltaL2_sp.col(p).sum()/nP;

      vector<Type> vonKVec  = deltaVonK_sp.col(p);
      vector<Type> L2Vec    = deltaL2_sp.col(p);

  //     // Within stock growth par deviation priors - essentially
  //     // a regularisation to stop the unobserved stocks from overfitting  
  //     // May be leading to undifferentiable posterior due to
  //     // the multiple local minima - think more on this
  //     // vonKnlp_p(p)          -= dnorm( vonKVec, deltaVonKbar_p(p), sigmavonK, true).sum();
  //     // L2nlp_p(p)            -= dnorm( L2Vec, deltaL2bar_p(p), sigmaL2, true).sum();

    }


  // Mortality prior
  // vector<Type>  Mnlp_s(nS);
  // vector<Type>  Mnlp_p(nP);
  // array<Type>   Mnlp_sp(nS,nP);
  // 0-initialise
  steepnessnlp_s.setZero();
  steepnessnlp_sp.setZero();
  Mnlp_s.setZero();
  Mnlp_sp.setZero();
  Mnlp_spx.setZero();
  Type sel_nlp = 0.;


  for( int s = 0; s < nS; s++ )
  {
    vector<Type> vonKVec  = deltaVonK_sp.transpose().col(s);
    vector<Type> L2Vec    = deltaL2_sp.transpose().col(s);
    vonKnlp_s(s)         -= dnorm( vonKVec, Type(0), Type(1), true).sum();
    L2nlp_s(s)           -= dnorm( L2Vec, Type(0), Type(1), true).sum();

    for( int p = 0; p < nP; p ++)
    {
      steepnessnlp_sp(s,p)  -= dnorm( epsSteep_sp(s,p), Type(0), Type(1), true);
      Mnlp_sp(s,p)          -= dnorm( epsM_sp(s,p), Type(0), Type(1), true);

      if( nX > 1 )
        for( int x = 0; x < nX; x++ )
        {
          vonKnlp_s(s)  -= dnorm( deltaVonK_sx(s,x), Type(0), Type(1), true);
          L2nlp_s(s)    -= dnorm( deltaL2_sx(s,x), Type(0), Type(1), true);
          Mnlp_sp(s,p)  -= dnorm( epsM_sx(s,x), Type(0), Type(1), true);

        }
    }

  }
  
  // add species level h prior
  steepnessnlp_s  -= dnorm( epsSteep_s, 0., Type(1), true );
  Mnlp_s          -= dnorm( epsM_s, 0., Type(1), true );
  
  // Now prior on complex mean M
  // Currently have sigmaM/sigmah doing double duty, replace
  // with another hyperparameter - we might want to estimate
  // the complex variance later
  Mnlp            -= dnorm( lnM, log(muM), sigmaM, true);
  steepnessnlp    -= dnorm( logitSteep, logit_muSteep, sigmah, true);
    
  // And penalties stock mean growth model
  // deviations
  L2nlp           -= dnorm( deltaL2bar_p, Type(0.), sigmaL2, true).sum();
  vonKnlp         -= dnorm( deltaVonKbar_p, Type(0.), sigmavonK, true).sum();

  // Add time-varying selectivity deviations
  sel_nlp -= dnorm( epsxSel50_vec, Type(0), sigmaSel, true).sum();
  sel_nlp -= dnorm( epsxSelStep_vec, Type(0), sigmaSel, true).sum();
  // Prior on xSel95_sf
  for(int s = 0; s < nS; s++)
    for( int f = 0; f < nF; f++)
    {
      sel_nlp -= dnorm( lnxSel50_sf(s,f), pmlnxSel50_sf(s,f), cvxSel, true);
      sel_nlp -= dnorm( lnxSelStep_sf(s,f), pmlnxSelStep_sf(s,f), cvxSel, true);

      for( int p = 0; p < nP; p++)
      {
        sel_nlp -= dnorm( epsxSel50_spf(s,p,f), Type(0), cvxSel, true);
        sel_nlp -= dnorm( epsxSelStep_spf(s,p,f), Type(0), cvxSel, true);  
      }
      
    }
  

  // VonB priors
  vector<Type>  L1nlp_s(nS);
  // Zero-init
  L1nlp_s.setZero();

  // Compute
  L1nlp_s   -= dnorm( lnL1_s, pmlnL1_s, cvL1, true);
  L2nlp_s   -= dnorm( log(L2_s), pmlnL2_s, cvL2, true);
  vonKnlp_s -= dnorm( lnvonK_s, pmlnVonK, cvVonK, true);


  array<Type> B0nlp_sp(nS,nP);
  B0nlp_sp.setZero();
  B0nlp_sp += lambdaB0 * log( B_spt.col(0) );
    
  

  
  // Observations
  // f += Ctnll_sp.sum();       // Catch
  f += ageLikeWt * lenCompsnll_spf.sum(); // Length compositions
  f += lenLikeWt * ageCompsnll_spf.sum(); // Age compositions
  f += idxLikeWt * CPUEnll_spf.sum();     // Survey CPUE
  f += sel_nlp;
  f += tauObsnlp_f.sum();
  // Growth model
  f += growthLikeWt * vonBnll_spf.sum();     // Growth model
  f += growthLikeWt * (L1nlp_s.sum() + L2nlp_s.sum() + vonKnlp_s.sum());
  f += growthLikeWt * ( L2nlp_p.sum() + vonKnlp_p.sum());
  f += growthLikeWt * ( L2nlp + vonKnlp );
  // Recruitment errors
  f += recnll_sp.sum();      // recruitment process errors
  // q groups
  f += qnlpSurv + qnlpSyn_s.sum() + qnlpSyn + qnlp_tv;
  
  // Biological parameters
  f += steepnessnlp_s.sum() + steepnessnlp_sp.sum() + steepnessnlp;
  f += Mnlp_s.sum() + Mnlp_sp.sum() + Mnlp;
  f += B0nlp_sp.sum();
  // joint_nlp += pop_nlp + spec_nlp;

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
  REPORT( h_sp );     
  REPORT( M_spx );     
  REPORT( L2_spx );  
  REPORT( vonK_spx );  
  REPORT( L1_spx ); 
  REPORT( sigmaLa_s );
  REPORT( sigmaLb_s );
  REPORT( sigmaR_sp );

  // Fishery/Survey model pars
  REPORT( q_spf );          // Catchability
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
  REPORT( sel_lfspt );      // Selectivity at length
  REPORT( sel_afsptx );      // Selectivity at age
  REPORT( Fbar_spf );       // Average fishing mortality
  REPORT( tauC_f );
  REPORT( tauD_f );

  // Model states
  REPORT( B_asptx );
  REPORT( N_asptx );
  REPORT( R_spt );
  REPORT( F_aspftx );
  REPORT( Z_aspxt );
  REPORT( C_aspftx );
  REPORT( Cw_aspftx );
  REPORT( predCw_spft );
  REPORT( predC_spft );
  REPORT( Bv_spft );
  REPORT( Nv_spft );
  REPORT( Nv_aspftx );
  REPORT( SB_spt );
  REPORT( B_spt );

  // Stock recruit parameters
  REPORT( Surv_aspx );
  REPORT( SSBpr_asp );
  REPORT( R0_sp );
  REPORT( phi_sp );
  REPORT( reca_sp );
  REPORT( recb_sp );

  // Growth and maturity
  REPORT( probLenAge_laspx );
  REPORT( probAgeLen_alspftx );
  REPORT( lenAge_aspx );
  REPORT( Wlen_ls );
  REPORT( matLen_ls );
  REPORT( matAge_asp );
  REPORT( xMat50_s );
  REPORT( xMat95_s );
  REPORT( meanWtAge_aspx );
  REPORT( ageAtLenResids_alspftx );
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
  REPORT( deltaVonK_sx );
  REPORT( deltaL2_sx );
  REPORT( deltaVonKbar_p );
  REPORT( deltaL2bar_p );

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
  REPORT( pmlnL2_s );
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
  REPORT( omegaRinit_asp );
  REPORT( omegaRmat_spt );

  // Species/stock effects
  REPORT( epsM_sx );
  REPORT( epsM_sp );
  REPORT( epsM_s );
  REPORT( epsSteep_sp );
  REPORT( epsSteep_s );


  // Likelihood values
  REPORT( joint_nlp );
  REPORT( recnll_sp );
  REPORT( vonBnll_spf );
  REPORT( CPUEnll_spf );
  REPORT( tau2Obs_spf );
  REPORT( tauObs_spf );
  REPORT( tauObsnlp_f );
  REPORT( validIdxObs_spf );
  REPORT( residCPUE_spft );
  REPORT( ssrIdx_spf );
  REPORT( ageCompsnll_spf );
  REPORT( etaSumSqAge_spf );
  REPORT( nResidsAge_spf );
  REPORT( nObsAge_spf );
  REPORT( tau2Age_spf );
  REPORT( ageRes_aspftx );
  REPORT( lenCompsnll_spf );
  REPORT( etaSumSqLen_spf );
  REPORT( nResidsLen_spf );
  REPORT( nObsLen_spf );
  REPORT( tau2Len_spf );
  REPORT( lenRes_lspftx );
  REPORT( steepnessnlp_sp );
  REPORT( steepnessnlp_s );
  REPORT( steepnessnlp );
  REPORT( Mnlp_sp );
  REPORT( Mnlp_s );
  REPORT( Mnlp );
  REPORT( qnlpSurv );
  REPORT( qnlpSyn );
  REPORT( qnlpSyn_s );
  REPORT( qnlp_tv );
  REPORT( tauObsnlp_f );
  REPORT( minAgeProp );
  REPORT( minLenProp );
  REPORT( L1nlp_s );
  REPORT( L2nlp_s );
  REPORT( L2nlp_p );
  REPORT( vonKnlp_s );
  REPORT( vonKnlp_p );
  REPORT( sel_nlp );
  REPORT( L2nlp );
  REPORT( vonKnlp );

  // Data
  REPORT( I_spft );
  REPORT( C_spft );
  REPORT( D_spft );
  REPORT( ALK_spalftx );
  REPORT( age_aspftx );
  REPORT( len_lspftx );
  REPORT( group_f );
  REPORT( A_s );
  REPORT( minA_s );
  REPORT( L_s );
  REPORT( minL_s );
  REPORT( lenD_s );
  REPORT( calcIndex_spf );
  REPORT( tvSelFleets );
  REPORT( tvqFleets );
  REPORT( regFfleets );

  // Predicted observations
  REPORT( lDist_lspftx_hat );
  REPORT( aDist_aspftx_hat );
  REPORT( I_spft_hat );

  // Reordered arrays
  REPORT( B_aspxt );
  REPORT( vB_aspfxt );
  REPORT( vB_spfxt );
  REPORT( vB_spft );
  REPORT( sel_aspfxt );
  REPORT( J_spft );

  // Echo switches
  REPORT( swRinit_s );
  REPORT( parSwitch );
  REPORT( calcIndex_spf );
  REPORT( tvSelFleets );
  REPORT( tvqFleets );
  REPORT( regFfleets );
  REPORT( idxLikeWt );
  REPORT( ageLikeWt );
  REPORT( lenLikeWt );
  REPORT( growthLikeWt );
  REPORT( tFirstRecDev_s );
  REPORT( tLastRecDev_s );
  REPORT( minAgeProp );
  REPORT( minLenProp );
  REPORT( minPAAL );
  REPORT( lambdaB0 );
  REPORT( minTimeIdx_spf );

  // REPORT(A);

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
  // ADREPORT(lntauObs_spf)        // Observation error
  // ADREPORT(lnxSel50_sf)   // length at 50% sel
  // ADREPORT(lnxSelStep_sf)   // length at 95% sel
  // ADREPORT(lnF_spft)         // Fishing mortality


  return( f );

}

