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
// Intended features:
// - Multistock and multispecies (DERPA)
//     - Species made up of a single stock have only species level 
//       priors
//     - Different initialisation times for each stock - wishlist
// - Multi-sector
//     - Fleets are fisheries by gear type (including surveys)
// - Multi-level RH priors on:
//     - Growth (vonB or F-W undecided)
//     - Fishing mortality (correlation in REs if estimated)
//     - Natural mortality (M0 prior and correlated deviations)
//     - Selectivity
//     - Catchability
//     - S-R Steepness
// - Age or Length observations are used (Multinomial or Dirichlet Multinomial)
//     - Age observations are used where available
//     - When ages not available (or effective sample size too low), lengths
//        are compared to a predicted length dist
//     - Currently, age-length distribution is stationary, could
//        augment to explicitly model growth-groups (like Atl Halibut)
// - Integrated growth model to predict length dist when age data missing
// - Discarding - need to talk to fisheries about discarding behaviour
//     - "Grading length" for fisheries with no size limit
// 
// ><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><><>><><>><><>><>

#include <TMB.hpp>                                // Links in the TMB libraries
#include <iostream>
// link personal library of TMB stock assessment functions
#include "/Users/sdnjohnson/Work/code/tmbFuns/stockAssessmentFuns.hpp"  

// objective function
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Call namespaces //
  using namespace density;

  /*data section*/
  // Data Structures
  DATA_ARRAY(I_pft);                // CPUE data
  DATA_ARRAY(C_pft);                // Catch data (biomass)
  DATA_ARRAY(D_pft);                // Discard data (biomass)
  DATA_ARRAY(ALK_alp);              // Age-length observations (in freq) by pop - time averaged
  DATA_ARRAY(age_afpt);             // Age observations for each population and fleet over time (-1 missing)
  DATA_ARRAY(length_lfpt);          // Length observations for each population and fleet over time (-1 missing)
  DATA_IVECTOR(s_p);                // species ID for each population e.g. (1,1,2,2,3)
  DATA_IVECTOR(type_f);             // Fleet type (0 survey, 1 fishery dependent)
  DATA_IVECTOR(A_s);                // +Group age by species
  DATA_IVECTOR(lenD_s);             // Discard length by species
  // DATA_MATRIX(intMatrix_pp);        // population interactions mtx (spawning?, adjacency)
  DATA_IVECTOR(nLenBins_s);         // Number of length bins by species
  DATA_ARRAY(lenBinHi_sl);          // UBs of length bins by species
  DATA_ARRAY(lenBinMids_sl);        // Middle of length bins by species

  DATA_INTEGER(maxL);               // highest number of length bins to aggregate into (taken over species)
  DATA_INTEGER(maxA);               // highest number of age classes (taken over species)

  // Model dimensions
  int nP = I_pft.dim(0);            // No. of stocks (populations)
  int nG = I_pft.dim(1);            // No. of fleets (surveys + trawl)  
  int nT = I_pft.dim(2);            // No of time steps
  int nS = A_s.size();              // No of species
  int nL = maxL;                    // max no of length bins (for creating state arrays)
  // int maxA = A_s.max();             // max no. of age classes

  // Model switches - fill these in when we know what we want
  DATA_IVECTOR(swRinit_p)           // pop-spec fished initialisation switch

  /*parameter section*/
  // Leading Parameters
  // Biological
  PARAMETER_VECTOR(lnB0_p);         // Biomass at unfished for ind. stocks
  PARAMETER_VECTOR(logitSteep_p);   // pop specific steepness
  PARAMETER_VECTOR(lnM_p);          // stock-specific M
  PARAMETER_ARRAY(lnRinit_p);       // Pop-specific initial R (non-eq init)
  PARAMETER_VECTOR(lnLinf_p);       // pop-specific vonB Linf parameter
  PARAMETER_VECTOR(lnvonK_p);       // pop-specific vonB K parameter
  PARAMETER_VECTOR(lnvont0_p);      // pop-specific vonB t0 parameter
  PARAMETER_ARRAY(LWa_s);           // species L-W conversion a par (units conv: cm -> kt)
  PARAMETER_ARRAY(LWb_s);           // species L-W conversion b par (length to volume)
  PARAMETER_VECTOR(lenMat50_s);     // length at 50% maturity
  PARAMETER_VECTOR(lenMat95_s);     // length at 95% maturity
  // Observation model
  // survey //
  PARAMETER_ARRAY(lnq_fp);          // fleet-Pop specific catchability
  PARAMETER_ARRAY(lntau_fp);        // fleet-species specific obs error sd
  // Fishery model
  // selectivity by fleet //
  PARAMETER_ARRAY(lnlenSel50_sf);   // length at 50% selectivity by species/fleet
  PARAMETER_ARRAY(lnlenSel95_sf);   // length at 95% selectivity by species/fleet
  // Fishing mortality
  PARAMETER_ARRAY(lnF_pft);         // Fishing mortality by population/fleet/time
  PARAMETER_VECTOR(lntauC_f);       // Catch observation SD by fleet
  PARAMETER_VECTOR(lntauD_f);       // Discards observation SD by fleet

  // Priors
  // Growth //
  PARAMETER_VECTOR(muLinf_s);       // species mean Linf
  PARAMETER_VECTOR(sigmaLinf_s);    // species Linf sd    
  PARAMETER_VECTOR(muvonK_s);       // species mean vonK  
  PARAMETER_VECTOR(sigmavonK_s);    // species vonK sd       
  PARAMETER_VECTOR(mut0_s);         // species mean t0  
  PARAMETER_VECTOR(sigmat0_s);      // species t0 sd       
  PARAMETER_VECTOR(lnsigmaL_s);     // age-length observation error SDs
  // selectivity //
  PARAMETER_VECTOR(mulenSel50_f);    // mean length at 50% selectivity by fleet
  PARAMETER_VECTOR(mulenSel95_f);    // mean length at 95% selectivity by fleet
  PARAMETER_VECTOR(sigmalenSel50_f); // SD in length at 50% selectivity by fleet
  PARAMETER_VECTOR(sigmalenSel95_f); // SD in length at 95% selectivity by fleet
  // catchability //
  PARAMETER_ARRAY(lnqbar_fs);       // fleet/species level catchability
  PARAMETER_ARRAY(lntauq_fs);       // fleet/species specific catchability SD
  PARAMETER_VECTOR(lnqbar_f);       // fleet level catchability (mean trawl survey efficiency?)
  PARAMETER_VECTOR(lntauq_f);       // fleet specific catchability SD
  // steepness //
  PARAMETER_VECTOR(logitSteep_s);   // species steepness
  PARAMETER_VECTOR(lnsigmaSteep_s); // species level steepness SD
  PARAMETER(logitSteep);            // complex steepness
  PARAMETER(lnsigmaSteep);          // complex level steepness SD
  // Natural Mortality //
  PARAMETER_VECTOR(lnM_s);          // Species Complex mean M
  PARAMETER_VECTOR(lnsigmaM_s);     // species M SD
  PARAMETER(ln_muM);                // Multispecies assemblage mean M
  PARAMETER(lnsigmaM);              // Assemblage M SD


  // Random Effects
  // recruitment //
  PARAMETER_ARRAY(omegaR_pt);       // stock specific recruitment errors 2:nT
  PARAMETER_ARRAY(omegaRinit_pa);   // stock-age specific recruitment initialisation errors
  PARAMETER_VECTOR(lnsigmaR_p);     // stock recruitment errors sd (sqrt cov matrix diag)
  PARAMETER_VECTOR(logitRCorr_chol);// stock recruitment errors corr chol factor off diag  
  PARAMETER_VECTOR(logitRgamma_p);  // stock-specific AR1 auto-corr on year effect (omega)
  // mortality deviations //
  PARAMETER_ARRAY(epsilonM_pt);     // M deviations by population, 2:nT
  PARAMETER_VECTOR(lnsigmaM_p);     // M deviation pop-specific SD
  PARAMETER_VECTOR(logitMgamma_p);  // AR1 auto-corr on M devs effect (omega)


  /* Procedure Section */
  // Exponentiate leading parameters
  // Biological //
  vector<Type>  B0_p          = exp(lnB0_p);
  vector<Type>  h_p           = Type(1.0) / ( Type(1.0) + exp( -logitSteep_p) ) ;
  vector<Type>  M_p           = exp(lnM_p);
  vector<Type>  Rinit_p       = exp(lnRinit_p);
  vector<Type>  Linf_p        = exp(lnLinf_p);
  vector<Type>  vonK_p        = exp(lnvonK_p);
  vector<Type>  vont0_p       = exp(lnvont0_p);
  vector<Type>  sigmaL_s      = exp(lnsigmaL_s);
  vector<Type>  sigmaR_p      = exp(lnsigmaR_p);
  // Observation model
  array<Type>   q_fp(nG,nP);
  array<Type>   tau_fp(nG,nP);
  q_fp.setZero();
  tau_fp.setZero();
  for( int f = 0; f < nG; f++ )
  {
    q_fp(f)         = exp( lnq_fp(f) );
    tau_fp(f)       = exp( lntau_fp(f) );
  } 
  // Fishery model
  // selectivity by fleet //
  vector<Type>  lenSel50_sf   = exp( lnlenSel50_sf);
  vector<Type>  lenSel95_sf   = exp( lnlenSel95_sf);
  // Fishing mortality
  array<Type>   F_pft(nP,nG,nT);
  F_pft.setZero();
  for( int p = 0; p < nP; p++ ) F_pft(p) = exp( lnF_pft(p));
  vector<Type>  tauC_f        = exp(lntauC_f);    
  vector<Type>  tauD_f        = exp(lntauD_f);    
  
  // Set up model state arrays
  array<Type> B_apt(maxA,nP,nT);      // Biomass at age-pop-time
  array<Type> N_apt(maxA,nP,nT);      // Numbers at age-pop-time
  array<Type> R_pt(nP,nT);            // Recruits by pop-time
  array<Type> F_afpt(maxA,nG,nP,nT);  // Fishing mortality by age, fleet, population and time
  array<Type> Z_apt(maxA,nG,nP,nT);   // Total mortality by age, population and time
  array<Type> C_afpt(maxA,nG,nP,nT);  // Total catch (in weight) by age, fleet, population and time
  array<Type> Cw_afpt(maxA,nP,nT);    // Total catch (in weight) by age, population and time

  // derived parameters
  // Stock recruitment //
  vector<Type>  R0_p(nP);             // eqbm recruitment
  vector<Type>  phi_p(nP);            // eqbm SSB per recruit
  vector<Type>  reca_p(nP);           // BH a parameter for pops
  vector<Type>  recb_p(nP);           // BH b parameter for pops

  // Growth //
  array<Type>   Wlen_ls(nL,nS);       // weight-at-length by species
  array<Type>   lenAge_ap(maxA,nP);   // Mean length-at-age by population

  // Maturity //
  array<Type>   matAge_ap(maxA,nP);   // Proportion mature at age by population (growth differences between pops)
  array<Type>   matLen_ls(maxL,nS);   // Proportion mature at Length by species

  // Prior Hyperparameters //
  // Steepness
  vector<Type>  h_s           = Type(1.0) / ( Type(1.0) + exp( -logitSteep_s) ) ;
  Type          muSteep       = Type(1.0) / ( Type(1.0) + exp( -logitSteep) ) ;
  // Natural mortality
  vector<Type>  M_s           = log(lnM_s);
  vector<Type>  sigmaM_s      = log(lnsigmaM_s);
  Type          muM           = log(ln_muM);
  Type          sigmaM        = log(lnsigmaM);


  // ---------------- Growth, Maturity, Stock-Recruitment ----------- //
  Wlen_ls.fill(-1.0);
  matLen_ls.fill(-1.0);
  for( int s = 0; s < nS; s++ )
  {
    // weight-at-length for each length bin midpoint
    for( int l = 0; l < nLenBins_s(s); l++ )
    {
      Wlen_ls(l,s) = LWa_s(s) * pow( lenBinMids_sl(s,l), LWb_s(s) );
    }
    // proportion mature-at-length
    for( int l = 0; l < maxL; l ++ )
    {
      Type len = lenBinMids_sl(s,l);
      if( len > 0 ) // Only calculate maturity if the length bin has a midpoint, -1 == NA
        matLen_ls(l,s) = 1 / (1 + exp( log(Type(19.0)) * (len - lenMat50_s(s)) / (lenMat95_s(s) - lenMat50_s(s)) ) );
    }
  }
  // Now produce the probability of length-at-age matrix for each stock
  // (taken from LIME by Rudd and Thorson, 2017)
  // In the same loop, we're going to compute the ssbpr parameter phi
  // for each stock
  array<Type> probLenAge_lap(nL,maxA,nP);
  array<Type> meanWtAge_ap(maxA,nP);
  lenAge_ap.setZero();
  probLenAge_lap.setZero();
  meanWtAge_ap.setZero();
  phi_p.setZero();
  matAge_ap.setZero();

  for( int p = 0; p < nP; p++ )
  {
    // Grab species ID, number of length bins and Linf CV
    int   specID    = s_p(p);
    int   nLbins    = nLenBins_s( specID );
    Type  lenAgeCV  = sigmaL_s( specID ) / Linf_p ( p );
    Type  sumProbs  = 0;
    for( int a = 0; a < A_s(specID); a ++ )
    {
      lenAge_ap(a,p) = Linf_p(p) * ( 1 - exp( -1. * vonK_p(p) * ( a + 1 - vont0_p ( p ) ) ) ); 
      for( int l = 0; l < nLbins; l++ )
      {
        if( l == 0 )
        {
          // Compute LH tail of the distribution
          probLenAge_lap(l,a,p) = pnorm( lenBinHi_sl(specID,l), lenAge_ap(a,p), lenAge_ap(a,p) * lenAgeCV );
          sumProbs += probLenAge_lap(l,a,p);
        }
        if( (l >= 1) & (l < nLbins - 1 ) )
        {
          // Compute interior bins
          probLenAge_lap(l,a,p) = pnorm( lenBinHi_sl(specID,l), lenAge_ap(a,p), lenAge_ap(a,p) * lenAgeCV );
          probLenAge_lap(l,a,p) -= pnorm( lenBinHi_sl(specID,l-1), lenAge_ap(a,p), lenAge_ap(a,p) * lenAgeCV );
          sumProbs += probLenAge_lap(l,a,p);
        }
        // Now the RH tail
        if( l == nLbins - 1 ) probLenAge_lap(l,a,p) = Type(1.0) - sumProbs;
        meanWtAge_ap(a,p) += probLenAge_lap(l,a,p) * Wlen_ls(l,specID);
        if( (matLen_ls(l,specID) > 0) )
          matAge_ap(a,p) += probLenAge_lap(l,a,p) * matLen_ls(l,specID);
      }

      // To compute ssbpr, we need to borrow N_apt for a moment 
      Type ssbpr_a = exp(-a * M_p(p)) * matAge_ap(a,p) * meanWtAge_ap(a,p);
      if ( a == A_s(specID) - 1 ) ssbpr_a /= (1. - exp( -M_p(p)) );
      phi_p(p) += ssbpr_a;
    }
  }



  // --------- Selectivity and Fishing Mortality -------- //
  array<Type> sel_lfp(nL,nG,nP);
  array<Type> sel_afp(maxA,nG,nP);
  sel_lfp.setZero();
  sel_afp.setZero();
  for( int p = 0; p < nP; p++ )
  {
    int specID  = s_p(p);
    int A       = A_s(specID);
    int nLbins  = nLenBins_s(specID);
    // Loop over fleets, create selectivities
    for( int f = 0; f < nG; f++ )
    {
      // selectivity-at-length
      for( int l = 0; l < nLbins; l++ )
      {
        sel_lfp(l,f,p) = 1;
        sel_lfp(l,f,p) /= ( 1 + exp( - log(Type(19.)) * ( lenBinMids_sl(specID,l) - lenSel50_sf(specID,f) ) / ( lenSel95_sf(specID,f) - lenSel50_sf(specID,f) ) ) );
      }
      // convert selectivity-at-length to selectivity-at-age using probability matrix
      for( int a = 0; a < A; a++ )
      {
        vector<Type> probLenAgea_l(nL) ;
        probLenAgea_l.setZero();
        probLenAgea_l = probLenAge_lap.col(p).col(a);
        sel_afp(a,f,p) = (probLenAgea_l*sel_lfp.col(p).col(f)).sum();
      }
    }
  }

  // --------- Population Dynamics ----------- //
  // First, compute eqbm recruitment
  R0_p = B0_p / phi_p;

  // Set all state arrays to zero
  B_apt.setZero();
  N_apt.setZero();
  R_pt.setZero();
  F_afpt.setZero();
  Z_apt.setZero();
  C_afpt.setZero();
  Cw_afpt.setZero();

  for( int p = 0; p < nP; p++ )
  {
    // Get species ID, set plus group age
    int specID = s_p(p);
    int A    = A_s(specID);

    for( int t = 0; t < nT; t++ )
    {
      // Initial time-step
      if( t == 0 )
      {
        // Initialise population either at unfished
        // or using the estimated fished initialisation
        Type InitRtmp;
        if( swRinit_p(p) == 0 ) InitRtmp = R0_p(p);
        else InitRtmp = Rinit_p(p);

        // Populate first year
        for( int a = 0; a < A; a++ )
        {
          // Produce recruits
          N_apt(a,p,t) = InitRtmp * exp(-Type(a) * M_p(p) ) * exp( omegaRinit_pa(p,a) - square(sigmaR_p(p)) / 2 );
          // Compute fishing mortality at age
          F_afpt.col(t).col(p)(a) = sel_afp.col(p)(a) * F_pft.col(t)(p);
          // Compute total mortality at age (for catch later)
          Z_apt(a,p,t) = M_p(p);
          for( int f = 0; f < nG; f++ ) Z_apt(a,p,t) += F_afpt(a,f,p,t);
        }
        // Correct plus-group numbers at eqbm if initialised at unfished
        if( swRinit_p(p) == 0 ) N_apt(A-1,p,t) /= (Type(1) - exp(-M_p(p)));
      }

      // Projections
      if( t > 0 )
      {
        // Generate recruitment
        Type SBt = B_apt.col(t-1).col(p).sum();
        N_apt(0,p,t) = Type(4) * R0_p(p) * SBt / ( B0_p(p) * (Type(1) - h_p(p) ) + (Type(5)*h_p(p) - Type(1)) * SBt );
        N_apt(0,p,t) *= exp( omegaR_pt(p,t) - square(sigmaR_p(p)) / 2 );

        // Now loop over ages and apply fishing mortality (no discarding yet)
        for( int a = 0; a < A-1; a ++ )
        {
          // Compute fishing mortality at age
          F_afpt.col(t).col(p)(a) = sel_afp.col(p)(a) * F_pft.col(t)(p);
          // Compute total mortality at age
          Z_apt(a,p,t) = M_p(p);
          for( int f = 0; f < nG; f++ ) Z_apt(a,p,t) += F_afpt(a,f,p,t);

          // Update numbers according to age
          if( a == 0 )    N_apt(0,p,t) = R_pt(p,t);
          if( a > 0 )     N_apt(a,p,t) = N_apt( a-1, p, t-1 ) *  exp( - M_p(p) - Z_apt(a,p,t-1) );
          if( a == A - 1) N_apt(a,p,t) += N_apt( a, p, t-1 ) *  exp( - M_p(p) - Z_apt(a,p,t-1) );
        }
      }

      // Save recruits in R_pt
      R_pt(p,t) = N_apt(0,p,t);

      // Loop over fleets and compute catch
      for( int f =0; f < nG; f++ )
      {
        C_afpt.col(t).col(p).col(f) = N_apt.col(t).col(p);
        C_afpt.col(t).col(p).col(f) *= (Type(1) - exp( Z_apt.col(t).col(p) ) * F_afpt.col(t).col(p).col(f) / Z_apt.col(t).col(p) );  
        Cw_afpt.col(t).col(p).col(f) = C_afpt.col(t).col(p).col(f) * meanWtAge_ap.col(p);
      }

      // Compute SSB
      B_apt.col(t).col(p) = N_apt.col(t).col(p) * meanWtAge_ap.col(p);
    }

  }

  
  // --------- Observation Model --------- //
  // Age/Length observations and CPUE
  array<Type> I_pft_hat(nP,nG,nT);
  array<Type> aDist_afpt_hat(maxA,nG,nP,nT);
  array<Type> lDist_lfpt_hat(nL,nG,nP,nT);
  I_pft_hat.setZero();
  for( int p = 0; p < nP; p++ )
  {
    int specID = s_p(p);
    int A      = A_s(specID);
    int L      = nLenBins_s(specID);
    for( int f = 0; f < nG; f++ )
    {
      for( int t = 0; t < nT; t++ )
      {
        // Predict CPUE for surveys
        Type expBio = 0.;
        expBio = ( N_apt.col(t).col(p) * meanWtAge_ap.col(p)*sel_afp.col(p).col(f) ).sum();
        I_pft_hat(p,f,t) = q_fp(f,p) * expBio;  

        // Create array of predicted age distributions - this should just be catch-at-age props
        aDist_afpt_hat.col(t).col(p).col(f) = C_afpt.col(t).col(p).col(f) / C_afpt.col(t).col(p).col(f).sum();

        // Create array of predicted length distributions
        // Need to do some thinking here... follow LIME to start with...
        // Probability of being harvested at age
        vector<Type> probHarvAge(A);
        probHarvAge.setZero();
        for( int a = 0; a < A; a ++ )
        {
           probHarvAge(a) = sel_afp(a,f,p) * N_apt (a,p,t) / N_apt.col(t).col(p).sum();
        }

        // Probability of sampling a given length bin
        vector<Type> probHarvLen(L);
        probHarvLen.setZero();
        for( int l = 0; l < L; l++ )
        {
          for( int a = 0; a < A; a++ )
          {
            probHarvLen(l) += probHarvAge(a)*probLenAge_lap(l,a,p);
          }
          lDist_lfpt_hat(l,f,p,t) = probHarvLen(l);
        }
        lDist_lfpt_hat.col(t).col(p).col(f) /= probHarvLen.sum();
      }

    }
  }

  // --------- Statistical Model ----------- //

  // Process model likelihood //
  // Recruitment deviations (initial and ongoing) - add correlation later
  // both auto-correlation and between-pop correlation can be added.
  vector<Type> recnll_p(nP);
  recnll_p.setZero();
  for( int p = 0; p < nP; p++ )
  {
    // First initial recruitment deviations
    int specID = s_p(p);
    int A = A_s(specID);
    for( int a = 0; a < A; A++ )
    {
      recnll_p(p) += Type(0.5) * (lnsigmaR_p(p) + square(omegaRinit_pa(p,a)/sigmaR_p(p)) );
    }
    for( int t = 0; t <  nT; t++ )
    {
      recnll_p(p) += Type(0.5) * (lnsigmaR_p(p) + square(omegaR_pt(p,t)/sigmaR_p(p)) ); 
    }
  }

  // vonB growth model likelihood //
  // Currently has lengths binned the same as the integration above, 
  // we should probably free this up to reduce bias caused by summary stats...
  // (ref review of integrated models (Maunder and Punt, 2013))
  // This could also be replaced with Cadigan's hierarchical
  // vonB model for redfish later (Cadigan 2016)...
  vector<Type> vonBnll_p(nP);
  vonBnll_p.setZero();
  for( int p = 0; p < nP; p++ )
  {
    int specID = s_p(p);
    int A      = A_s(specID);
    for( int a = 0; a < A; a++ )
    {
      for( int l = 0; l < nL; l++ )
      {
        vonBnll_p(p) += ALK_alp(a,l,p) * square( lenBinMids_sl(specID,l) - lenAge_ap(a,p) ) / 2 / square(sigmaL_s(specID) * lenAge_ap(a,p) );
      }
    }
  }

  // Observation likelihoods
  vector<Type> CPUEnll_p(nP);       // CPUE
  vector<Type> ageCompsnll_p(nP);   // Age observations
  vector<Type> lenCompsnll_p(nP);   // length observations
  vector<Type> Ctnll_p(nP);         // Catch
  vector<Type> Dtnll_p(nP);         // Discards
  CPUEnll_p.setZero();
  ageCompsnll_p.setZero();
  lenCompsnll_p.setZero();
  Ctnll_p.setZero();
  Dtnll_p.setZero();
  for( int p = 0; p < nP; p++ )
  {
    int specID = s_p(p);
    int A      = A_s(specID);
    // Loop over fleets
    for( int f = 0; f < nG; f++ )
    {
      // Now loop over time steps
      for( int t = 0; t < nT; t++ )
      {
        // Only use a year if the data exists and that fleet is used as a survey
        if( (I_pft(p,f,t) > 0) & (type_f(f) == 0) )
        {
          Type res = 0.;
          res = log(I_pft(p,f,t)) - log(I_pft_hat(p,f,t));
          CPUEnll_p(p) += Type(0.5) * ( pow(lntau_fp(f,p),2) + pow(res - lnq_fp(f,p), 2 ) / pow(lntau_fp(f,p),2) );
        }

        // Now for compositional data. First, check if ages are being used (or exist), if so
        // compute multinomial likelihood
        vector<Type> currAgeObs = age_afpt.col(t).col(p).col(f);
        vector<Type> predAgeDist = aDist_afpt_hat.col(t).col(p).col(f);
        if( max(currAgeObs) > 0 )
        {
          ageCompsnll_p(p) += dmultinom(currAgeObs, predAgeDist, true );
        }
        // If ages aren't being used, but lengths exist, then compute multinomial likelihood
        vector<Type> currLenObs = length_lfpt.col(t).col(f).col(p);
        vector<Type> predLenDist = lDist_lfpt_hat.col(t).col(p).col(f);
        if( max(currAgeObs) > 0 & max(currLenObs) > 0 )
        {
          lenCompsnll_p += dmultinom(currLenObs, predLenDist , true);
        }
        predAgeDist.setZero();
        predLenDist.setZero();
        // Now compute the catch and discards likelihoods
        Type predCt = Cw_afpt.col(t).col(p).col(f).sum();
        Ctnll_p(p) += Type(0.5)*( square(log(C_pft(p,f,t)) - log(predCt) ) ) / square(tauC_f(f)) ;
      }
    }
  }

  // Robin-Hood Priors //
  // Loop over populations, penalise pop-specific deviations from spec-specific values
  vector<Type> vonBnlp_p(nP);
  vector<Type> qnlp_p(nP);
  vector<Type> steepnessnlp_p(nP);
  vector<Type> Mnlp_p(nP);
  vonBnlp_p.setZero();
  qnlp_p.setZero();
  steepnessnlp_p.setZero();
  Mnlp_p.setZero();
  for( int p = 0; p < nP; p ++ )
  {
    int specID = s_p(p);
    // Growth
    // asymptotic length
    vonBnlp_p(p) += Type(0.5) * ( log(sigmaLinf_s(specID)) + square(Linf_p(p) - muLinf_s(specID) ) / square(sigmaLinf_s (specID)) );
    // growth coefficient
    vonBnlp_p(p) += Type(0.5) * ( log(sigmavonK_s(specID)) + square(vonK_p(p) - muvonK_s(specID) ) / square(sigmavonK_s (specID) ));
    // time at length 0
    vonBnlp_p(p) += Type(0.5) * ( log(sigmat0_s(specID)) + square(vont0_p(p) - mut0_s(specID) ) / square(sigmat0_s (specID) ));
    // Steepness
    Type sigmaSteep_spec = exp(lnsigmaSteep_s(specID));
    steepnessnlp_p(p) += Type(0.5) * ( lnsigmaSteep_s(specID) + square(h_p(p) - h_s(specID)) / square(sigmaSteep_spec) );
    // Natural mortality?
    Mnlp_p(p) += Type(0.5) * ( lnsigmaM_s(specID) + square(M_p(p) - M_s(specID)) / square(sigmaM_s(specID)) );

    // For selectivity and catchability, loop over fleets
    for( int f = 0; f < nG; f++ )
    {      
      // Catchability
      if( type_f(f) == 0) 
        qnlp_p(p) += Type(0.5) * ( lntauq_fs(f,specID) + square( lnq_fp(f,p) - lnqbar_fs(f,specID) ) / exp(2*lntauq_fs(f,specID)) );
    }
  }
  Type pop_nlp = 0.0;
  pop_nlp += vonBnlp_p.sum() + steepnessnlp_p.sum() + Mnlp_p.sum() + qnlp_p.sum();

  // Now loop over species and penalise by "complex level" means - will have to think
  // carefully about what parameters have a complex mean, and what don't
  // vector<Type> vonBnlp_p(nP);
  // vector<Type> selnlp_p(nP);
  vector<Type> qnlp_s(nS);
  vector<Type> selnlp_s(nS);
  vector<Type> steepnessnlp_s(nS);
  vector<Type> Mnlp_s(nS);
  // vonBnlp_p.setZero();
  selnlp_s.setZero();
  qnlp_s.setZero();
  steepnessnlp_s.setZero();
  Mnlp_s.setZero();
  for( int s = 0; s < nS; s++ )
  {
    // Steepness
    Type sigmaSteep = exp(lnsigmaSteep);
    steepnessnlp_s(s) += Type(0.5) * ( lnsigmaSteep + square(h_s(s) - muSteep) / square(sigmaSteep) );

    // Natural mortality
    Mnlp_s(s) += Type(0.5) * ( lnsigmaM + square(M_s(s) - muM) / square(sigmaM) );

    // Catchability
    for( int f = 0; f < nG; f++ )
    {
      // Catchability
      if( type_f(f) == 0) 
        qnlp_s(s) += Type(0.5) * ( lntauq_f(f) + square( lnqbar_fs(f,s) - lnqbar_f(f) ) / exp(2*lntauq_f(f)) );

      // Selectivity
      selnlp_s(s) += Type(0.5) * ( log(sigmalenSel50_f(f) ) + square( mulenSel50_f(f) - lenSel50_sf(s,f) ) / square(sigmalenSel50_f(f)) );
      selnlp_s(s) += Type(0.5) * ( log(sigmalenSel95_f(f) ) + square( mulenSel95_f(f) - lenSel95_sf(s,f) ) / square(sigmalenSel95_f(f)) );
    }
  }
  Type spec_nlp = 0.0;  
  spec_nlp += steepnessnlp_s.sum() + Mnlp_s.sum() + qnlp_s.sum();
  

  // Now take the sum of all the likelihoods, keep separate
  // in case we want to weight data later
  Type joint_nlp = 0.0;
  // Observations
  joint_nlp += Ctnll_p.sum();       // Catch
  joint_nlp += lenCompsnll_p.sum(); // Length compositions
  joint_nlp += ageCompsnll_p.sum(); // Age compositions
  joint_nlp += CPUEnll_p.sum();     // Survey CPUE
  // Growth model
  joint_nlp += vonBnll_p.sum();     // Growth model
  // Recruitment errors
  joint_nlp += recnll_p.sum();      // recruitment process errors
  joint_nlp += pop_nlp + spec_nlp;

  
  // Return quantities
  // Leading parameters
  REPORT(B0_p);    
  REPORT(h_p);     
  REPORT(M_p);     
  REPORT(Rinit_p); 
  REPORT(Linf_p);  
  REPORT(vonK_p);  
  REPORT(vont0_p); 
  REPORT(sigmaL_s);
  REPORT(sigmaR_p);
  // Fishery model pars
  REPORT(q_fp)          // Catchability
  REPORT(tau_fp)        // Observation error
  REPORT(lenSel50_sf)   // length at 50% sel by species
  REPORT(lenSel95_sf)   // length at 95% sel by specie
  REPORT(F_pft)         // Fishing mortality
  // Model states
  REPORT(B_apt);
  REPORT(N_apt);
  REPORT(R_pt);
  REPORT(F_afpt);
  REPORT(Z_apt);
  REPORT(C_afpt);
  REPORT(Cw_afpt);
  // Stock recruit parameters
  REPORT(R0_p);
  REPORT(phi_p);
  REPORT(reca_p);
  REPORT(recb_p);

  // Likelihood values
  REPORT(joint_nlp);
  REPORT(recnll_p);
  REPORT(vonBnll_p);
  REPORT(CPUEnll_p);
  REPORT(ageCompsnll_p);
  REPORT(Ctnll_p);

  // sd report for standard errors
  ADREPORT(B0_p);    
  ADREPORT(h_p);     
  ADREPORT(M_p);     
  ADREPORT(Rinit_p); 
  ADREPORT(Linf_p);  
  ADREPORT(vonK_p);  
  ADREPORT(vont0_p); 
  ADREPORT(sigmaL_s);
  ADREPORT(sigmaR_p);
  ADREPORT(q_fp)          // Catchability
  ADREPORT(tau_fp)        // Observation error
  ADREPORT(lenSel50_sf)   // length at 50% sel
  ADREPORT(lenSel95_sf)   // length at 95% sel
  ADREPORT(F_pft)         // Fishing mortality

  return( joint_nlp );

}

