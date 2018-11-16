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
  DATA_ARRAY(I_spft);               // CPUE data
  DATA_ARRAY(C_spft);               // Catch data (biomass)
  DATA_ARRAY(D_spft);               // Discard data (biomass)
  // DATA_ARRAY(ALK_spal);             // Age-length observations (in freq) by pop - time averaged
  DATA_ARRAY(age_aspft);            // Age observations for each population and fleet over time (-1 missing)
  DATA_ARRAY(len_lspft);            // Length observations for each population and fleet over time (-1 missing)
  DATA_IVECTOR(type_f);             // Fleet type (0 survey, 1 fishery dependent)
  DATA_IVECTOR(A_s);                // +Group age by species
  DATA_IVECTOR(L_s);                // +Group length by species
  DATA_IVECTOR(lenD_s);             // Discard length by species
  // DATA_MATRIX(intMatrix_pp);        // population interactions mtx (spawning?, adjacency)
  // DATA_IVECTOR(nLenBins_s);         // Number of length bins by species

  // Model dimensions
  int nS = I_spft.dim(0);           // No. of species
  int nP = I_spft.dim(1);           // No. of stocks in each species
  int nF = I_spft.dim(2);           // No. of fleets (surveys + trawl)  
  int nT = I_spft.dim(3);           // No of time steps
  int nL = len_lspft.dim(0);        // max no of length bins (for creating state arrays)
  int nA = age_aspft.dim(0);        // max no of age bins (for creating state arrays)

  // Model switches - fill these in when we know what we want
  DATA_MATRIX(swRinit_sp);          // pop-spec fished initialisation switch (0 == unfished, 1 == fished)
  DATA_INTEGER(parSwitch);          // parallel accumulator switch for faster run time

  /*parameter section*/
  // Leading Parameters
  // Biological
  PARAMETER_ARRAY(lnB0_sp);         // Biomass at unfished for ind. stocks
  PARAMETER_ARRAY(logitSteep_sp);   // pop specific steepness
  PARAMETER_ARRAY(lnM_sp);          // stock-specific M
  PARAMETER_ARRAY(lnRinit_sp);      // Pop-specific initial R (non-eq init)
  PARAMETER_ARRAY(lnLinf_sp);       // pop-specific vonB Linf parameter
  PARAMETER_ARRAY(lnvonK_sp);       // pop-specific vonB K parameter
  PARAMETER_ARRAY(lnL1_sp);         // pop-specific vonB Length-at-age 1
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
  PARAMETER_ARRAY(lnlenSel50_sf);   // Selectivity Alpha parameter by species/fleet
  PARAMETER_ARRAY(lnlenSel95_sf);   // Selectivity Beta parameter by species/fleet
  // Fishing mortality
  PARAMETER_ARRAY(lnF_spft);        // Fishing mortality by species/population/fleet/time
  PARAMETER_VECTOR(lntauC_f);       // Catch observation SD by fleet
  PARAMETER_VECTOR(lntauD_f);       // Discards observation SD by fleet

  // Priors
  // Growth //
  PARAMETER_VECTOR(muLinf_s);       // species mean Linf
  PARAMETER_VECTOR(sigmaLinf_s);    // species Linf sd    
  PARAMETER_VECTOR(muvonK_s);       // species mean vonK  
  PARAMETER_VECTOR(sigmavonK_s);    // species vonK sd       
  PARAMETER_VECTOR(muL1_s);         // species mean t0  
  PARAMETER_VECTOR(sigmaL1_s);      // species t0 sd       
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
  PARAMETER_ARRAY(omegaR_spt);      // species-stock specific recruitment errors 2:nT
  PARAMETER_ARRAY(omegaRinit_spa);  // stock-age specific recruitment initialisation errors
  PARAMETER_ARRAY(lnsigmaR_sp);    // stock recruitment errors sd (sqrt cov matrix diag)
  PARAMETER_VECTOR(logitRCorr_chol);// stock recruitment errors corr chol factor off diag  
  PARAMETER_ARRAY(logitRgamma_sp); // species-stock-specific AR1 auto-corr on year effect (omega)
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
  array<Type>  Rinit_sp(nS,nP);
  array<Type>  Linf_sp(nS,nP);
  array<Type>  vonK_sp(nS,nP);
  array<Type>  L1_sp(nS,nP);
  array<Type>  sigmaR_sp(nS,nP);

  for( int sIdx = 0; sIdx < nS; sIdx++ )
  {
    for( int pIdx = 0; pIdx < nP; pIdx++ )
    {
      B0_sp(sIdx,pIdx)      = exp(lnB0_sp(sIdx,pIdx));
      h_sp(sIdx,pIdx)       = .2 + Type(0.8) / ( Type(1.0) + exp( -logitSteep_sp(sIdx,pIdx)) );
      M_sp(sIdx,pIdx)       = exp(lnM_sp(sIdx,pIdx));
      Rinit_sp(sIdx,pIdx)   = exp(lnRinit_sp(sIdx,pIdx));
      Linf_sp(sIdx,pIdx)    = exp(lnLinf_sp(sIdx,pIdx));
      vonK_sp(sIdx,pIdx)    = exp(lnvonK_sp(sIdx,pIdx));
      L1_sp(sIdx,pIdx)      = exp(lnL1_sp(sIdx,pIdx));
      sigmaR_sp(sIdx,pIdx)  = exp(lnsigmaR_sp(sIdx,pIdx));
    }
  }

  vector<Type>  sigmaL_s       = exp(lnsigmaL_s);

  // Observation model
  array<Type>   q_spf(nS,nP,nF);
  array<Type>   tau_spf(nS,nP,nF);
  // selectivity by fleet //
  array<Type>  lenSel50_sf(nS,nF);
  array<Type>  lenSel95_sf(nS,nF);
  q_spf.setZero();
  tau_spf.setZero();
  for( int f = 0; f < nF; f++ )
    for( int s = 0; s < nS; s++ )
    {
      lenSel50_sf(s,f)   = exp( lnlenSel50_sf(s,f));
      lenSel95_sf(s,f)   = exp( lnlenSel95_sf(s,f));
      for( int p = 0; p < nP; p++ )
      {
        q_spf(s,p,f)    = exp( lnq_spf(s,p,f) );
        tau_spf(s,p,f)  = exp( lntau_spf(s,p,f) );    
      }
    }
  
  // Fishing mortality
  array<Type>   F_spft(nS,nP,nF,nT);
  F_spft.setZero();
  for( int s = 0; s < nS; s++ )
    for( int f = 0; f < nF; f++ )
      for( int p = 0; p < nP; p++ )
        for( int t = 0; t < nT; t++ )
          F_spft(s,p,f,t) = exp( lnF_spft(s,p,f,t));
  
  // Probably concentrate these out later....
  vector<Type>  tauC_f        = exp(lntauC_f);    
  vector<Type>  tauD_f        = exp(lntauD_f);    
  
  // Set up model state arrays
  array<Type> B_aspt(nA,nS,nP,nT);      // Biomass at age-pop-time
  array<Type> N_aspt(nA,nS,nP,nT);      // Numbers at age-pop-time
  array<Type> R_spt(nS,nP,nT);          // Recruits by pop-time
  array<Type> F_aspft(nA,nS,nP,nF,nT);  // Fishing mortality by age, fleet, population and time
  array<Type> Z_aspt(nA,nS,nP,nT);      // Total mortality by age, population and time
  array<Type> C_aspft(nA,nS,nP,nF,nT);  // Predicted catch-at-age (in numbers), fleet, population and time
  array<Type> Cw_aspft(nA,nS,nP,nF,nT); // Predicted catch-at-age (in weight), population and time
  array<Type> predCw_spft(nS,nP,nF,nT); // Predicted catch (in weight) by population and time
  array<Type> predC_spft(nS,nP,nF,nT);  // Predicted catch (in numbers) by population and time
  array<Type> B_spt(nS,nP,nT);          // Total biomass by species, pop, time
  array<Type> SB_spt(nS,nP,nT);         // Spawning biomass by species, pop, time
  array<Type> Bv_spft(nS,nP,nF,nT);     // Vulnerable biomass by species, pop, time


  // derived variables
  // Stock recruitment //
  array<Type>  R0_sp(nS,nP);            // eqbm recruitment
  array<Type>  phi_sp(nS,nP);           // eqbm SSB per recruit
  array<Type>  reca_sp(nS,nP);          // BH a parameter for pops
  array<Type>  recb_sp(nS,nP);          // BH b parameter for pops

  // Growth //
  array<Type>   Wlen_ls(nL,nS);         // weight-at-length by species
  array<Type>   lenAge_asp(nA,nS,nP);   // Mean length-at-age by population

  // Maturity //
  array<Type>   matAge_asp(nA,nS,nP);   // Proportion mature at age by population
  array<Type>   matLen_ls(nL,nS);       // Proportion mature at Length by species

  // Prior Hyperparameters //
  // Steepness
  vector<Type>  h_s           = .2 + Type(0.8) / ( Type(1.0) + exp( -logitSteep_s) ) ;
  Type          muSteep       = .2 + Type(0.8) / ( Type(1.0) + exp( -logitSteep) ) ;
  // Natural mortality
  vector<Type>  M_s           = exp(lnM_s);
  vector<Type>  sigmaM_s      = exp(lnsigmaM_s);
  Type          muM           = exp(ln_muM);
  Type          sigmaM        = exp(lnsigmaM);


  // ---------------- Growth, Maturity, Stock-Recruitment ----------- //
  Wlen_ls.fill(-1.0);
  matLen_ls.fill(-1.0);
  for( int s = 0; s < nS; s++ )
    for( int l = 0; l < L_s(s); l++ )
    {
      Type len = l+1;
      // weight-at-length for each length bin midpoint
      Wlen_ls(l,s) = LWa_s(s) * pow( l+1, LWb_s(s) );
      // proportion mature-at-length
      matLen_ls(l,s) = 1; 
      matLen_ls(l,s) /= (1 + exp( -1. * log(Type(19.0)) * (len - xMat50_s(s)) / (xMat95_s(s) - xMat50_s(s)) ) );
    }

  // Now produce the probability of length-at-age matrix for each stock
  // (taken from LIME by Rudd and Thorson, 2017)
  // In the same loop, we're going to compute the ssbpr parameter phi
  // for each stock
  array<Type> probLenAge_lasp(nL,nA,nS,nP);
  array<Type> meanWtAge_asp(nA,nS,nP);
  lenAge_asp.setZero();
  probLenAge_lasp.setZero();
  meanWtAge_asp.setZero();
  phi_sp.setZero();
  matAge_asp.setZero();

  for( int s = 0; s < nS; s++)
    for( int p = 0; p < nP; p++ )
    {
      // Calc Linf CV
      Type  lenAgeCV  = sigmaL_s( s ) / Linf_sp ( s, p );
      for( int a = 0; a < A_s(s); a ++ )
      {
        Type  sumProbs  = 0;
        lenAge_asp(a,s,p) = Linf_sp(s,p) + (L1_sp(s,p) - Linf_sp(s,p)) * exp( -1. * vonK_sp( s, p ) * ( a ) ) ; 
        for( int l = 0; l < L_s(s); l++ )
        {
          Type len = l + 1;
          Type lenHi = len + .5;
          Type lenLo = len - .5;
          if( l == 0 )
          {
            // Compute LH tail of the distribution
            probLenAge_lasp(l,a,s,p) = pnorm( lenHi, lenAge_asp(a,s,p), lenAge_asp(a,s,p) * lenAgeCV );
            sumProbs += probLenAge_lasp(l,a,s,p);
          }
          if( (l >= 1) & (l < L_s(s) - 1 ) )
          {
            // Compute interior bins
            probLenAge_lasp(l,a,s,p) = pnorm( lenHi, lenAge_asp(a,s,p), lenAge_asp(a,s,p) * lenAgeCV );
            probLenAge_lasp(l,a,s,p) -= pnorm( lenLo, lenAge_asp(a,s,p), lenAge_asp(a,s,p) * lenAgeCV );
            sumProbs += probLenAge_lasp(l,a,s,p);
          }
          // Now the RH tail
          if( l == L_s(s) - 1 ) 
            probLenAge_lasp(l,a,s,p) = Type(1.0) - sumProbs;

          // Calculate mean weight at age from probLenAge
          meanWtAge_asp(a,s,p) += probLenAge_lasp(l,a,s,p) * Wlen_ls(l,s);
          // Calculate maturity at age from probLenAge
          if( (matLen_ls(l,s) > 0) )
            matAge_asp(a,s,p) += probLenAge_lasp(l,a,s,p) * matLen_ls(l,s);
        }

        // To compute ssbpr, we need to borrow N_apt for a moment 
        Type ssbpr_a = exp(-a * M_sp(s,p)) * matAge_asp(a,s,p) * meanWtAge_asp(a,s,p);
        if ( a == A_s(s) - 1 ) ssbpr_a /= (1. - exp( -M_sp(s,p)) );
        phi_sp(s,p) += ssbpr_a;
      }
    }



  // --------- Selectivity and Fishing Mortality -------- //
  array<Type> sel_lfsp(nL,nF,nS,nP);
  array<Type> sel_afsp(nA,nF,nS,nP);
  sel_lfsp.setZero();
  sel_afsp.setZero();
  for( int s = 0; s < nS; s++ )
    for( int p = 0; p < nP; p++ )
    {
      // Loop over fleets, create selectivities
      for( int f = 0; f < nF; f++ )
      {
        // selectivity-at-length
        for( int l = 0; l < L_s(s); l++ )
        {
          sel_lfsp(l,f,s,p) = 1;
          sel_lfsp(l,f,s,p) /= ( 1 + exp( - log(Type(19.)) * ( l+1 - lenSel50_sf(s,f) ) / ( lenSel95_sf(s,f) - lenSel50_sf(s,f) ) ) );
        }
        // convert selectivity-at-length to selectivity-at-age using probability matrix
        for( int a = 0; a < A_s(s); a++ )
        {
          vector<Type> probLenAgea_l(nL) ;
          probLenAgea_l.setZero();
          probLenAgea_l = probLenAge_lasp.col(p).col(s).col(a);
          sel_afsp(a,f,s,p) = (probLenAgea_l*sel_lfsp.col(p).col(s).col(f)).sum();
        }
      }
    }

  // --------- Population Dynamics ----------- //
  // First, compute eqbm recruitment
  // Derived parameters
  for( int s = 0; s < nS; s++ )
    for( int p = 0; p < nP; p++ )
      R0_sp(s,p) = B0_sp(s,p) / phi_sp(s,p);

  // Set all state arrays to zero
  B_aspt.setZero();
  N_aspt.setZero();
  R_spt.setZero();
  F_aspft.setZero();
  Z_aspt.setZero();
  C_aspft.setZero();
  Cw_aspft.setZero();
  B_spt.setZero();
  SB_spt.setZero();
  Bv_spft.setZero();

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
          Type InitRtmp;
          if( swRinit_sp(s,p) == 0 ) InitRtmp = R0_sp(s,p);

          // Populate first year
          for( int a = 0; a < A; a++ )
          {
            // Produce recruits
            N_aspt(a,s,p,t) = InitRtmp * exp(-Type(a) * M_sp(s,p) );
            // Non-eqbm initialisation
            if(swRinit_sp(s,p) == 1)
              N_aspt(a,s,p,t) *= exp( omegaRinit_spa(s,p,a) - square(sigmaR_sp(s,p)) / 2 );
            // Compute total mortality at age (for catch later)
            Z_aspt(a,s,p,t) = M_sp(s,p);
            for( int f = 0; f < nF; f++ )
            {
              // Compute fishing mortality at age
              F_aspft(a,s,p,f,t) = sel_afsp(a,f,s,p) * F_spft(s,p,f,t);
              // Add to Z
              Z_aspt(a,s,p,t) += F_aspft(a,s,p,f,t);  
            }
          }
          // Correct plus-group numbers at eqbm if initialised at unfished
          if( swRinit_sp(s,p) == 0 ) 
            N_aspt(A-1,s,p,t) /= (Type(1) - exp(-M_sp(s,p)));
        }

        // time series history
        if( t > 0 )
        {
          // Generate recruitment
          Type SBt = (B_aspt.col(t-1).col(p).col(s) * matAge_asp.col(p).col(s) ).sum();
          N_aspt(0,s,p,t) = Type(4) * R0_sp(s,p) * SBt / ( B0_sp(s,p) * (Type(1) - h_sp(s,p) ) + (Type(5)*h_sp(s,p) - Type(1)) * SBt );
          N_aspt(0,s,p,t) *= exp( omegaR_spt(s,p,t-1) - square(sigmaR_sp(s,p)) / 2 );

          // Now loop over ages and apply fishing mortality (no discarding yet)
          for( int a = 0; a < A; a ++ )
          {
            // Compute fishing mortality at age
            // Compute total mortality at age
            Z_aspt(a,s,p,t) = M_sp(s,p);
            for( int f = 0; f < nF; f++ ) 
            {
              F_aspft(a,s,p,f,t) = sel_afsp(a,f,s,p) * F_spft(s,p,f,t);
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
          // Have to loop over ages here to avoid adding NaNs
          for(int a = 0; a < A; a++ )
          {
            C_aspft(a,s,p,f,t) = N_aspt(a,s,p,t);
            C_aspft(a,s,p,f,t) *= (Type(1) - exp( Z_aspt(a,s,p,t) ) * F_aspft(a,s,p,f,t) / Z_aspt(a,s,p,t) );  
            Cw_aspft(a,s,p,f,t) = C_aspft(a,s,p,f,t) * meanWtAge_asp(a,s,p);
            predCw_spft(s,p,f,t) += Cw_aspft(a,s,p,f,t);
            predC_spft(s,p,f,t) += C_aspft(a,s,p,f,t);

            // Calculate vulnerable biomass for this fleet
            Bv_spft(s,p,f,t) += N_aspt(a,s,p,t) * sel_afsp(a,f,s,p) * meanWtAge_asp(a,s,p);
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
  I_spft_hat.setZero();
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
          Type expBio = 0.;
          expBio = ( N_aspt.col(t).col(p).col(s) * meanWtAge_asp.col(p).col(s)*sel_afsp.col(p).col(s).col(f) ).sum();
          I_spft_hat(s,p,f,t) = q_spf(s,p,f) * expBio;  

          // Create array of predicted age distributions - this should just be catch-at-age props
          // Calculate total catch
          Type totCatch = 0.;
          for( int a = 0; a < A; a ++ )
            totCatch += C_aspft(a,s,p,f,t);
          // Convert catch-at-age to proportions-at-age
          for( int a = 0; a < A; a ++ )
            aDist_aspft_hat(a,s,p,f,t) = C_aspft(a,s,p,f,t) / totCatch;

          // Create array of predicted length distributions
          // Need to do some thinking here... follow LIME to start with...
          // Probability of being harvested at age
          vector<Type> probHarvAge(A);
          probHarvAge.setZero();
          for( int a = 0; a < A; a ++ )
          {
             probHarvAge(a) = sel_afsp(a,f,s,p) * N_aspt(a,s,p,t) / N_aspt.col(t).col(p).col(s).sum();
          }

          // Probability of sampling a given length bin
          vector<Type> probHarvLen(L);
          probHarvLen.setZero();
          // loop over length and age, calculate probability of harvest at length
          for( int l = 0; l < L; l++ )
          {
            for( int a = 0; a < A; a++ )
              probHarvLen(l) += probHarvAge(a)*probLenAge_lasp(l,a,s,p);

            // Save to length dists
            lDist_lspft_hat(l,s,p,f,t) = probHarvLen(l);
          }
          // renormalise
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
      // for( int a = 0; a < A; A++ )
      //   recnll_sp(s,p) += Type(0.5) * (2 * lnsigmaR_sp(s,p) + square(omegaRinit_spa(s,p,a)/sigmaR_sp(s,p)) );

      // then yearly recruitment deviations
      for( int t = 1; t <  nT; t++ )
        recnll_sp(s,p) += Type(0.5) * (2 * lnsigmaR_sp(s,p) + square(omegaR_spt(s,p,t-1)/sigmaR_sp(s,p)) ); 
    }

  // // vonB growth model likelihood //
  // // Currently has lengths binned the same as the integration above, 
  // // we should probably free this up to reduce bias caused by summary stats...
  // // (ref review of integrated models (Maunder and Punt, 2013))
  // // This could also be replaced with Cadigan's hierarchical
  // // vonB model for redfish later (Cadigan 2016)...
  // array<Type> vonBnll_sp(nS,nP);
  // vonBnll_sp.setZero();
  // for( int s = 0; s < nS; s++ )
  //   for( int p = 0; p < nP; p++ )
  //   {
  //     Type cvL = sigmaL_s(s) / Linf_sp(s,p);
  //     int A      = A_s(s);
  //     for( int a = 0; a < A; a++ )
  //       for( int l = 0; l < nL; l++ )
  //         vonBnll_sp(s,p) += ALK_spal(s,p,a,l) * 0.5 * ( 2 * log(cvL * lenAge_asp(a,s,p)) + square( (l+1 - lenAge_asp(a,s,p)) / (cvL * lenAge_asp(a,s,p) ) ) );
  //   }

  // Observation likelihoods
  array<Type> CPUEnll_sp(nS,nP);       // CPUE
  array<Type> ageCompsnll_sp(nS,nP);   // Age observations
  array<Type> lenCompsnll_sp(nS,nP);   // length observations
  array<Type> Ctnll_sp(nS,nP);         // Catch
  array<Type> Dtnll_sp(nS,nP);         // Discards
  CPUEnll_sp.setZero();
  ageCompsnll_sp.setZero();
  lenCompsnll_sp.setZero();
  Ctnll_sp.setZero();
  Dtnll_sp.setZero();
  for( int s = 0; s < nS; s++ )
    for( int p = 0; p < nP; p++ )
    {
      int A      = A_s(s);
      int L      = L_s(s);
      // Loop over fleets
      for( int f = 0; f < nF; f++ )
      {
        // tmp vectors to hold age and length observed and
        // predicted values
        vector<Type> fleetAgeObs(A);
        vector<Type> fleetAgePred(A);

        vector<Type> fleetLenObs(L);
        vector<Type> fleetLenPred(L);     

        // Set tmp vectors to zero 
        fleetAgeObs.setZero();
        fleetAgePred.setZero();
        fleetLenObs.setZero();
        fleetLenPred.setZero();

        // Now loop over time steps
        for( int t = 0; t < nT; t++ )
        {
          // Only use a year if the data exists and that fleet is used as a survey
          if( (I_spft(s,p,f,t) > 0) & (type_f(f) == 0) )
          {
            Type res = 0.;
            res = log(I_spft(s,p,f,t)) - log(I_spft_hat(s,p,f,t));
            CPUEnll_sp(s,p) += Type(0.5) * ( pow(lntau_spf(s,p,f),2) + pow(res - lnq_spf(s,p,f), 2 ) / pow(lntau_spf(s,p,f),2) );
          }

          // Now compute the catch and discards likelihoods
          if( C_spft(s,p,f,t) > 0 & predCw_spft(s,p,f,t) > 0 )
            Ctnll_sp(s,p) += 0.5* ( 2 * log(tauC_f(f)) +  square( (log(C_spft(s,p,f,t)) - log(predCw_spft(s,p,f,t)) ) / tauC_f(f)) ) ;

          // Loop and fill age obs and pred matrices
          for( int a = 0; a < A; a++ )
          {
            fleetAgeObs(a)  = age_aspft(a,s,p,f,t);
            fleetAgePred(a) = aDist_aspft_hat(a,s,p,f,t);
          }

          // Loop and fill length obs and pred matrices
          for( int l = 0; l < L; l++ )
          {
            fleetLenObs(l)  = len_lspft(l,s,p,f,t);
            fleetLenPred(l) = lDist_lspft_hat(l,s,p,f,t);
          }

          // Now for compositional data. First, check if ages are being used (or exist), if so
          // compute multinomial likelihood
          if( fleetAgeObs.sum() > 0 )
          {
            // Code a general comps likelhood function
            // to use here, and add a switch to the
            // function inputs to choose comps model
            ageCompsnll_sp(s,p) -= dmultinom(fleetAgeObs, fleetAgePred, true );
          }     

          // If ages aren't being used, but lengths exist, then compute multinomial likelihood
          if( fleetAgeObs.sum() <= 0 & fleetLenObs.sum() > 0 )
          {
            lenCompsnll_sp(s,p) -= dmultinom(fleetLenObs, fleetLenPred, true);
          }     
          
        }
        // Set tmp vectors to zero 
        fleetAgeObs.setZero();
        fleetAgePred.setZero();
        fleetLenObs.setZero();
        fleetLenPred.setZero();

      }
    }

  // Shared Hierarchical Prior Distributions //
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
  //   vonBnlp_p(p) += Type(0.5) * ( log(sigmaLinf_s(specID)) + square(Linf_p(p) - muLinf_s(specID) ) / square(sigmaLinf_s (specID)) );
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
  //     if( type_f(f) == 0) 
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
  //     if( type_f(f) == 0) 
  //       qnlp_s(s) += Type(0.5) * ( lntauq_f(f) + square( lnqbar_fs(f,s) - lnqbar_f(f) ) / exp(2*lntauq_f(f)) );

  //     // Selectivity
  //     selnlp_s(s) += Type(0.5) * ( log(sigmalenSel50_f(f) ) + square( mulenSel50_f(f) - lenSel50_sf(s,f) ) / square(sigmalenSel50_f(f)) );
  //     selnlp_s(s) += Type(0.5) * ( log(sigmalenSel95_f(f) ) + square( mulenSel95_f(f) - lenSel95_sf(s,f) ) / square(sigmalenSel95_f(f)) );
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

  
  // Observations
  joint_nlp += Ctnll_sp.sum();       // Catch
  joint_nlp += lenCompsnll_sp.sum(); // Length compositions
  joint_nlp += ageCompsnll_sp.sum(); // Age compositions
  joint_nlp += CPUEnll_sp.sum();     // Survey CPUE
  // Growth model
  // joint_nlp += vonBnll_sp.sum();     // Growth model
  // Recruitment errors
  joint_nlp += recnll_sp.sum();      // recruitment process errors
  // joint_nlp += pop_nlp + spec_nlp;

  
  // // Return quantities
  // Model dimensions
  REPORT(nS);
  REPORT(nP);
  REPORT(nF);
  REPORT(nT);
  REPORT(nL);
  REPORT(nA);

  // Leading parameters
  REPORT(B0_sp);    
  REPORT(h_sp);     
  REPORT(M_sp);     
  REPORT(Rinit_sp); 
  REPORT(Linf_sp);  
  REPORT(vonK_sp);  
  REPORT(L1_sp); 
  REPORT(sigmaL_s);
  REPORT(sigmaR_sp);
  REPORT(tauC_f);
  REPORT(tauD_f);
  // Fishery model pars
  REPORT(q_spf)           // Catchability
  REPORT(tau_spf)         // Observation error
  REPORT(lenSel50_sf)     // length at 50% sel by species
  REPORT(lenSel95_sf)     // length at 95% sel by specie
  REPORT(F_spft)          // Fishing mortality
  REPORT(sel_lfsp);
  REPORT(sel_afsp);
  // Model states
  REPORT(B_aspt);
  REPORT(N_aspt);
  REPORT(R_spt);
  REPORT(F_aspft);
  REPORT(Z_aspt);
  REPORT(C_aspft);
  REPORT(Cw_aspft);
  REPORT(predCw_spft);
  REPORT(predC_spft);
  REPORT(Bv_spft);
  REPORT(SB_spt);
  REPORT(B_spt);
  // Stock recruit parameters
  REPORT(R0_sp);
  REPORT(phi_sp);
  REPORT(reca_sp);
  REPORT(recb_sp);

  // Growth and maturity
  REPORT( probLenAge_lasp );
  REPORT( lenAge_asp );
  REPORT( Wlen_ls );
  REPORT( matLen_ls );
  REPORT( matAge_asp );

  // Likelihood values
  REPORT(joint_nlp);
  REPORT(recnll_sp);
  // REPORT(vonBnll_sp);
  REPORT(CPUEnll_sp);
  REPORT(ageCompsnll_sp);
  REPORT(lenCompsnll_sp);
  REPORT(Ctnll_sp);

  // Data
  REPORT( I_spft );
  REPORT( C_spft );
  REPORT( D_spft );
  // REPORT( ALK_spal );
  REPORT( age_aspft );
  REPORT( len_lspft );
  REPORT( type_f );
  REPORT( A_s );
  REPORT( L_s );
  REPORT( lenD_s );

  // Simulated observations
  REPORT( lDist_lspft_hat );
  REPORT( aDist_aspft_hat );
  REPORT( I_spft_hat );


  // sd report for standard errors
  ADREPORT(lnB0_sp);    
  // ADREPORT(logitSteep_sp);     
  // ADREPORT(lnM_sp);     
  // ADREPORT(lnRinit_sp); 
  // ADREPORT(lnLinf_sp);  
  // ADREPORT(lnvonK_sp);  
  // ADREPORT(lnL1_sp); 
  // ADREPORT(lnsigmaL_s);
  // ADREPORT(lnsigmaR_sp);
  // ADREPORT(lnq_spf)          // Catchability
  // ADREPORT(lntau_spf)        // Observation error
  // ADREPORT(lnlenSel50_sf)   // length at 50% sel
  // ADREPORT(lnlenSel95_sf)   // length at 95% sel
  // ADREPORT(lnF_spft)         // Fishing mortality


  return( joint_nlp );

}

