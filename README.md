# SCAL-RH

Statistical catch-at-length multi-species and multi-stock assessment model. Defined using shared hierarchical priors on biological and observation model parameters, with an integrated growth model to enable use of age and length composition data.


## ToDo

### Biological Data summary
Pull out and summarise CPUE, catch, and biological data by *major stat* area for each DERPA species (or make scale of data/assessment a parameter for the functions to do so). 

This will enable identification of whether a species exists in sufficient quantity in a stat area for assesment. - this will also require some post stratification of the synoptic survey, maybe this is a bad idea, and we manage at the level of the survey on this round (combine HS and WCHG). 

Create some visuals of age comps filled in with transformed length comps.

Arrange length compositions and age compositions in arrays for scalRH to selectivity observation models

Create time-averaged array of age/length observations for integrated growth model

Test whether maturity is length-based or age-based. Perhaps include a switch for this in scalRH

### Discarding
Is there a grading length for DERPA species? How strict is it? 

Discard induced mortality? Start with dTrawl = 0.8 (exploitation, not inst.)

Can probably get some grading lengths from the comm length dists, at least in the years that discarded fish were sampled.

### Model Priors
Define observation model and biological shared priors. Candidates are:
  a. spatial/multispecies:
    i. process error correlations
    ii. fishery dependent values (e.g. lenSel50 )
    iii. If not conditioning on catch, F values could be correlated
    iv. catchability within a survey/area
    v. natural mortality (M) deviations within an area - spatial effects? Time varying?
  b. intraspecies/across areas: 
    i. vonB parameters between areas within a species (sexually dimorphic growth?)
    ii. initial natural mortality (M0) within a species - basic demographic rate
