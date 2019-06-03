The posterior distributions of the parameters from the multi-state model go in this directory.
The following notation is used to denote model type:

* fp-p0 (1st order fractional polynomial model with p1 = 0 (i.e., a Weibull model))
* fp-p00 (2nd order fractional polynomial model with p1 = p2 = 0)
* fp-p1 order fractional polynomial model with p1 = 1 (i.e., a Gompertz model))
* fp-p01 (2nd order fractional polynomial model with p1 = 0 and p2 = 1)

The directory must contain files for the meta-analyses (models for absolute effects) 
and network meta-analyses (models for relative treatment effects).

# Meta-analysis files

First line fixed effects meta-analyses for gefitnib:

* ma-1L-fe-gef-fp-p0 
* ma-1L-fe-gef-fp-p00 
* ma-1L-fe-gef-fp-p1 
* ma-1L-fe-gef-fp-p01 

First line fixed effects meta-analyses for platinum-based doublet chemotherapy (PBDC):

* ma-1L-fe-pbdc-fp-p0 
* ma-1L-fe-pbdc-fp-p00 
* ma-1L-fe-pbdc-fp-p1 
* ma-1L-fe-pbdc-fp-p01 

Second line fixed effects meta-analysis for osimertinib among patients with a T790M mutation:

* ma-2L-fe-t790m-osi-fp-p0
* ma-2L-fe-t790m-osi-fp-p00 
* ma-2L-fe-t790m-osi-fp-p1 
* ma-2L-fe-t790m-osi-fp-p01 

# Network meta-analysis files

* nma-1L-fe-fp-p0
* nma-1L-fe-fp-p00
* nma-1L-fe-fp-p1
* nma-1L-fe-fp-p01
