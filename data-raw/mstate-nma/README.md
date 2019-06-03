The posterior distributions of the parameters from the multi-state model go in this directory. 
The file `../params_mstate_nma.R` uses this data to create the `params_mstate_nma` 
object that loads with the package.

The following notation is used to denote model type:

* `fp-p0` (1st order fractional polynomial model with p1 = 0 (i.e., a Weibull model))
* `fp-p00` (2nd order fractional polynomial model with p1 = p2 = 0)
* `fp-p1` (1st order fractional polynomial model with p1 = 1 (i.e., a Gompertz model))
* `fp-p01` (2nd order fractional polynomial model with p1 = 0 and p2 = 1)

The directory contains two types of files: (1) the parameter estimates and 
(2) lookup files.

# Parameter estimates
The directory must contain files for the meta-analyses (models for absolute effects) 
and network meta-analyses (models for relative treatment effects).

## Meta-analysis files
First line fixed effects meta-analyses for gefitnib:

* `ma-1L-fe-gef-fp-p0.csv` 
* `ma-1L-fe-gef-fp-p00.csv` 
* `ma-1L-fe-gef-fp-p1.csv` 
* `ma-1L-fe-gef-fp-p01.csv` 

Second line fixed effects meta-analyses for platinum-based doublet chemotherapy (PBDC):

* `ma-2L-fe-pbdc-fp-p0.csv` 
* `ma-2L-fe-pbdc-fp-p00.csv` 
* `ma-2L-fe-pbdc-fp-p1.csv` 
* `ma-2L-fe-pbdc-fp-p01.csv` 

Second line fixed effects meta-analysis for osimertinib among patients with a T790M mutation:

* `ma-2L-fe-t790m-osi-fp-p0.csv`
* `ma-2L-fe-t790m-osi-fp-p00.csv` 
* `ma-2L-fe-t790m-osi-fp-p1.csv` 
* `ma-2L-fe-t790m-osi-fp-p01.csv` 

## Network meta-analysis files
First line fixed effects network meta-analyses of treatments relative to gefitnib:

* `nma-1L-fe-fp-p0.csv`
* `nma-1L-fe-fp-p00.csv`
* `nma-1L-fe-fp-p1.csv`
* `nma-1L-fe-fp-p01.csv`

# Lookup files
There are two types of lookup files.

The first type links treatments in the model to the numeric identifiers in the NMA:

* `tx-lookup-1L.csv`
* `tx-lookup-2L.csv`

The second type summarizes the parameters of the NMA in terms of the parametric and flexible
parametric multi-state survival models from the `hesim` package (which is used to simulate
the economic model):

* `params-lookup-1L.xlsx`
* `params-lookup-2L-pbdc.xlsx`
* `params-lookup-2L-t790m-osi.xlsx`
