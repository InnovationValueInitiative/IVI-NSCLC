[![Travis-CI Build Status](https://travis-ci.org/InnovationValueInitiative/IVI-NSCLC.svg?branch=master)](https://travis-ci.org/InnovationValueInitiative/IVI-NSCLC)
[![Coverage Status](https://codecov.io/gh/InnovationValueInitiative/IVI-NSCLC/branch/master/graph/badge.svg)](https://codecov.io/gh/InnovationValueInitiative/IVI-NSCLC)

# IVI-NSCLC
`iviNSCLC` is an R package that runs the [Innovation and Value Initiative's (IVI's)](http://www.thevalueinitiative.org/)  non-small cell lung cancer (NSCLC) simulation model (the IVI-NSCLC model). The model simulates the costs, health outcomes, and risks associated with sequences of treatment including EGFR Tyrosine Kinase Inhibitors (TKIs), platinum-based doublet chemotherapy (PBDC), anti–vascular endothelial growth factor (anti-VEGF) therapy, and immune checkpoint inhibitors for patients with epidermal growth factor receptor (EGFR) positive NSCLC. 

# Installation
```r
# Install the development version from GitHub:
# install.packages("devtools")
devtools::install_github("InnovationValueInitiative/IVI-NSCLC")
```
# File structure
The componenst of the package `iviNSCLC` are distributed in eight directories:
1. `R`: contains all R related code that calls `hesim` package for perform simulations and estimates of health and economic values.
2. `data-raw`: Contains `R` code to generate tables used in the package which are saved into `data`.  Other input tables for adverse events, costs, productivity, and treatment information are included here too. 
3. `data`: Here are the tables created by the several `R` codes in `data-raw`. Posterior distributions derived from nma for various models, as well as costs, treatment, adverse evets parameter tables are also included here. 
4. `docs`: Articles, references and author information is stored here.
5. `man`: Manual and general documentation.
6. `pkgdown`: Website updates.
7. `tests`: R code to test of several components of the package.
8. `vignettes`: Files for tutorial of the package.

