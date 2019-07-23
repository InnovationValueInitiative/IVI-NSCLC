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
The components of the package `iviNSCLC` are distributed in eight directories:
1. `R`: contains all R related code for simulating disease progression, costs, and QALYs. The code relies heavily on the `hesim` package.
2. `data-raw`: Contains `R` scripts to generate tables used in the package which are saved into `data`. The input used in the NMA are included in the subdrectory `mstate-nma-data`. The output of the NMA is stored in the subdirectory `mstate-nma`.
3. `data`: Contains tables created by the  `R` scripts in `data-raw`. These are the parameter estimates that load with the package and can be viewed with `data(package = "iviNSCLC")`. 
4. `docs`: Articles, references, and author information are stored here.
5. `man`: Manual and general documentation.
6. `pkgdown`: The package [website](https://innovationvalueinitiative.github.io/IVI-NSCLC/).
7. `tests`: Unit tests for the package. 
8. `vignettes`: Files for the package tutorial.
