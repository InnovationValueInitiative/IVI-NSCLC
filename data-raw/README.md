This directory contains the R scripts and data used to generate inputs for the model.
The inputs are stored as datasets in the `data` directory and load with the package. 
The file `Make-data.R` runs all R scripts needed to create these datasets. 

Output from network meta-analyses (NMAs)---which are processed into forms suitable for simulating 
the economic model---are stored in the `mstate-nma` (multi-state NMA) and
`ae-nma` (Adverse event NMA) directories. Data used to run the NMAs are in the 
`mstate-nma-data` and `ae-nma-data` directories. 

The `figs` directory contains figures used to assess some of the analyses 
contained in the various R scripts in `Make-data.R`.
