# This script runs all files in the data-raw directory and creates all required
# datasets in the data directory. It should be run from inside the data-raw
# directory.

rm(list = ls())
unlink("../data/*") # deletes all files in data directory
source("params_mstate_nma.R")