rm(list = ls())
library("data.table")
library("readxl")
source("func.R")

# DRG payment interactive dataset: 
# https://data.cms.gov/Medicare-Inpatient/Inpatient-Prospective-Payment-System-IPPS-Provider/fm2n-hjj6
params_costs_ae <- data.table(read_excel("costs_ae.xlsx", sheet = "costs_ae"))
adj <- cpi_adj(params_costs_ae$year)
params_costs_ae[, mean := mean * adj]
params_costs_ae[, lower := lower * adj]
params_costs_ae[, upper := upper * adj]
params_costs_ae[, se := se * adj]
save(params_costs_ae, file = "../data/params_costs_ae.rda", compress = "bzip2")
