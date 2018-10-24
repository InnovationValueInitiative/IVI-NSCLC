rm(list = ls())
library("data.table")
library("readxl")

params_costs_ae <- data.table(read_excel("costs_ae.xlsx", sheet = "costs_ae"))
save(params_costs_ae, file = "../data/params_costs_ae.rda", compress = "bzip2")
