rm(list = ls())
library("data.table")
library("readxl")

params_costs_inpt <- data.table(read_excel("costs_inpt.xlsx", sheet = "costs_inpt"))
save(params_costs_inpt, file = "../data/params_costs_inpt.rda", compress = "bzip2")
