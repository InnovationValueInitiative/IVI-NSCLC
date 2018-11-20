rm(list = ls())
library("data.table")
library("readxl")

params_costs_op <- data.table(read_excel("costs_op.xlsx", sheet = "costs_op"))
save(params_costs_op, file = "../data/params_costs_op.rda", compress = "bzip2")

