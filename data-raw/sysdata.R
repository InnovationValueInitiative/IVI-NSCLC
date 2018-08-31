rm(list = ls())
library("data.table")
treatments <- fread("treatments.csv")
devtools::use_data(treatments, internal = TRUE)