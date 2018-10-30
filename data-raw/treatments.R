rm(list = ls())
library("data.table")
treatments <- fread("treatments.csv")
treatments[, source_approval_yr := NULL]
current_yr <- year(Sys.Date())
treatments[, yrs_since_approval := current_yr - approval_yr]
save(treatments, file = "../data/treatments.rda", compress = "bzip2")