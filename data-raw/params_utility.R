rm(list = ls())
library("data.table")
library("readxl")

# Load
params_state_utility <- data.table(read_excel("utility.xlsx", sheet = "state_utility"))
params_ae_disutility <- data.table(read_excel("utility.xlsx", sheet = "ae_disutility"))

# Save
params_utility <- list(state_utility = params_state_utility,
                       ae_disutility = params_ae_disutility)
setorderv(params_utility$ae_disutility, "ae_abb")
save(params_utility, file = "../data/params_utility.rda", compress = "bzip2")
