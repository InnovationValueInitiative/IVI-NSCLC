rm(list = ls())
library("data.table")
library("readxl")

# Load
mcda_3 <- data.table(read_excel("mcda-defaults.xlsx", sheet = "3-state"))
mcda_4 <- data.table(read_excel("mcda-defaults.xlsx", sheet = "4-state"))

# Save
mcda_defaults <- list(three_state = mcda_3,
                      four_state = mcda_4)
save(mcda_defaults, file = "../data/mcda_defaults.rda", compress = "bzip2")
