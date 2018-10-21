rm(list = ls())
library("data.table")
library("readxl")

# Load in excel tables
dosage <- data.table(read_excel("treatment-costs.xlsx", sheet = "dosage"))
acquisition_costs <- data.table(read_excel("treatment-costs.xlsx", sheet = "acquisition_costs"))
discounts <- data.table(read_excel("treatment-costs.xlsx", sheet = "discounts"))
administration_costs <- data.table(read_excel("treatment-costs.xlsx", sheet = "administration_costs"))
lookup <- data.table(read_excel("treatment-costs.xlsx", sheet = "lookup"))

# Modify
dosage[, duration_days := ifelse(duration_days == "Until progression", Inf, duration_days)]
dosage[, duration_days := as.numeric(duration_days)]

# Save
params_costs_tx <- list(dosage = dosage,
                        acquisition_costs = acquisition_costs,
                        discounts = discounts,
                        administration_costs = administration_costs,
                        lookup = lookup)
save(params_costs_tx, file = "../data/params_costs_tx.rda", compress = "bzip2")
