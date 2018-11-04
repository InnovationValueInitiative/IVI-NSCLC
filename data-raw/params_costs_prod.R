rm(list = ls())
library("data.table")
library("readxl")

# Load in excel tables
wages <- data.table(read_excel("productivity-costs.xlsx", sheet = "wages"))
temporary_disability <- unlist(read_excel("productivity-costs.xlsx", sheet = "temporary_disability"))
permanent_disability <- unlist(data.table(read_excel("productivity-costs.xlsx", sheet = "permanent_disability")))

# Modify
wages[, n_sum := sum(n), by = "gender"]
wages[, prop := n/n_sum, by = "gender"]
setnames(wages, "source_n", "source_prop")
wages <- wages[, c("gender", "employment_status", "prop", "weekly_wage",
                   "source_prop", "source_wage"), with = FALSE]


# Save
params_costs_prod <- list(wages = wages,
                          temporary_disability = temporary_disability,
                          permanent_disability = permanent_disability)
save(params_costs_prod, file = "../data/params_costs_prod.rda", compress = "bzip2")
