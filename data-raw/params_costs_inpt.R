rm(list = ls())
library("data.table")
library("iviNSCLC")

# Outpatient costs by health state
state_name <- iviNSCLC:::pkg_env$state_names_start1L_4[1:3]
mean <- c(30000, 20000, 22000)
se <- c(3000, 4000, 4200)
ref <- c(NA, NA, NA)

# Save
params_costs_inpt <- data.table(state_name = state_name,
                              mean = mean,
                              se = se,
                              ref = ref)
save(params_costs_inpt, file = "../data/params_costs_inpt.rda", compress = "bzip2")
