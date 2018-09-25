rm(list = ls())
library("data.table")
library("iviNSCLC")
rm(list = ls())
library("data.table")
library("iviNSCLC")

# Outpatient costs by health state
state_name <- iviNSCLC:::pkg_env$state_names_start1L_4[1:3]
mean <- c(1500, 2000, 2200)
se <- c(300, 400, 420)
ref <- c(NA, NA, NA)

# Save
params_costs_op <- data.table(state_name = state_name,
                              mean = mean,
                              se = se,
                              ref = ref)
save(params_costs_op, file = "../data/params_costs_op.rda", compress = "bzip2")
