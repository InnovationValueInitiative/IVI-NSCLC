rm(list = ls())
library("data.table")
library("iviNSCLC")

# Utility by health state
## Nafees 2008 et al.
nafees_util_prog_2L_mean <- .6532 - 0.1798
nafees_util_prog_2L_se <- sqrt(0.02223^2 + 0.02169^2)

## Create data table
utility_states_name <- iviNSCLC:::pkg_env$state_names_start1L_4[1:3]
utility_states_mean <- c(.78, .6532, nafees_util_prog_2L_mean)
utility_states_se <- c(.01, .02223, nafees_util_prog_2L_se)
utility_states_ref <- c(NA, "nafees2008health", "nafees2008health")
utility_states <- data.table(state_name = utility_states_name,
                             mean = utility_states_mean,
                             se = utility_states_se,
                             ref = utility_states_ref)

# Utility loss from adverse events
disutility_ae_name <- c("Diarrhea", "Dry skin",
                        "Elevated alanine transaminase",
                        "Elevated aspartate transaminase",
                        "Eye problems", "Paronychia",
                      "Pneumonitis", "Pruritus", "Rash", "Stomatitis")
disutility_ae_abb <- c("diarrhea", "dry_skin", "alt", "ast", "eye_problems",
                       "paronychia", "pneumonitis",
                       "pruritus", "rash", "stomatitis")
disutility_ae_mean <- c(.0468, 0, 0, 0, 0, 0, 0, 0, -.03248, 0)
disutility_ae_se <- c(.01553, 0, 0, 0, 0, 0, 0, 0, .01171, 0)
disutility_ae_ref <- c("nafees2008health", NA, NA, NA, NA, NA, NA, NA, "nafees2008health", NA)
disutility_ae <- data.table(ae_name = disutility_ae_name,
                            ae_abb = disutility_ae_abb,
                            mean = disutility_ae_mean,
                            se = disutility_ae_se,
                            ref = disutility_ae_ref)


# Save
params_utility <- list(state_utility = utility_states,
                       ae_disutility = disutility_ae)
setorderv(params_utility$ae_disutility, "ae_abb")
save(params_utility, file = "../data/params_utility.rda", compress = "bzip2")
