rm(list = ls())
library("data.table")
library("ggplot2")
library("xtable")
library("pracma") # For numerical integration with trapezoid rule
library("iviNSCLC")
theme_set(theme_minimal())
txt <- list() # List for storing numbers to use in text of model documentation
# Run this script out of the docs/model-doc directory. 
# setwd("docs/model-doc")

# Transition probabilities (i.e., multi-state NMA) -----------------------------
surv_mean <- function(x){
  # Compute mean survival times from survival curves by model and treatment.
  # Args:
  # x: A data table in the same format as iviNSCLC::mstate_nma_pfs and
  #    iviNSCLC::mstate_nma_os.
  by_vars <- c("model", "tx_name", "line", "mutation")
  x_list <- split(x, by = by_vars)
  by <- unique(x[, by_vars, with = FALSE])
  means <- sapply(x_list, function (y) pracma::trapz(y$month, y$mean))
  res <- data.table(model = by$model,
                    tx_name = by$tx_name,
                    mean = means) 
  return(res)
}

# PFS/OS curves
mstate_nma <- rbind(data.table(mstate_nma_pfs, outcome = "PFS"),
                    data.table(mstate_nma_os, outcome = "OS")) 

## First line
p <- ggplot(mstate_nma[line == 1], 
            aes(x = month, y = mean, col = tx_name, linetype = outcome)) +
     geom_line() +
     facet_wrap(~model) + 
     xlab("Month") + ylab("Proportion surviving") +
     scale_color_discrete(name = "") +
     scale_linetype_discrete(name = "") +
     theme(legend.position = "bottom")
ggsave("figs/surv-1L.pdf", p, width = 8, height = 8)

## Second line (PBDC)
p <- ggplot(mstate_nma[line == 2 & mutation == 0], 
            aes(x = month, y = mean, linetype = outcome)) +
     geom_line() +
     geom_ribbon(aes(ymin = l95, ymax = u95),
                alpha = 0.2) + 
     facet_wrap(~model) + 
     xlab("Month") + ylab("Proportion surviving") +
     scale_linetype_discrete(name = "") +
     theme(legend.position = "bottom")
ggsave("figs/surv-2L-pbdc.pdf", p, width = 7, height = 5)

## Second line (Osimertinib)
p <- ggplot(mstate_nma[line == 2 & mutation == 1], 
            aes(x = month, y = mean, linetype = outcome)) +
     geom_line() +
      geom_ribbon(aes(ymin = l95, ymax = u95),
                alpha = 0.2) +
     facet_wrap(~model) + 
     xlab("Month") + ylab("Proportion surviving") +
     scale_linetype_discrete(name = "") +
     theme(legend.position = "bottom")
ggsave("figs/surv-2L-t790m-osi.pdf", p, width = 7, height = 5)

# Median survival
surv_med_est <- hesim:::surv_quantile(mstate_nma, 
                                      surv_cols = c("mean", "l95", "u95"),
                                      t = "month", 
                                      by = c("outcome", "model", "tx_name", "line", "mutation"))
surv_med_est[, tx_name := ifelse(tx_name == "osimertinib" & line == 2,
                                 "osmiternib (2L)", tx_name)]
surv_med_est[, tx_name := ifelse(tx_name == "PBDC" & line == 2,
                                 "PBDC (2L)", tx_name)]
surv_med_est[, tx_name := factor(tx_name,
                                 levels = c("gefitinib", "erlotinib", "afatinib",
                                            "dacomitinib", "osimertinib",
                                            "osmiternib (2L)", "PBDC (2L)"))]
surv_med_est[, fill_var := ifelse(line == 1, "1L",
                                   ifelse(mutation == 0,
                                          "2L", "2L (T790M+)"))]
surv_med_est[, fill_var := factor(fill_var, levels = c("1L",  "2L (T790M+)", "2L"))]

## PFS
p <- ggplot(surv_med_est[outcome == "PFS"], 
       aes(x = tx_name, y = quantile_mean, fill = factor(fill_var),
           label = quantile_mean)) +
      geom_bar(stat = "identity", position = "dodge") + facet_wrap(~model) + 
      ylab("Median survival") + xlab("Month") +
      scale_fill_discrete(name = "") +
      geom_text(size = 3, position = position_stack(vjust = 0.4)) +
      geom_errorbar(aes(ymin = quantile_l95,
                        ymax = quantile_u95), 
                    width = .2,
                    col = "grey") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("figs/medsurv-pfs.pdf", p, width = 7, height = 5)

## OS
p <- ggplot(surv_med_est[outcome == "OS"], 
       aes(x = tx_name, y = quantile_mean, fill = factor(fill_var),
           label = quantile_mean)) +
      geom_bar(stat = "identity", position = "dodge") + facet_wrap(~model) + 
      ylab("Median survival") + xlab("Month") +
      scale_fill_discrete(name = "") +
      geom_text(size = 3, position = position_stack(vjust = 0.4)) +
      geom_errorbar(aes(ymin = quantile_l95,
                        ymax = quantile_u95), 
                    width = .2,
                    col = "grey") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("figs/medsurv-os.pdf", p, width = 7, height = 5)

# Mean survival
surv_mean(mstate_nma_pfs)
surv_mean(mstate_nma_os)

# Utility ----------------------------------------------------------------------
# Utility by health state
state_utility <- copy(params_utility$state_utility)
state_utility[, mean := formatC(mean, format = "f", digits = 4)]
state_utility[, se := formatC(se, format = "f", digits = 4)]
state_utility[, ref := paste0("\\citet{", ref, "}")]
print(xtable(state_utility), 
      include.rownames = FALSE, include.colnames = FALSE,
      only.contents = TRUE, sanitize.text.function = identity,
      file = "tables/state-utility.txt")

# Disutility by adverse event
ae_disutility <- params_utility$ae_disutility[ , .(ae_name, mean, se, ref)]
ae_disutility[, mean := formatC(mean, format = "f", digits = 4)]
ae_disutility[, se := formatC(se, format = "f", digits = 4)]
ae_disutility[, ref := paste0("\\citet{", ref, "}")]
print(xtable(ae_disutility), 
      include.rownames = FALSE, include.colnames = FALSE,
      only.contents = TRUE, sanitize.text.function = identity,
      file = "tables/ae-disutility.txt")

# Health care sector costs -----------------------------------------------------
# Treatment costs
## Drug dosage
dosage <- params_costs_tx$dosage[, .(agent_name, dosage)]
print(xtable(dosage), 
      include.rownames = FALSE, include.colnames = FALSE,
      only.contents = TRUE, sanitize.text.function = identity,
      file = "tables/dosage.txt")

## Drug acquisition costs
acq_costs <- params_costs_tx$acquisition_costs[, .(agent_name, strength, acquisition_cost)]
acq_costs <- acq_costs[agent_name %in% unique(dosage$agent_name)]
print(xtable(acq_costs), 
      include.rownames = FALSE, include.colnames = FALSE,
      only.contents = TRUE, sanitize.text.function = identity,
      file = "tables/acq_costs.txt")

## Drug administration
admin_costs <- params_costs_tx$administration_costs
print(xtable(admin_costs), 
      include.rownames = FALSE, include.colnames = FALSE,
      only.contents = TRUE, sanitize.text.function = identity,
      file = "tables/admin_costs.txt")

# Inpatient and outpatient costs
## Inpatient
inpt_costs = copy(params_costs_inpt)
inpt_costs[, ref := paste0("\\citet{", ref, "}")]
inpt_costs[state_name == "S1", ref := "None"]
inpt_costs[, mean := formatC(mean, format = "d", big.mark = ",")]
inpt_costs[, se := formatC(se, format = "d", big.mark = ",")]
print(xtable(inpt_costs), 
      include.rownames = FALSE, include.colnames = FALSE,
      only.contents = TRUE, sanitize.text.function = identity,
      file = "tables/inpt_costs.txt")

## Outpatient
op_costs = copy(params_costs_op)
op_costs[, ref := paste0("\\citet{", ref, "}")]
op_costs[state_name == "S1", ref := "None"]
op_costs[, mean := formatC(mean, format = "d", big.mark = ",")]
op_costs[, se := formatC(se, format = "d", big.mark = ",")]
print(xtable(op_costs), 
      include.rownames = FALSE, include.colnames = FALSE,
      only.contents = TRUE, sanitize.text.function = identity,
      file = "tables/op_costs.txt")

# Adverse event costs
ae_costs  <- params_costs_ae[, .(ae_name, mean, lower, upper, ref)]
ae_costs[, ref := ifelse(grepl("DRG", ref) == 0,
                         paste0("\\citet{", ref, "}"),
                                ref)]
ae_costs[, mean := formatC(mean, format = "d", big.mark = ",")]
ae_costs[, lower := formatC(lower, format = "d", big.mark = ",")]
ae_costs[, upper := formatC(upper, format = "d", big.mark = ",")]
print(xtable(ae_costs), 
      include.rownames = FALSE, include.colnames = FALSE,
      only.contents = TRUE, sanitize.text.function = identity,
      file = "tables/ae_costs.txt")

# Productivity -----------------------------------------------------------------
# Wages
wages <- params_costs_prod$wages[, .(gender, employment_status, prop, weekly_wage)]
wages[, gender := factor(gender, levels = c("female", "male"),
                         labels = c("Female", "Male"))]
wages[, employment_status := factor(employment_status,
                                    levels = c("full", "part", "unemployed"),
                                    labels = c("Full-time", "Part-time", "Unemployed"))]
setnames(wages, "prop", "percent")
wages[, percent := paste0(formatC(100 * percent, format = "f", digits = 1),
                          "\\%")]
wages[, weekly_wage := paste0("\\$",
                              formatC(weekly_wage, format = "d", big.mark = ","))]
print(xtable(wages), 
      include.rownames = FALSE, include.colnames = FALSE,
      only.contents = TRUE, sanitize.text.function = identity,
      file = "tables/wages.txt")

# Temporary disability
tmp_disability <- params_costs_prod$temporary_disability
txt$MissedDaysEst <- tmp_disability["missed_days_est"]
txt$MissedDaysLower <- tmp_disability["missed_days_lower"]
txt$MissedDaysUpper <- tmp_disability["missed_days_upper"]

# Permanent disability
perm_disability <- params_costs_prod$permanent_disability
txt$HoursReductionEst <- perm_disability["hours_reduction_est"]
txt$HoursReductionLower <- perm_disability["hours_reduction_lower"]
txt$HoursReductionUpper <- perm_disability["hours_reduction_upper"]

# Text for model documentation -------------------------------------------------
# convert statistics to data frame
txtstats <- data.frame(do.call(rbind, txt))
rownames(txtstats) <- gsub("txt.", "", rownames(txtstats))

# output to text file to input into latex
txtstats$def <-  "\\def"
names(txtstats)[1] <- "value"
txtstats$value <- as.character(txtstats$value)
txtstats <- data.frame(def = txtstats$def, name = rownames(txtstats), value =  txtstats$value)
txtstats$output <- paste(txtstats[, 1], " ", "\\", txtstats[, 2],
                         "{", txtstats[, 3], "}", sep = "")
fileConn <- file("output/txtstats.txt")
writeLines(txtstats$output, fileConn)
close(fileConn)