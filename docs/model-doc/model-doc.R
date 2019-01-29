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

# Population -------------------------------------------------------------------
txt$AgeMean <- round(attr(age_dist, "mean"), 2)
txt$AgeSd <- round(attr(age_dist, "sd"), 2)

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

# Hazard ratios
## 1L
p <- ggplot(mstate_nma_hr[transition == "Stable to progression" & 
                            month > 0], 
            aes(x = month, y = median, col = tx_name)) +
     geom_line() +
     facet_wrap(~model, ncol = 2) +
     xlab("Month") + ylab("Hazard ratio") +
     scale_colour_discrete(name = "") +
     theme(legend.position = "bottom") +
     coord_cartesian(ylim = c(0, 1.3))
ggsave("figs/hr-1L.pdf", p, width = 8, height = 8)

## 1L with CrI
hr_1L_cri_plot <- function(model_name){
 p <- ggplot(mstate_nma_hr[transition == "Stable to progression" &
                            model == model_name & month > 0], 
            aes(x = month, y = median)) +
     geom_line() +
     facet_wrap(~tx_name, ncol = 2) +
      geom_ribbon(aes(ymin = l95, ymax = u95),
                alpha = 0.2) +
     xlab("Month") + ylab("Hazard ratio") +
     scale_colour_discrete(name = "") +
     theme(legend.position = "bottom") +
     coord_cartesian(ylim = c(0, 1.3))
 return(p)
}
p <- hr_1L_cri_plot("Weibull")
ggsave("figs/hr-1L-weibull.pdf", p, width = 8, height = 8)
p <- hr_1L_cri_plot("Fractional polynomial (0, 0)")
ggsave("figs/hr-1L-fp-00.pdf", p, width = 8, height = 8)
p <- hr_1L_cri_plot("Fractional polynomial (0, 1)")
ggsave("figs/hr-1L-fp-01.pdf", p, width = 8, height = 8)
p <- hr_1L_cri_plot("Gompertz")
ggsave("figs/hr-1L-gompertz.pdf", p, width = 8, height = 8)

# Hazards
## Main body without CrIs
haz_plot <- function(line, tx_name, ylim){
  line_env <- line
  tx_name_env <- tx_name
  p <- ggplot(mstate_nma_hazard[line == line_env & 
                                tx_name == tx_name_env &
                                month > 0],
            aes(x = month, y = median, col = model)) +
     geom_line() +
     facet_wrap(~transition, ncol = 2) +
     xlab("Month") + ylab("Hazard") +
     scale_colour_discrete(name = "") +
     theme(legend.position = "bottom") + 
     coord_cartesian(ylim = ylim)
  return(p)
}
p <- haz_plot(line = 1, tx_name = "gefitinib", ylim = c(0, .3))
ggsave("figs/hazard-1L-gef.pdf", p, width = 8, height = 8)
p <- haz_plot(line = 2, tx_name = "osimertinib", ylim = c(0, .5))
ggsave("figs/hazard-2L-t790m-osi.pdf", p, width = 7, height = 5)
p <- haz_plot(line = 2, tx_name = "PBDC", ylim = c(0, .5))
ggsave("figs/hazard-2L-pbdc.pdf", p, width = 7, height = 5)

## With CrIs
haz_cri_plot <- function(line, tx_name, ylim, model){
  line_env <- line
  tx_name_env <- tx_name
  model_env <- model
  p <- ggplot(mstate_nma_hazard[line == line_env & 
                                tx_name == tx_name_env &
                                model == model_env &
                                month > 0],
            aes(x = month, y = median)) +
     geom_line() +
     geom_ribbon(aes(ymin = l95, ymax = u95),
                alpha = 0.2) +
     facet_wrap(~transition, ncol = 2) +
     xlab("Month") + ylab("Hazard") +
     scale_colour_discrete(name = "") +
     theme(legend.position = "bottom") + 
     coord_cartesian(ylim = ylim)
  return(p)
}

### 1L
p <- haz_cri_plot(line = 1, tx_name = "gefitinib", ylim = c(0, .4),
                  model = "Weibull")
ggsave("figs/hazard-1L-gef-weibull.pdf", p, width = 8, height = 8)
p <- haz_cri_plot(line = 1, tx_name = "gefitinib", ylim = c(0, .4),
                  model = "Fractional polynomial (0, 0)")
ggsave("figs/hazard-1L-gef-fp-00.pdf", p, width = 8, height = 8)
p <- haz_cri_plot(line = 1, tx_name = "gefitinib", ylim = c(0, .4),
                  model = "Fractional polynomial (0, 1)")
ggsave("figs/hazard-1L-gef-fp-01.pdf", p, width = 8, height = 8)
p <- haz_cri_plot(line = 1, tx_name = "gefitinib", ylim = c(0, .4),
                  model = "Gompertz")
ggsave("figs/hazard-1L-gef-gompertz.pdf", p, width = 8, height = 8)

### 2L
p <- haz_cri_plot(line = 2, tx_name = "osimertinib", ylim = c(0, .4),
                  model = "Weibull")
ggsave("figs/hazard-2L-t790m-osi-weibull.pdf", p, width = 7, height = 5)
p <- haz_cri_plot(line = 2, tx_name = "PBDC", ylim = c(0, .4),
                  model = "Weibull")
ggsave("figs/hazard-2L-pbdc-weibull.pdf", p, width = 7, height = 5)

p <- haz_cri_plot(line = 2, tx_name = "osimertinib", ylim = c(0, .4),
                  model = "Fractional polynomial (0, 0)")
ggsave("figs/hazard-2L-t790m-osi-fp-00.pdf", p, width = 7, height = 5)
p <- haz_cri_plot(line = 2, tx_name = "PBDC", ylim = c(0, .4),
                  model = "Fractional polynomial (0, 0)")
ggsave("figs/hazard-2L-pbdc-fp-00.pdf", p, width = 7, height = 5)

p <- haz_cri_plot(line = 2, tx_name = "osimertinib", ylim = c(0, .4),
                  model = "Fractional polynomial (0, 1)")
ggsave("figs/hazard-2L-t790m-osi-fp-01.pdf", p, width = 7, height = 5)
p <- haz_cri_plot(line = 2, tx_name = "PBDC", ylim = c(0, .4),
                  model = "Fractional polynomial (0, 1)")
ggsave("figs/hazard-2L-pbdc-fp-01.pdf", p, width = 7, height = 5)

p <- haz_cri_plot(line = 2, tx_name = "osimertinib", ylim = c(0, .4),
                  model = "Gompertz")
ggsave("figs/hazard-2L-t790m-osi-gompertz.pdf", p, width = 7, height = 5)
p <- haz_cri_plot(line = 2, tx_name = "PBDC", ylim = c(0, .4),
                  model = "Gompertz")
ggsave("figs/hazard-2L-pbdc-gompertz.pdf", p, width = 7, height = 5)

# PFS/OS curves
mstate_nma <- rbind(data.table(mstate_nma_pfs, outcome = "PFS"),
                    data.table(mstate_nma_os, outcome = "OS")) 

## 1L
p <- ggplot(mstate_nma[line == 1], 
            aes(x = month, y = median, col = tx_name, linetype = outcome)) +
     geom_line() +
     facet_wrap(~model) + 
     xlab("Month") + ylab("Proportion surviving") +
     scale_color_discrete(name = "") +
     scale_linetype_discrete(name = "") +
     theme(legend.position = "bottom")
ggsave("figs/surv-1L.pdf", p, width = 8, height = 8)

## 2L (PBDC)
p <- ggplot(mstate_nma[line == 2 & mutation == 0], 
            aes(x = month, y = median, linetype = outcome)) +
     geom_line() +
     geom_ribbon(aes(ymin = l95, ymax = u95),
                alpha = 0.2) + 
     facet_wrap(~model) + 
     xlab("Month") + ylab("Proportion surviving") +
     scale_linetype_discrete(name = "") +
     theme(legend.position = "bottom")
ggsave("figs/surv-2L-pbdc.pdf", p, width = 7, height = 5)

## 2L (osimertinib)
p <- ggplot(mstate_nma[line == 2 & mutation == 1], 
            aes(x = month, y = median, linetype = outcome)) +
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

# DIC
dic <- fread("output/dic.csv")
dic[, dic := formatC(dic, format = "d", big.mark = ",")]

## 1L NMA
tbl <- dic[line == 1 & effects == "Relative"]
tbl <- dcast(tbl, model ~ method + sd_tx_effect + pd_tx_effect,
                value.var = "dic")
setcolorder(tbl, c("model", 
                 "FE_none_none", "FE_none_constant",
                  "RE_none_none", "RE_none_constant"))
print(xtable(tbl), 
      include.rownames = FALSE, include.colnames = FALSE,
      only.contents = TRUE, sanitize.text.function = identity,
      file = "tables/dic-1L-nma.txt")

## 1L MA
tbl <- dic[line == 1 & effects == "Absolute", .(model, dic)]
print(xtable(tbl), 
      include.rownames = FALSE, include.colnames = FALSE,
      only.contents = TRUE, sanitize.text.function = identity,
      file = "tables/dic-1L-ma-gef.txt")

## 2L MA
tbl <- dic[line == 2 & effects == "Absolute",]
tbl <- dcast(tbl, model ~ method + mutation,
                value.var = "dic")
print(xtable(tbl), 
      include.rownames = FALSE, include.colnames = FALSE,
      only.contents = TRUE, sanitize.text.function = identity,
      file = "tables/dic-2L-ma.txt")

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
ae_disutility[, ref := ifelse(is.na(ref),
                              "N/A",
                              paste0("\\citet{", ref, "}"))]
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
#inpt_costs[state_name == "S1", ref := "None"]
inpt_costs[, mean := formatC(mean, format = "d", big.mark = ",")]
inpt_costs[, se := formatC(se, format = "d", big.mark = ",")]
print(xtable(inpt_costs), 
      include.rownames = FALSE, include.colnames = FALSE,
      only.contents = TRUE, sanitize.text.function = identity,
      file = "tables/inpt_costs.txt")

## Outpatient
op_costs = copy(params_costs_op)
op_costs[, ref := paste0("\\citet{", ref, "}")]
#op_costs[state_name == "S1", ref := "None"]
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
cols <- colnames(wages)[colnames(wages) != "gender"]
print(xtable(wages[gender == "Male", cols, with = FALSE]), 
      include.rownames = FALSE, include.colnames = FALSE,
      only.contents = TRUE, sanitize.text.function = identity,
      hline.after = NULL,
      file = "tables/wages_male.txt")
print(xtable(wages[gender == "Female", cols, with = FALSE]), 
      include.rownames = FALSE, include.colnames = FALSE,
      only.contents = TRUE, sanitize.text.function = identity,
      hline.after = NULL,
      file = "tables/wages_female.txt")

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