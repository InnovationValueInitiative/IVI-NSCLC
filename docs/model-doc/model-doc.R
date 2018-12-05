rm(list = ls())
library("data.table")
library("ggplot2")
library("xtable")
library("pracma") # For numerical integration with trapezoid rule
library("iviNSCLC")
theme_set(theme_minimal())
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
     facet_wrap(~model) + 
     xlab("Month") + ylab("Proportion surviving") +
     scale_linetype_discrete(name = "") +
     theme(legend.position = "bottom")
ggsave("figs/surv-2L-pbdc.pdf", p, width = 7, height = 5)

## Second line (Osimertinib)
p <- ggplot(mstate_nma[line == 2 & mutation == 1], 
            aes(x = month, y = mean, linetype = outcome)) +
     geom_line() +
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
dic[, keep := ifelse(analysis == "1L MA for gefitinib, FE" |
                     analysis == "1L NMA, FE" |
                     analysis == "2L MA for PBDC, FE" |
                     analysis == "2L MA for osimertinib (T790M+), FE",
                     1, 0)]
p <- ggplot(dic[keep == 1],
            aes(x = model, y = dic, label = dic)) +
     geom_point() +
     geom_text(size = 3, hjust = -.2, vjust = 0) +
     facet_wrap(~analysis, scales = "free_y") + 
     theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
     scale_y_continuous(breaks = pretty_breaks(n = 8)) + 
     xlab("") +
     ylab("Deviance information criterion") 
ggsave("figs/dic.pdf", p, width = 7.5, height = 7.5)

# Resource use and costs -------------------------------------------------------
# Drug dosage
dosage <- params_costs_tx$dosage[, .(agent_name, dosage)]
print(xtable(dosage), 
      include.rownames = FALSE, include.colnames = FALSE,
      only.contents = TRUE, sanitize.text.function = identity,
      file = "tables/dosage.txt")

# Drug acquisition costs
params_costs_tx$acquisition_costs