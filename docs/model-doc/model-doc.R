rm(list = ls())
library("data.table")
library("ggplot2")
library("pracma") # For numerical integration with trapezoid rule
library("iviNSCLC")
theme_set(theme_minimal())
# Run this script out of the docs/model-doc directory. 
# setwd("docs/model-doc")

# Functions --------------------------------------------------------------------
surv_mean <- function(x){
  # Compute mean survival times from survival curves by model and treatment.
  # Args:
  # x: A data table in the same format as iviNSCLC::mstate_nma_pfs and
  #    iviNSCLC::mstate_nma_os. Should be subset by line and mutation status.
  by_vars <- c("model", "tx_name")
  x_list <- split(x, by = by_vars)
  by <- unique(x[, by_vars, with = FALSE])
  means <- sapply(x_list, function (y) pracma::trapz(y$month, y$mean))
  res <- data.table(model = by$model,
                    tx_name = by$tx_name,
                    mean = means) 
  return(res)
}

# Transition probabilities (i.e., multi-state NMA) -----------------------------
# Progression free survival (1L)
dat <- mstate_nma_pfs[line == 1]

## Curves
p <- ggplot(dat, aes(x = month, y = mean, col = tx_name)) + geom_line() +
            facet_wrap(~model) +
            xlab("Month") + ylab("Progression-free survival") +
             scale_color_discrete(name = "") + theme(legend.position = "bottom")
ggsave("figs/pfs_1L.pdf", p, height = 8, width = 7)

## Median survival
surv_med_est <- hesim:::surv_quantile(dat, surv_cols = "mean", t = "month", by = c("model", "tx_name"))

## Mean survival
surv_mean_est <- surv_mean(dat)

# Overall survival (1L)
dat <- mstate_nma_os[line == 1]

## Curves
p <- ggplot(dat, aes(x = month, y = mean, col = tx_name)) + geom_line() +
            facet_wrap(~model) +
            xlab("Month") + ylab("Overall survival") +
             scale_color_discrete(name = "") + theme(legend.position = "bottom")
ggsave("figs/os_1L.pdf", p, height = 8, width = 7)

## Median survival
surv_med_est <- hesim:::surv_quantile(dat, surv_cols = "mean", t = "month", by = c("model", "tx_name"))

## Mean survival
surv_mean_est <- surv_mean(dat)
