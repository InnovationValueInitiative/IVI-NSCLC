rm(list = ls())
library("data.table")
library("ggplot2")
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

surv_results <- function(mstate_nma, 
                         outcome = c("PFS", "OS"),
                         line,
                         mutation = NA){
  outcome <- match.arg(outcome)
  line_env <- line
  dat <- mstate_nma[line == line_env]
  if (!is.na(mutation)){
    mutation_env <- mutation
    dat <- mstate_nma[mutation == mutation_env]
  }
  
  # Curves
  if (outcome == "PFS"){
    y_lab <- "Progression-free survival"
    filename <- paste0("figs/", "pfs-", line, "L.pdf")
  } else{
      y_lab <- "Overall survival"
      filename <- paste0("figs/", "os-", line, "L.pdf")
  }
  p <- ggplot(dat, aes(x = month, y = mean, col = tx_name)) + geom_line() +
            facet_wrap(~model) +
            xlab("Month") + ylab(y_lab) +
             scale_color_discrete(name = "") + theme(legend.position = "bottom")
  ggsave(filename, p, height = 8, width = 7)
  
  # Median survival
  surv_med_est <- hesim:::surv_quantile(dat, surv_cols = "mean", 
                                        t = "month", by = c("model", "tx_name"))
  
  # Mean survival
  surv_mean_est <- surv_mean(dat)
  
  # Return
  return(list(p_curves = p, 
              median = surv_med_est,
              mean = surv_mean_est))
  
}

pfs_1L <- surv_results(mstate_nma_pfs, outcome = "PFS", line = 1)
os_1L <- surv_results(mstate_nma_os, outcome = "OS", line = 1)