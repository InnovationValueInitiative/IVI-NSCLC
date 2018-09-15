#' Treatments
#'
#' Dataset of treatments for EGFR positive NSCLC.
#'
#' @format A data.table object with 1 row for each treatment and w columns:
#' \describe{
#'   \item{name}{Name of treatment.}
#'   \item{abb}{Abbreviation for the treatment.}
#' }
"treatments"

#' Multi-state NMA parameters
#'
#' The regression coefficients of the continuous time state transition models
#' estimated using the multi-state network meta-analysis. 
#'
#' @format A list of \code{\link[hesim]{params_surv}} objects from the \code{hesim}
#' package. The list contains regression coefficient estimates from the Weibull,
#' Gompertz, and 2nd order fractional polynomial survival distributions. The 
#' parameter (e.g., scale and shape for the Weibull distribution) of each survival 
#' distribution are predicted as a function of treatments, health state transitions, and
#' treatment history (for 2L treatments).
"params_mstate_nma"