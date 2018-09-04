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

#' Multi-state model paramters
#'
#' The parameters of the individual-level continuous time state transition model
#' estimated using the multi-state network meta-analysis. 
#'
#' @format A list containing parameter estimates from Weibull, Gompertz,
#' and 2nd order fractional polynomial survival distributions. Within each distribution,
#' parameters are contained in a list of lists with the outer list denoting a 
#' treatment line (1L, 2L, 2L+) and the inner list 
#' a \code{\link[hesim]{params_surv_list}} object from the \code{hesim} package.
"params_mstate"