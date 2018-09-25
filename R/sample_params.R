#' Sample parameters for PSA
#' 
#' Sample parameters from their joint probability distribution for a 
#' probabilistic sensitivity analysis (PSA).
#' @param n The number of random observations of the parameters to draw.
#' @param params_mstate_nma A list of \code{\link{params_surv}} objects, 
#' where each element in the list denotes a survival distribution. Should have
#' the same variable names as \code{\link{params_mstate_nma}}.
#' @return An object of class "sampled_params", which is a list with the
#' following components. 
#' @param params_utility Parameter estimates for health state utilities and
#' adverse event disutilities in the same format as \code{\link{params_utility}}.
#' \describe{
#' \item{mstate_nma}{A \code{\link{params_surv}} object where the number of 
#' rows in each \code{coef} element is equal to \code{n}.}
#' \item{utility}{A 3-dimensional array of matrices where each matrix represents
#' the utility values for a different treatment. The rows of each matrix are
#' random draws from a beta distribution and columns are health states.}
#' }
#' @examples
#' params <- sample_params(2)
#' names(params)
#' 
#' # Multi-state NMA
#' params$mstate_nma$weibull$coefs$shape
#' params$mstate_nma$weibull$dist
#' @export
sample_params <- function(n, params_mstate_nma = iviNSCLC::params_mstate_nma,
                          params_utility = iviNSCLC::params_utility){
  params <- list()
  params$mstate_nma <- sample_params_mstate_nma(n, params_mstate_nma)
  params$utility <- sample_params_utility(n, params_utility)
  
  class(params) <- "sampled_params"
  return(params)
}

sample_params_mstate_nma <- function(n, object){
  n_samples <- sapply(object, function(x) x$n_samples)
  if(!all(n_samples == n_samples[1])){
    msg <- paste0("The number of random draws from the posterior distribution in ",
                  "'params_mstate_nma' must be the same across probability distributions.")
    stop(msg, call. = FALSE)
  }
  if (n <= n_samples){
    sampled_rows <- sample.int(n_samples, n, replace = FALSE) 
  } else{
    warning("'n' is larger than the values of 'n_samples' in 'params_mstate_nma'.")
    sampled_rows <- sample.int(n_samples, n, replace = TRUE) 
  }  
  for (i in 1:length(object)){
    n_params <- length(object[[i]]$coefs)
    object[[i]]$n_samples <- n
    for (j in 1:n_params){
      object[[i]]$coefs[[j]] <- object[[i]]$coefs[[j]][sampled_rows, , drop = FALSE]
    }
  }  
  return(object)
}

sample_params_utility <- function(n, object){
  state_utility_params <- hesim::mom_beta(mean = object$state_utility$mean,
                                          sd = object$state_utility$se)
  state_utility_dist <- matrix(stats::rbeta(n * nrow(object$state_utility), 
                                            shape1 = state_utility_params$shape1,
                                            shape2 = state_utility_params$shape2),
                               nrow = n, byrow = TRUE)
  colnames(state_utility_dist) <- object$state_utility$state_name
  return(state_utility_dist)  
}
