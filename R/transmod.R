create_strategies <- function(struct, txseqs){
  start_line <- attributes(struct$txseqs)$start_line
  
  # Strategies table
  n_strategies <- length(struct$txseqs)
  strategy_id <- 1:n_strategies
  
  container <- matrix(NA, nrow = n_strategies, ncol = 5)
  colnames(container) <- c("abb_first", "abb_second_pos",
                           "abb_second_neg", "abb_second_plus_pos",
                           "abb_second_plus_neg")
  for (i in 1:n_strategies){
    txseqs_i <- unlist(struct$txseqs[[i]])
    if (txseqs_i[1] == "osimertinib" & struct$n_states == "four" & 
        start_line == "first"){
      msg <- paste0("There is no evidence to parameterize a model with four health states ", 
                    "when osimertinib is used as a first line treatment. ",
                    "Use a model with three health states for sequences beginning ",
                    "with osimertinib instead.")
      stop(msg, call. = FALSE)   
    }
    abb <- treatments$abb[match(txseqs_i, treatments$name)]
    names(abb) <- names(txseqs_i)
    container[i, "abb_first"] <- abb["first"]
    container[i, "abb_second_pos"] <- abb["second.pos"]
    container[i, "abb_second_neg"] <- abb["second.neg"]
    container[i, "abb_second_plus_pos"] <- abb["second_plus.pos"]
    container[i, "abb_second_plus_neg"] <- abb["second_plus.neg"]
  }  
  strategies <- data.table(strategy_id, container)
  return(strategies)
}

mutation_data <- function(txseqs, mutation_prob, patients){
  n_patients <- nrow(patients)
  n_strategies <- length(txseqs)
  if(length(mutation_prob) == 1){
    mutation_prob <- rep(mutation_prob, n_strategies)
  }
  n_mutations <- round(mutation_prob * n_patients)
  mutations <- sapply(n_mutations, function (x) c(rep(1, x), 
                                     rep(0, n_patients - x)))
  return(c(mutations))
}

dist_params <- function(dist = "weibull"){
  dist <- match.arg(dist)
  if (dist == "weibull"){
    return (c("scale", "shape"))
  }
}

transmod_vars <- function(struct, data){
  transition_id <- mutation <- tx_abb <- NULL # stop no visible binding CRAN warning
  
  start_line <- attributes(struct$txseqs)$start_line
  params <- dist_params(struct$dist)
  
  # Treatment variables by transition
  if (start_line == "first"){
      if (struct$n_states == "four"){
        data[, "tx_abb" := ifelse(transition_id %in% c(1, 2), abb_first, NA)]
        data[, "tx_abb" := ifelse(transition_id %in% c(3, 4, 5) & mutation == 1, abb_second_pos, tx_abb)]
        data[, "tx_abb" := ifelse(transition_id %in% c(3, 4, 5) & mutation == 0, abb_second_neg, tx_abb)]
        data[, "tx_hist" := ifelse(transition_id %in% c(3, 4, 5), abb_first, NA)]
      } else {
        data[, "tx_abb" := abb_first]
        data[, "tx_hist" := NA]
      }
  } else{
    stop("Functionality has not yet been added when starting at second line.")
  }

  
  
  # Treatments
  ## First line
  abb_first <- unique(data$abb_first)
  abb_first <- c("osi", abb_first[abb_first != "osi"])
  
  ## Second line
  ## Positive mutation
  abb_second_pos <- unique(data$abb_second_pos)
  abb_second_pos <- c("osi", abb_second_pos[abb_second_pos != "osi"])
  
  ## Negative mutation
  abb_second_neg <- unique(data$abb_second_neg)
  abb_second_neg <- c("pbdc", abb_second_neg[abb_second_neg != "pbdc"])  

  # Create dummy variables
  for (j in 1:length(params)){
    
    # (1) Model starting at first line
    if (start_line == "first"){
      
      # Loop through first line treatments
      for (i in 1:length(abb_first)){
        abb_first_i <- abb_first[i]
        if (abb_first_i =="osi"){
          data[, paste0("osi_s1p1_", params[j]) := ifelse(transition_id == 1, 1, 0)]
          data[, paste0("osi_s1d_", params[j]) := ifelse(transition_id == 2, 1, 0)]
          if (struct$n_states == "three"){
            data[, paste0("osi_p1d_", params[j]) := ifelse(transition_id == 3, 1, 0)]
          }
        } else {
          data[, paste0("d_", abb_first_i, "_s1p1_", params[j]) := ifelse(transition_id == 1 &
                                                                           tx_abb == abb_first_i, 
                                                                           1, 0)]
          data[, paste0("d_", abb_first_i, "_s1d_", params[j]) := ifelse(transition_id == 2 &
                                                                          tx_abb == abb_first_i, 
                                                                          1, 0)]    
          if (struct$n_states == "three"){
            data[, paste0("d_", abb_first_i, "_p1d_", params[j]) := ifelse(transition_id == 3 &
                                                                            tx_abb == abb_first_i, 
                                                                            1, 0)]              
          } 
        }
      } # end first line treatment loop      
      
      if (struct$n_states == "four") {
        # Loop through second line treatments for T790M positive patients
        data[, paste0("osi_s2p2_", params[j]) := ifelse(transition_id == 3 & mutation == 1, 1, 0)]
        data[, paste0("osi_s2d_", params[j]) := ifelse(transition_id == 4 & mutation == 1, 1, 0)]
        data[, paste0("osi_p2d_", params[j]) := ifelse(transition_id == 5 & mutation == 1, 1, 0)]
        
        # Loop through second line treatments for T790M negative patients
        for (i in 1:length(abb_second_neg)){    
          abb_second_neg_i <- abb_second_neg[i]
          if (abb_second_neg_i =="pbdc"){
            data[, paste0("pbdc_s2p2_", params[j]) := ifelse(transition_id == 3 & mutation == 0, 1, 0)]
            data[, paste0("pbdc_s2d_", params[j]) := ifelse(transition_id == 4 & mutation == 0, 1, 0)]
            data[, paste0("pbdc_p2d_", params[j]) := ifelse(transition_id == 5 & mutation == 0, 1, 0)]
          } else{
              data[, paste0("d_", abb_second_neg_i, "_s2p2_", params[j]) := ifelse(transition_id == 3 &
                                                                                    mutation == 0 &
                                                                                    tx_abb == abb_second_neg_i, 
                                                                                    1, 0)]
              data[, paste0("d_", abb_second_neg_i, "_s2d_", params[j]) := ifelse(transition_id == 4 &
                                                                                    mutation == 0 &
                                                                                    tx_abb == abb_second_neg_i, 
                                                                                    1, 0)]    
              data[, paste0("d_", abb_second_neg_i, "_p2d_", params[j]) := ifelse(transition_id == 5 &
                                                                                    mutation == 0 &
                                                                                    tx_abb == abb_second_neg_i, 
                                                                                    1, 0)]               
          }
        } # end T790M negative second line loop        
      } # end conditional requiring four-state model
    } # end model starting at first line (i.e., start_line = "first")
    
    
    # (2) Model starting at second line
    # if (start_line == "second"){
    #   # Second line T790M positive treatments 
    # }
    
  } # end parameter loop
  return(data)
}

#' Data for transition model
#' 
#' Create data used to simulate health state transitions with a 
#' continuous time state transition model (CTSTM) given the parameters from the
#' multi-state NMA (i.e., \code{\link{params_mstate_nma}}). The included variables are 
#' a function of the selected treatment sequences and the modeled patients.
#' @param struct A \code{\link{model_structure}} object.
#' @param trans_mat A transition matrix as returned by \code{\link{create_trans_mat}}.
#' @param patients A data table returned from \code{\link{create_patients}}.
#' @param mutation_prob The probability of a T790M mutation. A vector of either
#' of length 1 (the probability is constant across first line treatments in 
#' \code{txseqs}) or of length equal to the number of treatment sequences in
#' \code{txseqs}. 
#' @return An object of class "expanded_hesim_data" from the 
#' \href{https://innovationvalueinitiative.github.io/hesim/}{hesim} package, which
#' is a data table with one observation for each treatment strategy 
#' (i.e., treatment sequence), patient, and transition combination. The survival
#' distribution is stored as a "dist" attribute.  
#' @seealso \code{\link{create_transmod_params}}
#' @examples
#' ## Treatment sequences
#' txseq1 <- txseq(first = "erlotinib",
#'                 second = c("osimertinib", "PBDC"),
#'                 second_plus = c("PBDC + bevacizumab", "PBDC + bevacizumab"))
#' txseq2 <- txseq(first = "gefitinib",
#'                 second = c("osimertinib", "PBDC"),
#'                 second_plus = c("PBDC + bevacizumab", "PBDC + bevacizumab"))
#' txseqs <- txseq_list(seq1 = txseq1, seq2 = txseq2)
#'
#' # Patient population
#' pats <- create_patients(n = 2)
#'
#' ## Model structure
#' struct <- model_structure(txseqs, dist = "weibull")
#' tmat <- create_trans_mat(struct)
#'
#' ## Data
#' transmod_data <- create_transmod_data(struct, tmat, pats, mutation_prob = .45)
#' head(transmod_data)
#' @export
create_transmod_data <- function(struct, trans_mat, patients, mutation_prob = .45){
  strategies <- create_strategies(struct)
  hesim_data <- hesim::hesim_data(strategies = strategies,
                                  patients = patients)
  data <- hesim::expand(hesim_data, by = c("strategies", "patients"))
  
  # Add mutations
  mutation <- mutation_data(struct$txseqs, mutation_prob, patients)
  data[, "mutation" := mutation]
  
  # Expand for each transition
  n_transitions <- max(trans_mat, na.rm = TRUE)
  data <- data[rep(seq_len(nrow(data)), each = n_transitions)]
  data[, "transition_id" := 1:n_transitions]
  setcolorder(data, c("strategy_id", "patient_id", "transition_id"))
  setattr(data, "id_vars", c("strategy_id", "patient_id", "transition_id"))

  # Add treatment effect variables
  data <- transmod_vars(struct, data)
  data[, c("abb_first", "abb_second_pos", "abb_second_neg", 
           "abb_second_plus_pos", "abb_second_plus_neg") := NULL]
  setattr(data, "dist", struct$dist)
  return(data[, ])
}

sample_params_mstate_nma <- function(n, object){
  n_samples <- object$n_samples
  if (n <= n_samples){
    sampled_rows <- sample.int(n_samples, n, replace = FALSE) 
  } else{
    warning("'n' is larger than the values of 'n_samples' in 'params_mstate_nma'.")
    sampled_rows <- sample.int(n_samples, n, replace = TRUE) 
  }  
  n_params <- length(object$coefs)
  object$n_samples <- n
  for (j in 1:n_params){
    object$coefs[[j]] <- object$coefs[[j]][sampled_rows, , drop = FALSE]
  }
  return(object)
}

#' Parameters for transition model
#' 
#' Extract parameters from a multi-state NMA for use with the data table returned by
#' \code{\link{create_transmod_data}}, which are used to simulate health state 
#' transitions with a continuous time state transition model (CTSTM).
#' @param n The number of random observations of the parameters to draw.
#' @param data A data table of class "expanded_hesim_data" returned from 
#' \code{\link{create_transmod_data}}.
#' @param params_mstate_nma A list of \code{\link{params_surv}} objects, 
#' where each element in the list denotes a survival distribution. Should have
#' the same variable names as \code{\link{params_mstate_nma}}.
#' @param check_covs Logical indicating whether to check that all covariates in 
#' \code{data} are contained in \code{params}.
#' @param covs If \code{check_covs} is \code{TRUE}, then \code{data_covs}
#' cannot be \code{NULL} and must specify all of the covariates in \code{data}
#' that should be contained in \code{params}.
#' @details The "dist" attribute from \code{data} is used to select a survival
#' distribution from the \code{mstate_nma} element contained in \code{params}. The 
#' covariates for the selected survival distribution in \code{mstate_nma} 
#' that are also contained in \code{data} are extracted.
#' @return  A \code{\link[hesim]{params_surv}} objects from the 
#' \href{https://innovationvalueinitiative.github.io/hesim/}{hesim} package.
#' @seealso \code{\link{create_transmod_data}}
#' @examples
#' # Joint distribution of parameters
#' # Treatment sequences
#' txseq1 <- txseq(first = "erlotinib",
#'                 second = c("osimertinib", "PBDC"),
#'                 second_plus = c("PBDC + bevacizumab", "PBDC + bevacizumab"))
#' txseq2 <- txseq(first = "gefitinib",
#'                 second = c("osimertinib", "PBDC"),
#'                 second_plus = c("PBDC + bevacizumab", "PBDC + bevacizumab"))
#' txseqs <- txseq_list(seq1 = txseq1, seq2 = txseq2)
#'
#' # Patient population
#' pats <- create_patients(n = 2)
#'
#' # Model structure
#' struct <- model_structure(txseqs, dist = "weibull")
#' tmat <- create_trans_mat(struct)
#'
#' # Data and parameters for state transition model
#' transmod_data <- create_transmod_data(struct, tmat, pats, mutation_prob = .45)
#' transmod_params <- create_transmod_params(n = 2, transmod_data)
#' print(transmod_params)
#' @export
create_transmod_params <- function(n = 100,
                                   data,
                                   params_mstate_nma = iviNSCLC::params_mstate_nma,
                                   check_covs = FALSE,
                                   covs = NULL){
  if(!inherits(data, "expanded_hesim_data")){
    stop("'data' must be an object of class 'expanded_hesim_data'.")
  }
  if (check_covs){
    if (is.null(covs)){
      stop("If 'check_covs' = TRUE, then 'covs' cannot be NULL.")
    }
  }
  dist <- attributes(data)$dist
  for (i in 1:length(params_mstate_nma)){
    if(!inherits(params_mstate_nma[[i]], "params_surv")){
      stop("Each element of 'params_mstate_nma' must be an object of class 'params_surv'.")
    }
  }
  params <- params_mstate_nma[[dist]]
  params <- sample_params_mstate_nma(n, params)
  n_params <- length(params$coefs)
  for (i in 1:n_params){
    if (check_covs){
      param_name <- names(params$coefs)[i]
      vars <- colnames(params$coefs[[i]])
      covs_i <- covs[grep(param_name, covs)]
      if (!all(covs_i %in% vars)){
        stop("All variables from 'data' are non contained in 'params'.")
      }  
    }
    col_inds <- which(colnames(params$coefs[[i]]) %in% colnames(data))
    params$coefs[[i]] <- params$coefs[[i]][, col_inds]
  }
  return(params)
}