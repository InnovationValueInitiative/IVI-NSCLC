#' Adverse event probabilities
#' 
#' Obtain adverse event probabilities for treatment sequences given existing 
#' probabilities by treatment from a network meta-analysis.
#' @param n The number of random observations of the parameters to draw.
#' @param struct A \code{\link{model_structure}} object.
#' @param params_ae Parameter estimates of the probabilities of adverse 
#' events in the same format as \code{\link{params_ae_nma}}.  
#' @return A list of matrices where each matrix corresponds to a type of adverse
#' event. Rows of each matrix are posterior samples of the parameters and columns
#' are treatment strategies. The list is an object of class "ae_probs" with an
#' attribute "tx_abb" denoting the treatment used to estimate adverse events for
#'  each treatment sequence. 
#' @examples
#' txseq1 <- txseq(first = "erlotinib",
#'                 second = c("osimertinib", "PBDC"),
#'                 second_plus = c("PBDC + bevacizumab", "PBDC + bevacizumab"))
#' txseq2 <- txseq(first = "gefitinib",
#'                 second = c("osimertinib", "PBDC"),
#'                 second_plus = c("PBDC + bevacizumab", "PBDC + bevacizumab"))
#' txseqs <- txseq_list(seq1 = txseq1, seq2 = txseq2) 
#' struct <- model_structure(txseqs)
#' ae_probs <- ae_probs(n = 3, struct = struct)
#' print(ae_probs[1:3])
#' tidy_ae_probs <- tidy(ae_probs)
#' head(tidy_ae_probs)
#' @export
ae_probs <- function(n, struct, params_ae = iviNSCLC::params_ae_nma){
  check_is_class(struct, name = "struct", class = "model_structure") 
  
  # Treatment used for AE
  if (attr(struct$txseqs, "start_line") == "first"){
    tx <- sapply(struct$txseqs, function (x) x$first)
  } else if (attr(struct$txseqs, "start_line") == "second"){
      if (attr(struct$txseqs, "mutation") == "positive") {
        tx <- sapply(struct$txseqs, function (x) x$second["pos"]) 
      } else{
        tx <- sapply(struct$txseqs, function (x) x$second["neg"]) 
      }
  } else{
    stop("The starting line must either be 'first' or 'second'.")
  }
  tx_abb <- iviNSCLC::treatments$tx_abb[match(tx, iviNSCLC::treatments$tx_name)]  
  
  # Select columns and rows
  prob_vars <- paste0("prob_", tx_abb)
  for (i in 1:length(params_ae)){
    ## Sample rows for PSA based on posterior samples 
    n_samples <- nrow(params_ae[[i]])
    if (n <= n_samples){
      sampled_rows <- sample.int(n_samples, n, replace = FALSE) 
    } else if (n > n_samples) {
      warning("'n' is larger than the values of 'n_samples' in 'params_mstate_nma'.")
      sampled_rows <- sample.int(n_samples, n, replace = TRUE) 
    } 
    params_ae[[i]] <-  params_ae[[i]][sampled_rows, prob_vars, drop = FALSE] 
    colnames(params_ae[[i]]) <- names(struct$txseqs)
  }
  class(params_ae) <- "ae_probs"
  attr(params_ae, "tx_abb") <- tx_abb
  return(params_ae)
}

#' Tidy adverse event data
#' 
#' Create tidy data of adverse event probabilities.
#' @param object An "ae_probs" object as returned by \code{\link{ae_probs}}.
#' @param ae_lookup A \code{data.rame} or \code{data.table} in the same format
#' as \code{\link{adverse_events}} used to lookup long-form adverse event names
#' (i.e., \code{ae_name}) given values of \code{ae_abb} contained in \code{object}.
#' @param ... Further arguments passed to or from other methods. Currently unused. 
#' @seealso \code{\link{ae_probs}}
#' @export
tidy.ae_probs <- function(object, ae_lookup = iviNSCLC::adverse_events, ...){
  n_samples <- nrow(object[[1]])
  n_strategies <- ncol(object[[1]])
  strategy_id <- 1:ncol(object[[1]])
  ae_names <-  ae_lookup[match(names(object), ae_lookup$ae_abb)]$ae_name
  n_aes <- length(ae_names)
  probs <- unlist(object)
  tbl <- data.table(strategy_id = rep(rep(strategy_id, each = n_samples),
                                      times = n_aes),
                    ae_name = rep(ae_names, each = (n_samples * n_strategies)),
                    prob = probs)
  return(tbl[,])
}