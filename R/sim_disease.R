#' Create a line-specific transition matrix
#' 
#' Create a transition matrix in the same format as the \link[mstate]{mstate} 
#' conditional on the line of treatment.
#' @param line The line of treatment.
#' @return A matrix where rows/columns denote states and entries in the matrix
#' denote transitions between states.
#' @examples
#' txseq1 <- txseq(first = "erlotinib",
#'                 second = c("osimertinib", "PBDC"),
#'                 second_plus = c("PBDC + bevacizumab", "PBDC + bevacizumab"))
#' txseq2 <- txseq(first = "gefitinib",
#'                 second = c("osimertinib", "PBDC"),
#'                 second_plus = c("PBDC + bevacizumab", "PBDC + bevacizumab")) 
#' txseqs <- txseq_list(seq1 = txseq1, seq2 = txseq2)
#' tmat <- create_trans_mat(line = "2")
#' class(tmat)
#' print(tmat)
#' @export
create_trans_mat <- function(line = c("1", "2", "2+")){
  line <- match.arg(line)
  tmat <- matrix(NA, nrow = 5, ncol = 5) 
  colnames(tmat) <- rownames(tmat) <- c("S", "P1", "P2", "P2+", "D")
  if (line == "1"){
    tmat[1, ] <- c(NA, 1, NA, NA, 2)
  } else if (line == "2"){
    tmat[2, ] <- c(NA, NA, 1, NA, 2)
  } else if (line == "2+"){
    tmat[3, ] <- c(NA, NA, NA, 1, 2)
    tmat[4, ] <- c(NA, NA, NA, NA, 3)
  } else{
    stop("'line' must be 1, 2, or 2+",
         call. = FALSE)
  }
  return(tmat)
}

#' Create treatment strategies object
#'
#' Create a \code{\link{data.table}} containing information on treatment strategies used
#' at a particular treatment line. Of the same format as the \code{stratgies}
#' element in \link[hesim]{hesim_data} in the \code{hesim} package.
#' @param object The treatment sequences of interest. Must be objects of class
#' \code{\link{txseq_list}}.
#' @param line The line of treatment.
#' @param mutation If the line of treatmnet is > 1, whether the patient has
#' a T790M mutation.
#' @return A \code{\link{data.table}}.
#' @export
create_strategies <- function(object, line = c("1", "2", "2+"), 
                                         mutation = c("neg", "pos")){
  line <- switch(line,
                 "1" = "first",
                 "2" = "second",
                 "2+" = "second_plus")
  mutation <- match.arg(mutation)
  

  # Strategy ID and name
  if (line == "first"){
      name <- sapply(object, function (x) x[[line]])
  } else{
      name <- sapply(object, function (x) x[[line]])[mutation, ]
  }
  strategy_id <- seq_len(length(name))
  strategy_tbl <- data.table(strategy_id, name)

  # Treatment abbreviation
  abb_pos <- match(name, treatments$name)
  abb <- treatments$abb[abb_pos]
  unique_abb <- unique(abb)
  
  # Treatment effect variables
  d_vars_mat <- matrix(NA, nrow = length(abb), ncol = length(unique_abb))
  colnames(d_vars_mat) <- paste0("d_", unique_abb)
  for (i in 1:length(unique_abb)){
    d_vars_mat[, i] <- ifelse(abb == unique_abb[i], 1, 0)
  }
  strategy_tbl <- cbind(strategy_tbl, data.table(abb = abb, d_vars_mat))

  # Return
  return(strategy_tbl)

}



#' Simulate disease progression
#' 
#' Simulate disease progression for EGFR+ non-small cell lung cancer patients
#' given treatment sequences of interest. 
#' @param txseqs The treatment sequences of interest. Must be objects of class
#' \code{\link{txseq_list}}. 
#' @param params The posterior distribution of the parameters of the disease model. Must be an object of class
#' \code{\link[hesim]{params_surv}} from the \code{hesim} package.
#' @param n_patients The number of patients to model.
#' @return An object of class "indiv_ctstm_disprog". See the description for \code{\link[hesim]{IndivCtstmTrans}} 
#' from the \code{hesim} package.  
#' @export
sim_disease <- function(txseqs, params, n_patients){
  # Step 1: create_strategies(txseq, line, mutation)
  
  # Step 2: create_patients(n_patients)
  
  # Step 3: hesim_data(strategies, patients) -> expand_hesim_data()
  
  # Step 4: create_trans_mat(txseqs, line)
  
  # Step 5: extract_params(params)
  
  # Step 6: IndivCtstmTrans$new(data, params, tmat, start_state,
  #                             start_time, start_age)
  
  # Step 7: sim_disease()
  
  return(2)
}