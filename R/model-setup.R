#' Create health states
#' 
#' Generic function to create data table of health states.
#' @param object Object used to create health states.
#' @return The form of the returned value depends on \code{object}.
#' @keywords internal
#' @seealso \code{\link{create_states.txseq_list}}
#' @export
create_states <- function (object) {
  UseMethod("create_states")
}

create_states_txseq <- function(object, start_line){
  if (start_line == "first"){
    if (object$first != "osimertinib"){
      state_names <- pkg_env$state_names_start1L_7
      state_names_long <- pkg_env$state_names_long_start1L_7
    } else {
      state_names <- pkg_env$state_names_start1L_3
      state_names_long <- pkg_env$state_names_long_start1L_3  
    }
  } else { # Second line
    state_names <- pkg_env$state_names_start2L_3
    state_names_long <- pkg_env$state_names_long_start2L_3 
  }
  state_id <- seq_along(state_names)
  return(data.table(state_id = state_id,
                    state_name = state_names,
                    state_name_long = state_names_long)) 
}

#' Create health states tables
#' 
#' Create a list of data tables describing the health states for the model
#' structure associated with each modeled treatment sequence. 
#' @param object A \code{\link{txseq_list}} object.
#' @return
#' \describe{
#' \item{state_id}{The state ID number.}
#' \item{state_name}{The state name,}
#' \item{state_name_long}{A long-form state name.}
#' }
#' @examples
#' txseq1 <- txseq(first = "erlotinib",
#'                 second = c("osimertinib", "PBDC"),
#'                 second_plus = c("PBDC + bevacizumab", "PBDC + bevacizumab"))
#' txseq2 <- txseq(first = "osimertinib",
#'                 second = c("PBDC", "PBDC"),
#'                 second_plus = c("PBDC + bevacizumab", "PBDC + bevacizumab"))
#' txseqs <- txseq_list(seq1 = txseq1, seq2 = txseq2) 
#' create_states(txseqs)
#' txseqs <- txseq_list(seq1 = txseq1, seq2 = txseq2, start_line = "second") 
#' create_states(txseqs)
#' @export
create_states.txseq_list <- function(object){
  start_line <- attributes(object)$start_line
  l <- vector(mode = "list", length = length(object))
  names(l) <- names(object)
  for (i in 1:length(l)){
    l[[i]] <- create_states_txseq(object[[i]], start_line)
  }
  class(l) <- "states"
  return(l)
}

#' Create transition matrcies
#' 
#' Generic function to create transition matrices.
#' @param object Object used to create transition matrices.
#' @return The form of the returned value depends on \code{object}.
#' @keywords internal
#' @seealso \code{\link{create_trans_mats.txseq_list}}
#' @export
create_trans_mats <- function (object) {
  UseMethod("create_trans_mats")
}

create_trans_mat <- function(object, start_line){
  if (start_line == "first"){
    if (object$first != "osimertinib"){  
      tmat <- rbind(c(NA, 1, NA, NA, NA, NA, 2),
                    c(NA, NA, 3, NA, 4, NA, NA),
                    c(NA, NA, NA, 5, NA, NA, 6),
                    c(NA, NA, NA, NA, NA, NA, 7),
                    c(NA, NA, NA, NA, NA, 8, 9),
                    c(NA, NA, NA, NA, NA, NA, 10),
                    c(NA, NA, NA, NA, NA, NA, NA))    
      colnames(tmat) <- rownames(tmat) <- pkg_env$state_names_start1L_7        
    } else {
      tmat <- rbind(c(NA, 1, 2),
                    c(NA, NA, 3),
                    c(NA, NA, NA))
      colnames(tmat) <- rownames(tmat) <- pkg_env$state_names_start1L_3     
    }
  } else{
      tmat <- rbind(c(NA, 1, 2),
                    c(NA, NA, 3),
                    c(NA, NA, NA))
      colnames(tmat) <- rownames(tmat) <- pkg_env$state_names_start2L_3      
  }
  return(tmat)
}

#' Create transition matrices
#' 
#' Create a list of transition matrices describing patient transitions between 
#' health states for the model associated with each modeled treatment sequence. 
#' @param object A \code{\link{txseq_list}} object.
#' @return A list of transition matrices for each treatment sequence, each 
#' in the same format as in the \link[mstate]{mstate} package.
#' @examples
#' txseq1 <- txseq(first = "erlotinib",
#'                 second = c("osimertinib", "PBDC"),
#'                 second_plus = c("PBDC + bevacizumab", "PBDC + bevacizumab"))
#' txseq2 <- txseq(first = "osimertinib",
#'                 second = c("PBDC", "PBDC"),
#'                 second_plus = c("PBDC + bevacizumab", "PBDC + bevacizumab"))
#' txseqs <- txseq_list(seq1 = txseq1, seq2 = txseq2) 
#' create_trans_mats(txseqs)
#' txseqs <- txseq_list(seq1 = txseq1, seq2 = txseq2, start_line = "second") 
#' create_trans_mats(txseqs)
#' @export
create_trans_mats.txseq_list <- function(object){
  start_line <- attributes(object)$start_line
  l <- vector(mode = "list", length = length(object))
  names(l) <- names(object)
  for (i in 1:length(l)){
    l[[i]] <- create_trans_mat(object[[i]], start_line)
  }
  class(l) <- "trans_mats"
  return(l)
}

#' Create patient data table
#' 
#' Create a data table of patients to model.
#' @param n Number of patients to model.
#' @examples
#' create_patients(n = 10)
#' @export
create_patients <- function(n){
  patient_id <- 1:n
  return(data.table(patient_id = patient_id))
}

#' #' @export
#' create_strategies <- function (object, ...) {
#'   UseMethod("create_strategies")
#' }
#' 
#' #' Create treatment strategies 
#' #' 
#' #' Create a data table of treatment strategies as a function of the 
#' #' modeled treatment sequences.
#' #' @param object A \code{\link{txseqs_list}} object.
#' #' @examples
#' #' txseq1 <- txseq(first = "erlotinib",
#' #'                 second = c("osimertinib", "PBDC"),
#' #'                 second_plus = c("PBDC + bevacizumab", "PBDC + bevacizumab"))
#' #' txseq2 <- txseq(first = "osimertinib",
#' #'                 second = c("PBDC", "PBDC"),
#' #'                 second_plus = c("PBDC + bevacizumab", "PBDC + bevacizumab"))
#' #' txseqs <- txseq_list(seq1 = txseq1, seq2 = txseq2)  
#' #' create_strategies(txseqs, start_line = "first")
#' #' @export
#' create_strategies.txseq_list <- function(object){
#'   start_line <- attributes(object)$start_line
#'   n_strategies <- length(object)
#'   strategy_id <- 1:n_strategies
#'   container <- rep(NA, length(object))
#'   name_first <- container
#'   name_second_pos <- name_second_neg <- container
#'   name_second_plus_pos <- name_second_plus_neg <- container
#'   for (i in 1:n_strategies){
#'     name_first[i] <- object[[i]]$first
#'     name_second_pos[i] <- object[[i]]$second["pos"]
#'     name_second_neg[i] <- object[[i]]$second["neg"]
#'     name_second_plus_pos[i] <- object[[i]]$second_plus["pos"]    
#'     name_second_plus_neg[i] <- object[[i]]$second_plus["neg"]    
#'   }
#'   strategies <- data.table(strategy_id = strategy_id,
#'                            name_first = name_first,
#'                            name_second_pos = name_second_pos,
#'                            name_second_neg = name_second_neg,
#'                            name_second_plus_pos = name_second_plus_pos,
#'                            name_second_plus_neg = name_second_plus_neg)
#'   return(strategies)
#' }
