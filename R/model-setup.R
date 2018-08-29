#' Create health states table
#' 
#' Create a \code{\link{data.table}} containing information on health states for
#' use in \code{\link[hesim]{hesim_data}}.
#' @return An object of class "states", which is a \code{\link{data.table}} with 
#' the following columns:
#' \describe{
#' \item{state_id}{The state ID number.}
#' \item{state_name}{The state name,}
#' \item{state_name_long}{A long-form state name.}
#' }
#' @examples
#' txseq1 <- txseq(first = "erlotinib",
#'                 second = c("osimertinib", "PBDC"),
#'                 second_plus = c("PBDC + bevacizumab", "PBDC + bevacizumab"))
#' txseq2 <- txseq(first = "gefitinib",
#'                 second = c("osimertinib", "PBDC"),
#'                 second_plus = c("PBDC + bevacizumab", "PBDC + bevacizumab")) 
#' txseqs <- txseq_list(seq1 = txseq1, seq2 = txseq2)
#' states <- create_states(txseqs) 
#' class(states)
#' print(states)
#' @export
create_states <- function(object){
  UseMethod("create_states", object)
}

#' @rdname create_states
#' @export
create_states.txseq_list <- function(object){
  len <- length(object[[1]])
  n_states <- len + 2
  
  # Long progression name
  str_end <- paste0(seq_len(len), "L")
  str_end <- ifelse(str_end == "3L", "2L+", str_end)
  prog_states_long <- paste0("Progression on ", str_end)
  
  # Progression name
  prog_states <- paste0("P", seq_len(len))
  prog_states <- ifelse(prog_states == "P3", "P2+", prog_states)
  
  # Create table
  state_id <- seq_len(n_states)
  state_names <- c("S", prog_states, "D")
  state_names_long <- c("Stable", prog_states_long, "Death")
  tbl <- data.table(state_id = state_id,
                    state_name = state_names,
                    state_name_long = state_names_long)
  setattr(tbl, "class", c("states", "data.table", "data.frame"))
  return(tbl)
}

#' Create line-specific transitions matrices
#' 
#' Create transition matrices in the same format as the \link[mstate]{mstate} 
#' package for each treatment line.
#' @return An object of class "trans_mats", which is a list of matrices with length 
#' equal to the number of treatment lines.
#' @examples
#' txseq1 <- txseq(first = "erlotinib",
#'                 second = c("osimertinib", "PBDC"),
#'                 second_plus = c("PBDC + bevacizumab", "PBDC + bevacizumab"))
#' txseq2 <- txseq(first = "gefitinib",
#'                 second = c("osimertinib", "PBDC"),
#'                 second_plus = c("PBDC + bevacizumab", "PBDC + bevacizumab")) 
#' txseqs <- txseq_list(seq1 = txseq1, seq2 = txseq2)
#' trans_mats <- create_trans_mats(txseqs)
#' class(trans_mats)
#' print(trans_mats)
#' @export
create_trans_mats <- function(object){
  UseMethod("create_trans_mats", object)
}

#' @rdname create_trans_mats
#' @export
create_trans_mats.txseq_list <- function(object){
  len <- length(object[[1]])
  n_states <- len + 2
  tmat <- matrix(NA, nrow = n_states, ncol = n_states) 
  colnames(tmat) <- rownames(tmat) <- c("S", "P1", "P2", "P2+", "D")
  l <- vector(mode = "list", length = len)
  for (i in 1:length(l)){
    l[[i]] <- tmat
    if (i == 1){
      l[[i]][1, ] <- c(NA, 1, NA, NA, 2)
    } else if (i == 2){
      l[[i]][2, ] <- c(NA, NA, 1, NA, 2)
    } else if (i == 3){
      l[[i]][3, ] <- c(NA, NA, NA, 1, 2)
      l[[i]][4, ] <- c(NA, NA, NA, NA, 3)
    }
  }
  class(l) <- c("trans_mats")
  return(l)
}