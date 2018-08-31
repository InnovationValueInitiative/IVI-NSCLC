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

#' Create a line-specific transition matrix
#' 
#' Create a transition matrix in the same format as the \link[mstate]{mstate} 
#' conditional on the line of treatment.
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
#' tmat <- create_trans_mat(txseqs)
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
#' Create a list of data tables containing information on treatment strategies where 
#' each table denotes a different treatment line. Of the same format as the \code{stratgies}
#' element in \link[hesim]{hesim_data} in the \code{hesim} package. 
#' @param txseq_list The treatment sequences of interest. Must be objects of class
#' \code{\link{txseq_list}}. 
#' @return A list of \code{\link{data.table}}
#' @export
create_strategies <- function(txseq_list ){
  return(2)
}