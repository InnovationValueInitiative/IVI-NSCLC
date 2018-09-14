#' First line treatment options
#' 
#' First line treatments for EGFR+ NSCLC patients.
#' @return A vector of possible first line treatments.
#' @export
#' @examples
#' tx_1L()
tx_1L <- function(){
  treatments <- c("erlotinib",
                  "gefitinib",
                  "afatinib",
                  "dacomitinib",
                  "osimertinib")
  return(treatments)
}

#' Second line treatment options
#'
#' Possible second line treatment for EGFR+ NSCLC patients conditional on the first
#' line treatment.
#' @param first A first line treatment.
#' @return A list of two elements:
#' \describe{
#' \item{pos}{Possible treatments with a T790M mutation.}
#' \item{neg}{Possible treatments without a T790M mutation.}
#' }
#' @export
#' @examples
#' first <- tx_1L()[2]
#' print(first)
#' tx_2L(first)
tx_2L <- function(first){
  if (length(first) > 1){
    stop("Length of 'first' must be 1.",
         call. = FALSE)
  }
  tki_old <- c(pkg_env$tki_1gen, pkg_env$tki_2gen)
  tx <- list()
  if(first %in% tki_old){
    tx$pos <- "osimertinib"
    tx$neg <- c("PBDC", paste0("PBDC + ", pkg_env$anti_vegf))
    } else if (first == "osimertinib") {
      tx$pos <- tx$neg <- c("PBDC", paste0("PBDC + ", pkg_env$anti_vegf), tki_old)
    } else{
      stop(paste0(first, " is not a first line treatment."),
           call. = FALSE)
    }
  return(tx)
}

#' Second line plus treatment options
#'
#' Treatment after second line for EGFR+ NSCLC patients conditional on possible second line treatments.
#' @param second A vector of possible treatments after second line. Must be of length 2 with the 
#' first element corresponding to use with a T790M mutation and the second election corresponding 
#' to treatment without a T790M mutation.
#' @return A list of two elements:
#' \describe{
#' \item{pos}{Possible treatments with a T790M mutation.}
#' \item{neg}{Possible treatments without a T790M mutation.}
#' }
#' @export
#' @examples
#' first <- tx_1L()[3]
#' print(first)
#' second_opts <- tx_2L(first)
#' print(second_opts)
#' tx_2LP(c(second_opts$pos[1], second_opts$neg[1]))
tx_2LP <- function(second){
  if (length(second) != 2){
    stop("Length of 'second' must be 2.",
         call. = FALSE)
  }  
  tx <- list()  
  tx$pos <- tx_2LP_work(second[1])
  tx$neg <- tx_2LP_work(second[2])
  return(tx)
}


tx_2LP_work <- function(second){
  all_tx <- c("PBDC", 
            pkg_env$ici,
            paste0("PBDC + ", pkg_env$anti_vegf),
            paste0("PBDC + ", pkg_env$ici),
            paste0("PBDC + ", pkg_env$anti_vegf, " + ", pkg_env$ici)
            )
  if (second == "osimertinib"){
    tx <- all_tx
  } else if (second == "PBDC"){
    tx <- all_tx[all_tx != "PBDC"]
  } else if (second %in% paste0("PBDC + ", pkg_env$anti_vegf)){
    tx <- all_tx[!all_tx %in% c("PBDC",  paste0("PBDC + ", pkg_env$anti_vegf))]
  } else if (second %in% c(pkg_env$tki_1gen, pkg_env$tki_2gen)){
    tx <- all_tx
  } else {
    stop(paste0(second, " is not a second line treatment"),
         call. = FALSE)
  }
  return(tx)
}

#' A treatment sequence
#' 
#' Create a treatment sequence for simulation modeling. 
#' @param first The treatment to be used at first line. Must be a character of length 1.
#' @param second A vector of length 2 where the first element is the second line treatment 
#' to be used if a patient develops a T790M mutation after first line treatment and the second element is 
#' the second line treatment to be used if the patient did not.
#' @param second_plus A vector of length 2 where the first element is the post second line treatment 
#' to be used if a patient develops a T790M mutation after first line treatment and the second element is 
#' the post second line treatment to be used if the patient did not.
#' @return An object of class "txseq", which is a list of treatments used at first line (\code{first}),
#'  second line (\code{second}), and post second line (\code{second_plus}).
#' @export
#' @examples
#' txseq <- txseq(first = "erlotinib",
#'                second = c("osimertinib", "PBDC"),
#'                second_plus = c("PBDC + bevacizumab", "PBDC + bevacizumab"))
#' txseq$first
#' txseq$second
#' txseq$second_plus 
txseq <- function(first, second, second_plus){
  if(!is.character(first) | length(first) != 1){
    stop("'first' must be a character of length 1",
         call. = FALSE)
  }  
  if(!is.vector(second, mode = "character") | length(second) != 2){
    stop("'second' must be a character vector of length 2",
         call. = FALSE)
  }
  if(!is.vector(second_plus, mode = "character") | length(second_plus) != 2){
    stop("'second_plus' must be a character vector of length 2",
         call. = FALSE)
  }  
  
  object <- new_txseq(first, second, second_plus)
  return(check(object))
}

new_txseq <- function(first, second, second_plus){
  names(second) <- c("pos", "neg")
  names(second_plus) <- c("pos", "neg")
  object <- list(first = first,
             second = second,
             second_plus = second_plus)
  class(object) <- "txseq"
  return(object)
}

check.txseq <- function(x){
  
  # Check 1st line
  if(!x$first %in% tx_1L()){
    stop(paste0(x$first, " cannot be selected as a 1st line treatment."),
                call. = FALSE)
  }
  
  # Check 2nd line
  tx_2L_opts <- tx_2L(x$first)
  if (!x$second[1] %in% tx_2L_opts$pos){
    stop(paste0(x$second[1], " cannot be selected as a 2nd line treatment ",
                "for patients with a T790M mutation."),
          call. = FALSE)    
  }
  if (!x$second[2] %in% tx_2L_opts$neg){
      stop(paste0(x$second[2], " cannot be selected as a 2nd line treatment ",
                  "for patients without a T790M mutation"),
                call. = FALSE)
  }  
  
  # Check 2nd line plus
  tx_2LP_opts <- tx_2LP(x$second)  
  if (!x$second_plus[1] %in% tx_2LP_opts$pos){
    stop(paste0(x$second_plus[1], " cannot be selected as a post 2nd line treatment ",
                "for patients with a T790M mutation."),
          call. = FALSE)    
  }
  if (!x$second_plus[2] %in% tx_2LP_opts$neg){
    stop(paste0(x$second_plus[2], " cannot be selected as a post 2nd line treatment ",
                "for patients without a T790M mutation."),
          call. = FALSE)    
  }  
  
  return(x)
}

#' A list of treatment sequences
#' 
#' Create a list of objects of class "txseq". 
#' @param ... Objects to form a list.
#' @param start_line The starting line of treatmnet that is being modeled. When
#' modeling second line treatment, the first line must be specified
#' in order to characterize a treatment history.
#' @return An object of class "txseq_list", which is a list of objects of class. 
#' \code{start_line} is stored as an attribute.
#' \code{\link{txseq}}. 
#' @export
#' @examples
#' txseq1 <- txseq(first = "erlotinib",
#'                second = c("osimertinib", "PBDC"),
#'                second_plus = c("PBDC + bevacizumab", "PBDC + bevacizumab"))
#' txseq2 <- txseq(first = "gefitinib",
#'                second = c("osimertinib", "PBDC"),
#'                second_plus = c("PBDC + bevacizumab", "PBDC + bevacizumab")) 
#' txseqs <- txseq_list(seq1 = txseq1, seq2 = txseq2) 
#' class(txseqs)        
#' print(txseqs$seq1)
#' attributes(txseqs)
txseq_list <- function(..., start_line = c("first", "second")){
  start_line <- match.arg(start_line)
  objects <- new_txseq_list(..., start_line = start_line)
  return(check(objects))
}

new_txseq_list <- function(..., start_line){
  objects <- list(...)
  if(length(objects) == 1 & inherits(objects[[1]], "list")){
      objects <- objects[[1]]
  }
  class(objects) <- "txseq_list"
  attr(objects, "start_line") <- start_line
  return(objects)
}

check.txseq_list <- function(x){
  len1 <- length(x[[1]])
  for (i in 1:length(x)){
    if(!inherits(x[[i]], "txseq")){
      stop("Each element in ... must of of class 'txseq'.",
           call. = FALSE)
    }
    if (length(x[[i]]) != len1){
      stop("Each sequence must be the same length.",
           call. = FALSE)
    }
  } 
  return(x)
}