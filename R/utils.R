#' Input validation for class objects
#' 
#' \code{check} is a generic function for validating the inputs of class objects.
#' @param object object to check.
#' @param ... Further arguments passed to or from other methods. Currently unused.
#' 
#' @return If validation is successful, returns the object in question; otherwise,
#' informs the user that an error has occurred.  
#' @keywords internal
check <- function (object, ...) {
  UseMethod("check")
}

check_is_class <- function(object, name, class){
  if (!inherits(object, class)){
      stop(paste0("'", name, "' must be of class '", class, "'"),
         call. = FALSE)
  }  
}

#' Tidy data
#' 
#' A generic function for creating tidy data.
#' @param object An object to create tidy data from.
#' @param ... Further arguments passed to or from other methods. 
#' @export
#' @keywords internal
#' @seealso \code{\link{tidy.ae_probs}}
tidy <- function(object, ...){
  UseMethod("tidy", object)
}

tx_by_state <- function(struct){
  line <- state_id <- NULL
  
  n_txseqs <- length(struct$txseqs)
  txseq_dt <- vector(mode = "list", length = n_txseqs)
  for (i in 1:n_txseqs){
    txseq_i <- unlist(struct$txseqs[[i]])
    txseq_dt[[i]] <- data.table(strategy_id = i,
                                tx_name = c(txseq_i[1], txseq_i),
                                line = rep(c("first", "second", "second_plus"), 
                                           each = 2),
                                mutation = rep(c(0, 1), 3),
                                grp_id = rep(c(1, 2), 3)) 
  }
  txseq_dt <- rbindlist(txseq_dt)
  if (attr(struct$txseqs, "start_line") == "second"){
    txseq_dt <- txseq_dt[line != "first"]
  }
  txseq_dt[, state_id := .GRP, by = "line"]
  return(txseq_dt[,])
}