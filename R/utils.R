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