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
