#' Package environment
#'
#' An environment to store global variables for the package. 
pkg_env <- new.env(parent = emptyenv())
pkg_env$tki_1gen <- c("erlotinib", "gefitinib")
pkg_env$tki_2gen <- c("afatinib", "dacomitinib")
pkg_env$anti_vegf <- "bevacizumab"
pkg_env$ici <- c("nivolumab", "pembrolizumab", "atezolizumab")

