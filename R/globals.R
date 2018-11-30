#' Package environment
#'
#' An environment to store global variables for the package. 
pkg_env <- new.env(parent = emptyenv())

# Treatements
pkg_env$tki_1gen <- c("erlotinib", "gefitinib")
pkg_env$tki_2gen <- c("afatinib", "dacomitinib")
pkg_env$anti_vegf <- "bevacizumab"
pkg_env$ici <- c("nivolumab", "pembrolizumab", "atezolizumab")

# State names
## Start at first line
### Four state model
pkg_env$state_names_start1L_4 <- c("S1", "P1/S2", "P2", "D")
pkg_env$state_names_long_start1L_4 <- c("Stable with 1L",
                                        "Progression with 1L -> Stable with 2L",
                                        "Progression with 2L",
                                        "Death")

### Three state model
pkg_env$state_names_start1L_3 <- c("S1", "P1/S2", "D") 
pkg_env$state_names_long_start1L_3 <- c("Stable with 1L",
                                        "Progression with 1L",
                                        "Death")  

## Start at second line
pkg_env$state_names_start2L_3 <- c("P1/S2", "P2", "D") 
pkg_env$state_names_long_start2L_3 <- c("Stable with 2L",
                                        "Progression with 2L", 
                                        "Death")  

