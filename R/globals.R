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
### 1st/2nd generation TKIs with 7 states
pkg_env$state_names_start1L_7 <- c("S1", "P1", "S2pos", "P2pos", "S2neg", "P2neg", "D") 
pkg_env$state_names_long_start1L_7 <- c("Stable with 1L", "Progression with 1L",
                                     "Stable with 2L (T790M positive)", 
                                     "Progression with 2L (T790M positive)",
                                     "Stable with 2L (T790M negative)", 
                                     "Progression with 2L (T790M negative)",
                                     "Death")

### Osimertinib
pkg_env$state_names_start1L_3 <- c("S1", "P1", "D") 
pkg_env$state_names_long_start1L_3 <- c("Stable with 1L",
                                          "Progression with 1L", "Death")  

## Start at second line
pkg_env$state_names_start2L_3 <- c("S2", "P2", "D") 
pkg_env$state_names_long_start2L_3 <- c("Stable with 2L",
                                     "Progression with 2L", "Death")  



pkg_env$state_names7 <- c("S1", "P1", "S2pos", "P2pos", "S2neg", "P2neg", "D") 
pkg_env$state_names_long7 <- c("Stable with 1L", "Progression with 1L",
                               "Stable with 2L (T790M positive)", 
                               "Progression with 2L (T790M positive)",
                                "Stable with 2L (T790M negative)", 
                               "Progression with 2L (T790M negative)",
                               "Death")  
pkg_env$state_names3 <- c("S2", "P2", "D") 
pkg_env$state_names_long3 <- c("Stable with 2L",
                                "Progression with 2L", "Death")  

