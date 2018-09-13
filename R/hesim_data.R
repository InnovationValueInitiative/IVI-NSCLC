#' Create health states table
#' 
#' Create a \code{\link{data.table}} containing information on health states for
#' use in \code{\link[hesim]{hesim_data}} as a function of whether the patient
#' is beginning first or second line treatment..
#' @param start_line The line of treatment that a patient is starting.
#' @return An object of class "states", which is a \code{\link{data.table}} with 
#' the following columns:
#' \describe{
#' \item{state_id}{The state ID number.}
#' \item{state_name}{The state name,}
#' \item{state_name_long}{A long-form state name.}
#' }
#' @examples
#' create_states(start_line = "first")
#' create_states(start_line = "second")
#' @export
create_states <- function(start_line = c("first", "second")){
  start_line <- match.arg(start_line)
  if (start_line == "first"){
    state_names <- c("S1", "P1 -> S2", "P2", "D") 
    state_names_long <- c("Stable with 1L", "Progression with 1L -> Stable with 2L",
                        "Progression with 2L", "Death")    
  } else {
    state_names <- c("P1 -> S2", "P2", "D") 
    state_names_long <- c("Progression with 1L -> Stable with 2L",
                        "Progression with 2L", "Death")    
  }
  state_id <- seq_along(state_names)
  return(data.table(state_id = state_id,
                    state_name = state_names,
                    state_name_long = state_names_long))
}