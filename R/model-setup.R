#' Model structure
#' 
#' Define the structure of the oncology model.
#' @param n_states Number of modeled health states.
#' @param start_line The starting line of treatmnet that is being modeled. When
#' modeling second line treatment, the first line must be specified
#' in order to characterize a treatment history. 
#' @return A list
#' @examples
#' model_structure()
#' model_structure(n_states = "three", start_line = "first")
#' model_structure(n_states = "three", start_line = "second")
#' @export
model_structure <- function(n_states = c("four", "three"), 
                            start_line = c("first", "second")) {
  start_line <- match.arg(start_line)
  n_states <- match.arg(n_states)
  if (n_states == "four" & start_line == "second"){
    stop("If start_line == 'second', then n_states must be 'three'.",
         call. = FALSE)
  }
  l <- list(n_states = n_states,
            start_line = start_line)
  class(l) <- "model_structure"
  return(l)
}

check_is_model_structure <- function(object){
  if (!inherits(object, "model_structure")){
      stop("'object' must be of class 'model_stucture'",
         call. = FALSE)
  }  
}

#' Create health states tables
#' 
#' Create a data table describing the health states for the model
#' structure. 
#' @param object A \code{\link{model_structure}} object.
#' @return
#' \describe{
#' \item{state_id}{The state ID number.}
#' \item{state_name}{The state name,}
#' \item{state_name_long}{A long-form state name.}
#' }
#' @examples
#' struct <- model_structure(n_states = "four", start_line = "first")
#' create_states(struct)
#' 
#' struct <- model_structure(n_states = "three", start_line = "first")
#' create_states(struct)
#' 
#' struct <- model_structure(n_states = "three", start_line = "second")
#' create_states(struct)
#' @export
create_states <- function(object){
  check_is_model_structure(object)
  start_line <- object$start_line
  if (start_line == "first"){
    if (object$n_states == "four"){
      state_names <- pkg_env$state_names_start1L_4
      state_names_long <- pkg_env$state_names_long_start1L_4
    } else {
      state_names <- pkg_env$state_names_start1L_3
      state_names_long <- pkg_env$state_names_long_start1L_3  
    }
  } else { # Second line
    state_names <- pkg_env$state_names_start2L_3
    state_names_long <- pkg_env$state_names_long_start2L_3 
  }  
  state_id <- seq_along(state_names)
  return(data.table(state_id = state_id,
                    state_name = state_names,
                    state_name_long = state_names_long))   
}

#' Create transition matrix
#' 
#' Create a transition matrix describing patient transitions between 
#' health states.
#' @param object A \code{\link{model_structure}} object.
#' @return A transition matrix of the same format as in the \link[mstate]{mstate} 
#' package.
#' @examples
#' struct <- model_structure(n_states = "four", start_line = "first")
#' create_trans_mat(struct)
#' 
#' struct <- model_structure(n_states = "three", start_line = "first")
#' create_trans_mat(struct)
#' 
#' struct <- model_structure(n_states = "three", start_line = "second")
#' create_trans_mat(struct)
#' @export
create_trans_mat <- function(object){
  check_is_model_structure(object)
  if (object$start_line == "first"){
    if (object$n_states == "four"){  
      tmat <- rbind(c(NA, 1, NA, 2),
                    c(NA, NA, 3, 4),
                    c(NA, NA, NA, 5),
                    c(NA, NA, NA, NA))
      colnames(tmat) <- rownames(tmat) <- pkg_env$state_names_start1L_4        
    } else {
      tmat <- rbind(c(NA, 1, 2),
                    c(NA, NA, 3),
                    c(NA, NA, NA))
      colnames(tmat) <- rownames(tmat) <- pkg_env$state_names_start1L_3     
    }
  } else{
      tmat <- rbind(c(NA, 1, 2),
                    c(NA, NA, 3),
                    c(NA, NA, NA))
      colnames(tmat) <- rownames(tmat) <- pkg_env$state_names_start2L_3      
  }
  return(tmat)
}

#' Create patient data table
#' 
#' Create a data table of patients to model.
#' @param n Number of patients to model.
#' @examples
#' create_patients(n = 10)
#' @export
create_patients <- function(n){
  patient_id <- 1:n
  return(data.table(patient_id = patient_id))
}
