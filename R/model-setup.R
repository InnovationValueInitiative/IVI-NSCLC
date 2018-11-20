#' Model structure
#' 
#' Define the structure of the oncology model.
#' @param txseqs A \code{\link{txseq_list}} object.
#' @param n_states Number of modeled health states.
#' @param dist Parametric distribution used to model health state transitions.
#' Options are \code{"weibull"} (Weibull), \code{"gompertz"} (gompertz),
#' \code{"fracpoly1"}(2nd order fractional polynomial with \eqn{p_1 = 0} and \eqn{p_2 = 0}),
#' and \code{"fracpoly2"}(2nd order fractional polynomial with \eqn{p_1 = 0} and \eqn{p_2 = 1}).
#' @return A list containing the elements \code{txseqs}, \code{n_states} and
#' \code{dist}.
#' @examples
#' txseq1 <- txseq(first = "erlotinib",
#'                 second = c("osimertinib", "PBDC"),
#'                second_plus = c("PBDC + bevacizumab", "PBDC + bevacizumab"))
#' txseq2 <- txseq(first = "gefitinib",
#'                second = c("osimertinib", "PBDC"),
#'                second_plus = c("PBDC + bevacizumab", "PBDC + bevacizumab")) 
#' txseqs <- txseq_list(seq1 = txseq1, seq2 = txseq2)
#' 
#' # Model with 4 health states
#' struct <- model_structure(txseqs)
#' names(struct)
#' class(struct$txseqs)
#' struct$n_states
#' struct$dist
#' 
#' # Model with 3 health states
#' struct <- model_structure(txseqs, n_states = "three", dist = "weibull")
#' struct$n_states
#' @export
model_structure <- function(txseqs,
                            n_states = c("four", "three"),
                            dist = c("weibull", "gompertz", "fracpoly1",
                                     "fracpoly2")) {
  dist <- match.arg(dist)
  if (dist != "weibull") dist <- "weibull" # NEED TO CHANGE THIS!
  n_states <- match.arg(n_states)
  check_is_class(txseqs, "txseqs", "txseq_list")
  if (n_states == "four" & attributes(txseqs)$start_line == "second"){
    stop("If the model starts at second line, then n_states must be 'three'.",
         call. = FALSE)
  }  
  if (n_states == "four" & "osimertinib" %in% sapply(txseqs, function (x) x$first)){
      msg <- paste0("There is no evidence to parameterize a model with four health states ", 
                    "when osimertinib is used as a first line treatment. ",
                    "Use a model with three health states for sequences beginning ",
                    "with osimertinib instead.")
      stop(msg, call. = FALSE) 
  }    
  l <- list(txseqs = txseqs,
            n_states = n_states,
            dist = dist)
  class(l) <- "model_structure"
  return(l)
}

#' Create health states tables
#' 
#' Create a data table describing the health states for the model.
#' @param object A \code{\link{model_structure}} object.
#' @return
#' A data table with the following columns:
#' \describe{
#' \item{state_id}{The state ID number.}
#' \item{state_name}{The state name,}
#' \item{state_name_long}{A long-form state name.}
#' }
#' @examples
#' txseq1 <- txseq(first = "erlotinib",
#'                second = c("osimertinib", "PBDC"),
#'                second_plus = c("PBDC + bevacizumab", "PBDC + bevacizumab"))
#' txseq2 <- txseq(first = "gefitinib",
#'                second = c("osimertinib", "PBDC"),
#'                second_plus = c("PBDC + bevacizumab", "PBDC + bevacizumab")) 
#' txseqs <- txseq_list(seq1 = txseq1, seq2 = txseq2) 
#' 
#' struct1 <- model_structure(txseqs, n_states = "four")
#' create_states(struct1)
#' 
#' struct2 <- model_structure(txseqs, n_states = "three")
#' create_states(struct2)
#' 
#' txseqs <- txseq_list(seq1 = txseq1, seq2 = txseq2, 
#'                      start_line = "second", mutation = "negative")
#' struct3 <- model_structure(txseqs, n_states = "three")
#' create_states(struct3)
#' @export
create_states <- function(object){
  check_is_class(object, "object", "model_structure")
  start_line <- attributes(object$txseqs)$start_line
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
#' txseq1 <- txseq(first = "erlotinib",
#'                second = c("osimertinib", "PBDC"),
#'                second_plus = c("PBDC + bevacizumab", "PBDC + bevacizumab"))
#' txseq2 <- txseq(first = "gefitinib",
#'                second = c("osimertinib", "PBDC"),
#'                second_plus = c("PBDC + bevacizumab", "PBDC + bevacizumab")) 
#' txseqs <- txseq_list(seq1 = txseq1, seq2 = txseq2) 
#' 
#' struct1 <- model_structure(txseqs, n_states = "four")
#' create_trans_mat(struct1)
#'  
#' struct2 <- model_structure(txseqs, n_states = "three")
#' create_trans_mat(struct2)
#'  
#' txseqs <- txseq_list(seq1 = txseq1, seq2 = txseq2, 
#'                      start_line = "second", mutation = "positive")
#' struct3 <- model_structure(txseqs, n_states = "three")
#' create_trans_mat(struct3)
#' @export
create_trans_mat <- function(object){
  check_is_class(object, "object", "model_structure")
  if (attributes(object$txseqs)$start_line == "first"){
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
#' @param female_prop The proportion of patients that are female.
#' @param age_mean Mean age. Based on sources cited in \code{\link{age_dist}}.
#' @param age_sd Standard deviation of age. Based on sources cited in \code{\link{age_dist}}.
#' @param mutation_prob The probability of a T790M mutation. The default value
#' is based on Table 3 from the article by Ma et al. cited below.
#' @examples
#' create_patients(n = 10)
#' @references 
#' Ma C, Wei S, Song Y. T790M and acquired resistance of EGFR TKI: a literature 
#' review of clinical reports. Journal of thoracic disease. 2011 Mar;3(1):10.
#' @return An object of class "patients", which is a \code{data.table} 
#' containing each modeled patient. Columns are:
#' \describe{
#' \item{patient_id}{An integer from 1 to \code{n} denoting a unique patient.}
#' \item{mutation}{1 if a patient has a T790M mutation and 0 otherwise.}
#' \item{female}{1 if a patient is female and 0 otherwise.}
#' }
#' 
#' @export
create_patients <- function(n, female_prop = .45, 
                            age_mean = 70.39, age_sd = 11.68,
                            mutation_prob = .52){
  patient_id <- 1:n
  
  # Age
  age <- truncnorm::rtruncnorm(n, a = 0, b = 100, mean = age_mean, sd = age_sd)
  
  # Gender
  female <- stats::rbinom(n, 1, female_prop)
  
  # Mutations
  mutation <- stats::rbinom(n, 1, mutation_prob)
  
  # Create dataset
  object <- data.table(patient_id = patient_id,
                       female = female,
                       age = age,
                       mutation = mutation)
  setattr(object, "class", c("patients", "data.table", "data.frame"))
  return(object)
}
