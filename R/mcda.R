#' Linear partial value function
#' 
#' Convert performance for criterion on original scale to a score on a common scale
#' using a linear partial value function.
#' @param x Performance of the criterion on the original scale.
#' @param x_min Minimum value of criterion on original scale. If higher performance
#' is better, then the minimum value is the lowest possible value; otherwise, it
#' is the highest.
#' @param x_max Maximum value of criterion on original scale. If higher performance
#' is better, then the maximum value is the highest possible value; otherwise, it
#' is the lowest.
#' @param score_min Minimum value of score on common scale.
#' @param score_max Maximum value of score on common scale.
#' @return The value of the score on the common scale.
lpvf <- function(x, x_min, x_max, score_min = 0, score_max = 100){
  score <- (x - x_min) * (score_max - score_min) / (x_max - x_min) + score_min
  if (x_min == x_max){
    score <- rep(0, length(score))
  }
  return(score)
}

 
#' Multi criteria decision analysis
#' 
#' Conduct a multi criteria decision analysis (MCDA) and compute scores for
#' competing treatment strategies using output from a probabilistic sensitivity 
#' analysis (PSA).
#' @param x A \code{data.frame} or \code{data.table} of simulation output 
#' characterizing the probability distribution of model outcomes. 
#' @param sample Character name of column from \code{x} denoting a randomly sampled 
#' parameter set from the PSA.
#' @param strategy Character name of column from \code{x} denoting treatment strategy.
#' @param criteria A vector of character names of columns from \code{x} denoting
#' the criteria to use in the MCDA.
#' @param criteria_min A vector of minimum values for each criterion. If \code{NULL}, 
#' then the minimum value is computed automatically.
#' @param criteria_max A vector of maximum values for each criterion. If \code{NULL}, 
#' then the maximum value is computed automatically.
#' @param optimal A character vector denoting whether the optimal value of each criteria
#' is \code{"low"} or \code{"high"}. If an element is \code{"low"}, then 
#' lower performance on that criterion is better, and,
#' if an element is \code{"high"}, then higher performance on that criterion is better. 
#' Must be specified if either \code{criteria_min} or 
#' \code{criteria_min} is \code{NULL}.
#' @param weights Weights to apply to each criteria. Internally normalized to
#' sum to 1.
#' @param score_min Minimum of total value score. Default is 0.
#' @param score_max Maximum of total value score. Default is 100.
#' @examples
#' n_samples <- 5
#' strategies <- c("Strategy 1", "Strategy 2")
#' outcome1 <- c(rnorm(n_samples, mean = 10, sd = 5),
#'               rnorm(n_samples, mean = 8, sd = 4))
#' outcome2 <- c(rnorm(n_samples, mean = 1500, sd = 90),
#'               rnorm(n_samples, mean = 1000, sd = 100))
#' outcomes <- data.frame(sample = rep(1:n_samples, length(strategies)),
#'                        strategy_id = rep(strategies, each = n_samples),
#'                        criteria1 = outcome1,
#'                        criteria2 = outcome2)
#'
#' # Performance matrix
#' performance_mat <- performance_matrix(outcomes, 
#'                                       strategy = "strategy_id", 
#'                                       criteria = c("criteria1", "criteria2"),
#'                                       rownames = c("Criteria 1", "Criteria 2"), 
#'                                       colnames = strategies)
#' print(performance_mat)
#'                                                                                                 
#' # MCDA                        
#' weights <- c(.7, .3)
#' mcda <- mcda(outcomes, sample = "sample", strategy = "strategy_id",
#'              criteria = c("criteria1", "criteria2"),
#'              weights = weights,
#'              optimal = c("low", "high"))
#' names(mcda)
#' 
#' # Scores on common scale
#' print(mcda$scores)
#' 
#' # "Total value"
#' print(mcda$total_value)
#' 
#' # "Total value" decomposed by criteria
#' print(mcda$weighted_scores)
#' 
#' # Probability of ranking
#' print(mcda$prob_rank)
#' @export
#' @seealso \code{\link{performance_matrix}}, \code{\link{lpvf_plot_data}}
mcda <- function(x, sample, strategy, criteria, 
                 criteria_min = NULL,
                 criteria_max = NULL,
                 optimal = NULL,
                 weights,
                 score_min = 0, score_max = 100){
  x <- data.table(x)
  n_samples <- length(unique(x[[sample]]))
  
  # Compute minimum and maximums for criteria
  if (is.null(criteria_min) | is.null(criteria_max)){
    if (is.null(optimal)){
      stop(paste0("If 'criteria_min' = NULL or 'criteria_max' = NULL ",
                  "then 'optimal' cannot be NULL."),
           call. = FALSE)
    }
    which_low <- which(optimal == "low")
    which_high <- which(optimal == "high") 
    criteria_low <- criteria[which_low]
    criteria_high <- criteria[which_high]
  }
  
  if (is.null(criteria_min)){
    criteria_min <- rep(NA, length(criteria))
    for (i in 1:length(criteria)){
      if (optimal[i] == "low"){
        criteria_min[i] <- max(x[[criteria[i]]])
      } else{
        criteria_min[i] <- min(x[[criteria[i]]])
      }
    }
  }
  if (is.null(criteria_max)){
   criteria_max <- rep(NA, length(criteria))
    for (i in 1:length(criteria)){
      if (optimal[i] == "low"){
        criteria_max[i] <- min(x[[criteria[i]]])
      } else{
        criteria_max[i] <- max(x[[criteria[i]]])
      }
    }
  }

  # Compute scores for each criteria on common scale
  scores <- matrix(NA, nrow = nrow(x), ncol = length(criteria))
  colnames(scores) <- criteria
  for (i in 1:length(criteria)){
    scores[, i] <- lpvf(x[[criteria[i]]], 
                                  x_min = criteria_min[i], x_max = criteria_max[i],
                                  score_min = score_min, score_max = score_max)
  }  
  scores_dt <- data.table(sample = x[[sample]],
                          strategy = x[[strategy]],
                          scores)
  
  # Compute weighted scores
  weights <- weights/sum(weights)
  weighted_scores <- t(t(scores) * weights)
  total_value <- rowSums(weighted_scores) 
  weighted_scores_dt = data.table(sample = x[[sample]],
                                  strategy = x[[strategy]],
                                  weighted_scores)
  weighted_scores_dt <- melt(weighted_scores_dt, 
                             id.vars = c("sample", "strategy"),
                             variable.name = "criteria",
                             value.name = "weighted_score")
  total_value_dt <- data.table(sample = x[[sample]],
                               strategy = x[[strategy]],
                               score = total_value)
  
  # Probability of ranking
  total_value_dt[, "rank" := frankv(-get("score"), ties.method = "random"), by = "sample"]
  prank_dt <- copy(total_value_dt)
  prank_dt <- setkey(prank_dt, strategy, rank)[CJ(strategy, rank, unique = TRUE), 
                                                .N, by = .EACHI]
  prank_dt[, "prob" := get("N")/n_samples]
  
  # Return
  setnames(scores_dt, c("sample", "strategy"), c(sample, strategy))
  setnames(weighted_scores_dt, c("sample", "strategy"), c(sample, strategy))
  setnames(total_value_dt, c("sample", "strategy"), c(sample, strategy))
  setnames(prank_dt, c("strategy"), c(strategy))
  res <- list(scores = scores_dt, 
              weighted_scores = weighted_scores_dt,
              total_value = total_value_dt,
              prob_rank = prank_dt)
}

#' Performance matrix
#' 
#' Create a performance matrix, which compares mean values of different
#' criteria across treatment strategies.
#' @param x A \code{data.frame} or \code{data.table} of simulation output 
#' characterizing the probability distribution of model outcomes. 
#' @param strategy Character name of column from \code{x} denoting treatment strategy.
#' @param criteria A vector of character names of columns from \code{x} denoting
#' the criteria to use in the MCDA.
#' @param cri If \code{TRUE}, credible intervals are computed; otherwise they are not.
#' @param digits Number of digits to use.
#' @param rownames Row names for returned table.
#' @param colnames Column names for returned table.
#' @return A matrix where each row is a criteria and a each column is a treatment
#' strategy.
#' @export
#' @seealso \code{\link{mcda}}
performance_matrix <- function(x, strategy, criteria, cri = TRUE, digits = 2,
                               rownames = NULL,
                               colnames = NULL){
  x <- data.table(x)
  
  # Format table based on number digtis to left of decimal place
  nchars <- sapply(c(x[1, criteria, with = FALSE]), function(x) nchar(trunc(x)))
  format_tbl <- function(x, nchars){
    x_str <- matrix(NA, nrow = nrow(x), ncol = ncol(x))
    for (i in 1:ncol(x)){
      if (nchars[i] <= 3){
        x_str[, i] <- formatC(x[, i], format = "f", digits = 2)
      } else{
        x_str[, i] <- formatC(x[, i], format = "d", big.mark = ",")
      }
    }  
    return(x_str)
  }
  
  # means
  means_tbl <- x[, lapply(.SD, mean), by = strategy, .SDcols = criteria] 
  means_mat <- as.matrix(means_tbl[, criteria, with = FALSE])
  means_mat_str <- format_tbl(means_mat, nchars)
  
  # credible intervals
  if (cri){
    lower <- x[, lapply(.SD, stats::quantile, .025), by = strategy, .SDcols = criteria] 
    lower_mat <- as.matrix(lower[, criteria, with = FALSE])
    lower_mat_str <- format_tbl(lower_mat, nchars)
    upper <- x[, lapply(.SD, stats::quantile, .975), by = strategy, .SDcols = criteria] 
    upper_mat <- as.matrix(upper[, criteria, with = FALSE])
    upper_mat_str <- format_tbl(upper_mat, nchars)
    tbl <- matrix(paste0(means_mat_str, " (", lower_mat_str, ", ", upper_mat_str, ")"),
                     nrow = nrow(means_mat_str), ncol = ncol(means_mat_str))
  } else{
    tbl <- means_mat_str
  }
  
  tbl <- t(tbl)
  if (is.null(rownames)){
    rownames(tbl) <- criteria
  } else{
    rownames(tbl) <- rownames
  }
  if (is.null(colnames)){
    colnames(tbl) <- means_tbl$strategy_id
  } else{
    colnames(tbl) <- colnames
  }  
  return(tbl)
}

#' Data to plot linear partial value function
#' 
#' Generate data used to plot a linear partial value function for a 
#' particular criteria given model outcomes.
#' @param x Model outcomes for a particular criteria on the original scale.
#' @param optimal If \code{"low"}, then lower values of the outcome are better, and,
#' if \code{"high"}, then higher values of the outcome are better.
#' @param length_out Number of points between minimum and maximum values of 
#' \code{x}. 
#' @return A \code{data.table} containing \code{x} and \code{y} coordinates for
#' a line plot.
#' @examples 
#' outcome <- rnorm(10, mean = 100, sd = 11)
#' plot_data <- lpvf_plot_data(outcome, optimal = "high")
#' print(plot_data)
#' @export
#' @seealso \code{\link{mcda}}
lpvf_plot_data <- function(x, optimal = c("low", "high"),
                           length_out = 1000){
  optimal <- match.arg(optimal)
  if (optimal == "low"){
    max_x <- min(x)
    min_x <- max(x)
  } else{
    min_x <- min(x)
    max_x <- max(x) 
  }
  if (min_x < max_x){
      x_data <-  seq(min_x, max_x, length.out = length_out)
      range_x <- c(min_x, max_x)
  } else {
      x_data <- seq(max_x, min_x, length.out = length_out)
      range_x <- c(max_x, min_x)
  }
  data <- data.table(x = x_data,
                    y = lpvf(x = x_data, x_min = min_x, x_max = max_x))
}

#' MCDA treatment attribute performance
#' 
#' Compute performance for MCDA criteria related to treatment attributes.
#' @param struct A \code{\link{model_structure}} object.
#' @param patients A data table returned from \code{\link{create_patients}}.
#' @param econmod An economic model of class \code{"IndivCtstm"}. Disease progression
#' must have been previously simulated (i.e., \code{$disprog_} cannot be \code{NULL}.)
#' @param treatments An object in the same format as \code{\link{treatments}}.
#' @return A \code{data.table} with four columns:
#' \describe{
#' \item{sample}{\code{sample} from \code{econmod$disprog_}}
#' \item{strategy_id}{\code{strategy_id} from \code{econmod$disprog_}}
#' \item{route}{Route of administration, weighted by the time each
#' simulated patient uses each therapy in a treatment sequence. Oral administration
#' is given a value of 1 and intravenous administration is given a value of 0.}
#' \item{yrs_since_approval}{Years since FDA approval, weighted by the time each
#' simulated patient uses each therapy in a treatment sequence.}
#' }
#' @seealso See the example in the \href{https://innovationvalueinitiative.github.io/IVI-NSCLC/articles/tutorial.html}{tutorial}.
#' @export
txattr_performance <- function(struct, patients, econmod, treatments = iviNSCLC::treatments){
  # Prevent CRAN warnings
  time <- time_stop <- time_start <- NULL
  route <- yrs_since_approval <- NULL
  weight <- weighted_route <- weighted_yrs_since_approval <- NULL
  
  # Errors
  check_is_class(struct, name = "struct", class = "model_structure") 
  check_is_class(patients, name = "patients", class = "patients")  
  check_is_class(econmod, name = "econmod", class = "IndivCtstm")   
  
  # Commpute performance
  txseq_dt <- tx_by_state(struct)
  txseq_dt <- merge(txseq_dt, 
                    treatments[, c("tx_name", "route", "yrs_since_approval"),
                               with = FALSE])
  txseq_dt[, route := ifelse(route == "IV", 0, 1)]
  
  disprog <- copy(econmod$disprog_)
  disprog[, time := time_stop - time_start]
  lys <- disprog[, lapply(.SD, sum),
                 .SDcols = "time",
                  by = c("sample", "strategy_id", "patient_id", "from", "to")]
  lys[, ("mutation") := patients[match(lys$patient_id, patients$patient_id)]$mutation]
  setnames(lys, "from", "state_id")
  lys <- merge(lys, 
               txseq_dt[, c("strategy_id", "state_id", "mutation", "route", 
                            "yrs_since_approval"), with = FALSE], 
               by = c("strategy_id", "state_id", "mutation"),
               sort = FALSE)
  lys[, weight := time/sum(time), by = c("sample", "strategy_id", "patient_id")]
  lys[, weight := ifelse(time == 0, 1, weight)]
  lys[, weighted_route := route * weight]
  lys[, weighted_yrs_since_approval := yrs_since_approval * weight]
  lys <- lys[, lapply(.SD, sum),
             .SDcols = c("weighted_route", "weighted_yrs_since_approval"),
             by = c("sample", "strategy_id", "patient_id")]
  lys <- lys[, lapply(.SD, mean),
             .SDcols = c("weighted_route", "weighted_yrs_since_approval"),
             by = c("sample", "strategy_id")]
  setnames(lys, c("weighted_route", "weighted_yrs_since_approval"), 
           c("route", "yrs_since_approval"))
  return(lys[,])
}