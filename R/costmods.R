#' Create cost models
#' 
#' Create a cost model for four cost categories: (i) treatment costs (i.e., drug
#' acquisition and administration costs), (ii) inpatient medical costs, 
#' (iii) outpatient medical costs, and (iv) costs due to adverse events.
#' @param n The number of random observations of the parameters to draw.
#' @param struct A \code{\link{model_structure}} object.
#' @param patients A data table returned from \code{\link{create_patients}}.
#' @param ae_probs An "ae_probs" object as returned by \code{\link{ae_probs}}.
#' @param params_costs_tx Parameter estimates for treatment costs (i.e.,
#' acquisition and administration costs) in the same format as 
#' \code{\link{params_costs_tx}}.
#' @param params_costs_op Parameter estimates for outpatient medical costs
#' in the same format as \code{\link{params_costs_op}}.
#' @param params_costs_inpt Parameter estimates for inpatient medical costs
#' in the same format as \code{\link{params_costs_inpt}}.
#' @param params_costs_ae Parameter estimates for adverse event costs
#' in the same format as \code{\link{params_costs_ae}}.
#' @examples
#' # Treatment sequences
#' txseq1 <- txseq(first = "erlotinib",
#'                 second = c("osimertinib", "PBDC"),
#'                 second_plus = c("PBDC + bevacizumab", "PBDC + bevacizumab"))
#' txseq2 <- txseq(first = "gefitinib",
#'                 second = c("osimertinib", "PBDC"),
#'                 second_plus = c("PBDC + bevacizumab", "PBDC + bevacizumab"))
#' txseqs <- txseq_list(seq1 = txseq1, seq2 = txseq2)
#'
#' # Patient population
#' pats <- create_patients(n = 2)
#'
#' # Model structure
#' struct <- model_structure(txseqs, dist = "weibull")
#'
#' ## Cost models
#' n_samples <- 2
#' ae_probs <- ae_probs(n = n_samples, struct = struct)
#' costmods <- create_costmods(n = 2, struct = struct, patients = pats,
#'                             ae_probs = ae_probs)
#' @return A list of objects of class "StateVals" from the 
#' \href{https://innovationvalueinitiative.github.io/hesim/}{hesim} package.
#' @export
create_costmods <- function(n = 100, struct, patients,
                            ae_probs,
                            params_costs_tx = iviNSCLC::params_costs_tx,
                            params_costs_op = iviNSCLC::params_costs_op,
                            params_costs_inpt = iviNSCLC::params_costs_inpt,
                            params_costs_ae = iviNSCLC::params_costs_ae
                            ){
  costmods <- list()
  costmods_tx <- create_costmod_tx(n, struct, patients, params_costs_tx)
  costmods$tx_ac <- costmods_tx$tx_ac
  costmods$tx_admin <- costmods_tx$tx_admin
  costmods$op <- create_costmod_default(n, struct, patients,
                                        params_costs_op)
  costmods$inpt <- create_costmod_default(n, struct, patients,
                                          params_costs_inpt)
  costmods$ae <- create_costmod_ae(n, struct, patients,
                                   ae_probs,
                                   params_costs_ae)
  return(costmods)
}

create_costmod_default <- function(n = 100, 
                                    struct, patients, 
                                    params){
  
  strategies <- data.table(strategy_id = 1:length(struct$txseqs))
  hesim_dat <- hesim::hesim_data(strategies = strategies,
                                  patients = patients)
  
  states <- create_states(struct)[get("state_name") != "D"]
  tbl <- merge(states, params, by = "state_name")  
  tbl <- hesim::stateval_tbl(tbl, dist = "gamma", hesim_data = hesim_dat)
  
  mod <- hesim::create_StateVals(tbl, n = n)
  
  return(mod)  
}

#' Annualized treatment costs
#' 
#' Compute annualized drug acquisition and administration costs.
#' @param x An object in the same format as \code{\link{params_costs_tx}}.
#' @examples 
#' annualized_tx_costs(iviNSCLC::params_costs_tx)
#' @return An object of class "annualized_tx_costs", which is a
#'  \code{data.table} of annualized treatment costs. Contains the columns:
#' \describe{
#' \item{agent}{See description from  the \code{dosage} element from \code{\link{params_costs_tx}}.}
#' \item{treatment}{See description from the \code{dosage} element from \code{\link{params_costs_tx}}.}
#' \item{dosage}{See description from the \code{dosage} element from \code{\link{params_costs_tx}}.}
#' \item{dose}{See description from the \code{dosage} element from \code{\link{params_costs_tx}}.}
#' \item{unit}{See description from the \code{dosage} element from \code{\link{params_costs_tx}}.}
#' \item{units_per_day}{See description from the \code{dosage} element from \code{\link{params_costs_tx}}.}
#' \item{duration_days}{See description from the \code{dosage} element from \code{\link{params_costs_tx}}.}
#' \item{discount_lower}{See description from the \code{discount} element from \code{\link{params_costs_tx}}.}
#' \item{discount_upper}{See description from the \code{discount} element from \code{\link{params_costs_tx}}.}
#' \item{acquisition_costs}{Annualized drug acquisition costs.}
#' \item{administration_costs}{Annualized drug administration costs.}
#' }
#' 
#' @seealso \code{\link{params_costs_tx}}
#' @export
annualized_tx_costs <- function(x = iviNSCLC::params_costs_tx){
  # Acquisition costs
  ## First unit
  tbl <- merge(x$dosage, 
                    x$acquisition_costs[, c("agent_name", "strength", 
                                            "acquisition_cost")], 
                    by.x = c("agent_name", "strength1"),
                    by.y = c("agent_name", "strength"))
  setnames(tbl,"acquisition_cost", "acquisition_cost1")
  
  ## Second unit
  tbl <- merge(tbl, 
               x$acquisition_costs[, c("agent_name", "strength",
                                       "acquisition_cost")], 
               by.x = c("agent_name", "strength2"),
               by.y = c("agent_name", "strength"),
               all.x = TRUE)
  setnames(tbl,"acquisition_cost", "acquisition_cost2")
  
  ## Combine first and second units
  tbl[, ("acquisition_cost2") := ifelse(is.na(get("acquisition_cost2")), 
                                        0, 
                                        get("acquisition_cost2"))]
  tbl[, ("acquisition_unit_cost") := get("quantity1") * get("acquisition_cost1") +
                                     get("quantity2") * get("acquisition_cost2")]  
  
  # Administration costs
  tbl <- merge(tbl, 
               x$administration_costs[, c("agent_name", "administration_cost")], 
               by = "agent_name",
               all.x = TRUE)
  setnames(tbl, "administration_cost", "administration_unit_cost")
  
  # Annualized costs
  tbl[, ("acquisition_costs") := get("units_per_day") * get("acquisition_unit_cost") * 365.25]
  tbl[, ("administration_costs") := get("units_per_day") * get("administration_unit_cost") * 365.25]
  
  # Set NULL
  tbl[, c("strength1", "strength2", "quantity1", "quantity2", 
          "acquisition_cost1", "acquisition_cost2", "acquisition_unit_cost",
          "administration_unit_cost", "source") := NULL]
  setorderv(tbl, c( "agent_name"))
  setcolorder(tbl, c("agent_name"))
  setattr(tbl, "class", c("annualized_tx_costs", "data.table", "data.frame"))
  return (tbl[,])
}
  
create_costmod_tx <- function(n = 100, 
                               struct, patients,
                               params){
  
  # 1. Treatments by strategy, health state, and mutation status
  txseq_dt <- tx_by_state(struct)

  # 2. Compute annualized costs
  annualized_costs <- annualized_tx_costs(params)
  annualized_costs <- annualized_costs[, c("agent_name", "duration_days", 
                                           "acquisition_costs",
                                           "administration_costs"), 
                                       with = FALSE]
  
  ## Subset to relevant treatments
  tx_names <- unique(txseq_dt$tx_name)
  lookup <- melt(params$lookup[get("tx_name") %in% tx_names],
                 id.vars = c("tx_name"),
                 value.name = "agent_name")
  lookup <- lookup[!is.na(get("agent_name")), c("tx_name", "agent_name")]
  annualized_costs <- merge(lookup,
                            annualized_costs,
                            by = "agent_name")
  
  ## Compute distribution of annualized costs. Relevant for acquisition costs
  ## due to uncertainty in discount rates.
  annualized_costs <- merge(annualized_costs, params$discounts,
                            by = "agent_name")
  discount_dist <- stats::runif(n * nrow(annualized_costs),
                                annualized_costs$discount_lower,
                                annualized_costs$discount_upper)
  annualized_costs_dist <- data.table(sample = rep(1:n, each = nrow(annualized_costs)),
                                      agent_name = annualized_costs$agent_name,
                                      tx_name = annualized_costs$tx_name,
                                      duration_days = annualized_costs$duration_days,
                                      discount = discount_dist,
                                      acquisition_costs = annualized_costs$acquisition_costs,
                                      administration_costs = annualized_costs$administration_costs)  
  annualized_costs_dist[, ("acquisition_costs") := (1 - get("discount")) * get("acquisition_costs")]
  
  ## Aggregate across agents to treatment level
  annualized_costs_dist <- annualized_costs_dist[, list(acquisition_costs = sum(get("acquisition_costs")),
                                                        administration_costs = sum(get("administration_costs"))),
                                                  by = c("sample", "tx_name", "duration_days")]
  setorderv(annualized_costs_dist, c("sample", "tx_name", "duration_days"), order = c(1, 1, -1))
  annualized_costs_dist[, ("cum_acquisition_costs") := cumsum(get("acquisition_costs")),
                        by = c("sample", "tx_name")]
  annualized_costs_dist[, ("cum_administration_costs") := cumsum(get("administration_costs")),
                        by = c("sample", "tx_name")]  
  annualized_costs_dist[, c("acquisition_costs", "administration_costs") := NULL]
  setnames(annualized_costs_dist, 
           c("cum_acquisition_costs", "cum_administration_costs"),
           c("acquisition_costs", "administration_costs"))
  setorderv(annualized_costs_dist, c("sample", "tx_name", "duration_days"))  
  
  ## Create time intervals and rectangularize dataset
  ### Time Intervals
  annualized_costs_dist[, ("time_interval") := 1:.N, by = c("sample", "tx_name")]
  annualized_costs_dist[, ("time_start") := ifelse(get("time_interval") == 1, 0, shift(get("duration_days"))),
                        by = "sample"]
  max_time <- annualized_costs_dist[, list(max_time = max(get("duration_days"))),
                                    by = "tx_name"]
  annualized_costs_dist[, c("time_interval", "duration_days") := NULL] 
  
  ### Rectangularize by time interval
  annualized_costs_dist_rect <- expand.grid(sample = unique(annualized_costs_dist$sample),
                                           tx_name = unique(annualized_costs_dist$tx_name),
                                           time_start = unique(annualized_costs_dist$time_start)) 
  annualized_costs_dist_rect <- data.table(annualized_costs_dist_rect)
  annualized_costs_dist_rect <- merge(annualized_costs_dist_rect, 
                                      annualized_costs_dist[, c("sample", "tx_name", "time_start", 
                                                              "acquisition_costs", "administration_costs"), with = FALSE],
                                    by = c("sample", "tx_name", "time_start"),
                                    all.x = TRUE)
  annualized_costs_dist_rect <- merge(annualized_costs_dist_rect, 
                                     max_time,
                                    by = c("tx_name"),
                                    all.x = TRUE)  
  annualized_costs_dist_rect[, ("lag_acquisition_costs") := shift(get("acquisition_costs")),
                             by = c("sample","tx_name")]
  annualized_costs_dist_rect[, ("lag_administration_costs") := shift(get("administration_costs")),
                             by = c("sample","tx_name")]  
  annualized_costs_dist_rect[, ("acquisition_costs") := ifelse(is.na(get("acquisition_costs")),
                                                               ifelse(get("time_start") < get("max_time"), get("lag_acquisition_costs"), 0),
                                                               get("acquisition_costs"))]
  annualized_costs_dist_rect[, ("administration_costs") := ifelse(is.na(get("administration_costs")),
                                                                  ifelse(get("time_start") < get("max_time"), get("lag_administration_costs"), 0),
                                                                  get("administration_costs"))]
  annualized_costs_dist_rect[, c("lag_acquisition_costs", "lag_administration_costs") := NULL] 
  
  # 3. Add annualized costs to treatment table
  txseq_dt <- merge(txseq_dt, annualized_costs_dist_rect,
                    by = "tx_name",
                    all.x = TRUE,
                    allow.cartesian = TRUE)

  # 4. Create StateVals models
  patients2 <- copy(patients)
  patients2[, ("grp_id") := ifelse(get("mutation") == 0, 1, 2)]
  hesim_dat <- hesim::hesim_data(strategies = data.table(strategy_id = 1:length(struct$txseqs)),
                                patients = patients2)
  
  
  mods <- vector(mode = "list", length = 2)
  names(mods) <- c("tx_ac", "tx_admin")
  
  ## Acquisition costs
  setnames(txseq_dt, "acquisition_costs", "value")
  stateval_tbl <- hesim::stateval_tbl(txseq_dt, dist = "custom",
                                      hesim_data = hesim_dat)
  mods[[1]] <- hesim::create_StateVals(stateval_tbl, n = n)
  setnames(txseq_dt, "value", "acquisition_costs")
  
  ## Administration costs
  setnames(txseq_dt, "administration_costs", "value")
  stateval_tbl <- hesim::stateval_tbl(txseq_dt, dist = "custom",
                                      hesim_data = hesim_dat)
  mods[[2]] <- hesim::create_StateVals(stateval_tbl, n = n)  
  
  # Return
  return(mods)
}

create_costmod_ae <- function(n = 100, 
                              struct, patients, 
                              ae_probs,
                              params_costs_ae){
  
  
  # Probability distribution for adverse event costs
  costs_ae_dist <- matrix(stats::rnorm(n * nrow(params_costs_ae),
                                      params_costs_ae$mean, params_costs_ae$se),
                                      nrow = n, byrow = TRUE)
  colnames(costs_ae_dist) <- params_costs_ae$ae_abb
  
  # Compute costs weighted by adverse event probabilities
  indices <- match(params_costs_ae$ae_abb, names(ae_probs))
  if (any(is.na(indices))){
    stop(paste0("The adverse event abbreviations in 'params_costs_ae' do not match ",
                "the names of the adverse events in 'ae_probs'."))
  }
  expected_costs <- vector(mode = "list", length = length(ae_probs))
  names(expected_costs) <- names(ae_probs)
  for (i in 1:length(ae_probs)){
    name_i <- names(ae_probs)[i]
    expected_costs[[i]] <- ae_probs[[i]] * costs_ae_dist[, name_i]
  }
  expected_costs <- Reduce('+', expected_costs)
  
  # Create model
  strategy_id <- 1:length(struct$txseqs)
  hesim_dat <- hesim::hesim_data(strategies = data.table(strategy_id = strategy_id),
                                 patients = patients,
                                 states = create_states(struct))
  ## stateval_tbl
  tbl1 <- data.table(strategy_id = rep(strategy_id, each = n),
                     sample = rep(1:n, times = length(strategy_id)),
                     value = c(expected_costs),
                     time_start = 0)
  tbl2 <- data.table(strategy_id = tbl1$strategy_id,
                     sample = tbl1$sample,
                     value = 0,
                     time_start = 1/12)
  tbl <- rbind(tbl1, tbl2)
  tbl <- hesim::stateval_tbl(tbl, dist = "custom", hesim_data = hesim_dat)
  mod <- hesim::create_StateVals(tbl, n = n)
  return(mod)  
}