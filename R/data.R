#' Treatments
#'
#' Dataset of treatments for EGFR positive NSCLC.
#'
#' @format A data.table object with 1 row for each available treatment and with
#' columns:
#' \describe{
#'   \item{tx_name}{Name of treatment.}
#'   \item{tx_abb}{Abbreviation for the treatment.}
#'   \item{route}{Route of administration.}
#'   \item{approval_data}{Date that treatment was approved in the US by the Federal Drug Administration (FDA)}
#'   \item{years_since_approval}{Years (from last update of the model) since the treatment was approved.}
#' }
#' @examples 
#' print(treatments)
"treatments"

#' Multi-state NMA parameters
#'
#' Posterior distributions of the regression coefficients of the continuous 
#' time state transition model estimated using the multi-state network meta-analysis. 
#'
#' @format A list of \code{\link[hesim]{params_surv}} objects from the 
#' \href{https://innovationvalueinitiative.github.io/hesim/}{hesim} package
#' package. The list contains regression coefficient estimates from the Weibull,
#' Gompertz, and 2nd order fractional polynomial survival distributions. The 
#' parameter (e.g., scale and shape for the Weibull distribution) of each survival 
#' distribution are predicted as a function of treatments, health state transitions, and
#' treatment history (for 2L treatments).
"params_mstate_nma"

#' Adverse event NMA parameters
#'
#' The posterior distribution of the probability of an adverse event by available
#' first and second line treatments. Based on separate models by adverse event 
#' with a binomial likelihood and a logit link.
#'
#' @format A list where each element is a matrix for a distinct adverse event.
#' Each column is a distinct treatment and each row is a random draw from
#' the posterior distribution. Each matrix contains a "imputed" attribute, which
#' is vector of 0's and 1's with length equal to the number of columns, and a 1
#' denotes whether values of a column have been imputed because of insufficient
#' clinical trial data.
#' 
#'  When there is missing data for a TKI for a given adverse 
#' event, it is assumed to have the same adverse event probabilities as the TKI 
#' (among the TKI's with data) with the highest mean adverse event probability. When there
#' is missing data for a combination therapy containing PBDC, it is assumed to have
#' the same adverse event probabilities as PBDC. Finally, when there is not data for
#' any PBDC based therapy, all PBDC based therapies are assumed to have the same
#' adverse event probabilities as the TKI with the highest mean adverse event probability.
#' @examples 
#' names(params_ae_nma)
#' head(params_ae_nma[[1]])
#' attr(params_ae_nma[[1]], "imputed")
"params_ae_nma"

#' Utility parameters
#'
#' Utility estimates by health state and disutilities by adverse event.
#'
#' @format A list containing the following elements:
#' \itemize{
#' \item{state_utility}{ Utility estimates by health state.}
#' \item{as_disutility}{ Disutility by health state.}
#' }
#' 
#' @section State utility:
#' The \code{state_utility} element is a data table with the following columns:
#' \describe{
#' \item{state_name}{Name of the health state.}
#' \item{mean}{Mean utility.}
#' \item{se}{Standard error of utility.}
#' \item{ref}{BibTeX reference for estimate.}
#' }
#' 
#' @section Adverse event disutilities:
#' The \code{as_disutility} element is a data table with the following columns:
#' \describe{
#' \item{ae_name}{Name of the adverse event.}
#' \item{ae_abb}{Abbreviation for the adverse event.}
#' \item{mean}{Mean disutility.}
#' \item{se}{Standard error of disutility.}
#' \item{ref}{BibTeX reference for estimate.}
#' }
#' @examples 
#' print(params_utility)
"params_utility"

#' Treatment cost parameters
#'
#' Information on dosage, strength, pricing, and administration costs needed
#' to estimate costs of treatments sequences.
#'
#' @format A list with the following elements:
#' 
#' @section dosage:
#' A \code{data.table} containing information on dosage for each
#' agent used for treatments that can be included in a treatment sequence. The variables
#' \code{strength1}, \code{strength2}, \code{quantity1}, \code{quantity2}, and \code{units_per_day} can
#' be modified to change annualized costs. 
#' \describe{
#' \item{agent_name}{Name of the agent used in the treatment; combination therapies will have
#' more than one agents.}
#' \item{dosage}{Amount, number. and frequency of doses over a specified time period.}
#' \item{dose}{Specified amount of medication taken at one time. By default, doses
#' based on body surface area (BSA) and weight are based on mean BSA and weight,
#' respectively.}
#' \item{unit}{The unit used to compute costs.}
#' \item{strength1}{The amount of the drug in a single dosage form (e.g., tablet or vial).}
#' \item{strength2}{If a unit consists of multiple dosage forms, then the amount of the drug
#' in the dosage form not included in \code{strength1}.}
#' \item{quantity1}{Quantity of dosage form specified in \code{strength1} in a unit.}
#' \item{quantity2}{Quantity of dosage form specified in \code{strength2} in a unit.}
#' \item{units_per_day}{Number of units of the treatment (based on dosage) used per day.}
#' \item{duration_days}{Number of days to use the treatment. If used until progression, then 
#' must equal \code{Inf}.}
#' \item{source}{Source for dosage information.}
#' }
#' 
#' @section acquisition_costs:
#' A \code{data.table} of acquisition costs by \code{agent} and \code{strength}.
#' Default costs are based on wholesale acquisition costs (WACs). The variable 
#' \code{acquisition_cost} can be modified to change annualized acquisition costs.
#' \describe{
#' \item{agent_name}{Name of the agent}
#' \item{strength}{The amount of the agent in a given unit of the dosage form.}
#' \item{acquisition_cost}{The acquisition cost of an \code{agent} for a given \code{strength}.}
#' \item{source}{Source used for \code{acquisition cost}.}
#' }
#' 
#' @section discounts:
#' \describe{
#'  \item{agent_name}{Name of the agent.}
#'  \item{discount_lower}{Lower bound for discount and rebates as a 
#'  proportion of \code{acquisition_cost} in the \code{acquisition_costs} table.
#'   Assumed to follow a uniform distribution in the probabilistic sensitivity analysis.}
#'  \item{discount_upper}{Upper bound for discount and rebates as a 
#'  proportion of \code{acquisition_cost} in the \code{acquisition_costs} table.
#'   Assumed to follow a uniform distribution in the probabilistic sensitivity analysis.}
#' }
#' 
#' @section administration_costs:
#' A \code{data.table} of administration costs. The variable \code{costs}
#' can be modified to change annualized administration costs.
#' \describe{
#' \item{agent_name}{Name of the agent.}
#' \item{administration_cost}{The costs of administering an agent. }
#' \item{source}{Source used for \code{administration cost}.}
#' }
#' 
#' @section lookup:
#' A lookup \code{data.table} used to match treatments from the \code{\link{treatments}}
#' table to agents comprising a particular treatment.
#' \describe{
#' \item{tx_name}{Name of treatment.}
#' \item{agent_name1}{Name of first agent used for treatment.}
#' \item{agent_name2}{Name of second agent used for treatment.}
#' \item{agent_name3}{Name of third agent used for treatment.}
#' \item{agent_name4}{Name of fourth agent used for treatment.}
#' } 
#' 
#' 
#' @examples 
#' print(params_costs_tx)
"params_costs_tx"

#' Outpatient cost parameters
#'
#' Mean outpatient costs by health state.
#'
#' @format A data table with the following columns:
#' \describe{
#' \item{state_name}{Name of the health state.}
#' \item{mean}{Mean outpatient costs.}
#' \item{se}{Standard error of outpatient costs.}
#' \item{ref}{BibTeX reference for estimate.}
#' }
#' @examples 
#' print(params_costs_op)
"params_costs_op"

#' Inpatient cost parameters
#'
#' Mean inpatient costs by health state.
#'
#' @format A data table with the following columns:
#' \describe{
#' \item{state_name}{Name of the health state.}
#' \item{mean}{Mean outpatient costs.}
#' \item{se}{Standard error of outpatient costs.}
#' \item{ref}{BibTeX reference for estimate.}
#' }
#' @examples 
#' print(params_costs_inpt)
"params_costs_inpt"

#' Adverse event cost parameters
#'
#' Mean costs by adverse event.
#'
#' @format A data table with the following columns:
#' \describe{
#' \item{ae_name}{Name of the adverse event.}
#' \item{ae_abb}{Abbreviation for the adverse event.}
#' \item{mean}{Mean costs.}
#' \item{lower}{2.5\% quantile for mean costs. Currently assumed to be +20\%
#' of mean costs given lack of data for estimating standard errors.}
#' \item{upper}{97.5\% quantile for mean costs. Currently assumed to be -20\%
#' of mean costs given lack of data for estimating standard errors.}
#' \item{se}{Standard error of costs.}
#' \item{ref}{Reference for mean costs.}
#' }
#' @examples 
#' print(params_costs_ae)
"params_costs_ae"