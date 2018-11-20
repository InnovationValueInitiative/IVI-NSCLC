## ---- out.width = "800px", echo = FALSE----------------------------------
knitr::include_graphics("treatment-strategies.png")

## ---- out.width = "800px", echo = FALSE----------------------------------
knitr::include_graphics("modstruct4_1L.png")

## ---- out.width = "800px", echo = FALSE----------------------------------
knitr::include_graphics("modstruct3_1L.png")

## ---- message=FALSE, warning=FALSE---------------------------------------
library("iviNSCLC")
library("hesim")
library("data.table")
library("ggplot2")
library("gridExtra")
library("scales")
library("knitr")
ggplot2::theme_set(theme_minimal())

## ------------------------------------------------------------------------
pats <- create_patients(n = 1000)

## ------------------------------------------------------------------------
txseq1 <- txseq(first = "erlotinib",
                second = c("osimertinib", "PBDC"),
                second_plus = c("PBDC + bevacizumab", "PBDC + bevacizumab"))
txseq2 <- txseq(first = "gefitinib",
                second = c("osimertinib", "PBDC"),
                second_plus = c("PBDC + bevacizumab", "PBDC + bevacizumab"))
txseqs <- txseq_list("Sequence 1" = txseq1, "Sequence 2" = txseq2, start_line = "first",
                     mutation = "unknown")

## ------------------------------------------------------------------------
# Convenience function to add factor names to data table
# for plotting
strategy_factor <- function(x, rev = FALSE){ 
  strategy_names <- names(txseqs)
  if (rev == FALSE){
    x[, strategy_name := factor(strategy_id, levels = 1:length(txseqs), 
                               labels = strategy_names)]
  } else{
    x[, strategy_name := factor(strategy_id, levels = rev(1:length(txseqs)), 
                                 labels = rev(strategy_names))]
  }
}

## ------------------------------------------------------------------------
txseqs_2L_pos <- txseq_list("Sequence 1" = txseq1, "Sequence 2" = txseq2, start_line = "second",
                              mutation = "positive")
txseqs_2L_neg <- txseq_list("Sequence 1" = txseq1, "Sequence 2" = txseq2, start_line = "second",
                              mutation = "negative")

## ------------------------------------------------------------------------
struct <- model_structure(txseqs, dist = "weibull")

# Health states
states <- create_states(struct)
print(states)

# Transition matrix
tmat <- create_trans_mat(struct)
print(tmat)

## ------------------------------------------------------------------------
# 1st line model
struct_3_1L <- model_structure(txseqs, dist = "weibull",
                            n_states = "three")
print(create_states(struct_3_1L))
print(create_trans_mat(struct_3_1L))

# 2nd line model (T790M positive)
struct_3_2L <- model_structure(txseqs_2L_pos, dist = "weibull",
                              n_states = "three")
print(create_states(struct_3_2L))
print(create_trans_mat(struct_3_2L))

## ------------------------------------------------------------------------
tmat_3_1L <- create_trans_mat(struct_3_1L)
tmat_3_2L <- create_trans_mat(struct_3_2L)
print(tmat_3_1L)
print(tmat_3_2L)

## ------------------------------------------------------------------------
state_factor <- function(x){
  x[, state_name := factor(state_id, levels = 1:nrow(states), 
                           labels = states$state_name)]
}

## ------------------------------------------------------------------------
n_samples <- 100

## ------------------------------------------------------------------------
# Input data
transmod_data <- create_transmod_data(struct, tmat, pats)
## Print first 5 rows and and 10 covariates from data
print(transmod_data[1:5, 1:10]) 

# Parameters
transmod_params <- create_transmod_params(n = n_samples, data = transmod_data)
## Print first 5 samples from the probability distribution and 4 covariates (which
## match those in 'transmod_data')
transmod_params$coefs$scale[1:5, 1:4]

# Health state transition model
transmod <- hesim::create_IndivCtstmTrans(transmod_params, transmod_data, tmat,
                                          start_age = pats$age)
class(transmod)

## ----ae_probs, fig.height = 4--------------------------------------------
ae_probs <- ae_probs(n_samples, struct)

## ------------------------------------------------------------------------
utilmod <- create_utilmod(n = n_samples, struct = struct, patients = pats,
                          ae_probs = ae_probs)

## ------------------------------------------------------------------------
costmods <- create_costmods(n = n_samples, struct = struct, patients = pats, ae_probs = ae_probs)

## ------------------------------------------------------------------------
econmod <- hesim::IndivCtstm$new(trans_model = transmod,
                                 utility_model = utilmod,
                                 cost_models = costmods)

## ------------------------------------------------------------------------
econmod$sim_disease(max_t = 20)
econmod$sim_stateprobs(t = seq(0, 20 , 1/26)) # Biweekly probabilities

## ----stateprobs----------------------------------------------------------
# Plot of state probabilities
stprobs <- econmod$stateprobs_[, .(prob_mean = mean(prob),
                                     prob_lower = quantile(prob, .025),
                                     prob_upper = quantile(prob, .975)),
                                 by = c("strategy_id", "state_id", "t")]
state_factor(stprobs)
strategy_factor(stprobs)
ggplot(stprobs[t < 2], aes(x = t, y = prob_mean, col = strategy_name)) +
  geom_line() + facet_wrap(~state_name) + 
  xlab("Years") + ylab("Probability in health state") +
  scale_color_discrete(name = "") 

## ----pfs_os_curves-------------------------------------------------------
curves <- econmod$stateprobs_[state_id == 1 | state_id == 4]
curves[, prob := ifelse(state_id == 4, 1 - prob, prob)] # OS is 1 - pr(death)
curves[, lab := factor(state_id, levels = c(1, 4), labels = c("PFS", "OS"))]
curves <- curves[, .(prob_mean = mean(prob),
                            prob_lower = quantile(prob, .025),
                            prob_upper = quantile(prob, .975)),
                          by = c("strategy_id", "lab", "t")]
strategy_factor(curves)
ggplot(curves[t < 2], aes(x = t, y = prob_mean, col = strategy_name)) +
  geom_line() + facet_wrap(~lab) + 
  xlab("Years") + ylab("Survival probability") +
  scale_color_discrete(name = "") 

## ----pfs_os_quantiles, fig.height = 3------------------------------------
quantiles <- hesim::surv_quantile(curves, probs = c(.025, .25, .5, .75, .975),
                                  t = "t",
                                  surv_cols = c("prob_mean", "prob_lower", "prob_upper"),
                                  by = c("strategy_id", "lab"))
strategy_factor(quantiles, rev = TRUE)
ggplot(quantiles[prob == .5], aes(x = strategy_name, y = quantile_prob_mean)) + 
  facet_wrap(~lab, nrow = 2) +
  geom_bar(stat = "identity", fill = "#d9230f") +
  geom_errorbar(aes(ymin = quantile_prob_lower,
                    ymax = quantile_prob_upper), width = .2) +
  scale_fill_discrete(name = "") +
  xlab("") + ylab("Median years") + 
  coord_flip() 

## ----ae_probs_plot, fig.height = 6.5-------------------------------------
tidy_ae_probs <- tidy(ae_probs)
plot_data <- tidy_ae_probs[, .(prob_mean = mean(prob),
                               prob_lower = quantile(prob, .025),
                                prob_upper = quantile(prob, .975)),
                          by = c("strategy_id", "ae_name")]
strategy_factor(plot_data, rev = TRUE)
plot_data[, ae_name := ifelse(ae_name == "Elevated alanine transaminase",
                              "Elevated \n alanine \n transaminase",
                              ae_name)]
plot_data[, ae_name := ifelse(ae_name == "Elevated aspartate transaminase",
                              "Elevated \n aspartate \n transaminase",
                              ae_name)]
ggplot(plot_data, aes(x = rev(strategy_name), y = prob_mean)) + 
  facet_wrap(~factor(ae_name), ncol = 3) +
  geom_bar(stat = "identity", fill = "#d9230f") +
  geom_errorbar(aes(ymin = prob_lower,
                    ymax = prob_upper), width = .2) +
  scale_fill_discrete(name = "") +
  scale_y_continuous(breaks = pretty_breaks(n = 3)) +
  xlab("") + ylab("Probability of adverse event") + 
  coord_flip() +
  theme(panel.spacing = unit(2, "lines"))

## ----mean_lys, fig.height = 1.5------------------------------------------
# Simulate
econmod$sim_qalys(dr = c(0, .03))

# Plot
lys <- econmod$qalys_[dr == 0, .(mean_lys = mean(lys)),
                            by = c("strategy_id", "state_id")]
state_factor(lys)
strategy_factor(lys, rev = TRUE)
ggplot(lys, aes(x = strategy_name, y = mean_lys, fill = state_name)) + 
  geom_bar(stat = "identity") +
  scale_fill_discrete(name = "") +
  xlab("") + ylab("Mean life-years") + 
  coord_flip() 

## ----costs, fig.width = 4------------------------------------------------
# Simulate
econmod$sim_costs(dr = .03)
ce_sim <- econmod$summarize()

# Plot
costs <- ce_sim$costs[dr == .03 , .(costs_mean = mean(costs),
                           costs_lower = quantile(costs, .025),
                           costs_upper = quantile(costs, .975)),
                        by = c("strategy_id", "category")]
strategy_factor(costs)
ggplot(costs[category == "total"], 
       aes(x = strategy_name, y = costs_mean)) + 
  geom_bar(stat = "identity", fill = "#d9230f") +
  geom_errorbar(aes(ymin = costs_lower, 
                    ymax = costs_upper), width=.2) +
  scale_fill_discrete(name = "") +
  xlab("") + ylab("Costs") + 
  scale_y_continuous(label = dollar_format()) 

## ----costs_cat, fig.width = 5--------------------------------------------
costs[, category_name := factor(category, 
                                levels = c("ae", "tx_ac", "tx_admin", "inpt", "op", "total"),
                                labels = c("Adverse event", "Drug acquisition", "Drug administration",
                                           "Inpatient", "Outptient", "Total"))]
ggplot(costs[category != "total"], 
       aes(x = strategy_name, y = costs_mean, fill = category_name)) + 
  geom_bar(stat = "identity") +
  scale_fill_discrete(name = "") +
  xlab("") + ylab("Costs") + 
  scale_y_continuous(label = dollar_format()) 

## ----incr_costs, fig.width = 4-------------------------------------------
# Compute incremental costs for each possible reference treatment sequence
incr_costs <- vector(mode = "list", length = length(txseqs))
incr_costs_i <- copy(ce_sim$costs)
for (i in 1:length(txseqs)){
  incr_costs_i[, costs_comparator := costs[i], 
                           by = c("category", "dr", "sample")]
  incr_costs_i[, icosts := costs - costs_comparator]
  incr_costs[[i]] <- incr_costs_i[, .(icosts_mean = mean(icosts),
                                   icosts_lower = quantile(icosts, .025),
                                   icosts_upper = quantile(icosts, .975)),
                               by = c("category", "dr", "strategy_id")]
  strategy_factor(incr_costs[[i]])
}

# Plot incremental costs with treatment sequence 1 as the reference treatment
ggplot(incr_costs[[1]][dr == .03 & category == "total" & 
                    strategy_id != 1], 
       aes(x = strategy_name, y = icosts_mean)) + 
  geom_bar(stat = "identity", fill = "#d9230f") +
  geom_errorbar(aes(ymin = icosts_lower, 
                    ymax = icosts_upper), width=.2) +
  scale_fill_discrete(name = "") +
  scale_x_discrete(drop = FALSE) +
  xlab("") + ylab("Incremental costs") + 
  scale_y_continuous(label = dollar_format())

## ----compute_prod_costs--------------------------------------------------
prodcosts <- sim_prod_costs(econmod, patients = pats)
ce_sim2 <- add_prod_costs(ce_sim, prod_costs = prodcosts)

## ----plot_prod_costs_cat, fig.width = 5----------------------------------
costs2 <- ce_sim2$costs[dr == .03 , .(costs_mean = mean(costs)),
                        by = c("strategy_id", "category")]
strategy_factor(costs2)
costs2[, category_name := factor(category, 
                                levels = c("ae", "tx_ac", "tx_admin", "prod", "inpt", "op", "total"),
                                labels = c("Adverse event", "Drug acquisition", "Drug administration",
                                           "Productivity", "Inpatient", "Outptient", "Total"))]
ggplot(costs2[category != "total"], 
       aes(x = strategy_name, y = costs_mean, fill = category_name)) + 
  geom_bar(stat = "identity") +
  scale_fill_discrete(name = "") +
  xlab("") + ylab("Costs") + 
  scale_y_continuous(label = dollar_format()) 

## ----icea----------------------------------------------------------------
icea <- hesim::icea(ce_sim, dr = .03)
icea_pw <- hesim::icea_pw(ce_sim, comparator = 1, dr = .03)

## ----icer----------------------------------------------------------------
icer <- hesim::icer_tbl(icea_pw,
                        k = 100000, # WTP per QALY 
                        cri = TRUE,
                        rownames = c("Incremental QALYs", "Incremental costs ($)", 
                                     "Incremental NMB ($)", "ICER ($ per QALY)", 
                                     "Conclusion"),
                        colnames = names(txseqs))
knitr::kable(icer)

## ----ceplane-------------------------------------------------------------
ylim <- max(icea_pw$delta[, ic]) * 1.1
xlim <- max(icea_pw$delta[, ie]) * 1.1
strategy_factor(icea_pw$delta)
ggplot(icea_pw$delta, aes(x = ie, y = ic, col = strategy_name)) + 
  geom_jitter(size = .5)  + 
  xlab("Incremental QALYs") + ylab("Incremental costs") +
  scale_y_continuous(label = dollar, limits = c(-ylim, ylim)) +
  scale_x_continuous(limits = c(-xlim, xlim)) +
 scale_colour_discrete(name = "") +
  geom_abline(slope = 100000, linetype = "dashed") +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0)

## ----ceac----------------------------------------------------------------
strategy_factor(icea_pw$ceac)
ggplot(icea_pw$ceac, aes(x = k, y = prob, col = strategy_name)) +
  geom_line() + xlab("Willingness to pay") +
  ylab("Probability cost-effective") +
  scale_x_continuous(label = scales::dollar) +
  scale_colour_discrete(name = "")

## ----evpi----------------------------------------------------------------
ggplot(icea$evpi, aes(x = k, y = evpi)) +
  geom_line() + xlab("Willingness to pay") +
  ylab("Expected value of perfect information") +
  scale_x_continuous(label = scales::dollar) +
  scale_y_continuous(label = scales::dollar)

## ----performance-matrix1-------------------------------------------------
hc_costs <- ce_sim$costs[!category %in% c("prod", "total"), 
                           .(costs = sum(costs)),
                             by = c("sample", "strategy_id")]
txattr <- txattr_performance(struct = struct, patients = pats, econmod = econmod)
outcomes <- data.table(sample = rep(1:n_samples, each = length(txseqs)),
                       strategy_id = rep(1:length(txseqs), times = n_samples),
                       lys_s1 = econmod$qalys_[state_id == 1 & dr == 0.03, lys],
                       lys_p1 = econmod$qalys_[state_id == 2 & dr == 0.03, lys],
                       lys_p2 = econmod$qalys_[state_id == 3 & dr == 0.03, lys],
                       hc_costs = hc_costs$costs,
                       route = txattr$route,
                       yrs_since_approval = txattr$yrs_since_approval)
criteria_vars <- c("lys_s1", "lys_p1", "lys_p2", "hc_costs", "route", "yrs_since_approval")
criteria_names <- c("Life-years stable disease", "Life-years progressed 1L", 
                    "Life-years progressed 2L", "Health care sector costs",
                    "Route of administration", "Years since FDA approval")
optimal <- c("high", "high", "high", "low", "high", "high")
performance_mat_orig <- performance_matrix(outcomes, 
                                           strategy = "strategy_id", 
                                           criteria = criteria_vars,
                                           rownames = criteria_names, 
                                           colnames = names(txseqs))
knitr::kable(performance_mat_orig)

## ----lpvf----------------------------------------------------------------
plot_data <- lpvf_plot_data(outcomes$hc_costs, optimal = "low")
ggplot(plot_data, aes(x = x, y = y)) + geom_line() +
  xlab("Health care sector costs") + ylab("Score") + 
  scale_x_continuous(label = dollar_format()) 

## ----mcda----------------------------------------------------------------
weights <- c(.2, .2, .2, .3, .05, .05)
mcda <- mcda(outcomes, sample = "sample", strategy = "strategy_id",
             criteria = criteria_vars,
             weights = weights,
             optimal = optimal)

## ----performance-matrix--------------------------------------------------
performance_mat_common <- performance_matrix(mcda$scores, 
                                             strategy = "strategy_id", 
                                             criteria = criteria_vars,
                                             rownames = criteria_names, 
                                             colnames = names(txseqs))
knitr::kable(performance_mat_common)

## ----mcda-total-value, fig.width = 4-------------------------------------
# Total value
total_value <- mcda$total_value[, .(mean = mean(score),
                                    lower = quantile(score, .025),
                                    upper = quantile(score, .975)),
                                by = c("strategy_id")]
strategy_factor(total_value)
ggplot(total_value, 
       aes(x = strategy_name, y = mean)) + 
  geom_bar(stat = "identity", fill = "#d9230f") +
  geom_errorbar(aes(ymin = lower, 
                    ymax = upper), width=.2) +
  scale_fill_discrete(name = "") +
  xlab("") + ylab("Total value") 

## ----mcda-total-value-decompose, fig.width = 5---------------------------
# Weighted scores
weighted_scores <- mcda$weighted_scores[, .(mean = mean(weighted_score)),
                                        by = c("strategy_id", "criteria")]
weighted_scores[, criteria := factor(criteria, 
                                     levels = criteria_vars,
                                labels = criteria_names)]
strategy_factor(weighted_scores)
ggplot(weighted_scores, 
       aes(x = strategy_name, y = mean, fill = criteria)) + 
  geom_bar(stat = "identity") +
  scale_fill_discrete(name = "") +
  xlab("") + ylab("Total value") 

## ----mcda-prob-rank------------------------------------------------------
prank <- mcda$prob_rank
strategy_factor(prank)
prank_levels <- sort(unique(prank$rank))
prank[, f1_rank := factor(rank)]
prank[, f2_rank := factor(rank,
                          levels = prank_levels,
                          labels = paste0("Rank = ", prank_levels))]

# Rank on x-axis
p1 <- ggplot(prank, 
       aes(x = f1_rank, y = prob, fill = strategy_name)) + 
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_discrete(name = "") +
  xlab("Rank") + ylab("Probability") +
  theme(legend.position = "bottom")

# Treatment sequence on x-axis
p2 <- ggplot(prank, 
       aes(x = strategy_name, y = prob, fill = f2_rank)) + 
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_discrete(name = "") +
  xlab("") + ylab("Probability") +
  theme(legend.position = "bottom")

# Combine
grid.arrange(p1, p2, ncol = 2)

## ----voh-----------------------------------------------------------------
voh <- value_of_hope(econmod, comparator = 1)

## ----voh-plot-qalys------------------------------------------------------
plot_data <- copy(voh)
plot_data[, iqalys_voh := iqalys + voh]
plot_data <- melt(plot_data, id.vars = c("strategy_id"),
                  measure.vars = c("iqalys", "iqalys_voh"),
                  value.name = "qalys")
plot_data[, variable := ifelse(variable == "iqalys", 
                               "Standard",
                               "With value of hope")]
strategy_factor(plot_data)
ggplot(plot_data, 
       aes(x = strategy_name, y = qalys, fill = variable)) + 
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_discrete(name = "") +
  geom_hline(yintercept = 0) +
  xlab("") + ylab("Incremental QALYs")

## ----voh-plot-nmb--------------------------------------------------------
k <- 100000
incr_costs_mean <- incr_costs[[1]][category == "total" & dr == .03, 
                                   .(strategy_id, icosts_mean)]
plot_data[, icosts := incr_costs_mean$icosts_mean[match(strategy_id,
                                            incr_costs_mean$strategy_id)]]
plot_data[, inmb := k * qalys - icosts]
ggplot(plot_data, 
       aes(x = strategy_name, y = inmb, fill = variable)) + 
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_discrete(name = "") +
  geom_hline(yintercept = 0) +
  xlab("") + ylab("Incremental NMB") + 
  scale_y_continuous(label = scales::dollar)

