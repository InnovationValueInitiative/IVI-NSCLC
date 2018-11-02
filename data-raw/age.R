rm(list = ls())
library("data.table")
library("ggplot2")
theme_set(theme_bw())

# Age distribution table
# Source: https://seer.cancer.gov/archive/csr/1975_2010/results_merged/topic_med_age.pdf
age_dist <- fread("age-distribution.csv")
age_dist[, age_top := ifelse(is.infinite(age_top), 100, age_top)]
age_dist[, age_cat := paste0(age_bot, " - ", age_top)]
age_dist[, age_cat := ifelse(age_bot == 85, "85+", age_cat)]
age_dist[, age_cat := factor(age_cat)]
age_dist[, prop := percentage/100]
age_dist[, age_mid := (age_top - age_bot)/2 + age_bot]
age_mean <- weighted.mean(x = age_dist$age_mid, w = age_dist$prop)
print(paste0("Mean age: ", age_mean))
age_var <-  sum(age_dist$prop * (age_dist$age_mid - age_mean)^2)
age_sd <- sqrt(age_var)
attr(age_dist, "mean") <- age_mean
attr(age_dist, "sd") <- age_sd

# A barplot to examine the distribution
p <- ggplot(age_dist, aes(x = age_cat, y = prop)) +
  geom_bar(stat = "identity", fill = "lightblue", color = "blue") + xlab("Age category") +
  ylab("Proportion")
ggsave("figs/age-barplot.pdf", p, width = 5, height = 7)

# Can we replicate the distribution using a normal distribution?
sim_data <- vector(mode = "list", length = nrow(age_dist))
for (i in 1:length(sim_data)){
  sim_data[[i]] <- runif(1000 * age_dist[i]$prop,
                        age_dist[i]$age_bot,
                        age_dist[i]$age_top)
}
sim_data <- unlist(sim_data)
sim_data <- data.table(age = c(sim_data))
                   
p <- ggplot(sim_data, aes(x = age)) +
  geom_histogram(aes(y = ..density..), fill = "lightblue", color = "blue") +
  stat_function(data = age_dist,
                mapping = aes(x = age_mid),
                fun = dnorm,
                args = list(mean = age_mean, sd = age_sd),
                col = "red")
ggsave("figs/age-density.pdf", p, width = 5, height = 7)


# Save 
setcolorder(age_dist, c("age_cat", "age_bot", "age_top", "age_mid", "percentage", "prop"))
age_dist[, percentage := NULL]
save(age_dist, file = "../data/age_dist.rda", compress = "bzip2")