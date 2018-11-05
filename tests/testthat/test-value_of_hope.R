context("value_of_hope.R unit tests")
library("data.table")
library("flexsurv")
library("hesim")
rm(list = ls())


test_that("value_of_hope", {
  econmod <- example_IndivCtstm(n_samples = 2, n_patients = 3)

  # Works
  voh <- value_of_hope(econmod, comparator = 1)
  expect_true(inherits(voh, "data.table"))
  
  # Errors
  expect_error(value_of_hope(econmod = 2, comparator = 1))

})

utility_fun1 <- function(x, r){
    x^r
}

utility_fun2 <- function(x, r){
  if (r !=1){
    return((x^(1-r) - 1)/(1 - r))
  } else{
    return(log(x))
  }
}

r <- seq(0, 3, by = .5)
x <- seq(0, 10, by = .1)
u1 <- u2 <-  matrix(NA, nrow = length(x), ncol = length(r_vec))
colnames(u1) <- colnames(u2) <- r
for (i in 1:length(r_vec)){
  u1[, i] <- utility_fun1(x, r[i])
  u2[, i] <- utility_fun2(x, r[i])
}
u1_dt <- data.table(u1)
u1_dt$x <- x
u1_dt$utility_fun <- "Shafrin (2017)"
u2_dt <- data.table(u2)
u2_dt$x <- x
u2_dt$utility_fun <- "Standard"
u_dt <- rbind(u1_dt, u2_dt)
u_dt <- melt(u_dt, id.vars = c("x", "utility_fun"),
             variable.name = "r")
library("ggplot2")
ggplot(u_dt, aes(x = x, y = value, col = utility_fun)) + geom_line() + facet_wrap(~r, scales = "free") 

ggplot(u_dt[utility_fun == "Standard"], aes(x = x, y = value, col = utility_fun)) + geom_line() + facet_wrap(~r, scales = "free") 


