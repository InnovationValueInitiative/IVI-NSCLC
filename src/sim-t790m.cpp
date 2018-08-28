#include <Rcpp.h>
using namespace Rcpp;

// This file is to create a function to simulate a patient's T790m mutation status.

// [[Rcpp::export]]
double test_fun(double x){
  return 2 * x;
}