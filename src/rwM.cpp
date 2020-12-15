#include <Rcpp.h>
using namespace Rcpp;

//' @title rwMetropolis
//' @description Implement a random walk Metropolis sampler for generating the standard Laplace distribution
//' @param sigma variance of normal distrubtion
//' @param x0 Initial value
//' @param N length of chain
//' @return dataframe include number of candidate points rejected and chain
//' @examples
//' \dontrun{
//' N = 2000
//' sigma = c(.05, .5, 2, 16)
//' x0 = 25
//' rw1 = rw.Metropolis(sigma[1],x0,N)
//' }
//' @export
// [[Rcpp::export]]
NumericVector rwMetropolis (double sigma, double x0, int N) {
  NumericVector x(N);
  x[0] = x0; 
  NumericVector u = runif(N);
  for (int i = 1; i < N;i++ ) {
    NumericVector y = rnorm(1, x[i-1], sigma);
    if (u[i] <= (exp(-abs(y[0])) / exp(-abs(x[i-1])))){
      x[i] = y[0];
    }
    else { 
      x[i] = x[i-1]; 
    }
  }
  return(x);
} 

