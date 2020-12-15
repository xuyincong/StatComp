#' @title rwMetropolis
#' @description Implement a random walk Metropolis sampler for generating the standard Laplace distribution
#' @param sigma variance of normal distrubtion
#' @param x0 Initial value
#' @param N length of chain
#' @return dataframe include number of candidate points rejected and chain
#' @examples
#' \dontrun{
#' N = 2000
#' sigma = c(.05, .5, 2, 16)
#' x0 = 25
#' rw1 = rw.Metropolis(sigma[1],x0,N)
#' }
#' @importFrom stats runif
#' @importFrom stats rnorm
#' @export
rw.Metropolis = function(sigma, x0, N){
  x = numeric(N)
  x[1] = x0
  u = runif(N)
  k = 0
  for (i in 2:N) {
    y = rnorm(1, x[i-1], sigma)
    if (u[i] <= (lap_f(y) / lap_f(x[i-1]))) x[i] = y 
    else {
      x[i] = x[i-1]
      k = k+1
    }
  }
  return(list(x = x, k = k))
}
lap_f = function(x) exp(-abs(x))


#' @title GelmanRubin
#' @description Use the Gelman-Rubin method to monitor convergence of the chain
#' @param psi psi is the statistic 
#' @return G-R statistic
#' @importFrom stats var
#' @export
Gelman.Rubin <- function(psi) {
  # psi[i,j] is the statistic psi(X[i,1:j])
  # for chain in i-th row of X
  psi <- as.matrix(psi)
  n <- ncol(psi)
  k <- nrow(psi)
  psi.means <- rowMeans(psi) #row means
  B <- n * var(psi.means) #between variance est.
  psi.w <- apply(psi, 1, "var") #within variances
  W <- mean(psi.w) #within est.
  v.hat <- W*(n-1)/n + (B/n) #upper variance est.
  r.hat <- v.hat / W #G-R statistic
  return(r.hat)
}

#' @title get roots of 11.4 
#' @description for k = 4 : 25, 100, 500, 1000, where t(k)  is a Student t random variable with k degrees of freedom.
#' @param k degree of freedom
#' @return root
#' @importFrom stats uniroot
#' @importFrom stats pt
#' @export
solve = function(k){
  output = uniroot(function(a){S(a,k)-S(a,k-1)},lower=1,upper=2)
  output$root
}

S = function(a,k){
  ck = sqrt(a^2*k/(k+1-a^2))
  pt(ck,df=k,lower.tail=FALSE)
}

#' @title Mle for A-B-O blood type problem  
#' @description calculate the corresponding log-maximum likelihood values (for observed data)
#' @param x p is probability of A and q is probability of B
#' @param x1 p1 is probability of A and q1 is probability of B
#' @param n.A number of n.A
#' @param n.B number of n.B
#' @param nAB number of nAB
#' @param noo number of noo
#' @return mle
#' @export
eval_f0 = function(x,x1,n.A=444,n.B=132,noo=361,nAB=63) {
  
  r1 = 1-sum(x1)
  nAA = n.A*x1[1]^2/(x1[1]^2+2*x1[1]*r1)
  nBB = n.B*x1[2]^2/(x1[2]^2+2*x1[2]*r1)
  r = 1-sum(x)
  return(-2*nAA*log(x[1])-2*nBB*log(x[2])-2*noo*log(r)-
           (n.A-nAA)*log(2*x[1]*r)-(n.B-nBB)*log(2*x[2]*r)-nAB*log(2*x[1]*x[2]))
}

#' @title mylapply
#' @description Implement a combination of Map() and vapply() to create an lapply() variant that iterates in parallel over all of its inputs and stores its outputs in a vector (or a matrix).
#' @param X data
#' @param FUN function
#' @param FUN.VALUE FUN.VALUE in vapply
#' @param simplify turn out into array 
#' @return function
#' @export
mylapply <- function(X, FUN, FUN.VALUE, simplify = FALSE){
  out <- Map(function(x) vapply(x, FUN, FUN.VALUE), X)
  if(simplify == TRUE) return(simplify2array(out))
  unlist(out)
}

