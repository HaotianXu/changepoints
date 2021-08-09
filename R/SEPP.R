#' @title Simulate a (stable) SEPP model (without change point)
#' @description TO DO
#' @param intercept  A \code{numeric} scalar representing the intercept of the model.
#' @param M          A \code{integer} scalar representing dimensionality.
#' @param A          A \code{numeric} M-by_M matrix representing the coefficient matrix.
#' @param threshold  A \code{numeric} scalar representing the upper bound for each coordinate of X_t (for stability).
#' @param n          A \code{integer} scalar representing sample size.
#' @param vzero      A \code{numeric} vector representing the observation at time 0.
#' @param ...        Additional arguments.
#' @return  A p-by-n matrix.
#' @export
#' @author 
#' @examples
#' p = 20
#' sigma = 1
#' n = 100
#' A = matrix(rnorm(p*p), nrow = p)
#' simu.change.VAR1(sigma, p, n, A)
simu.SEPP = function(intercept, M, A, threshold, n, vzero = NULL, ...){
  #X is the data matrix with horizontal axis being time
  if(is.null(vzero)){
    vthreshold = rpois(M, intercept)
  }else{
      vthreshold = vzero
      vthreshold[which(vzero>threshold)] = threshold
  }
  X = matrix(0, ncol = n, nrow = M)
  X[,1] = rpois(M, lambda=exp(intercept + A %*% as.matrix(vthreshold)))
  for(t in 2:n){
    X.temp = X[,t-1]
    X.temp[which(X[,t-1]>threshold)] = threshold
    X[,t] = rpois(M, lambda = exp(intercept + A %*% X.temp))
  }
  return(X)
}