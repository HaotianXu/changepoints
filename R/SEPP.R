#' @title Simulate a (stable) SEPP model (without change point).
#' @description Simulate a (stable) SEPP model (without change point).
#' @param intercept  A \code{numeric} scalar representing the intercept of the model.
#' @param M          An \code{integer} scalar representing dimensionality.
#' @param A          A \code{numeric} M-by_M matrix representing the coefficient matrix.
#' @param threshold  A \code{numeric} scalar representing the upper bound for each coordinate of X_t (for stability).
#' @param n          An \code{integer} scalar representing sample size.
#' @param vzero      A \code{numeric} vector representing the observation at time 0.
#' @param ...        Additional arguments.
#' @return  A p-by-n matrix.
#' @export
#' @author  Daren Wang & Haotian Xu
#' @references Wang, D., Yu, Y., & Willett, R. (2020). Detecting Abrupt Changes in High-Dimensional Self-Exciting Poisson Processes. arXiv preprint arXiv:2006.03572.
#' @examples
#' M = 30
#' n = 300
#' threshold = 6 #thresholding makes the process stable
#' intercept = 1/2 #intercept of the model. Assume to be known as in the existing literature
#' v1 = (2*(seq(1,30,1)%%2)-1)*(0.05+0.05*6)
#' v2 = -v1
#' AA=matrix(0,nrow = 30,ncol=28)
#' A=cbind(v1,v2,AA)
#' A2=cbind(v2,v1,AA)
#' DATA1 = simu.SEPP(intercept, M, A, threshold, n, vzero = NULL)
#' DATA2 = simu.SEPP(intercept, M, A2, threshold, n, vzero = DATA1[,n])
#' DATA = cbind(DATA1, DATA2)
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


#' @title Dynamic programming for SEPP change points detection through l0 penalty.
#' @description Perform dynamic programming for SEPP change points detection through l0 penalty.
#' @param DATA      A \code{numeric} matrix of observations.
#' @param gamma     A \code{numeric} scalar of the tuning parameter associated with the l0 penalty.
#' @param delta     A positive \code{integer} scalar of minimum spacing.
#' @param lambda    A \code{numeric} scalar of tuning parameter for lasso penalty.
#' @param threshold A \code{numeric} scalar representing the upper bound for each coordinate of X_t (for stability).
#' @param ...      Additional arguments.
#' @return A vector of the best partition.
#' @export
#' @author Daren Wang & Haotian Xu
#' @references Wang, D., Yu, Y., & Willett, R. (2020). Detecting Abrupt Changes in High-Dimensional Self-Exciting Poisson Processes. arXiv preprint arXiv:2006.03572.
#' @examples
#' parti = DP.SEPP(gamma = 1, delta = 5, delta2 = 600, DATA, intercept = 1/2, lambda = 100, threshold = 6)$partition
#' cpt_hat = part2local(parti)
DP.SEPP = function(gamma, delta, delta2, DATA, intercept, lambda, threshold, ...){
  M = nrow(DATA)
  N = ncol(DATA)
  bestvalue = rep(0,N+1)
  partition = rep(0,N)
  bestvalue[1] = -gamma
  for(r in 1:N){
    bestvalue[r+1] = Inf
    for(l in 1:r){
      b = bestvalue[l] + gamma + error.seg.SEPP(l, r, DATA, intercept, lambda, threshold, delta, delta2)$dist
      if(b < bestvalue[r+1]){
        bestvalue[r+1] = b
        partition[r] = l-1
      }
    }
  }
  r = N
  l = partition[r]
  while(r > 0){
    r = l
    l = partition[r]
  }
  return(list(partition = partition))
}


#' @title Internal Function: Prediction error in terms of poisson log-likelihood function for the lasso estimator of transition matrix[12] 
#' @param s         A \code{integer} scalar of starting index.
#' @param e         A \code{integer} scalar of ending index.
#' @param DATA      A \code{integer} matrix of observations with horizontal axis being time.
#' @param intercept A \code{numeric} scalar representing the intercept of the model. Assume to be known as in the existing literature
#' @param lambda    A \code{numeric} scalar of lasso penalty.
#' @param threshold A \code{numeric} scalar representing the upper bound for each coordinate of X_t (for stability).
#' @param delta     A \code{integer} scalar of minimum spacing.
#' @param delta2    A \code{integer} scalar representing the maximal of the change point spacing (for reducing computation cost).
#' @return  A \code{numeric} scalar of prediction error in terms of poisson log-likelihood function.
#' @noRd
error.seg.SEPP = function(s, e, DATA, intercept, lambda, threshold, delta, delta2, ...){
  M = nrow(DATA)
  estimate = matrix(0, nrow = M, ncol = M)
  n.temp = e - s + 1
  d = rep(NA, M)
  DATA.x = DATA
  DATA.x[which(DATA.x > threshold)] = threshold
  if(n.temp - 1 > delta & n.temp - 1 < delta2){
    for(m in 1:M){
      pen <- penalized(DATA[m, (s+1):e] / exp(intercept), 
                       penalized = t(DATA.x[, s:(e-1)]), unpenalized = ~0,
                       lambda1 = lambda/(n.temp)^(1/2), lambda2 = 10, model = c("poisson"), trace = F, maxiter = 500) # lambda is rescaled
      estimate[m,] = coefficients(pen, "all")
      d[m] = sum(abs(sapply((s+1):e, function(x) exp(intercept + estimate[m,] %*% DATA.x[,x-1]) - DATA[m,x] * (intercept + estimate[m,] %*% DATA.x[,x-1]))))
    }
  }else{
    estimate = NA
    d = Inf
  }
  result = list(dist = sum(d), transition.hat = estimate)
  return(result)
}
