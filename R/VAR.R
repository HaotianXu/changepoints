#' @title Simulate a VAR1 model (without change point)
#' @description TO DO
#' @param sigma      A \code{numeric} scalar representing error standard deviation.
#' @param p          A \code{integer} scalar representing dimensionality.
#' @param n          A \code{integer} scalar representing sample size.
#' @param A          A \code{numeric} p-by-p matrix representing transition matrix of VAR1 model.
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
simu.change.VAR1= function(sigma, p, n, A, vzero = NULL, ...){
  X = matrix(0, nrow = p, ncol = n)
  if(is.null(vzero)){
    X[,1] = rnorm(p, mean = 0, sd = sigma)
  }else if(length(vzero) != p){
    stop("If the observation at time 0 (vzero) is specified, it should be a p-dim vector.")
  }else{
    X[,1] = vzero
  }
  for (t in 2:n){
    X[,t] = rnorm(p, mean = A%*%X[,t-1], sd = sigma)
  }
  return(X)
}


#' @title Internal Function: Prediction error in squared Frobenius norm for the lasso estimator of trasition matrix[11] 
#' @param s         A \code{integer} scalar of starting index.
#' @param e         A \code{integer} scalar of ending index.
#' @param X_futu    A \code{numeric} matrix of time series at one step ahead.
#' @param X_curr    A \code{numeric} matrix of time series at current step.
#' @param lambda    A \code{numeric} scalar of lasso penalty.
#' @param delta     A \code{integer} scalar of minimum spacing.
#' @return  A \code{numeric} scalar of prediction error in Frobenius norm.
#' @noRd
error.pred.seg.VAR1 = function(s, e, X_futu, X_curr, lambda, delta){
  #options(warn = -1)
  n = ncol(X_curr)
  p = nrow(X_curr)
  if(e > n | s > e | s < 1){
    stop("s and e are not correctly specified.")
  }
  if(e-s > delta){
    estimate = matrix(0, nrow=p, ncol=p)
    for(m in 1:p){
      # obtain Lasso estimator of each row of transition matrix
      out = glmnet(x = t(X_curr[,s:e]), y = X_curr[m,s:e], family=c("gaussian"),
                   alpha = 1, lambda = lambda/sqrt(e-s))#,intercept=F)
      estimate[m,] = as.vector(out$beta)
    }
    #norm(estimate-A1, type="F")
    #norm(A1, type="F")
    X_futu_hat = estimate%*%X_curr[,s:e]
    d = norm(X_futu_hat - X_futu[,s:e], type="F")
  }else{
    d = Inf
  }
  return(d^2)
}


# #dp matrix  function
# dp.matrix.function=function(X.train, Y.train, lambda.lasso , delta){
#   N = ncol(X.train)
#   p = nrow(X.train)
#   dp.matrix = matrix(Inf, N, N) 
#   for (s in 1:N){
#     for (e in 1:N){
#       if (e-s > delta){
#         #print(c(i,j))
#         dp.matrix[s,e] = error.pred.seg.VAR1(s, e, X_futu = Y.train, X_curr = X.train, lambda.lasso)
#       }
#     }
#   }
#   return(dp.matrix)
# }


#' @title Dynamic programming for VAR1 change points detection by l0 penalty
#' @description TO DO
#' @param DATA      A \code{numeric} matrix of observations.
#' @param gamma     A \code{numeric} scalar of the tuning parameter associated with the l0 penalty.
#' @param delta     A strictly \code{integer} scalar of minimum spacing.
#' @param lambda    A \code{numeric} scalar of tuning parameter for lasso penalty.
#' @param ...      Additional arguments.
#' @return TO DO.
#' @export
#' @author Haotian Xu
#' @examples
#' p = 20
#' sigma = 1
#' n = 100
#' A = matrix(rnorm(p*p), nrow = p)
#' DATA = simu.change.VAR1(sigma, p, n, A)
#' parti = DP.VAR1(DATA, gamma = 1, delta = 5, lambda = 1)$partition
#' part2local(parti)
DP.VAR1 = function(DATA, gamma, delta, lambda, ...){
  N = ncol(DATA)
  p = nrow(DATA)
  X_curr = DATA[,1:(N-1)]
  X_futu = DATA[,2:N]
  bestvalue = rep(0,N)
  partition = rep(0,N-1)
  bestvalue[1] = -gamma
  for(r in 1:(N-1)){
    bestvalue[r+1] = Inf
    for(l in 1:r){
      b = bestvalue[l] + gamma + error.pred.seg.VAR1(l, r, X_futu, X_curr, lambda, delta)
      if(b < bestvalue[r+1]){
        bestvalue[r+1] = b
        partition[r] = l-1
      }
    }
  }
  r = N-1
  l = partition[r]
  while(r > 0){
    r = l
    l = partition[r]
  }
  return(list(partition = partition))
}
