#' @title Dynamic programming for regression change points detection by l0 penalty
#' @description TO DO
#' @param gamma     A \code{numeric} scalar of the tuning parameter associated with the l0 penalty.
#' @param delta     A strictly \code{integer} scalar of minimum spacing.
#' @param y         A \code{numeric} vector of observations.
#' @param X         A \code{numeric} matrix of covariates. If missing, then the mean change of y is considered.
#' @param ...      Additional arguments.
#' @return TO DO.
#' @export
#' @author Haotian Xu
#' @examples
#' data = simu.change.regression(10, c(10, 30, 40, 70, 90), 30, 100, 1, 9)
#' DP.regression(2, 5, data$y, X = data$X, lambda = 1)
DP.regression = function(gamma, delta, y, X, lambda, ...){
  N = length(y)
  bestvalue = rep(0,N+1)
  partition = rep(0,N)
  yhat = rep(NA, N)
  if(length(dim(X)) != 2 | dim(X)[2] != N){
    stop("X should be a p-by-n design matrix")
  }
  p = dim(X)[1]
  bestvalue[1] = -gamma*log(max(N,p))
  for(r in 1:N){
    bestvalue[r+1] = Inf
    for(l in 1:r){
      b = bestvalue[l] + gamma*log(max(N,p)) + error.pred.seg.regression(l, r, y, X, lambda, delta)
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


#' @title Simulate sparse regression model with changepoints in coefficients (Model1 in [10]) under the setting of Simulations 4.2 [10]
#' @description TO DO
#' @param d0         A \code{numeric} scalar of number of nonzero coefficients.
#' @param cpt.ture   A \code{integer} vector of true changepoints (sorted in strictly increasing order).
#' @param p          A \code{integer} scalar of dimensionality.
#' @param n          A \code{integer} scalar of sample size.
#' @param sigma      A \code{numeric} scalar of error standard deviation.
#' @param kappa      A \code{numeric} scalar of minimum jump size in terms of the l2 norm.
#' @param ...        Additional arguments.
#' @return A \code{list} with the structure:
#' \itemize{
#'  \item cpt.true    A vector of true changepoints (sorted in strictly increasing order).
#'  \item X           A p-by-n matrix of (random) design matrix.
#'  \item y           A n-dim vector of response variable.
#'  \item betafullmat A p-by-n matrix of coefficients.
#'  \item ...         Additional parameters.
#' }
#' @export
#' @author 
#' @examples
#' TO DO
simu.change.regression = function(d0, cpt.true, p, n, sigma, kappa){
  #seed = 10
  if(d0 >= p){
    stop("d0 should be strictly smaller than p")
  }
  if(sigma <= 0){
    stop("sigma should be strictly larger than 0")
  }
  if(kappa <= 0){
    stop("kappa should be strictly larger than 0")
  }
  no.cpt = length(cpt.true)
  if(is.unsorted(cpt.true, strictly = TRUE) | min(cpt.true) <= 1 | max(cpt.true >= n) | no.cpt > n-2){
    stop("cpt.true is not correctly specified")
  }
  X = matrix(rnorm(p*n,0,1), p, n)
  y = matrix(0, n, 1)
  nonzero.element.loc = c(1:d0)
  cpt = c(0, cpt.true, n)
  beta = matrix(0, p, no.cpt+1)
  betafullmat = matrix(0, p, n)
  for (i in 1:(no.cpt+1)){
    if (i%%2 == 1){
      beta[nonzero.element.loc, i] = kappa/(2*sqrt(d0))
    }
    else{
      beta[nonzero.element.loc, i] = -kappa/(2*sqrt(d0))
    }
    y[(1+cpt[i]):cpt[i+1],] = rnorm(cpt[i+1] - cpt[i], t(X[,(1+cpt[i]):cpt[i+1]]) %*% beta[,i], sigma)
    for (j in (1+cpt[i]):cpt[i+1]){
      betafullmat[,j] = beta[,i] 
    }
  }
  List = list(cpt.true = cpt.true, X = X, y = y, betafullmat = betafullmat)
  return(List)
}


#' @title Internal Function: Prediction error in squared l2 norm for the lasso [10] 
#' @param s         A \code{integer} scalar of starting index.
#' @param e         A \code{integer} scalar of ending index.
#' @param y         A \code{numeric} vector of response variable.
#' @param X         A \code{numeric} matrix of covariates.
#' @param lambda    A \code{numeric} scalar of lasso penalty.
#' @param delta     A \code{integer} scalar of minimum spacing.
#' @return  A \code{numeric} scalar of prediction error in l2 norm.
#' @noRd
error.pred.seg.regression = function(s, e, y, X, lambda, delta){
  n = ncol(X)
  p = nrow(X)
  if(e > n | s > e | s < 1){
    stop("s and e are not correctly specified.")
  }
  if (e-s > delta){
    fit = glmnet(x = t(X[,s:e]), y = y[s:e], family = c("gaussian"),
                 alpha = 1, lambda = lambda*sqrt(max(log(max(n,p)), e-s))*sqrt(log(max(n,p)))/(e-s),intercept=F)
    coef_est = as.vector(fit$beta)
    yhat = t(X[,s:e])%*%coef_est
    d = norm(y[s:e] - yhat, type = "2")
  }
  else{
    d = Inf
  }
  return(d^2)
}


#' @title Local refinement [10] 
#' @description TO DO
#' @param cpt.init  A \code{integer} vector of initial changepoints estimation (sorted in strictly increasing order).
#' @param y         A \code{numeric} vector of response variable.
#' @param X         A \code{numeric} matrix of covariates.
#' @param zeta      A \code{numeric} scalar of lasso penalty.
#' @param w         A \code{numeric} scalar of weight for interpolation.
#' @param ...       Additional arguments.
#' @return  A \code{numeric} scalar of prediction error in l2 norm.
#' @export
#' @author 
#' @examples
#' data = simu.change.regression(10, c(10, 30, 40, 70, 90), 30, 100, 1, 9)
#' cpt.init = part2local(DP.regression(2, 5, data$y, X = data$X, lambda = 2)$partition)$cpt
#' local.refine(cpt.init, data$y, X = data$X, 1, 1/3)
local.refine = function(cpt.init, y, X, zeta, w = 1/3){
  n = ncol(X)
  cpt.init.ext = c(0, cpt.init, n)
  cpt.init.numb = length(cpt.init)
  cpt.refined = rep(0, cpt.init.numb+1)
  for (k in 1:cpt.init.numb){
    s = w*cpt.refined[k] + (1-w)*cpt.init.ext[k+1]
    e = (1-w)*cpt.init.ext[k+1] + w*cpt.init.ext[k+2]
    lower = ceiling(s) + 1
    upper = floor(e) - 1
    b = sapply(lower:upper, function(eta)obj.func.lr(s, e, eta, y, X, zeta))
    cpt.refined[k+1] = ceiling(s) + which.min(b)
  }
  return(cpt.refined[-1])
}


#' @title Internal Function: An objective function to select the best splitting location in the local refinement, see eq(4) in [10]
#' @param s         A \code{integer} scalar of starting index.
#' @param e         A \code{integer} scalar of ending index.
#' @param y         A \code{numeric} vector of response variable.
#' @param X         A \code{numeric} matrix of covariates.
#' @param zeta      A \code{numeric} scalar of tuning parameter for the group lasso.
#' @noRd
obj.func.lr = function(s, e, eta, y, X, zeta){
  n = ncol(X)
  p = nrow(X)
  group = rep(1:p, 2)
  X.convert = X.glasso.converter(X[,(ceiling(s)):(floor(e))], eta, ceiling(s))
  y.convert = y[(ceiling(s)):(floor(e)),]
  lambda.LR = zeta*sqrt(log(max(n, p)))
  auxfit = gglasso(x = X.convert, y = y.convert, group = group, loss="ls",
                   lambda = lambda.LR/(floor(e)-ceiling(s)+1), intercept = FALSE, eps = 0.001)
  coef = as.vector(auxfit$beta)
  coef1 = coef[1:p]
  coef2 = coef[(p+1):(2*p)]
  btemp = norm(y.convert - X.convert %*% coef, type = "2")^2 + lambda.LR*sum(sqrt(coef1^2 + coef2^2))
  return(btemp)
}


#' @title Internal Function: Convert a p-by-n design submatrix X with partial consecutive observations into a n-by-(2p) matrix, which fits the group lasso, see eq(4) in [10]
#' @param  X         A \code{numeric} matrix of covariates with partial consecutive observations.
#' @param  eta       A \code{integer} scalar of splitting index.
#' @param  s_ceil    A \code{numeric} matrix of covariates.
#' @return A n-by-(2p) matrix
#' @noRd
X.glasso.converter=function(X, eta, s_ceil){
  n = ncol(X)
  xx1 = xx2 = t(X)
  t = eta - s_ceil + 1
  xx1[(t+1):n,] = 0
  xx2[1:t,] = 0
  xx = cbind(xx1/sqrt(t-1), xx2/sqrt(n-t+1))
  return(xx)
}