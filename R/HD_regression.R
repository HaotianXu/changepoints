#' @title Dynamic programming algorithm for regression change points detection through l0 penalty.
#' @description     Perform dynamic programming algorithm for regression change points detection through l0 penalty.
#' @param y         A \code{numeric} vector of observations.
#' @param X         A \code{numeric} matrix of covariates.
#' @param gamma     A positive \code{numeric} scalar of the tuning parameter associated with the l0 penalty.
#' @param lambda    A positive \code{numeric} scalar of tuning parameter for the lasso penalty.
#' @param delta     A positive \code{integer} scalar of minimum spacing.
#' @param ...       Additional arguments.
#' @return          A vector of the best partition.
#' @export
#' @author
#' @examples
#' data = simu.change.regression(10, c(10, 30, 40, 70, 90), 30, 100, 1, 9)
#' temp = DP.regression(data$y, X = data$X, gamma = 2, lambda = 1, delta = 5)
#' part2local(temp$partition)
DP.regression = function(y, X, gamma, lambda, delta, ...){
  N = length(y)
  bestvalue = rep(0,N+1)
  partition = rep(0,N)
  #yhat = rep(NA, N)
  if(length(dim(X)) != 2 | dim(X)[2] != N){
    stop("X should be a p-by-n design matrix")
  }
  p = dim(X)[1]
  bestvalue[1] = -gamma*log(max(N,p))
  for(r in 1:N){
    bestvalue[r+1] = Inf
    for(l in 1:r){
      b = bestvalue[l] + gamma*log(max(N,p)) + error.pred.seg.regression(l, r, y, X, lambda, delta)$MSE
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


#' @title Simulate sparse regression model with changepoints in coefficients.
#' @description      Simulate sparse regression model with changepoints in coefficients under the setting of Simulations 4.2 [10].
#' @param d0         A \code{numeric} scalar of number of nonzero coefficients.
#' @param cpt.ture   An \code{integer} vector of true changepoints (sorted in strictly increasing order).
#' @param p          An \code{integer} scalar of dimensionality.
#' @param n          An \code{integer} scalar of sample size.
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
#' data = simu.change.regression(10, c(10, 30, 40, 70, 90), 30, 100, 1, 9)
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
#' @param s         An \code{integer} scalar of starting index.
#' @param e         An \code{integer} scalar of ending index.
#' @param y         A \code{numeric} vector of response variable.
#' @param X         A \code{numeric} matrix of covariates.
#' @param lambda    A \code{numeric} scalar of tuning parameter for lasso penalty.
#' @param delta     A \code{integer} scalar of minimum spacing.
#' @return    A \code{list} with the structure:
#' \itemize{
#'  \item MSE       A \code{numeric} scalar of prediction error in l2 norm.
#'  \item beta.hat  A p-dim vector of estimated coefficients.
#'  \item ...         Additional parameters.
#' }
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
    coef_est = NA
    d = Inf
  }
  result = list(MSE = d^2, beta.hat = coef_est)
  return(result)
}


#' @title Local refinement for regression change points detection.
#' @description     Perform local refinement for regression change points detection.
#' @param cpt.init  An \code{integer} vector of initial changepoints estimation (sorted in strictly increasing order).
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
#' cpt.init = part2local(DP.regression(data$y, X = data$X, gamma = 2, lambda = 2, delta = 5)$partition)
#' local.refine.regression(cpt.init, data$y, X = data$X, 1, 1/3)
local.refine.regression = function(cpt.init, y, X, zeta, w = 1/3){
  n = ncol(X)
  cpt.init.ext = c(0, cpt.init, n)
  cpt.init.numb = length(cpt.init)
  cpt.refined = rep(0, cpt.init.numb+1)
  for (k in 1:cpt.init.numb){
    s = w*cpt.init.ext[k] + (1-w)*cpt.init.ext[k+1]
    e = (1-w)*cpt.init.ext[k+1] + w*cpt.init.ext[k+2]
    lower = ceiling(s) + 1
    upper = floor(e) - 1
    b = sapply(lower:upper, function(eta)obj.func.lr.regression(s, e, eta, y, X, zeta))
    cpt.refined[k+1] = ceiling(s) + which.min(b)
  }
  return(cpt.refined[-1])
}


#' @title Internal Function: An objective function to select the best splitting location in the local refinement, see eq(4) in [10]
#' @param s.inter   A \code{numeric} scalar of interpolated starting index.
#' @param e.inter   A \code{numeric} scalar of interpolated ending index.
#' @param y         A \code{numeric} vector of response variable.
#' @param X         A \code{numeric} matrix of covariates.
#' @param zeta      A \code{numeric} scalar of tuning parameter for the group lasso.
#' @noRd
obj.func.lr.regression = function(s.inter, e.inter, eta, y, X, zeta){
  n = ncol(X)
  p = nrow(X)
  group = rep(1:p, 2)
  X.convert = X.glasso.converter.regression(X[,(ceiling(s.inter)):(floor(e.inter))], eta, ceiling(s.inter))
  y.convert = y[(ceiling(s.inter)):(floor(e.inter))]
  lambda.LR = zeta*sqrt(log(max(n, p)))
  auxfit = gglasso(x = X.convert, y = y.convert, group = group, loss="ls",
                   lambda = lambda.LR/(floor(e.inter)-ceiling(s.inter)+1), intercept = FALSE, eps = 0.001)
  coef = as.vector(auxfit$beta)
  coef1 = coef[1:p]
  coef2 = coef[(p+1):(2*p)]
  btemp = norm(y.convert - X.convert %*% coef, type = "2")^2 + lambda.LR*sum(sqrt(coef1^2 + coef2^2))
  return(btemp)
}


#' @title Internal Function: Convert a p-by-n design submatrix X with partial consecutive observations into a n-by-(2p) matrix, which fits the group lasso, see eq(4) in [10]
#' @param  X         A \code{numeric} matrix of covariates with partial consecutive observations.
#' @param  eta       A \code{integer} scalar of splitting index.
#' @param  s_ceil    A \code{integer} scalar of starting index.
#' @return A n-by-(2p) matrix
#' @noRd
X.glasso.converter.regression = function(X, eta, s_ceil){
  n = ncol(X)
  xx1 = xx2 = t(X)
  t = eta - s_ceil + 1
  xx1[(t+1):n,] = 0
  xx2[1:t,] = 0
  xx = cbind(xx1/sqrt(t-1), xx2/sqrt(n-t))
  return(xx)
}



#' @title Internal function: Cross-validation of dynamic programming algorithm for regression change points detection through l0 penalty.
#' @description     Perform cross-validation of dynamic programming algorithm for regression change points detection through l0 penalty.
#' @param y         A \code{numeric} vector of observations.
#' @param X         A \code{numeric} matrix of covariates.
#' @param gamma     A \code{numeric} scalar of the tuning parameter associated with the l0 penalty.
#' @param lambda    A \code{numeric} scalar of tuning parameter for the lasso penalty.
#' @param delta     A strictly \code{integer} scalar of minimum spacing.
#' @param ...      Additional arguments.
#' @return  A \code{list} with the structure:
#' \itemize{
#'  \item{cpt_hat}: A list of vector of estimated change points locations (sorted in strictly increasing order).
#'  \item{K_hat}: A list of scalar of number of estimated change points.
#'  \item{test_error}: A list of vector of testing errors.
#'  \item{train_error}: A list of vector of training errors.
#' } 
#' @noRd
CV.DP.regression = function(y, X, gamma, lambda, delta, ...){
  N = ncol(X)
  even_indexes = seq(2, N, 2)
  odd_indexes = seq(1, N, 2)
  train.X = X[,odd_indexes]
  train.y = y[odd_indexes]
  validation.X = X[,even_indexes]
  validation.y = y[even_indexes]
  init_cpt_train = part2local(DP.regression(train.y, train.X, gamma, lambda, delta)$partition)
  init_cpt_train.long = c(0, init_cpt_train, ncol(train.X))
  diff.point = diff(init_cpt_train.long)
  if (length(which(diff.point == 1)) > 0){
    print(paste("gamma =", gamma,",", "lambda =", lambda, ".","Warning: Consecutive points detected. Try a larger gamma."))
    init_cpt = odd_indexes[init_cpt_train]
    len = length(init_cpt)
    result = list(cpt_hat = init_cpt, K_hat = len, test_error = Inf, train_error = Inf)
  }
  else{
    init_cpt = odd_indexes[init_cpt_train]
    len = length(init_cpt)
    init_cpt_long = c(init_cpt_train, floor(N/2))
    interval = matrix(0, nrow = len+1, ncol = 2)
    interval[1,] = c(1, init_cpt_long[1])
    if(len > 0){
      for(j in 2:(1+len)){
        interval[j,] = c(init_cpt_long[j-1]+1, init_cpt_long[j])
      }
    }
    p = nrow(train.X)
    trainmat = sapply(1:(len+1), function(index) error.pred.seg.regression(interval[index,1], interval[index,2], train.y, train.X, lambda, delta))
    betamat = matrix(0, nrow = p, ncol = len+1)
    training_loss = matrix(0, nrow = 1, ncol = len+1)
    for(col in 1:(len+1)){
      betamat[,col] = as.numeric(trainmat[2,col]$beta.hat)
      training_loss[,col] = as.numeric(trainmat[1,col]$MSE)
    }
    validationmat = sapply(1:(len+1), function(index) error.test.regression(interval[index,1], interval[index,2], validation.y, validation.X, betamat[,index]))
    result = list(cpt_hat = init_cpt, K_hat = len, test_error = sum(validationmat), train_error = sum(training_loss))
  }
  return(result)
}


#' @title Internal Function: compute testing error for regression
#' @param  lower     A \code{integer} scalar of starting index.
#' @param  upper     A \code{integer} scalar of ending index.
#' @param y         A \code{numeric} vector of observations.
#' @param X         A \code{numeric} matrix of covariates.
#' @return A numeric scalar of testing error in squared l2 norm.
#' @noRd
error.test.regression = function(lower, upper, y, X, beta.hat){
  res = norm(y[lower:upper] - t(X[,lower:upper])%*%beta.hat, type = "2")^2
  return(res)
} 


#' @title Grid search based on Cross-Validation of Dynamic Programming for regression change points detection via l0 penalty
#' @description TO DO
#' @param y             A \code{numeric} vector of observations.
#' @param X             A \code{numeric} matrix of covariates.
#' @param gamma.set     A \code{numeric} vector of candidate tuning parameter associated with the l0 penalty.
#' @param lambda.set    A \code{numeric} vector of candidate tuning parameter for the lasso penalty.
#' @param delta         A strictly \code{integer} scalar of minimum spacing.
#' @param ...           Additional arguments.
#' @return  A \code{list} with the structure:
#' \itemize{
#'  \item{cpt_hat}: A list of vector of estimated change points locations (sorted in strictly increasing order).
#'  \item{K_hat}: A list of scalar of number of estimated change points.
#'  \item{test_error}: A list of vector of testing errors.
#'  \item{train_error}: A list of vector of training errors.
#' } 
#' @export
#' @author Daren Wang
#' @examples
#' data = simu.change.regression(10, c(10, 30, 40, 70, 90), 30, 100, 1, 9)
#' CV.search.DP.regression(data$y, data$X, gamma.set = 1:5, lambda.set = 1:5, delta = 5)
#' cpt.init = part2local(DP.regression(data$y, X = data$X, gamma = 4, lambda = 2, delta = 5)$partition)
#' local.refine.regression(cpt.init, data$y, X = data$X, 1, 1/3)
CV.search.DP.regression = function(y, X, gamma.set, lambda.set, delta){
  output = sapply(1:length(lambda.set), function(i) sapply(1:length(gamma.set), 
                                                           function(j) CV.DP.regression(y, X, gamma.set[j], lambda.set[i], delta)))
  print(output)
  cpt_hat = output[seq(1,4*length(gamma.set),4),]## estimated change points
  K_hat = output[seq(2,4*length(gamma.set),4),]## number of estimated change points
  test_error = output[seq(3,4*length(gamma.set),4),]## validation loss
  train_error = output[seq(4,4*length(gamma.set),4),]## training loss                                                      
  result = list(cpt_hat = cpt_hat, K_hat = K_hat, test_error = test_error, train_error = train_error)
  return(result)
}
