#' @title Simulate from a VAR1 model (without change point).
#' @description Simulate data of size n and dimension p from a VAR1 model (without change point) with Gaussian i.i.d. error terms.
#' @param sigma      A \code{numeric} scalar representing the standard deviation of error terms.
#' @param p          An \code{integer} scalar representing dimension.
#' @param n          An \code{integer} scalar representing sample size.
#' @param A          A \code{numeric} p-by-p matrix representing the transition matrix of the VAR1 model.
#' @param vzero      A \code{numeric} vector representing the observation at time 0.
#' @param ...        Additional arguments.
#' @return  A p-by-n matrix.
#' @export
#' @author  Daren Wang
#' @examples
#' p = 20
#' sigma = 1
#' n = 100
#' A = matrix(rnorm(p*p), nrow = p)
#' simu.VAR1(sigma, p, n, A)
simu.VAR1= function(sigma, p, n, A, vzero = NULL, ...){
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


#' @export
error_pred_seg_VAR1 <- function(X_futu, X_curr, s, e, alpha, lambda, delta) {
  .Call('_changepoints_rcpp_error_pred_seg_VAR1', PACKAGE = 'changepoints', X_futu, X_curr, s, e, alpha, lambda, delta)
}


#' @title Internal Function: Prediction error in squared Frobenius norm for the lasso estimator of transition matrix[11] 
#' @param s         A \code{integer} scalar of starting index.
#' @param e         A \code{integer} scalar of ending index.
#' @param X_futu    A \code{numeric} matrix of time series at one step ahead.
#' @param X_curr    A \code{numeric} matrix of time series at current step.
#' @param lambda    A \code{numeric} scalar of lasso penalty.
#' @param delta     A \code{integer} scalar of minimum spacing.
#' @return  A \code{numeric} scalar of prediction error in Frobenius norm.
#' @noRd
error.pred.seg.VAR1 = function(s, e, X_futu, X_curr, lambda, delta){
  n = ncol(X_curr)
  p = nrow(X_curr)
  if(e > n | s > e | s < 1){
    stop("s and e are not correctly specified.")
  }
  if(e-s > delta){
    estimate = as.matrix(do.call(cbind, sapply(1:p, function(m) glmnet(x = t(X_curr[,s:e]), y = X_curr[m,s:e], family=c("gaussian"), alpha = 1, lambda = lambda/sqrt(e-s))$beta)))
    X_futu_hat = t(estimate) %*% X_curr[,s:e]
    d = norm(X_futu_hat - X_futu[,s:e], type = "F")
  }else{
    estimate = NA
    d = Inf
  }
  result = list(MSE = d^2, transition.hat = estimate)
  return(result)
}


#' @export
DP_VAR1 <- function(X_futu, X_curr, alpha, gamma, lambda, delta) {
  .Call('_changepoints_rcpp_DP_VAR1', PACKAGE = 'changepoints', X_futu, X_curr, alpha, gamma, lambda, delta)
}

#' @title Dynamic programming for VAR1 change points detection through l0 penalty.
#' @description Perform dynamic programming for VAR1 change points detection through l0 penalty.
#' @param X_futu    A \code{numeric} matrix of time series at one step ahead.
#' @param X_curr    A \code{numeric} matrix of time series at current step.
#' @param gamma     A \code{numeric} scalar of the tuning parameter associated with the l0 penalty.
#' @param lambda    A \code{numeric} scalar of tuning parameter for lasso penalty.
#' @param delta     A strictly \code{integer} scalar of minimum spacing.
#' @param ...      Additional arguments.
#' @return partition: A vector of the best partition.
#' @export
#' @author Daren Wang, Haotian Xu
#' @examples
#' p = 10
#' sigma = 1
#' n = 5
#' v1 = 2*(seq(1,p,1)%%2) - 1
#' v2 = -v1
#' AA = matrix(0, nrow = p, ncol = p-2)
#' A1=cbind(v1,v2,AA)
#' A2=cbind(v2,v1,AA)
#' A3=A1
#' data = simu.VAR1(sigma, p, 2*n+1, A1)
#' data = cbind(data, simu.VAR1(sigma, p, 2*n, A2, vzero=c(data[,ncol(data)])))
#' data = cbind(data, simu.VAR1(sigma, p, 2*n, A3, vzero=c(data[,ncol(data)])))
#' N = ncol(data)
#' X_curr = data[,1:(N-1)]
#' X_futu = data[,2:N]
#' parti = DP.VAR1(X_futu, X_curr, gamma = 1, lambda = 1, delta = 5)$partition
#' part2local(parti)
DP.VAR1 = function(X_futu, X_curr, gamma, lambda, delta, ...){
  p = nrow(X_futu)
  N = ncol(X_futu) + 1
  bestvalue = rep(0,N)
  partition = rep(0,N-1)
  bestvalue[1] = -gamma
  for(r in 1:(N-1)){
    bestvalue[r+1] = Inf
    for(l in 1:r){
      b = bestvalue[l] + gamma + error.pred.seg.VAR1(l, r, X_futu, X_curr, lambda, delta)$MSE
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



#' @title Local refinement for VAR1 change points detection [11] 
#' @description TO DO
#' @param cpt.init   A \code{integer} vector of initial changepoints estimation (sorted in strictly increasing order).
#' @param DATA       A \code{numeric} matrix of observations with with horizontal axis being time, and vertical axis being dimensions.
#' @param zeta.group A \code{numeric} scalar of lasso penalty.
#' @param w          A \code{numeric} scalar of weight for interpolation.
#' @param ...       Additional arguments.
#' @return  A \code{numeric} scalar of prediction error in l2 norm.
#' @export
#' @author 
#' @examples
#' #TODO data = simu.change.regression(10, c(10, 30, 40, 70, 90), 30, 100, 1, 9)
#' cpt.init = part2local(DP.regression(2, 5, data$y, X = data$X, lambda = 2)$partition)$cpt
#' local.refine.regression(cpt.init, data$y, X = data$X, 1, 1/3)
local.refine.VAR1 = function(cpt.init, DATA, zeta.group, w = 1/3, ...){
  N = ncol(DATA)
  p = nrow(DATA)
  X_curr = DATA[,1:(N-1)]
  X_futu = DATA[,2:N]
  cpt.init.ext = c(0, cpt.init, N)
  cpt.init.numb = length(cpt.init)
  cpt.refined = rep(0, cpt.init.numb+1)
  for (k in 1:cpt.init.numb){
    s.inter = w*cpt.init.ext[k] + (1-w)*cpt.init.ext[k+1]
    e.inter = (1-w)*cpt.init.ext[k+1] + w*cpt.init.ext[k+2]
    lower = ceiling(s.inter) + 1
    upper = floor(e.inter) - 1
    b = sapply(lower:upper, function(eta)obj.func.lr.VAR1(ceiling(s.inter), floor(e.inter), eta, X_futu, X_curr, zeta.group))
    cpt.refined[k+1] = ceiling(s.inter) + which.min(b)
  }
  return(cpt.refined[-1])
}


#' @title Internal Function: An objective function to select the best splitting location in the local refinement
#' @param s.inter    A \code{numeric} scalar of interpolated starting index.
#' @param e.inter    A \code{numeric} scalar of ending index.
#' @param eta        A \code{integer} scalar of splitting index.
#' @param X_futu     A \code{numeric} matrix of time series at one step ahead.
#' @param X_curr     A \code{numeric} matrix of time series at current step.
#' @param zeta.group A \code{numeric} scalar of tuning parameter for the group lasso.
#' @noRd
obj.func.lr.VAR1 = function(s.inter.ceil, e.inter.floor, eta, X_futu, X_curr, zeta.group){
  n = ncol(X_futu)
  p = nrow(X_futu)
  btemp = rep(NA, p)
  group = rep(1:p, 2)
  X.convert = X.glasso.converter.VAR1(X_curr[,s.inter.ceil:e.inter.floor], eta, s.inter.ceil)
  for(m in 1:p){
    y.convert = X_futu[m, s.inter.ceil:e.inter.floor]
    auxfit = gglasso(x = X.convert, y = y.convert, group = group, loss="ls",
                     lambda = zeta.group/(e.inter.floor-s.inter.ceil+1), intercept = FALSE, eps = 0.001)
    coef = as.vector(auxfit$beta)
    btemp[m] = norm(y.convert - X.convert %*% coef, type = "2")^2
  }
  return(sum(btemp))
}


#' @title Internal Function: Convert a p-by-n submatrix X with partial consecutive observations into a n-by-(2p) matrix, which fits the group lasso, see eq(7) in [11]
#' @param  X         A \code{numeric} matrix of partial consecutive observations.
#' @param  eta       A \code{integer} scalar of splitting index.
#' @param  s_ceil    A \code{integer} scalar of starting index.
#' @return   A n-by-(2p) matrix
#' @noRd
X.glasso.converter.VAR1 = function(X, eta, s_ceil){
  n = ncol(X)
  xx1 = xx2 = t(X)
  t = eta - s_ceil + 1
  xx1[(t+1):n,] = 0
  xx2[1:t,] = 0
  xx = cbind(xx1/sqrt(t-1), xx2/sqrt(n-t))
  return(xx)
}



#' @title Internal function: Cross-Validation of Dynamic Programming algorithm for VAR1 change points detection via l0 penalty
#' @param DATA      A \code{numeric} matrix of observations with with horizontal axis being time, and vertical axis being dimensions.
#' @param gamma     A positive \code{numeric} scalar of the tuning parameter associated with the l0 penalty.
#' @param lambda    A positive \code{numeric} scalar of tuning parameter for the lasso penalty.
#' @param delta     A positive \code{integer} scalar of minimum spacing.
#' @param ...      Additional arguments.
#' @noRd
CV.DP.VAR1 = function(DATA, gamma, lambda, delta, ...){
  DATA.temp = DATA
  if (ncol(DATA)%%2 == 0){
    DATA.temp = DATA[,2:ncol(DATA)]
  }
  N = ncol(DATA.temp)
  p = nrow(DATA.temp)
  X_curr = DATA.temp[,1:(N-1)]
  X_futu = DATA.temp[,2:N]
  X_curr.train = X_curr[,seq(1,N-1,2)]
  X_curr.test = X_curr[,seq(2,N-1,2)]
  X_futu.train = X_futu[,seq(1,N-1,2)]
  X_futu.test = X_futu[,seq(2,N-1,2)]
  init_cpt_train = part2local(DP.VAR1(X_futu.train, X_curr.train, gamma, lambda, delta)$partition)
  init_cpt = 2*init_cpt_train
  len = length(init_cpt)
  init_cpt_long = c(init_cpt_train, ncol(X_curr.train))
  interval = matrix(0, nrow = len+1, ncol = 2)
  interval[1,] = c(1, init_cpt_long[1])
  if(len > 0){
    for(j in 2:(1+len)){
      interval[j,] = c(init_cpt_long[j-1]+1, init_cpt_long[j])
    }
  }
  trainmat = sapply(1:(len+1), function(index) error.pred.seg.VAR1(interval[index,1], interval[index,2], X_futu.train, X_curr.train, lambda, delta))
  transition.list = vector("list", len+1)
  training_loss = matrix(0, nrow = 1, ncol = len+1)
  for(col in 1:(len+1)){
    transition.list[[col]] = trainmat[2,col]$transition.hat
    training_loss[,col] = as.numeric(trainmat[1,col]$MSE)
  }
  validationmat = sapply(1:(len+1), function(index) error.test.VAR1(interval[index,1], interval[index,2], X_futu.test, X_curr.test, transition.list[[index]]))
  result = list(cpt_hat = init_cpt, K_hat = len, test_error = sum(validationmat), train_error = sum(training_loss))
  return(result)
}


#' @title Internal Function: compute testing error for VAR1
#' @param  lower     A \code{integer} scalar of starting index.
#' @param  upper     A \code{integer} scalar of ending index.
#' @param  X_futu    A \code{numeric} matrix of observations.
#' @param  X_curr    A \code{numeric} matrix of covariates.
#' @param  transition.hat A \code{numeric} matrix of transition matrix estimator.
#' @return A numeric scalar of testing error in squared Frobenius norm.
#' @noRd
error.test.VAR1 = function(lower, upper, X_futu, X_curr, transition.hat){
  res = norm(X_futu[lower:upper] - transition.hat%*%X_curr[,lower:upper], type = "F")^2
  return(res)
}


#' @title Grid search based on Cross-Validation of Dynamic Programming for regression change points detection via l0 penalty
#' @description TO DO
#' @param DATA          A \code{numeric} matrix of observations with with horizontal axis being time, and vertical axis being dimensions.
#' @param gamma.set     A \code{numeric} vector of candidate tuning parameter associated with the l0 penalty.
#' @param lambda.set    A \code{numeric} vector of candidate tuning parameter for the lasso penalty.
#' @param delta         A strictly \code{integer} scalar of minimum spacing.
#' @param ...           Additional arguments.
#' @return Row: lambda.set; column: lambda.set
#' @export
#' @author
#' @examples
#' TO DO
CV.search.DP.VAR1 = function(DATA, gamma.set, lambda.set, delta, ...){
  output = sapply(1:length(lambda.set), function(i) sapply(1:length(gamma.set), 
                                                           function(j) CV.DP.VAR1(DATA, gamma.set[j], lambda.set[i], delta)))
  cpt_hat = output[seq(1,4*length(gamma.set),4),]## estimated change points
  K_hat = output[seq(2,4*length(gamma.set),4),]## number of estimated change points
  test_error = output[seq(3,4*length(gamma.set),4),]## validation loss
  train_error = output[seq(4,4*length(gamma.set),4),]## training loss                                                      
  result = list(cpt_hat = cpt_hat, K_hat = K_hat, test_error = test_error, train_error = train_error)
  return(result)
}


#' @title Internal function: Cross-Validation of local refinement for VAR1 change points detection via l0 penalty
#' @param DATA      A \code{numeric} matrix of observations with with horizontal axis being time, and vertical axis being dimensions.
#' @param gamma     A positive \code{numeric} scalar of the tuning parameter associated with the l0 penalty.
#' @param lambda    A positive \code{numeric} scalar of tuning parameter for the lasso penalty.
#' @param delta     A positive \code{integer} scalar of minimum spacing.
#' @param ...      Additional arguments.
#' @noRd
CV.lr.VAR1 = function(DATA, gamma, lambda, delta, ...){
  DATA.temp = DATA
  if (ncol(DATA)%%2 == 0){
    DATA.temp = DATA[,2:ncol(DATA)]
  }
  N = ncol(DATA.temp)
  p = nrow(DATA.temp)
  X_curr = DATA.temp[,1:(N-1)]
  X_futu = DATA.temp[,2:N]
  X_curr.train = X_curr[,seq(1,N-1,2)]
  X_curr.test = X_curr[,seq(2,N-1,2)]
  X_futu.train = X_futu[,seq(1,N-1,2)]
  X_futu.test = X_futu[,seq(2,N-1,2)]
  init_cpt_train = part2local(DP.VAR1(X_futu.train, X_curr.train, gamma, lambda, delta)$partition)
  init_cpt = 2*init_cpt_train
  len = length(init_cpt)
  init_cpt_long = c(init_cpt_train, ncol(X_curr.train))
  interval = matrix(0, nrow = len+1, ncol = 2)
  interval[1,] = c(1, init_cpt_long[1])
  if(len > 0){
    for(j in 2:(1+len)){
      interval[j,] = c(init_cpt_long[j-1]+1, init_cpt_long[j])
    }
  }
  trainmat = sapply(1:(len+1), function(index) error.pred.seg.VAR1(interval[index,1], interval[index,2], X_futu.train, X_curr.train, lambda, delta))
  transition.list = vector("list", len+1)
  training_loss = matrix(0, nrow = 1, ncol = len+1)
  for(col in 1:(len+1)){
    transition.list[[col]] = trainmat[2,col]$transition.hat
    training_loss[,col] = as.numeric(trainmat[1,col]$MSE)
  }
  validationmat = sapply(1:(len+1), function(index) error.test.VAR1(interval[index,1], interval[index,2], X_futu.test, X_curr.test, transition.list[[index]]))
  result = list(cpt_hat = init_cpt, K_hat = len, test_error = sum(validationmat), train_error = sum(training_loss))
  return(result)
}



glasso.error.test = function(s, e, X.train, Y.train, X.test, Y.test, delta.local, zeta_group_set){
  len = length(zeta_group_set)
  estimate_temp = rep(0, len)
  test_error_temp = rep(0, len)
  for(ll in 1:len){
    estimate_temp[ll] = find.one.change.grouplasso(s, e, X.train, Y.train, delta.local, zeta_group_set[ll])
    test_error_temp[ll] = sum(sapply(1:p, function(m) 
                              test.res.glasso(estimate_temp[ll] - s + 2 , y = Y.train[m,s:e], X= X.train[,s:e], 
                                              y.test = Y.test[m,s:e], X.test = X.test[,s:e], zeta_group_set[ll])))
  }
  return(lambda.group.list[which.min(test.errors.temp)])
}

obj.func.lr.VAR1()

test.res.glasso = function(eta, X.train, y.train, X.test, y.test, zeta.group){
  group = rep(1:p, each=2)
  X.convert = X.glasso.converter.VAR1(X.train, eta, 1)
  out = gglasso(x = X.convert, y = y.train, group=group, loss="ls",
                lambda = zeta.group/ncol(X.train), intercept = FALSE, eps = 0.001)
  coef = as.vector(out$beta)
  res = sum((y.test - X.glasso.converter.VAR1(X.test, eta, 1) %*% coef)^2)
  return(res)
}
