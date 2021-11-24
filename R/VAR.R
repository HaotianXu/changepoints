#' @title Simulate from a VAR1 model (without change point).
#' @description Simulate data of size n and dimension p from a VAR1 model (without change point) with Gaussian i.i.d. error terms.
#' @param sigma      A \code{numeric} scalar representing the standard deviation of error terms.
#' @param p          An \code{integer} scalar representing dimension.
#' @param n          An \code{integer} scalar representing sample size.
#' @param A          A \code{numeric} p-by-p matrix representing the transition matrix of the VAR1 model.
#' @param vzero      A \code{numeric} vector representing the observation at time 0. If \code{vzero = NULL},it is generated following the distribution of the error terms. 
#' @param ...        Additional arguments.
#' @return  A p-by-n matrix.
#' @export
#' @author  Daren Wang
#' @examples
#' p = 20
#' sigma = 1
#' n = 100
#' A = matrix(rnorm(p*p), nrow = p)*0.1  # transition matrix
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


#' @title Internal Function: Prediction error in squared Frobenius norm for the lasso estimator of transition matrix.
#' @param X_futu    A \code{numeric} matrix of time series at one step ahead.
#' @param X_curr    A \code{numeric} matrix of time series at current step.
#' @param s         A \code{integer} scalar of starting index.
#' @param e         A \code{integer} scalar of ending index.
#' @param lambda    A \code{numeric} scalar of lasso penalty.
#' @param delta     A \code{integer} scalar of minimum spacing.
#' @param eps       A \code{numeric} scalar of precision level for convergence.
#' @return    A \code{list} with the following structure:
#'  \item{MSE}{A \code{numeric} scalar of prediction error in Frobenius norm}
#'  \item{tran_hat}{A p-by-p matrix of estimated transition matrix}
#' @return  A \code{numeric} scalar of prediction error in Frobenius norm.
#' @noRd
error.pred.seg.VAR1 <- function(X_futu, X_curr, s, e, lambda, delta, eps) {
  .Call('_changepoints_rcpp_error_pred_seg_VAR1', PACKAGE = 'changepoints', X_futu, X_curr, s, e, lambda, delta, eps)
}


#' @title Dynamic programming for VAR1 change points detection through \eqn{l_0} penalty.
#' @description Perform dynamic programming for VAR1 change points detection through \eqn{l_0} penalty.
#' @param X_futu    A \code{numeric} matrix of time series at one step ahead, with horizontal axis being time.
#' @param X_curr    A \code{numeric} matrix of time series at current step, with horizontal axis being time.
#' @param gamma     A \code{numeric} scalar of the tuning parameter associated with \eqn{l_0} penalty.
#' @param lambda    A \code{numeric} scalar of the tuning parameter for lasso penalty.
#' @param delta     A strictly \code{integer} scalar of minimum spacing.
#' @param eps       A \code{numeric} scalar of precision level for convergence of lasso.
#' @param ...       Additional arguments.
#' @return  partition: A vector of the best partition.
#' @export
#' @author Daren Wang, Haotian Xu
#' @examples
#' p = 20
#' sigma = 1
#' n = 20
#' v1 = 2*(seq(1,p,1)%%2) - 1
#' v2 = -v1
#' AA = matrix(0, nrow = p, ncol = p-2)
#' A1 = cbind(v1,v2,AA)*0.1
#' A2 = cbind(v2,v1,AA)*0.1
#' A3 = A1
#' data = simu.VAR1(sigma, p, 2*n+1, A1)
#' data = cbind(data, simu.VAR1(sigma, p, 2*n, A2, vzero=c(data[,ncol(data)])))
#' data = cbind(data, simu.VAR1(sigma, p, 2*n, A3, vzero=c(data[,ncol(data)])))
#' N = ncol(data)
#' X_curr = data[,1:(N-1)]
#' X_futu = data[,2:N]
#' parti = DP.VAR1(X_futu, X_curr, gamma = 1, lambda = 1, delta = 5)$partition
#' part2local(parti)
DP.VAR1 <- function(X_futu, X_curr, gamma, lambda, delta, eps = 0.001, ...) {
  .Call('_changepoints_rcpp_DP_VAR1', PACKAGE = 'changepoints', X_futu, X_curr, gamma, lambda, delta, eps)
}


#' @title Internal function: Cross-Validation of Dynamic Programming algorithm for VAR1 change points detection via \eqn{l_0} penalty
#' @param DATA      A \code{numeric} matrix of observations with with horizontal axis being time, and vertical axis being dimensions.
#' @param gamma     A positive \code{numeric} scalar of the tuning parameter associated with the \eqn{l_0} penalty.
#' @param lambda    A positive \code{numeric} scalar of the tuning parameter for lasso penalty.
#' @param delta     A positive \code{integer} scalar of minimum spacing.
#' @param eps       A \code{numeric} scalar of precision level for convergence of lasso.
#' @noRd
CV.DP.VAR1 = function(DATA, gamma, lambda, delta, eps = 0.001){
  DATA.temp = DATA
  if (ncol(DATA)%%2 != 0){
    DATA.temp = DATA[,2:ncol(DATA)]
  }
  N = ncol(DATA.temp)
  p = nrow(DATA.temp)
  X_train = DATA.temp[,seq(1,N,2)]
  X_test = DATA.temp[,seq(2,N,2)]
  X_curr.train = X_train[,1:(N/2-1)]
  X_futu.train = X_train[,2:(N/2)]
  X_curr.test = X_test[,1:(N/2-1)]
  X_futu.test = X_test[,2:(N/2)]
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
  trainmat = sapply(1:(len+1), function(index) error.pred.seg.VAR1(X_futu.train, X_curr.train, interval[index,1], interval[index,2], lambda, delta, eps))
  transition.list = vector("list", len+1)
  training_loss = matrix(0, nrow = 1, ncol = len+1)
  for(col in 1:(len+1)){
    transition.list[[col]] = trainmat[2,col]$tran_hat
    training_loss[,col] = as.numeric(trainmat[1,col]$MSE)
  }
  validationmat = sapply(1:(len+1), function(index) error.test.VAR1(X_futu.test, X_curr.test, interval[index,1], interval[index,2], transition.list[[index]]))
  result = list(cpt_hat = init_cpt, K_hat = len, test_error = sum(validationmat), train_error = sum(training_loss))
  return(result)
}


#' @title Internal Function: compute testing error for VAR1
#' @param  X_futu    A \code{numeric} matrix of time series at one step ahead, with horizontal axis being time.
#' @param  X_curr    A \code{numeric} matrix of time series at current step, with horizontal axis being time.
#' @param  lower     A \code{integer} scalar of starting index.
#' @param  upper     A \code{integer} scalar of ending index.
#' @param  tran_hat  A \code{numeric} matrix of transition matrix estimator.
#' @return A numeric scalar of testing error in squared Frobenius norm.
#' @noRd
error.test.VAR1 = function(X_futu, X_curr, lower, upper, tran_hat){
  if(sum(dim(tran_hat) == c(0,0))){
    res = Inf
  }else{
    res = norm(X_futu[,lower:upper] - tran_hat%*%X_curr[,lower:upper], type = "F")^2
  }
  return(res)
}


#' @title Grid search based on cross-validation of dynamic programming for VAR change points detection via \eqn{l_0} penalty.
#' @description  Perform grid search based on cross-validation of dynamic programming for VAR change points detection.
#' @param DATA          A \code{numeric} matrix of observations with horizontal axis being time, and vertical axis being dimensions.
#' @param gamma_set     A \code{numeric} vector of candidate tuning parameters associated with the \eqn{l_0} penalty.
#' @param lambda_set    A \code{numeric} vector of candidate tuning parameters for lasso penalty.
#' @param delta         A strictly \code{integer} scalar of minimum spacing.
#' @param eps           A \code{numeric} scalar of precision level for convergence of lasso.
#' @param ...           Additional arguments.
#' @return  A \code{list} with the following structure:
#'  \item{cpt_hat}{A list of vector of estimated change points}
#'  \item{K_hat}{A list of scalar of number of estimated change points}
#'  \item{test_error}{A list of vector of testing errors (each row corresponding to each gamma, and each column corresponding to each lambda)}
#'  \item{train_error}{A list of vector of training errors}
#' @export
#' @author   Daren Wang & Haotian Xu
#' @examples
#' set.seed(123)
#' p = 20
#' sigma = 1
#' n = 20
#' v1 = 2*(seq(1,p,1)%%2) - 1
#' v2 = -v1
#' AA = matrix(0, nrow = p, ncol = p-2)
#' A1 = cbind(v1,v2,AA)*0.1
#' A2 = cbind(v2,v1,AA)*0.1
#' A3 = A1
#' cpt_true = c(40, 80)
#' data = simu.VAR1(sigma, p, 2*n+1, A1)
#' data = cbind(data, simu.VAR1(sigma, p, 2*n, A2, vzero=c(data[,ncol(data)])))
#' data = cbind(data, simu.VAR1(sigma, p, 2*n, A3, vzero=c(data[,ncol(data)])))
#' gamma_set = c(0.1, 0.5, 1)
#' lambda_set = c(0.1, 1, 3.2)
#' temp = CV.search.DP.VAR1(data, gamma_set, lambda_set, delta = 5)
#' temp$test_error # test error result
#' # find the indices of gamma.set and lambda.set which minimizes the test error
#' min_idx = as.vector(arrayInd(which.min(temp$test_error), dim(temp$test_error)))
#' cpt_init = unlist(temp$cpt_hat[min_idx[1], min_idx[2]])
#' Hausdorff.dist(cpt_init, cpt_true)
CV.search.DP.VAR1 = function(DATA, gamma_set, lambda_set, delta, eps = 0.001, ...){
  output = sapply(1:length(lambda_set), function(i) sapply(1:length(gamma_set), 
                                                           function(j) CV.DP.VAR1(DATA, gamma_set[j], lambda_set[i], delta, eps)))
  cpt_hat = output[seq(1,4*length(gamma_set),4),]## estimated change points
  K_hat = output[seq(2,4*length(gamma_set),4),]## number of estimated change points
  test_error = output[seq(3,4*length(gamma_set),4),]## validation loss
  train_error = output[seq(4,4*length(gamma_set),4),]## training loss                                                      
  result = list(cpt_hat = cpt_hat, K_hat = K_hat, test_error = test_error, train_error = train_error)
  return(result)
}




#' @title Local refinement for VAR1 change points detection.
#' @description Perform local refinement for VAR1 change points detection.
#' @param cpt_init   A \code{integer} vector of initial change points estimation (sorted in strictly increasing order).
#' @param DATA       A \code{numeric} matrix of observations with with horizontal axis being time, and vertical axis being dimensions.
#' @param zeta       A \code{numeric} scalar of lasso penalty.
#' @param ...       Additional arguments.
#' @return  An \code{integer} vector of locally refined change points estimation.
#' @export
#' @author  Daren Wang & Haotian Xu
#' @examples
#' set.seed(123)
#' p = 20
#' sigma = 1
#' n = 20
#' v1 = 2*(seq(1,p,1)%%2) - 1
#' v2 = -v1
#' AA = matrix(0, nrow = p, ncol = p-2)
#' A1 = cbind(v1,v2,AA)*0.1
#' A2 = cbind(v2,v1,AA)*0.1
#' A3 = A1
#' data = simu.VAR1(sigma, p, 2*n+1, A1)
#' data = cbind(data, simu.VAR1(sigma, p, 2*n, A2, vzero=c(data[,ncol(data)])))
#' data = cbind(data, simu.VAR1(sigma, p, 2*n, A3, vzero=c(data[,ncol(data)])))
#' gamma_set = c(0.1, 0.5, 1)
#' lambda_set = c(0.1, 1, 3.2)
#' temp = CV.search.DP.VAR1(data, gamma_set, lambda_set, delta = 5)
#' temp$test_error # test error result
#' # find the indices of gamma.set and lambda.set which minimizes the test error
#' min_idx = as.vector(arrayInd(which.min(temp$test_error), dim(temp$test_error))) 
#' cpt_init = unlist(temp$cpt_hat[min_idx[1], min_idx[2]])
#' local.refine.VAR1(cpt_init, data, zeta = 0.5)
local.refine.VAR1 = function(cpt_init, DATA, zeta, ...){
  w = 0.9
  N = ncol(DATA)
  p = nrow(DATA)
  X_curr = DATA[,1:(N-1)]
  X_futu = DATA[,2:N]
  cpt_init_ext = c(0, cpt_init, N)
  cpt_init_numb = length(cpt_init)
  cpt_refined = rep(0, cpt_init_numb+1)
  for (k in 1:cpt_init_numb){
    s_inter = w*cpt_init_ext[k] + (1-w)*cpt_init_ext[k+1]
    e_inter = (1-w)*cpt_init_ext[k+1] + w*cpt_init_ext[k+2]
    lower = ceiling(s_inter) + 1
    upper = floor(e_inter) - 1
    b = sapply(lower:upper, function(eta)obj.LR.VAR1(ceiling(s_inter), floor(e_inter), eta, X_futu, X_curr, zeta))
    cpt_refined[k+1] = which.min(b) + lower - 1
  }
  return(cpt_refined[-1])
}


#' @title Internal Function: An objective function to select the best splitting location in the local refinement
#' @param s.inter    A \code{numeric} scalar of interpolated starting index.
#' @param e.inter    A \code{numeric} scalar of ending index.
#' @param eta        A \code{integer} scalar of splitting index.
#' @param X_futu     A \code{numeric} matrix of time series at one step ahead.
#' @param X_curr     A \code{numeric} matrix of time series at current step.
#' @param zeta.group A \code{numeric} scalar of tuning parameter for the group lasso.
#' @noRd
obj.LR.VAR1 = function(s.inter.ceil, e.inter.floor, eta, X_futu, X_curr, zeta.group){
  n = ncol(X_futu)
  p = nrow(X_futu)
  len.inter = e.inter.floor-s.inter.ceil+1
  btemp = rep(NA, p)
  group = rep(1:p, 2)
  X.convert = X.glasso.converter.VAR1(X_curr[,s.inter.ceil:e.inter.floor], eta, s.inter.ceil)
  Y.convert = X_futu[, s.inter.ceil:e.inter.floor]
  for(m in 1:p){
    y.convert = Y.convert[m,]
    auxfit = gglasso(x = X.convert, y = y.convert, group = group, intercept = FALSE, loss="ls",
             lambda = zeta.group/len.inter)
    coef = as.vector(auxfit$beta)
    btemp[m] = sum((y.convert - X.convert %*% coef)^2) + zeta.group*sqrt(sum(coef^2))
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
  p = nrow(X)
  xx1 = xx2 = t(X)
  t = eta - s_ceil + 1
  xx1[(t+1):n,] = 0
  xx2[1:t,] = 0
  xx = cbind(xx1/sqrt(t-1), xx2/sqrt(n-t))
  index=c()
  for(pp in 1:p){
    index = c(index, pp, pp+p)
  }
  xxx = xx[,index]
  return(xxx)
}


#' @title Local refinement for VAR1 change points detection.
#' @description Perform local refinement for VAR1 change points detection.
#' @param cpt_init        A \code{integer} vector of initial change points estimation (sorted in strictly increasing order).
#' @param DATA            A \code{numeric} matrix of observations with with horizontal axis being time, and vertical axis being dimensions.
#' @param zeta_set        A \code{numeric} vector of candidate tuning parameters for group lasso penalty.
#' @param delta_local     A strictly \code{integer} scalar of minimum spacing for group lasso.
#' @param ...             Additional arguments.
#' @return  A \code{list} with the following structure:
#'  \item{cpt_hat}{A vector of estimated change point locations (sorted in strictly increasing order)}
#'  \item{zeta}{A scalar of selected zeta by cross-validation}
#' @export
#' @author  Daren Wang & Haotian Xu
#' @examples
#' set.seed(123)
#' p = 20
#' sigma = 1
#' n = 20
#' v1 = 2*(seq(1,p,1)%%2) - 1
#' v2 = -v1
#' AA = matrix(0, nrow = p, ncol = p-2)
#' A1 = cbind(v1,v2,AA)*0.1
#' A2 = cbind(v2,v1,AA)*0.1
#' A3 = A1
#' data = simu.VAR1(sigma, p, 2*n+1, A1)
#' data = cbind(data, simu.VAR1(sigma, p, 2*n, A2, vzero=c(data[,ncol(data)])))
#' data = cbind(data, simu.VAR1(sigma, p, 2*n, A3, vzero=c(data[,ncol(data)])))
#' gamma_set = c(0.1, 0.5, 1)
#' lambda_set = c(0.1, 1, 3.2)
#' temp = CV.search.DP.VAR1(data, gamma_set, lambda_set, delta = 5)
#' temp$test_error # test error result
#' # find the indices of gamma.set and lambda.set which minimizes the test error
#' min_idx = as.vector(arrayInd(which.min(temp$test_error), dim(temp$test_error))) 
#' cpt_init = unlist(temp$cpt_hat[min_idx[1], min_idx[2]])
#' temp_zeta = local.refine.CV.VAR1(cpt_init, data, c(0.1, 0.5, 1), delta_local = 10)
#' temp_zeta$cpt_hat
local.refine.CV.VAR1 = function(cpt_init, DATA, zeta_set, delta_local, ...){
  w = 0.9
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
  cpt_hat_train_ext = c(1, round(cpt_init/2), floor(ncol(X_curr)/2))
  KK = length(cpt_init)
  if(KK > 0){
    for(kk in 1:KK){
      # find the zeta which gives the smallest testing error 
      zeta_temp = glasso.error.test(s = cpt_hat_train_ext[kk], e = cpt_hat_train_ext[kk+2], X_futu.train, X_curr.train, X_futu.test, X_curr.test, delta_local, zeta_set)
      s.inter = w*cpt_hat_train_ext[kk] + (1-w)*cpt_hat_train_ext[kk+1]
      e.inter = (1-w)*cpt_hat_train_ext[kk+1] + w*cpt_hat_train_ext[kk+2]
      temp_estimate = find.one.change.grouplasso.VAR1(s = ceiling(2*s.inter), e = floor(2*e.inter), X_futu, X_curr, delta_local, zeta.group = zeta_temp)
      cpt_hat_train_ext[kk+1] = round(temp_estimate/2)
    }
  }
  return(list(cpt_hat = 2*cpt_hat_train_ext[-c(1, length(cpt_hat_train_ext))], zeta = zeta_temp))
}


# return best zeta
#' @noRd
glasso.error.test = function(s, e, Y.train, X.train, Y.test, X.test, delta.local, zeta_group_set){
  p = nrow(Y.train)
  len = length(zeta_group_set)
  estimate_temp = rep(0, len)
  test_error_temp = rep(0, len)
  for(ll in 1:len){
    estimate_temp[ll] = find.one.change.grouplasso.VAR1(s, e, Y.train, X.train, delta.local, zeta_group_set[ll])
    test_error_temp[ll] = sum(sapply(1:p, function(m)
                              test.res.glasso(p, estimate_temp[ll] - s + 2 , Y.train[m,s:e], X.train[,s:e],
                                              Y.test[m,s:e], X.test[,s:e], zeta_group_set[ll])))
  }
  return(zeta_group_set[which.min(test_error_temp)])
}

# return the best split location
#' @noRd
find.one.change.grouplasso.VAR1 = function(s, e, Y.train, X.train, delta.local, zeta.group){
  estimate = (s+e)/2
  if(e-s > 2*delta.local){
    can.vec = c((s+delta.local):(e-delta.local))
    #can.vec = can.vec[which(can.vec%%2==0)]
    res.seq = sapply(can.vec, function(t) 
                     obj.LR.VAR1(s, e, t, Y.train, X.train, zeta.group)) 
    estimate = can.vec[which.min(res.seq)]
  }
  return(estimate)
}

#' @noRd
test.res.glasso = function(p, eta, y.train, X.train, y.test, X.test, zeta.group){
  group = rep(1:p, each=2)
  X.convert = X.glasso.converter.VAR1(X.train, eta, 1)
  X.test.convert = X.glasso.converter.VAR1(X.test, eta, 1)
  auxfit = gglasso(x = X.convert, y = y.train, group = group, intercept = FALSE, loss="ls",
                   lambda = zeta.group/ncol(X.train))
  coef = as.vector(auxfit$beta)
  res = sum((y.test - X.test.convert %*% coef)^2) + zeta.group*sqrt(sum(coef^2))
  return(res)
}
