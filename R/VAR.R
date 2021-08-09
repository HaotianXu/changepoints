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
    d = norm(X_futu_hat - X_futu[,s:e], type = "F")
  }else{
    estimate = NA
    d = Inf
  }
  result = list(MSE = d^2, transition.hat = estimate)
  return(result)
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
#' n = 30
#' v1 = 2*(seq(1,p,1)%%2) - 1
#' v2 = -v1
#' AA = matrix(0, nrow = p, ncol = p-2)
#' A1=cbind(v1,v2,AA)
#' A2=cbind(v2,v1,AA)
#' A3=A1
#' data = simu.VAR1(sigma, p, 2*n+1, A1)
#' data = cbind(simu.VAR1(sigma, p, 2*n, A2, vzero=c(data[,ncol(data)])))
#' data = cbind(simu.VAR1(sigma, p, 2*n, A3, vzero=c(data[,ncol(data)])))
#' N = ncol(data)
#' X_curr = data[,1:(N-1)]
#' X_futu = data[,2:N]
#' parti = DP.VAR1(gamma = 1, delta = 5, X_futu, X_curr, lambda = 1)$partition
#' localization = part2local(parti)
DP.VAR1 = function(gamma, delta, X_futu, X_curr, lambda, ...){
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
#' @param y          A \code{numeric} vector of response variable.
#' @param X          A \code{numeric} matrix of covariates.
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
local.refine.VAR1 = function(cpt.init, DATA, zeta.group, w = 1/3){
  N = ncol(DATA)
  p = nrow(DATA)
  X_curr = DATA[,1:(N-1)]
  X_futu = DATA[,2:N]
  cpt.init.ext = c(0, cpt.init, n)
  cpt.init.numb = length(cpt.init)
  cpt.refined = rep(0, cpt.init.numb+1)
  for (k in 1:cpt.init.numb){
    s.inter = w*cpt.refined[k] + (1-w)*cpt.init.ext[k+1]
    e.inter = (1-w)*cpt.init.ext[k+1] + w*cpt.init.ext[k+2]
    lower = ceiling(s.inter) + 1
    upper = floor(e.inter) - 1
    b = sapply(lower:upper, function(eta)obj.func.lr.VAR1(s.inter, e.inter, eta, X_futu, X_curr, zeta.group))
    cpt.refined[k+1] = ceiling(s.inter) + which.min(b)
  }
  return(cpt.refined[-1])
}


#' @title Internal Function: An objective function to select the best splitting location in the local refinement
#' @param s.inter    A \code{numeric} scalar of interpolated starting index.
#' @param e.inter    A \code{numeric} scalar of ending index.
#' @param y          A \code{numeric} vector of response variable.
#' @param X          A \code{numeric} matrix of covariates.
#' @param zeta.group A \code{numeric} scalar of tuning parameter for the group lasso.
#' @noRd
obj.func.lr.VAR1 = function(s.inter, e.inter, eta, X_futu, X_curr, zeta.group){
  n = ncol(X_futu)
  p = nrow(X_futu)
  btemp = rep(NA, p)
  group = rep(1:p, 2)
  X.convert = X.glasso.converter.VAR1(X_curr[,(ceiling(s.inter)):(floor(e.inter))], eta, ceiling(s.inter))
  for(m in 1:p){
    y.convert = X_futu[m, (ceiling(s.inter)):(floor(e.inter))]
    auxfit = gglasso(x = X.convert, y = y.convert, group = group, loss="ls",
                     lambda = zeta.group/(floor(e.inter)-ceiling(s.inter)+1), intercept = FALSE, eps = 0.001)
    coef = as.vector(auxfit$beta)
    coef1 = coef[1:p]
    coef2 = coef[(p+1):(2*p)]
    btemp[m] = norm(y.convert - X.convert %*% coef, type = "2")^2 + zeta.group*sum(sqrt(coef1^2 + coef2^2))
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



#' @title Cross-Validation of Dynamic Programming algorithm for regression change points detection by l0 penalty
#' @description TO DO
#' @param gamma     A \code{numeric} scalar of the tuning parameter associated with the l0 penalty.
#' @param delta     A strictly \code{integer} scalar of minimum spacing.
#' @param y         A \code{numeric} vector of observations.
#' @param X         A \code{numeric} matrix of covariates.
#' @param lambda    A \code{numeric} scalar of tuning parameter for the lasso penalty.
#' @param ...      Additional arguments.
#' @return TO DO.
#' @export
#' @author
#' @examples
#' TO DO
CV.DP.VAR1 = function(DATA, gamma, delta, lambda, ...){
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
  init_cpt_train = part2local(DP.VAR1(gamma, delta, X_futu.train, X_curr.train, lambda)$partition)


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
    transition.list[[col]] = as.numeric(trainmat[2,col]$transition.hat)
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