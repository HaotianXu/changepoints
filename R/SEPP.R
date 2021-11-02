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
#' M = 50 #M is dimension
#' n = 100
#' s = 30 # s is sparsity
#' factor = 0.12 #large factor gives exact recovery
#' threshold = 4 #thresholding makes the process stable
#' intercept = 1/2 #intercept of the model. Assume to be known as in the existing literature
#' A1 = matrix(0,M,M)
#' diag(A1[,-1]) = 1
#' diag(A1) = 1
#' diag(A1[-1,]) = -1
#' A1 = A1*factor
#' A1[(s+1):M,(s+1):M] = 0
#' A2 = matrix(0,M,M)
#' diag(A2[,-1]) = 1
#' diag(A2) = -1
#' diag(A2[-1,]) = 1
#' A2 = A2*factor
#' A2[(s+1):M,(s+1):M] = 0
#' A3 = matrix(0,M,M)
#' diag(A3[,-1]) = 1
#' diag(A3) = 1
#' diag(A3[-1,]) = -1
#' A3 = A3*factor
#' A3[(s+1):M,(s+1):M] = 0
#' data1 = simu.SEPP(intercept, M, A1, threshold, n, vzero = NULL)
#' data2 = simu.SEPP(intercept, M, A2, threshold, n, vzero = data1[,n])
#' data3 = simu.SEPP(intercept, M, A3, threshold, n, vzero = data2[,n])
#' data = cbind(data1, data2, data3)
#' dim(data)
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
#' @param DATA      A \code{numeric} matrix of observations with horizontal axis being time.
#' @param gamma     A \code{numeric} scalar of the tuning parameter associated with the l0 penalty.
#' @param lambda    A \code{numeric} scalar of tuning parameter for lasso penalty.
#' @param delta     An \code{integer} scalar of minimum spacing.
#' @param delta2    An \code{integer} scalar representing the maximal of the change point spacing (for reducing computation cost).
#' @param intercept A \code{numeric} scalar representing the intercept of the model, which is assumed to be known.
#' @param threshold A \code{numeric} scalar representing the upper bound for each coordinate of X_t (for stability).
#' @param ...      Additional arguments.
#' @return A vector of the best partition.
#' @export
#' @author Daren Wang & Haotian Xu
#' @references Wang, D., Yu, Y., & Willett, R. (2020). Detecting Abrupt Changes in High-Dimensional Self-Exciting Poisson Processes. arXiv preprint arXiv:2006.03572.
#' @examples
#' gamma = n/50
#' delta = n/2-1
#' delta2 = 1.5*n
#' intercept = 1/2
#' threshold = 4
#' parti = DP.SEPP(DATA, gamma = gamma, lambda = 500, delta, delta2, intercept, threshold)$partition
#' cpt_hat = part2local(parti)
DP.SEPP = function(DATA, gamma, lambda, delta, delta2, intercept, threshold, ...){
  M = nrow(DATA)
  N = ncol(DATA)
  bestvalue = rep(0,N+1)
  partition = rep(0,N)
  bestvalue[1] = -gamma
  for(r in 1:N){
    bestvalue[r+1] = Inf
    for(l in 1:r){
      b = bestvalue[l] + gamma + error.seg.SEPP(l, r, DATA, lambda, delta, delta2, intercept, threshold)$dist
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
#' @param lambda    A \code{numeric} scalar of lasso penalty.
#' @param delta     A \code{integer} scalar of minimum spacing.
#' @param delta2    A \code{integer} scalar representing the maximal of the change point spacing (for reducing computation cost).
#' @param intercept A \code{numeric} scalar representing the intercept of the model, which is assumed to be known.
#' @param threshold A \code{numeric} scalar representing the upper bound for each coordinate of X_t (for stability).
#' @return  A \code{numeric} scalar of prediction error in terms of poisson log-likelihood function.
#' @noRd
error.seg.SEPP = function(s, e, DATA, lambda, delta, delta2, intercept, threshold){
  print(c(s,e))
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













#' #' @title Internal function: Cross-Validation of Dynamic Programming algorithm for VAR1 change points detection via l0 penalty
#' #' @param DATA      A \code{numeric} matrix of observations with with horizontal axis being time, and vertical axis being dimensions.
#' #' @param gamma     A positive \code{numeric} scalar of the tuning parameter associated with the l0 penalty.
#' #' @param lambda    A positive \code{numeric} scalar of tuning parameter for the lasso penalty.
#' #' @param delta     A positive \code{integer} scalar of minimum spacing.
#' #' @param ...      Additional arguments.
#' #' @noRd
#' CV.DP.SEPP = function(DATA, gamma, lambda, delta, delta2, intercept, threshold, ...){
#'   DATA.temp = DATA
#'   if (ncol(DATA)%%2 != 0){
#'     DATA.temp = DATA[,2:ncol(DATA)]
#'   }
#'   N = ncol(DATA.temp)
#'   p = nrow(DATA.temp)
#'   DATA.train = DATA.temp[,seq(1,N,2)]
#'   DATA.test = DATA.temp[,seq(2,N,2)]
#'   init_cpt_train = part2local(DP.SEPP(DATA.train, gamma, lambda, delta, delta2, intercept, threshold)$partition)
#'   init_cpt = 2*init_cpt_train
#'   len = length(init_cpt)
#'   init_cpt_long = c(init_cpt_train, ncol(DATA.train))
#'   interval = matrix(0, nrow = len+1, ncol = 2)
#'   interval[1,] = c(1, init_cpt_long[1])
#'   if(len > 0){
#'     for(j in 2:(1+len)){
#'       interval[j,] = c(init_cpt_long[j-1]+1, init_cpt_long[j])
#'     }
#'   }
#'   trainmat = sapply(1:(len+1), function(index) error.seg.SEPP(interval[index,1], interval[index,2], DATA.train, lambda, delta, delta2, intercept, threshold))
#'   transition.list = vector("list", len+1)
#'   training_loss = matrix(0, nrow = 1, ncol = len+1)
#'   for(col in 1:(len+1)){
#'     transition.list[[col]] = trainmat[2,col]$transition.hat
#'     training_loss[,col] = as.numeric(trainmat[1,col]$MSE)
#'   }
#'   validationmat = sapply(1:(len+1), function(index) error.test.VAR1(interval[index,1], interval[index,2], X_futu.test, X_curr.test, transition.list[[index]]))
#'   result = list(cpt_hat = init_cpt, K_hat = len, test_error = sum(validationmat), train_error = sum(training_loss))
#'   return(result)
#' }
#' 
#' 
#' #' @title Internal Function: compute testing error for VAR1
#' #' @param  lower     A \code{integer} scalar of starting index.
#' #' @param  upper     A \code{integer} scalar of ending index.
#' #' @param  X_futu    A \code{numeric} matrix of observations.
#' #' @param  X_curr    A \code{numeric} matrix of covariates.
#' #' @param  transition.hat A \code{numeric} matrix of transition matrix estimator.
#' #' @return A numeric scalar of testing error in squared Frobenius norm.
#' #' @noRd
#' error.test.SEPP = function(lower, upper, X_futu, X_curr, transition.hat){
#'   res = norm(X_futu[lower:upper] - transition.hat%*%X_curr[,lower:upper], type = "F")^2
#'   return(res)
#' }
#' 
#' 
#' #' @title Grid search based on Cross-Validation of Dynamic Programming for regression change points detection via l0 penalty
#' #' @description TO DO
#' #' @param DATA          A \code{numeric} matrix of observations with with horizontal axis being time, and vertical axis being dimensions.
#' #' @param gamma.set     A \code{numeric} vector of candidate tuning parameter associated with the l0 penalty.
#' #' @param lambda.set    A \code{numeric} vector of candidate tuning parameter for the lasso penalty.
#' #' @param delta         A strictly \code{integer} scalar of minimum spacing.
#' #' @param ...           Additional arguments.
#' #' @return Row: lambda.set; column: gamma.set
#' #' @export
#' #' @author
#' #' @examples
#' #' TO DO
#' CV.search.DP.VAR1 = function(DATA, gamma.set, lambda.set, delta, ...){
#'   output = sapply(1:length(lambda.set), function(i) sapply(1:length(gamma.set), 
#'                                                            function(j) CV.DP.VAR1(DATA, gamma.set[j], lambda.set[i], delta)))
#'   cpt_hat = output[seq(1,4*length(gamma.set),4),]## estimated change points
#'   K_hat = output[seq(2,4*length(gamma.set),4),]## number of estimated change points
#'   test_error = output[seq(3,4*length(gamma.set),4),]## validation loss
#'   train_error = output[seq(4,4*length(gamma.set),4),]## training loss                                                      
#'   result = list(cpt_hat = cpt_hat, K_hat = K_hat, test_error = test_error, train_error = train_error)
#'   return(result)
#' }
#' 
