#' @title Simulate the design matrix following AR1
#' @noRd
simu_ar_data <- function(ar1_rho, n, Cov_X, burnin=50){
  p <- dim(Cov_X)[1]
  X_data <- matrix(0,n+burnin,p)
  for(time_index in 2:(n+burnin)){
    # print(time_index)
    X_data[time_index,] <- ar1_rho*X_data[time_index-1,]+mvrnorm(n=1,mu=rep(0,p),Sigma=Cov_X)
  }
  return(X_data[-c(1:burnin),])
}

#' @title Simulate the design matrix following MA1
#' @noRd
simu_ma_data <- function(ma1_theta, n, Cov_X){
  p <- dim(Cov_X)[1]
  error_data0 <- mvrnorm(n=1,mu=rep(0,p),Sigma=Cov_X)
  error_data1 <- mvrnorm(n=1,mu=rep(0,p),Sigma=Cov_X)
  X_data <- c()
  for(time_index in 1:n){
    X_data <- rbind(X_data, error_data1+ma1_theta*error_data0)
    error_data0 <- error_data1
    error_data1 <- mvrnorm(n=1,mu=rep(0,p),Sigma=Cov_X)
  }
  return(X_data)
}


#' @title Simulate a sparse regression model with change points in coefficients.
#' @description      Simulate a sparse regression model with change points in coefficients under temporal dependence.
#' @param d0         A \code{numeric} scalar stands for the number of nonzero coefficients.
#' @param cpt_true   An \code{integer} vector contains true change points (sorted in strictly increasing order).
#' @param p          An \code{integer} scalar stands for the dimensionality.
#' @param n          An \code{integer} scalar stands for the sample size.
#' @param sigma      A \code{numeric} scalar stands for error standard deviation.
#' @param kappa      A \code{numeric} scalar stands for the minimum jump size of coefficient vector in \eqn{l_2} norm.
#' @param cov_type   A \code{character} string stands for the type of covariance matrix of covariates. 'I': Identity; 'T': Toeplitz; 'E': Equal-correlation.
#' @param mod_X      A \code{character} string stands for the time series model followed by the covariates. 'IID': IID multivariate Gaussian; 'AR': Multivariate AR1 with rho = 0.3; Multivariate MA1 theta = 0.3.
#' @param mod_e      A \code{character} string stands for the time series model followed by the errors 'IID': IID univariate Gaussian; 'AR': Univariate AR1 with rho = 0.3; Univariate MA1 theta = 0.3.
#' @return A \code{list} with the following structure:
#'  \item{cpt_true}{A vector of true changepoints (sorted in strictly increasing order).}
#'  \item{X}{An n-by-p design matrix.}
#'  \item{y}{An n-dim vector of response variable.}
#'  \item{betafullmat}{A p-by-n matrix of coefficients.}
#' @export
#' @author Daren Wang, Zifeng Zhao & Haotian Xu
#' @references Rinaldo, Wang, Wen, Willett and Yu (2020) <arxiv:2010.10410>; Xu, Wang, Zhao and Yu (2022) <arXiv:2207.12453>.
#' @examples
#' d0 = 10
#' p = 30
#' n = 100
#' cpt_true = c(10, 30, 40, 70, 90)
#' data = simu.change.regression(d0, cpt_true, p, n, sigma = 1, kappa = 9)
simu.change.regression = function(d0, cpt_true, p, n, sigma, kappa, cov_type = 'I', mod_X = 'IID', mod_e = 'IID'){
  if(d0 >= p){
    stop("d0 should be strictly smaller than p")
  }
  if(sigma <= 0){
    stop("sigma should be strictly larger than 0")
  }
  if(kappa <= 0){
    stop("kappa should be strictly larger than 0")
  }
  if(!is.element(cov_type, c('I', 'T', 'E'))){
    stop("cov_type should be one element of c('I', 'T', 'E')")
  }
  if(!is.element(mod_X, c('IID', 'AR', 'MA'))){
    stop("mod_X should be one element of c('IID', 'AR', 'MA')")
  }
  if(!is.element(mod_e, c('IID', 'AR', 'MA'))){
    stop("mod_e should be one element of c('IID', 'AR', 'MA')")
  }
  no.cpt = length(cpt_true)
  if(is.unsorted(cpt_true, strictly = TRUE) | min(cpt_true) <= 1 | max(cpt_true >= n) | no.cpt > n-2){
    stop("cpt_true is not correctly specified")
  }
  ### covariance matrix of X
  if(cov_type == 'I'){
    cov_X <- diag(p)
  }else if(cov_type == 'T'){
    cov_X <- matrix(NA, nrow = p, ncol = p)
    for(i in 1:p){
      for(j in 1:p){
        cov_X[i,j] = 0.5^(abs(i-j))
      }
    }
  }else if(cov_type == 'E'){
    cov_X <- matrix(0.3,p,p)
    diag(cov_X) <- 1
  }
  cov_X <- cov_X
  ### time series model of X
  if(mod_X == "IID"){
    X = mvrnorm(n = n, mu = rep(0,p), Sigma = cov_X)
  }else if(mod_X == "AR"){
    X = simu_ar_data(ar1_rho = 0.3, n = n, cov_X)
  }else if(mod_X == "MA"){
    X = simu_ma_data(ma1_theta = 0.3, n = n, cov_X)
  }
  ### time series model of errors
  if(mod_e == "IID"){
    err = rnorm(n, mean = 0, sd = sigma)
  }else if(mod_X == "AR"){
    err = as.numeric(arima.sim(list(order=c(1,0,0), ar=0.3), n = n, sd = sqrt(0.91)*sigma))
  }else if(mod_X == "MA"){
    err = as.numeric(arima.sim(list(order=c(0,0,1), ma=0.3), n = n, sd = sigma/sqrt(1.09)))
  }
  y = matrix(0, n, 1)
  nonzero.element.loc = c(1:d0)
  cpt = c(0, cpt_true, n)
  beta = matrix(0, p, no.cpt+1)
  betafullmat = matrix(0, p, n)
  for (i in 1:(no.cpt+1)){
    if (i%%2 == 1){
      beta[nonzero.element.loc, i] = kappa/(2*sqrt(d0))
    }
    else{
      beta[nonzero.element.loc, i] = -kappa/(2*sqrt(d0))
    }
    y[(1+cpt[i]):cpt[i+1],] = X[(1+cpt[i]):cpt[i+1],] %*% beta[,i] + err[(1+cpt[i]):cpt[i+1]]
    for (j in (1+cpt[i]):cpt[i+1]){
      betafullmat[,j] = beta[,i] 
    }
  }
  List = list(cpt_true = cpt_true, X = X, y = as.vector(y), betafullmat = betafullmat)
  return(List)
}



#' @title Dynamic programming algorithm for regression change points localisation with \eqn{l_0} penalisation.
#' @description     Perform dynamic programming algorithm for regression change points localisation.
#' @param y         A \code{numeric} vector of response variable.
#' @param X         A \code{numeric} matrix of covariates with vertical axis being time.
#' @param gamma     A positive \code{numeric} scalar stands for tuning parameter associated with \eqn{l_0} penalty.
#' @param lambda    A positive \code{numeric} scalar stands for tuning parameter associated with the lasso penalty.
#' @param delta     A positive \code{integer} scalar stands for minimum spacing.
#' @param eps       A \code{numeric} scalar of precision level for convergence of lasso.
#' @return An object of \code{\link[base]{class}} "DP", which is a \code{list} with the following structure:
#'  \item{partition}{A vector of the best partition.}
#'  \item{cpt}{A vector of change points estimation.}
#' @export
#' @author Daren Wang & Haotian Xu
#' @references Rinaldo, Wang, Wen, Willett and Yu (2020) <arxiv:2010.10410>
#' @examples
#' d0 = 10
#' p = 20
#' n = 100
#' cpt_true = c(30, 70)
#' data = simu.change.regression(d0, cpt_true, p, n, sigma = 1, kappa = 9)
#' temp = DP.regression(y = data$y, X = data$X, gamma = 2, lambda = 1, delta = 5)
#' cpt_hat = temp$cpt
#' @export
DP.regression <- function(y, X, gamma, lambda, delta, eps = 0.001) {
  DP_result = .Call('_changepoints_rcpp_DP_regression', PACKAGE = 'changepoints', y, X, gamma, lambda, delta, eps)
  result = append(DP_result, list(cpt = part2local(DP_result$partition)))
  class(result) = "DP"
  return(result)
}




#' @title Internal Function: Prediction error in squared \eqn{l_2} norm for the lasso.
#' @param y         A \code{numeric} vector of response variable.
#' @param X         A \code{numeric} matrix of covariates with vertical axis being time.
#' @param s         An \code{integer} scalar of starting index.
#' @param e         An \code{integer} scalar of ending index.
#' @param lambda    A \code{numeric} scalar of tuning parameter for lasso penalty.
#' @param delta     A \code{integer} scalar of minimum spacing.
#' @param eps       A \code{numeric} scalar of precision level for convergence of lasso.
#' @return    A \code{list} with the following structure:
#'  \item{MSE}{A \code{numeric} scalar of prediction error in \eqn{l_2} norm}
#'  \item{beta_hat}{A p-dim vector of estimated coefficients}
#' @noRd
error.pred.seg.regression <- function(y, X, s, e, lambda, delta, eps = 0.001) {
  .Call('_changepoints_rcpp_error_pred_seg_regression', PACKAGE = 'changepoints', y, X, s, e, lambda, delta, eps)
}



#' @title Internal function: Cross-validation of dynamic programming algorithm for regression change points localisation through \eqn{l_0} penalisation.
#' @description     Perform cross-validation of dynamic programming algorithm for regression change points.
#' @param y         A \code{numeric} vector of response variable.
#' @param X         A \code{numeric} matrix of covariates with vertical axis being time.
#' @param gamma     A \code{numeric} scalar of tuning parameter associated with the \eqn{l_0} penalty.
#' @param lambda    A \code{numeric} scalar of tuning parameter for the lasso penalty.
#' @param delta     A strictly \code{integer} scalar of minimum spacing.
#' @param eps       A \code{numeric} scalar of precision level for convergence of lasso.
#' @return  A \code{list} with the following structure:
#'  \item{cpt_hat}{A vector of estimated change points locations (sorted in strictly increasing order)}
#'  \item{K_hat}{A scalar of number of estimated change points}
#'  \item{test_error}{A list of vector of testing errors in squared \eqn{l_2} norm}
#'  \item{train_error}{A list of vector of training errors in squared \eqn{l_2} norm}
#' @noRd
CV.DP.regression = function(y, X, gamma, lambda, delta, eps = 0.001){
  N = nrow(X)
  even_indexes = seq(2, N, 2)
  odd_indexes = seq(1, N, 2)
  train.X = X[odd_indexes,]
  train.y = y[odd_indexes]
  validation.X = X[even_indexes,]
  validation.y = y[even_indexes]
  init_cpt_train = DP.regression(train.y, train.X, gamma, lambda, delta, eps)$cpt
  init_cpt_train.long = c(0, init_cpt_train, nrow(train.X))
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
    p = ncol(train.X)
    trainmat = sapply(1:(len+1), function(index) error.pred.seg.regression(train.y, train.X, interval[index,1], interval[index,2], lambda, delta, eps))
    betamat = matrix(0, nrow = p, ncol = len+1)
    training_loss = matrix(0, nrow = 1, ncol = len+1)
    for(col in 1:(len+1)){
      betamat[,col] = as.numeric(trainmat[2,col]$beta_hat)
      training_loss[,col] = as.numeric(trainmat[1,col]$MSE)
    }
    validationmat = sapply(1:(len+1), function(index) error.test.regression(validation.y, validation.X, interval[index,1], interval[index,2], betamat[,index]))
    result = list(cpt_hat = init_cpt, K_hat = len, test_error = sum(validationmat), train_error = sum(training_loss))
  }
  return(result)
}


#' @title Internal Function: compute testing error for regression.
#' @param lower     A \code{integer} scalar of starting index.
#' @param upper     A \code{integer} scalar of ending index.
#' @param y         A \code{numeric} vector of response variable.
#' @param X         A \code{numeric} matrix of covariates with vertical axis being time.
#' @return A numeric scalar of testing error in squared l2 norm.
#' @noRd
error.test.regression = function(y, X, lower, upper, beta_hat){
  res = norm(y[lower:upper] - X[lower:upper,]%*%beta_hat, type = "2")^2
  return(res)
} 


#' @title Grid search based on cross-validation of dynamic programming for regression change points localisation with \eqn{l_0} penalisation.
#' @description Perform grid search to select tuning parameters gamma (for \eqn{l_0} penalty of DP) and lambda (for lasso penalty) based on cross-validation.
#' @param y             A \code{numeric} vector of response variable.
#' @param X             A \code{numeric} matrix of covariates with vertical axis being time.
#' @param gamma_set     A \code{numeric} vector of candidate tuning parameters associated with \eqn{l_0} penalty of DP.
#' @param lambda_set    A \code{numeric} vector of candidate tuning parameters for lasso penalty.
#' @param delta         A strictly positive \code{integer} scalar of minimum spacing.
#' @param eps           A \code{numeric} scalar of precision level for convergence of lasso.
#' @return  A \code{list} with the following structure:
#'  \item{cpt_hat}{A list of vector of estimated change points}
#'  \item{K_hat}{A list of scalar of number of estimated change points}
#'  \item{test_error}{A list of vector of testing errors (each row corresponding to each gamma, and each column corresponding to each lambda)}
#'  \item{train_error}{A list of vector of training errors}
#' @export
#' @author Daren Wang
#' @examples
#' d0 = 10
#' p = 20
#' n = 100
#' cpt_true = c(30, 70)
#' data = simu.change.regression(d0, cpt_true, p, n, sigma = 1, kappa = 9)
#' gamma_set = c(0.01, 0.1, 1)
#' lambda_set = c(0.01, 0.1, 1, 3)
#' temp = CV.search.DP.regression(y = data$y, X = data$X, gamma_set, lambda_set, delta = 2)
#' temp$test_error # test error result
#' # find the indices of gamma_set and lambda_set which minimizes the test error
#' min_idx = as.vector(arrayInd(which.min(temp$test_error), dim(temp$test_error))) 
#' gamma_set[min_idx[1]]
#' lambda_set[min_idx[2]]
#' cpt_init = unlist(temp$cpt_hat[min_idx[1], min_idx[2]])
#' @references Rinaldo, Wang, Wen, Willett and Yu (2020) <arxiv:2010.10410>
CV.search.DP.regression = function(y, X, gamma_set, lambda_set, delta, eps = 0.001){
  output = sapply(1:length(lambda_set), function(i) sapply(1:length(gamma_set), 
                                                           function(j) CV.DP.regression(y, X, gamma_set[j], lambda_set[i], delta)))
  cpt_hat = output[seq(1,4*length(gamma_set),4),]## estimated change points
  K_hat = output[seq(2,4*length(gamma_set),4),]## number of estimated change points
  test_error = output[seq(3,4*length(gamma_set),4),]## validation loss
  train_error = output[seq(4,4*length(gamma_set),4),]## training loss                                                      
  result = list(cpt_hat = cpt_hat, K_hat = K_hat, test_error = test_error, train_error = train_error)
  return(result)
}





#' @title Local refinement for regression change points localisation.
#' @description     Perform local refinement for regression change points localisation.
#' @param cpt_init  An \code{integer} vector of initial changepoints estimation (sorted in strictly increasing order).
#' @param y         A \code{numeric} vector of response variable.
#' @param X         A \code{numeric} matrix of covariates with vertical axis being time..
#' @param zeta      A \code{numeric} scalar of tuning parameter for the group lasso.
#' @return  A vector of locally refined change points estimation.
#' @export
#' @author Daren Wang & Haotian Xu
#' @references Rinaldo, A., Wang, D., Wen, Q., Willett, R., & Yu, Y. (2021, March). Localizing changes in high-dimensional regression models. In International Conference on Artificial Intelligence and Statistics (pp. 2089-2097). PMLR.
#' @examples
#' d0 = 10
#' p = 20
#' n = 100
#' cpt_true = c(30, 70)
#' data = simu.change.regression(d0, cpt_true, p, n, sigma = 1, kappa = 9)
#' gamma_set = c(0.01, 0.1, 1)
#' lambda_set = c(0.01, 0.1, 1, 3)
#' temp = CV.search.DP.regression(y = data$y, X = data$X, gamma_set, lambda_set, delta = 2)
#' temp$test_error # test error result
#' # find the indices of gamma_set and lambda_set which minimizes the test error
#' min_idx = as.vector(arrayInd(which.min(temp$test_error), dim(temp$test_error)))
#' gamma_set[min_idx[1]]
#' lambda_set[min_idx[2]]
#' cpt_init = unlist(temp$cpt_hat[min_idx[1], min_idx[2]])
#' local.refine.regression(cpt_init, data$y, X = data$X, zeta = 0.5)
#' @references Rinaldo, Wang, Wen, Willett and Yu (2020) <arxiv:2010.10410>
local.refine.regression = function(cpt_init, y, X, zeta){
  w = 0.9
  n = nrow(X)
  cpt_init_ext = c(0, cpt_init, n)
  cpt_init_numb = length(cpt_init)
  cpt_refined = rep(0, cpt_init_numb+1)
  for (k in 1:cpt_init_numb){
    s = w*cpt_init_ext[k] + (1-w)*cpt_init_ext[k+1]
    e = (1-w)*cpt_init_ext[k+1] + w*cpt_init_ext[k+2]
    lower = ceiling(s) + 1
    upper = floor(e) - 1
    b = sapply(lower:upper, function(eta)obj.LR.regression(ceiling(s), floor(e), eta, y, X, zeta))
    cpt_refined[k+1] = ceiling(s) + which.min(b)
  }
  return(cpt_refined[-1])
}


#' @title Internal Function: An objective function for selecting the best splitting location in the local refinement.
#' @description See equation (4) of the reference.
#' @param s_extra   An \code{integer} scalar of extrapolated starting index.
#' @param e_extra   An \code{integer} scalar of extrapolated ending index.
#' @param y         A \code{numeric} vector of response variable.
#' @param X         A \code{numeric} matrix of covariates with vertical axis being time.
#' @param zeta      A \code{numeric} scalar of tuning parameter for the group lasso.
#' @noRd
obj.LR.regression = function(s_extra, e_extra, eta, y, X, zeta){
  n = nrow(X)
  p = ncol(X)
  group = rep(1:p, 2)
  X_convert = X.glasso.converter.regression(X[s_extra:e_extra,], eta, s_extra)
  y_convert = y[s_extra:e_extra]
  lambda_LR = zeta*sqrt(log(max(n, p)))
  auxfit = gglasso(x = X_convert, y = y_convert, group = group, loss="ls",
                   lambda = lambda_LR/(e_extra-s_extra+1), intercept = FALSE)
  coef = as.vector(auxfit$beta)
  coef1 = coef[1:p]
  coef2 = coef[(p+1):(2*p)]
  btemp = norm(y_convert - X_convert %*% coef, type = "2")^2 + lambda_LR*sum(sqrt(coef1^2 + coef2^2))
  return(btemp)
}


#' @title Internal Function: Compute prediction error based on different zeta.
#' @param y       A \code{numeric} vector of response variable.
#' @param X       A \code{numeric} matrix of covariates with vertical axis being time.
#' @param lower   An \code{integer} scalar of starting index.
#' @param upper   An \code{integer} scalar of ending index.
#' @param zeta    A \code{numeric} scalar of tuning parameter for group lasso.
#' @noRd
distance.CV.LR = function(y, X, lower, upper, zeta){
  n = nrow(X)
  p = ncol(X)
  lambda_LR = zeta*sqrt(log(max(n,p)))
  fit = glmnet(x = X[lower:upper,], y = y[lower:upper], lambda = lambda_LR)
  yhat = X[lower:upper,] %*% as.vector(fit$beta)
  d = norm(y[lower:upper] - yhat, type = "2")
  result = list("MSE" = d^2, "beta" = as.vector(fit$beta))
  return(result)
}                                   


#' @title Internal function: Cross-validation of local refinement for regression.
#' @description     Perform cross-validation of local refinement for regression.
#' @param y         A \code{numeric} vector of observations.
#' @param X         A \code{numeric} matrix of covariates.
#' @param gamma     A \code{numeric} scalar of the tuning parameter associated with the l0 penalty.
#' @param lambda    A \code{numeric} scalar of tuning parameter for the lasso penalty.
#' @param zeta      A \code{numeric} scalar of tuning parameter for the group lasso.
#' @param delta     A strictly positive \code{integer} scalar of minimum spacing.
#' @param eps           A \code{numeric} scalar of precision level for convergence of lasso.
#' @return  A \code{list} with the structure:
#' \itemize{
#'  \item{cpt_hat}: A list of vector of estimated change points locations (sorted in strictly increasing order).
#'  \item{K_hat}: A list of scalar of number of estimated change points.
#'  \item{test_error}: A list of vector of testing errors.
#'  \item{train_error}: A list of vector of training errors.
#' } 
#' @noRd
CV.DP.LR.regression = function(y, X, gamma, lambda, zeta, delta, eps = 0.001){
  n = nrow(X)
  even_indexes = seq(2,n,2)
  odd_indexes = seq(1,n,2)
  train.X = X[odd_indexes,]
  train.y = y[odd_indexes]
  validation.X = X[even_indexes,]
  validation.y = y[even_indexes]
  init_cpt_train = DP.regression(train.y, train.X, gamma, lambda, delta, eps)$cpt
  if(length(init_cpt_train) != 0){
    init_cp_dp = odd_indexes[init_cpt_train]
    init_cp = local.refine.regression(init_cp_dp, y, X, zeta)
  }
  else{
    init_cp_dp = c()
    init_cp = c()
  }
  len = length(init_cp)
  init_cp_train = (1+init_cp_dp)/2
  init_cp_long = c(init_cp_train, n/2)
  interval = matrix(0, nrow = len + 1, ncol = 2)
  interval[1,] = c(1, init_cp_long[1])
  if (len > 0){
    for (j in 2:(1+len)){
      interval[j,] = c(init_cp_long[j-1]+1, init_cp_long[j])
    }
  }
  p = ncol(train.X)
  trainmat = sapply(1:(len+1), function(index) distance.CV.LR(train.y, train.X, interval[index,1], interval[index,2], zeta))
  betamat = matrix(0, nrow = p, ncol = len+1)
  training_loss = matrix(0, nrow = 1, ncol = len+1)                
  for(col in 1:(len+1)){
    betamat[,col] = as.numeric(trainmat[2,col]$beta)
    training_loss[,col] = as.numeric(trainmat[1,col]$MSE)
  }      
  validationmat = sapply(1:(len+1),function(index) error.test.regression(validation.y, validation.X, interval[index,1], interval[index,2], betamat[,index]))                       
  result = list(cpt_hat = init_cp, K_hat = len, test_error = sum(validationmat), train_error = sum(training_loss))                       
  return(result)
}   



#' @title Internal function: Grid search based on Cross-Validation (only gamma and lambda) of local refinement for regression.
#' @param y             A \code{numeric} vector of response variable.
#' @param X             A \code{numeric} matrix of covariates.
#' @param gamma.set     A \code{numeric} vector of candidate tuning parameter associated with the l0 penalty.
#' @param lambda.set    A \code{numeric} vector of candidate tuning parameter for the lasso penalty.
#' @param zeta          A \code{numeric} scalar of tuning parameter for the group lasso.
#' @param delta         A strictly positive \code{integer} scalar of minimum spacing.
#' @param eps           A \code{numeric} scalar of precision level for convergence of lasso.
#' @return  A \code{list} with the following structure:
#'  \item{cpt_hat}{A list of vector of estimated changepoints (sorted in strictly increasing order)}
#'  \item{K_hat}{A list of scalar of number of estimated changepoints}
#'  \item{test_error}{A list of vector of testing errors (each row corresponding to each gamma, and each column corresponding to each lambda)}
#'  \item{train_error}{A list of vector of training errors}
#' @noRd
CV.search.DP.LR.gl = function(y, X, gamma.set, lambda.set, zeta, delta, eps = 0.001){
  output = sapply(1:length(lambda.set), function(i) sapply(1:length(gamma.set), 
                                                           function(j) CV.DP.LR.regression(y, X, gamma.set[j], lambda.set[i], zeta, delta, eps)))
  cpt_hat = output[seq(1,4*length(gamma.set),4),]## estimated change points
  K_hat = output[seq(2,4*length(gamma.set),4),]## number of estimated change points
  test_error = output[seq(3,4*length(gamma.set),4),]## validation loss
  train_error = output[seq(4,4*length(gamma.set),4),]## training loss                                                         
  result = list(cpt_hat = cpt_hat, K_hat = K_hat, test_error = test_error, train_error = train_error)
  return(result)
}                           



#' @title Grid search based on Cross-Validation of all tuning parameters (gamma, lambda and zeta) for regression.
#' @description Perform grid search based on Cross-Validation of all tuning parameters (gamma, lambda and zeta)
#' @param y             A \code{numeric} vector of response variable.
#' @param X             A \code{numeric} matrix of covariates with vertical axis being time.
#' @param gamma_set     A \code{numeric} vector of candidate tuning parameter associated with the l0 penalty.
#' @param lambda_set    A \code{numeric} vector of candidate tuning parameter for the lasso penalty.
#' @param zeta_set      A \code{numeric} vector of candidate tuning parameter for the group lasso.
#' @param delta         A strictly positive \code{integer} scalar of minimum spacing.
#' @param eps           A \code{numeric} scalar of precision level for convergence of lasso.
#' @return  A \code{list} with the following structure:
#'  \item{cpt_hat}{A list of vector of estimated changepoints (sorted in strictly increasing order)}
#'  \item{K_hat}{A list of scalar of number of estimated changepoints}
#'  \item{test_error}{A list of vector of testing errors (each row corresponding to each gamma, and each column corresponding to each lambda)}
#'  \item{train_error}{A list of vector of training errors}
#' @export
#' @author Daren Wang & Haotian Xu
#' @examples
#' set.seed(123)
#' d0 = 8
#' p = 15
#' n = 100
#' cpt_true = c(30, 70)
#' data = simu.change.regression(d0, cpt_true, p, n, sigma = 1, kappa = 9)
#' gamma_set = c(0.01, 0.1)
#' lambda_set = c(0.01, 0.1)
#' temp = CV.search.DP.regression(y = data$y, X = data$X, gamma_set, lambda_set, delta = 2)
#' temp$test_error # test error result
#' # find the indices of gamma_set and lambda_set which minimizes the test error
#' min_idx = as.vector(arrayInd(which.min(temp$test_error), dim(temp$test_error))) 
#' gamma_set[min_idx[1]]
#' lambda_set[min_idx[2]]
#' cpt_init = unlist(temp$cpt_hat[min_idx[1], min_idx[2]])
#' zeta_set = c(0.1, 1)
#' temp_zeta = CV.search.DP.LR.regression(data$y, data$X, gamma_set[min_idx[1]],
#'                   lambda_set[min_idx[2]], zeta_set, delta = 2, eps = 0.001)
#' min_zeta_idx = which.min(unlist(temp_zeta$test_error))
#' cpt_LR = local.refine.regression(cpt_init, data$y, X = data$X, zeta = zeta_set[min_zeta_idx])
#' Hausdorff.dist(cpt_init, cpt_true)
#' Hausdorff.dist(cpt_LR, cpt_true)
#' @references Rinaldo, Wang, Wen, Willett and Yu (2020) <arxiv:2010.10410>
CV.search.DP.LR.regression = function(y, X, gamma_set, lambda_set, zeta_set, delta, eps = 0.001){
  cpt_hat = vector("list", length(zeta_set))
  K_hat = vector("list", length(zeta_set))
  test_error = vector("list", length(zeta_set))
  train_error = vector("list", length(zeta_set))
  for(ii in 1:length(zeta_set)){
    temp = CV.search.DP.LR.gl(y, X, gamma_set, lambda_set, zeta_set[ii], delta, eps = eps)
    cpt_hat[[ii]] = temp$cpt_hat
    K_hat[[ii]] = temp$K_hat
    test_error[[ii]] = temp$test_error
    train_error[[ii]] = temp$train_error
  }
  return(list(cpt_hat = cpt_hat, K_hat = K_hat, test_error = test_error, train_error = train_error))                  
}                         



#' @title Internal Function: Convert a p-by-n design submatrix X with partial consecutive observations into a n-by-(2p) matrix, which fits the group lasso.
#' @description See equation (4) of the reference.
#' @param  X         A \code{numeric} matrix of covariates with partial consecutive observations.
#' @param  eta       A \code{integer} scalar of splitting index.
#' @param  s_ceil    A \code{integer} scalar of starting index.
#' @return A n-by-(2p) matrix
#' @noRd
X.glasso.converter.regression = function(X, eta, s_ceil){
  n = nrow(X)
  xx1 = xx2 = X
  t = eta - s_ceil + 1
  xx1[(t+1):n,] = 0
  xx2[1:t,] = 0
  xx = cbind(xx1/sqrt(t-1), xx2/sqrt(n-t))
  return(xx)
}





#' @title Dynamic programming with dynamic update algorithm for regression change points localisation with \eqn{l_0} penalisation.
#' @description     Perform DPDU algorithm for regression change points localisation.
#' @param y         A \code{numeric} vector of response variable.
#' @param X         A \code{numeric} matrix of covariates with vertical axis being time.
#' @param lambda    A positive \code{numeric} scalar of tuning parameter for lasso penalty.
#' @param zeta      A positive \code{integer} scalar of tuning parameter associated with \eqn{l_0} penalty (minimum interval size).
#' @param eps       A \code{numeric} scalar of precision level for convergence of lasso.
#' @return An object of \code{\link[base]{class}} "DP", which is a \code{list} with the following structure:
#'  \item{partition}{A vector of the best partition.}
#'  \item{cpt}{A vector of change points estimation.}
#' @export
#' @author Haotian Xu
#' @references Xu, Wang, Zhao and Yu (2022) <arXiv:2207.12453>.
#' @examples
#' d0 = 10
#' p = 20
#' n = 100
#' cpt_true = c(30, 70)
#' data = simu.change.regression(d0, cpt_true, p, n, sigma = 1, kappa = 9)
#' temp = DPDU.regression(y = data$y, X = data$X, lambda = 1, zeta = 20)
#' cpt_hat = temp$cpt
#' @export
DPDU.regression <- function(y, X, lambda, zeta, eps = 0.001) {
  DPDU_result = .Call('_changepoints_rcpp_DPDU_regression', PACKAGE = 'changepoints', y, X, lambda, zeta, eps)
  cpt_est = part2local(DPDU_result$partition)
  if(length(cpt_est) >= 1){
    if(min(cpt_est) < zeta){
      cpt_est = cpt_est[-1]
    }
  }
  if(length(cpt_est) >= 1){
    if(length(y) - max(cpt_est) < zeta){
      cpt_est = cpt_est[-1]
    }
  }
  result = append(DPDU_result, list(cpt = cpt_est))
  class(result) = "DP"
  return(result)
}


#' @noRd
lassoDPDU_error <- function(y, X, beta_hat) {
  .Call('_changepoints_rcpp_lassoDPDU_error', PACKAGE = 'changepoints', y, X, beta_hat)
}


#' @title Internal function: Cross-validation for DPDU.
#' @description     Perform cross-validation of dynamic programming algorithm for regression change points.
#' @param y         A \code{numeric} vector of response variable.
#' @param X         A \code{numeric} matrix of covariates with vertical axis being time.
#' @param lambda    A positive \code{numeric} scalar of tuning parameter for lasso penalty.
#' @param zeta      A positive \code{integer} scalar of tuning parameter associated with \eqn{l_0} penalty (minimum interval size).
#' @param eps       A \code{numeric} scalar of precision level for convergence of lasso.
#' @return  A \code{list} with the following structure:
#'  \item{cpt_hat}{A vector of estimated change points locations (sorted in strictly increasing order)}
#'  \item{K_hat}{A scalar of number of estimated change points}
#'  \item{test_error}{A scalar of testing errors in squared \eqn{l_2} norm}
#'  \item{train_error}{A scalar of training errors in squared \eqn{l_2} norm}
#' @noRd
CV.DPDU.regression = function(y, X, lambda, zeta, eps = 0.001){
  N = nrow(X)
  p = ncol(X)
  even_indexes = seq(2, N, 2)
  odd_indexes = seq(1, N, 2)
  train.X = X[odd_indexes,]
  train.y = y[odd_indexes]
  validation.X = X[even_indexes,]
  validation.y = y[even_indexes]
  init_train = DPDU.regression(train.y, train.X, lambda, zeta, eps)
  init_train_cpt = init_train$cpt
  if(length(init_train_cpt) >= 1){
    init_train_cpt_long = c(0, init_train_cpt, nrow(train.X))
    init_train_beta = init_train$beta_mat[,c(init_train_cpt, nrow(train.X))]
    len = length(init_train_cpt_long)-1
    train_error = 0
    test_error = 0
    init_test_cpt_long = c(0, init_train_cpt, nrow(validation.X))
    for(i in 1:len){
      train_error = train_error + lassoDPDU_error(train.y[(init_train_cpt_long[i]+1):init_train_cpt_long[i+1]], cbind(rep(1, nrow(train.X)), train.X)[(init_train_cpt_long[i]+1):init_train_cpt_long[i+1],], init_train_beta[,i])
      test_error = test_error + lassoDPDU_error(validation.y[(init_test_cpt_long[i]+1):init_test_cpt_long[i+1]], cbind(rep(1, nrow(validation.X)), validation.X)[(init_test_cpt_long[i]+1):init_test_cpt_long[i+1],], init_train_beta[,i])
    }
    init_cpt = odd_indexes[init_train_cpt]
    K_hat = len-1
  }else{
    init_cpt = init_train_cpt
    K_hat = 0
    init_train_beta = init_train$beta_mat[,1]
    train_error = lassoDPDU_error(train.y, cbind(rep(1, nrow(train.X)), train.X), init_train_beta)
    test_error = lassoDPDU_error(validation.y, cbind(rep(1, nrow(validation.X)), validation.X), init_train_beta)
  }
  result = list(cpt_hat = init_cpt, K_hat, test_error = test_error, train_error = train_error, beta_hat = init_train_beta)
  return(result)
}


#' @title Grid search based on cross-validation of dynamic programming for regression change points localisation with \eqn{l_0} penalisation.
#' @description Perform grid search to select tuning parameters gamma (for \eqn{l_0} penalty of DP) and lambda (for lasso penalty) based on cross-validation.
#' @param y             A \code{numeric} vector of response variable.
#' @param X             A \code{numeric} matrix of covariates with vertical axis being time.
#' @param lambda_set    A \code{numeric} vector of candidate tuning parameters for lasso penalty.
#' @param zeta_set      An \code{integer} vector of tuning parameter associated with \eqn{l_0} penalty (minimum interval size).
#' @param eps           A \code{numeric} scalar of precision level for convergence of lasso.
#' @return  A \code{list} with the following structure:
#'  \item{cpt_hat}{A list of vectors of estimated change points}
#'  \item{K_hat}{A list of scalars of number of estimated change points}
#'  \item{test_error}{A matrix of testing errors (each row corresponding to each gamma, and each column corresponding to each lambda)}
#'  \item{train_error}{A matrix of training errors}
#'  \item{beta_hat}{A list of matrices of estimated regression coefficients}
#' @export
#' @author Haotian Xu
#' @examples
#' d0 = 5
#' p = 30
#' n = 200
#' cpt_true = 100
#' data = simu.change.regression(d0, cpt_true, p, n, sigma = 1, kappa = 9)
#' lambda_set = c(0.01, 0.1, 1, 2)
#' zeta_set = c(10, 15, 20)
#' temp = CV.search.DPDU.regression(y = data$y, X = data$X, lambda_set, zeta_set)
#' temp$test_error # test error result
#' # find the indices of lambda_set and zeta_set which minimizes the test error
#' min_idx = as.vector(arrayInd(which.min(temp$test_error), dim(temp$test_error))) 
#' lambda_set[min_idx[2]]
#' zeta_set[min_idx[1]]
#' cpt_init = unlist(temp$cpt_hat[min_idx[1], min_idx[2]])
#' beta_hat = matrix(unlist(temp$beta_hat[min_idx[1], min_idx[2]]), ncol = length(cpt_init)+1)
#' @references Xu, Wang, Zhao and Yu (2022) <arXiv:2207.12453>.
CV.search.DPDU.regression = function(y, X, lambda_set, zeta_set, eps = 0.001){
  output = sapply(1:length(lambda_set), function(i) sapply(1:length(zeta_set), 
                                                           function(j) CV.DPDU.regression(y, X, lambda_set[i], zeta_set[j])))
  cpt_hat = output[seq(1,5*length(zeta_set),5),]## estimated change points
  K_hat = output[seq(2,5*length(zeta_set),5),]## number of estimated change points
  test_error = output[seq(3,5*length(zeta_set),5),]## validation loss
  train_error = output[seq(4,5*length(zeta_set),5),]## training loss
  beta_hat = output[seq(5,5*length(zeta_set),5),]
  result = list(cpt_hat = cpt_hat, K_hat = K_hat, test_error = test_error, train_error = train_error, beta_hat = beta_hat)
  return(result)
}


#' @title Local refinement for DPDU regression change points localisation.
#' @description     Perform local refinement for regression change points localisation.
#' @param cpt_init  An \code{integer} vector of initial changepoints estimation (sorted in strictly increasing order).
#' @param beta_hat  A \code{numeric} (px(K_hat+1))matrix of estimated regression coefficients.
#' @param y         A \code{numeric} vector of response variable.
#' @param X         A \code{numeric} matrix of covariates with vertical axis being time..
#' @param w         A \code{numeric} scalar in (0,1) representing the weight for interval truncation.
#' @return  A vector of locally refined change points estimation.
#' @export
#' @author Haotian Xu
#' @references Xu, Wang, Zhao and Yu (2022) <arXiv:2207.12453>.
#' @examples
#' d0 = 5
#' p = 30
#' n = 200
#' cpt_true = 100
#' data = simu.change.regression(d0, cpt_true, p, n, sigma = 1, kappa = 9)
#' lambda_set = c(0.01, 0.1, 1, 2)
#' zeta_set = c(10, 15, 20)
#' temp = CV.search.DPDU.regression(y = data$y, X = data$X, lambda_set, zeta_set)
#' temp$test_error # test error result
#' # find the indices of lambda_set and zeta_set which minimizes the test error
#' min_idx = as.vector(arrayInd(which.min(temp$test_error), dim(temp$test_error))) 
#' lambda_set[min_idx[2]]
#' zeta_set[min_idx[1]]
#' cpt_init = unlist(temp$cpt_hat[min_idx[1], min_idx[2]])
#' beta_hat = matrix(unlist(temp$beta_hat[min_idx[1], min_idx[2]]), ncol = length(cpt_init)+1)
#' cpt_refined = local.refine.DPDU.regression(cpt_init, beta_hat, data$y, data$X, w = 0.9)
#' @references Xu, Wang, Zhao and Yu (2022) <arXiv:2207.12453>.
local.refine.DPDU.regression = function(cpt_init, beta_hat, y, X, w = 0.9){
  n = nrow(X)
  cpt_init_long = c(0, cpt_init, n)
  cpt_init_numb = length(cpt_init)
  cpt_refined = rep(0, cpt_init_numb+1)
  for (k in 1:cpt_init_numb){
    s = w*cpt_init_long[k] + (1-w)*cpt_init_long[k+1]
    e = (1-w)*cpt_init_long[k+1] + w*cpt_init_long[k+2]
    lower = ceiling(s) + 2
    upper = floor(e) - 2
    b = sapply(lower:upper, function(eta)(lassoDPDU_error(y[ceiling(s):eta], cbind(rep(1, n), X)[ceiling(s):eta,], beta_hat[,k]) + lassoDPDU_error(y[(eta+1):floor(e)], cbind(rep(1, n), X)[(eta+1):floor(e),], beta_hat[,k+1])))
    cpt_refined[k+1] = ceiling(s) + which.min(b)
  }
  return(cpt_refined[-1])
}



#' @title Interval trimming based on initial change point localisation.
#' @description     Performing the interval trimming for local refinement.
#' @param n         An \code{integer} scalar corresponding to the sample size.
#' @param cpt_init  An \code{integer} vector of initial changepoints estimation (sorted in strictly increasing order).
#' @param w         A \code{numeric} scalar in (0,1) representing the weight for interval truncation.
#' @return  A matrix with each row be a trimmed interval.
#' @export
#' @author Haotian Xu
#' @references Xu, Wang, Zhao and Yu (2022) <arXiv:2207.12453>.
trim_interval = function(n, cpt_init, w = 0.9){
  cpt_init_long = c(0, cpt_init, n)
  cpt_init_numb = length(cpt_init)
  interval_mat = matrix(NA, nrow = cpt_init_numb, ncol = 2)
  for (k in 1:cpt_init_numb){
    interval_mat[k,1] = w*cpt_init_long[k] + (1-w)*cpt_init_long[k+1]
    interval_mat[k,2] = (1-w)*cpt_init_long[k+1] + w*cpt_init_long[k+2]
  }
  return(interval_mat)
}





#' @title Long-run variance estimation for regression settings with change points.
#' @description     Estimating long-run variance for regression settings with change points.
#' @param cpt_init   An \code{integer} vector of initial changepoints estimation (sorted in strictly increasing order).
#' @param beta_hat   A \code{numeric} (px(K_hat+1))matrix of estimated regression coefficients.
#' @param y          A \code{numeric} vector of response variable.
#' @param X          A \code{numeric} matrix of covariates with vertical axis being time.
#' @param w          A \code{numeric} scalar in (0,1) representing the weight for interval truncation.
#' @param block_size An \code{integer} scalar corresponding to the block size S in the paper.
#' @return  A vector of long-run variance estimators associated with all local refined intervals.
#' @export
#' @author Haotian Xu
#' @references Xu, Wang, Zhao and Yu (2022) <arXiv:2207.12453>.
#' @examples
#' d0 = 5
#' p = 10
#' n = 200
#' cpt_true = c(70, 140)
#' data = simu.change.regression(d0, cpt_true, p, n, sigma = 1, kappa = 9, cov_type = "T", mod_X = "MA", mod_e = "AR")
#' lambda_set = c(0.1, 0.5, 1, 2)
#' zeta_set = c(10, 15, 20)
#' temp = CV.search.DPDU.regression(y = data$y, X = data$X, lambda_set, zeta_set)
#' temp$test_error # test error result
#' # find the indices of lambda_set and zeta_set which minimizes the test error
#' min_idx = as.vector(arrayInd(which.min(temp$test_error), dim(temp$test_error))) 
#' lambda_set[min_idx[2]]
#' zeta_set[min_idx[1]]
#' cpt_init = unlist(temp$cpt_hat[min_idx[1], min_idx[2]])
#' beta_hat = matrix(unlist(temp$beta_hat[min_idx[1], min_idx[2]]), ncol = length(cpt_init)+1)
#' interval_refine = trim_interval(n, cpt_init)
#' # choose S
#' block_size = ceiling(sqrt(min(floor(interval_refine[,2]) - ceiling(interval_refine[,1])))/2)
#' LRV_est = LRV.regression(cpt_init, beta_hat, data$y, data$X, w = 0.9, block_size)
#' @references Xu, Wang, Zhao and Yu (2022) <arXiv:2207.12453>.
LRV.regression = function(cpt_init, beta_hat, y, X, w = 0.9, block_size){
  n = nrow(X)
  cpt_init_long = c(0, cpt_init, n)
  cpt_init_numb = length(cpt_init)
  X_full = cbind(rep(1, n), X)
  lrv_hat = rep(NA, cpt_init_numb)
  for (k in 1:cpt_init_numb){
    kappa2_hat = sum((beta_hat[,k] - beta_hat[,k+1])^2)
    s = w*cpt_init_long[k] + (1-w)*cpt_init_long[k+1]
    e = (1-w)*cpt_init_long[k+1] + w*cpt_init_long[k+2]
    z_vec = rep(NA, floor(e)-ceiling(s)+1)
    for (t in ceiling(s):floor(e)){
      z_vec[t-ceiling(s)+1] = as.numeric(2*y[t] - crossprod(X_full[t,], beta_hat[,k]+beta_hat[,k+1]))*crossprod(X_full[t,], beta_hat[,k+1]-beta_hat[,k])
    }
    pair_numb = floor((floor(e)-ceiling(s)+1)/(2*block_size))
    z_mat1 = matrix(z_vec[1:(block_size*pair_numb*2)], nrow = block_size)
    z_mat1_colsum = apply(z_mat1, 2, sum)
    z_mat2 = matrix(rev(z_vec)[1:(block_size*pair_numb*2)], nrow = block_size)
    z_mat2_colsum = apply(z_mat2, 2, sum)
    lrv_hat[k] = (mean((z_mat1_colsum[2*(1:pair_numb)-1] - z_mat1_colsum[2*(1:pair_numb)])^2/(2*block_size)) + mean((z_mat2_colsum[2*(1:pair_numb)-1] - z_mat2_colsum[2*(1:pair_numb)])^2/(2*block_size)))/(2*kappa2_hat)
  }
  return(lrv_hat)
}


#' @title Internal Function: simulate a two-sided Brownian motion with drift.
#' @param  n         An \code{integer} scalar representing the length of one side.
#' @param  drift     A \code{numeric} scalar of drift coefficient.
#' @param  LRV    A \code{integer} scalar of LRV.
#' @return A (2n+1) vector
#' @noRd
simu.2BM_Drift = function(n, drift, LRV){
  z_vec = rnorm(2*n)
  w_vec = c(rev(cumsum(z_vec[n:1])/sqrt(1:n)), 0, cumsum(z_vec[(n+1):(2*n)])/sqrt(1:n))
  v_vec = drift*abs(seq(-n, n)) + sqrt(LRV)*w_vec
  return(v_vec)
}


#' @title Confidence interval construction of change points for regression settings with change points.
#' @description     Construct element-wise confidence interval for change points.
#' @param cpt_init   An \code{integer} vector of initial changepoints estimation (sorted in strictly increasing order).
#' @param cpt_LR     An \code{integer} vector of refined changepoints estimation (sorted in strictly increasing order).
#' @param beta_hat   A \code{numeric} (px(K_hat+1))matrix of estimated regression coefficients.
#' @param y          A \code{numeric} vector of response variable.
#' @param X          A \code{numeric} matrix of covariates with vertical axis being time.
#' @param w          A \code{numeric} scalar in (0,1) representing the weight for interval truncation.
#' @param B          An \code{integer} scalar corresponding to the number of simulated two-sided Brownian motion with drift.
#' @param M          An \code{integer} scalar corresponding to the length for each side of the limiting distribution, i.e. the two-sided Brownian motion with drift.
#' @param alpha_vec  An \code{numeric} vector in (0,1) representing the vector of significance levels.
#' @param rounding   A \code{boolean} scalar representing if the confidence intervals need to be rounded into integer intervals.
#' @return  An length(cpt_init)-2-length(alpha_vec) array of confidence intervals.
#' @export
#' @author Haotian Xu
#' @references Xu, Wang, Zhao and Yu (2022) <arXiv:2207.12453>.
#' @examples
#' d0 = 5
#' p = 10
#' n = 200
#' cpt_true = c(70, 140)
#' data = simu.change.regression(d0, cpt_true, p, n, sigma = 1, kappa = 9, cov_type = "T", mod_X = "MA", mod_e = "AR")
#' lambda_set = c(0.1, 0.5, 1, 2)
#' zeta_set = c(10, 15, 20)
#' temp = CV.search.DPDU.regression(y = data$y, X = data$X, lambda_set, zeta_set)
#' temp$test_error # test error result
#' # find the indices of lambda_set and zeta_set which minimizes the test error
#' min_idx = as.vector(arrayInd(which.min(temp$test_error), dim(temp$test_error))) 
#' lambda_set[min_idx[2]]
#' zeta_set[min_idx[1]]
#' cpt_init = unlist(temp$cpt_hat[min_idx[1], min_idx[2]])
#' beta_hat = matrix(unlist(temp$beta_hat[min_idx[1], min_idx[2]]), ncol = length(cpt_init)+1)
#' cpt_LR = local.refine.DPDU.regression(cpt_init, beta_hat, data$y, data$X, w = 0.9)
#' alpha_vec = c(0.01, 0.05, 0.1)
#' CI.regression(cpt_init, cpt_LR, beta_hat, data$y, data$X, w = 0.9, B = 1000, M = n, alpha_vec)
#' @references Xu, Wang, Zhao and Yu (2022) <arXiv:2207.12453>.
CI.regression = function(cpt_init, cpt_LR, beta_hat, y, X, w = 0.9, B = 1000, M, alpha_vec, rounding = TRUE){
  if(length(cpt_init) != length(cpt_LR)){
    stop("The initial and the refined change point estimators should have the same length.")
  }
  n = length(y)
  CI_array = array(NA, c(length(cpt_init), 2, length(alpha_vec)))
  interval_refine_mat = trim_interval(n = n, cpt_init, w = w)
  block_size = ceiling((min(floor(interval_refine_mat[,2]) - ceiling(interval_refine_mat[,1])))^(2/5)/2) # choose S
  LRV_hat = LRV.regression(cpt_init, beta_hat, y, X, w = w, block_size)
  kappa2_hat = apply(beta_hat[,-dim(beta_hat)[2]] - beta_hat[,-1], MARGIN = 2, crossprod)
  drift_hat = diag(t(beta_hat[,-dim(beta_hat)[2]] - beta_hat[,-1]) %*% (t(cbind(rep(1,n), X))%*%cbind(rep(1,n), X)) %*% (beta_hat[,-dim(beta_hat)[2]] - beta_hat[,-1]) / (n*kappa2_hat))
  for(i in 1:length(cpt_init)){
    u_vec = rep(NA, B)
    for(b in 1:B){
      set.seed(12345+b+10000*i)
      u_vec[b] = seq(-M, M)[which.min(simu.2BM_Drift(M, drift_hat[i], LRV_hat[i]))]
    }
    for(j in 1:length(alpha_vec)){
      alpha = alpha_vec[j]
      CI_array[i,,j] = quantile(u_vec, probs = c(alpha/2, 1-alpha/2))/kappa2_hat[i] + cpt_LR[i]
      if(rounding == TRUE){
        CI_array[i,1,j] = floor(CI_array[i,1,j])
        CI_array[i,2,j] = ceiling(CI_array[i,2,j])
      }
    }
  }
  return(CI_array)
}