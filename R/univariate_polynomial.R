#' @title Internal Function: Generate a polynomial basis matrix of order r.
#' @param n          An \code{integer} scalar of sample size.
#' @param s          An \code{integer} scalar of starting index.
#' @param e          An \code{integer} scalar of ending index.
#' @param r          An \code{integer} scalar representing the order of polynomials.
#' @return  A (e-s+1)-by-(r+1) \code{numeric} basis matrix.
#' @noRd
basis.poly <- function(n, s, e, r) {
  .Call('_changepoints_rcpp_basis_poly', PACKAGE = 'changepoints', n, s, e, r)
}


#' @title Dynamic programming algorithm for univariate polynomials change points detection. 
#' @description Perform dynamic programming algorithm for univariate polynomials change points detection.
#' @param y         A \code{numeric} vector of observations.
#' @param r         An \code{integer} scalar order of polynomials.
#' @param gamma     A \code{numeric} scalar of the tuning parameter associated with the \eqn{l_0} penalty.
#' @param delta     A strictly \code{integer} scalar of minimum spacing.
#' @param ...       Additional arguments.
#' @return A \code{list} with the following structure:
#'  \item{partition}{A vector of the best partition}
#'  \item{yhat}{A vector of mean estimation for corresponding to the best partition}
#' @export
#' @author  Haotian Xu
#' @examples
#' set.seed(0)
#' cpt_true = c(20, 50, 170)
#' y = rnorm(300) + c(rep(0,20),rep(2,30),rep(0,120),rep(2,130))
#' plot.ts(y)
#' temp = DP.poly(y, r = 2, gamma = 15, delta = 5)
#' part2local(temp$partition)
DP.poly <- function(y, r, gamma, delta, ...) {
  .Call('_changepoints_rcpp_DP_poly', PACKAGE = 'changepoints', y, r, gamma, delta)
}
# DP.poly = function(y, r, gamma, delta, ...){
#   n = length(y)
#   bestvalue = rep(0,n+1)
#   partition = rep(0,n)
#   yhat = rep(NA, n)
#   bestvalue[1] = -gamma
#   for(i in 1:n){
#     bestvalue[i+1] = Inf
#     for(l in 1:i){
#       if(i - l > delta){
#         u_mat = basis.poly(n, l, i, r)
#         proj_mat = u_mat %*% solve(t(u_mat)%*%u_mat) %*% t(u_mat)
#         b = bestvalue[l] + gamma + (norm((diag(1, i-l+1) - proj_mat) %*% y[l:i], type = "2"))^2
#       }else{
#         b = Inf
#       }
#       if(b < bestvalue[i+1]){
#         bestvalue[i+1] = b
#         partition[i] = l-1
#       }
#     }
#   }
#   i = n
#   l = partition[i]
#   while(i > 0){
#     u_mat = basis.poly(n, l+1, i, r)
#     proj_mat = u_mat %*% solve(t(u_mat)%*%u_mat) %*% t(u_mat)
#     yhat[(l+1):i] = proj_mat %*% y[(l+1):i]
#     i = l
#     l = partition[i]
#   }
#   return(list(partition = partition, yhat = yhat))
# }


#' @title Internal function: Cross-Validation of Dynamic Programming algorithm for univariate polynomials change points detection.
#' @description     Perform cross-validation by sample splitting. Using the sample with odd indices as training data to estimate the changepoints, then computing sample estimation for each segment within two consecutive changepoints, and computing the validation error based on the sample with even indices.
#' @param y         A \code{numeric} vector of observations.
#' @param r         An \code{integer} scalar order of polynomials.
#' @param gamma     A \code{numeric} scalar of the tuning parameter associated with the l0 penalty.
#' @param delta     A positive \code{integer} scalar of minimum spacing.
#' @param ...       Additional arguments.
#' @return  A \code{list} with the following structure:
#'  \item{cpt_hat}{A vector of estimated change points locations (sorted in strictly increasing order)}
#'  \item{K_hat}{A scalar of number of estimated change points}
#'  \item{test_error}{A vector of testing errors}
#'  \item{train_error}{A vector of training errors}
#' @noRd
CV.DP.poly = function(y, r, gamma, delta, ...){
  N = length(y)
  even_indexes = seq(2, N, 2)
  odd_indexes = seq(1, N, 2)
  train.y = y[odd_indexes]
  validation.y = y[even_indexes]
  temp = DP.poly(train.y, r, gamma, delta)
  init_cpt_train = part2local(temp$partition)
  y_hat_train = temp$yhat
  init_cpt_train.long = c(0, init_cpt_train, length(train.y))
  diff.point = diff(init_cpt_train.long)
  if (length(which(diff.point == 1)) > 0){
    print(paste("gamma =", gamma, ".", "Warning: Consecutive points detected. Try a larger gamma."))
    init_cpt = odd_indexes[init_cpt_train]
    len = length(init_cpt)
    result = list(cpt_hat = init_cpt, K_hat = len, test_error = Inf, train_error = Inf)
  }
  else{
    init_cpt = odd_indexes[init_cpt_train]
    len = length(init_cpt)
    train_error = norm(train.y - y_hat_train, type = "2")
    test_error = norm(validation.y - y_hat_train[1:length(validation.y)], type = "2")
    result = list(cpt_hat = init_cpt, K_hat = len, test_error, train_error)
  }
  return(result)
}


#' @title Grid search for dynamic programming to select the tuning parameter through Cross-Validation.
#' @description Perform grid search for dynamic programming to select the tuning parameter through Cross-Validation.
#' @param gamma_set     A \code{numeric} vector of candidate tuning parameter associated with the l0 penalty.
#' @param y             A \code{numeric} vector of observations.
#' @param r             An \code{integer} scalar order of polynomials.
#' @param delta         A positive \code{integer} scalar of minimum spacing.
#' @param ...           Additional arguments.
#' @return  A \code{list} with the following structure:
#'  \item{cpt_hat}{A list of vector of estimated change points locations (sorted in strictly increasing order)}
#'  \item{K_hat}{A list of scalar of number of estimated change points}
#'  \item{test_error}{A list of vector of testing errors}
#'  \item{train_error}{A list of vector of training errors}
#' @export
#' @author  Haotian Xu
#' @examples
#' set.seed(0)
#' cpt_true = c(20, 50, 170)
#' y = rnorm(300) + c(rep(0,20),rep(2,30),rep(0,120),rep(2,130))
#' plot.ts(y)
#' gamma_set = 3:9
#' DP_result = CV.search.DP.poly(y, r = 2, gamma_set, delta = 5)
#' min_idx = which.min(DP_result$test_error)
#' cpt_init = unlist(DP_result$cpt_hat[min_idx])
#' local.refine.poly(cpt_init, y, r = 2, delta_lr = 5)
CV.search.DP.poly = function(y, r, gamma_set, delta, ...){
  output = sapply(1:length(gamma_set), function(j) CV.DP.poly(y, r, gamma_set[j], delta))
  print(output)
  cpt_hat = output[1,]## estimated change points
  K_hat = output[2,]## number of estimated change points
  test_error = output[3,]## validation loss
  train_error = output[4,]## training loss                                                      
  result = list(cpt_hat = cpt_hat, K_hat = K_hat, test_error = test_error, train_error = train_error)
  return(result)
}


#' @title Local refinement for univariate polynomials change point detection.
#' @description     Perform local refinement for univariate polynomials change point detection.
#' @param cpt_init  An \code{integer} vector of initial change points estimation (sorted in strictly increasing order).
#' @param y         A \code{numeric} vector of univariate time series.
#' @param r         An \code{integer} scalar order of polynomials.
#' @param delta_lr A positive \code{integer} scalar of minimum spacing for local refinement.
#' @param ...       Additional arguments.
#' @return  An \code{integer} vector of locally refined change point estimation.
#' @export
#' @author  Haotian Xu
#' @examples
#' set.seed(0)
#' cpt_true = c(20, 50, 170)
#' y = rnorm(300) + c(rep(0,20),rep(2,30),rep(0,120),rep(2,130))
#' plot.ts(y)
#' gamma_set = 3:9
#' DP_result = CV.search.DP.poly(y, r = 2, gamma_set, delta = 5)
#' min_idx = which.min(DP_result$test_error)
#' cpt_init = unlist(DP_result$cpt_hat[min_idx])
#' local.refine.poly(cpt_init, y, r = 2, delta_lr = 5)
local.refine.poly = function(cpt_init, y, r, delta_lr, ...){
  w = 0.9
  n = length(y)
  cpt_init_ext = c(0, cpt_init, n)
  cpt_init_numb = length(cpt_init)
  cpt_refined = rep(0, cpt_init_numb+1)
  for (k in 1:cpt_init_numb){
    s = w*cpt_init_ext[k] + (1-w)*cpt_init_ext[k+1]
    e = (1-w)*cpt_init_ext[k+1] + w*cpt_init_ext[k+2]
    lower = ceiling(s) + 1
    upper = floor(e) - 1
    b = sapply(lower:upper, function(eta) obj.func.lr.poly(y, s, e, eta, r, delta_lr))
    cpt_refined[k+1] = ceiling(s) + which.min(b)
  }
  return(cpt_refined[-1])
}


#' @title Internal Function: An objective function to select the best splitting location in the local refinement, see eq(4) in [10]
#' @param y         A \code{numeric} vector of observations.
#' @param s.inter   A \code{numeric} scalar of interpolated starting index.
#' @param e.inter   A \code{numeric} scalar of interpolated ending index.
#' @param eta       An \code{integer} scalar between s.inter and e.inter.
#' @noRd
obj.func.lr.poly = function(y, s.inter, e.inter, eta, r, delta_lr){
  n = length(y)
  s_star = ceiling(s.inter)
  e_star = floor(e.inter)
  if((eta - s_star < 2*delta_lr) | (e_star - eta + 1 < 2*delta_lr)){
    btemp = Inf
  }else{
    u_mat1 = basis.poly(n, s_star, eta-1, r)
    proj_mat1 = u_mat1 %*% solve(t(u_mat1)%*%u_mat1) %*% t(u_mat1)
    u_mat2 = basis.poly(n, eta, e_star, r)
    proj_mat2 = u_mat2 %*% solve(t(u_mat2)%*%u_mat2) %*% t(u_mat2)
    btemp = norm(y[s_star:(eta-1)] - proj_mat1 %*% y[s_star:(eta-1)], type = "2")^2 + norm(y[eta:e_star] - proj_mat2 %*% y[eta:e_star], type = "2")^2
  }
  return(btemp)
}

