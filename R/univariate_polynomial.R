#' @title Internal function: Coefficients reparametrization for polynomials with different bases
#' @description Compute the transformed coefficients for polynomials with different bases (currently, only the linear, quadratic and cubic functions are considered). 
#' @param coef_vec  A \code{numeric} vector of coefficients for polynomials associated with cpt1.
#' @param cpt1      An \code{integer} scalar of the first change point.
#' @param cpt2      An \code{integer} scalar of the second change point.
#' @param n         An \code{integer} scalar of sample size.
#' @return A vector of transformed coefficients for polynomials associated with cpt2.
#' @noRd
#' @references Yu and Chatterjee (2020) <arXiv:2007.09910>
#' @author  Haotian Xu
coef.repara = function(coef_vec, cpt1, cpt2, n){
  if(length(coef_vec) == 2){
    A = matrix(c(1, (cpt2-cpt1)/n, 0, 1), nrow = 2, byrow = TRUE)
  }else if(length(coef_vec) == 3){
    A = matrix(c(1, (cpt2-cpt1)/n, ((cpt2-cpt1)/n)^2, 0, 1, 2*(cpt2-cpt1)/n, 0, 0, 1), nrow = 3, byrow = TRUE)
  }else if(length(coef_vec) == 4){
    A = matrix(c(1, (cpt2-cpt1)/n, ((cpt2-cpt1)/n)^2, ((cpt2-cpt1)/n)^3, 0, 1, 2*(cpt2-cpt1)/n, 3*((cpt2-cpt1)/n)^2, 0, 0, 1, 3*((cpt2-cpt1)/n), 0, 0, 0, 1), nrow = 4, byrow = TRUE)
  }else{
    stop("Currently, only the linear, quadratic and cubic functions are considered.")
  }
  return(A %*% coef_vec)
}



#' @title Generate univariate data from piecewise polynomials of degree at most r.
#' @description Generate univariate data from piecewise polynomials (currently, only the linear, quadratic functions and cubic functions are considered). 
#' @param init_coef_vec  A (r+1)-dim \code{numeric} vector of coefficients for the first segment.
#' @param cpt_vec        A K-dim \code{integer} vector of change points.
#' @param kappa_mat      A (r+1)xK \code{numeric} matrix where the i-th column represents the jump sizes for coefficients associated with the i-th change point.
#' @param n              An \code{integer} scalar of sample size.
#' @param sigma          A \code{numeric} scalar of standard deviation of error terms.
#' @return A vector of data generated from piecewise polynomials.
#' @export
#' @references Yu and Chatterjee (2020) <arXiv:2007.09910>.
#' @author  Haotian Xu
#' @examples
#' r = 2
#' init_coef_vec = c(-2, 2, 9)
#' cpt_true = c(100, 200)
#' n = 300
#' sigma = 1
#' kappa_mat = cbind(c(3, 9, -27), c(-3, 9, -27))
#' plot.ts(gen.piece.poly(init_coef_vec, cpt_true, kappa_mat, n, sigma), ylab = "y")
gen.piece.poly = function(init_coef_vec, cpt_vec, kappa_mat, n, sigma){
  r = length(init_coef_vec) - 1
  cpt_ext = c(cpt_vec, n)
  if(any(dim(kappa_mat) != c(length(init_coef_vec), length(cpt_vec)))){
    stop("kappa_mat is not correct")
  }
  result = sapply((1:cpt_vec[1])/n, function(x){sum((x - cpt_vec[1]/n)^(0:r) * init_coef_vec)})
  coef_vec = init_coef_vec + kappa_mat[,1]
  for(i in 1:length(cpt_vec)){
    result = c(result, sapply(((cpt_ext[i]+1):cpt_ext[i+1])/n, function(x){sum((x - cpt_vec[i]/n)^(0:r) * coef_vec)}))
    if(i <= length(cpt_vec)-1){
      coef_vec = coef.repara(coef_vec, cpt_vec[i], cpt_vec[i+1], n) + kappa_mat[,i+1]
    }
  }
  result = result + rnorm(n, mean = 0, sd = sigma)
  return(result)
}


#' @title Mean function of piecewise polynomials.
#' @description Compute mean function of piecewise polynomials (currently, only the linear, quadratic functions and cubic functions are considered). 
#' @param init_coef_vec  A \code{numeric} vector of coefficients for the first segment.
#' @param cpt_vec        An \code{integer} vector of change points.
#' @param kappa_mat      A \code{numeric} matrix where the i-th column represents the jump sizes for coefficients associated with the i-th change point.
#' @param n              An \code{integer} scalar of sample size.
#' @return A vector of mean function of piecewise polynomials.
#' @export
#' @references Yu and Chatterjee (2020) <arXiv:2007.09910>
#' @author  Haotian Xu
#' @examples
#' r = 2
#' init_coef_vec = c(-2, 2, 9)
#' cpt_true = c(100, 200)
#' n = 300
#' kappa_mat = cbind(c(3, 9, -27), c(-3, 9, -27))
#' plot.ts(gen.piece.poly.noiseless(init_coef_vec, cpt_true, kappa_mat, n), 
#'         ylab = "Values of piecewise polynomials")
gen.piece.poly.noiseless = function(init_coef_vec, cpt_vec, kappa_mat, n){
  r = length(init_coef_vec) - 1
  cpt_ext = c(cpt_vec, n)
  if(any(dim(kappa_mat) != c(length(init_coef_vec), length(cpt_vec)))){
    stop("kappa_mat is not correct")
  }
  result = sapply((1:cpt_vec[1])/n, function(x){sum((x - cpt_vec[1]/n)^(0:r) * init_coef_vec)})
  coef_vec = init_coef_vec + kappa_mat[,1]
  for(i in 1:length(cpt_vec)){
    result = c(result, sapply(((cpt_ext[i]+1):cpt_ext[i+1])/n, function(x){sum((x - cpt_vec[i]/n)^(0:r) * coef_vec)}))
    if(i <= length(cpt_vec)-1){
      coef_vec = coef.repara(coef_vec, cpt_vec[i], cpt_vec[i+1], n) + kappa_mat[,i+1]
    }
  }
  result = result# + rnorm(n, mean = 0, sd = sigma)
  return(result)
}



#' @title Internal Function: Generate a polynomial basis matrix of order r.
#' @param n          An \code{integer} scalar of sample size.
#' @param s          An \code{integer} scalar of starting index.
#' @param e          An \code{integer} scalar of ending index.
#' @param r          An \code{integer} scalar representing the order of polynomials.
#' @return  A (e-s+1)-by-(r+1) \code{numeric} basis matrix with (i,j)-th entry being \eqn{[(s+i-1)/n]^j}.
#' @noRd
#' @author  Haotian Xu
#' @examples 
#' poly_mat = basis.poly(100, 1, 50, 2)
basis.poly <- function(n, s, e, r){
  .Call('_changepoints_rcpp_basis_poly', PACKAGE = 'changepoints', n, s, e, r)
}


#' @title Dynamic programming algorithm for univariate polynomials change points detection. 
#' @description Perform dynamic programming algorithm for univariate polynomials change points detection.
#' @param y         A \code{numeric} vector of observations.
#' @param r         An \code{integer} scalar order of polynomials.
#' @param gamma     A \code{numeric} scalar of the tuning parameter associated with the \eqn{l_0} penalty.
#' @param delta     A strictly \code{integer} scalar of minimum spacing.
#' @return A \code{list} with the following structure:
#' @return An object of \code{\link[base]{class}} "DP", which is a \code{list} with the following structure:
#'  \item{partition}{A vector of the best partition.}
#'  \item{yhat}{A vector of mean estimation for corresponding to the best partition.}
#'  \item{cpt}{A vector of change points estimation.}
#' @export
#' @author  Haotian Xu
#' @references Yu and Chatterjee (2020) <arXiv:2007.09910>
#' @examples
#' set.seed(0)
#' cpt_true = c(20, 50, 170)
#' y = rnorm(300) + c(rep(0,20),rep(2,30),rep(0,120),rep(2,130))
#' plot.ts(y)
#' temp = DP.poly(y, r = 2, gamma = 15, delta = 5)
#' temp$cpt
DP.poly <- function(y, r, gamma, delta) {
  DP_result = .Call('_changepoints_rcpp_DP_poly', PACKAGE = 'changepoints', y, r, gamma, delta)
  result = append(DP_result, list(cpt = part2local(DP_result$partition)))
  class(result) = "DP"
  return(result)
}



#' @title Internal function: Cross-Validation of Dynamic Programming algorithm for univariate polynomials change points detection.
#' @description     Perform cross-validation by sample splitting. Using the sample with odd indices as training data to estimate the changepoints, then computing sample estimation for each segment within two consecutive changepoints, and computing the validation error based on the sample with even indices.
#' @param y         A \code{numeric} vector of observations.
#' @param r         An \code{integer} scalar order of polynomials.
#' @param gamma     A \code{numeric} scalar of the tuning parameter associated with the l0 penalty.
#' @param delta     A positive \code{integer} scalar of minimum spacing.
#' @return  A \code{list} with the following structure:
#'  \item{cpt_hat}{A vector of estimated change points locations (sorted in strictly increasing order)}
#'  \item{K_hat}{A scalar of number of estimated change points}
#'  \item{test_error}{A vector of testing errors}
#'  \item{train_error}{A vector of training errors}
#' @noRd
CV.DP.poly = function(y, r, gamma, delta){
  N = length(y)
  even_indexes = seq(2, N, 2)
  odd_indexes = seq(1, N, 2)
  train.y = y[odd_indexes]
  validation.y = y[even_indexes]
  temp = DP.poly(train.y, r, gamma, delta)
  init_cpt_train = temp$cpt
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
#' @references Yu and Chatterjee (2020) <arXiv:2007.09910>
CV.search.DP.poly = function(y, r, gamma_set, delta){
  output = sapply(1:length(gamma_set), function(j) CV.DP.poly(y, r, gamma_set[j], delta))
  #print(output)
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
#' @references Yu and Chatterjee (2020) <arXiv:2007.09910>
local.refine.poly = function(cpt_init, y, r, delta_lr){
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

