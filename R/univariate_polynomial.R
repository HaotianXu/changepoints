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
#' @description TO DO
#' @param y         A \code{numeric} vector of observations.
#' @param r         An \code{integer} scalar order of polynomials.
#' @param gamma     A \code{numeric} scalar of the tuning parameter associated with the l0 penalty.
#' @param delta     A strictly \code{integer} scalar of minimum spacing.
#' @param ...      Additional arguments.
#' @return TO DO.
#' @export
#' @author
#' @examples
#' data = simu.change.regression(10, c(10, 30, 40, 70, 90), 30, 100, 1, 9)
#' DP.regression(2, 5, data$y, X = data$X, lambda = 1)
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


#Two-step estimation
local.refine.poly = function(cpt.init, y, w = 1/2){
  n = length(y)
  cpt_init_ext = c(0, cpt_init, n)
  cpt_init_numb = length(cpt_init)
  cpt_refined = rep(0, cpt_init_numb+1)
  for (k in 1:cpt_init_numb){
    s = w*cpt_init_ext[k] + (1-w)*cpt_init_ext[k+1]
    e = (1-w)*cpt_init_ext[k+1] + w*cpt_init_ext[k+2]
    lower = ceiling(s) + 1
    upper = floor(e) - 1
    b = sapply(lower:upper, function(eta) obj.func.lr.poly(y, s, e, eta))
    cpt.refined[k+1] = ceiling(s) + which.min(b)
  }
  return(cpt.refined[-1])
}


#' @title Internal Function: An objective function to select the best splitting location in the local refinement, see eq(4) in [10]
#' @param y         A \code{numeric} vector of observations.
#' @param s.inter   A \code{numeric} scalar of interpolated starting index.
#' @param e.inter   A \code{numeric} scalar of interpolated ending index.
#' @param eta       An \code{integer} scalar between s.inter and e.inter.
#' @noRd
obj.func.lr.poly = function(y, s.inter, e.inter, eta){
  n = length(y)
  s_star = ceiling(s.inter)
  e_star = floor(e.inter)
  u_mat1 = basis.poly(n, s_star, eta-1, r)
  proj_mat1 = u_mat1 %*% solve(t(u_mat1)%*%u_mat1) %*% t(u_mat1)
  u_mat2 = basis.poly(n, eta, e_star, r)
  proj_mat2 = u_mat2 %*% solve(t(u_mat2)%*%u_mat2) %*% t(u_mat2)
  btemp = norm(y[s_star:(eta-1)] - proj_mat1 %*% y[s_star:(eta-1)], type = "2")^2 + norm(y[eta:e_star] - proj_mat2 %*% y[eta:e_star], type = "2")^2
  return(btemp)
}

