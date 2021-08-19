#' @title Internal Function: Compute the Covariance CUSUM statistic.
#' @param X         A \code{numeric} matrix of observations with with horizontal axis being time, and vertical axis being dimensions.
#' @param s         A \code{integer} scalar of starting index.
#' @param e         A \code{integer} scalar of ending index.
#' @param t         A \code{integer} scalar of splitting index.
#' @return  A \code{numeric} matrix of the Covariance CUSUM statistic.
#' @noRd
CUSUM.cov = function(X, s, e, t){
  n_st = t - s
  n_se = e - s
  n_te = e - t
  cov_st = cov(t(X[, (s+1):t]))
  cov_te = cov(t(X[, (t+1):e]))
  result = sqrt(n_st * n_te / n_se) * (cov_st - cov_te)
  return(result)
}


#' @title Binary Segmentation for covariance change points detection through Operator Norm.
#' @description TO DO
#' @param X         A \code{numeric} matrix of observations with with horizontal axis being time, and vertical axis being dimensions.
#' @param s         A \code{integer} scalar of starting index.
#' @param e         A \code{integer} scalar of ending index.
#' @param level     A parameter using to track the level a change point is detected. Should be fixed as 0.
#' @param ...      Additional arguments.
#' @return  A \code{list} with the structure:
#' \itemize{
#'  \item S           A vector of estimated changepoints (sorted in strictly increasing order).
#'  \item Dval        A vector of values of CUSUM statistic based on KS distance.
#'  \item Level       A vector representing the levels at which each change point is detected.
#'  \item Parent      A matrix with the starting indices on the first row and the ending indices on the second row.
#'  \item ...         Additional parameters.
#' } 
#' @export
#' @author Haotian Xu
#' @examples
#' p = 10
#' A1 = gen.cov.mat(p, 1, 2)
#' A2 = gen.cov.mat(p, 2, 1)
#' A3 = gen.cov.mat(p, 3, 2)
#' X = cbind(t(mvrnorm(100, mu = rep(0, p), A1)), t(mvrnorm(150, mu = rep(0, p), A2)), t(mvrnorm(200, mu = rep(0, p), A3)))
#' temp = BS.cov(X, 1, 450)
#' BS.threshold(temp, 10)
BS.cov = function(X, s, e, level = 0, ...){
  p = dim(X)[1]
  n = dim(X)[2]
  delta = 2*p*log(n) + 1
  S = NULL
  Dval = NULL
  Level = NULL
  Parent = NULL
  if(e-s <= delta){
    return(list(S = S, Dval = Dval, Level = Level, Parent = Parent))
  }else{
    level = level + 1
    parent = matrix(c(s, e), nrow = 2)
    s_star = ceiling(s + p*log(n))
    e_star = floor(e - p*log(n))
    a = rep(0, e_star-s_star+1)
    for(t in s_star:e_star){
      a[t-s_star+1] = norm(CUSUM.cov(X, s, e, t), type = "2")
    }
    best_value = max(a)
    best_t = which.max(a) + s_star - 1
    temp1 = BS.cov(X, s, best_t-1, level)
    temp2 = BS.cov(X, best_t, e, level)
    S = c(temp1$S, best_t, temp2$S)
    Dval = c(temp1$Dval, best_value, temp2$Dval)
    Level = c(temp1$Level, level, temp2$Level)
    Parent = cbind(temp1$Parent, parent, temp2$Parent)
    return(list(S = S, Dval = Dval, Level = Level, Parent = Parent))
  }
}