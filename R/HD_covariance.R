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
#' @description Perform binary Segmentation for covariance change points detection through Operator Norm.
#' @param X         A \code{numeric} matrix of observations with with horizontal axis being time, and vertical axis being dimensions.
#' @param s         A \code{integer} scalar of starting index.
#' @param e         A \code{integer} scalar of ending index.
#' @param level     A parameter for tracking the level at which a change point is detected. Should be fixed as 0.
#' @param ...      Additional arguments.
#' @return  A \code{list} with the structure:
#' \itemize{
#'  \item S:           A vector of estimated changepoints (sorted in strictly increasing order).
#'  \item Dval:        A vector of values of CUSUM statistic based on KS distance.
#'  \item Level:       A vector representing the levels at which each change point is detected.
#'  \item Parent:      A matrix with the starting indices on the first row and the ending indices on the second row.
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
    result = list(S = S, Dval = Dval, Level = Level, Parent = Parent)
    class(result) = "BS"
    return(result)
  }
}



#' @title Internal Function: Principal Component Estimation.
#' @param X         A \code{numeric} matrix of observations with with horizontal axis being time, and vertical axis being dimensions.
#' @param Alpha     A \code{integer} vector of starting indices of random intervals.
#' @param Beta      A \code{integer} vector of ending indices of random intervals.
#' @return  A p-by-M matrix of the Covariance CUSUM statistic.
#' @noRd
PC.cov = function(X, Alpha, Beta){
  p = dim(X)[1]
  n = dim(X)[2]
  delta = 2*p*log(n)+1
  M = length(Alpha)
  u_mat = matrix(0, nrow = p, ncol = M)
  for(m in 1:M){
    if(Beta[m] - Alpha[m] > delta){
      values = sapply(ceiling(Alpha[m]+p*log(n)):floor(Beta[m]-p*log(n)), function(x) norm(CUSUM.cov(X, Alpha[m], Beta[m], x), type = "2"))
      best_t = which.max(values)
      u_mat[,m] = (svd(CUSUM.cov(X, Alpha[m], Beta[m], best_t))$u)[,1]
    }
  }
}


#' @title Wild binary segmentation for covariance change points detection through Independent Projection.
#' @description  Perform wild binary segmentation for covariance change points detection through Independent Projection
#' @param X         A \code{numeric} vector of observations.
#' @param s         A \code{integer} scalar of starting index.
#' @param e         A \code{integer} scalar of ending index.
#' @param Alpha     A \code{integer} vector of starting indices of random intervals.
#' @param Beta      A \code{integer} vector of ending indices of random intervals.
#' @param delta     A positive \code{integer} scalar of minimum spacing.
#' @param level     A parameter for tracking the level at which a change point is detected. Should be fixed as 0.
#' @param ...      Additional arguments.
#' @return  A \code{list} with the structure:
#' \itemize{
#'  \item S:           A vector of estimated changepoints (sorted in strictly increasing order).
#'  \item Dval:        A vector of values of CUSUM statistic based on KS distance.
#'  \item Level:       A vector representing the levels at which each change point is detected.
#'  \item Parent:      A matrix with the starting indices on the first row and the ending indices on the second row.
#' } 
#' @export
#' @author Haotian Xu
#' @examples
#' y = c(rnorm(100, 0, 1), rnorm(100, 0, 10), rnorm(100, 0, 40))
#' M = 120
#' Alpha = sample.int(size = M, n = 300, replace = TRUE)
#' Beta = sample.int(size = M, n = T, replace = TRUE)
#' for(j in 1:M){
#'   aux =  Alpha[j]
#'   aux2 = Beta[j]
#'   Alpha[j] = min(aux, aux2)
#'   Beta[j] = max(aux, aux2)
#' }
#' temp = WBS.univar(y, 1, 300, Alpha, Beta, 5)
#' plot.ts(y)
#' points(x = tail(temp$S[order(temp$Dval)], 4), y = Y[,tail(temp$S[order(temp$Dval)],4)], col = "red")
#' BS.threshold(temp, 1.5)
WBSIP.cov = function(X, X_prime, s, e, Alpha, Beta, delta, level = 0){
  S = NULL
  Dval = NULL
  Level = NULL
  Parent = NULL
  n = dim(X)[2]
  M = length(Alpha)
  u_mat = PC.cov(X_prime, Alpha, Beta)
  y_mat = matrix(NA, nrow = M, ncol = e-s+1)
  for(i in s:e){
    for(m in 1:M){
      y_mat[m, i-s+1] = (t(u_mat[,m]) %*% X[,i])^2
    }
  }
  a = rep(NA, M)
  b = rep(NA, M)
  Alpha_new = ceiling(pmax(Alpha, s) + delta)
  Beta_new = floor(pmin(Beta, e) - delta)
  for(m in 1:M){
    if(Beta_new[m] - Alpha_new[m] >= 2*log(n)){
      s_star = ceiling(Alpha_new[m]+log(n))
      e_star = floor(Beta_new[m]-log(n))
      temp = rep(0, e_star - s_star + 1)
      for(t in s_star:e_star){
        temp[t-s_star+1] = sqrt((t-Alpha_new[m]) * (Beta_new[m]-t) / (Beta_new[m]-Alpha_new[m])) * abs(mean(y_mat[m,(Alpha_new[m]+1):t]) - mean(y_mat[m,(t+1):Beta_new[m]]))
      }
      best_value = max(temp)
      best_t = which.max(temp) + s_star - 1
      a[m] = best_value
      b[m] = best_t
    }else{
      a[m] = -1
    }
  }
  if(all(is.na(b))){
    return(list(S = S, Dval = Dval, Level = Level, Parent = Parent))
  }else{
    level = level + 1
    parent = matrix(c(s, e), nrow = 2)
    m_star = which.max(a)
  }
  temp1 = WBSIP.cov(X, X_prime, s, b[m_star]-1, Alpha, Beta, delta, level)
  temp2 = WBSIP.cov(X, X_prime, b[m_star], e, Alpha, Beta, delta, level)
  S = c(temp1$S, b[m_star], temp2$S)
  Dval = c(temp1$Dval, a[m_star], temp2$Dval)
  Level = c(temp1$Level, level, temp2$Level)
  Parent = cbind(temp1$Parent, parent, temp2$Parent)
  result = list(S = S, Dval = Dval, Level = Level, Parent = Parent)
  class(result) = "BS"
  return(result)
}

