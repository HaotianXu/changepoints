#' @title Internal Function: Compute the CUSUM statistic based on KS distance (multivariate).
#' @param Y         A \code{numeric} matrix of observations with with horizontal axis being time, and vertical axis being dimension.
#' @param s         A \code{integer} scalar of starting index.
#' @param e         A \code{integer} scalar of ending index.
#' @param t         A \code{integer} scalar of splitting index.
#' @param h         A \code{integer} scalar of bandwidth.
#' @return  A \code{numeric} scalar of the CUSUM statistic based on KS distance.
#' @noRd
CUSUM.KS.multivariate = function(Y, s, e, t, h){
  p = dim(Y)[1]
  n = dim(Y)[2]
  n_st = t - s
  n_se = e - s
  n_te = e - t
  aux = Y[,(s+1):t]
  aux = aux[which(is.na(aux)==FALSE)]
  temp1 = kde(t(aux), gridsize = 30, eval.points = t(Y), H = h*diag(p))
  aux = Y[,(t+1):e]
  aux = aux[which(is.na(aux)==FALSE)]
  temp2 = kde(t(aux), gridsize = 30, eval.points = t(Y), H = h*diag(p))
  result = sqrt(n_st * n_te / n_se) * max(abs(temp1$estimate - temp2$estimate))
  return(result)
}


#' @title Wild binary segmentation for multivariate nonparametric change points detection
#' @description TO DO
#' @param Y         A \code{numeric} matrix of observations with horizontal axis being time, and vertical axis being multiple observations on each time point.
#' @param s         A \code{integer} scalar of starting index.
#' @param e         A \code{integer} scalar of ending index.
#' @param Alpha     A \code{integer} vector of starting indices of random intervals.
#' @param Beta      A \code{integer} vector of ending indices of random intervals.
#' @param h         A \code{numeric} vector bandwidth parameter
#' @param level     Should be fixed as 0.
#' @param ...       Additional arguments.
#' @return     A \code{list} with the structure:
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
#' new_MWBS
#' 
MNWBS = function(Y, s, e, Alpha, Beta, h, level = 0, ...){
  Alpha_new = pmax(Alpha, s)
  Beta_new = pmin(Beta, e)
  idx = which(Beta_new - Alpha_new > 2 * max(p,h^(-p)))
  Alpha_new = Alpha_new[idx]
  Beta_new = Beta_new[idx]
  M = length(Alpha_new)
  # xi = 1/8
  # Alpha_new2 = Alpha_new
  # Beta_new2  = Beta_new
  # Alpha_new = ceiling((1-xi)*Alpha_new2+  xi*Beta_new2)
  # Beta_new = ceiling((1-xi)*Beta_new2 +  xi*Alpha_new2)
  # idx = which(Beta_new - Alpha_new > 1)
  # Alpha_new = Alpha_new[idx]
  # Beta_new = Beta_new[idx]
  # M = length(Alpha_new)
  S = NULL
  Dval = NULL
  Level = NULL
  Parent = NULL
  if(M == 0){
    return(list(S = S, Dval = Dval, Level = Level, Parent = Parent))
  }else{
    level = level + 1
    parent = matrix(c(s, e), nrow = 2)
    a = rep(0, M)
    b = rep(0, M)
    for(m in 1:M){
      temp = rep(0, Beta_new[m] - Alpha_new[m] - 2*max(h^{-p},p) + 1)
      for(t in (Alpha_new[m]+max(h^{-p},p)):(Beta_new[m]-max(h^{-p},p))){
        temp[t-(Alpha_new[m])] = CUSUM.KS.multivariate(Y, Alpha_new[m], Beta_new[m], t, h)
      }
      best_value = max(temp)
      best_t = which.max(temp) + Alpha_new[m]
      a[m] = best_value
      b[m] = best_t
    }
    m_star = which.max(a)
  }
  temp1 = MNWBS(Y, s, b[m_star]-1, Alpha, Beta, h, level)
  temp2 = MNWBS(Y, b[m_star], e, Alpha, Beta, h, level)
  S = c(temp1$S, b[m_star], temp2$S)
  Dval = c(temp1$Dval, a[m_star], temp2$Dval)
  Level = c(temp1$Level, level, temp2$Level)
  Parent = cbind(temp1$Parent, parent, temp2$Parent)
  return(list(S = S, Dval = Dval, Level = Level, Parent = Parent))
}