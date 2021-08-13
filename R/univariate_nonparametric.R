#' @title Standard binary segmentation for univariate nonparametric change points detection
#' @description TO DO
#' @param Y         A \code{numeric} matrix of observations with horizontal axis being time, and vertical axis being multiple observations on each time point.
#' @param s         A \code{integer} scalar of starting index.
#' @param e         A \code{integer} scalar of ending index.
#' @param N         A \code{integer} vector representing number of multiple observations on each time point.
#' @param delta     A positive \code{integer} scalar of minimum spacing.
#' @param level     Should be fixed as 0.
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
#' Y = t(as.matrix(c(rnorm(100, 0, 1), rnorm(100, 0, 10), rnorm(100, 0, 40))))
#' N = rep(1, 300)
#' temp = NBS(Y, 1, 300, N, 5)
#' plot.ts(t(Y))
#' points(x = tail(temp$S[order(temp$Dval)],4), y = Y[,tail(temp$S[order(temp$Dval)],4)], col = "red")
NBS = function(Y, s, e, N, delta, level = 0, ...){
  S = NULL
  Dval = NULL
  Level = NULL
  Parent = NULL
  if(e-s <= delta){
    return(list(S = S, Dval = Dval, Level = Level, Parent = Parent))
  }else{
    level = level + 1
    parent = matrix(c(s, e), nrow = 2)
    a = rep(0, e-s-1)
    for(t in (s+1):(e-1)){
      a[t-s] = CUSUM.KS(Y, s, e, t, N)
    }
    best_value = max(a)
    best_t = which.max(a) + s
    temp1 = NBS(Y, s, best_t-1, N, delta, level)
    temp2 = NBS(Y, best_t, e, N, delta, level)
    S = c(temp1$S, best_t, temp2$S)
    Dval = c(temp1$Dval, best_value, temp2$Dval)
    Level = c(temp1$Level, level, temp2$Level)
    Parent = cbind(temp1$Parent, parent, temp2$Parent)
    return(list(S = S, Dval = Dval, Level = Level, Parent = Parent))
  }
}



#' @title Internal Function: Compute value of CUSUM statistic based on KS distance.
#' @param Y         A \code{numeric} matrix of observations with with horizontal axis being time, and vertical axis being multiple observations on each time point.
#' @param s         A \code{integer} scalar of starting index.
#' @param e         A \code{integer} scalar of ending index.
#' @param t         A \code{integer} scalar of splitting index.
#' @param N         A \code{integer} vector representing number of multiple observations on each time point.
#' @return  A \code{numeric} scalar of value of CUSUM statistic based on KS distance.
#' @noRd
CUSUM.KS = function(Y, s, e, t, N){
  n_st = sum(N[(s+1):t])
  n_se = sum(N[(s+1):e])
  n_te = sum(N[(t+1):e])
  aux = Y[, (s+1):t]
  aux = aux[which(is.na(aux)==FALSE)]
  temp = ecdf(aux)
  vec_y = Y[, s:e]
  vec_y = vec_y[which(is.na(vec_y)==FALSE)]
  Fhat_st = temp(vec_y)# temp(grid)
  aux = Y[, (t+1):e]
  aux = aux[which(is.na(aux)==FALSE)]
  temp = ecdf(aux)
  Fhat_te = temp(vec_y)# temp(grid)
  result = sqrt(n_st * n_te / n_se) * max(abs(Fhat_te - Fhat_st))
  return(result)
}



#' @title Wild binary segmentation for univariate nonparametric change points detection
#' @description TO DO
#' @param Y         A \code{numeric} matrix of observations with horizontal axis being time, and vertical axis being multiple observations on each time point.
#' @param s         A \code{integer} scalar of starting index.
#' @param e         A \code{integer} scalar of ending index.
#' @param Alpha     A \code{integer} vector of starting indices of random intervals.
#' @param Beta      A \code{integer} vector of ending indices of random intervals.
#' @param N         A \code{integer} vector representing number of multiple observations on each time point.
#' @param delta     A positive \code{integer} scalar of minimum spacing.
#' @param level     Should be fixed as 0.
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
#' Y = t(as.matrix(c(rnorm(100, 0, 1), rnorm(100, 0, 10), rnorm(100, 0, 40))))
#' M =   120
#' Alpha = sample.int(size = M, n = 300, replace = TRUE)
#' Beta = sample.int(size = M, n = T, replace = TRUE)
#' for(j in 1:M){
#'   aux =  Alpha[j]
#'   aux2 = Beta[j]
#'   Alpha[j] = min(aux, aux2)
#'   Beta[j] = max(aux, aux2)
#' }
#' temp = NWBS(Y, 1, 300, Alpha, Beta, N, 5)
#' plot.ts(t(Y))
#' points(x = tail(temp$S[order(temp$Dval)], 4), y = Y[,tail(temp$S[order(temp$Dval)],4)], col = "red")
NWBS = function(Y, s, e, Alpha, Beta, N, delta, level = 0){ 
  Alpha_new = pmax(Alpha, s)
  Beta_new = pmin(Beta, e)
  idx = which(Beta_new - Alpha_new > delta)
  Alpha_new = Alpha_new[idx]
  Beta_new = Beta_new[idx]
  M = length(Alpha_new)
  # xi = 1/8
  # Alpha_new2 = Alpha_new
  # Beta_new2  = Beta_new
  # Alpha_new = ceiling((1-xi)*Alpha_new2+  xi*Beta_new2)
  # Beta_new = ceiling((1-xi)*Beta_new2 +  xi*Alpha_new2)
  # idx = which(Beta_new - Alpha_new > delta)
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
      temp = rep(0, Beta_new[m] - Alpha_new[m] - 1)
      for(t in (Alpha_new[m]+1):(Beta_new[m]-1)){
        temp[t-(Alpha_new[m])] = CUSUM.KS(Y, Alpha_new[m], Beta_new[m], t, N)
      }
      best_value = max(temp)
      best_t = which.max(temp) + Alpha_new[m]
      a[m] = best_value
      b[m] = best_t
    }
    m_star = which.max(a)
  }
  temp1 = NWBS(Y, s, b[m_star]-1, Alpha, Beta, N, delta, level)
  temp2 = NWBS(Y, b[m_star], e, Alpha, Beta, N, delta, level)
  S = c(temp1$S, b[m_star], temp2$S)
  Dval = c(temp1$Dval, a[m_star], temp2$Dval)
  Level = c(temp1$Level, level, temp2$Level)
  Parent = cbind(temp1$Parent, parent, temp2$Parent)
  return(list(S = S, Dval = Dval, Level = Level, Parent = Parent))
}


