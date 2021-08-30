#' @title Standard binary segmentation for univariate nonparametric change points detection.
#' @description Perform standard binary segmentation for univariate nonparametric change points detection.
#' @param Y         A \code{numeric} matrix of observations with horizontal axis being time, and vertical axis being multiple observations on each time point.
#' @param s         A \code{integer} scalar of starting index.
#' @param e         A \code{integer} scalar of ending index.
#' @param N         A \code{integer} vector representing number of multiple observations on each time point.
#' @param delta     A positive \code{integer} scalar of minimum spacing.
#' @param level     Should be fixed as 0.
#' @param ...       Additional arguments.
#' @return  A \code{list} with the structure:
#' \itemize{
#'  \item S:           A vector of estimated changepoints (sorted in strictly increasing order).
#'  \item Dval:        A vector of values of CUSUM statistic based on KS distance.
#'  \item Level:       A vector representing the levels at which each change point is detected.
#'  \item Parent:      A matrix with the starting indices on the first row and the ending indices on the second row.
#' } 
#' @export
#' @author Oscar Hernan Madrid Padilla, Haotian Xu
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
  if(e-s <= 2*delta){
    return(list(S = S, Dval = Dval, Level = Level, Parent = Parent))
  }else{
    level = level + 1
    parent = matrix(c(s, e), nrow = 2)
    a = rep(0, e-s-2*delta+1)
    for(t in (s+delta):(e-delta)){
      a[t-s-delta+1] = CUSUM.KS(Y, s, e, t, N)
    }
    best_value = max(a)
    best_t = which.max(a) + s + delta - 1
    temp1 = NBS(Y, s, best_t-1, N, delta, level)
    temp2 = NBS(Y, best_t, e, N, delta, level)
    S = c(temp1$S, best_t, temp2$S)
    Dval = c(temp1$Dval, best_value, temp2$Dval)
    Level = c(temp1$Level, level, temp2$Level)
    Parent = cbind(temp1$Parent, parent, temp2$Parent)
    result = list(S = S, Dval = Dval, Level = Level, Parent = Parent)
    class(result) = "BS"
    return(result)
  }
}


#' @export
NBS.tune = function(Y, W, N, len_tau, delta){
  #n = max(N)
  len_time = ncol(Y)
  temp1 = NBS(W, 1, len_time, N, delta, level = 0)  
  Dval = temp1$Dval
  aux = sort(Dval, decreasing = TRUE)
  tau_grid = rev(aux[1:min(len_tau,length(Dval))])
  tau_grid = c(tau_grid, 10)
  B_list = c()
  for(j in 1:length(tau_grid)){
    aux = threshold.BS(temp1, tau_grid[j])$change_points[,1]
    if(length(aux) == 0){
      break;
    }
    B_list[[j]] = sort(aux)
  }
  B_list = unique(B_list)
  O_set = B_list[[1]]
  lambda = log(sum(N))/1.5#2.5#1.5#2#2.555#
  for(m in 1:(length(B_list)-1)){
    eta = min(setdiff(O_set, B_list[[m+1]]))
    B_set_ext = c(1, B_list[[m+1]], len_time)
    k = which.max(B_set_ext < eta)
    eta1 = B_set_ext[k]
    eta2 = B_set_ext[k+1]
    z_hat = as.vector(Y[,eta1:eta2])[which.max(CUSUM.KS(Y, eta1, eta2, eta, N, vector = TRUE))]
    if(error.ECDF(Y, eta1, eta, N, z_hat) + error.ECDF(Y, eta+1, eta2, N, z_hat) - error.ECDF(Y, eta1, eta2, N, z_hat) > lambda){
      O_set = B_list[[m+1]]
    }else{
      break;
    }
  }
  return(O_set)
}



#' @title Internal Function: Compute the CUSUM statistic based on KS distance.
#' @param Y         A \code{numeric} matrix of observations with with horizontal axis being time, and vertical axis being multiple observations on each time point.
#' @param s         A \code{integer} scalar of starting index.
#' @param e         A \code{integer} scalar of ending index.
#' @param t         A \code{integer} scalar of splitting index.
#' @param N         A \code{integer} vector representing number of multiple observations on each time point.
#' @param vector    If TRUE, return a CUSUM vector of empirical distribution function evaluated at the vectorized Y[,s:e]; otherwise, return a scalar representing the maximum of the CUSUM vector.
#' @return  A \code{numeric} scalar of the CUSUM statistic based on KS distance.
#' @noRd
CUSUM.KS = function(Y, s, e, t, N, vector = FALSE){
  n_st = sum(N[s:t])
  n_se = sum(N[s:e])
  n_te = sum(N[(t+1):e])
  aux = as.vector(Y[,s:t])
  aux = aux[which(is.na(aux)==FALSE)]
  temp = ecdf(aux)
  vec_y = as.vector(Y[,s:e])
  vec_y = vec_y[which(is.na(vec_y)==FALSE)]
  Fhat_st = temp(vec_y)# temp(grid)
  aux = as.vector(Y[,(t+1):e])
  aux = aux[which(is.na(aux)==FALSE)]
  temp = ecdf(aux)
  Fhat_te = temp(vec_y)# temp(grid)
  if(vector == TRUE){
    result = sqrt(n_st * n_te / n_se) * abs(Fhat_te - Fhat_st)
  }else{
    result = sqrt(n_st * n_te / n_se) * max(abs(Fhat_te - Fhat_st)) 
  }
  return(result)
}


#' @noRd
error.ECDF = function(Y, s, e, N, z){
  aux = as.vector(Y[,s:e])
  aux = aux[which(is.na(aux)==FALSE)]
  temp = ecdf(aux)
  Fhat = temp(z)
  result = sum(((Y[,s:e] <= z) - Fhat)^2)
  return(result)
}


#' @title Wild binary segmentation for univariate nonparametric change points detection.
#' @description Perform wild binary segmentation for univariate nonparametric change points detection.
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
#'  \item S:           A vector of estimated changepoints (sorted in strictly increasing order).
#'  \item Dval:        A vector of values of CUSUM statistic based on KS distance.
#'  \item Level:       A vector representing the levels at which each change point is detected.
#'  \item Parent:      A matrix with the starting indices on the first row and the ending indices on the second row.
#' } 
#' @export
#' @author Oscar Hernan Madrid Padilla, Haotian Xu
#' @examples
#' Y = t(as.matrix(c(rnorm(100, 0, 1), rnorm(100, 0, 10), rnorm(100, 0, 40))))
#' M = 120
#' intervals = WBS.intervals(M = M, lower = 1, upper = ncol(Y))
#' temp = NWBS(Y, 1, 300, intervals$Alpha, intervals$Beta, N, 5)
#' plot.ts(t(Y))
#' points(x = tail(temp$S[order(temp$Dval)], 4), y = Y[,tail(temp$S[order(temp$Dval)],4)], col = "red")
#' threshold.BS(temp, 1)
NWBS = function(Y, s, e, Alpha, Beta, N, delta, level = 0, ...){ 
  Alpha_new = pmax(Alpha, s)
  Beta_new = pmin(Beta, e)
  idx = which(Beta_new - Alpha_new > 2*delta)
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
      temp = rep(0, Beta_new[m] - Alpha_new[m] - 2*delta + 1)
      for(t in (Alpha_new[m]+delta):(Beta_new[m]-delta)){
        temp[t-(Alpha_new[m]+delta)+1] = CUSUM.KS(Y, Alpha_new[m], Beta_new[m], t, N)
      }
      best_value = max(temp)
      best_t = which.max(temp) + Alpha_new[m] + delta - 1
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
  result = list(S = S, Dval = Dval, Level = Level, Parent = Parent)
  class(result) = "BS"
  return(result)
}



#' @export
NWBS.tune = function(Y, W, Alpha, Beta, N, len_tau = 20, delta){
  #n = max(N)
  len_time = ncol(Y)
  temp1 = NWBS(W, 1, len_time, Alpha, Beta, N, delta, level = 0)  
  Dval = temp1$Dval
  aux = sort(Dval, decreasing = TRUE)
  tau_grid = rev(aux[1:min(len_tau,length(Dval))]-10^{-5})
  tau_grid = c(tau_grid, 10)
  B_list =  c()
  for(j in 1:length(tau_grid)){
    aux = threshold.BS(temp1, tau_grid[j])$change_points[,1]
    if(length(aux) == 0){
      break;
    }
    B_list[[j]] = sort(aux)
  }
  B_list = unique(B_list)
  O_set = B_list[[1]]
  lambda = log(sum(N))/1.5#2.5#1.5#2#2.555#
  for(m in 1:(length(B_list)-1)){
    eta = min(setdiff(O_set, B_list[[m+1]]))
    k = which.max(B_list[[m+1]] < eta)
    eta1 = B_list[[m+1]][k]
    eta2 = B_list[[m+1]][k+1]
    z_hat = as.vector(Y[,eta1:eta2])[which.max(CUSUM.KS(Y, eta1, eta2, eta, N, vector = TRUE))]
    if(error.ECDF(Y, eta1, eta, N, z_hat) + error.ECDF(Y, eta+1, eta2, N, z_hat) - error.ECDF(Y, eta1, eta2, N, z_hat) > lambda){
      O_set = B_list[[m+1]]
    }else{
      break;
    }
  }
  return(O_set)
}
