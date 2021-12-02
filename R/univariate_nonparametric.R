#' @title Standard binary segmentation for univariate nonparametric change points detection.
#' @description Perform standard binary segmentation for univariate nonparametric change points detection.
#' @param Y         A \code{numeric} matrix of observations with horizontal axis being time, and vertical axis being multiple observations on each time point.
#' @param s         A \code{integer} scalar of starting index.
#' @param e         A \code{integer} scalar of ending index.
#' @param N         A \code{integer} vector representing number of multiple observations on each time point.
#' @param delta     A positive \code{integer} scalar of minimum spacing.
#' @param level     Should be fixed as 0.
#' @return  A \code{list} with the following structure:
#'  \item{S}{A vector of estimated change points (sorted in strictly increasing order)}
#'  \item{Dval}{A vector of values of CUSUM statistic based on KS distance}
#'  \item{Level}{A vector representing the levels at which each change point is detected}
#'  \item{Parent}{A matrix with the starting indices on the first row and the ending indices on the second row}
#' @export
#' @author Oscar Hernan Madrid Padilla & Haotian Xu
#' @examples
#' Y = t(as.matrix(c(rnorm(100, 0, 1), rnorm(100, 0, 10), rnorm(100, 0, 40))))
#' N = rep(1, 300)
#' temp = BS.uni.nonpar(Y, 1, 300, N, 5)
#' plot.ts(t(Y))
#' points(x = tail(temp$S[order(temp$Dval)],4), y = Y[,tail(temp$S[order(temp$Dval)],4)], col = "red")
#' @references Padilla, Yu, Wang and Rinaldo (2021) <doi:10.1214/21-EJS1809>
BS.uni.nonpar = function(Y, s, e, N, delta, level = 0){
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
    temp1 = BS.uni.nonpar(Y, s, best_t-1, N, delta, level)
    temp2 = BS.uni.nonpar(Y, best_t, e, N, delta, level)
    S = c(temp1$S, best_t, temp2$S)
    Dval = c(temp1$Dval, best_value, temp2$Dval)
    Level = c(temp1$Level, level, temp2$Level)
    Parent = cbind(temp1$Parent, parent, temp2$Parent)
    result = list(S = S, Dval = Dval, Level = Level, Parent = Parent)
    class(result) = "BS"
    return(result)
  }
}

#' @title Standard binary segmentation for univariate nonparametric change points detection with tuning parameter selection.
#' @description Perform standard binary segmentation with tuning parameter selection based on sample splitting.
#' @param Y         A \code{numeric} matrix of observations with horizontal axis being time, and vertical axis being multiple observations on each time point.
#' @param W         A copy of the matrix Y, it can be Y itself.
#' @param N         A \code{integer} vector representing number of multiple observations on each time point.
#' @param delta     A positive \code{integer} scalar of minimum spacing.
#' @return  A vector of estimated change points (sorted in strictly increasing order).
#' @export
#' @author Oscar Hernan Madrid Padilla & Haotian Xu
#' @examples
#' Y = t(as.matrix(c(rnorm(100, 0, 1), rnorm(100, 0, 10), rnorm(100, 0, 40))))
#' N = rep(1, 300)
#' temp = BS.uni.nonpar.CPD(Y, Y, N, 5)
#' @references Padilla, Yu, Wang and Rinaldo (2021) <doi:10.1214/21-EJS1809>
BS.uni.nonpar.CPD = function(Y, W, N, delta){
  obs_num = ncol(Y)
  temp1 = BS.uni.nonpar(W, 1, obs_num, N, delta, level = 0)  
  Dval = temp1$Dval
  aux = sort(Dval, decreasing = TRUE)
  len_tau = 30 # number of candidate values considered to tune the threshold
  tau_grid = rev(aux[1:min(len_tau,length(Dval))]) - 10^{-30}
  tau_grid = c(tau_grid, 10)
  B_list = c()
  for(j in 1:length(tau_grid)){
    aux = thresholdBS(temp1, tau_grid[j])$cpt_hat[,1]
    if(length(aux) == 0){
      break
    }
    B_list[[j]] = sort(aux)
  }
  B_list = unique(B_list)
  if(length(B_list) == 0){
    return(NULL)
  }
  if(length(B_list[[1]]) == 0){
    return(B_list[[1]])
  }
  lambda = log(sum(N))/1.5#2.5#1.5#2#2.555#
  for(j in 1:(length(B_list))){
    B2 = B_list[[j]]
    if(j < length(B_list)){
      B1 = B_list[[j+1]]
    }else if(j == length(B_list)){
      B1 = NULL
    }
    temp = setdiff(B2, B1)
    st = -10^15
    for(l in 1:length(temp)){
      eta = temp[l]
      if(length(B1) == 0){
        eta1 = 1
        eta2 = obs_num
      }else if(length(B1) > 0){
        for(k in 1:length(B1)){
          if(B1[k] > eta){
            break
          }
        }
        if(B1[k] > eta){
          eta2 = B1[k]
          if(k == 1)
            eta1 = 1
          if(k > 1)
            eta1 = B1[k-1] + 1
        }
        if(B1[k] < eta){
          eta1 = B1[k] + 1
          eta2 = obs_num
        }
      }######## if length(B1) > 0
      st_aux = CUSUM.KS(Y, eta1, eta2, eta, N)^2
      if(st_aux > st){
        st = st_aux
      }
    }###  close for defining  eta1 and eta2
    if(st > lambda){
      return(B2)
    }
  }
  return(B1)
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
#' @return  A \code{list} with the following structure:
#'  \item{S}{A vector of estimated change points (sorted in strictly increasing order)}
#'  \item{Dval}{A vector of values of CUSUM statistic based on KS distance}
#'  \item{Level}{A vector representing the levels at which each change point is detected}
#'  \item{Parent}{A matrix with the starting indices on the first row and the ending indices on the second row}
#' @export
#' @author Oscar Hernan Madrid Padilla, Haotian Xu
#' @examples
#' Y = t(as.matrix(c(rnorm(100, 0, 1), rnorm(100, 0, 10), rnorm(100, 0, 40))))
#' N = rep(1, 300)
#' M = 20
#' intervals = WBS.intervals(M = M, lower = 1, upper = ncol(Y))
#' temp = WBS.uni.nonpar(Y, 1, 300, intervals$Alpha, intervals$Beta, N, 5)
#' plot.ts(t(Y))
#' points(x = tail(temp$S[order(temp$Dval)], 4), y = Y[,tail(temp$S[order(temp$Dval)],4)], col = "red")
#' thresholdBS(temp, 2)
#' @references Padilla, Yu, Wang and Rinaldo (2021) <doi:10.1214/21-EJS1809>
WBS.uni.nonpar = function(Y, s, e, Alpha, Beta, N, delta, level = 0){ 
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
  temp1 = WBS.uni.nonpar(Y, s, b[m_star]-1, Alpha, Beta, N, delta, level)
  temp2 = WBS.uni.nonpar(Y, b[m_star], e, Alpha, Beta, N, delta, level)
  S = c(temp1$S, b[m_star], temp2$S)
  Dval = c(temp1$Dval, a[m_star], temp2$Dval)
  Level = c(temp1$Level, level, temp2$Level)
  Parent = cbind(temp1$Parent, parent, temp2$Parent)
  result = list(S = S, Dval = Dval, Level = Level, Parent = Parent)
  class(result) = "BS"
  return(result)
}



#' @title Wild binary segmentation for univariate nonparametric change points detection with tuning parameter selection.
#' @description Perform wild binary segmentation with tuning parameter selection based on sample splitting.
#' @param Y         A \code{numeric} matrix of observations with horizontal axis being time, and vertical axis being multiple observations on each time point.
#' @param W         A copy of the matrix Y, it can be Y itself.
#' @param Alpha     A \code{integer} vector of starting indices of random intervals.
#' @param Beta      A \code{integer} vector of ending indices of random intervals.
#' @param N         A \code{integer} vector representing number of multiple observations on each time point.
#' @param delta     A positive \code{integer} scalar of minimum spacing.
#' @return  A vector of estimated change points (sorted in strictly increasing order).
#' @export
#' @author Oscar Hernan Madrid Padilla & Haotian Xu
#' @examples
#' Y = t(as.matrix(c(rnorm(100, 0, 1), rnorm(100, 0, 10), rnorm(100, 0, 40))))
#' N = rep(1, 300)
#' M = 20
#' intervals = WBS.intervals(M = M, lower = 1, upper = ncol(Y))
#' temp = WBS.uni.nonpar.CPD(Y, Y, intervals$Alpha, intervals$Beta, N, 5)
#' @references Padilla, Yu, Wang and Rinaldo (2021) <doi:10.1214/21-EJS1809>
WBS.uni.nonpar.CPD = function(Y, W, Alpha, Beta, N, delta){
  obs_num = ncol(Y)
  temp1 = WBS.uni.nonpar(W, 1, obs_num, Alpha, Beta, N, delta, level = 0)  
  Dval = temp1$Dval
  aux = sort(Dval, decreasing = TRUE)
  len_tau = 30
  tau_grid = rev(aux[1:min(len_tau,length(Dval))]) - 10^{-30}
  tau_grid = c(tau_grid, 10)
  B_list = c()
  for(j in 1:length(tau_grid)){
    aux = thresholdBS(temp1, tau_grid[j])$cpt_hat[,1]
    if(length(aux) == 0){
      break
    }
    B_list[[j]] = sort(aux)
  }
  B_list = unique(B_list)
  if(length(B_list) == 0){
    return(NULL)
  }
  if(length(B_list[[1]]) == 0){
    return(B_list[[1]])
  }
  lambda = log(sum(N))/1.5#2.5#1.5#2#2.555#
  for(j in 1:(length(B_list))){
    B2 = B_list[[j]]
    if(j < length(B_list)){
      B1 = B_list[[j+1]]
    }else if(j == length(B_list)){
      B1 = NULL
    }
    temp = setdiff(B2, B1)
    st = -10^15
    for(l in 1:length(temp)){
      eta = temp[l]
      if(length(B1) == 0){
        eta1 = 1
        eta2 = obs_num
      }else if(length(B1) > 0){
        for(k in 1:length(B1)){
          if(B1[k] > eta){
            break
          }
        }
        if(B1[k] > eta){
          eta2 = B1[k]
          if(k == 1)
            eta1 = 1
          if(k > 1)
            eta1 = B1[k-1] + 1
        }
        if(B1[k] < eta){
          eta1 = B1[k] + 1
          eta2 = obs_num
        }
      }######## if length(B1) > 0
      st_aux = CUSUM.KS(Y, eta1, eta2, eta, N)^2
      if(st_aux > st){
        st = st_aux
      }
    }###  close for defining  eta1 and eta2
    if(st > lambda){
      return(B2)
    }
  }
  return(B1)
}

