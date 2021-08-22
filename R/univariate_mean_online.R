#' @title Online change point detection via CUSUM (only one change point).
#' @description TO DO
#' @param y_vec       A \code{numeric} vector of observations.
#' @param tau_vec     A \code{numeric} vector of thresholds at time t>= 2.
#' @param ...         Additional arguments.
#' @return  An \code{integer} scalar of estimated change point location.
#' @export
#' @author Haotian Xu
#' @examples
#' TO DO
online.univar.one = function(y_vec, tau_vec, ...){
  if(length(y_vec) - length(tau_vec) != 1){
    stop("tau_vec should be the vector of thresholds at time t >= 2.")
  }
  t = 1
  FLAG = 0
  while(FLAG == 0 & t <= length(y_vec)){
    t = t + 1
    cusum_vec = sapply(1:(t-1), function(s) sqrt((t-s)*s/t) * abs(mean(y_vec[1:s]) - mean(y_vec[(s+1):t])))
    FLAG = 1 - prod(cusum_vec <= tau_vec[t])
  }
  return(t)
}


#' @title Online change point detection via CUSUM 2 (only one change point)
#' @description TO DO
#' @param y_vec       A \code{numeric} vector of observations.
#' @param gamma       A \code{integer} scalar of interval length (>= 2)
#' @param tau_gamma   A \code{numeric} scalar of threshold.
#' @param ...         Additional arguments.
#' @return  An \code{integer} scalar of estimated change point location.
#' @export
#' @author Haotian Xu
#' @examples
#' TO DO
online.univar.one2 = function(y_vec, gamma, tau_gamma, ...){
  t = 1
  FLAG = 0
  while(FLAG == 0 & t <= length(y_vec)){
    t = t + 1
    e = max(t-gamma, 0)
    cusum_vec = sapply((e+1):(t-1), function(s) sqrt((t-s)*(s-e)/(t-e)) * abs(mean(y_vec[(e+1):s]) - mean(y_vec[(s+1):t])))
    FLAG = 1 - prod(cusum_vec <= tau_gamma)
  }
  return(t)
}



#' @title Online change point detection via CUSUM 3 (only one change point).
#' @description TO DO
#' @param y_vec       A \code{numeric} vector of observations.
#' @param tau_vec     A \code{numeric} vector of thresholds at time t>= 1.
#' @param ...         Additional arguments.
#' @return  An \code{integer} scalar of estimated change point location.
#' @export
#' @author Haotian Xu
#' @examples
#' TO DO
online.univar.one3 = function(y_vec, tau_vec, ...){
  if(length(y_vec) != length(tau_vec)){
    stop("y_vec and tau_vec should have the same length.")
  }
  t = 1
  FLAG = 0
  while(FLAG == 0 & t <= length(y_vec)){
    t = t + 1
    J = floor(log2(t))
    j = 0
    while(j < J & FLAG == 0){
      j = j + 1
      s_j = t - 2^(j-1)
      cusum = sqrt((t-s_j)*s_j/t) * abs(mean(y_vec[1:s_j]) - mean(y_vec[(s_j+1):t]))
      FLAG = (cusum > tau_vec[t])
    }
  }
  return(t)
}


online.univar.multi = function(y_vec, tau_vec, ...){
  if(length(y_vec) != length(tau_vec)){
    stop("y_vec and tau_vec should have the same length.")
  }
  t = 1
  FLAG = 0
  while(FLAG == 0 & t <= length(y_vec)){
    t = t + 1
    J = floor(log2(t))
    j = 0
    while(j < J & FLAG == 0){
      j = j + 1
      s_j = t - 2^(j-1)
      cusum = sqrt((t-s_j)*s_j/t) * abs(mean(y_vec[1:s_j]) - mean(y_vec[(s_j+1):t]))
      FLAG = (cusum > tau_vec[t])
    }
  }
  return(t)
}
