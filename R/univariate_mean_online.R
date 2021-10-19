#' @title Online change point detection via CUSUM (single change point).
#' @description Perform online change point detection via CUSUM (single change point).
#' @param y_vec       A \code{numeric} vector of observations.
#' @param b_vec       A \code{numeric} vector of thresholds b_t with t >= 2.
#' @param train_vec   A \code{numeric} vector of training data from a pre-change distribution, which is only needed to when b_t is NULL in order to calibrate b_t.
#' @param ...         Additional arguments.
#' @return  An \code{integer} scalar of estimated change point location.
#' @export
#' @author Haotian Xu
#' @examples
#' TO DO
online.univar.one = function(y_vec, b_vec, train_vec = NULL, alpha = 0.05, permu_num = 100, ...){
  if(!is.null(b_vec)){
    if(length(y_vec) - length(b_vec) != 1){
      stop("b_vec should be the vector of thresholds b_t with t >= 2.")
    }
  }else{
    obs_train = length(train_vec)
    obs_data = length(y_vec)
    if(obs_train > obs_data){
      # only use part of train_vec if it's sample size is bigger than that of y_vec
      train_vec = train_vec[(obs_train-obs_data+1):obs_train]
      obs_train = obs_data
    }
    score = matrix(NA, permu_num, obs_train-1)
    C_vec = matrix(NA, permu_num)
    trend = sapply(2:obs_train, function(t) sqrt(log(t/alpha)))
    for(sim in 1:permu_num){
      idx_permu = sample(1:obs_train)
      train_permu = train_vec[idx_permu]
      for(t in 2:obs_train){
        score[sim,t-1] = max(sapply(1:(t-1), function(s) sqrt((t-s)*s/t) * abs(mean(train_permu[1:s]) - mean(train_permu[(s+1):t]))))
      }
      C_vec[sim] = max(score[sim,]/trend)
    }
    b_vec = quantile(C_vec, 1-alpha) * sapply(2:obs_data, function(t) sqrt(log(t/alpha)))
  }
  t = 1
  FLAG = 0
  while(FLAG == 0 & t <= length(y_vec)){
    t = t + 1
    cusum_vec = sapply(1:(t-1), function(s) sqrt((t-s)*s/t) * abs(mean(y_vec[1:s]) - mean(y_vec[(s+1):t])))
    FLAG = 1 - prod(cusum_vec <= b_vec[t-1])
  }
  return(t)
}


#' @title Online change point detection via CUSUM (single change point, type 2).
#' @description  Perform online change point detection via CUSUM (single change point, type 2).
#' @param y_vec       A \code{numeric} vector of observations.
#' @param gamma       A \code{integer} scalar of interval length (>= 2).
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



#' @title Online change point detection via CUSUM (single change point, type 3).
#' @description Perform online change point detection via CUSUM (single change point, type 3).
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


#' @title Online change point detection via CUSUM (multiple change points).
#' @description Perform Online change point detection via CUSUM (multiple change points).
#' @param y_vec       A \code{numeric} vector of observations.
#' @param tau_vec     A \code{numeric} vector of thresholds at time t>= 2.
#' @param ...         Additional arguments.
#' @return  An \code{integer} scalar of estimated change point location.
#' @export
#' @author Haotian Xu
#' @examples
#' TO DO
online.univar.multi = function(y_vec, tau_vec, ...){
  if(length(y_vec) - length(tau_vec) != 1){
    stop("tau_vec should be the vector of thresholds at time t >= 2.")
  }
  cpt = NULL
  e = 0
  t = 1
  FLAG = 0
  while(t < length(y_vec)-1){
    t = t + 1
    cusum_vec = sapply((e+1):(t-1), function(s) sqrt((t-s)*(s-e)/(t-e)) * abs(mean(y_vec[(e+1):s]) - mean(y_vec[(s+1):t])))
    FLAG = 1 - prod(cusum_vec <= tau_vec[t])
    if(FLAG == 1){
      cpt = c(cpt, t)
      FLAG = 0
      e = t
    }
  }
  return(cpt)
}
