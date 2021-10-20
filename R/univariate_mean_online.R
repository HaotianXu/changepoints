#' @title Online change point detection with controlled false alarm rate or average run length.
#' @description Perform online change point detection with controlled false alarm rate or average run length.
#' @param y_vec       A \code{numeric} vector of observations.
#' @param b_vec       A \code{numeric} vector of thresholds b_t with t >= 2.
#' @param train_vec   A \code{numeric} vector of training data from a pre-change distribution, which is only needed to when b_t is NULL in order to calibrate b_t.
#' @param alpha       A \code{numeric} scalar of desired false alarm rate.
#' @param gamma       An \code{integer} scalar of desired average run length.
#' @param permu_num   An \code{integer} scalar of number of random permutation for calibration.
#' @param ...         Additional arguments.
#' @return  A \code{list} with the structure:
#' \itemize{
#'  \item cpt_hat:     An \code{integer} scalar of estimated change point location.
#'  \item b_vec:       A \code{numeric} vector of thresholds b_t with t >= 2.
#' } 
#' @return  
#' @export
#' @author Haotian Xu
#' @examples
#' y_vec = rnorm(300) + c(rep(0, 150), rep(1, 150))
#' train_vec = rnorm(200)
#' temp1 = online.univar.one(y_vec = y_vec, train_vec = train_vec, alpha = 0.05, permu_num = 100)
#' temp1$cpt_hat
#' temp1$b_vec
#' temp2 = online.univar.one(y_vec = y_vec, train_vec = train_vec, gamma = 300, permu_num = 100)
#' temp2$cpt_hat
#' temp2$b_vec
online.univar.one = function(y_vec, b_vec = NULL, train_vec = NULL, alpha = NULL, gamma = NULL, permu_num = NULL, ...){
  if(!is.null(b_vec)){
    if(length(y_vec) - length(b_vec) != 1){
      stop("b_vec should be the vector of thresholds b_t with t >= 2.")
    }
  }else{
    if(is.null(train_vec)){
      stop("Given b_vec is missing, train_vec should be provided to calibrate b_vec.")
    }
    if(is.null(alpha)+is.null(gamma)!=1){
      stop("Given b_vec is missing, one and only one of the parameters alpha and gamma should be provided.")
    }
    if(!is.null(alpha)){
      obs_train = length(train_vec)
      obs_data = length(y_vec)
      if(obs_train > obs_data){
        # only use part of train_vec if it's sample size is bigger than that of y_vec
        train_vec = train_vec[(obs_train-obs_data+1):obs_train]
        obs_train = obs_data
      }
      score = matrix(NA, permu_num, obs_train-1)
      C_vec = rep(NA, permu_num)
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
    }else if(!is.null(gamma)){
      obs_train = length(train_vec)
      if(gamma > obs_train){
        gamma = obs_train
        warning(paste0("gamma is set to be ", obs_train, ". To allow larger value of gamma, please increase the length of train_vec."))
      }
      obs_data = length(y_vec)
      score = matrix(NA, permu_num, obs_train-1)
      C_mat = matrix(NA, permu_num, obs_train-1)
      b = sqrt(log(2^(1/3)*(gamma+1)))
      for(sim in 1:permu_num){
        idx_permu = sample(1:obs_train)
        train_permu = train_vec[idx_permu]
        for(t in 2:obs_train){
          score[sim,t-1] = max(sapply(1:(t-1), function(s) sqrt((t-s)*s/t) * abs(mean(train_permu[1:s]) - mean(train_permu[(s+1):t]))))
        }
        C_mat[sim,] = score[sim,]/b
      }
      C_grid = seq(min(C_mat), max(C_mat), length = 300)
      alarm = rep(0, length(C_grid))
      for(j in 1:length(alarm)){
        aux = apply(C_mat, 1,function(x){min(c(which(x > C_grid[j]), obs_train))})
        alarm[j] = mean(aux)
      }
      ind = which.min(abs(alarm-gamma))
      b_vec = rep(C_grid[ind] * b, obs_data-1)
    }
  }
  t = 1
  FLAG = 0
  while(FLAG == 0 & t <= length(y_vec)){
    t = t + 1
    cusum_vec = sapply(1:(t-1), function(s) sqrt((t-s)*s/t) * abs(mean(y_vec[1:s]) - mean(y_vec[(s+1):t])))
    FLAG = 1 - prod(cusum_vec <= b_vec[t-1])
  }
  return(list(cpt_hat = t, b_vec = b_vec))
}


#' @title Online change point detection with controlled average run length.
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
