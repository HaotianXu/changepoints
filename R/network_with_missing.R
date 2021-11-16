#' @export
softImpute.network.missing <- function(data_incomplete_list, eta_list, lambda, rho_star, it_max = 10000) {
  .Call('_changepoints_rcpp_soft_impute', PACKAGE = 'changepoints', data_incomplete_list, eta_list, lambda, rho, it_max)
}

#' @export
lambda.network.missing <- function(s, e, t, alpha, rho, m, p, C_lambda) {
  .Call('_changepoints_rcpp_lambda', PACKAGE = 'changepoints', s, e, t, alpha, rho, m, p, C_lambda)
}

#' @export
CUSUM.network.missing <- function(data_incomplete_list, eta_list, s, e, t, alpha, rho, m, C_lambda, delta) {
  .Call('_changepoints_rcpp_CUSUM', PACKAGE = 'changepoints', data_incomplete_list, eta_list, s, e, t, alpha, rho, m, C_lambda, delta)
}

#' @export
threshold.network.missing = function(s, t, rank, pi_lb, p, rho, pi_ub, alpha){
  return(sqrt(rank*rho*p*pi_ub*log(s/alpha)/(pi_lb^2*s)) + sqrt(rank*rho*p*pi_ub*log(t/alpha)/(pi_lb^2*(t-s))))
}


#' @title Online change point detection for network data with missing values.
#' @description  Perform online change point detection for network data by controlling the false alarm rate at level alpha or controlling the average run length gamma. The default choice of the tuning parameters tau1, tau2 and tau3 are used (see Section 4.1 of the reference).
#' @param data_mat1  A \code{numeric} matrix of observations with with horizontal axis being time, and with each column be the vectorized adjacency matrix.
#' @param data_mat2  A \code{numeric} matrix of observations with with horizontal axis being time, and with each column be the vectorized adjacency matrix (data_mat1 and data_mat2 are independent and have the same dimensions ).
#' @param b_vec      A \code{numeric} vector of thresholds b_t with t >= 2.
#' @param train_mat  A \code{numeric} matrix of training data from a pre-change distribution(no change point), which is only needed to when b_vec is NULL in order to calibrate b_t.
#' @param alpha      A \code{numeric} scalar in (0,1) representing the level.
#' @param gamma      An \code{integer} scalar of desired average run length.
#' @param permu_num  An \code{integer} scalar of number of random permutation for calibration.
#' @param ...        Additional arguments.
#' @return  A \code{list} with the following structure:
#'  \item{cpt}{Estimated change point}
#'  \item{score}{A \code{numeric} vector of computed cumsum statistics}
#'  \item{b_vec}{A \code{numeric} vector of thresholds b_t with t >= 2}
#' @export
#' @author  Haotian Xu
#' @references Yu Y, Padilla, O, Wang D, Rinaldo A. Optimal network online change point localisation. arXiv preprint arXiv:2101.05477.
#' @examples
#' p = 12 # number of nodes
#' rho = 0.5 # sparsity parameter
#' block_num = 3 # number of groups for SBM
#' train_obs_num = 200 # sample size for each segment
#' conn1_mat = rho * matrix(c(0.6,1,0.6,1,0.6,0.5,0.6,0.5,0.6), nrow = 3) # connectivity matrix 
#' set.seed(1)
#' can_vec = sample(1:p, replace = F) # randomly assign nodes into groups
#' sbm = simu.SBM(conn1_mat, can_vec, train_obs_num, symm = TRUE, self = TRUE)
#' train_mat = sbm$obs_mat
#' train_list = lapply(1:ncol(train_mat), function(t) lowertri2mat(train_mat[,t], p, diag = TRUE))
#' pi_mat = matrix(0.9, p, p)
#' train_eta_list = lapply(1:length(train_list), function(t) gen.missing(pi_mat, symm = TRUE))
#' train_miss_list = lapply(1:length(train_list), function(t) train_eta_list[[t]] * train_list[[t]])
#' pi_lb_hat = quantile(Reduce("+", train_eta_list)/train_obs_num, 0.05) # estimator of pi_lb
#' pi_ub_hat = quantile(Reduce("+", train_eta_list)/train_obs_num, 0.95) # estimator of pi_ub
#' C_lambda = 2/3
#' graphon_miss_impute = softImpute.network.missing(train_miss_list, train_eta_list, lambda.network.missing(1, length(train_miss_list), length(train_miss_list), 0.05, rho = 0.509, m = pi_ub_hat, p, C_lambda), 1)
#' graphon_miss_hat = graphon_miss_impute$u %*% diag(as.numeric(graphon_miss_impute$d)) %*% t(graphon_miss_impute$v)
#' rho_hat = quantile(graphon_miss_hat, 0.95)
#' rank_hat = sum(graphon_miss_impute$d != 0)
#' alpha_grid = c(0.05, 0.01)
#' permu_num = 100
#' threshold_len = 300
#' temp = calibrate.network.missing(train_miss_list, train_eta_list, threshold_len, alpha_grid, permu_num, pi_lb_hat, pi_ub_hat, rho_hat, rank_hat, C_lambda, delta = 5)
calibrate.network.missing = function(train_miss_list, train_eta_list, threshold_len, alpha_grid, permu_num, pi_lb_hat, pi_ub_hat, rho_hat, rank_hat, C_lambda, delta = 5, ...){
  burnin_idx = ceiling(log2(2*delta))
  train_obs_num = length(train_miss_list)
  p = ncol(train_miss_list[[1]])
  if(train_obs_num != length(train_eta_list)){
    stop("train_miss_list and train_eta_list should have the same length.")
  }
  if(any(alpha_grid <= 0) | any(alpha_grid >= 1)){
    stop("All the elements of alpha_grid should be in (0,1).")
  }
  if(train_obs_num <= max(4*delta+1, 100)){
    stop("The length of train_miss_list is too short, please increase.")
  }
  scores_si = array(0, c(train_obs_num, floor(log2(train_obs_num))-burnin_idx+2, length(alpha_grid), permu_num))
  print("Start the calibration step:")
  pb = txtProgressBar(min = 0, max = permu_num, style = 3)
  counter = 0
  for(sim in 1:permu_num){
    idx_permu = sample(1:train_obs_num)
    for(t in 2:train_obs_num){
      if(t >= 4*delta+1){
        m = floor(log2(t)) - 1
        if(t - 2^(m+1) >= 2*delta){
          m = m + 1
        }
        if(m < burnin_idx){
          N_grid = c(2*delta+1)
        }else{
          N_grid = c(2*delta+1, 2^(burnin_idx:m))
        }
        for(ind_alpha in 1:length(alpha_grid)){
          alpha = alpha_grid[ind_alpha]
          scores_si[t, 1:length(N_grid), ind_alpha, sim] = sapply(N_grid, function(s) CUSUM.network.missing(train_miss_list[idx_permu], train_eta_list[idx_permu], 1, t, s, alpha, rho_hat, pi_ub_hat, C_lambda, delta = delta))
        }
      }
    }
    counter = counter + 1
    setTxtProgressBar(pb, counter)
  }
  print("Finish the calibration.")
  thresholds = array(1, c(train_obs_num, floor(log2(train_obs_num))-burnin_idx+2, length(alpha_grid)))
  for(t in 2:train_obs_num){
    if(t >= 4*delta+1){
      m = floor(log2(t)) - 1
      if(t - 2^(m+1) >= 2*delta){
        m = m + 1
      }
      if(m < burnin_idx){
        N_grid = c(2*delta+1)
      }else{
        N_grid = c(2*delta+1, 2^(burnin_idx:m))
      }
      for(ind_alpha in 1:length(alpha_grid)){
        alpha = alpha_grid[ind_alpha]
        thresholds[t,1:length(N_grid),ind_alpha] = sapply(N_grid, function(s) threshold.network.missing(s, t, rank_hat, pi_lb_hat, p, rho_hat, pi_ub_hat, alpha))
      }
    }
  }
  C_mat = matrix(0, length(alpha_grid), permu_num)
  for(ind_alpha in 1:length(alpha_grid)){
    for(sim in 1:permu_num){
      C_mat[ind_alpha, sim] = max(scores_si[,,ind_alpha,sim]/thresholds[,,ind_alpha])
    }
  }
  C_vec = rep(0, length(alpha_grid))
  for(ind_alpha in 1:length(alpha_grid)){
    C_vec[ind_alpha] = quantile(C_mat[ind_alpha,], 1-alpha_grid[ind_alpha])
  }
  thresholds_si_array = array(1, c(threshold_len, floor(log2(threshold_len))-burnin_idx+2, length(alpha_grid)))
  for(t in 2:threshold_len){
    if(t >= 4*delta+1){
      m = floor(log2(t)) - 1
      if(t - 2^(m+1) >= 2*delta){
        m = m + 1
      }
      if(m < burnin_idx){
        N_grid = c(2*delta+1)
      }else{
        N_grid = c(2*delta+1, 2^(burnin_idx:m))
      }
      for(ind_alpha in 1:length(alpha_grid)){
        alpha = alpha_grid[ind_alpha]
        thresholds_si_array[t,1:length(N_grid),ind_alpha] = C_vec[ind_alpha] * sapply(N_grid, function(s) threshold.network.missing(s, t, rank_hat, pi_lb_hat, p, rho_hat, pi_ub_hat, alpha))
      }
    }
  }
  return(list(C_lambda = C_lambda, rho_hat = rho_hat, rank_hat = rank_hat, pi_lb_hat = pi_lb_hat, pi_ub_hat = pi_ub_hat, thresholds_si_array = thresholds_si_array))
}
