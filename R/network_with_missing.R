#' @title Internal function: Estimate graphon matrix by soft-impute for independent adjacency matrices with missing values
#' @param data_incomplete_list  A \code{list} of adjacency matrices (with entries being 0 or 1) with missing values being coercing to 0.
#' @param eta_list              A \code{list} of matrices associated with data_incomplete_list, each matrix indicates the missing entries in corresponding adjacency matrix.
#' @param lambda                A \code{numeric} scalar of thresholding parameter for leading singular value in the soft-impute algorithm.
#' @param a                     A \code{numeric} scalar of truncation parameter in the soft-impute algorithm.
#' @param it_max                An \code{integer} scalar of maximum iteration for the soft-impute algorithm.
#' @return  Estimated graphon matrix
#' @export
#' @author  Haotian Xu
#' @references Dubey, Xu and Yu (2021) <arxiv:2110.06450>
softImpute.network.missing <- function(data_incomplete_list, eta_list, lambda, a, it_max = 10000) {
  .Call('_changepoints_rcpp_soft_impute', PACKAGE = 'changepoints', data_incomplete_list, eta_list, lambda, a, it_max)
}


#' @title Function to compute the thresholding parameter for leading singular value in the soft-impute algorithm (see Theorem 2 in the reference)
#' @param s          An \code{integer} scalar of the starting index.
#' @param e          An \code{integer} scalar of the ending index.
#' @param t          An \code{integer} scalar of the splitting index.
#' @param alpha      A \code{numeric} scalar in (0,1) representing the desired false alarm rate.
#' @param rho        A \code{numeric} scalar of the sparsity parameter.
#' @param pi_ub      A \code{numeric} scalar of the upper bound of the missing probability.
#' @param p          An \code{integer} scalar of the dimensionality of the graphon matrix.
#' @param C_lambda   A \code{numeric} scalar of an absolute constant, which is set to be 2/3 by default.
#' @export
#' @return  The default thresholding parameter for leading singular value in the soft-impute algorithm
#' @references Dubey, Xu and Yu (2021) <arxiv:2110.06450>
lambda.network.missing <- function(s, e, t, alpha, rho, pi_ub, p, C_lambda) {
  .Call('_changepoints_rcpp_lambda', PACKAGE = 'changepoints', s, e, t, alpha, rho, pi_ub, p, C_lambda)
}

#' @title Internal function: CUSUM statistic based on soft-imput estimators
#' @noRd
#' @references Dubey, Xu and Yu (2021) <arxiv:2110.06450>
CUSUM.network.missing <- function(data_incomplete_list, eta_list, s, e, t, alpha, rho, m, C_lambda, delta) {
  .Call('_changepoints_rcpp_CUSUM', PACKAGE = 'changepoints', data_incomplete_list, eta_list, s, e, t, alpha, rho, m, C_lambda, delta)
}

#' @title Internal function: Function to compute the threshold for online change point detection (see Theorem 2 in the reference)
#' @noRd
#' @references Dubey, Xu and Yu (2021) <arxiv:2110.06450>
threshold.network.missing = function(s, t, rank, pi_lb, p, rho, pi_ub, alpha){
  return(sqrt(rank*rho*p*pi_ub*log(s/alpha)/(pi_lb^2*s)) + sqrt(rank*rho*p*pi_ub*log(t/alpha)/(pi_lb^2*(t-s))))
}


#' @title Calibrate step for online change point detection for network data with missing values.
#' @description  Calibrate step for online change point detection for network data by controlling the false alarm rate at level alpha.
#' @param train_miss_list  A \code{list} of adjacency matrices (with entries being 0 or 1) with missing values being coercing to 0.
#' @param train_eta_list   A \code{list} of matrices associated with data_incomplete_list, each matrix indicates the missing entries in corresponding adjacency matrix.
#' @param threshold_len    An \code{integer} scalar of the length of tuned thresholds.
#' @param alpha_grid       A \code{numeric} vector in (0,1) representing the desired false alarm rate.
#' @param permu_num        An \code{integer} scalar of number of random permutation for calibration.
#' @param pi_lb_hat        A \code{numeric} scalar of the lower bound of the missing probability.
#' @param pi_ub_hat        A \code{numeric} scalar of the upper bound of the missing probability.
#' @param rho_hat          A \code{numeric} scalar of the sparsity parameter.
#' @param rank_hat         An \code{integer} scalar of the rank of the underlying graphon matrix.
#' @param C_lambda         A \code{numeric} scalar of an absolute constant, which is set to be 2/3 by default.
#' @param delta            An \code{integer} scalar of minimum spacing.
#' @return  A \code{list} with the following structure:
#'  \item{C_lambda}{The absolute constant}
#'  \item{rho_hat}{the (estimated) sparsity parameter}
#'  \item{rank_hat}{the (estimated) rank of underlying graphon matrix}
#'  \item{pi_lb_hat}{the (estimated) lower bound of the missing probability}
#'  \item{pi_ub_hat}{the (estimated) upper bound of the missing probability}
#'  \item{thresholds_array}{A \code{numeric} array of calibrated threshold}
#' @export
#' @author  Haotian Xu
#' @references Dubey, Xu and Yu (2021) <arxiv:2110.06450>
#' @examples
#' p = 6 # number of nodes
#' rho = 0.5 # sparsity parameter
#' block_num = 3 # number of groups for SBM
#' train_obs_num = 150 # sample size for each segment
#' conn1_mat = rho * matrix(c(0.6,1,0.6,1,0.6,0.5,0.6,0.5,0.6), nrow = 3) # connectivity matrix 
#' set.seed(1)
#' can_vec = sample(1:p, replace = FALSE) # randomly assign nodes into groups
#' sbm = simu.SBM(conn1_mat, can_vec, train_obs_num, symm = TRUE, self = TRUE)
#' train_mat = sbm$obs_mat
#' train_list = lapply(1:ncol(train_mat), function(t) lowertri2mat(train_mat[,t], p, diag = TRUE))
#' pi_mat = matrix(0.9, p, p)
#' train_eta_list = lapply(1:length(train_list), function(t) gen.missing(pi_mat, symm = TRUE))
#' train_miss_list = lapply(1:length(train_list), function(t) train_eta_list[[t]] * train_list[[t]])
#' pi_lb_hat = quantile(Reduce("+", train_eta_list)/train_obs_num, 0.05) # estimator of pi_lb
#' pi_ub_hat = quantile(Reduce("+", train_eta_list)/train_obs_num, 0.95) # estimator of pi_ub
#' C_lambda = 2/3
#' lambda = lambda.network.missing(1, length(train_miss_list), length(train_miss_list), 0.05, 
#'                                 rho = 0.509, pi_ub = pi_ub_hat, p, C_lambda)
#' graphon_miss_impute = softImpute.network.missing(train_miss_list, train_eta_list, lambda, 1)
#' graphon_miss_hat = graphon_miss_impute$u %*% diag(as.numeric(graphon_miss_impute$d)) %*% 
#'                    t(graphon_miss_impute$v)
#' rho_hat = quantile(graphon_miss_hat, 0.95)
#' rank_hat = sum(graphon_miss_impute$d != 0)
#' alpha_grid = c(0.05, 0.01)
#' permu_num = 20
#' threshold_len = 100
#' temp = calibrate.online.network.missing(train_miss_list, train_eta_list, threshold_len, alpha_grid, 
#'                    permu_num, pi_lb_hat, pi_ub_hat, rho_hat, rank_hat, C_lambda, delta = 5)
#' @seealso \code{\link{online.network.missing}} for detecting online change point.
calibrate.online.network.missing = function(train_miss_list, train_eta_list, threshold_len, alpha_grid, permu_num, pi_lb_hat, pi_ub_hat, rho_hat, rank_hat, C_lambda = 2/3, delta = 5){
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
  thresholds_array = array(1, c(threshold_len, floor(log2(threshold_len))-burnin_idx+2, length(alpha_grid)))
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
        thresholds_array[t,1:length(N_grid),ind_alpha] = C_vec[ind_alpha] * sapply(N_grid, function(s) threshold.network.missing(s, t, rank_hat, pi_lb_hat, p, rho_hat, pi_ub_hat, alpha))
      }
    }
  }
  return(list(C_lambda = C_lambda, rho_hat = rho_hat, rank_hat = rank_hat, pi_lb_hat = pi_lb_hat, pi_ub_hat = pi_ub_hat, thresholds_array = thresholds_array))
}


#' @title Online change point detection for network data with missing values.
#' @description  Perform online change point detection for network with missing values by controlling the false alarm rate at level alpha.
#' @param data_incomplete_list  A \code{list} of adjacency matrices (with entries being 0 or 1) with missing values being coercing to 0.
#' @param eta_list              A \code{list} of matrices associated with data_incomplete_list, each matrix indicates the missing entries in corresponding adjacency matrix.
#' @param alpha_grid            A \code{numeric} vector in (0,1) representing the desired false alarm rate.
#' @param thresholds_array      A \code{numeric} array of calibrated thresholds.
#' @param pi_ub_hat             A \code{numeric} scalar of the upper bound of the missing probability.
#' @param rho_hat               A \code{numeric} scalar of the sparsity parameter.
#' @param C_lambda              A \code{numeric} scalar of an absolute constant, which is set to be 2/3 by default.
#' @param delta                 An \code{integer} scalar of minimum spacing.
#' @return  Online change point estimator.
#' @export
#' @author  Haotian Xu
#' @references Dubey, Xu and Yu (2021) <arxiv:2110.06450>
#' @seealso \code{\link{calibrate.online.network.missing}} for calibrating thresholds.
online.network.missing = function(data_incomplete_list, eta_list, alpha_grid, thresholds_array, rho_hat, pi_ub_hat, C_lambda = 2/3, delta = 5){
  burnin_idx = ceiling(log2(2*delta))
  cpt_hat = rep(NA, length(alpha_grid))
  for(ind_alpha in 1:length(alpha_grid)){
    alpha = alpha_grid[ind_alpha]
    t = 4*delta
    FLAG = 0
    while(t < length(data_incomplete_list) & FLAG == 0){
      t = t + 1
      print(t)
      m = floor(log2(t)) - 1
      if(t - 2^(m+1) >= 2*delta){
        m = m + 1
      }
      if(m < burnin_idx){
        N_grid = c(2*delta+1)
      }else{
        N_grid = c(2*delta+1, 2^(burnin_idx:m))
      }
      j = 1
      while(j <= length(N_grid) & FLAG == 0){
        FLAG = rcpp_CUSUM(data_incomplete_list, eta_list, 1, t, N_grid[j], alpha, rho_hat, pi_ub_hat, C_lambda, delta) > thresholds_array[t,j,ind_alpha]
        j = j + 1
      }
    }
    cpt_hat[ind_alpha] = t
  }
  return(cpt_hat)
}