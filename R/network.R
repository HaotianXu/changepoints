#' @title Simulate a Stochastic Block Model (without change point).
#' @description  Simulate a Stochastic Block Model (without change point). The generated data is a matrix with each column corresponding to the vectorized adjacency (sub)matrix at a time point. For example, if the network matrix is required to be symmetric and without self-loop, only the strictly lower diagonal entries are considered.
#' @param connec_mat  A \code{numeric} symmetric matrix representing the connectivity matrix (each entry takes value in [0,1]).
#' @param can_vec     A \code{integer} p-dim vector of node indices. can_vec is then divided into subvectors corresponding to blocks.
#' @param n           A \code{integer} scalar representing the number of observations.
#' @param symm        A \code{logic} scalar indicating if adjacency matrices are required to be symmetric.
#' @param self        A \code{logic} scalar indicating if adjacency matrices are required to have self-loop.
#' @param ...        Additional arguments.
#' @return  A \code{list} with the structure:
#' \itemize{
#'  \item obs_mat:       A matrix, with each column be the vectorized adjacency (sub)matrix. For example, if "symm = TRUE" and "self = FALSE", only the strictly lower triangular matrix is considered.
#'  \item graphon_mat:   Underlying graphon matrix.
#' } 
#' @export
#' @author 
#' @examples
#' d = 100 # number of nodes
#' rho = 0.5 # sparsity parameter
#' block_num = 3 # number of groups for SBM
#' n = 150 # sample size for each segment
#' conn1_mat = rho * matrix(c(0.6,1,0.6,1,0.6,0.5,0.6,0.5,0.6), nrow = 3) # connectivity matrix for the first and the third segments
#' conn2_mat = rho * matrix(c(0.6,0.5,0.6,0.5,0.6,1,0.6,1,0.6), nrow = 3) # connectivity matrix for the second segment
#' set.seed(1)
#' can_vec = sample(1:d, replace = F) # randomly assign nodes into groups
#' sbm = simu.SBM(conn1_mat, can_vec, train_obs_num, symm = TRUE, self = TRUE)
simu.SBM = function(connec_mat, can_vec, n, symm = FALSE, self = TRUE, ...){
  block_num = dim(connec_mat)[1]
  p = length(can_vec)
  SBM_mean_mat = matrix(0, nrow = p, ncol = p)
  aa = p / block_num
  for(i in 1:block_num){
    for(j in 1:block_num){
      can_temp1 = rep(0, p)
      can_temp1[can_vec[(1+(i-1)*aa):(i*aa)]] = 1
      can_temp2 = rep(0, p)
      can_temp2[can_vec[(1+(j-1)*aa):(j*aa)]] = 1
      SBM_mean_mat = SBM_mean_mat + connec_mat[i,j] * can_temp1 %*% t(can_temp2)
    }
  }
  if((symm == FALSE) & (self == TRUE)){
    SBM_mean_vec = as.vector(SBM_mean_mat)
  }else if((symm == FALSE) & (self == FALSE)){
    diag(SBM_mean_mat) = 0
    SBM_mean_vec = as.vector(SBM_mean_mat)
  }else if((symm == TRUE) & (self == TRUE)){
    SBM_mean_vec = SBM_mean_mat[lower.tri(SBM_mean_mat, diag = TRUE)]
  }else{
    diag(SBM_mean_mat) = 0
    SBM_mean_vec = SBM_mean_mat[lower.tri(SBM_mean_mat, diag = F)]
  }
  obs_mat = t(sapply(SBM_mean_vec, function(x) rbinom(n, 1, x)))
  return(list(obs_mat = obs_mat, graphon_mat = SBM_mean_mat))
}


#' @title Internal Function: Compute value of CUSUM statistic (multivariate).
#' @param data_mat  A \code{numeric} matrix of observations with with horizontal axis being time.
#' @param s         A \code{integer} scalar of starting index.
#' @param e         A \code{integer} scalar of ending index.
#' @param t         A \code{integer} scalar of splitting index.
#' @return  A \code{numeric} vector of value of CUSUM statistic.
#' @noRd
CUSUM.vec = function(data_mat, s, e, t){
  n_st = t - s
  n_se = e - s
  n_te = e - t
  p = dim(data_mat)[1]
  if(t-s<3 | e-t<2){
    result_vec = rep(0, p)
  }else{
    result_vec = sqrt(n_te/(n_se*n_st)) * rowSums(data_mat[,(s+1):t]) - sqrt(n_st/(n_se*n_te)) * rowSums(data_mat[,(t+1):e])
  }
  return(result_vec)
}


#' @title Internal Function: Compute inner product of two CUSUM vectors based on two independent samples.
#' @param data_mat1  A \code{numeric} matrix of observations with with horizontal axis being time, and vertical axis being dimension.
#' @param data_mat2  A \code{numeric} matrix of observations with with horizontal axis being time, and vertical axis being dimension (data_mat1 and data_mat2 are independent and have the same dimensions ).
#' @param s          A \code{integer} scalar of starting index.
#' @param e          A \code{integer} scalar of ending index.
#' @param t          A \code{integer} scalar of splitting index.
#' @return  A \code{numeric} scalar representing the inner product of two CUSUM vectors based on two independent samples.
#' @noRd
CUSUM.innerprod = function(data_mat1, data_mat2, s, e, t){
  return(sum(CUSUM.vec(data_mat1, s, e, t) * CUSUM.vec(data_mat2, s, e, t)))
}



#' @title Wild binary segmentation for network change points detection.
#' @description  Perform wild binary segmentation for network change points detection.
#' @param data_mat1  A \code{numeric} matrix of observations with with horizontal axis being time, and with each column be the vectorized adjacency matrix.
#' @param data_mat2  A \code{numeric} matrix of observations with with horizontal axis being time, and with each column be the vectorized adjacency matrix (data_mat1 and data_mat2 are independent and have the same dimensions ).
#' @param s          A \code{integer} scalar of starting index.
#' @param e          A \code{integer} scalar of ending index.
#' @param Alpha      A \code{integer} vector of starting indices of random intervals.
#' @param Beta       A \code{integer} vector of ending indices of random intervals.
#' @param delta      A positive \code{integer} scalar of minimum spacing.
#' @param level      Should be fixed as 0.
#' @param ...        Additional arguments.
#' @return  A \code{list} with the structure:
#' \itemize{
#'  \item S:           A vector of estimated changepoints (sorted in strictly increasing order).
#'  \item Dval:        A vector of values of CUSUM statistic based on KS distance.
#'  \item Level:       A vector representing the levels at which each change point is detected.
#'  \item Parent:      A matrix with the starting indices on the first row and the ending indices on the second row.
#' } 
#' @export
#' @author  Daren Wang & Haotian Xu
#' @references Wang D, Yu Y, Rinaldo A. Optimal change point detection and localization in sparse dynamic networks. The Annals of Statistics. 2021 Feb;49(1):203-32.
#' @examples
#' y = c(rnorm(100, 0, 1), rnorm(100, 10, 10), rnorm(100, 40, 10))
WBS.network = function(data_mat1, data_mat2, s, e, Alpha, Beta, delta, level = 0, ...){
  Alpha_new = pmax(Alpha, s)
  Beta_new = pmin(Beta, e)
  xi = 1/64 #using the shrunk constant given in the paper
  Alpha_new2 = Alpha_new
  Beta_new2 = Beta_new
  Alpha_new = ceiling((1-xi)*Alpha_new2 + xi*Beta_new2)
  Beta_new = ceiling((1-xi)*Beta_new2 + xi*Alpha_new2)
  idx = which(Beta_new - Alpha_new > 2*delta)
  Alpha_new = Alpha_new[idx]
  Beta_new = Beta_new[idx]
  M = length(Alpha_new)
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
      for(t in (Alpha_new[m]+delta):(Beta_new[m]-delta)){
        temp[t-(Alpha_new[m]+delta)+1] = changepoints:::CUSUM.innerprod(data_mat1,data_mat2, Alpha_new[m], Beta_new[m], t)
      }
      best_value = max(temp)
      best_t = which.max(temp) + Alpha_new[m] + delta - 1
      a[m] = best_value
      b[m] = best_t
    }
    m_star = which.max(a)
  }
  temp1 = WBS.network(data_mat1, data_mat2, s, b[m_star]-1, Alpha, Beta, delta, level)
  temp2 = WBS.network(data_mat1, data_mat2, b[m_star], e, Alpha, Beta, delta, level)
  S = c(temp1$S, b[m_star], temp2$S)
  Dval = c(temp1$Dval, a[m_star], temp2$Dval)
  Level = c(temp1$Level, level, temp2$Level)
  Parent = cbind(temp1$Parent, parent, temp2$Parent)
  result = list(S = S, Dval = Dval, Level = Level, Parent = Parent)
  class(result) = "BS"
  return(result)
}



#' @title Local refinement for network change points detection.
#' @description Perform local refinement for network change points detection.
#' @param cpt_init   A \code{integer} vector of initial change points estimation (sorted in strictly increasing order).
#' @param data_mat1  A \code{numeric} matrix of observations with with horizontal axis being time, and with each column be the vectorized adjacency matrix.
#' @param data_mat2  A \code{numeric} matrix of observations with with horizontal axis being time, and with each column be the vectorized adjacency matrix (data_mat1 and data_mat2 are independent and have the same dimensions ).
#' @param self       A \code{logic} scalar indicating if adjacency matrices are required to have self-loop.
#' @param tau2       A \code{numeric} scalar corresponding to the first parameter of the USVT.
#' @param tau3       A \code{numeric} scalar corresponding to the second parameter of the USVT.
#' @param ...       Additional arguments.
#' @return  A \code{numeric} vector of locally refined change point locations.
#' @export
#' @author 
#' @examples
#' TO DO
local.refine.network = function(cpt_init, data_mat1, data_mat2, self = FALSE, tau2, tau3 = Inf, ...){
  obs_num = ncol(data_mat1)
  cpt_init_ext = c(0, cpt_init, obs_num)
  K = length(cpt_init)
  cpt_refined = rep(0, K+1)
  for(k in 1:K){
    s_inter = ceiling(0.5*cpt_init_ext[k] + 0.5*cpt_init_ext[k+1])
    e_inter = floor(0.5*cpt_init_ext[k+1] + 0.5*cpt_init_ext[k+2])
    Delta_tilde = sqrt((e_inter-cpt_init_ext[k+1])*(cpt_init_ext[k+1]-s_inter)/(e_inter-s_inter))
    A_tilde = sapply((s_inter + 1):(e_inter - 1), function(eta) CUSUM.vec(data_mat1, s_inter, e_inter, eta))
    cpt_refined[k+1] = s_inter + which.max(sapply((s_inter + 1):(e_inter - 1), function(x) USVT.norm(A_tilde[,x-s_inter], CUSUM.vec(data_mat2, s_inter, e_inter, cpt_init_ext[k+1]), self, tau2, Delta_tilde*tau3)))
  }
  return(cpt_refined[-1])
}


#' @title Internal Function: Compute the Universal Singular Value Thresholding (USVT) of a symmetric matrix.
#' @param cusum_vec  A \code{numeric} vector corresponding to a cusum vector .
#' @param self       A \code{logic} scalar indicating if adjacency matrices are required to have self-loop.
#' @param tau2       A positive \code{numeric} scalar corresponding to the threshold for singular values of input matrix.
#' @param tau3       A positive \code{numeric} scalar corresponding to the threshold for entries of output matrix.
#' @return  A \code{numeric} matrix.
#' @noRd
USVT = function(cusum_vec, self = FALSE, tau2, tau3){
  if(self == TRUE){
    p = sqrt(2*length(cusum_vec) + 1/4) - 1/2
    cusum_mat = lowertri2mat(cusum_vec, p, diag = self)
  }else{
    p = 1/2 + sqrt(2*length(cusum_vec) + 1/4) #obtain p
    cusum_mat = lowertri2mat(cusum_vec, p, diag = self)
  }
  result_mat = matrix(0, nrow = p, ncol = p)
  re = eigen(cusum_mat, symmetric = TRUE)
  kk1 = length(which(re$values > tau2))
  kk2 = length(which(re$values < (-1)*tau2)) 
  if(kk1 + kk2 == 0){
    return(result_mat)
  }else{
    if(kk1 > 0){
      for(i in 1:kk1){
        result_mat = result_mat + re$values[i] * re$vectors[,i] %*% t(re$vectors[,i])
      }
    }
    if(kk2 > 0){
      for(i in 1:kk2){
        result_mat = result_mat + re$values[p-i+1] * re$vectors[,p-i+1] %*% t(re$vectors[,p-i+1])
      }
    }
  }
  result_mat[result_mat > tau3] = tau3
  result_mat[result_mat < (-1)*tau3] = (-1)*tau3
  return(result_mat)
}


#' @title Internal Function: Compute the Frobenius norm of USVT matrix.
#' @param cusum_vec1  A \code{numeric} vector corresponding to a cusum vector computed based on data_mat1.
#' @param cusum_vec2  A \code{numeric} vector corresponding to a cusum vector computed based on data_mat2.
#' @param self        A \code{logic} scalar indicating if adjacency matrices are required to have self-loop.
#' @param tau2        A positive \code{numeric} scalar corresponding to the threshold for singular values of input matrix.
#' @param tau3        A positive \code{numeric} scalar corresponding to the threshold for entries of output matrix.
#' @return  A \code{numeric} scalar.
#' @noRd
USVT.norm = function(cusum_vec1, cusum_vec2, self = FALSE, tau2, tau3 = Inf){
  if(self == TRUE){
    p = sqrt(2*length(cusum_vec1) + 1/4) - 1/2
    cusum_mat1 = lowertri2mat(cusum_vec1, p, diag = self)
    cusum_mat2 = lowertri2mat(cusum_vec2, p, diag = self)
  }else{
    p = 1/2 + sqrt(2*length(cusum_vec1) + 1/4) #obtain p
    cusum_mat1 = lowertri2mat(cusum_vec1, p, diag = self)
    cusum_mat2 = lowertri2mat(cusum_vec2, p, diag = self)
  }
  Theta_mat2 = USVT(cusum_mat2, tau2, tau3)
  return(sum(cusum_mat1*Theta_mat2))
}




online.network = function(data_mat1, data_mat2, b_vec, tau1_mat, tau2_mat, c, alpha){
  p = sqrt(nrow(data_mat2))
  t = 1
  FLAG = 0
  while(FLAG == 0){
    t = t + 1
    J = floor(log2(t))
    j = 0
    while(j < J & FLAG == 0){
      s_j = t - 2^j
      B_tilde_vec = as.vector(USVT(matrix(CUSUM.vec(data_mat2, 0, s_j, t), p, p), tau1_mat[s_j, t], tau2_mat[s_j, t]))
      B_tilde_norm = sqrt(sum(B_tilde_vec^2))
      FLAG = (sum(CUSUM.vec(data_mat1, 0, s_j, t) * B_tilde_vec)/B_tilde_norm > b_vec[t-1]) * (B_tilde_norm > c*sqrt(log(t/alpha))) 
      j = j + 1
    }
  }
  return(t)
}


online.network.variant = function(data_mat1, data_mat2, b_vec, tau1_mat, tau2_mat, c, gamma){
  p = sqrt(nrow(data_mat2))
  t = 1
  FLAG = 0
  while(FLAG == 0){
    t = t + 1
    J = floor(log2(t))
    j = 0
    while(j < J & FLAG == 0){
      s_j = t - 2^j
      B_tilde_vec = as.vector(USVT(matrix(CUSUM.vec(data_mat2, 0, s_j, t), p, p), tau1_mat[s_j, t], tau2_mat[s_j, t]))
      B_tilde_norm = sqrt(sum(B_tilde_vec^2))
      FLAG = (sum(CUSUM.vec(data_mat1, 0, s_j, t) * B_tilde_vec)/B_tilde_norm > b_vec[t-1]) * (B_tilde_norm > c*sqrt(log(gamma))) 
      j = j + 1
    }
  }
  return(t)
}



