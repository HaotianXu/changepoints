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
#' p = 100 # number of nodes
#' rho = 0.5 # sparsity parameter
#' block_num = 3 # number of groups for SBM
#' n = 150 # sample size for each segment
#' conn1_mat = rho * matrix(c(0.6,1,0.6,1,0.6,0.5,0.6,0.5,0.6), nrow = 3) # connectivity matrix for the first and the third segments
#' conn2_mat = rho * matrix(c(0.6,0.5,0.6,0.5,0.6,1,0.6,1,0.6), nrow = 3) # connectivity matrix for the second segment
#' set.seed(1)
#' can_vec = sample(1:p, replace = F) # randomly assign nodes into groups
#' sbm1 = simu.SBM(conn1_mat, can_vec, n, symm = TRUE, self = TRUE)
#' sbm2 = simu.SBM(conn2_mat, can_vec, n, symm = TRUE, self = TRUE)
#' data_mat = cbind(sbm1$obs_mat, sbm2$obs_mat)
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
  n_st = t - s + 1
  n_se = e - s + 1
  n_te = e - t
  p = dim(data_mat)[1]
  if(t-s<3 | e-t<2){
    result_vec = rep(0, p)
  }else{
    result_vec = sqrt(n_te/(n_se*n_st)) * rowSums(data_mat[,s:t]) - sqrt(n_st/(n_se*n_te)) * rowSums(data_mat[,(t+1):e])
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
#' p = 100 # number of nodes
#' rho = 0.5 # sparsity parameter
#' block_num = 3 # number of groups for SBM
#' n = 150 # sample size for each segment
#' conn1_mat = rho * matrix(c(0.6,1,0.6,1,0.6,0.5,0.6,0.5,0.6), nrow = 3) # connectivity matrix for the first and the third segments
#' conn2_mat = rho * matrix(c(0.6,0.5,0.6,0.5,0.6,1,0.6,1,0.6), nrow = 3) # connectivity matrix for the second segment
#' set.seed(1)
#' can_vec = sample(1:p, replace = F) # randomly assign nodes into groups
#' sbm1 = simu.SBM(conn1_mat, can_vec, n, symm = TRUE, self = TRUE)
#' sbm2 = simu.SBM(conn2_mat, can_vec, n, symm = TRUE, self = TRUE)
#' data_mat = cbind(sbm1$obs_mat, sbm2$obs_mat)
#' data_mat1 = data_mat[,seq(1,ncol(data_mat),2)]
#' data_mat2 = data_mat[,seq(2,ncol(data_mat),2)]
#' M = 120
#' intervals = WBS.intervals(M = M, lower = 1, upper = ncol(data_mat1))
#' temp = WBS.network(data_mat1, data_mat2, 1, ncol(data_mat1), intervals$Alpha, intervals$Beta, delta = 5)
#' rho_hat = quantile(rowMeans(data_mat), 0.95)
#' tau = p*rho_hat*(log(n))^2/20 # default threshold given in the paper
#' cpt_init = unlist(threshold.BS(temp, tau)$change_points[1])
#' cpt_refined = local.refine.network(cpt_init, data_mat1, data_mat2, self = TRUE, tau2 = p*rho_hat/3, tau3 = Inf)
#' cpt_WBS = 2*cpt_init
#' cpt_refin = 2*cpt_refined
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
      temp = rep(0, Beta_new[m] - Alpha_new[m] - 2*delta + 1)
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
#' @param tau2       A positive \code{numeric} scalar for USVT corresponding to the threshold for singular values of input matrix.
#' @param tau3       A positive \code{numeric} scalar for USVT corresponding to the threshold for entries of output matrix.
#' @param ...       Additional arguments.
#' @return  A \code{numeric} vector of locally refined change point locations.
#' @export
#' @author  Daren Wang & Haotian Xu
#' @references Wang D, Yu Y, Rinaldo A. Optimal change point detection and localization in sparse dynamic networks. The Annals of Statistics. 2021 Feb;49(1):203-32.
#' @examples
#' p = 100 # number of nodes
#' rho = 0.5 # sparsity parameter
#' block_num = 3 # number of groups for SBM
#' n = 150 # sample size for each segment
#' conn1_mat = rho * matrix(c(0.6,1,0.6,1,0.6,0.5,0.6,0.5,0.6), nrow = 3) # connectivity matrix for the first and the third segments
#' conn2_mat = rho * matrix(c(0.6,0.5,0.6,0.5,0.6,1,0.6,1,0.6), nrow = 3) # connectivity matrix for the second segment
#' set.seed(1)
#' can_vec = sample(1:p, replace = F) # randomly assign nodes into groups
#' sbm1 = simu.SBM(conn1_mat, can_vec, n, symm = TRUE, self = TRUE)
#' sbm2 = simu.SBM(conn2_mat, can_vec, n, symm = TRUE, self = TRUE)
#' data_mat = cbind(sbm1$obs_mat, sbm2$obs_mat)
#' data_mat1 = data_mat[,seq(1,ncol(data_mat),2)]
#' data_mat2 = data_mat[,seq(2,ncol(data_mat),2)]
#' M = 120
#' intervals = WBS.intervals(M = M, lower = 1, upper = ncol(data_mat1))
#' temp = WBS.network(data_mat1, data_mat2, 1, ncol(data_mat1), intervals$Alpha, intervals$Beta, delta = 5)
#' rho_hat = quantile(rowMeans(data_mat), 0.95)
#' tau = p*rho_hat*(log(n))^2/20 # default threshold given in the paper
#' cpt_init = unlist(threshold.BS(temp, tau)$change_points[1])
#' cpt_refined = local.refine.network(cpt_init, data_mat1, data_mat2, self = TRUE, tau2 = p*rho_hat/3, tau3 = Inf)
#' cpt_WBS = 2*cpt_init
#' cpt_refin = 2*cpt_refined
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
  if(p != round(p)){
    stop("Either cusum_vec or self is not correctly specified.")
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
  }else{
    p = 1/2 + sqrt(2*length(cusum_vec1) + 1/4) #obtain p
    cusum_mat1 = lowertri2mat(cusum_vec1, p, diag = self)
  }
  Theta_mat2 = USVT(cusum_vec2, self = FALSE, tau2, tau3)
  return(sum(cusum_mat1*Theta_mat2))
}




#' @title Internal Function: Compute value of CUSUM statistic (multivariate) at the current time t and s being the splitting time point. The default choice of the tuning parameters tau1, tau2 and tau3 are used (see Section 4.1 of the reference).
#' @param data_mat1  A \code{numeric} matrix of observations with with horizontal axis being time, and with each column be the vectorized adjacency matrix.
#' @param data_mat2  A \code{numeric} matrix of observations with with horizontal axis being time, and with each column be the vectorized adjacency matrix (data_mat1 and data_mat2 are independent and have the same dimensions ).
#' @param self       A \code{logic} scalar indicating if adjacency matrices are required to have self-loop.
#' @param t          A \code{integer} scalar of current time index.
#' @param s          A \code{integer} scalar of splitting index.
#' @param rho        A \code{numeric} scalar of sparsity parameter of the network data.
#' @param alpha      A \code{numeric} scalar in (0,1) representing the level.
#' @param gamma     An \code{integer} scalar of desired average run length.
#' @return  A \code{numeric} scalar of value of CUSUM statistic.
#' @noRd
data.split.statistic = function(data_mat1, data_mat2, self = FALSE, t, s, rho, alpha = NULL, gamma = NULL){
  C = 32.1*(2^.25)*exp(2)
  if(self == TRUE){
    p = sqrt(2*nrow(data_mat1) + 1/4) - 1/2
  }else{
    p = 1/2 + sqrt(2*nrow(data_mat1) + 1/4) #obtain p
  }
  if(is.null(alpha)+is.null(gamma)!=1){
    stop("Either alpha or gamma should be provided.")
  }
  if(!is.null(alpha)){
    tau1 = ((C/100)*sqrt(p*rho) + sqrt(2*log(  t*(t+1)*log(t)/(alpha*log(2))   )    )      )/5
    tau2 =  rho*sqrt(  (t-s)*s/t)/10
    tau3 = (C/100)*sqrt( log(t/alpha)   )/50
  }else{
    tau1 =((C/100)*sqrt(p*rho)  + sqrt(  2*log((2*(gamma+1)*(gamma+1)*log(gamma+1) )/log(2))  ))/5
    tau2 =  rho*sqrt(  (t-s)*s/t)/10
    tau3 = (C/100)*sqrt( log(gamma)   )/50
  }

  vec1 = CUSUM.vec(data_mat1, 1, t, s)
  vec2 = CUSUM.vec(data_mat2, 1, t, s)  
  
  mat1 = lowertri2mat(vec1, p, diag = self)
  mat3 = USVT(vec2, self = FALSE, tau1, tau2)
  frob_norm = sum(as.vector(mat3)^2)
  if(frob_norm < tau3){
    return(0)
  }
  aux = sum(as.vector(mat1) * as.vector(mat3))/frob_norm
  return(aux)
}



#' @title Online changepoint detection for network data by controlling the false alarm rate at level alpha.
#' @description  Perform online changepoint detection for network data by controlling the false alarm rate at level alpha or controlling the average run length gamma. The default choice of the tuning parameters tau1, tau2 and tau3 are used (see Section 4.1 of the reference). 
#' @param data_mat1  A \code{numeric} matrix of observations with with horizontal axis being time, and with each column be the vectorized adjacency matrix.
#' @param data_mat2  A \code{numeric} matrix of observations with with horizontal axis being time, and with each column be the vectorized adjacency matrix (data_mat1 and data_mat2 are independent and have the same dimensions ).
#' @param b_vec      A \code{numeric} vector of thresholds b_t with t >= 2.
#' @param train_mat  A \code{numeric} matrix of training data from a pre-change distribution(no changepoint), which is only needed to when b_vec is NULL in order to calibrate b_t.
#' @param alpha      A \code{numeric} scalar in (0,1) representing the level.
#' @param gamma      An \code{integer} scalar of desired average run length.
#' @param permu_num  An \code{integer} scalar of number of random permutation for calibration.
#' @param ...        Additional arguments.
#' @return  A \code{list} with the structure:
#' \itemize{
#'  \item t:           Estimated changepoint.
#'  \item score:       A \code{numeric} vector of computed cumsum statistics.
#'  \item b_vec        A \code{numeric} vector of thresholds b_t with t >= 2.
#' } 
#' @export
#' @author  Oscar Hernan Madrid Padilla & Haotian Xu
#' @references Yu Y, Padilla, O, Wang D, Rinaldo A. Optimal network online change point localisation. arXiv preprint arXiv:2101.05477.
#' @examples
#' p = 100 # number of nodes
#' rho = 0.5 # sparsity parameter
#' block_num = 3 # number of groups for SBM
#' n = 150 # sample size for each segment
#' conn1_mat = rho * matrix(c(0.6,1,0.6,1,0.6,0.5,0.6,0.5,0.6), nrow = 3) # connectivity matrix for the first and the third segments
#' conn2_mat = rho * matrix(c(0.6,0.5,0.6,0.5,0.6,1,0.6,1,0.6), nrow = 3) # connectivity matrix for the second segment
#' set.seed(1)
#' can_vec = sample(1:p, replace = F) # randomly assign nodes into groups
#' sbm1 = simu.SBM(conn1_mat, can_vec, n, symm = TRUE, self = TRUE)
#' sbm2 = simu.SBM(conn2_mat, can_vec, n, symm = TRUE, self = TRUE)
#' data_mat = cbind(sbm1$obs_mat, sbm2$obs_mat)
#' data_mat1 = data_mat[,seq(1,ncol(data_mat),2)]
#' data_mat2 = data_mat[,seq(2,ncol(data_mat),2)]
#' train_mat = simu.SBM(conn1_mat, can_vec, n = 200, symm = TRUE, self = TRUE)$obs_mat
#' temp = online.network(data_mat1, data_mat2, self = TRUE, b_vec = NULL, train_mat, alpha = 0.05, gamma = NULL, permu_num = 100)
#' cpt_hat = 2 * temp$t
#' temp2 = online.network(data_mat1, data_mat2, self = TRUE, b_vec = NULL, train_mat, alpha = NULL, gamma = 200, permu_num = 100)
#' cpt_hat2 = 2 * temp2$t
online.network = function(data_mat1, data_mat2, self = TRUE, b_vec = NULL, train_mat = NULL, alpha = NULL, gamma = NULL, permu_num = NULL, ...){
  n = ncol(data_mat1)
  if(self == TRUE){
    p = sqrt(2*nrow(data_mat1) + 1/4) - 1/2
  }else{
    p = 1/2 + sqrt(2*nrow(data_mat1) + 1/4) #obtain p
  }
  if(is.null(alpha)+is.null(gamma)!=1){
    stop("Either alpha or gamma should be provided.")
  }
  if(!is.null(b_vec)){
    rho_hat = quantile(rowMeans(cbind(data_mat1, data_mat2)), 0.95)
    if(n - length(b_vec) != 1){
      stop("b_vec should be the vector of thresholds b_t with t >= 2.")
    }
  }else{
    if(is.null(train_mat)){
      stop("Given b_vec is missing, train_mat should be provided to calibrate b_vec.")
    }
    rho_hat = quantile(rowMeans(cbind(data_mat1, data_mat2, train_mat)), 0.95)
    train_mat1 = train_mat[,seq(1,ncol(train_mat),2)]
    train_mat2 = train_mat[,seq(2,ncol(train_mat),2)]
    obs_train = min(ncol(train_mat1), ncol(train_mat2))
    if(obs_train > n){
      # only use part of train_mat1 if it's sample size is bigger than that of data_mat1
      train_mat1 = train_mat1[,(obs_train-n+1):obs_train]
      train_mat2 = train_mat2[,(obs_train-n+1):obs_train]
      obs_train = n
    }
    scores = matrix(0, permu_num, obs_train-1)
    if(!is.null(alpha)){
      C_vec = rep(NA, permu_num)
      trend = sapply(2:obs_train, function(t) sqrt(rho_hat*log(t/alpha)))
      print("Start calibrating the thresholds b_t:")
      pb = txtProgressBar(min = 0, max = permu_num, style = 3)
      counter = 0
      for(sim in 1:permu_num){
        idx_permu_odd = sample(1:obs_train)
        idx_permu_even = sample(1:obs_train)
        train_permu_odd = train_mat1[,idx_permu_odd]
        train_permu_even = train_mat2[,idx_permu_even]
        
        for(t in 2:obs_train){
          if(t>10){
            m = floor(log(t)/log(2))-1
            N_grid = as.matrix(2^{1:m})
            aux = apply(N_grid, 1, function(par){data.split.statistic(train_permu_odd, train_permu_even, self = FALSE, t, par, rho = rho_hat, alpha = alpha)})
            scores[sim, t-1] = max(aux)
          }
        }
        C_vec[sim] = max(scores[sim,]/trend)
        counter = counter + 1
        setTxtProgressBar(pb, counter)
      }
      b_vec = quantile(C_vec, 1-alpha) * sapply(2:n, function(t) sqrt(rho_hat*log(t/alpha)))
    }else if(!is.null(gamma)){
      if(gamma > obs_train){
        gamma = obs_train
        warning(paste0("gamma is set to be ", obs_train, ". To allow larger value of gamma, please increase the sample size of train_mat."))
      }
      C_mat = matrix(NA, permu_num, obs_train-1)
      trend = sqrt(rho_hat*log(gamma))
      print("Start calibrating the thresholds b_t:")
      pb = txtProgressBar(min = 0, max = permu_num, style = 3)
      counter = 0
      for(sim in 1:permu_num){
        idx_permu_odd = sample(1:obs_train)
        idx_permu_even = sample(1:obs_train)
        train_permu_odd = train_mat1[,idx_permu_odd]
        train_permu_even = train_mat2[,idx_permu_even]
        
        for(t in 2:obs_train){
          if(t>10){
            m = floor(log(t)/log(2))-1
            N_grid = as.matrix(2^{1:m})
            aux = apply(N_grid, 1, function(par){data.split.statistic(train_permu_odd, train_permu_even, self = FALSE, t, par, rho = rho_hat, gamma = gamma)})
            scores[sim, t-1] = max(aux)
          }
        }
        C_mat[sim,] = scores[sim,]/trend
        counter = counter + 1
        setTxtProgressBar(pb, counter)
      }
      C_grid = seq(min(C_mat), max(C_mat), length = 300)
      alarm = rep(0, length(C_grid))
      for(j in 1:length(alarm)){
        aux = apply(C_mat, 1,function(x){min(c(which(x > C_grid[j]), obs_train))})
        alarm[j] = mean(aux)
      }
      ind = which.min(abs(alarm-gamma))
      b_vec = rep(C_grid[ind] * trend, n-1)
    }
    print("Finish calibration.")
  }
  score = rep(0, n)
  for(t in 2:n){
    if(t > 3){
      m = floor(log(t)/log(2))-1
      N_grid = as.matrix(2^{1:m})
      aux = apply(N_grid, 1, function(par){data.split.statistic(data_mat1, data_mat2, self = FALSE, t, par, rho, alpha, gamma)})
      score[t] = max(aux)
    }
    if(score[t] > b_vec[t]){
      break
    }
  }
  #ind =  which(score/b >1)
  if(score[t] > b_vec[t]){
    return(list(t = t, score = score, b_vec = b_vec))
  }
  return(list(t = Inf, score = score, b_vec = b_vec))
}
