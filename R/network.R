#' @title Simulate a Stochastic Block Model (of symmetric type and without change point).
#' @description  Simulate a Stochastic Block Model (without change point). The data generated are lower diagonal as the matrices are symmetry and 0 on the diagonal.
#' @param connec_mat  A \code{numeric} symmetric matrix representing the connectivity matrix (entries in [0,1]).
#' @param can_vec     A \code{integer} p-dim vector representing the candidate vector.
#' @param n           A \code{integer} scalar representing the number of observations.
#' @param ...        Additional arguments.
#' @return  A (p*(p-1)/2)-by-n matrix, with each column be the vectorized adjacency matrix.
#' @export
#' @author 
#' @examples
#' TO DO
#' 
simu.SBM = function(connec_mat, can_vec, n, ...){
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
  SBM_mean_vec = SBM_mean_mat[lower.tri(SBM_mean_mat, diag = F)]
  obs_mat = t(sapply(SBM_mean_vec, function(x) rbinom(n, 1, x)))
  return(obs_mat)
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



#' @title Binary segmentation for network change points detection.
#' @description  Perform binary segmentation for network change points detection.
#' @param data_mat1  A \code{numeric} matrix of observations with with horizontal axis being time, and with each column be the vectorized adjacency matrix.
#' @param data_mat2  A \code{numeric} matrix of observations with with horizontal axis being time, and with each column be the vectorized adjacency matrix (data_mat1 and data_mat2 are independent and have the same dimensions ).
#' @param s          A \code{integer} scalar of starting index.
#' @param e          A \code{integer} scalar of ending index.
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
#' @author
#' @examples
#' y = c(rnorm(100, 0, 1), rnorm(100, 10, 10), rnorm(100, 40, 10))
BS.network = function(data_mat1, data_mat2, s, e, delta, level = 0, ...){
  S = NULL
  Dval = NULL
  Level = NULL
  Parent = NULL
  if(e-s <= delta){
    return(list(S = S, Dval = Dval, Level = Level, Parent = Parent))
  }else{
    level = level + 1
    parent = matrix(c(s, e), nrow = 2)
    a = sapply(seq((s+1), (e-1), 1), function(x) changepoints:::CUSUM.innerprod(data_mat1,data_mat2, s, e, x))
    best_value = max(a)
    best_t = which.max(a) + s
    temp1 = BS.network(data_mat1, data_mat2, s, best_t-1, delta, level)
    temp2 = BS.network(data_mat1, data_mat2, best_t, e, delta, level)
    S = c(temp1$S, best_t, temp2$S)
    Dval = c(temp1$Dval, best_value, temp2$Dval)
    Level = c(temp1$Level, level, temp2$Level)
    Parent = cbind(temp1$Parent, parent, temp2$Parent)
    result = list(S = S, Dval = Dval, Level = Level, Parent = Parent)
    class(result) = "BS"
    return(result)
  }
}



#' @title Local refinement for network change points detection.
#' @description Perform local refinement for network change points detection.
#' @param cpt_init  A \code{integer} vector of initial change points estimation (sorted in strictly increasing order).
#' @param DATA      A \code{numeric} matrix of observations.
#' @param tau1      A \code{numeric} scalar corresponding to the first parameter of the USVT.
#' @param tau2      A \code{numeric} scalar corresponding to the second parameter of the USVT.
#' @param w         A \code{numeric} scalar of weight for interpolation of starting and ending indices.
#' @param ...       Additional arguments.
#' @return  A \code{numeric} vector of locally refined change point locations.
#' @export
#' @author 
#' @examples
#' TO DO
local.refine.network = function(cpt_init, DATA, tau1, tau2 = Inf, w = 1/3, ...){
  n = ncol(DATA)
  p = nrow(DATA)
  cpt_init_ext = c(0, cpt_init, n)
  K = length(cpt_init)
  cpt_refined = rep(0, K+1)
  for (k in 1:K){
    s_inter = w*cpt_init_ext[k] + (1-w)*cpt_init_ext[k+1]
    e_inter = (1-w)*cpt_init_ext[k+1] + w*cpt_init_ext[k+2]
    lower = ceiling(s_inter) + 1
    upper = floor(e_inter) - 1
    b_mat = sapply(lower:upper, function(eta) CUSUM.vec(DATA, s_inter, e_inter, eta))
    cpt_refined[k+1] = ceiling(s_inter) + which.max(sapply(lower:upper, function(x) USVT.norm(b_mat[,x-lower+1], tau1, tau2)))
  }
  return(cpt_refined[-1])
}


#' @title Internal Function: Compute the Universal Singular Value Thresholding of a symmetric matrix.
#' @param symm_mat   A \code{numeric} matrix of observations with with horizontal axis being time, and vertical axis being dimension.
#' @param tau1       A positive \code{numeric} scalar corresponding to the threshold for singular values of input matrix.
#' @param tau2       A positive \code{numeric} scalar corresponding to the threshold for entries of output matrix.
#' @return  A \code{numeric} matrix.
#' @noRd
USVT = function(symm_mat, tau1, tau2){
  p = dim(symm_mat)[1]
  result_mat = matrix(0, nrow = p, ncol = p)
  re = eigen(symm_mat, symmetric = TRUE)
  kk1 = length(which(re$values > tau1))
  kk2 = length(which(re$values < (-1)*tau1)) 
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
  result_mat[result_mat > tau2] = tau2
  result_mat[result_mat < (-1)*tau2] = (-1)*tau2
  return(result_mat)
}


#' @title Internal Function: Compute the Frobenius norm of USVT matrix (only take the lower triangular part of the matrix to reduce computation complexity).
#' @param cusum_vec  A \code{numeric} p*(p-1)/2-dim cusum vector.
#' @param tau1       A positive \code{numeric} scalar corresponding to the threshold for singular values of input matrix.
#' @param tau2       A positive \code{numeric} scalar corresponding to the threshold for entries of output matrix.
#' @return  A \code{numeric} scalar.
#' @noRd
USVT.norm = function(cusum_vec, tau1, tau2 = Inf){
  p = 1/2 + sqrt(2*length(cusum_vec) + 1/4) #obtain p
  cusum_mat = matrix(0, nrow = p, ncol=p)
  cusum_mat[gen.lower.coordinate(p)] = cusum_vec
  cusum_mat = cusum_mat + t(cusum_mat)
  #tau1=sqrt(p*rho)/2
  cusum_mat_temp = USVT(cusum_mat, tau1, tau2)
  return(norm(cusum_mat_temp, type="F"))
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
      FLAG = (sum(CUSUM.vec(data_mat2, 0, s_j, t) * B_tilde_vec)/B_tilde_norm > b_vec[t-1]) * (B_tilde_norm > c*sqrt(log(t/alpha))) 
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
      FLAG = (sum(CUSUM.vec(data_mat2, 0, s_j, t) * B_tilde_vec)/B_tilde_norm > b_vec[t-1]) * (B_tilde_norm > c*sqrt(log(gamma))) 
      j = j + 1
    }
  }
  return(t)
}



