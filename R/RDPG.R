#' @title Simulate a dot product graph (without change point).
#' @description  Simulate a dot product graph (without change point). The generated data is a matrix with each column corresponding to the vectorized adjacency (sub)matrix at a time point. For example, if the network matrix is required to be symmetric and without self-loop, only the strictly lower diagonal entries are considered.
#' @param x_mat       A \code{numeric} matrix representing the latent positions with horizontal axis being latent dimensions and vertical axis being nodes (each entry takes value in \eqn{[0,1]}).
#' @param n           A \code{integer} scalar representing the number of observations.
#' @param symm        A \code{logic} scalar indicating if adjacency matrices are required to be symmetric.
#' @param self        A \code{logic} scalar indicating if adjacency matrices are required to have self-loop.
#' @return  A \code{list} with the following structure:
#'  \item{obs_mat}{A matrix, with each column be the vectorized adjacency (sub)matrix. For example, if "symm = TRUE" and "self = FALSE", only the strictly lower triangular matrix is considered.}
#'  \item{graphon_mat}{Underlying graphon matrix.}
#' @export
#' @author Haotian Xu
#' @examples
#' p = 20 # number of nodes
#' n = 50 # sample size for each segment
#' lat_dim_num = 5 # number of latent dimensions
#' set.seed(1)
#' x_mat = matrix(runif(p*lat_dim_num), nrow = p, ncol = lat_dim_num)
#' x_tilde_mat = matrix(runif(p*lat_dim_num), nrow = p, ncol = lat_dim_num)
#' y_mat = rbind(x_tilde_mat[1:floor(p/4),], x_mat[(floor(p/4)+1):p,])
#' rdpg1 = simu.RDPG(x_mat, n, symm = TRUE, self = FALSE)
#' rdpg2 = simu.RDPG(y_mat, n, symm = TRUE, self = FALSE)
#' data1_mat = rdpg1$obs_mat
#' data2_mat = rdpg2$obs_mat
#' data_mat = cbind(data1_mat, data2_mat)
simu.RDPG = function(x_mat, n, symm = TRUE, self = FALSE){
  if(self == TRUE){
    stop("this function does not allow self-loop.")
  }
  p = nrow(x_mat)
  x_l2_vec = apply(x_mat, 1, function(x) (sum(x^2))^(1/2))
  ratio_x_mat = (x_mat %*% t(x_mat))/(x_l2_vec %*% t(x_l2_vec))
  diag(ratio_x_mat) = 0
  if((symm == TRUE) & (self == FALSE)){
    obs_mat = matrix(NA, p*(p-1)/2, n)
    ratio_vec = ratio_x_mat[lower.tri(ratio_x_mat, diag = F)]
    for(t in 1:n){
      obs_mat[,t] = rbinom(rep(1,length(ratio_vec)), rep(1,length(ratio_vec)), ratio_vec)
    }
  }else if((symm == FALSE) & (self == FALSE)){
    obs_mat = matrix(NA, p*p, n)
    ratio_vec = as.vector(ratio_x_mat)
    for(t in 1:n){
      obs_mat[,t] = rbinom(rep(1,length(ratio_vec)), rep(1,length(ratio_vec)), ratio_vec)
    }
  }
  return(list(obs_mat = obs_mat, graphon_mat = ratio_x_mat))
}


#' @title Wild binary segmentation for dependent dynamic random dot product graph models.
#' @description Perform wild binary segmentation for dependent dynamic random dot product graph models.
#' @param data_mat  A \code{numeric} matrix of observations with horizontal axis being time, and vertical axis being vectorized adjacency matrix.
#' @param lowerdiag A \code{logic} scalar. TRUE, if each row of data_mat is the vectorization of the strictly lower diagonal elements in an adjacency matrix. FALSE, if each row of data_mat is the vectorization of all elements in an adjacency matrix.
#' @param d         A \code{numeric} scalar of the number of leading singular values of an adjacency matrix considered in the scaled PCA algorithm.
#' @param Alpha     A \code{integer} vector of starting indices of random intervals.
#' @param Beta      A \code{integer} vector of ending indices of random intervals.
#' @param delta     A positive \code{integer} scalar of minimum spacing.
#' @return  An object of \code{\link[base]{class}} "BS", which is a \code{list} with the following structure:
#'  \item{S}{A vector of estimated change point locations (sorted in strictly increasing order).}
#'  \item{Dval}{A vector of values of CUSUM statistic.}
#'  \item{Level}{A vector representing the levels at which each change point is detected.}
#'  \item{Parent}{A matrix with the starting indices on the first row and the ending indices on the second row.}
#' @export
#' @author Oscar Hernan Madrid Padilla, Haotian Xu
#' @references Padilla, Yu and Priebe (2019) <arxiv:1911.07494>.
#' @seealso \code{\link{thresholdBS}} for obtaining change points estimation, \code{\link{tuneBSnonparRDPG}} for a tuning version.
#' @examples
#' ### generate data 
#' p = 20 # number of nodes
#' n = 50 # sample size for each segment
#' lat_dim_num = 5 # number of latent dimensions
#' set.seed(1)
#' x_mat = matrix(runif(p*lat_dim_num), nrow = p, ncol = lat_dim_num)
#' x_tilde_mat = matrix(runif(p*lat_dim_num), nrow = p, ncol = lat_dim_num)
#' y_mat = rbind(x_tilde_mat[1:floor(p/4),], x_mat[(floor(p/4)+1):p,])
#' rdpg1 = simu.RDPG(x_mat, n, symm = TRUE, self = FALSE)
#' rdpg2 = simu.RDPG(y_mat, n, symm = TRUE, self = FALSE)
#' data1_mat = rdpg1$obs_mat
#' data2_mat = rdpg2$obs_mat
#' data_mat = cbind(data1_mat, data2_mat)
#' ### detect change points
#' M = 30 # number of random intervals for WBS
#' d = 10 # parameter for scaled PCA algorithm
#' delta = 5
#' intervals = WBS.intervals(M = M, lower = 1, upper = ncol(data_mat))
#' WBS_result = WBS.nonpar.RDPG(data_mat, lowerdiag = TRUE, d, 
#'              Alpha = intervals$Alpha, Beta = intervals$Beta, delta)
WBS.nonpar.RDPG = function(data_mat, lowerdiag = FALSE, d, Alpha, Beta, delta){
  obs_num = ncol(data_mat)
  if(lowerdiag == TRUE){
    n = 1/2 + sqrt(2*nrow(data_mat) + 1/4) #obtain the number of nodes
    data_mat = apply(data_mat, MARGIN = 2, function(x) lowertri2mat(x, n, diag = FALSE))
  }else{
    n = sqrt(nrow(data_mat))
  }
  xhat = lapply(1:obs_num, function(i){scaledPCA(matrix(data_mat[,i], n, n),d)})
  Y_mat = matrix(0, floor(n/2), obs_num)
  for(t in 1:obs_num){
    phat = drop(xhat[[t]] %*% t(xhat[[t]]))
    ind = sample(1:n, n, replace = FALSE)
    #aux = phat[ind[2*(1:floor(n/2))], ind[2*(1:floor(n/2)) -1]  ]
    for(i in 1:floor(n/2)){
      Y_mat[i,t] = phat[ind[2*i], ind[2*i-1]]
    }
  }
  result = WBS.uni.nonpar(Y_mat, s = 1, e = obs_num, Alpha, Beta, N = rep(nrow(Y_mat), obs_num), delta)
  return(result)
}


#' @title Change points detection for dependent dynamic random dot product graph models.
#' @description Perform Change points detection for dependent dynamic random dot product graph models. The tuning parameter tau for WBS is automatically selected based on the BIC-type scores defined in Equation (2.4) in Zou et al. (2014).
#' @param BS_object A "BS" object produced by \code{WBS.nonpar.RDPG}.
#' @param data_mat  A \code{numeric} matrix of observations with horizontal axis being time, and vertical axis being vectorized adjacency matrix.
#' @param lowerdiag A \code{logic} scalar. TRUE, if each row of data_mat is the vectorization of the strictly lower diagonal elements in an adjacency matrix. FALSE, if each row of data_mat is the vectorization of all elements in an adjacency matrix.
#' @param d         A \code{numeric} scalar of the number of leading singular values of an adjacency matrix considered in the scaled PCA algorithm.
#' @return  A \code{numeric} vector of estimated change points.
#' @export
#' @author Oscar Hernan Madrid Padilla & Haotian Xu
#' @references Padilla, Yu and Priebe (2019) <arxiv:1911.07494>.
#' @seealso \code{\link{WBS.nonpar.RDPG}}.
#' @examples
#' ### generate data 
#' p = 20 # number of nodes
#' n = 50 # sample size for each segment
#' lat_dim_num = 5 # number of latent dimensions
#' set.seed(1)
#' x_mat = matrix(runif(p*lat_dim_num), nrow = p, ncol = lat_dim_num)
#' x_tilde_mat = matrix(runif(p*lat_dim_num), nrow = p, ncol = lat_dim_num)
#' y_mat = rbind(x_tilde_mat[1:floor(p/4),], x_mat[(floor(p/4)+1):p,])
#' rdpg1 = simu.RDPG(x_mat, n, symm = TRUE, self = FALSE)
#' rdpg2 = simu.RDPG(y_mat, n, symm = TRUE, self = FALSE)
#' data1_mat = rdpg1$obs_mat
#' data2_mat = rdpg2$obs_mat
#' data_mat = cbind(data1_mat, data2_mat)
#' ### detect change points
#' M = 20 # number of random intervals for WBS
#' d = 10 # parameter for scaled PCA algorithm
#' delta = 5
#' intervals = WBS.intervals(M = M, lower = 1, upper = ncol(data_mat))
#' WBS_result = WBS.nonpar.RDPG(data_mat, lowerdiag = TRUE, d, 
#'              Alpha = intervals$Alpha, Beta = intervals$Beta, delta)
#' cpt_hat = tuneBSnonparRDPG(WBS_result, data_mat, lowerdiag = TRUE, d)
tuneBSnonparRDPG = function(BS_object, data_mat, lowerdiag = FALSE, d){
  UseMethod("tuneBSnonparRDPG", BS_object)
}

#' @export
tuneBSnonparRDPG.BS = function(BS_object, data_mat, lowerdiag = FALSE, d){
  obs_num = ncol(data_mat)
  if(lowerdiag == TRUE){
    n = 1/2 + sqrt(2*nrow(data_mat) + 1/4) #obtain the number of nodes
    data_mat = apply(data_mat, MARGIN = 2, function(x) lowertri2mat(x, n, diag = FALSE))
  }else{
    n = sqrt(nrow(data_mat))
  }
  xhat = lapply(1:obs_num, function(i){scaledPCA(matrix(data_mat[,i], n, n),d)})
  Y_mat = matrix(0, floor(n/2), obs_num)
  for(t in 1:obs_num){
    phat = drop(xhat[[t]] %*% t(xhat[[t]]))
    ind = sample(1:n, n, replace = FALSE)
    #aux = phat[ind[2*(1:floor(n/2))], ind[2*(1:floor(n/2)) -1]  ]
    for(i in 1:floor(n/2)){
      Y_mat[i,t] = phat[ind[2*i], ind[2*i-1]]
    }
  }
  Dval = BS_object$Dval
  aux = sort(Dval, decreasing = TRUE)
  tau_grid = rev(aux[1:50]-10^{-4})
  tau_grid = tau_grid[which(is.na(tau_grid)==FALSE)]
  tau_grid = c(tau_grid,10) 
  S = c()
  for(j in 1:length(tau_grid)){
    aux = unlist(thresholdBS(BS_object, tau_grid[j])$cpt_hat[,1])
    if(length(aux) == 0){
      break
    }
    S[[j]] = sort(aux)
  }
  S = unique(S)
  score = rep(0, length(S)+1)
  for(i in 1:(n/2)){
    for(j in 1:length(S)){
      score[j] = score[j]+ BICtype.obj(as.matrix(Y_mat[i,]), S[[j]])
    }
    score[j+1] = score[j+1]+ BICtype.obj(as.matrix(Y_mat[i,]), NULL)
  }
  best_ind = which.min(score)
  if(best_ind == length(score)){
    return(NULL)
  }
  return(S[[best_ind]])
}


#' @title scaled PCA
#' @noRd
scaledPCA = function(A_mat, d){
  temp = svd(A_mat)
  xhat_mat =  temp$u[,1:d] %*% diag(sqrt(temp$d[1:d]))
  return(xhat_mat)
}


#' @title Internal Function: Compute BIC-type score for each row of Y_mat.
#' @noRd
BICtype.obj = function(y, S){
  K = length(S)
  obs_num = dim(y)[1]
  cost = 0
  sorted_y = sort(y)
  for(k in 0:(K)){
    if(k == 0){
      s = 1
      if(K > 0){
        e =  S[1]
      }
      if(K == 0){
        e = obs_num
      }
    }
    ########
    ##3
    if(k == K && k > 0){
      s = S[K]+1
      e = obs_num
    }
    ##############3
    if(0 < k && k < K){
      s = S[k]+1
      e = S[k+1]
    }
    ######################33
    aux = as.vector(y[s:e,])
    aux = aux[which(is.na(aux)==FALSE)]
    
    Fhat = ecdf(aux)
    F_hat_z = Fhat(sorted_y)
    
    ny = length(sorted_y)
    a = 1:ny
    a = a*(ny-a)
    
    ind = which(F_hat_z==0)
    if(length(ind)>0){
      F_hat_z = F_hat_z[-ind]
      a = a[-ind] 
    }
    ind = which(F_hat_z==1)
    if(length(ind)>0){
      F_hat_z = F_hat_z[-ind]
      a = a[-ind] 
    }
    cost = cost + ny*(e-s+1)*sum((F_hat_z*log(F_hat_z)+(1-F_hat_z)*log(1-F_hat_z))/a) 
  }
  cost = -cost + 0.2*length(S)*(log(ny))^{2.1}
  return(cost)
}
