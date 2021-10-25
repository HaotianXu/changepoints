



# function to simulate RDPG model
# @export
# n = 40 # sample size for each segment
# p = 10 # dimension of the network
# block_num = 4 # number of blocks
# d = 5 # number of leading eigenvalues considered for CUSUM
# block_size = floor(p/block_num)
# # generate P_mat
# P_mat = matrix(0.3, p, p)
# for(i in 1:(block_num-1)){
#   P_mat[(1+(i-1)*block_size):(i*block_size), (1+(i-1)*block_size):(i*block_size)] = 0.5
# }
# P_mat[(1+(block_num-1)*block_size):p, (1+(block_num-1)*block_size):p] = 0.5
# diag(P_mat) = 0
# # generate Q_mat
# Q_mat = matrix(0.2, p, p)
# for(i in 1:(block_num-1)){
#   Q_mat[(1+(i-1)*block_size):(i*block_size), (1+(i-1)*block_size):(i*block_size)] = 0.45
# }
# Q_mat[(1+(block_num-1)*block_size):p, (1+(block_num-1)*block_size):p] = 0.45
# diag(Q_mat) = 0
# data1_array = simu.RDPG.SBM(P_mat, rho, n, symm = TRUE, self = FALSE, Azero_mat = NULL)
# data2_array = simu.RDPG.SBM(Q_mat, rho, n, symm = TRUE, self = FALSE, Azero_mat = data1_array[,,n])
# data3_array = simu.RDPG.SBM(P_mat, rho, n, symm = TRUE, self = FALSE, Azero_mat = data2_array[,,n])
# library(abind)
# data_array = abind(data1_array, data2_array, data3_array, along = 3)
# simu.RDPG.SBM = function(P_mat, rho, n, symm = TRUE, self = FALSE, Azero_mat = NULL){
#   p = ncol(P_mat)
#   obs_array = array(NA, dim = c(p,p,n))
#   if(is.null(Azero_mat)){
#     A_mat = matrix(rbinom(matrix(1,p,p),matrix(1,p,p),P_mat),p,p)
#   }else if(any(dim(Azero_mat) != c(p,p))){
#     stop("If the observation (adjacency matrix) at time 0 (Azero) is specified, it should be a pxp matrix.")
#   }else{
#     A_mat = Azero_mat
#   }
#   if((symm == TRUE) & (self == FALSE)){
#     diag(A_mat) = 0
#     A_mat[upper.tri(A_mat)] = t(A_mat)[upper.tri(A_mat)]
#     obs_array[,,1] = A_mat
#     for(t in 2:n){
#       aux1 = P_mat + (1-P_mat)*rho
#       aux2 = P_mat*(1-rho)
#       aux1 = matrix(rbinom(matrix(1,p,p),matrix(1,p,p),aux1),p,p)
#       aux2 = matrix(rbinom(matrix(1,p,p),matrix(1,p,p),aux2),p,p)
#       A_mat = aux1*A_mat + aux2*(1-A_mat)
#       A_mat[upper.tri(A_mat)] = t(A_mat)[upper.tri(A_mat)]
#       obs_array[,,t] = A_mat
#     }
#   }else if((symm == TRUE) & (self == TRUE)){
#     A_mat[upper.tri(A_mat)] = t(A_mat)[upper.tri(A_mat)]
#     obs_array[,,1] = A_mat
#     for(t in 2:n){
#       aux1 = P_mat + (1-P_mat)*rho
#       aux2 = P_mat*(1-rho)
#       aux1 = matrix(rbinom(matrix(1,p,p),matrix(1,p,p),aux1),p,p)
#       aux2 = matrix(rbinom(matrix(1,p,p),matrix(1,p,p),aux2),p,p)
#       A_mat = aux1*A_mat + aux2*(1-A_mat)
#       A_mat[upper.tri(A_mat)] = t(A_mat)[upper.tri(A_mat)]
#       obs_array[,,t] = A_mat
#     }
#   }else if((symm == FALSE) & (self == FALSE)){
#     diag(A_mat) = 0
#     obs_array[,,1] = A_mat
#     for(t in 2:n){
#       aux1 = P_mat + (1-P_mat)*rho
#       aux2 = P_mat*(1-rho)
#       aux1 = matrix(rbinom(matrix(1,p,p),matrix(1,p,p),aux1),p,p)
#       aux2 = matrix(rbinom(matrix(1,p,p),matrix(1,p,p),aux2),p,p)
#       A_mat = aux1*A_mat + aux2*(1-A_mat)
#       obs_array[,,t] = A_mat
#     }
#   }else{
#     obs_array[,,1] = A_mat
#     for(t in 2:n){
#       aux1 = P_mat + (1-P_mat)*rho
#       aux2 = P_mat*(1-rho)
#       aux1 = matrix(rbinom(matrix(1,p,p),matrix(1,p,p),aux1),p,p)
#       aux2 = matrix(rbinom(matrix(1,p,p),matrix(1,p,p),aux2),p,p)
#       A_mat = aux1*A_mat + aux2*(1-A_mat)
#       obs_array[,,t] = A_mat
#     }
#   }
#   return(obs_array)
# }
# 
# # data1_array = simu.RDPG.SBM(P_mat, rho, n, symm = TRUE, self = FALSE, Azero_mat = NULL)
# # data2_array = simu.RDPG.SBM(Q_mat, rho, n, symm = TRUE, self = FALSE, Azero_mat = data1_array[,,n])
# # data3_array = simu.RDPG.SBM(P_mat, rho, n, symm = TRUE, self = FALSE, Azero_mat = data2_array[,,n])
# # library(abind)
# # data = abind(data1_array, data2_array, data3_array, along = 3)
# # 
# # data[,t] = as.vector(A)
# # temp = svd(A)
# # xhat[t,,] =  temp$u[,1:d] %*%  diag( sqrt(temp$d[1:d]) ) 
# 
# 
# 
# #
# p = 100 # number of nodes
# lat_dim_num = 5
# x_mat = matrix(runif(p*lat_dim_num), nrow = p, ncol = lat_dim_num)
# x_tilde_mat = matrix(runif(p*lat_dim_num), nrow = p, ncol = lat_dim_num)
# 
# 
# x_l2_vec = apply(x_mat, 1, function(x) (sum(x^2))^(1/2))
# x_l2prod_mat = x_l2_vec %*% t(x_l2_vec)
# x_cross_mat = x_mat %*% t(x_mat)
# ratio_x_mat = x_cross_mat/x_l2prod_mat
# diag(ratio_x_mat) = 1
# 
# y_mat = rbind(x_tilde_mat[1:floor(p/4),], x_mat[(floor(p/4)+1):p,])
# y_l2_vec = apply(y_mat, 1, function(x) (sum(x^2))^(1/2))
# y_l2prod_mat = y_l2_vec %*% t(y_l2_vec)
# y_cross_mat = y_mat %*% t(y_mat)
# ratio_y_mat = y_cross_mat/y_l2prod_mat
# diag(ratio_y_mat) = 1
# 
# A1_mat = matrix(rbinom(matrix(1,p,p), matrix(1,p,p), ratio_x_mat), p, p)
# A2_mat = matrix(rbinom(matrix(1,p,p), matrix(1,p,p), ratio_y_mat), p, p)















#function(data_array)

# xhat_list = lapply(1:dim(data_array)[3], function(t) scaledPCA(data_array[,,t], d))
# yhat_mat = do.call(cbind, lapply(1:length(xhat_list), function(t) diag(tcrossprod(xhat_list[[t]])[1:floor(p/2), (1+floor(p/2)):(2*floor(p/2))])))
# M = 120
# intervals = WBS.intervals(M = M, lower = 1, upper = dim(data_array)[3])
# Alpha = intervals$Alpha
# Beta = intervals$Beta
# temp = NWBS(yhat_mat, 1, 120, Alpha, Beta, N = rep(5, 120), 5)
# threshold.BS(temp, 2)




#scaled PCA
#' @noRd
scaledPCA = function(A_mat, d){
  temp = svd(A_mat)
  xhat_mat =  temp$u[,1:d] %*% diag(sqrt(temp$d[1:d]))
  return(xhat_mat)
}


#' @title Change point detection for dependent dynamic random dot product graph models.
#' @description Perform Change point detection for dependent dynamic random dot product graph models.
#' @param data_mat  A \code{numeric} matrix of observations with horizontal axis being time, and vertical axis being vectorized adjacency matrix.
#' @param d         A \code{numeric} scalar used in the scaledPCA algorithm, which corresponds to the number of leading singular values of an adjacency matrix considered.
#' @param Alpha     A \code{integer} vector of starting indices of random intervals.
#' @param Beta      A \code{integer} vector of ending indices of random intervals.
#' @param delta     A positive \code{integer} scalar of minimum spacing.
#' @param ...      Additional arguments.
#' @return  A \code{numeric} vector of estimated changepoint locations.
#' @export
#' @author Oscar Hernan Madrid Padilla, Haotian Xu
#' @examples
#' ### generate data 
#' d = 10
#' n = 100
#' M = 120
#' delta = 5
#' obs_num = 150
#' rho_a = 0.9
#' v = c(floor(obs_num/3)+1, 2*floor(obs_num/3)+1)
#' data_mat = matrix(0, n^2, obs_num)
#' for(t in 1:obs_num){
#'   if(t == 1 || t == v[2]+1){
#'     P = matrix(0.3,n,n)
#'     P[1:floor(n/4), 1:floor(n/4)] = 0.5
#'     P[(1+floor(n/4)):(2*floor(n/4)),(1+floor(n/4)):(2*floor(n/4))] = 0.5
#'     P[(1+2*floor(n/4)):(3*floor(n/4)),(1+2*floor(n/4)):(3*floor(n/4))] = 0.5
#'     P[(1+3*floor(n/4)):n,(1+3*floor(n/4)):n] = 0.5
#'     diag(P) = 0
#'     A = matrix(rbinom(matrix(1,n,n),matrix(1,n,n),P),n,n)
#'     aux = drop(A)
#'     aux[lower.tri(aux)] = t(aux)[lower.tri(aux)]
#'     diag(aux) = 0
#'     data_mat[,t] = drop(matrix(aux,n^2,1))
#'   }
#'   if((t > 1 && t <= v[1]) || (t > v[2]+1)){
#'     aux1 = P + (1-P)*rho_a
#'     aux2 = P*(1-rho_a)
#'     aux1 = matrix(rbinom(matrix(1,n,n),matrix(1,n,n),aux1),n,n)
#'     aux2 = matrix(rbinom(matrix(1,n,n),matrix(1,n,n),aux2),n,n)
#'     A =  aux1*A + aux2*(1-A)
#'     aux = drop(A)
#'     aux[lower.tri(aux)] = t(aux)[lower.tri(aux)]  
#'     diag(aux) = 0
#'     data_mat[,t] = drop(matrix(aux,n^2,1))
#'   }
#'   if(t == v[1]+1){
#'     Q = matrix(0.2,n,n)
#'     Q[1:floor(n/4), 1:floor(n/4)] = 0.45
#'     Q[(1+floor(n/4)):(2*floor(n/4)),(1+floor(n/4)):(2*floor(n/4)) ] = 0.45
#'     Q[(1+2*floor(n/4)):(3*floor(n/4)),(1+2*floor(n/4)):(3*floor(n/4)) ] = 0.45
#'     Q[(1+3*floor(n/4)):n,(1+3*floor(n/4)):n ] = 0.45
#'     diag(Q) = 0
#'     A = matrix(rbinom(matrix(1,n,n),matrix(1,n,n),Q),n,n)
#'     aux = drop(A)
#'     aux[lower.tri(aux)] = t(aux)[lower.tri(aux)]  
#'     diag(aux) = 0
#'     data_mat[,t] = drop(matrix(aux,n^2,1))
#'   }
#'   if(t > v[1]+1 && t <= v[2]){
#'     aux1 = Q + (1-Q)*rho_a
#'     aux2 = Q*(1-rho_a)
#'     aux1 = matrix(rbinom(matrix(1,n,n),matrix(1,n,n),aux1),n,n)
#'     aux2 = matrix(rbinom(matrix(1,n,n),matrix(1,n,n),aux2),n,n)
#'     A = aux1*A + aux2*(1-A)
#'     aux = drop(A)
#'     aux[lower.tri(aux)] = t(aux)[lower.tri(aux)]  
#'     diag(aux) = 0
#'     data_mat[,t] = drop(matrix(aux,n^2,1))
#'   }
#' }
#' intervals = WBS.intervals(M = M, lower = 1, upper = obs_num)
#' cpt_hat = NonPar.RDPG.CPD(data_mat, d, Alpha = intervals$Alpha, Beta = intervals$Beta, delta)
#' 
NonPar.RDPG.CPD = function(data_mat, d, Alpha, Beta, delta){
  obs_num = ncol(data_mat)
  n = sqrt(nrow(data_mat))
  grid = seq(min(Y_mat),max(Y_mat), length = length(Alpha))
  
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
  temp1 = NWBS(Y_mat, s = 1, e = obs_num, Alpha, Beta, N = rep(nrow(Y_mat), obs_num), delta)
  Dval = temp1$Dval
  aux = sort(Dval, decreasing = TRUE)
  tau_grid = rev(aux[1:50]-10^{-4})
  tau_grid = tau_grid[which(is.na(tau_grid)==FALSE)]
  tau_grid = c(tau_grid,10) 
  S = c()
  for(j in 1:length(tau_grid)){
    aux = threshold.BS(temp1, tau_grid[j])$change_points[,1]
    if(length(aux) == 0)
      break;
    S[[j]] = sort(aux)
  }
  S = unique(S)
  score = rep(0, length(S)+1)
  for(i in 1:(n/2)){
    for(j in 1:length(S)){
      score[j] = score[j]+ BICtype.obj(as.matrix(Y_mat[i,]), S[[j]], grid)
    }
    score[j+1] = score[j+1]+ BICtype.obj(as.matrix(Y_mat[i,]), NULL, grid)
  }
  best_ind = which.min(score)
  if(best_ind == length(score)){
    return(NULL)
  }
  return(S[[best_ind]])
}


#' @title Internal Function: Compute BIC-type score for each row of Y_mat.
#' @noRd
BICtype.obj = function(y, S, grid){
  K = length(S)
  obs_num = dim(y)[1]
  val = rep(0, length(grid)) 
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
