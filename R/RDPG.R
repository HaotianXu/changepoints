
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


# function to simulate RDPG model
#' @export
simu.RDPG.SBM = function(P_mat, rho, n, symm = TRUE, self = FALSE, Azero_mat = NULL){
  p = ncol(P_mat)
  obs_array = array(NA, dim = c(p,p,n))
  if(is.null(Azero_mat)){
    A_mat = matrix(rbinom(matrix(1,p,p),matrix(1,p,p),P_mat),p,p)
  }else if(any(dim(Azero_mat) != c(p,p))){
    stop("If the observation (adjacency matrix) at time 0 (Azero) is specified, it should be a pxp matrix.")
  }else{
    A_mat = Azero_mat
  }
  if((symm == TRUE) & (self == FALSE)){
    diag(A_mat) = 0
    A_mat[upper.tri(A_mat)] = t(A_mat)[upper.tri(A_mat)]
    obs_array[,,1] = A_mat
    for(t in 2:n){
      aux1 = P_mat + (1-P_mat)*rho
      aux2 = P_mat*(1-rho)
      aux1 = matrix(rbinom(matrix(1,p,p),matrix(1,p,p),aux1),p,p)
      aux2 = matrix(rbinom(matrix(1,p,p),matrix(1,p,p),aux2),p,p)
      A_mat = aux1*A_mat + aux2*(1-A_mat)
      A_mat[upper.tri(A_mat)] = t(A_mat)[upper.tri(A_mat)]
      obs_array[,,t] = A_mat
    }
  }else if((symm == TRUE) & (self == TRUE)){
    A_mat[upper.tri(A_mat)] = t(A_mat)[upper.tri(A_mat)]
    obs_array[,,1] = A_mat
    for(t in 2:n){
      aux1 = P_mat + (1-P_mat)*rho
      aux2 = P_mat*(1-rho)
      aux1 = matrix(rbinom(matrix(1,p,p),matrix(1,p,p),aux1),p,p)
      aux2 = matrix(rbinom(matrix(1,p,p),matrix(1,p,p),aux2),p,p)
      A_mat = aux1*A_mat + aux2*(1-A_mat)
      A_mat[upper.tri(A_mat)] = t(A_mat)[upper.tri(A_mat)]
      obs_array[,,t] = A_mat
    }
  }else if((symm == FALSE) & (self == FALSE)){
    diag(A_mat) = 0
    obs_array[,,1] = A_mat
    for(t in 2:n){
      aux1 = P_mat + (1-P_mat)*rho
      aux2 = P_mat*(1-rho)
      aux1 = matrix(rbinom(matrix(1,p,p),matrix(1,p,p),aux1),p,p)
      aux2 = matrix(rbinom(matrix(1,p,p),matrix(1,p,p),aux2),p,p)
      A_mat = aux1*A_mat + aux2*(1-A_mat)
      obs_array[,,t] = A_mat
    }
  }else{
    obs_array[,,1] = A_mat
    for(t in 2:n){
      aux1 = P_mat + (1-P_mat)*rho
      aux2 = P_mat*(1-rho)
      aux1 = matrix(rbinom(matrix(1,p,p),matrix(1,p,p),aux1),p,p)
      aux2 = matrix(rbinom(matrix(1,p,p),matrix(1,p,p),aux2),p,p)
      A_mat = aux1*A_mat + aux2*(1-A_mat)
      obs_array[,,t] = A_mat
    }
  }
  return(obs_array)
}


# data1_array = simu.RDPG.SBM(P_mat, rho, n, symm = TRUE, self = FALSE, Azero_mat = NULL)
# data2_array = simu.RDPG.SBM(Q_mat, rho, n, symm = TRUE, self = FALSE, Azero_mat = data1_array[,,n])
# data3_array = simu.RDPG.SBM(P_mat, rho, n, symm = TRUE, self = FALSE, Azero_mat = data2_array[,,n])
# library(abind)
# data = abind(data1_array, data2_array, data3_array, along = 3)
# 
# data[,t] = as.vector(A)
# temp = svd(A)
# xhat[t,,] =  temp$u[,1:d] %*%  diag( sqrt(temp$d[1:d]) ) 

