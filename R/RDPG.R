



# function to simulate RDPG model
#' @export
#' @example
#' p = 10 # dimension of the network
#' block_num = 4 # number of blocks
#' d = 5 # number of leading eigenvalues considered for CUSUM
#' block_size = floor(p/block_num)
#' # generate P_mat
#' P_mat = matrix(0.3, p, p)
#' for(i in 1:(block_num-1)){
#'   P_mat[(1+(i-1)*block_size):(i*block_size), (1+(i-1)*block_size):(i*block_size)] = 0.5
#' }
#' P_mat[(1+(block_num-1)*block_size):p, (1+(block_num-1)*block_size):p] = 0.5
#' diag(P_mat) = 0
#' # generate Q_mat
#' Q_mat = matrix(0.2, p, p)
#' for(i in 1:(block_num-1)){
#'   Q_mat[(1+(i-1)*block_size):(i*block_size), (1+(i-1)*block_size):(i*block_size)] = 0.45
#' }
#' Q_mat[(1+(block_num-1)*block_size):p, (1+(block_num-1)*block_size):p] = 0.45
#' diag(Q_mat) = 0
#' data1_array = simu.RDPG.SBM(P_mat, rho, n, symm = TRUE, self = FALSE, Azero_mat = NULL)
#' data2_array = simu.RDPG.SBM(Q_mat, rho, n, symm = TRUE, self = FALSE, Azero_mat = data1_array[,,n])
#' data3_array = simu.RDPG.SBM(P_mat, rho, n, symm = TRUE, self = FALSE, Azero_mat = data2_array[,,n])
#' library(abind)
#' data_array = abind(data1_array, data2_array, data3_array, along = 3)
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


#scaled PCA
# assume A_mat is symmetric
scaledPCA = function(A_mat, d){
  temp = svd(A_mat)
  xhat_mat =  temp$u[,1:d] %*% diag(sqrt(temp$d[1:d]))
  return(xhat_mat)
}


#function(data_array)

xhat_list = lapply(1:dim(data_array)[3], function(t) scaledPCA(data_array[,,t], d))

yhat_mat = do.call(cbind, lapply(1:length(xhat_list), function(t) diag(tcrossprod(xhat_list[[t]])[1:floor(p/2), (1+floor(p/2)):(2*floor(p/2))])))

NWBS_estimator = function(y, grid, Ntau =15, gam, N, alpha, beta){
  n = max(N)
  T = dim(y)[1]
  #grid = 
  #temp1 = BS_path(y, gam,1,T,0,NULL,NULL,1, N)  
  temp1 = new_WBS(y, gam,1,T,0,NULL,NULL,1,alpha,beta,N)
  #temp1 =  new_WBS(y[, (1+n/2):n], gam,1,T,0,NULL,NULL,1,alpha,beta,N)
  #n = max(N/2)
  #new_WBS(y1, gam=1,1,T,0,NULL,NULL,1, N)  
  Dval = temp1$Dval
  p1 =  parent(temp1)
  aux = sort(Dval,decreasing = TRUE)
  tau_grid = rev(aux[1:50]-10^{-4})
  tau_grid =  tau_grid[which(is.na(tau_grid)==FALSE)] ### *
  tau_grid = c(tau_grid,10) 
  
  S =  c() 
  for( j in 1:length(tau_grid))
  {
    aux = new_BS_threshold(temp1,tau_grid[j],p1)
    
    if(length(aux)==0)  ##*
      break;###*
    
    S[[j]] = sort(aux)
  }
  
  T= dim(y)[1]
  S = unique(S)
  score  = rep(0,length(S)+1)
  #n = dim(y)[2]
  
  for(i in 1:n)
  {
    for(j in 1:length(S))
    {
      
      score[j] = score[j]+ L_obj_v4(as.matrix(y[,i]),as.matrix(y[,i]),S[[j]])
    }
    score[j+1] = score[j+1]+ L_obj_v4(as.matrix(y[,i]),as.matrix(y[,i]),NULL)
    
  }
  
  
  best_ind = which.min(score)
  
  if(best_ind==length(score))
  {
    return(NULL)
  }
  return( S[[best_ind]])
}

L_obj_v4 = function(y1,y2,S){
  K =  length(S)
  # # 
  T =  dim(y1)[1]
  val = rep(0,length(grid)) 
  cost = 0
  sorted_y2 = sort(y2)
  # 
  for(k in 0:(K)){
    if(k== 0){
      s =  1
      if(K>0){
        e =  S[1]
      }
      if(K==0){
        e = T
      }
    }
    ########
    ##3
    if(k   == K && k > 0){
      s =  S[K]+1
      e =  T
    }
    ##############3
    if( 0 < k &&  k < K){
      s= S[k]+1
      e= S[k+1]
    }
    ######################33
    
    aux =  as.vector(y2[s:e,])
    aux =  aux[which(is.na(aux)==FALSE)]
    
    Fhat = ecdf(  aux )
    F_hat_z = Fhat(sorted_y2)
    
    ny =  length(sorted_y2)
    a = 1:ny
    a = a*(ny-a)
    
    ind =  which(F_hat_z==0)
    if(length(ind)>0){
      F_hat_z = F_hat_z[-ind]
      a = a[-ind] 
    }
    ind =  which(F_hat_z==1)
    if(length(ind)>0){
      F_hat_z = F_hat_z[-ind]
      a = a[-ind] 
    }
    cost = cost +   ny*(e-s+1)*sum( (F_hat_z*log(F_hat_z)+(1-F_hat_z)*log(1-F_hat_z))/a     ) 
  }
  #0.5*(1/log(10))
  cost = -cost +  0.2*length(S)*(log(ny))^{2.1}
  return(cost)
}



