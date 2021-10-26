#' @title Change point detection for dependent dynamic random dot product graph models.
#' @description Perform Change point detection for dependent dynamic random dot product graph models. The tuning parameter tau for WBS is automatically selected based on the BIC-type scores defined in Equation (2.4) in Zou et al. (2014).
#' @param data_mat  A \code{numeric} matrix of observations with horizontal axis being time, and vertical axis being vectorized adjacency matrix.
#' @param d         A \code{numeric} scalar used in the scaledPCA algorithm, which corresponds to the number of leading singular values of an adjacency matrix considered.
#' @param Alpha     A \code{integer} vector of starting indices of random intervals.
#' @param Beta      A \code{integer} vector of ending indices of random intervals.
#' @param delta     A positive \code{integer} scalar of minimum spacing.
#' @param ...      Additional arguments.
#' @return  A \code{numeric} vector of estimated changepoint locations.
#' @export
#' @author Oscar Hernan Madrid Padilla, Haotian Xu
#' @references Zou et al. (2014), Nonparametric maximum likelihood approach to multiple change-point problems,  Ann. Statist. 42(3): 970-1002 (June 2014). DOI: 10.1214/14-AOS1210
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
