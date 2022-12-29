#' @title Internal Function: Compute the CUSUM statistic based on KS distance (multivariate).
#' @param Y         A \code{numeric} matrix of observations with with horizontal axis being time, and vertical axis being dimension.
#' @param s         A \code{integer} scalar of starting index.
#' @param e         A \code{integer} scalar of ending index.
#' @param t         A \code{integer} scalar of splitting index.
#' @param h         A \code{integer} scalar of bandwidth parameter.
#' @return  A \code{numeric} scalar of the CUSUM statistic based on KS distance.
#' @noRd
CUSUM.KS.multivariate = function(Y, W, s, e, t, h){
  p = dim(Y)[1]
  n_st = t - s + 1
  n_se = e - s + 1
  n_te = e - t
  aux = Y[,s:t]
  temp1 = kde.eval(t(aux), eval.points = t(W), H = h*diag(p))
  aux = Y[,(t+1):e]
  temp2 = kde.eval(t(aux), eval.points = t(W), H = h*diag(p))
  result = sqrt(n_st * n_te / n_se) * max(abs(temp1 - temp2))
  return(result)
}


#' @title Wild binary segmentation for multivariate nonparametric change points detection.
#' @description Perform wild binary segmentation for multivariate nonparametric change points detection.
#' @param Y         A \code{numeric} matrix of observations with with horizontal axis being time, and vertical axis being dimension.
#' @param W         A copy of the matrix Y, it can be Y itself.
#' @param s         A \code{integer} scalar of starting index.
#' @param e         A \code{integer} scalar of ending index.
#' @param Alpha     A \code{integer} vector of starting indices of random intervals.
#' @param Beta      A \code{integer} vector of ending indices of random intervals.
#' @param h         A \code{numeric} scalar of bandwidth parameter.
#' @param delta     A \code{integer} scalar of minimum spacing.
#' @param level     Should be fixed as 0.
#' @return   An object of \code{\link[base]{class}} "BS", which is a \code{list} with the following structure:
#'  \item{S}{A vector of estimated change points (sorted in strictly increasing order).}
#'  \item{Dval}{A vector of values of CUSUM statistic based on KS distance.}
#'  \item{Level}{A vector representing the levels at which each change point is detected.}
#'  \item{Parent}{A matrix with the starting indices on the first row and the ending indices on the second row.}
#' @export
#' @author Oscar Hernan Madrid Padilla & Haotian Xu
#' @references Padilla, Yu, Wang and Rinaldo (2019) <arxiv:1910.13289>.
#' @seealso \code{\link{thresholdBS}} for obtain change points estimation, \code{\link{tuneBSmultinonpar}} for a tuning version.
#' @examples
#' n = 70
#' v = c(floor(n/3), 2*floor(n/3)) # location of change points
#' p = 4
#' Y = matrix(0, p, n) # matrix for data
#' mu0 = rep(0, p) # mean of the data
#' mu1 = rep(0, p)
#' mu1[1:floor(p/2)] = 2
#' Sigma0 = diag(p) #Covariance matrices of the data
#' Sigma1 = diag(p)*2
#' # Generate data
#' for(t in 1:n){
#'   if(t < v[1] || t > v[2]){
#'      Y[,t] = MASS::mvrnorm(n = 1, mu0, Sigma0)
#'   }
#'   if(t >= v[1] && t < v[2]){
#'      Y[,t] = MASS::mvrnorm(n = 1, mu1, Sigma1)
#'   }
#' }## close for generate data
#' M = 10
#' intervals = WBS.intervals(M = M, lower = 1, upper = ncol(Y)) #Random intervals
#' K_max = 30
#' h = 5*(K_max*log(n)/n)^{1/p} # bandwith
#' temp = WBS.multi.nonpar(Y, Y, 1, ncol(Y), intervals$Alpha, intervals$Beta, h, delta = 10)
#' result = thresholdBS(temp, median(temp$Dval))
WBS.multi.nonpar = function(Y, W, s, e, Alpha, Beta, h, delta, level = 0){
  print(paste0("WBS at level: ", level))
  Alpha_new = pmax(Alpha, s)
  Beta_new = pmin(Beta, e)
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
      s_star = Alpha_new[m] + delta
      e_star = Beta_new[m] - delta
      temp = rep(0, e_star - s_star + 1)
      for(t in s_star:e_star){
        temp[t-s_star+1] = CUSUM.KS.multivariate(Y, W, Alpha_new[m], Beta_new[m], t, h)
      }
      best_value = max(temp)
      best_t = which.max(temp) + s_star - 1
      a[m] = best_value
      b[m] = best_t
    }
    m_star = which.max(a)
  }
  temp1 = WBS.multi.nonpar(Y, W, s, b[m_star]-1, Alpha, Beta, h, delta, level)
  temp2 = WBS.multi.nonpar(Y, W, b[m_star], e, Alpha, Beta, h, delta, level)
  S = c(temp1$S, b[m_star], temp2$S)
  Dval = c(temp1$Dval, a[m_star], temp2$Dval)
  Level = c(temp1$Level, level, temp2$Level)
  Parent = cbind(temp1$Parent, parent, temp2$Parent)
  result = list(S = S, Dval = Dval, Level = Level, Parent = Parent)
  class(result) = "BS"
  return(result)
}



#' A function to compute change points based on the multivariate nonparametic method with tuning parameter selected by FDR control.
#' @param BS_object A "BS" object produced by \code{WBS.multi.nonpar}.
#' @param Y         A \code{numeric} matrix of observations with with horizontal axis being time, and vertical axis being dimension.
#' @return A vector of estimated change points.
#' @author Oscar Hernan Madrid Padilla & Haotian Xu
#' @references Padilla, Yu, Wang and Rinaldo (2019) <arxiv:1910.13289>.
#' @seealso \code{\link{WBS.multi.nonpar}}.
#' @export
#' @examples 
#' n = 70
#' v = c(floor(n/3), 2*floor(n/3)) # location of change points
#' p = 4
#' Y = matrix(0, p, n) # matrix for data
#' mu0 = rep(0, p) # mean of the data
#' mu1 = rep(0, p)
#' mu1[1:floor(p/2)] = 2
#' Sigma0 = diag(p) #Covariance matrices of the data
#' Sigma1 = diag(p)*2
#' # Generate data
#' for(t in 1:n){
#'   if(t < v[1] || t > v[2]){
#'      Y[,t] = MASS::mvrnorm(n = 1, mu0, Sigma0)
#'   }
#'   if(t >= v[1] && t < v[2]){
#'      Y[,t] = MASS::mvrnorm(n = 1, mu1, Sigma1)
#'   }
#' }## close for generate data
#' M = 8
#' intervals = WBS.intervals(M = M, lower = 1, upper = ncol(Y)) #Random intervals
#' K_max = 30
#' h = 5*(K_max*log(n)/n)^{1/p} # bandwith
#' temp = WBS.multi.nonpar(Y, Y, 1, ncol(Y), intervals$Alpha, intervals$Beta, h, delta = 10)
#' S = tuneBSmultinonpar(temp, Y)
tuneBSmultinonpar = function(BS_object, Y){
  UseMethod("tuneBSmultinonpar", BS_object)
}

#' @export
tuneBSmultinonpar.BS = function(BS_object, Y){
  p = nrow(Y)
  obs_num = ncol(Y)
  Dval = BS_object$Dval
  aux = sort(Dval, decreasing = TRUE)
  len_tau = 30
  tau_grid = rev(aux[1:min(len_tau,length(Dval))]) - 10^{-30}
  B_list = c()
  for(j in 1:length(tau_grid)){
    aux = thresholdBS(BS_object, tau_grid[j])$cpt_hat[,1]
    if(length(aux) == 0){
      break
    }
    B_list[[j]] = sort(aux)
  }
  B_list = unique(B_list)
  if(length(B_list) == 0){
    return(NULL)
  }
  if(length(B_list[[1]]) == 0){
    return(B_list[[1]])
  }
  lambda = 1.8^2
  v_n = 100
  min_pval = 10
  for(j in 1:length(B_list)){
    B2 = B_list[[j]]
    if(j < length(B_list)){
      B1 = B_list[[j+1]]
    }else if(j == length(B_list)){
      B1 = NULL
    }
    temp = setdiff(B2,B1)
    for(l in 1:length(temp)){
      eta = temp[l]
      if(length(B1) == 0){
        eta1 = 1
        eta2 = obs_num
      }else if(length(B1) > 0){
        for(k in 1:length(B1)){
          if(B1[k] > eta){
            break
          }
        }
        if(B1[k] > eta){
          eta2 = B1[k]
          if(k == 1)
            eta1 = 1
          if(k > 1)
            eta1 = B1[k-1]+1
        }
        if(B1[k] < eta){
          eta1 = B1[k] + 1
          eta2 = obs_num
        }
      }######## if length(B1) > 0
      pval = rep(0, length(v_n))
      for(ind_v in 1:v_n){
        vec = runif(p)
        vec = vec/sqrt(sum(vec^2))
        if(ind_v < p+1){
          vec = rep(0, p)
          vec[ind_v] = 1
        }
        aux = CUSUM.KS(t(t(Y)%*%vec), eta1, eta2, eta, rep(1,obs_num))  
        pval[ind_v] = exp(-2*aux^2)
        st_aux = aux
      }## for ind_v
      pval_adj = p.adjust(pval, method ="fdr")
      if(min_pval > min(pval_adj)){
        min_pval = min(pval_adj)
      }
    }### for l
    if(min_pval < 0.0005){
      break
    }
  }##### j 
  if(min_pval > .1){
    est = NULL
    return(est)
  }
  est = B_list[[j]]
  return(est)
}
