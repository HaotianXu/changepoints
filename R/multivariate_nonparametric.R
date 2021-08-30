#' @title Internal Function: Compute the CUSUM statistic based on KS distance (multivariate).
#' @param Y         A \code{numeric} matrix of observations with with horizontal axis being time, and vertical axis being dimension.
#' @param s         A \code{integer} scalar of starting index.
#' @param e         A \code{integer} scalar of ending index.
#' @param t         A \code{integer} scalar of splitting index.
#' @param h         A \code{integer} scalar of bandwidth parameter.
#' @return  A \code{numeric} scalar of the CUSUM statistic based on KS distance.
#' @noRd
CUSUM.KS.multivariate = function(Y, s, e, t, h){
  p = dim(Y)[1]
  n = dim(Y)[2]
  n_st = t - s + 1
  n_se = e - s + 1
  n_te = e - t
  aux = Y[,s:t]
  aux = aux[which(is.na(aux)==FALSE)]
  temp1 = kde(t(aux), gridsize = 30, eval.points = t(Y), H = h*diag(p))
  aux = Y[,(t+1):e]
  aux = aux[which(is.na(aux)==FALSE)]
  temp2 = kde(t(aux), gridsize = 30, eval.points = t(Y), H = h*diag(p))
  result = sqrt(n_st * n_te / n_se) * max(abs(temp1$estimate - temp2$estimate))
  return(result)
}


#' @title Wild binary segmentation for multivariate nonparametric change points detection
#' @description TO DO
#' @param Y         A \code{numeric} matrix of observations with with horizontal axis being time, and vertical axis being dimension.
#' @param s         A \code{integer} scalar of starting index.
#' @param e         A \code{integer} scalar of ending index.
#' @param Alpha     A \code{integer} vector of starting indices of random intervals.
#' @param Beta      A \code{integer} vector of ending indices of random intervals.
#' @param h         A \code{numeric} scalar of bandwidth parameter.
#' @param level     Should be fixed as 0.
#' @param ...       Additional arguments.
#' @return     A \code{list} with the structure:
#' \itemize{
#'  \item S:           A vector of estimated changepoints (sorted in strictly increasing order).
#'  \item Dval:        A vector of values of CUSUM statistic based on KS distance.
#'  \item Level:       A vector representing the levels at which each change point is detected.
#'  \item Parent:      A matrix with the starting indices on the first row and the ending indices on the second row.
#' } 
#' @export
#' @author Haotian Xu
#' @examples
#' new_MWBS
#' 
MNWBS = function(Y, s, e, Alpha, Beta, h, delta = ceiling(max(p, h^(-p))), level = 0, ...){
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
        temp[t-s_star+1] = CUSUM.KS.multivariate(Y, Alpha_new[m], Beta_new[m], t, h)
      }
      best_value = max(temp)
      best_t = which.max(temp) + s_star - 1
      a[m] = best_value
      b[m] = best_t
    }
    m_star = which.max(a)
  }
  temp1 = MNWBS(Y, s, b[m_star]-1, Alpha, Beta, h, level)
  temp2 = MNWBS(Y, b[m_star], e, Alpha, Beta, h, level)
  S = c(temp1$S, b[m_star], temp2$S)
  Dval = c(temp1$Dval, a[m_star], temp2$Dval)
  Level = c(temp1$Level, level, temp2$Level)
  Parent = cbind(temp1$Parent, parent, temp2$Parent)
  result = list(S = S, Dval = Dval, Level = Level, Parent = Parent)
  class(result) = "BS"
  return(result)
}



#' @export
MNWBS.tune = function(Y, W, Alpha, Beta, h, len_tau, delta){
  len_time = ncol(Y)
  temp1 = MNWBS(W, 1, len_time, Alpha, Beta, h, delta, level = 0)  
  Dval = temp1$Dval
  aux = sort(Dval, decreasing = TRUE)
  tau_grid = rev(aux[1:min(len_tau,length(Dval))])
  B_list = c()
  for(j in 1:length(tau_grid)){
    aux = threshold.BS(temp1, tau_grid[j])$change_points[,1]
    if(length(aux) == 0){
      break;
    }
    B_list[[j]] = sort(aux)
  }
  B_list = unique(B_list)
  O_set = B_list[[1]]
  lambda = log(sum(N))/1.5#2.5#1.5#2#2.555#
  v_n = 100
  for(m in 1:(length(B_list)-1)){
    eta = min(setdiff(O_set, B_list[[m+1]]))
    B_set_ext = c(1, B_list[[m+1]], len_time)
    k = which.max(B_list[[m+1]] < eta)
    eta1 = B_list[[m+1]][k]
    eta2 = B_list[[m+1]][k+1]
    pval = rep(0, length(v_n))
    vec_s = matrix(0, v_n, len_time)
    for(ind_v in 1:v_n){
      vec = runif(len_time)
      vec = vec/sqrt(sum(vec^2))
      vec_s[ind_v,] = vec 
      if(ind_v < len_time + 1){
        vec = rep(0, len_time)
        vec[ind_v] = 1
      }
      aux = CUSUM.KS(Y%*%vec, eta1, eta2, eta, rep(1,len_time))  
      pval[ind_v] = exp(-2*aux^2)
      st_aux = aux
    }
    pval_adj = p.adjust(pval, method ="fdr")
    min_pval = min(pval_adj)
    if(min_pval < 0.0005){
      break;
    }
    if(min_pval > 0.1){
      O_set = NULL
      return(O_set)
    }
    O_set = B_list[[m+1]]
  }
  return(O_set)
}
