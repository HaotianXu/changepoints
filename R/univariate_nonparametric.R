#' @title Standard binary segmentation for univariate nonparametric change points detection
#' @description TO DO
#' @param Y         A \code{numeric} matrix of observations with horizontal axis being time, and vertical axis being multiple observations on each time point.
#' @param delta     A strictly \code{integer} scalar of minimum spacing.
#' @param s         A \code{integer} scalar of starting index.
#' @param e         A \code{integer} scalar of ending index.
#' @param N         A \code{integer} vector representing number of multiple observations on each time point.
#' @param level     Should be fixed as 0.
#' @param ...      Additional arguments.
#' @return  A \code{list} with the structure:
#' \itemize{
#'  \item S           A vector of estimated changepoints (sorted in strictly increasing order).
#'  \item Dval        A vector of values of CUSUM statistic based on KS distance.
#'  \item Level       A vector representing the levels at which each change point is detected.
#'  \item ...         Additional parameters.
#' } 
#' @export
#' @author Haotian Xu
#' @examples
#' Y = t(as.matrix(c(rnorm(100, 0, 1), rnorm(100, 0, 10), rnorm(100, 0, 40))))
#' temp = NBS(Y, 5, 1, 300, N)
#' plot.ts(t(Y))
#' points(x = tail(temp$S[order(temp$Dval)],4), y = Y[,tail(temp$S[order(temp$Dval)],4)], col = "red")
NBS = function(Y, delta, s, e, N, level = 0, ...){
  S = NULL
  Dval = NULL
  if(e-s <= delta){
    return(list(S = S, Dval = Dval, Level = Level))
  }else{
    level = level + 1
    a = rep(0, e-s-1)
    for(t in (s+1):(e-1)){
      a[t-s] = Delta_se_t(Y, s, e, t, N)
    }
    best_value = max(a)
    best_t = which.max(a) + s
    temp1 = NBS(Y, delta, s, best_t-1, N)
    temp2 = NBS(Y, delta, best_t, e, N)
    S = c(temp1$S, best_t, temp2$S)
    Dval = c(temp1$Dval, best_value, temp2$Dval)
    Level = c(temp1$Level, level, temp2$Level)
    return(list(S = S, Dval = Dval, Level = Level))
  }
}



#' @title Internal Function: Compute value of CUSUM statistic based on KS distance.
#' @param Y         A \code{numeric} matrix of observations with with horizontal axis being time, and vertical axis being multiple observations on each time point.
#' @param s         A \code{integer} scalar of starting index.
#' @param e         A \code{integer} scalar of ending index.
#' @param t         A \code{integer} scalar of splitting index.
#' @param N         A \code{integer} vector representing number of multiple observations on each time point.
#' @return  A \code{numeric} scalar of value of CUSUM statistic based on KS distance.
#' @noRd
Delta_se_t = function(Y, s, e, t, N){
  n_st = sum(N[s:t])  #n*(t-s+1)
  n_se = sum(N[s:e])  #n*(e-s+1)
  n_te = sum(N[(t+1):e]) #n*(e-(t+1) +1)
  
  aux = Y[, s:t]
  aux = aux[which(is.na(aux)==FALSE)]
  temp = ecdf(aux)
  vec_y = Y[, s:e]
  vec_y = vec_y[which(is.na(vec_y)==FALSE)]
  Fhat_st = temp(vec_y)# temp(grid)
  
  aux = Y[, (t+1):e]
  aux = aux[which(is.na(aux)==FALSE)]
  temp = ecdf(aux)
  Fhat_te = temp(vec_y)# temp(grid)
  
  temp = sqrt(n_st * n_te / n_se) * max(abs(Fhat_te - Fhat_st))
  return(temp)
}


new_BS_threshold = function(temp, tau, p){
  idx = which(temp$Dval > tau)
  S_hat = c()
  if(length(idx) == 0){
    return(NULL)
  }
  for(j in 1:length(idx)){
    if(p[idx[j]] == 0){
      S_hat = c(S_hat, temp$S[idx[j]])
    }
    if(p[ind[j]] > 0 && min(abs(S_hat - temp$S[p[ind[j]]])) == 0){
      S_hat = c(S_hat, temp$S[ind[j]])
    }
  }
  return(Shat)
}




NBS_full =  function(y, z, gam, N){
  n = max(N)
  T = dim(y)[1]
  temp1 = new_BS(z, gam,1,T,0,NULL,NULL,1, N)
  Dval = temp1$Dval
  p1 =  parent(temp1)
  aux = sort(Dval,decreasing = TRUE)
  tau_grid = rev(aux[1:min(20,length(Dval))]-10^{-5})
  tau_grid =  tau_grid[which(is.na(tau_grid)==FALSE)] ### *
  tau_grid = c(tau_grid,10)

  S =  c()
  for( j in 1:length(tau_grid))
  {
    aux = new_BS_threshold(temp1,tau_grid[j],p1)
    if(length(aux)==0)
    {
      break;
    }
    S[[j]] = sort(aux)
  }
  #}
  T= dim(y)[1]
  S = unique(S)
  if(length(S)==0)
  {
    return(NULL)
  }
  lamb =log(sum(N))/1.5#2.5#1.5#2#2.555#
  for(j in 1:length(S))#)
  {
    if(length(S[[j]])==0)
    {
      j = j+1;

      if(j>length(S))
        break;
    }

    B2  =  S[[j]]
    if(j==length(S))
    {
      B1 = NULL
    }
    if(j< length(S))
    {
      B1 = S[[j+1]]
    }
    temp = setdiff(B2,B1)

    st =  -10^15
    #Delta_se_t(z,eta1,eta2,eta,N)^2
    for(l in 1:length(temp))
    {
      eta =  temp[l]

      if(j == length(S))
      {
        eta1 = 1
        eta2 = T
      }
      if(j < length(S))
      {
        for(k in 1:length(S[[j+1]]))
        {
          if(S[[j+1]][k]> eta  )
            break;
        }
        if(S[[j+1]][k]> eta )
        {
          eta2 = S[[j+1]][k]

          if(k ==1)
            eta1 = 1

          if(k > 1)
            eta1 = S[[j+1]][k-1]+1
        }
        if(S[[j+1]][k]< eta )
        {
          eta1 = S[[j+1]][k]+1
          eta2 = T
        }
      }
      st_aux = Delta_se_t(y,eta1,eta2,eta,N)^2
      # print(st_au)
      if(st_aux> st)
      {
        st = st_aux
      }
    }###  close for defining  eta1 and eta2


    # print(c1 - c2 - Delta_se_t(z,eta1+1,eta2,eta,N)^2 + lamb)
    if(st >   lamb)
    {
      #B1 = B2
      return(B2)
    }
    # print(st)
  }
  #c1 - c2 - Delta_se_t(z,eta1+1,eta2,eta,N)^2 + lamb

  return(B1)
  #

}