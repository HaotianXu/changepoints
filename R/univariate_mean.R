#' @title Dynamic programming for univariate mean change points detection by l0 penalty (RCPP)
#' @description TO DO
#' @param gamma     A \code{numeric} scalar of the tuning parameter associated with the l0 penalty.
#' @param delta     A strictly \code{integer} scalar of minimum spacing.
#' @param y         A \code{numeric} vector of observations.
#' @param ...      Additional arguments.
#' @return TO DO.
#' @export
#' @author Haotian Xu
#' @examples
#' y = rnorm(300)+ c(rep(0,130),rep(-1,20),rep(1,20),rep(0,130))
#' D_P_univar(1, 5, y)
D_P_univar <- function(gamma, delta, y, ...) {
  .Call('_changepoints_rcpp_D_P_univar', PACKAGE = 'changepoints', gamma, delta, y)
}

#' @title Dynamic programming for univariate mean change points detection by l0 penalty
#' @description TO DO
#' @param gamma     A \code{numeric} scalar of the tuning parameter associated with the l0 penalty.
#' @param delta     A strictly \code{integer} scalar of minimum spacing.
#' @param y         A \code{numeric} vector of observations.
#' @param ...      Additional arguments.
#' @return TO DO.
#' @export
#' @author Haotian Xu
#' @examples
#' y = c(rep(0, 100), rep(1, 100))
#' DP.univar(1, 4, y)
DP.univar = function(gamma, delta, y, ...){
  N = length(y)
  bestvalue = rep(0,N+1)
  partition = rep(0,N)
  yhat = rep(NA, N)
  bestvalue[1] = -gamma
  for(r in 1:N){
    bestvalue[r+1] = Inf
    for(l in 1:r){
      if(r - l > delta){
        b = bestvalue[l] + gamma + (norm(y[l:r] - mean(y[l:r]), type = "2"))^2 
      }else{
        b = Inf
      }
      if(b < bestvalue[r+1]){
        bestvalue[r+1] = b
        partition[r] = l-1
      }
    }
  }
  r = N
  l = partition[r]
  while(r > 0){
    yhat[(l+1):r] = rep(mean(y[(l+1):r]), r-l)
    r = l
    l = partition[r]
  }
  return(list(partition = partition, yhat = yhat))
}