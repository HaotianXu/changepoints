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
#' y = rnorm(300) + c(rep(0,130),rep(-1,20),rep(1,20),rep(0,130))
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


#' @title Internal Function: Cross-Validation of Dynamic Programming algorithm for univariate mean change points detection by l0 penalty
#' @description TO DO
#' @param gamma     A \code{numeric} scalar of the tuning parameter associated with the l0 penalty.
#' @param delta     A strictly \code{integer} scalar of minimum spacing.
#' @param y         A \code{numeric} vector of observations.
#' @param ...      Additional arguments.
#' @return TO DO.
#' @export
#' @author
#' @examples
#' TO DO
CV.DP.univar = function(gamma, delta, y, ...){
  N = length(y)
  even_indexes = seq(2, N, 2)
  odd_indexes = seq(1, N, 2)
  train.y = y[odd_indexes]
  validation.y = y[even_indexes]
  init_cpt_train = part2local(D_P_univar(gamma, delta, train.y)$partition)
  init_cpt_train.long = c(0, init_cpt_train, length(train.y))
  diff.point = diff(init_cpt_train.long)
  if (length(which(diff.point == 1)) > 0){
    print(paste("gamma =", gamma, ".", "Warning: Consecutive points detected. Try a larger gamma."))
    init_cpt = odd_indexes[init_cpt_train]
    len = length(init_cpt)
    result = list(cpt_hat = init_cpt, K_hat = len, test_error = Inf, train_error = Inf)
  }
  else{
    init_cpt = odd_indexes[init_cpt_train]
    len = length(init_cpt)
    init_cpt_long = c(init_cpt_train, N/2)
    interval = matrix(0, nrow = len+1, ncol = 2)
    interval[1,] = c(1, init_cpt_long[1])
    if(len > 0){
      for(j in 2:(1+len)){
        interval[j,] = c(init_cpt_long[j-1]+1, init_cpt_long[j])
      }
    }
    y_hat_train = sapply(1:(len+1), function(index) mean(train.y[(interval[index,1]):(interval[index,2])]))
    training_loss = sapply(1:(len+1), function(index) norm(train.y[(interval[index,1]):(interval[index,2])] - y_hat_train[index], type = "2")^2)
    validationmat = sapply(1:(len+1), function(index) norm(validation.y[(interval[index,1]):(interval[index,2])] - y_hat_train[index], type = "2")^2)
    result = list(cpt_hat = init_cpt, K_hat = len, test_error = sum(validationmat), train_error = sum(training_loss))
  }
  return(result)
}


#' @title Perform grid search based on Cross-Validation of Dynamic Programming algorithm for univariate mean change points detection by l0 penalty
#' @description TO DO
#' @param gamma.set     A \code{numeric} vector of candidate tuning parameter associated with the l0 penalty.
#' @param y             A \code{numeric} vector of observations.
#' @param delta         A strictly \code{integer} scalar of minimum spacing.
#' @param ...           Additional arguments.
#' @return TO DO.
#' @export
#' @author
#' @examples
#' TO DO
CV.search.DP.regression = function(gamma.set, delta, y, ...){
  output = sapply(1:length(gamma.set), function(j) CV.DP.univar(gamma.set[j], delta, y))
  print(output)
  cpt_hat = output[1,]## estimated change points
  K_hat = output[2,]## number of estimated change points
  test_error = output[3,]## validation loss
  train_error = output[4,]## training loss                                                      
  result = list(cpt_hat = cpt_hat, K_hat = K_hat, test_error = test_error, train_error = train_error)
  return(result)
}
