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
CV.search.DP.univar = function(gamma.set, delta, y, ...){
  output = sapply(1:length(gamma.set), function(j) CV.DP.univar(gamma.set[j], delta, y))
  print(output)
  cpt_hat = output[1,]## estimated change points
  K_hat = output[2,]## number of estimated change points
  test_error = output[3,]## validation loss
  train_error = output[4,]## training loss                                                      
  result = list(cpt_hat = cpt_hat, K_hat = K_hat, test_error = test_error, train_error = train_error)
  return(result)
}


#' @title Standard binary segmentation for univariate mean change points detection
#' @description TO DO
#' @param y         A \code{numeric} vector of observations.
#' @param s         A \code{integer} scalar of starting index.
#' @param e         A \code{integer} scalar of ending index.
#' @param delta     A strictly \code{integer} scalar of minimum spacing.
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
#' y = c(rnorm(100, 0, 1), rnorm(100, 10, 10), rnorm(100, 40, 10))
#' temp = BS.univar(y, 1, 300, 5)
#' plot.ts(y)
#' points(x = tail(temp$S[order(temp$Dval)],20), y =rep(0, 20), col = "red")
BS.univar = function(y, s, e, delta, level = 0, ...){
  S = NULL
  Dval = NULL
  Level = NULL
  if(e-s <= delta){
    return(list(S = S, Dval = Dval, Level = Level))
  }else{
    level = level + 1
    a = rep(0, e-s-1)
    for(t in (s+1):(e-1)){
      a[t-s] = sqrt((t-s+1) * (e-t) / (e-s+1)) * abs(mean(y[s:t]) - mean(y[(t+1):e]))
    }
    best_value = max(a)
    best_t = which.max(a) + s
    temp1 = BS.univar(y, s, best_t-1, delta, level)
    temp2 = BS.univar(y, best_t, e, delta, level)
    S = c(temp1$S, best_t, temp2$S)
    Dval = c(temp1$Dval, best_value, temp2$Dval)
    Level = c(temp1$Level, level, temp2$Level)
    return(list(S = S, Dval = Dval, Level = Level))
  }
}
