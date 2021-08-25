#' @title Dynamic programming for univariate mean change points detection.
#' @description     Perform dynamic programming for univariate mean change points detection through l0 penalty
#' @param y         A \code{numeric} vector of observations.
#' @param gamma     A \code{numeric} scalar of the tuning parameter associated with the l0 penalty.
#' @param delta     A positive \code{integer} scalar of minimum spacing.
#' @param ...       Additional arguments.
#' @return A \code{list} with the structure:
#' \itemize{
#'  \item{partition}{A vector of the best partition.}
#'  \item{yhat}{A vector of mean estimation for corresponding to the best partition.}
#' }
#' @export
#' @author Haotian Xu
#' @examples
#' y = rnorm(300) + c(rep(0,130),rep(-1,20),rep(1,20),rep(0,130))
#' DP.univar(y, 1, 5)
DP.univar <- function(y, gamma, delta, ...) {
  .Call('_changepoints_rcpp_DP_univar', PACKAGE = 'changepoints', y, gamma, delta)
}

# DP.univar = function(y, gamma, delta, ...){
#   N = length(y)
#   bestvalue = rep(0,N+1)
#   partition = rep(0,N)
#   yhat = rep(NA, N)
#   bestvalue[1] = -gamma
#   for(r in 1:N){
#     bestvalue[r+1] = Inf
#     for(l in 1:r){
#       if(r - l > delta){
#         b = bestvalue[l] + gamma + (norm(y[l:r] - mean(y[l:r]), type = "2"))^2 
#       }else{
#         b = Inf
#       }
#       if(b < bestvalue[r+1]){
#         bestvalue[r+1] = b
#         partition[r] = l-1
#       }
#     }
#   }
#   r = N
#   l = partition[r]
#   while(r > 0){
#     yhat[(l+1):r] = rep(mean(y[(l+1):r]), r-l)
#     r = l
#     l = partition[r]
#   }
#   return(list(partition = partition, yhat = yhat))
# }


#' @title Internal function: Cross-Validation of Dynamic Programming algorithm for univariate mean change points detection.
#' @description     Perform cross-validation for Dynamic Programming algorithm for univariate mean change points detection through l0 penalty.
#' @param y         A \code{numeric} vector of observations.
#' @param gamma     A \code{numeric} scalar of the tuning parameter associated with the l0 penalty.
#' @param delta     A positive \code{integer} scalar of minimum spacing.
#' @param ...       Additional arguments.
#' @return  A \code{list} with the structure:
#' \itemize{
#'  \item{cpt_hat}: A vector of estimated change points locations (sorted in strictly increasing order).
#'  \item{K_hat}: A scalar of number of estimated change points.
#'  \item{test_error}: A vector of testing errors.
#'  \item{train_error}: A vector of training errors.
#' } 
#' @noRd
CV.DP.univar = function(y, gamma, delta, ...){
  N = length(y)
  even_indexes = seq(2, N, 2)
  odd_indexes = seq(1, N, 2)
  train.y = y[odd_indexes]
  validation.y = y[even_indexes]
  init_cpt_train = part2local(DP.univar(train.y, gamma, delta)$partition)
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


#' @title Grid search for cross-validation of dynamic programming for univariate mean change points detection.
#' @description Perform grid search for Cross-Validation of Dynamic Programming for univariate mean change points detection through l0 penalty
#' @param gamma.set     A \code{numeric} vector of candidate tuning parameter associated with the l0 penalty.
#' @param y             A \code{numeric} vector of observations.
#' @param delta         A positive \code{integer} scalar of minimum spacing.
#' @param ...           Additional arguments.
#' @return  A \code{list} with the structure:
#' \itemize{
#'  \item{cpt_hat}: A list of vector of estimated change points locations (sorted in strictly increasing order).
#'  \item{K_hat}: A list of scalar of number of estimated change points.
#'  \item{test_error}: A list of vector of testing errors.
#'  \item{train_error}: A list of vector of training errors.
#' } 
#' @export
#' @author  Daren Wang and Haotian Xu
#' @examples
#' y = rnorm(300) + c(rep(0,130),rep(-1,20),rep(1,20),rep(0,130))
#' CV.search.DP.univar(y, gamma.set = 3:6, delta = 5)
CV.search.DP.univar = function(y, gamma.set, delta, ...){
  output = sapply(1:length(gamma.set), function(j) CV.DP.univar(y, gamma.set[j], delta))
  print(output)
  cpt_hat = output[1,]## estimated change points
  K_hat = output[2,]## number of estimated change points
  test_error = output[3,]## validation loss
  train_error = output[4,]## training loss                                                      
  result = list(cpt_hat = cpt_hat, K_hat = K_hat, test_error = test_error, train_error = train_error)
  return(result)
}


#' @title Standard binary segmentation for univariate mean change points detection.
#' @description     Perform standard binary segmentation for univariate mean change points detection
#' @param y         A \code{numeric} vector of observations.
#' @param s         A \code{integer} scalar of starting index.
#' @param e         A \code{integer} scalar of ending index.
#' @param delta     A positive \code{numeric} scalar of minimum spacing.
#' @param level     Should be fixed as 0.
#' @param ...      Additional arguments.
#' @return  A \code{list} with the structure:
#' \itemize{
#'  \item{S}: A vector of estimated change point locations (sorted in strictly increasing order).
#'  \item{Dval}: A vector of values of CUSUM statistic based on KS distance.
#'  \item{Level}: A vector representing the levels at which each change point is detected.
#'  \item{Parent}: A matrix with the starting indices on the first row and the ending indices on the second row.
#' } 
#' @export
#' @author Haotian Xu
#' @examples
#' y = c(rnorm(100, 0, 1), rnorm(100, 10, 10), rnorm(100, 40, 10))
#' temp = BS.univar(y, 1, length(y), delta = 5)
#' plot.ts(y)
#' points(x = tail(temp$S[order(temp$Dval)],4),
#'        y = y[tail(temp$S[order(temp$Dval)],4)], col = "red")
#' threshold.BS(temp, 20)
BS.univar = function(y, s, e, delta = 2, level = 0, ...){
  S = NULL
  Dval = NULL
  Level = NULL
  Parent = NULL
  if(e-s <= delta){
    return(list(S = S, Dval = Dval, Level = Level, Parent = Parent))
  }else{
    level = level + 1
    parent = matrix(c(s, e), nrow = 2)
    a = rep(0, e-s-1)
    for(t in (s+1):(e-1)){
      a[t-s] = sqrt((t-s) * (e-t) / (e-s)) * abs(mean(y[(s+1):t]) - mean(y[(t+1):e]))
    }
    best_value = max(a)
    best_t = which.max(a) + s
    temp1 = BS.univar(y, s, best_t-1, delta, level)
    temp2 = BS.univar(y, best_t, e, delta, level)
    S = c(temp1$S, best_t, temp2$S)
    Dval = c(temp1$Dval, best_value, temp2$Dval)
    Level = c(temp1$Level, level, temp2$Level)
    Parent = cbind(temp1$Parent, parent, temp2$Parent)
    result = list(S = S, Dval = Dval, Level = Level, Parent = Parent)
    class(result) = "BS"
    return(result)
  }
}


#' @title Wild binary segmentation for univariate mean change points detection.
#' @description     Perform wild binary segmentation for univariate mean change points detection.
#' @param y         A \code{numeric} vector of observations.
#' @param s         A \code{integer} scalar of starting index.
#' @param e         A \code{integer} scalar of ending index.
#' @param Alpha     A \code{integer} vector of starting indices of random intervals.
#' @param Beta      A \code{integer} vector of ending indices of random intervals.
#' @param delta     A positive \code{integer} scalar of minimum spacing.
#' @param level     Should be fixed as 0.
#' @param ...      Additional arguments.
#' @return  A \code{list} with the structure:
#' \itemize{
#'  \item{S}: A vector of estimated change point locations (sorted in strictly increasing order).
#'  \item{Dval}: A vector of values of CUSUM statistic based on KS distance.
#'  \item{Level}: A vector representing the levels at which each change point is detected.
#'  \item{Parent}: A matrix with the starting indices on the first row and the ending indices on the second row.
#' }
#' @export
#' @author Haotian Xu
#' @examples
#' y = c(rnorm(100, 0, 1), rnorm(100, 0, 10), rnorm(100, 0, 40))
#' intervals = WBS.intervals(M = 120, lower = 1, upper = length(y))
#' temp = WBS.univar(y, 1, length(y), intervals$Alpha, intervals$Beta, delta = 5)
#' plot.ts(y)
#' points(x = tail(temp$S[order(temp$Dval)], 4),
#'        y = y[tail(temp$S[order(temp$Dval)],4)], col = "red")
#' threshold.BS(temp, 100)
WBS.univar = function(y, s, e, Alpha, Beta, delta = 2, level = 0){ 
  Alpha_new = pmax(Alpha, s)
  Beta_new = pmin(Beta, e)
  idx = which(Beta_new - Alpha_new > delta)
  Alpha_new = Alpha_new[idx]
  Beta_new = Beta_new[idx]
  M = length(Alpha_new)
  # xi = 1/8
  # Alpha_new2 = Alpha_new
  # Beta_new2  = Beta_new
  # Alpha_new = ceiling((1-xi)*Alpha_new2+  xi*Beta_new2)
  # Beta_new = ceiling((1-xi)*Beta_new2 +  xi*Alpha_new2)
  # idx = which(Beta_new - Alpha_new > delta)
  # Alpha_new = Alpha_new[idx]
  # Beta_new = Beta_new[idx]
  # M = length(Alpha_new)
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
      temp = rep(0, Beta_new[m] - Alpha_new[m] - 1)
      for(t in (Alpha_new[m]+1):(Beta_new[m]-1)){
        temp[t-(Alpha_new[m])] = sqrt((t-Alpha_new[m]) * (Beta_new[m]-t) / (Beta_new[m]-Alpha_new[m])) * abs(mean(y[(Alpha_new[m]+1):t]) - mean(y[(t+1):Beta_new[m]]))
      }
      best_value = max(temp)
      best_t = which.max(temp) + Alpha_new[m]
      a[m] = best_value
      b[m] = best_t
    }
    m_star = which.max(a)
  }
  temp1 = WBS.univar(y, s, b[m_star]-1, Alpha, Beta, delta, level)
  temp2 = WBS.univar(y, b[m_star], e, Alpha, Beta, delta, level)
  S = c(temp1$S, b[m_star], temp2$S)
  Dval = c(temp1$Dval, a[m_star], temp2$Dval)
  Level = c(temp1$Level, level, temp2$Level)
  Parent = cbind(temp1$Parent, parent, temp2$Parent)
  result = list(S = S, Dval = Dval, Level = Level, Parent = Parent)
  class(result) = "BS"
  return(result)
}



