#' @title Dynamic programming for univariate mean change points detection through \eqn{l_0} penalty.
#' @description     Perform dynamic programming for univariate mean change points detection.
#' @param y         A \code{numeric} vector of observations.
#' @param gamma     A \code{numeric} scalar of the tuning parameter associated with \eqn{l_0} penalty.
#' @param delta     A positive \code{integer} scalar of minimum spacing.
#' @return An object of \code{\link[base]{class}} "DP", which is a \code{list} with the following structure:
#'  \item{partition}{A vector of the best partition.}
#'  \item{yhat}{A vector of mean estimation for corresponding to the best partition.}
#'  \item{cpt}{A vector of change points estimation.}
#' @export
#' @author Haotian Xu
#' @references Wang, Yu and Rinaldo (2020) <doi:10.1214/20-EJS1710>
#' @examples
#' set.seed(123)
#' cpt_true = c(20, 50, 170)
#' y = rnorm(300) + c(rep(0,20),rep(1,30),rep(0,120),rep(1,130))
#' DP_result = DP.univar(y, gamma = 5, delta = 5)
#' cpt_hat = DP_result$cpt
#' Hausdorff.dist(cpt_hat, cpt_true)
DP.univar <- function(y, gamma, delta){
  DP_result = .Call('_changepoints_rcpp_DP_univar', PACKAGE = 'changepoints', y, gamma, delta)
  result = append(DP_result, list(cpt = part2local(DP_result$partition)))
  class(result) = "DP"
  return(result)
}


#' @title Internal function: cross-validation of dynamic programming algorithm for univariate mean change points detection.
#' @description     Perform cross-validation by sample splitting. Using the sample with odd indices as training data to estimate the change points, then computing sample mean for each segment within two consecutive change points, and computing the validation error based on the sample with even indices.
#' @param y         A \code{numeric} vector of observations.
#' @param gamma     A \code{numeric} scalar of the tuning parameter associated with the \eqn{l_0} penalty.
#' @param delta     A positive \code{integer} scalar of minimum spacing.
#' @return  A \code{list} with the following structure:
#'  \item{cpt_hat}{A vector of estimated change points locations (sorted in strictly increasing order).}
#'  \item{K_hat}{A scalar of number of estimated change points.}
#'  \item{test_error}{A vector of testing errors.}
#'  \item{train_error}{A vector of training errors.}
#' @noRd
CV.DP.univar = function(y, gamma, delta){
  N = length(y)
  even_indexes = seq(2, N, 2)
  odd_indexes = seq(1, N, 2)
  train_y = y[odd_indexes]
  validation_y = y[even_indexes]
  init_cpt_train = DP.univar(train_y, gamma, delta)$cpt
  init_cpt_train_ext = c(0, init_cpt_train, length(train_y))
  init_cpt = odd_indexes[init_cpt_train]
  len = length(init_cpt)
  init_cpt_ext = c(init_cpt_train, N/2)
  interval = matrix(0, nrow = len+1, ncol = 2)
  interval[1,] = c(1, init_cpt_ext[1])
  if(len > 0){
    for(j in 2:(1+len)){
      interval[j,] = c(init_cpt_ext[j-1]+1, init_cpt_ext[j])
    }
  }
  y_hat_train = sapply(1:(len+1), function(index) mean(train_y[(interval[index,1]):(interval[index,2])]))
  training_loss = sapply(1:(len+1), function(index) norm(train_y[(interval[index,1]):(interval[index,2])] - y_hat_train[index], type = "2")^2)
  validationmat = sapply(1:(len+1), function(index) norm(validation_y[(interval[index,1]):(interval[index,2])] - y_hat_train[index], type = "2")^2)
  result = list(cpt_hat = init_cpt, K_hat = len, test_error = sum(validationmat), train_error = sum(training_loss))
  return(result)
}


#' @title Grid search for dynamic programming to select the tuning parameter through Cross-Validation.
#' @description Perform grid search for dynamic programming to select the tuning parameter through Cross-Validation.
#' @param gamma_set     A \code{numeric} vector of candidate tuning parameter associated with the \eqn{l_0} penalty.
#' @param y             A \code{numeric} vector of observations.
#' @param delta         A positive \code{integer} scalar of minimum spacing.
#' @return  A \code{list} with the following structure:
#'  \item{cpt_hat}{A list of vector of estimated change points (sorted in strictly increasing order).}
#'  \item{K_hat}{A list of scalar of number of estimated change points.}
#'  \item{test_error}{A list of vector of testing errors.}
#'  \item{train_error}{A list of vector of training errors.}
#' @export
#' @author  Daren Wang & Haotian Xu
#' @references Wang, Yu and Rinaldo (2020) <doi:10.1214/20-EJS1710>
#' @examples
#' set.seed(0)
#' cpt_true = c(20, 50, 170)
#' y = rnorm(300) + c(rep(0,20),rep(2,30),rep(0,120),rep(2,130))
#' gamma_set = 1:5
#' DP_result = CV.search.DP.univar(y, gamma_set, delta = 5)
#' min_idx = which.min(DP_result$test_error)
#' cpt_hat = unlist(DP_result$cpt_hat[min_idx])
#' Hausdorff.dist(cpt_hat, cpt_true)
CV.search.DP.univar = function(y, gamma_set, delta){
  output = sapply(1:length(gamma_set), function(j) CV.DP.univar(y, gamma_set[j], delta))
  #print(output)
  cpt_hat = output[1,]## estimated change points
  K_hat = output[2,]## number of estimated change points
  test_error = output[3,]## validation loss
  train_error = output[4,]## training loss                                                      
  result = list(cpt_hat = cpt_hat, K_hat = K_hat, test_error = test_error, train_error = train_error)
  return(result)
}


#' @title Local refinement of an initial estimator for univariate mean change points detection.
#' @description     Perform local refinement for univariate mean change points detection.
#' @param cpt_init  An \code{integer} vector of initial change points estimation (sorted in strictly increasing order).
#' @param y         A \code{numeric} vector of univariate time series.
#' @return  An \code{integer} vector of locally refined change point estimation.
#' @export
#' @author Haotian Xu
#' @references Wang, Yu and Rinaldo (2020) <doi:10.1214/20-EJS1710>.
#' @examples
#' set.seed(0)
#' cpt_true = c(20, 50, 170)
#' y = rnorm(300) + c(rep(0,20),rep(2,30),rep(0,120),rep(2,130))
#' gamma_set = 1:5
#' DP_result = CV.search.DP.univar(y, gamma_set, delta = 5)
#' min_idx = which.min(DP_result$test_error)
#' cpt_hat = unlist(DP_result$cpt_hat[min_idx])
#' Hausdorff.dist(cpt_hat, cpt_true)
#' cpt_LR = local.refine.univar(cpt_hat, y)
#' Hausdorff.dist(cpt_LR, cpt_true)
local.refine.univar = function(cpt_init, y){
  w = 0.9
  n = length(y)
  cpt_init_ext = c(0, cpt_init, n)
  cpt_init_numb = length(cpt_init)
  cpt_refined = rep(0, cpt_init_numb+1)
  for (k in 1:cpt_init_numb){
    s = w*cpt_init_ext[k] + (1-w)*cpt_init_ext[k+1]
    e = (1-w)*cpt_init_ext[k+1] + w*cpt_init_ext[k+2]
    lower = ceiling(s) + 1
    upper = floor(e) - 1
    b = sapply(lower:upper, function(eta) sum((y[ceiling(s):eta] - mean(y[ceiling(s):eta]))^2) + sum((y[(eta+1):floor(e)] - mean(y[(eta+1):floor(e)]))^2))
    cpt_refined[k+1] = ceiling(s) + which.min(b)
  }
  return(cpt_refined[-1])
}


#' @title Standard binary segmentation for univariate mean change points detection.
#' @description     Perform standard binary segmentation for univariate mean change points detection.
#' @param y         A \code{numeric} vector of observations.
#' @param s         A \code{integer} scalar of starting index.
#' @param e         A \code{integer} scalar of ending index.
#' @param delta     A positive \code{numeric} scalar of minimum spacing.
#' @param level     Should be fixed as 0.
#' @return  An object of \code{\link[base]{class}} "BS", which is a \code{list} with the following structure:
#'  \item{S}{A vector of estimated change point locations (sorted in strictly increasing order).}
#'  \item{Dval}{A vector of values of CUSUM statistic.}
#'  \item{Level}{A vector representing the levels at which each change point is detected.}
#'  \item{Parent}{A matrix with the starting indices on the first row and the ending indices on the second row.}
#' @export
#' @author  Haotian Xu
#' @references Wang, Yu and Rinaldo (2020) <doi:10.1214/20-EJS1710>.
#' @seealso \code{\link{thresholdBS}} for obtaining change points estimation, \code{\link{tuneBSunivar}} for a tuning version.
#' @examples
#' set.seed(0)
#' cpt_true = c(20, 50, 170)
#' y = rnorm(300) + c(rep(0,20),rep(2,30),rep(0,120),rep(2,130))
#' temp = BS.univar(y, 1, length(y), delta = 5)
#' plot.ts(y)
#' points(x = tail(temp$S[order(temp$Dval)],4),
#'        y = y[tail(temp$S[order(temp$Dval)],4)], col = "red")
#' BS_result = thresholdBS(temp, tau = 4)
#' BS_result
#' print(BS_result$BS_tree, "value")
#' plot(BS_result$BS_tree)
#' print(BS_result$BS_tree_trimmed, "value")
#' plot(BS_result$BS_tree_trimmed)
#' cpt_hat = sort(BS_result$cpt_hat[,1]) # the threshold tau is specified to be 4
#' Hausdorff.dist(cpt_hat, cpt_true)
#' cpt_LR = local.refine.univar(cpt_hat, y)
#' Hausdorff.dist(cpt_LR, cpt_true)
BS.univar = function(y, s, e, delta = 2, level = 0){
  S = NULL
  Dval = NULL
  Level = NULL
  Parent = NULL
  if(e-s <= 2*delta){
    return(list(S = S, Dval = Dval, Level = Level, Parent = Parent))
  }else{
    level = level + 1
    parent = matrix(c(s, e), nrow = 2)
    a = rep(0, e-s-2*delta+1)
    for(t in (s+delta):(e-delta)){
      a[t-s-delta+1] = sqrt((t-s) * (e-t) / (e-s)) * abs(mean(y[(s+1):t]) - mean(y[(t+1):e]))
    }
    best_value = max(a)
    best_t = which.max(a) + s + delta - 1
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
#' @return  An object of \code{\link[base]{class}} "BS", which is a \code{list} with the following structure:
#'  \item{S}{A vector of estimated change point locations (sorted in strictly increasing order).}
#'  \item{Dval}{A vector of values of CUSUM statistic.}
#'  \item{Level}{A vector representing the levels at which each change point is detected.}
#'  \item{Parent}{A matrix with the starting indices on the first row and the ending indices on the second row.}
#' @export
#' @author Haotian Xu
#' @references Wang, Yu and Rinaldo (2020) <doi:10.1214/20-EJS1710>.
#' @seealso \code{\link{thresholdBS}} for obtaining change points estimation, \code{\link{tuneBSunivar}} for a tuning version.
#' @examples
#' set.seed(0)
#' cpt_true = c(20, 50, 170)
#' y = rnorm(300) + c(rep(0,20),rep(2,30),rep(0,120),rep(2,130))
#' intervals = WBS.intervals(M = 300, lower = 1, upper = length(y))
#' temp = WBS.univar(y, 1, length(y), intervals$Alpha, intervals$Beta, delta = 5)
#' plot.ts(y)
#' points(x = tail(temp$S[order(temp$Dval)],4),
#'        y = y[tail(temp$S[order(temp$Dval)],4)], col = "red")
#' WBS_result = thresholdBS(temp, tau = 4)
#' print(WBS_result$BS_tree, "value")
#' plot(WBS_result$BS_tree)
#' print(WBS_result$BS_tree_trimmed, "value")
#' plot(WBS_result$BS_tree_trimmed)
#' cpt_hat = sort(WBS_result$cpt_hat[,1]) # the threshold tau is specified to be 4
#' Hausdorff.dist(cpt_hat, cpt_true)
#' cpt_LR = local.refine.univar(cpt_hat, y)
#' Hausdorff.dist(cpt_LR, cpt_true)
WBS.univar = function(y, s, e, Alpha, Beta, delta = 2, level = 0){ 
  Alpha_new = pmax(Alpha, s)
  Beta_new = pmin(Beta, e)
  idx = which(Beta_new - Alpha_new > 2*delta)
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
      temp = rep(0, Beta_new[m] - Alpha_new[m] - 2*delta + 1)
      for(t in (Alpha_new[m]+delta):(Beta_new[m]-delta)){
        temp[t-(Alpha_new[m]+delta)+1] = sqrt((t-Alpha_new[m]) * (Beta_new[m]-t) / (Beta_new[m]-Alpha_new[m])) * abs(mean(y[(Alpha_new[m]+1):t]) - mean(y[(t+1):Beta_new[m]]))
      }
      best_value = max(temp)
      best_t = which.max(temp) + Alpha_new[m] + delta - 1
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


#' @title Univariate mean change points detection based on standard or wild binary segmentation with tuning parameter selected by sSIC.
#' @description Perform univariate mean change points detection based on standard or wild binary segmentation. The threshold parameter tau for WBS is automatically selected based on the sSIC score defined in Equation (4) in Fryzlewicz (2014).
#' @param BS_object A "BS" object produced by \code{BS.univar} or \code{WBS.univar}.
#' @param y         A \code{numeric} vector of observations.
#' @return  A \code{list} with the following structure:
#'  \item{cpt}{A vector of estimated change point locations (sorted in strictly increasing order).}
#'  \item{tau}{A scalar of selected threshold tau based on sSIC.}
#' @export
#' @author Daren Wang & Haotian Xu
#' @references Wang, Yu and Rinaldo (2020) <doi:10.1214/20-EJS1710>;
#'             Fryzlewicz (2014), Wild binary segmentation for multiple change-point detection, <DOI: 10.1214/14-AOS1245>.
#' @seealso \code{\link{BS.univar}} and \code{\link{WBS.univar}}.
#' @examples
#' set.seed(0)
#' cpt_true = c(20, 50, 170)
#' y = rnorm(300) + c(rep(0,20),rep(2,30),rep(0,120),rep(2,130))
#' ## change points detection by BS
#' temp1 = BS.univar(y, 1, length(y), delta = 5)
#' BS_result = tuneBSunivar(temp1, y)
#' cpt_BS = BS_result$cpt
#' Hausdorff.dist(cpt_BS, cpt_true)
#' cpt_BS_LR = local.refine.univar(cpt_BS, y)
#' Hausdorff.dist(cpt_BS_LR, cpt_true)
#' ## change points detection by WBS
#' intervals = WBS.intervals(M = 300, lower = 1, upper = length(y))
#' temp2 = WBS.univar(y, 1, length(y), intervals$Alpha, intervals$Beta, delta = 5)
#' WBS_result = tuneBSunivar(temp2, y)
#' cpt_WBS = WBS_result$cpt
#' Hausdorff.dist(cpt_WBS, cpt_true)
#' cpt_WBS_LR = local.refine.univar(cpt_WBS, y)
#' Hausdorff.dist(cpt_WBS_LR, cpt_true)
tuneBSunivar = function(BS_object, y){
  UseMethod("tuneBSunivar", BS_object)
}

#' @export
tuneBSunivar.BS = function(BS_object, y){
  Dval = BS_object$Dval
  aux = sort(Dval, decreasing = TRUE)
  tau_grid = rev(aux[1:100]-10^{-4})
  tau_grid = tau_grid[which(is.na(tau_grid)==FALSE)]
  tau_grid = c(tau_grid,10)
  S = c()
  for(j in 1:length(tau_grid)){
    aux = thresholdBS(BS_object, tau_grid[j])$cpt_hat[,1]
    if(length(aux) == 0)
      break;
    S[[j]] = sort(aux)
  }
  S = unique(S)
  score = rep(0, length(S))
  for(j in 1:length(S)){
    score[j] = sSIC.obj(y, S[[j]])
  }
  best_ind = which.min(score)
  return(list(cpt = S[[best_ind]], tau = tau_grid[best_ind]))
}



#' @title Strengthened Schwarz information criterion (sSIC)
#' @references Fryzlewicz (2014), Wild binary segmentation for multiple change-point detection, <DOI: 10.1214/14-AOS1245>
#' @noRd
sSIC.obj = function(y, S){
  K = length(S)
  obs_num = length(y)
  S_ext = c(0, S, obs_num)
  f_hat = NULL
  for(i in 1:(K+1)){
    f_hat = c(f_hat, rep(mean(y[(S_ext[i]+1):S_ext[i+1]]), S_ext[i+1]-S_ext[i]))
  }
  sigma_hat = sum((y - f_hat)^2)/obs_num
  ssic = obs_num/2*log(sigma_hat) + K*(log(obs_num))^(1.01)
  return(ssic)
}


