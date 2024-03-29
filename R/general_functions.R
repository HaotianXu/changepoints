#' @title Partition to localization for dynamic programming.
#' @description     Find change point locations from the best partition produced by dynamic programming altorithm.
#' @param parti_vec    A \code{integer} vector of best partition produced by dynamic programming algorithm.
#' @return  A vector of change point locations.
#' @noRd
part2local = function(parti_vec){
  N = length(parti_vec)
  localization = c()
  r = N
  l = parti_vec[r]
  localization = c(l, localization)
  while(r > 0){
    r = l
    l = parti_vec[r]
    localization = c(l, localization)
  }
  return(localization[-1])
}


#' @title Bidirectional Hausdorff distance.
#' @description     Compute the bidirectional Hausdorff distance between two sets.
#' @param vec1      A \code{integer} vector forms a subset of {1, 2, ..., n}.
#' @param vec2      A \code{integer} vector forms a subset of {1, 2, ..., n}.
#' @return  An integer scalar of bidirectional Hausdorff distance.
#' @export
#' @author  Daren Wang
#' @examples
#' vec1 = sample.int(1000, size = 50)
#' vec2 = sample.int(2000, size = 100)
#' Hausdorff.dist(vec1, vec2)
Hausdorff.dist = function(vec1, vec2){
  vec = c(vec1, vec2)
  if(!all(c(vec == floor(vec)), vec >= 0)){
    stop("vec1 and vec2 should be subsets of {0, 1, ...}")
  }
  dist = matrix(0, nrow = length(vec1), ncol = length(vec2))
  for (i in 1:nrow(dist)){
    for (j in 1:ncol(dist)){
      dist[i,j] = abs(vec1[i] - vec2[j])
    }
  }
  dH = max(max(apply(dist, 2, function(x) min(x))), max(apply(dist, 1, function(x) min(x))))
  return(dH)
}


#' @title Thresholding a BS object with threshold value tau.
#' @description         Given a BS object, perform thresholding to find the change point locations.
#' @param BS_object     A \code{BS} object.
#' @param tau           A positive \code{numeric} scalar of thresholding value.
#' @return  A \code{list} with the following structure:
#'  \item{BS_tree_trimmed}{BS_tree with change points which do not satisfy the thresholding criteria removed}
#'  \item{cpt_hat}{A matrix contains change point locations, values of corresponding statistic, and levels at which each change point is detected}
#' @export
#' @author Haotian Xu
#' @seealso \code{\link{BS.univar}}, \code{\link{BS.uni.nonpar}}, \code{\link{BS.cov}}, \code{\link{WBS.univar}}, \code{\link{WBS.uni.nonpar}}, \code{\link{WBS.multi.nonpar}}, \code{\link{WBS.network}}, \code{\link{WBSIP.cov}}
#' @examples
#' y = c(rnorm(100, 0, 1), rnorm(100, 10, 10), rnorm(100, 40, 10))
#' temp = BS.univar(y, 1, 300, 5)
#' plot.ts(y)
#' points(x = tail(temp$S[order(temp$Dval)],4), y = y[tail(temp$S[order(temp$Dval)],4)], col = "red")
#' thresholdBS(temp, 20)
thresholdBS <- function(BS_object, tau){
  UseMethod("thresholdBS", BS_object)
}


#' @export
thresholdBS.BS = function(BS_object, tau){
  if(tau <= 0){
    stop("The threshold tau should be a positive value.")
  }
  level_unique = unique(BS_object$Level[order(BS_object$Level)])
  level_length = length(level_unique)
  BS_tree = vector("list", level_length)
  BS_tree[[1]] = data.frame(current = 1, parent = NA, location = BS_object$S[order(BS_object$Level)][1], value = BS_object$Dval[order(BS_object$Level)][1])
  for(i in 2:level_length){
    idx_curr = cumsum(table(BS_object$Level))[i-1] + 1:table(BS_object$Level)[i]
    idx_prev = cumsum(table(BS_object$Level))[i-1] + 1 - table(BS_object$Level)[i-1]:1
    interval_prev = as.matrix(BS_object$Parent[,order(BS_object$Level)][,idx_prev])
    e_curr = BS_object$Parent[,order(BS_object$Level)][2,idx_curr]
    BS_tree[[i]] = data.frame(current = 1:length(idx_curr),
                              parent = sapply(e_curr, function(x) which(rbind(interval_prev[1,] <= x & interval_prev[2,] >= x))), 
                              location = BS_object$S[order(BS_object$Level)][idx_curr],
                              value = BS_object$Dval[order(BS_object$Level)][idx_curr])
  }
  BS_tree_new = BS_tree
  BS_tree_new[[1]] = BS_tree[[1]][,3:4]
  BS_tree_new[[1]]$location = paste0("N",BS_tree_new[[1]]$location)
  for(j in 2:level_length){
    BS_tree_new[[j]]$parent = sapply(BS_tree_new[[j]]$parent, function(x){BS_tree_new[[j-1]]$location[x]})
    BS_tree_new[[j]]$location = paste0(BS_tree_new[[j]]$parent, "$N", BS_tree_new[[j]]$location)
    BS_tree_new[[j]] = BS_tree_new[[j]][,3:4]
  }
  binary_tree = list()
  binary_tree$name = "Binary Segmentation Tree"
  for(j in 1:level_length){
    for(k in 1:nrow(BS_tree_new[[j]])){
      eval(parse(text=paste0("binary_tree$",BS_tree_new[[j]]$location[k],"$value","<-",BS_tree_new[[j]]$value[k])))
    }
  }
  BS_tree_node = data.tree::as.Node(binary_tree)
  
  BS_tree_trimmed = BS_tree
  for(i in 1:level_length){
    idx_remove = BS_tree_trimmed[[i]]$current[BS_tree_trimmed[[i]]$value <= tau]
    BS_tree_trimmed[[i]] = BS_tree_trimmed[[i]][BS_tree_trimmed[[i]]$value > tau,]
    if(length(idx_remove) > 0){
      idx_remove_parent = idx_remove
      k = i+1
      while(length(idx_remove_parent) > 0 & k <= level_length){
        temp = one.step.trim(idx_remove_parent, BS_tree_trimmed[[k]])
        BS_tree_trimmed[[k]] = temp$data_children_trimmed
        idx_remove_parent = temp$idx_remove_children
        k = k + 1
      }
    }
  }
  points_at_level = sapply(BS_tree_trimmed, dim)[1,]
  level = unlist(sapply(1:length(points_at_level), function(x) rep(x, points_at_level[x])))
  change_points = do.call(rbind, BS_tree_trimmed)[,c(3,4)]
  change_points$level = level
  rownames(change_points) = c()
  BS_tree_trimmed = BS_tree_trimmed[points_at_level != 0]
  if(length(BS_tree_trimmed) == 0){
    return(list(BS_tree = BS_tree_node, BS_tree_trimmed = NULL, cpt_hat = NULL))
  }
  BS_tree_trimmed_new = BS_tree_trimmed
  BS_tree_trimmed_new[[1]] = BS_tree_trimmed[[1]][,3:4]
  BS_tree_trimmed_new[[1]]$location = paste0("N",BS_tree_trimmed_new[[1]]$location)
  if(length(BS_tree_trimmed) == 1){
    return(list(BS_tree_trimmed = BS_tree_trimmed[[1]], cpt_hat = change_points))
  }
  for(j in 2:sum(points_at_level != 0)){
    BS_tree_trimmed_new[[j]]$parent = sapply(BS_tree_trimmed_new[[j]]$parent, function(x){BS_tree_trimmed_new[[j-1]]$location[rownames(BS_tree_trimmed_new[[j-1]]) == x]})
    BS_tree_trimmed_new[[j]]$location = paste0(BS_tree_trimmed_new[[j]]$parent, "$N", BS_tree_trimmed_new[[j]]$location)
    BS_tree_trimmed_new[[j]] = BS_tree_trimmed_new[[j]][,3:4]
  }
  binary_tree_trimmed = list()
  binary_tree_trimmed$name = paste0("Binary Segmentation Tree Trimmed with tau = ", tau)
  for(j in 1:sum(points_at_level != 0)){
    for(k in 1:nrow(BS_tree_trimmed_new[[j]])){
      eval(parse(text=paste0("binary_tree_trimmed$",BS_tree_trimmed_new[[j]]$location[k],"$value","<-",BS_tree_trimmed_new[[j]]$value[k])))
    }
  }
  BS_tree_trimmed_node = data.tree::as.Node(binary_tree_trimmed)
  return(list(BS_tree_trimmed = BS_tree_trimmed_node, cpt_hat = change_points))
}



#' @title Internal Function for BS thresholding: Given a set of indices in the current step, search for all indices of children in the next step.
#' @param idx_remove_parent     A \code{integer} vector of indices in the current step.
#' @param data_children         A \code{data.frame} including "current" indices, "parent" indices, "location" indices and "value" of CUSUM at the next step.
#' @return  A \code{list} with the following structure:
#'  \item{idx_remove_children}{A vector of indices of children in the next step}
#'  \item{data_children_trimmed}{A data.frame being data_children with the observations (idx_remove_children) removed}
#' @noRd
one.step.trim = function(idx_remove_parent, data_children){
  idx_remove_children = NULL
  for(j in idx_remove_parent){
    idx_remove_children = c(idx_remove_children, data_children$current[data_children$parent == j])
  }
  if(length(idx_remove_children) == 0){
    data_children_trimmed = data_children
  }else{
    data_children_trimmed = data_children[!(data_children$current %in% idx_remove_children),]
  }
  return(list(idx_remove_children = idx_remove_children, data_children_trimmed = data_children_trimmed))
}



#' @title Generate random intervals for WBS.
#' @description    Generate random intervals for WBS.
#' @param M        A positive \code{integer} scalar of number of random intervals.
#' @param lower    A positive \code{integer} scalar of lower bound of random intervals.
#' @param upper    A positive \code{integer} scalar of upper bound of random intervals.
#' @return         A \code{list} with the following structure:
#'  \item{Alpha}{A M-dim vector representing the starting indices of random intervals}
#'  \item{Beta}{A M-dim vector representing the ending indices of random intervals}
#' @export
#' @author    Oscar Hernan Madrid Padilla
#' @examples
#' WBS.intervals(120, lower = 1, upper = 300)
WBS.intervals = function(M, lower = 1, upper){
  if(lower >= upper){
    stop("Integer 'lower' should be strictly smaller than integer 'upper'.")
  }
  Alpha = sample(x = lower:upper, size = M, replace = TRUE)
  Beta = sample(x = lower:upper, size = M, replace = TRUE)
  for(j in 1:M){
    aux = Alpha[j]
    aux2 = Beta[j]
    Alpha[j] = min(aux, aux2)
    Beta[j] = max(aux, aux2)
  }
  return(list(Alpha = Alpha, Beta = Beta))
}



#' @title Transform a vector containing lower diagonal entries into a symmetric matrix of dimension p.
#' @param lowertri_vec  A \code{numeric} vector containing lower diagonal entries.
#' @param p             A \code{integer} scalar of dimensionality.
#' @param diag          A \code{logic} scalar indicating if the diagonal entries are contained in lowertri_vec.
#' @return   A \code{numeric} p x p symmetric matrix.
#' @export
#' @author    Haotian Xu
#' @examples
#' A = matrix(1:16, 4, 4)
#' B = lowertri2mat(A[lower.tri(A)], 4, diag = FALSE)
#' C = lowertri2mat(A[lower.tri(A, diag = TRUE)], 4, diag = TRUE)
lowertri2mat = function(lowertri_vec, p, diag = FALSE){
  aux_mat = matrix(0, nrow = p, ncol = p)
  if(diag == FALSE){
    aux_mat[lower.tri(aux_mat)] = lowertri_vec
    aux_mat[upper.tri(aux_mat)] = t(aux_mat)[upper.tri(aux_mat)]
  }else{
    aux_mat[lower.tri(aux_mat, diag = TRUE)] = lowertri_vec
    aux_mat[upper.tri(aux_mat)] = t(aux_mat)[upper.tri(aux_mat)]
  }
  return(aux_mat)
}




#' @title Generate population covariance matrix with dimension p.
#' @param p       A \code{integer} scalar of dimensionality.
#' @param sigma2  A positive \code{numeric} scalar representing the variance of each entry.
#' @param type    Specify the type of a covariance matrix: Diagonal structure ("diagonal"); Equal correlation structure ("equal"); Power decay structure ("power").
#' @return    A \code{numeric} p-by-p matrix.
#' @export
#' @author  Haotian Xu
#' @examples
#' gen.cov.mat(p = 5, sigma2 = 1, type = "diagonal")
gen.cov.mat = function(p, sigma2, type){
  if(type == "diagonal"){
    Sigma = diag(1, p)
  }else if(type == "equal"){
    Sigma = matrix(0.5, nrow = p, ncol = p) + diag(0.5, p)
  }else if(type == "power"){
    Sigma = matrix(NA, nrow = p, ncol = p)
    for(i in 1:p){
      for(j in 1:p){
        Sigma[i,j] = 0.5^(abs(i-j))
      }
    }
  }
  return(Sigma)
}


#' @title Function to generate a matrix with values 0 or 1, where 0 indicating the entry is missing
#' @param pi_mat   A \code{numeric} pxp matrix, for each entry, the value representing the probability of missing.
#' @param symm     A \code{logic} scalar indicating if the output matrix needs to be symmetric.
#' @return   A \code{numeric} p x p matrix.
#' @export
#' @author    Haotian Xu
#' @examples
#' p = 5
#' pi_mat = matrix(0.9, p, p)
#' eta_mat = gen.missing(pi_mat, symm = TRUE)
gen.missing = function(pi_mat, symm = TRUE){
  if(ncol(pi_mat) != nrow(pi_mat)){
    stop("pi_mat should be a square matrix.")
  }
  if((symm == TRUE) & (isSymmetric(pi_mat) == FALSE)){
    stop("If symm is TRUE, pi_mat should be a symmetric matrix.")
  }
  p = ncol(pi_mat)
  temp_mat = matrix(rbinom(matrix(1,p,p), matrix(1,p,p), pi_mat), p, p)
  if(symm == TRUE){
    temp_mat[upper.tri(temp_mat)] = t(temp_mat)[upper.tri(temp_mat)]
  }
  return(temp_mat)
}



#' @noRd
lasso_standardized_seq <- function(Xtilde, Ytilde, lambda_seq, eps = 0.0001, ...){
  .Call('_changepoints_rcpp_lasso_standardized_seq', PACKAGE = 'changepoints', Xtilde, Ytilde, lambda_seq, eps)
}


#' @noRd
lasso_standardized_obj <- function(Xtilde, Ytilde, beta, lambda){
  .Call('_changepoints_rcpp_lasso_standardized_obj', PACKAGE = 'changepoints', Xtilde, Ytilde, beta, lambda)
}
  

#' @noRd
soft_threshold_scalar <- function(x, lambda){
  .Call('_changepoints_rcpp_soft_threshold_scalar', PACKAGE = 'changepoints', x, lambda)
}


#' @noRd
lasso_seq <- function(X, Y, lambda_seq, eps = 0.0001, ...){
  .Call('_changepoints_rcpp_lasso_seq', PACKAGE = 'changepoints', X, Y, lambda_seq, eps)
}


#' @noRd
lasso_standardized <- function(Xtilde, Ytilde, beta_start, lambda, eps = 0.0001){
  .Call('_changepoints_rcpp_lasso_standardized', PACKAGE = 'changepoints', Xtilde, Ytilde, beta_start, lambda, eps)
}

#' @noRd
standardizeXY <- function(X, Y){
  .Call('_changepoints_rcpp_standardizeXY', PACKAGE = 'changepoints', X, Y)
}

#' @noRd
rcpp_error_pred_seg_VAR1 <- function(X_futu, X_curr, s, e, lambda, delta, eps = 0.0001){
  .Call('_changepoints_rcpp_error_pred_seg_VAR1', PACKAGE = 'changepoints', X_futu, X_curr, s, e, lambda, delta, eps)
}



#' @export
print.BS = function(x, ...){
  level_unique = unique(x$Level[order(x$Level)])
  level_length = length(level_unique)
  BS_tree = vector("list", level_length)
  BS_tree[[1]] = data.frame(current = 1, parent = NA, location = x$S[order(x$Level)][1], value = x$Dval[order(x$Level)][1])
  for(i in 2:level_length){
    idx_curr = cumsum(table(x$Level))[i-1] + 1:table(x$Level)[i]
    idx_prev = cumsum(table(x$Level))[i-1] + 1 - table(x$Level)[i-1]:1
    interval_prev = as.matrix(x$Parent[,order(x$Level)][,idx_prev])
    e_curr = x$Parent[,order(x$Level)][2,idx_curr]
    BS_tree[[i]] = data.frame(current = 1:length(idx_curr),
                              parent = sapply(e_curr, function(x) which(rbind(interval_prev[1,] <= x & interval_prev[2,] >= x))), 
                              location = x$S[order(x$Level)][idx_curr],
                              value = x$Dval[order(x$Level)][idx_curr])
  }
  BS_tree_new = BS_tree
  BS_tree_new[[1]] = BS_tree[[1]][,3:4]
  BS_tree_new[[1]]$location = paste0("N",BS_tree_new[[1]]$location)
  for(j in 2:level_length){
    BS_tree_new[[j]]$parent = sapply(BS_tree_new[[j]]$parent, function(x){BS_tree_new[[j-1]]$location[x]})
    BS_tree_new[[j]]$location = paste0(BS_tree_new[[j]]$parent, "$N", BS_tree_new[[j]]$location)
    BS_tree_new[[j]] = BS_tree_new[[j]][,3:4]
  }
  binary_tree = list()
  binary_tree$name = "Binary Segmentation Tree"
  for(j in 1:level_length){
    for(k in 1:nrow(BS_tree_new[[j]])){
      eval(parse(text=paste0("binary_tree$",BS_tree_new[[j]]$location[k],"$value","<-",BS_tree_new[[j]]$value[k])))
    }
  }
  BS_tree_node = data.tree::as.Node(binary_tree)
  print(BS_tree_node, "value")
}


#' @export
print.DP = function(x, ...){
  cat("Estimated change points by DP:", x$cpt, "\n")
}


#' @title Element-wise adaptive Huber mean estimator.
#' @description Computes the element-wise adaptive Huber mean estimator.
#' @param x         A \code{numeric} vector of observations.
#' @param tau       A \code{numeric} scalar corresponding to the robustification parameter (larger than 0).
#' @return          A \code{numeric} scalar corresponding to the adaptive Huber mean estimator.
#' @export
#' @author Haotian Xu
#' @examples
#' set.seed(123)
#' y = rnorm(100)
#' mean(y)
#' huber_mean(y, 1.345)
huber_mean <- function(x, tau){
  if (tau <= 0){
    stop("tau should be larger than 0.")
  }
  .Call('_changepoints_rcpp_huber_mean', PACKAGE = 'changepoints', x, tau)
}


