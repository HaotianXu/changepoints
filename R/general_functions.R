#' @title Partition to localization.
#' @description     Find change point locations from the best partition produced by dynamic programming.
#' @param parti_vec    A \code{numeric} vector of observations.
#' @param ...          Additional arguments.
#' @return  A vector of change point locations.
#' @export
#' @author Haotian Xu
#' @examples
#' y = c(rep(0, 100), rep(1, 100))
#' part2local(DP.univar(y, 1, 4)$partition)
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
#' @param ...       Additional arguments.
#' @return  An integer scalar.
#' @export
#' @author
#' @examples
#' vec1 = sample.int(1000, size = 50)
#' vec2 = sample.int(2000, size = 100)
#' Hausdorff.dist(vec1, vec2)
Hausdorff.dist = function(vec1, vec2, ...){
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


#' @title Thresholding a BS object.
#' @description         Given a BS object, perform thresholding to find the change point locations.
#' @param BS_object     A \code{numeric} vector of observations.
#' @param tau           A positive \code{numeric} scalar of thresholding value.
#' @param ...           Additional arguments.
#' @return  A \code{list} with the structure:
#' \itemize{
#'  \item{BS_tree}: {A list of data.frame containing "current" indices, "parent" indices, "location" indices and "value" of CUSUM at all steps.}
#'  \item{BS_tree_trimmed}: {BS_tree with change points which do not satisfy the thresholding criteria removed.}
#'  \item{change_points}: {A matrix contains change point locations, values of corresponding statistic, and levels at which each change point is detected.}
#' } 
#' @export
#' @author Haotian Xu
#' @examples
#' y = c(rnorm(100, 0, 1), rnorm(100, 10, 10), rnorm(100, 40, 10))
#' temp = BS.univar(y, 1, 300, 5)
#' plot.ts(y)
#' points(x = tail(temp$S[order(temp$Dval)],4), y = y[tail(temp$S[order(temp$Dval)],4)], col = "red")
#' threshold.BS(temp, 20)
threshold.BS = function(BS_object, tau, ...){
  level_unique = unique(BS_object$Level[order(BS_object$Level)])
  level_length = length(level_unique)
  BS_tree = vector("list", level_length)
  BS_tree[[1]] = data.frame(current = 1, parent = 1, location = BS_object$S[order(BS_object$Level)][1], value = BS_object$Dval[order(BS_object$Level)][1])
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
  return(list(BS_tree = BS_tree, BS_tree_trimmed = BS_tree_trimmed, change_points = change_points))
}



#' @title Internal Function for BS thresholding: Given a set of indices in the current step, search for all indices of children in the next step.
#' @param idx_remove_parent     A \code{integer} vector of indices in the current step.
#' @param data_children         A \code{data.frame} including "current" indices, "parent" indices, "location" indices and "value" of CUSUM at the next step.
#' @return  A \code{list} with the structure:
#' \itemize{
#'  \item{idx_remove_children}: {A vector of indices of children in the next step.}
#'  \item{data_children_trimmed}: {A data.frame being data_children with the observations (idx_remove_children) removed.}
#' } 
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
#' @param ...      Additional arguments.
#' @return         A \code{list} with the structure:
#' \itemize{
#'  \item{Alpha}: {A M-dim vector representing the starting indices of random intervals.}
#'  \item{Beta}: {A M-dim vector representing the ending indices of random intervals.}
#' } 
#' @export
#' @author    Oscar Hernan Madrid Padilla
#' @examples
#' WBS.intervals(120, lower = 1, upper = 300)
WBS.intervals = function(M, lower = 1, upper, ...){
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




#' @title Generate coordinates of lower triangular matrix of dimension p.
#' @param p   A \code{integer} scalar of dimensionality.
#' @return    A \code{integer} (p*(p-1)/2)-dim vector representing the indices of the lower triangular entries (indices correspond to the vectorization by stacking columns).
#' @export
#' @author    Oscar Hernan Madrid Padilla
#' @examples
#' gen.lower.coordinate(5)
gen.lower.coordinate=function(p){
  mat=matrix(1:p^2, nrow = p)
  return(mat[lower.tri(mat, diag = F)])  
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