#' @title Partition to localization
#' @description TO DO
#' @param x         A \code{numeric} vector of observations.
#' @param ...      Additional arguments.
#' @return TO DO.
#' @export
#' @author Haotian Xu
#' @examples
#' y = c(rep(0, 100), rep(1, 100))
#' part2local(D_P_univar(1, 4, y)$partition)
part2local = function(x){
  N = length(x)
  localization = c()
  r = N
  l = x[r]
  localization = c(l, localization)
  while(r > 0){
    r = l
    l = x[r]
    localization = c(l, localization)
  }
  return(localization[-1])
}


#' @title Bidirectional Hausdorff distance
#' @description TO DO
#' @param vec1      A \code{integer} vector forms a subset of {1, 2, ..., n}.
#' @param vec2      A \code{integer} vector forms a subset of {1, 2, ..., n}.
#' @param ...      Additional arguments.
#' @return  An integer scalar.
#' @export
#' @author
#' @examples
#' vec1 = sample.int(1000, size = 50)
#' vec2 = sample.int(2000, size = 100)
#' Hausdorff.dist(vec1, vec2)
Hausdorff.dist = function(vec1, vec2, ...){
  vec = c(vec1, vec2)
  if(!all(c(vec == floor(vec)), vec >= 1)){
    stop("vec1 and vec2 should be subsets of {1, 2, ..., n}")
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


#' @title Thresholding the standard binary segmentation for univariate mean change points detection
#' @description TO DO
#' @param BS_result     A \code{numeric} vector of observations.
#' @param tau           A \code{numeric} scalar of thresholding value.
#' @param ...      Additional arguments.
#' @return  A \code{list} with the structure:
#' \itemize{
#'  \item BS_tree             A list of data.frame containing "current" indices, "parent" indices, "location" indices and "value" of CUSUM at all steps.
#'  \item BS_tree_trimmed     BS_tree with change points which do not satisfy the thresholding criteria removed.
#'  \item ...         Additional parameters.
#' } 
#' @export
#' @author Haotian Xu
#' @examples
#' y = c(rnorm(100, 0, 1), rnorm(100, 10, 10), rnorm(100, 40, 10))
#' temp = BS.univar(y, 1, 300, 5)
#' plot.ts(y)
#' points(x = tail(temp$S[order(temp$Dval)],4), y = y[tail(temp$S[order(temp$Dval)],4)], col = "red")
#' BS.threshold(temp, 10)
BS.threshold = function(BS_result, tau, ...){
  level_unique = unique(BS_result$Level[order(BS_result$Level)])
  level_length = length(level_unique)
  BS_tree = vector("list", level_length)
  BS_tree[[1]] = data.frame(current = 1, parent = 1, location = BS_result$S[order(BS_result$Level)][1], value = BS_result$Dval[order(BS_result$Level)][1])
  for(i in 2:level_length){
    idx_curr = cumsum(table(BS_result$Level))[i-1] + 1:table(BS_result$Level)[i]
    idx_prev = cumsum(table(BS_result$Level))[i-1] + 1 - table(BS_result$Level)[i-1]:1
    interval_prev = as.matrix(BS_result$Parent[,order(BS_result$Level)][,idx_prev])
    e_curr = BS_result$Parent[,order(BS_result$Level)][2,idx_curr]
    BS_tree[[i]] = data.frame(current = 1:length(idx_curr),
                              parent = sapply(e_curr, function(x) which(rbind(interval_prev[1,] <= x & interval_prev[2,] >= x))), 
                              location = BS_result$S[order(BS_result$Level)][idx_curr],
                              value = BS_result$Dval[order(BS_result$Level)][idx_curr])
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
#'  \item idx_remove_children       A vector of indices of children in the next step.
#'  \item data_children_trimmed     A data.frame being data_children with the observations (idx_remove_children) removed.
#'  \item ...         Additional parameters.
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



#' @title Generate coordinates of lower triangular matrix of dimension p.
#' @param p          A \code{integer} scalar of dimensionality.
#' @return    A \code{integer} (p*(p-1)/2)-dim vector representing the indices of the lower triangular entries (indices correspond to the vectorization by stacking columns).
#' @export
#' @author
#' @examples
#' TO DO
gen.lower.coordinate=function(p){
  mat=matrix(1:p^2, nrow = p)
  return(mat[lower.tri(mat, diag = F)])  
}