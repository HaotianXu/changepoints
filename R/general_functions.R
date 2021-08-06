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