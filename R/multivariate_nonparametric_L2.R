###############################################################################
## Multivariate normal mixture - density values
## Rcpp version of the function dmvnorm.mixt from R package "ks"
## Parameters
## x - points to compute density at 
## mus - vertical stack of means (row vectors)
## Sigmas - vertical stack of covariance matrices 
## props - vector of mixing proportions (sum(props) be 1)
## Returns
## Density values from the normal mixture (at x)
###############################################################################
#' @noRd
rcpp_dmvnorm_mixt <- function(evalpoints, mus, sigmas, props){
  result = .Call('_changepoints_rcpp_dmvnorm_mixt', PACKAGE = 'changepoints', evalpoints, mus, sigmas, props)
  return(result)
}

#' @noRd
rcpp_biweiker_mixt <- function(evalpoints, mus, bw, props){
  result = .Call('_changepoints_rcpp_biweiker_mixt', PACKAGE = 'changepoints', evalpoints, mus, bw, props)
  return(result)
}

#' @noRd
rcpp_triweiker_mixt <- function(evalpoints, mus, bw, props){
  result = .Call('_changepoints_rcpp_triweiker_mixt', PACKAGE = 'changepoints', evalpoints, mus, bw, props)
  return(result)
}

#' @noRd
rcpp_epanker_mixt <- function(evalpoints, mus, bw, props){
  result = .Call('_changepoints_rcpp_epanker_mixt', PACKAGE = 'changepoints', evalpoints, mus, bw, props)
  return(result)
}

#' @title Multivariate kernel density estimation based on Gaussian kernel.
#' @description Perform multivariate kernel density estimation and evaluated the estimated densities at specified points.
#' @param x            A \code{numeric} matrix of observations with horizontal axis being dimension, and vertical axis being time.
#' @param H            A \code{numeric} (symmetric and positive definite) matrix of bandwidth parameters.
#' @param eval.points  A \code{numeric} matrix of evaluated data points with horizontal axis being dimension, and vertical axis being time..
#' @return   A numeric vector of evaluated densities.
#' @export
#' @author Haotian Xu
#' @examples
#' n = 100
#' p = 10
#' x = matrix(rnorm(n*p), nrow = n)
#' h = 2*(1/n)^{1/(4+p)} # bandwith
#' kde.eval(x, h*diag(p), x)
kde.eval <- function(x, H, eval.points){
  if(is.vector(eval.points)){
    eval.points = matrix(eval.points, nrow = 1)
  }
  n <- nrow(x)
  d <- ncol(x)
  Hs <- replicate(n, H, simplify=FALSE) 
  Hs <- do.call(rbind, Hs)
  fhat <- as.numeric(rcpp_dmvnorm_mixt(evalpoints = eval.points, mus=x, sigmas=Hs, props=rep(1/n, n)))
  return(fhat)
}


#' @title Multivariate kernel density estimation based on Biweight kernel.
#' @description Perform multivariate kernel density estimation and evaluated the estimated densities at specified points.
#' @param x            A \code{numeric} matrix of observations with horizontal axis being dimension, and vertical axis being time.
#' @param h            A \code{numeric} scalar of bandwidth parameter.
#' @param eval.points  A \code{numeric} matrix of evaluated data points with horizontal axis being dimension, and vertical axis being time..
#' @return   A numeric vector of evaluated densities.
#' @export
#' @author Haotian Xu
#' @examples
#' n = 100
#' p = 10
#' x = matrix(rnorm(n*p), nrow = n)
#' h = 2*(1/n)^{1/(4+p)} # bandwith
#' kde.biwei.eval(x, h, x)
kde.biwei.eval <- function(x, bw, eval.points){
  if(is.vector(eval.points)){
    eval.points = matrix(eval.points, nrow = 1)
  }
  n <- nrow(x)
  fhat <- as.numeric(rcpp_biweiker_mixt(evalpoints = eval.points, mus=x, bw=bw, props=rep(1/n, n)))
  return(fhat)
}

#' @title Multivariate kernel density estimation based on Triweight kernel.
#' @description Perform multivariate kernel density estimation and evaluated the estimated densities at specified points.
#' @param x            A \code{numeric} matrix of observations with horizontal axis being dimension, and vertical axis being time.
#' @param h            A \code{numeric} scalar of bandwidth parameter.
#' @param eval.points  A \code{numeric} matrix of evaluated data points with horizontal axis being dimension, and vertical axis being time..
#' @return   A numeric vector of evaluated densities.
#' @export
#' @author Haotian Xu
#' @examples
#' n = 100
#' p = 10
#' x = matrix(rnorm(n*p), nrow = n)
#' h = 2*(1/n)^{1/(4+p)} # bandwith
#' kde.triwei.eval(x, h, x)
kde.triwei.eval <- function(x, bw, eval.points){
  if(is.vector(eval.points)){
    eval.points = matrix(eval.points, nrow = 1)
  }
  n <- nrow(x)
  fhat <- as.numeric(rcpp_triweiker_mixt(evalpoints = eval.points, mus=x, bw=bw, props=rep(1/n, n)))
  return(fhat)
}


#' @title Multivariate kernel density estimation based on Epanechnikov kernel.
#' @description Perform multivariate kernel density estimation and evaluated the estimated densities at specified points.
#' @param x            A \code{numeric} matrix of observations with horizontal axis being dimension, and vertical axis being time.
#' @param h            A \code{numeric} scalar of bandwidth parameter.
#' @param eval.points  A \code{numeric} matrix of evaluated data points with horizontal axis being dimension, and vertical axis being time..
#' @return   A numeric vector of evaluated densities.
#' @export
#' @author Haotian Xu
#' @examples
#' n = 100
#' p = 10
#' x = matrix(rnorm(n*p), nrow = n)
#' h = 2*(1/n)^{1/(4+p)} # bandwith
#' kde.epan.eval(x, h, x)
kde.epan.eval <- function(x, bw, eval.points){
  if(is.vector(eval.points)){
    eval.points = matrix(eval.points, nrow = 1)
  }
  n <- nrow(x)
  fhat <- as.numeric(rcpp_epanker_mixt(evalpoints = eval.points, mus=x, bw=bw, props=rep(1/n, n)))
  return(fhat)
}


#' @title Internal Function: Compute the CUSUM statistic based on L2 distance (multivariate).
#' @param Y             A \code{numeric} matrix of observations with with horizontal axis being time, and vertical axis being dimension.
#' @param s             A \code{integer} scalar of starting index.
#' @param e             A \code{integer} scalar of ending index.
#' @param t             A \code{integer} scalar of splitting index.
#' @param h             A \code{numeric} scalar of bandwidth parameter.
#' @return  A \code{numeric} scalar of the CUSUM statistic based on L2 distance.
#' @noRd
CUSUM.L2.multivariate = function(Y, s, e, t, h){
  p = dim(Y)[1]
  n = dim(Y)[2]
  n_st = t - s + 1
  n_se = e - s + 1
  n_te = e - t
  aux = Y[,s:t]
  temp1 = kde.eval(t(aux), eval.points = t(Y), H = h*diag(p))
  aux = Y[,(t+1):e]
  temp2 = kde.eval(t(aux), eval.points = t(Y), H = h*diag(p))
  result = sqrt(n_st * n_te / n_se) * sqrt(sum((temp1 - temp2)^2))
  return(result)
}


#' @title Internal Function: Compute the error in L2 distance for local refinement.
#' @param Y             A \code{numeric} matrix of observations with with horizontal axis being time, and vertical axis being dimension.
#' @param s             A \code{integer} scalar of starting index.
#' @param e             A \code{integer} scalar of ending index.
#' @param h             A \code{numeric} scalar of bandwidth parameter.
#' @return  A \code{numeric} scalar of error in L2 distance.
#' @noRd
error.L2.multivariate = function(Y, s, e, h){
  p = dim(Y)[1]
  temp_mean = kde.eval(t(Y[,s:e]), eval.points = t(Y), H = h*diag(p))
  error = 0
  for(i in 1:(e-s+1)){
    error = error + sum((kde.eval(t(Y[,s+i-1]), eval.points = t(Y), H = h*diag(p)) - temp_mean)^2)
  }
  return(error)
}



#' @title Wild binary segmentation for multivariate nonparametric change points detection based on L2 distance.
#' @description Perform wild binary segmentation for multivariate nonparametric change points detection based on L2 distance.
#' @param Y         A \code{numeric} matrix of observations with with horizontal axis being time, and vertical axis being dimension.
#' @param s         A \code{integer} scalar of starting index.
#' @param e         A \code{integer} scalar of ending index.
#' @param Alpha     A \code{integer} vector of starting indices of random intervals.
#' @param Beta      A \code{integer} vector of ending indices of random intervals.
#' @param h         A \code{numeric} scalar of bandwidth parameter.
#' @param delta     A \code{integer} scalar of minimum spacing.
#' @param level     Should be fixed as 0.
#' @return   An object of \code{\link[base]{class}} "BS", which is a \code{list} with the following structure:
#'  \item{S}{A vector of estimated change points (sorted in strictly increasing order).}
#'  \item{Dval}{A vector of values of CUSUM statistic based on KS distance.}
#'  \item{Level}{A vector representing the levels at which each change point is detected.}
#'  \item{Parent}{A matrix with the starting indices on the first row and the ending indices on the second row.}
#' @export
#' @author Haotian Xu
#' @seealso \code{\link{thresholdBS}} for obtain change points estimation, \code{\link{tuneBSmultinonpar}} for a tuning version.
#' @examples
#' n = 150
#' v = c(floor(n/3), 2*floor(n/3)) # location of change points
#' r = 2
#' p = 5
#' Y = matrix(0, p, n) # matrix for data
#' mu0 = rep(0, p) # mean of the data
#' mu1 = rep(0, p)
#' mu1[1:floor(p/2)] = 2
#' Sigma0 = diag(p) #Covariance matrices of the data
#' Sigma1 = diag(p)*2
#' # Generate data
#' for(t in 1:n){
#'   if(t < v[1] || t > v[2]){
#'      Y[,t] = MASS::mvrnorm(n = 1, mu0, Sigma0)
#'   }
#'   if(t >= v[1] && t < v[2]){
#'      Y[,t] = MASS::mvrnorm(n = 1, mu1, Sigma1)
#'   }
#' }## close for generate data
#' h = 2*(1/n)^{1/(2*r+p)} # bandwith
#' M = 50
#' intervals = WBS.intervals(M = M, lower = 1, upper = ncol(Y)) #Random intervals
#' temp = WBS.multi.nonpar.L2(Y, 1, ncol(Y), intervals$Alpha, intervals$Beta, h, delta = 15)
#' cpt_init = tuneBSmultinonpar(temp, Y)
WBS.multi.nonpar.L2 = function(Y, s, e, Alpha, Beta, h, delta, level = 0){
  print(paste0("WBS at level: ", level))
  Alpha_new = pmax(Alpha, s)
  Beta_new = pmin(Beta, e)
  idx = which(Beta_new - Alpha_new > 2*delta)
  Alpha_new = Alpha_new[idx]
  Beta_new = Beta_new[idx]
  M = length(Alpha_new)
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
      s_star = Alpha_new[m] + delta
      e_star = Beta_new[m] - delta
      temp = rep(0, e_star - s_star + 1)
      for(t in s_star:e_star){
        temp[t-s_star+1] = CUSUM.L2.multivariate(Y, Alpha_new[m], Beta_new[m], t, h)
      }
      best_value = max(temp)
      best_t = which.max(temp) + s_star - 1
      a[m] = best_value
      b[m] = best_t
    }
    m_star = which.max(a)
  }
  temp1 = WBS.multi.nonpar.L2(Y, s, b[m_star]-1, Alpha, Beta, h, delta, level)
  temp2 = WBS.multi.nonpar.L2(Y, b[m_star], e, Alpha, Beta, h, delta, level)
  S = c(temp1$S, b[m_star], temp2$S)
  Dval = c(temp1$Dval, a[m_star], temp2$Dval)
  Level = c(temp1$Level, level, temp2$Level)
  Parent = cbind(temp1$Parent, parent, temp2$Parent)
  result = list(S = S, Dval = Dval, Level = Level, Parent = Parent)
  class(result) = "BS"
  return(result)
}


#' @title Local refinement for multivariate nonparametric change points localisation based on L2 distance.
#' @description     Perform local refinement for multivariate nonparametric change points localisation based on L2 distance.
#' @param cpt_init  An \code{integer} vector of initial change points estimation (sorted in strictly increasing order).
#' @param Y         A \code{numeric} matrix of observations with with horizontal axis being time, and vertical axis being dimension.
#' @param kappa_hat A \code{numeric} vector of jump sizes estimator.
#' @param r         An \code{integer} scalar of smoothness parameter of the underlying Holder space.
#' @param w         A \code{numeric} scalar in (0,1) representing the weight for interval truncation.
#' @param c_kappa   A \code{numeric} scalar to be multiplied by kappa estimator as the bandwidth in the local refinment.
#' @return  A vector of locally refined change points estimation.
#' @export
#' @author Haotian Xu
#' @examples
#' n = 150
#' v = c(floor(n/3), 2*floor(n/3)) # location of change points
#' r = 2
#' p = 6
#' Y = matrix(0, p, n) # matrix for data
#' mu0 = rep(0, p) # mean of the data
#' mu1 = rep(0, p)
#' mu1[1:floor(p/2)] = 2
#' Sigma0 = diag(p) #Covariance matrices of the data
#' Sigma1 = diag(p)
#' # Generate data
#' for(t in 1:n){
#'   if(t <= v[1] || t > v[2]){
#'      Y[,t] = MASS::mvrnorm(n = 1, mu0, Sigma0)
#'   }
#'   if(t > v[1] && t <= v[2]){
#'      Y[,t] = MASS::mvrnorm(n = 1, mu1, Sigma1)
#'   }
#' }## close for generate data
#' M = 100
#' intervals = WBS.intervals(M = M, lower = 1, upper = ncol(Y)) #Random intervals
#' h = 2*(1/n)^{1/(2*r+p)} # bandwith
#' temp = WBS.multi.nonpar.L2(Y, 1, ncol(Y), intervals$Alpha, intervals$Beta, h, delta = 15)
#' cpt_init = tuneBSmultinonpar(temp, Y)
#' kappa_hat = kappa.multi.nonpar.L2(cpt_init, Y, h_kappa = 0.01)
#' local.refine.multi.nonpar.L2(cpt_init, Y, kappa_hat, r = 2, w = 0.9, c_kappa = 2)
local.refine.multi.nonpar.L2 = function(cpt_init, Y, kappa_hat, r = 2, w = 0.9, c_kappa = 10){
  p = dim(Y)[1]
  n = dim(Y)[2]
  cpt_init_long = c(0, cpt_init, n)
  cpt_init_numb = length(cpt_init)
  cpt_refined = rep(0, cpt_init_numb+1)
  h_refine = rep(NA, cpt_init_numb)
  minmax_mat = apply(Y, MARGIN = 1, function(x) c(min(x), max(x)))
  for (k in 1:cpt_init_numb){
    #aux1 = Y[,(cpt_init_long[k]+1):(cpt_init_long[k+1])]
    #aux2 = Y[,(cpt_init_long[k+1]+1):(cpt_init_long[k+2])]
    #integrateFCT = function(x){(kde.eval(t(aux1), eval.points = x, H = h_init*diag(p)) - kde.eval(t(aux2), eval.points = x, H = h_init*diag(p)))^2}
    #temp_integrate = cubature::cubintegrate(integrateFCT, lower = minmax_mat[1,], upper = minmax_mat[2,], method = "suave", maxEval = 1000)$integral
    #kappa_hat[k] = sqrt(temp_integrate)/sqrt((cpt_init_long[k+2]-cpt_init_long[k+1])*(cpt_init_long[k+1]-cpt_init_long[k])/(cpt_init_long[k+2]-cpt_init_long[k]))
    h_refine[k] = c_kappa*(kappa_hat[k])^(1/r)
    s = w*cpt_init_long[k] + (1-w)*cpt_init_long[k+1]
    e = (1-w)*cpt_init_long[k+1] + w*cpt_init_long[k+2]
    lower = ceiling(s) + 2
    upper = floor(e) - 2
    b = sapply(lower:upper, function(eta)(error.L2.multivariate(Y, ceiling(s), eta, h_refine[k]) + error.L2.multivariate(Y, (eta+1), floor(e), h_refine[k])))
    cpt_refined[k+1] = ceiling(s) + which.min(b)
  }
  return(list(cpt_refined = cpt_refined[-1]+1, kappa_hat = kappa_hat, h_refine = h_refine))
}



#' @export
jumpsize.multinonpar.L2 = function(cpt_init, Y, h_kappa = 0.01, c_fac){
  p = dim(Y)[1]
  n = dim(Y)[2]
  cpt_init_long = c(0, cpt_init, n)
  cpt_init_numb = length(cpt_init)
  kappa_hat = rep(NA, cpt_init_numb)
  minmax_mat = apply(Y, MARGIN = 1, function(x) c(min(x), max(x)))
  for (k in 1:cpt_init_numb){
    aux1 = Y[,(cpt_init_long[k]+1):(cpt_init_long[k+1])]
    aux2 = Y[,(cpt_init_long[k+1]+1):(cpt_init_long[k+2])]
    integrateFCT = function(x){(kde.eval(t(aux1), eval.points = x, H = h_kappa*diag(p)) - kde.eval(t(aux2), eval.points = x, H = h_kappa*diag(p)))^2}
    temp_integrate = cubature::cubintegrate(integrateFCT, lower = minmax_mat[1,], upper = minmax_mat[2,], method = "suave", maxEval = 1000)$integral
    kappa_hat[k] = c_fac*sqrt(temp_integrate)
  }
  return(kappa_hat)
}


#' @export
LRV.multinonpar.L2 = function(cpt_init, kappa_hat, Y, r = 2, w = 0.9, c_lrv, block_size){
  n = ncol(Y)
  p = nrow(Y)
  cpt_init_long = c(0, cpt_init, n)
  cpt_init_numb = length(cpt_init)
  lrv_hat = rep(NA, cpt_init_numb)
  minmax_mat = apply(Y, MARGIN = 1, function(x) c(min(x), max(x)))
  volume = prod(minmax_mat[2,] - minmax_mat[1,])
  for (k in 1:cpt_init_numb){
    kappahat = kappa_hat[k]
    h2 = c_lrv[k]*kappahat^(1/r)
    s = w*cpt_init_long[k] + (1-w)*cpt_init_long[k+1]
    e = (1-w)*cpt_init_long[k+1] + w*cpt_init_long[k+2]
    mean1 = kde.eval(t(Y[,s:cpt_init_long[k+1]]), eval.points = t(Y), H = h2*diag(p))
    mean2 = kde.eval(t(Y[,(cpt_init_long[k+1]+1):e]), eval.points = t(Y), H = h2*diag(p))
    z_vec = kappahat^(p/(2*r)-1)*volume*apply(cbind(sapply(1:(cpt_init_long[k+1]-s+1), function(i){(kde.eval(t(Y[,s+i-1]), eval.points = t(Y), H = h2*diag(p)) - mean1)*(mean1 - mean2)}),
                                                    sapply((cpt_init_long[k+1]-s+2):(e-s+1), function(i){(kde.eval(t(Y[,s+i-1]), eval.points = t(Y), H = h2*diag(p)) - mean2)*(mean1 - mean2)})), MARGIN = 2, mean)
    numb = floor((floor(e)-ceiling(s)+1)/(block_size))
    z_mat1 = matrix(z_vec[1:(block_size*numb)], nrow = block_size)
    z_mat1_colsum = apply(z_mat1, 2, sum)
    z_mat2 = matrix(rev(z_vec)[1:(block_size*numb)], nrow = block_size)
    z_mat2_colsum = apply(z_mat2, 2, sum)
    lrv_hat[k] = (mean((z_mat1_colsum)^2/block_size) + mean((z_mat2_colsum)^2/block_size))/2
  }
  return(lrv_hat)
}


##############################################################################################
#' @title Internal Function: Compute the CUSUM statistic based on L2 distance (multivariate).
#' @param Y             A \code{numeric} matrix of observations with with horizontal axis being time, and vertical axis being dimension.
#' @param s             A \code{integer} scalar of starting index.
#' @param e             A \code{integer} scalar of ending index.
#' @param t             A \code{integer} scalar of splitting index.
#' @param h             A \code{numeric} scalar of bandwidth parameter.
#' @return  A \code{numeric} scalar of the CUSUM statistic based on L2 distance.
#' @noRd
CUSUM.L2.multivariate.epan = function(Y, s, e, t, h){
  p = dim(Y)[1]
  n = dim(Y)[2]
  n_st = t - s + 1
  n_se = e - s + 1
  n_te = e - t
  aux = Y[,s:t]
  temp1 = kde.epan.eval(t(aux), bw = h, eval.points = t(Y))
  aux = Y[,(t+1):e]
  temp2 = kde.epan.eval(t(aux), bw = h, eval.points = t(Y))
  result = sqrt(n_st * n_te / n_se) * sqrt(sum((temp1 - temp2)^2))
  return(result)
}


#' @title Internal Function: Compute the error in L2 distance for local refinement.
#' @param Y             A \code{numeric} matrix of observations with with horizontal axis being time, and vertical axis being dimension.
#' @param s             A \code{integer} scalar of starting index.
#' @param e             A \code{integer} scalar of ending index.
#' @param h             A \code{numeric} scalar of bandwidth parameter.
#' @return  A \code{numeric} scalar of error in L2 distance.
#' @noRd
error.L2.multivariate.epan = function(Y, s, e, h){
  p = dim(Y)[1]
  temp_mean = kde.epan.eval(t(Y[,s:e]), bw = h, eval.points = t(Y))
  error = 0
  for(i in 1:(e-s+1)){
    error = error + sum((kde.epan.eval(t(Y[,s+i-1]), bw = h, eval.points = t(Y)) - temp_mean)^2)
  }
  return(error)
}



#' @title Wild binary segmentation for multivariate nonparametric change points detection based on L2 distance.
#' @description Perform wild binary segmentation for multivariate nonparametric change points detection based on L2 distance.
#' @param Y         A \code{numeric} matrix of observations with with horizontal axis being time, and vertical axis being dimension.
#' @param s         A \code{integer} scalar of starting index.
#' @param e         A \code{integer} scalar of ending index.
#' @param Alpha     A \code{integer} vector of starting indices of random intervals.
#' @param Beta      A \code{integer} vector of ending indices of random intervals.
#' @param h         A \code{numeric} scalar of bandwidth parameter.
#' @param delta     A \code{integer} scalar of minimum spacing.
#' @param level     Should be fixed as 0.
#' @return   An object of \code{\link[base]{class}} "BS", which is a \code{list} with the following structure:
#'  \item{S}{A vector of estimated change points (sorted in strictly increasing order).}
#'  \item{Dval}{A vector of values of CUSUM statistic based on KS distance.}
#'  \item{Level}{A vector representing the levels at which each change point is detected.}
#'  \item{Parent}{A matrix with the starting indices on the first row and the ending indices on the second row.}
#' @export
#' @author Haotian Xu
#' @seealso \code{\link{thresholdBS}} for obtain change points estimation, \code{\link{tuneBSmultinonpar}} for a tuning version.
#' @examples
#' n = 150
#' v = c(floor(n/3), 2*floor(n/3)) # location of change points
#' r = 2
#' p = 5
#' Y = matrix(0, p, n) # matrix for data
#' mu0 = rep(0, p) # mean of the data
#' mu1 = rep(0, p)
#' mu1[1:floor(p/2)] = 2
#' Sigma0 = diag(p) #Covariance matrices of the data
#' Sigma1 = diag(p)*2
#' # Generate data
#' for(t in 1:n){
#'   if(t < v[1] || t > v[2]){
#'      Y[,t] = MASS::mvrnorm(n = 1, mu0, Sigma0)
#'   }
#'   if(t >= v[1] && t < v[2]){
#'      Y[,t] = MASS::mvrnorm(n = 1, mu1, Sigma1)
#'   }
#' }## close for generate data
#' h = 2*(1/n)^{1/(2*r+p)} # bandwith
#' M = 50
#' intervals = WBS.intervals(M = M, lower = 1, upper = ncol(Y)) #Random intervals
#' temp = WBS.multi.nonpar.L2.epan(Y, 1, ncol(Y), intervals$Alpha, intervals$Beta, h, delta = 15)
#' cpt_init = tuneBSmultinonpar(temp, Y)
WBS.multi.nonpar.L2.epan = function(Y, s, e, Alpha, Beta, h, delta, level = 0){
  print(paste0("WBS at level: ", level))
  Alpha_new = pmax(Alpha, s)
  Beta_new = pmin(Beta, e)
  idx = which(Beta_new - Alpha_new > 2*delta)
  Alpha_new = Alpha_new[idx]
  Beta_new = Beta_new[idx]
  M = length(Alpha_new)
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
      s_star = Alpha_new[m] + delta
      e_star = Beta_new[m] - delta
      temp = rep(0, e_star - s_star + 1)
      for(t in s_star:e_star){
        temp[t-s_star+1] = CUSUM.L2.multivariate.epan(Y, Alpha_new[m], Beta_new[m], t, h)
      }
      best_value = max(temp)
      best_t = which.max(temp) + s_star - 1
      a[m] = best_value
      b[m] = best_t
    }
    m_star = which.max(a)
  }
  temp1 = WBS.multi.nonpar.L2.epan(Y, s, b[m_star]-1, Alpha, Beta, h, delta, level)
  temp2 = WBS.multi.nonpar.L2.epan(Y, b[m_star], e, Alpha, Beta, h, delta, level)
  S = c(temp1$S, b[m_star], temp2$S)
  Dval = c(temp1$Dval, a[m_star], temp2$Dval)
  Level = c(temp1$Level, level, temp2$Level)
  Parent = cbind(temp1$Parent, parent, temp2$Parent)
  result = list(S = S, Dval = Dval, Level = Level, Parent = Parent)
  class(result) = "BS"
  return(result)
}


#' @title Local refinement for multivariate nonparametric change points localisation based on L2 distance.
#' @description     Perform local refinement for multivariate nonparametric change points localisation based on L2 distance.
#' @param cpt_init  An \code{integer} vector of initial change points estimation (sorted in strictly increasing order).
#' @param Y         A \code{numeric} matrix of observations with with horizontal axis being time, and vertical axis being dimension.
#' @param kappa_hat A \code{numeric} vector of jump sizes estimator.
#' @param r         An \code{integer} scalar of smoothness parameter of the underlying Holder space.
#' @param w         A \code{numeric} scalar in (0,1) representing the weight for interval truncation.
#' @param c_kappa   A \code{numeric} scalar to be multiplied by kappa estimator as the bandwidth in the local refinment.
#' @return  A vector of locally refined change points estimation.
#' @export
#' @author Haotian Xu
#' @examples
#' n = 150
#' v = c(floor(n/3), 2*floor(n/3)) # location of change points
#' r = 2
#' p = 6
#' Y = matrix(0, p, n) # matrix for data
#' mu0 = rep(0, p) # mean of the data
#' mu1 = rep(0, p)
#' mu1[1:floor(p/2)] = 2
#' Sigma0 = diag(p) #Covariance matrices of the data
#' Sigma1 = diag(p)
#' # Generate data
#' for(t in 1:n){
#'   if(t <= v[1] || t > v[2]){
#'      Y[,t] = MASS::mvrnorm(n = 1, mu0, Sigma0)
#'   }
#'   if(t > v[1] && t <= v[2]){
#'      Y[,t] = MASS::mvrnorm(n = 1, mu1, Sigma1)
#'   }
#' }## close for generate data
#' M = 100
#' intervals = WBS.intervals(M = M, lower = 1, upper = ncol(Y)) #Random intervals
#' h = 2*(1/n)^{1/(2*r+p)} # bandwith
#' temp = WBS.multi.nonpar.L2(Y, 1, ncol(Y), intervals$Alpha, intervals$Beta, h, delta = 15)
#' cpt_init = tuneBSmultinonpar(temp, Y)
#' kappa_hat = kappa.multi.nonpar.L2(cpt_init, Y, h_kappa = 0.01)
#' local.refine.multi.nonpar.L2(cpt_init, Y, kappa_hat, r = 2, w = 0.9, c_kappa = 2)
local.refine.multi.nonpar.L2.epan = function(cpt_init, Y, kappa_hat, r = 2, w = 0.9, c_kappa = 10){
  p = dim(Y)[1]
  n = dim(Y)[2]
  cpt_init_long = c(0, cpt_init, n)
  cpt_init_numb = length(cpt_init)
  cpt_refined = rep(0, cpt_init_numb+1)
  h_refine = rep(NA, cpt_init_numb)
  minmax_mat = apply(Y, MARGIN = 1, function(x) c(min(x), max(x)))
  for (k in 1:cpt_init_numb){
    #aux1 = Y[,(cpt_init_long[k]+1):(cpt_init_long[k+1])]
    #aux2 = Y[,(cpt_init_long[k+1]+1):(cpt_init_long[k+2])]
    #integrateFCT = function(x){(kde.eval(t(aux1), eval.points = x, H = h_init*diag(p)) - kde.eval(t(aux2), eval.points = x, H = h_init*diag(p)))^2}
    #temp_integrate = cubature::cubintegrate(integrateFCT, lower = minmax_mat[1,], upper = minmax_mat[2,], method = "suave", maxEval = 1000)$integral
    #kappa_hat[k] = sqrt(temp_integrate)/sqrt((cpt_init_long[k+2]-cpt_init_long[k+1])*(cpt_init_long[k+1]-cpt_init_long[k])/(cpt_init_long[k+2]-cpt_init_long[k]))
    h_refine[k] = c_kappa*(kappa_hat[k])^(1/r)
    s = w*cpt_init_long[k] + (1-w)*cpt_init_long[k+1]
    e = (1-w)*cpt_init_long[k+1] + w*cpt_init_long[k+2]
    lower = ceiling(s) + 2
    upper = floor(e) - 2
    b = sapply(lower:upper, function(eta)(error.L2.multivariate.epan(Y, ceiling(s), eta, h_refine[k]) + error.L2.multivariate.epan(Y, (eta+1), floor(e), h_refine[k])))
    cpt_refined[k+1] = ceiling(s) + which.min(b)
  }
  return(list(cpt_refined = cpt_refined[-1]+1, kappa_hat = kappa_hat, h_refine = h_refine))
}


#' @export
jumpsize.multinonpar.L2.epan = function(cpt_init, Y, h_kappa = 0.01, c_fac){
  p = dim(Y)[1]
  n = dim(Y)[2]
  cpt_init_long = c(0, cpt_init, n)
  cpt_init_numb = length(cpt_init)
  kappa_hat = rep(NA, cpt_init_numb)
  minmax_mat = apply(Y, MARGIN = 1, function(x) c(min(x), max(x)))
  for (k in 1:cpt_init_numb){
    aux1 = Y[,(cpt_init_long[k]+1):(cpt_init_long[k+1])]
    aux2 = Y[,(cpt_init_long[k+1]+1):(cpt_init_long[k+2])]
    integrateFCT = function(x){(kde.epan.eval(t(aux1), bw = h_kappa, eval.points = x) - kde.epan.eval(t(aux2), bw = h_kappa, eval.points = x))^2}
    temp_integrate = cubature::cubintegrate(integrateFCT, lower = minmax_mat[1,], upper = minmax_mat[2,], method = "suave", maxEval = 1000)$integral
    kappa_hat[k] = c_fac*sqrt(temp_integrate)
  }
  return(kappa_hat)
}



##############################################################################################
#' @title Internal Function: Compute the CUSUM statistic based on L2 distance (multivariate).
#' @param Y             A \code{numeric} matrix of observations with with horizontal axis being time, and vertical axis being dimension.
#' @param s             A \code{integer} scalar of starting index.
#' @param e             A \code{integer} scalar of ending index.
#' @param t             A \code{integer} scalar of splitting index.
#' @param h             A \code{numeric} scalar of bandwidth parameter.
#' @return  A \code{numeric} scalar of the CUSUM statistic based on L2 distance.
#' @noRd
CUSUM.L2.multivariate.biwei = function(Y, s, e, t, h){
  p = dim(Y)[1]
  n = dim(Y)[2]
  n_st = t - s + 1
  n_se = e - s + 1
  n_te = e - t
  aux = Y[,s:t]
  temp1 = kde.biwei.eval(t(aux), bw = h, eval.points = t(Y))
  aux = Y[,(t+1):e]
  temp2 = kde.biwei.eval(t(aux), bw = h, eval.points = t(Y))
  result = sqrt(n_st * n_te / n_se) * sqrt(sum((temp1 - temp2)^2))
  return(result)
}


#' @title Internal Function: Compute the error in L2 distance for local refinement.
#' @param Y             A \code{numeric} matrix of observations with with horizontal axis being time, and vertical axis being dimension.
#' @param s             A \code{integer} scalar of starting index.
#' @param e             A \code{integer} scalar of ending index.
#' @param h             A \code{numeric} scalar of bandwidth parameter.
#' @return  A \code{numeric} scalar of error in L2 distance.
#' @noRd
error.L2.multivariate.biwei = function(Y, s, e, h){
  p = dim(Y)[1]
  temp_mean = kde.biwei.eval(t(Y[,s:e]), bw = h, eval.points = t(Y))
  error = 0
  for(i in 1:(e-s+1)){
    error = error + sum((kde.biwei.eval(t(Y[,s+i-1]), bw = h, eval.points = t(Y)) - temp_mean)^2)
  }
  return(error)
}



#' @title Wild binary segmentation for multivariate nonparametric change points detection based on L2 distance.
#' @description Perform wild binary segmentation for multivariate nonparametric change points detection based on L2 distance.
#' @param Y         A \code{numeric} matrix of observations with with horizontal axis being time, and vertical axis being dimension.
#' @param s         A \code{integer} scalar of starting index.
#' @param e         A \code{integer} scalar of ending index.
#' @param Alpha     A \code{integer} vector of starting indices of random intervals.
#' @param Beta      A \code{integer} vector of ending indices of random intervals.
#' @param h         A \code{numeric} scalar of bandwidth parameter.
#' @param delta     A \code{integer} scalar of minimum spacing.
#' @param level     Should be fixed as 0.
#' @return   An object of \code{\link[base]{class}} "BS", which is a \code{list} with the following structure:
#'  \item{S}{A vector of estimated change points (sorted in strictly increasing order).}
#'  \item{Dval}{A vector of values of CUSUM statistic based on KS distance.}
#'  \item{Level}{A vector representing the levels at which each change point is detected.}
#'  \item{Parent}{A matrix with the starting indices on the first row and the ending indices on the second row.}
#' @export
#' @author Haotian Xu
#' @seealso \code{\link{thresholdBS}} for obtain change points estimation, \code{\link{tuneBSmultinonpar}} for a tuning version.
#' @examples
#' n = 150
#' v = c(floor(n/3), 2*floor(n/3)) # location of change points
#' r = 2
#' p = 5
#' Y = matrix(0, p, n) # matrix for data
#' mu0 = rep(0, p) # mean of the data
#' mu1 = rep(0, p)
#' mu1[1:floor(p/2)] = 2
#' Sigma0 = diag(p) #Covariance matrices of the data
#' Sigma1 = diag(p)*2
#' # Generate data
#' for(t in 1:n){
#'   if(t < v[1] || t > v[2]){
#'      Y[,t] = MASS::mvrnorm(n = 1, mu0, Sigma0)
#'   }
#'   if(t >= v[1] && t < v[2]){
#'      Y[,t] = MASS::mvrnorm(n = 1, mu1, Sigma1)
#'   }
#' }## close for generate data
#' h = 2*(1/n)^{1/(2*r+p)} # bandwith
#' M = 50
#' intervals = WBS.intervals(M = M, lower = 1, upper = ncol(Y)) #Random intervals
#' temp = WBS.multi.nonpar.L2.biwei(Y, 1, ncol(Y), intervals$Alpha, intervals$Beta, h, delta = 15)
#' cpt_init = tuneBSmultinonpar(temp, Y)
WBS.multi.nonpar.L2.biwei = function(Y, s, e, Alpha, Beta, h, delta, level = 0){
  print(paste0("WBS at level: ", level))
  Alpha_new = pmax(Alpha, s)
  Beta_new = pmin(Beta, e)
  idx = which(Beta_new - Alpha_new > 2*delta)
  Alpha_new = Alpha_new[idx]
  Beta_new = Beta_new[idx]
  M = length(Alpha_new)
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
      s_star = Alpha_new[m] + delta
      e_star = Beta_new[m] - delta
      temp = rep(0, e_star - s_star + 1)
      for(t in s_star:e_star){
        temp[t-s_star+1] = CUSUM.L2.multivariate.biwei(Y, Alpha_new[m], Beta_new[m], t, h)
      }
      best_value = max(temp)
      best_t = which.max(temp) + s_star - 1
      a[m] = best_value
      b[m] = best_t
    }
    m_star = which.max(a)
  }
  temp1 = WBS.multi.nonpar.L2.biwei(Y, s, b[m_star]-1, Alpha, Beta, h, delta, level)
  temp2 = WBS.multi.nonpar.L2.biwei(Y, b[m_star], e, Alpha, Beta, h, delta, level)
  S = c(temp1$S, b[m_star], temp2$S)
  Dval = c(temp1$Dval, a[m_star], temp2$Dval)
  Level = c(temp1$Level, level, temp2$Level)
  Parent = cbind(temp1$Parent, parent, temp2$Parent)
  result = list(S = S, Dval = Dval, Level = Level, Parent = Parent)
  class(result) = "BS"
  return(result)
}


#' @title Local refinement for multivariate nonparametric change points localisation based on L2 distance.
#' @description     Perform local refinement for multivariate nonparametric change points localisation based on L2 distance.
#' @param cpt_init  An \code{integer} vector of initial change points estimation (sorted in strictly increasing order).
#' @param Y         A \code{numeric} matrix of observations with with horizontal axis being time, and vertical axis being dimension.
#' @param kappa_hat A \code{numeric} vector of jump sizes estimator.
#' @param r         An \code{integer} scalar of smoothness parameter of the underlying Holder space.
#' @param w         A \code{numeric} scalar in (0,1) representing the weight for interval truncation.
#' @param c_kappa   A \code{numeric} scalar to be multiplied by kappa estimator as the bandwidth in the local refinment.
#' @return  A vector of locally refined change points estimation.
#' @export
#' @author Haotian Xu
#' @examples
#' n = 150
#' v = c(floor(n/3), 2*floor(n/3)) # location of change points
#' r = 2
#' p = 6
#' Y = matrix(0, p, n) # matrix for data
#' mu0 = rep(0, p) # mean of the data
#' mu1 = rep(0, p)
#' mu1[1:floor(p/2)] = 2
#' Sigma0 = diag(p) #Covariance matrices of the data
#' Sigma1 = diag(p)
#' # Generate data
#' for(t in 1:n){
#'   if(t <= v[1] || t > v[2]){
#'      Y[,t] = MASS::mvrnorm(n = 1, mu0, Sigma0)
#'   }
#'   if(t > v[1] && t <= v[2]){
#'      Y[,t] = MASS::mvrnorm(n = 1, mu1, Sigma1)
#'   }
#' }## close for generate data
#' M = 100
#' intervals = WBS.intervals(M = M, lower = 1, upper = ncol(Y)) #Random intervals
#' h = 2*(1/n)^{1/(2*r+p)} # bandwith
#' temp = WBS.multi.nonpar.L2(Y, 1, ncol(Y), intervals$Alpha, intervals$Beta, h, delta = 15)
#' cpt_init = tuneBSmultinonpar(temp, Y)
#' kappa_hat = kappa.multi.nonpar.L2(cpt_init, Y, h_kappa = 0.01)
#' local.refine.multi.nonpar.L2(cpt_init, Y, kappa_hat, r = 2, w = 0.9, c_kappa = 2)
local.refine.multi.nonpar.L2.biwei = function(cpt_init, Y, kappa_hat, r = 2, w = 0.9, c_kappa = 10){
  p = dim(Y)[1]
  n = dim(Y)[2]
  cpt_init_long = c(0, cpt_init, n)
  cpt_init_numb = length(cpt_init)
  cpt_refined = rep(0, cpt_init_numb+1)
  h_refine = rep(NA, cpt_init_numb)
  minmax_mat = apply(Y, MARGIN = 1, function(x) c(min(x), max(x)))
  for (k in 1:cpt_init_numb){
    #aux1 = Y[,(cpt_init_long[k]+1):(cpt_init_long[k+1])]
    #aux2 = Y[,(cpt_init_long[k+1]+1):(cpt_init_long[k+2])]
    #integrateFCT = function(x){(kde.eval(t(aux1), eval.points = x, H = h_init*diag(p)) - kde.eval(t(aux2), eval.points = x, H = h_init*diag(p)))^2}
    #temp_integrate = cubature::cubintegrate(integrateFCT, lower = minmax_mat[1,], upper = minmax_mat[2,], method = "suave", maxEval = 1000)$integral
    #kappa_hat[k] = sqrt(temp_integrate)/sqrt((cpt_init_long[k+2]-cpt_init_long[k+1])*(cpt_init_long[k+1]-cpt_init_long[k])/(cpt_init_long[k+2]-cpt_init_long[k]))
    h_refine[k] = c_kappa*(kappa_hat[k])^(1/r)
    s = w*cpt_init_long[k] + (1-w)*cpt_init_long[k+1]
    e = (1-w)*cpt_init_long[k+1] + w*cpt_init_long[k+2]
    lower = ceiling(s) + 2
    upper = floor(e) - 2
    b = sapply(lower:upper, function(eta)(error.L2.multivariate.biwei(Y, ceiling(s), eta, h_refine[k]) + error.L2.multivariate.biwei(Y, (eta+1), floor(e), h_refine[k])))
    cpt_refined[k+1] = ceiling(s) + which.min(b)
  }
  return(list(cpt_refined = cpt_refined[-1]+1, kappa_hat = kappa_hat, h_refine = h_refine))
}


#' @export
jumpsize.multinonpar.L2.biwei = function(cpt_init, Y, h_kappa = 0.01, c_fac){
  p = dim(Y)[1]
  n = dim(Y)[2]
  cpt_init_long = c(0, cpt_init, n)
  cpt_init_numb = length(cpt_init)
  kappa_hat = rep(NA, cpt_init_numb)
  minmax_mat = apply(Y, MARGIN = 1, function(x) c(min(x), max(x)))
  for (k in 1:cpt_init_numb){
    aux1 = Y[,(cpt_init_long[k]+1):(cpt_init_long[k+1])]
    aux2 = Y[,(cpt_init_long[k+1]+1):(cpt_init_long[k+2])]
    integrateFCT = function(x){(kde.biwei.eval(t(aux1), bw = h_kappa, eval.points = x) - kde.biwei.eval(t(aux2), bw = h_kappa, eval.points = x))^2}
    temp_integrate = cubature::cubintegrate(integrateFCT, lower = minmax_mat[1,], upper = minmax_mat[2,], method = "suave", maxEval = 1000)$integral
    kappa_hat[k] = c_fac*sqrt(temp_integrate)
  }
  return(kappa_hat)
}



##############################################################################################
#' @title Internal Function: Compute the CUSUM statistic based on L2 distance (multivariate).
#' @param Y             A \code{numeric} matrix of observations with with horizontal axis being time, and vertical axis being dimension.
#' @param s             A \code{integer} scalar of starting index.
#' @param e             A \code{integer} scalar of ending index.
#' @param t             A \code{integer} scalar of splitting index.
#' @param h             A \code{numeric} scalar of bandwidth parameter.
#' @return  A \code{numeric} scalar of the CUSUM statistic based on L2 distance.
#' @noRd
CUSUM.L2.multivariate.triwei = function(Y, s, e, t, h){
  p = dim(Y)[1]
  n = dim(Y)[2]
  n_st = t - s + 1
  n_se = e - s + 1
  n_te = e - t
  aux = Y[,s:t]
  temp1 = kde.triwei.eval(t(aux), bw = h, eval.points = t(Y))
  aux = Y[,(t+1):e]
  temp2 = kde.triwei.eval(t(aux), bw = h, eval.points = t(Y))
  result = sqrt(n_st * n_te / n_se) * sqrt(sum((temp1 - temp2)^2))
  return(result)
}


#' @title Internal Function: Compute the error in L2 distance for local refinement.
#' @param Y             A \code{numeric} matrix of observations with with horizontal axis being time, and vertical axis being dimension.
#' @param s             A \code{integer} scalar of starting index.
#' @param e             A \code{integer} scalar of ending index.
#' @param h             A \code{numeric} scalar of bandwidth parameter.
#' @return  A \code{numeric} scalar of error in L2 distance.
#' @noRd
error.L2.multivariate.triwei = function(Y, s, e, h){
  p = dim(Y)[1]
  temp_mean = kde.triwei.eval(t(Y[,s:e]), bw = h, eval.points = t(Y))
  error = 0
  for(i in 1:(e-s+1)){
    error = error + sum((kde.triwei.eval(t(Y[,s+i-1]), bw = h, eval.points = t(Y)) - temp_mean)^2)
  }
  return(error)
}



#' @title Wild binary segmentation for multivariate nonparametric change points detection based on L2 distance.
#' @description Perform wild binary segmentation for multivariate nonparametric change points detection based on L2 distance.
#' @param Y         A \code{numeric} matrix of observations with with horizontal axis being time, and vertical axis being dimension.
#' @param s         A \code{integer} scalar of starting index.
#' @param e         A \code{integer} scalar of ending index.
#' @param Alpha     A \code{integer} vector of starting indices of random intervals.
#' @param Beta      A \code{integer} vector of ending indices of random intervals.
#' @param h         A \code{numeric} scalar of bandwidth parameter.
#' @param delta     A \code{integer} scalar of minimum spacing.
#' @param level     Should be fixed as 0.
#' @return   An object of \code{\link[base]{class}} "BS", which is a \code{list} with the following structure:
#'  \item{S}{A vector of estimated change points (sorted in strictly increasing order).}
#'  \item{Dval}{A vector of values of CUSUM statistic based on KS distance.}
#'  \item{Level}{A vector representing the levels at which each change point is detected.}
#'  \item{Parent}{A matrix with the starting indices on the first row and the ending indices on the second row.}
#' @export
#' @author Haotian Xu
#' @seealso \code{\link{thresholdBS}} for obtain change points estimation, \code{\link{tuneBSmultinonpar}} for a tuning version.
#' @examples
#' n = 150
#' v = c(floor(n/3), 2*floor(n/3)) # location of change points
#' r = 2
#' p = 5
#' Y = matrix(0, p, n) # matrix for data
#' mu0 = rep(0, p) # mean of the data
#' mu1 = rep(0, p)
#' mu1[1:floor(p/2)] = 2
#' Sigma0 = diag(p) #Covariance matrices of the data
#' Sigma1 = diag(p)*2
#' # Generate data
#' for(t in 1:n){
#'   if(t < v[1] || t > v[2]){
#'      Y[,t] = MASS::mvrnorm(n = 1, mu0, Sigma0)
#'   }
#'   if(t >= v[1] && t < v[2]){
#'      Y[,t] = MASS::mvrnorm(n = 1, mu1, Sigma1)
#'   }
#' }## close for generate data
#' h = 2*(1/n)^{1/(2*r+p)} # bandwith
#' M = 50
#' intervals = WBS.intervals(M = M, lower = 1, upper = ncol(Y)) #Random intervals
#' temp = WBS.multi.nonpar.L2.triwei(Y, 1, ncol(Y), intervals$Alpha, intervals$Beta, h, delta = 15)
#' cpt_init = tuneBSmultinonpar(temp, Y)
WBS.multi.nonpar.L2.triwei = function(Y, s, e, Alpha, Beta, h, delta, level = 0){
  print(paste0("WBS at level: ", level))
  Alpha_new = pmax(Alpha, s)
  Beta_new = pmin(Beta, e)
  idx = which(Beta_new - Alpha_new > 2*delta)
  Alpha_new = Alpha_new[idx]
  Beta_new = Beta_new[idx]
  M = length(Alpha_new)
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
      s_star = Alpha_new[m] + delta
      e_star = Beta_new[m] - delta
      temp = rep(0, e_star - s_star + 1)
      for(t in s_star:e_star){
        temp[t-s_star+1] = CUSUM.L2.multivariate.triwei(Y, Alpha_new[m], Beta_new[m], t, h)
      }
      best_value = max(temp)
      best_t = which.max(temp) + s_star - 1
      a[m] = best_value
      b[m] = best_t
    }
    m_star = which.max(a)
  }
  temp1 = WBS.multi.nonpar.L2.triwei(Y, s, b[m_star]-1, Alpha, Beta, h, delta, level)
  temp2 = WBS.multi.nonpar.L2.triwei(Y, b[m_star], e, Alpha, Beta, h, delta, level)
  S = c(temp1$S, b[m_star], temp2$S)
  Dval = c(temp1$Dval, a[m_star], temp2$Dval)
  Level = c(temp1$Level, level, temp2$Level)
  Parent = cbind(temp1$Parent, parent, temp2$Parent)
  result = list(S = S, Dval = Dval, Level = Level, Parent = Parent)
  class(result) = "BS"
  return(result)
}


#' @title Local refinement for multivariate nonparametric change points localisation based on L2 distance.
#' @description     Perform local refinement for multivariate nonparametric change points localisation based on L2 distance.
#' @param cpt_init  An \code{integer} vector of initial change points estimation (sorted in strictly increasing order).
#' @param Y         A \code{numeric} matrix of observations with with horizontal axis being time, and vertical axis being dimension.
#' @param kappa_hat A \code{numeric} vector of jump sizes estimator.
#' @param r         An \code{integer} scalar of smoothness parameter of the underlying Holder space.
#' @param w         A \code{numeric} scalar in (0,1) representing the weight for interval truncation.
#' @param c_kappa   A \code{numeric} scalar to be multiplied by kappa estimator as the bandwidth in the local refinment.
#' @return  A vector of locally refined change points estimation.
#' @export
#' @author Haotian Xu
#' @examples
#' n = 150
#' v = c(floor(n/3), 2*floor(n/3)) # location of change points
#' r = 2
#' p = 6
#' Y = matrix(0, p, n) # matrix for data
#' mu0 = rep(0, p) # mean of the data
#' mu1 = rep(0, p)
#' mu1[1:floor(p/2)] = 2
#' Sigma0 = diag(p) #Covariance matrices of the data
#' Sigma1 = diag(p)
#' # Generate data
#' for(t in 1:n){
#'   if(t <= v[1] || t > v[2]){
#'      Y[,t] = MASS::mvrnorm(n = 1, mu0, Sigma0)
#'   }
#'   if(t > v[1] && t <= v[2]){
#'      Y[,t] = MASS::mvrnorm(n = 1, mu1, Sigma1)
#'   }
#' }## close for generate data
#' M = 100
#' intervals = WBS.intervals(M = M, lower = 1, upper = ncol(Y)) #Random intervals
#' h = 2*(1/n)^{1/(2*r+p)} # bandwith
#' temp = WBS.multi.nonpar.L2(Y, 1, ncol(Y), intervals$Alpha, intervals$Beta, h, delta = 15)
#' cpt_init = tuneBSmultinonpar(temp, Y)
#' kappa_hat = kappa.multi.nonpar.L2(cpt_init, Y, h_kappa = 0.01)
#' local.refine.multi.nonpar.L2(cpt_init, Y, kappa_hat, r = 2, w = 0.9, c_kappa = 2)
local.refine.multi.nonpar.L2.triwei = function(cpt_init, Y, kappa_hat, r = 2, w = 0.9, c_kappa = 10){
  p = dim(Y)[1]
  n = dim(Y)[2]
  cpt_init_long = c(0, cpt_init, n)
  cpt_init_numb = length(cpt_init)
  cpt_refined = rep(0, cpt_init_numb+1)
  h_refine = rep(NA, cpt_init_numb)
  minmax_mat = apply(Y, MARGIN = 1, function(x) c(min(x), max(x)))
  for (k in 1:cpt_init_numb){
    #aux1 = Y[,(cpt_init_long[k]+1):(cpt_init_long[k+1])]
    #aux2 = Y[,(cpt_init_long[k+1]+1):(cpt_init_long[k+2])]
    #integrateFCT = function(x){(kde.eval(t(aux1), eval.points = x, H = h_init*diag(p)) - kde.eval(t(aux2), eval.points = x, H = h_init*diag(p)))^2}
    #temp_integrate = cubature::cubintegrate(integrateFCT, lower = minmax_mat[1,], upper = minmax_mat[2,], method = "suave", maxEval = 1000)$integral
    #kappa_hat[k] = sqrt(temp_integrate)/sqrt((cpt_init_long[k+2]-cpt_init_long[k+1])*(cpt_init_long[k+1]-cpt_init_long[k])/(cpt_init_long[k+2]-cpt_init_long[k]))
    h_refine[k] = c_kappa*(kappa_hat[k])^(1/r)
    s = w*cpt_init_long[k] + (1-w)*cpt_init_long[k+1]
    e = (1-w)*cpt_init_long[k+1] + w*cpt_init_long[k+2]
    lower = ceiling(s) + 2
    upper = floor(e) - 2
    b = sapply(lower:upper, function(eta)(error.L2.multivariate.triwei(Y, ceiling(s), eta, h_refine[k]) + error.L2.multivariate.triwei(Y, (eta+1), floor(e), h_refine[k])))
    cpt_refined[k+1] = ceiling(s) + which.min(b)
  }
  return(list(cpt_refined = cpt_refined[-1]+1, kappa_hat = kappa_hat, h_refine = h_refine))
}


#' @export
jumpsize.multinonpar.L2.triwei = function(cpt_init, Y, h_kappa = 0.01, c_fac){
  p = dim(Y)[1]
  n = dim(Y)[2]
  cpt_init_long = c(0, cpt_init, n)
  cpt_init_numb = length(cpt_init)
  kappa_hat = rep(NA, cpt_init_numb)
  minmax_mat = apply(Y, MARGIN = 1, function(x) c(min(x), max(x)))
  for (k in 1:cpt_init_numb){
    aux1 = Y[,(cpt_init_long[k]+1):(cpt_init_long[k+1])]
    aux2 = Y[,(cpt_init_long[k+1]+1):(cpt_init_long[k+2])]
    integrateFCT = function(x){(kde.triwei.eval(t(aux1), bw = h_kappa, eval.points = x) - kde.triwei.eval(t(aux2), bw = h_kappa, eval.points = x))^2}
    temp_integrate = cubature::cubintegrate(integrateFCT, lower = minmax_mat[1,], upper = minmax_mat[2,], method = "suave", maxEval = 1000)$integral
    kappa_hat[k] = c_fac*sqrt(temp_integrate)
  }
  return(kappa_hat)
}