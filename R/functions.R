#' @title Dynamic programming for l0 change points detection (RCPP)
#' @description TO DO
#' @param y         A \code{numeric} vector of observations.
#' @param gamma       A \code{numeric} scalar corresponding to the tuning parameter associated with the l0 penalty.
#' @param ...      Additional arguments.
#' @return TO DO.
#' @export
#' @author Haotian Xu
#' @examples
#' y = rnorm(300)+ c(rep(0,130),rep(-1,20),rep(1,20),rep(0,130))
#' D_P(y, 1)
D_P <- function(y, gamma, ...) {
  .Call('_changepoints_rcpp_D_P', PACKAGE = 'changepoints', y, gamma)
}

#' @title Dynamic programming for l0 change points detection
#' @description TO DO
#' @param gamma     A \code{numeric} scalar of the tuning parameter associated with the l0 penalty.
#' @param delta     A \code{numeric} scalar of minimum spacing.
#' @param y         A \code{numeric} vector of observations.
#' @param X         A \code{numeric} matrix of covariates. If missing, then the mean change of y is considered.
#' @param ...      Additional arguments.
#' @return TO DO.
#' @export
#' @author Haotian Xu
#' @examples
#' y = c(rep(0, 100), rep(1, 100))
#' DP(y, 1)
DP = function(gamma, delta, y, X = NULL, lambda = NULL, ...){
  N = length(y)
  bestvalue = rep(0,N+1)
  partition = rep(0,N)
  yhat = rep(NA, N)
  if(is.null(X) & is.null(lambda)){
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
      for(t in (l+1):r){
        yhat[t] = mean(y[(l+1):r])
      }
      r = l
      l = partition[r]
    }
    return(list(partition = partition, yhat = yhat))
  }else{
    if(is.null(X) + is.null(lambda) == 1){
      stop("In regression settings, X and lambda should both be specified")
    }
    if(length(dim(X)) != 2 | dim(X)[2] != N){
      stop("In regression settings, X should be a p-by-n design matrix")
    }
    p = dim(X)[1]
    bestvalue[1] = -gamma*log(max(N,p))
    for(r in 1:N){
      bestvalue[r+1] = Inf
      for(l in 1:r){
        b = bestvalue[l] + gamma*log(max(N,p)) + distanceR(l, r, y, X, lambda, delta)
        if(b < bestvalue[r+1]){
          bestvalue[r+1] = b
          partition[r] = l-1
        }
      }
    }
    r = N
    l = partition[r]
    while(r > 0){
      r = l
      l = partition[r]
    }
    return(list(partition = partition))
  }
}

#' @title Partition to localization
#' @description TO DO
#' @param x         A \code{numeric} vector of observations.
#' @param ...      Additional arguments.
#' @return TO DO.
#' @export
#' @author Haotian Xu
#' @examples
#' TO DO
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
  return(list(cpt = localization[-1], no.cpt = length(localization[-1])))
}


#' @title Simulate Model1 in [10] and follows Simulations 4.2
#' @description TO DO
#' @param d0         A \code{numeric} scalar of number of nonzero coefficients.
#' @param cpt.ture   A \code{integer} vector of true changepoints (sorted in strictly increasing order).
#' @param p          A \code{integer} scalar of dimensionality.
#' @param n          A \code{integer} scalar of sample size.
#' @param sigma      A \code{numeric} scalar of error standard deviation.
#' @param kappa      A \code{numeric} scalar of minimum jump size in terms of the l2 norm.
#' @param ...        Additional arguments.
#' @return A \code{list} with the structure:
#' \itemize{
#'  \item cpt.true    A vector of true changepoints (sorted in strictly increasing order).
#'  \item X           A p-by-n matrix of (random) design matrix.
#'  \item y           A n-dim vector of response variable.
#'  \item betafullmat A p-by-n matrix of coefficients.
#'  \item ...         Additional parameters.
#' }
#' @export
#' @author 
#' @examples
#' TO DO
data.generate.k.fluctuate2 = function(d0, cpt.true, p, n, sigma, kappa){
  #seed = 10
  if(d0 >= p){
    stop("d0 should be strictly smaller than p")
  }
  if(sigma <= 0){
    stop("sigma should be strictly larger than 0")
  }
  if(kappa <= 0){
    stop("kappa should be strictly larger than 0")
  }
  no.cpt = length(cpt.true)
  if(is.unsorted(cpt.true, strictly = TRUE) | min(cpt.true) <= 1 | max(cpt.true >= n) | no.cpt > n-2){
    stop("cpt.true is not correctly specified")
  }
  X = matrix(rnorm(p*n,0,1), p, n)
  y = matrix(0, n, 1)
  nonzero.element.loc = c(1:d0)
  cpt = c(0, cpt.true, n)
  beta = matrix(0, p, no.cpt+1)
  betafullmat = matrix(0, p, n)
  for (i in 1:(no.cpt+1)){
    if (i%%2 == 1){
      beta[nonzero.element.loc, i] = kappa/(2*sqrt(d0))
    }
    else{
      beta[nonzero.element.loc, i] = -kappa/(2*sqrt(d0))
    }
    y[(1+cpt[i]):cpt[i+1],] = rnorm(cpt[i+1] - cpt[i], t(X[,(1+cpt[i]):cpt[i+1]]) %*% beta[,i], sigma)
    for (j in (1+cpt[i]):cpt[i+1]){
      betafullmat[,j] = beta[,i] 
    }
  }
  List = list(cpt.true = cpt.true, X = X, y = y, betafullmat = betafullmat)
  return(List)
}


#' @title Prediction error in l2 norm for the lasso [10] 
#' @description TO DO
#' @param s         A \code{integer} scalar of starting index.
#' @param e         A \code{integer} scalar of ending index.
#' @param y         A \code{numeric} vector of response variable.
#' @param X         A \code{numeric} matrix of covariates.
#' @param lambda    A \code{numeric} scalar of lasso penalty.
#' @param delta     A \code{integer} scalar of minimum spacing.
#' @param ...        Additional arguments.
#' @return  A \code{numeric} scalar of prediction error in l2 norm.
#' @export
#' @author 
#' @examples
#' TO DO
distanceR = function(s, e, y, X, lambda, delta, ...){
  n = ncol(X)
  p = nrow(X)
  if (abs(s-e) > delta){
    fit = glmnet(x = t(X[,s:e]), y = y[s:e,], family = c("gaussian"),
                 alpha = 1, lambda = lambda*sqrt(max(log(max(n,p)), e-s))*sqrt(log(max(n,p)))/(e-s),intercept=F)
    coef_est = as.vector(fit$beta)
    yhat = t(X[,s:e])%*%coef_est
    d = norm(y[s:e] - yhat, type = "2")
  }
  else{
    d = Inf
  }
  return(d^2)
}