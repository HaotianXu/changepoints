% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multivariate_nonparametric_L2.R
\name{local.refine.multi.nonpar.L2.biwei}
\alias{local.refine.multi.nonpar.L2.biwei}
\title{Local refinement for multivariate nonparametric change points localisation based on L2 distance.}
\usage{
local.refine.multi.nonpar.L2.biwei(
  cpt_init,
  Y,
  kappa_hat,
  r = 2,
  w = 0.9,
  c_kappa = 10
)
}
\arguments{
\item{cpt_init}{An \code{integer} vector of initial change points estimation (sorted in strictly increasing order).}

\item{Y}{A \code{numeric} matrix of observations with with horizontal axis being time, and vertical axis being dimension.}

\item{kappa_hat}{A \code{numeric} vector of jump sizes estimator.}

\item{r}{An \code{integer} scalar of smoothness parameter of the underlying Holder space.}

\item{w}{A \code{numeric} scalar in (0,1) representing the weight for interval truncation.}

\item{c_kappa}{A \code{numeric} scalar to be multiplied by kappa estimator as the bandwidth in the local refinment.}
}
\value{
A vector of locally refined change points estimation.
}
\description{
Perform local refinement for multivariate nonparametric change points localisation based on L2 distance.
}
\examples{
n = 150
v = c(floor(n/3), 2*floor(n/3)) # location of change points
r = 2
p = 6
Y = matrix(0, p, n) # matrix for data
mu0 = rep(0, p) # mean of the data
mu1 = rep(0, p)
mu1[1:floor(p/2)] = 2
Sigma0 = diag(p) #Covariance matrices of the data
Sigma1 = diag(p)
# Generate data
for(t in 1:n){
  if(t <= v[1] || t > v[2]){
     Y[,t] = MASS::mvrnorm(n = 1, mu0, Sigma0)
  }
  if(t > v[1] && t <= v[2]){
     Y[,t] = MASS::mvrnorm(n = 1, mu1, Sigma1)
  }
}## close for generate data
M = 100
intervals = WBS.intervals(M = M, lower = 1, upper = ncol(Y)) #Random intervals
h = 2*(1/n)^{1/(2*r+p)} # bandwith
temp = WBS.multi.nonpar.L2(Y, 1, ncol(Y), intervals$Alpha, intervals$Beta, h, delta = 15)
cpt_init = tuneBSmultinonpar(temp, Y)
kappa_hat = kappa.multi.nonpar.L2(cpt_init, Y, h_kappa = 0.01)
local.refine.multi.nonpar.L2(cpt_init, Y, kappa_hat, r = 2, w = 0.9, c_kappa = 2)
}
\author{
Haotian Xu
}
