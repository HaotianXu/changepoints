% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multivariate_nonparametric.R
\name{tuneBSmultinonpar}
\alias{tuneBSmultinonpar}
\title{A function to compute change points based on the multivariate nonparametic method with tuning parameter selected by FDR control.}
\usage{
tuneBSmultinonpar(BS_object, Y)
}
\arguments{
\item{BS_object}{A "BS" object produced by \code{WBS.multi.nonpar}.}

\item{Y}{A \code{numeric} matrix of observations with with horizontal axis being time, and vertical axis being dimension.}
}
\value{
A vector of estimated change points.
}
\description{
A function to compute change points based on the multivariate nonparametic method with tuning parameter selected by FDR control.
}
\examples{
n = 70
v = c(floor(n/3), 2*floor(n/3)) # location of change points
p = 4
Y = matrix(0, p, n) # matrix for data
mu0 = rep(0, p) # mean of the data
mu1 = rep(0, p)
mu1[1:floor(p/2)] = 2
Sigma0 = diag(p) #Covariance matrices of the data
Sigma1 = diag(p)*2
# Generate data
for(t in 1:n){
  if(t < v[1] || t > v[2]){
     Y[,t] = MASS::mvrnorm(n = 1, mu0, Sigma0)
  }
  if(t >= v[1] && t < v[2]){
     Y[,t] = MASS::mvrnorm(n = 1, mu1, Sigma1)
  }
}## close for generate data
M = 8
intervals = WBS.intervals(M = M, lower = 1, upper = ncol(Y)) #Random intervals
K_max = 30
h = 5*(K_max*log(n)/n)^{1/p} # bandwith
temp = WBS.multi.nonpar(Y, Y, 1, ncol(Y), intervals$Alpha, intervals$Beta, h, delta = 10)
S = tuneBSmultinonpar(temp, Y)
}
\references{
Padilla, Yu, Wang and Rinaldo (2019) <arxiv:1910.13289>.
}
\seealso{
\code{\link{WBS.multi.nonpar}}.
}
\author{
Oscar Hernan Madrid Padilla & Haotian Xu
}
