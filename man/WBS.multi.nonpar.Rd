% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multivariate_nonparametric.R
\name{WBS.multi.nonpar}
\alias{WBS.multi.nonpar}
\title{Wild binary segmentation for multivariate nonparametric change points detection.}
\usage{
WBS.multi.nonpar(Y, W, s, e, Alpha, Beta, h, delta, level = 0)
}
\arguments{
\item{Y}{A \code{numeric} matrix of observations with with horizontal axis being time, and vertical axis being dimension.}

\item{W}{A copy of the matrix Y, it can be Y itself.}

\item{s}{A \code{integer} scalar of starting index.}

\item{e}{A \code{integer} scalar of ending index.}

\item{Alpha}{A \code{integer} vector of starting indices of random intervals.}

\item{Beta}{A \code{integer} vector of ending indices of random intervals.}

\item{h}{A \code{numeric} scalar of bandwidth parameter.}

\item{delta}{A \code{integer} scalar of minimum spacing.}

\item{level}{Should be fixed as 0.}
}
\value{
An object of \code{\link[base]{class}} "BS", which is a \code{list} with the following structure:
 \item{S}{A vector of estimated change points (sorted in strictly increasing order).}
 \item{Dval}{A vector of values of CUSUM statistic based on KS distance.}
 \item{Level}{A vector representing the levels at which each change point is detected.}
 \item{Parent}{A matrix with the starting indices on the first row and the ending indices on the second row.}
}
\description{
Perform wild binary segmentation for multivariate nonparametric change points detection.
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
M = 10
intervals = WBS.intervals(M = M, lower = 1, upper = ncol(Y)) #Random intervals
K_max = 30
h = 5*(K_max*log(n)/n)^{1/p} # bandwith
temp = WBS.multi.nonpar(Y, Y, 1, ncol(Y), intervals$Alpha, intervals$Beta, h, delta = 10)
result = thresholdBS(temp, median(temp$Dval))
}
\references{
Padilla, Yu, Wang and Rinaldo (2019) <arxiv:1910.13289>.
}
\seealso{
\code{\link{thresholdBS}} for obtain change points estimation, \code{\link{tuneBSmultinonpar}} for a tuning version.
}
\author{
Oscar Hernan Madrid Padilla & Haotian Xu
}
