% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HD_covariance.R
\name{WBSIP.cov}
\alias{WBSIP.cov}
\title{Wild binary segmentation for covariance change points detection through Independent Projection.}
\usage{
WBSIP.cov(X, X_prime, s, e, Alpha, Beta, delta, level = 0)
}
\arguments{
\item{X}{A \code{numeric} vector of observations.}

\item{X_prime}{A \code{numeric} vector of observations which are independent copy of X.}

\item{s}{A \code{integer} scalar of starting index.}

\item{e}{A \code{integer} scalar of ending index.}

\item{Alpha}{A \code{integer} vector of starting indices of random intervals.}

\item{Beta}{A \code{integer} vector of ending indices of random intervals.}

\item{delta}{A positive \code{integer} scalar of minimum spacing.}

\item{level}{A parameter for tracking the level at which a change point is detected. Should be fixed as 0.}
}
\value{
An object of \code{\link[base]{class}} "BS", which is a \code{list} with the following structure:
 \item{S}{A vector of estimated change points (sorted in strictly increasing order)}
 \item{Dval}{A vector of values of CUSUM statistic based on KS distance}
 \item{Level}{A vector representing the levels at which each change point is detected}
 \item{Parent}{A matrix with the starting indices on the first row and the ending indices on the second row}
}
\description{
Perform wild binary segmentation for covariance change points detection through Independent Projection
}
\examples{
p = 10
A1 = gen.cov.mat(p, 1, "equal")
A2 = gen.cov.mat(p, 3, "power")
A3 = A1
set.seed(1234)
X = cbind(t(MASS::mvrnorm(50, mu = rep(0, p), A1)), 
          t(MASS::mvrnorm(50, mu = rep(0, p), A2)), 
          t(MASS::mvrnorm(50, mu = rep(0, p), A3)))
X_prime = cbind(t(MASS::mvrnorm(50, mu = rep(0, p), A1)), 
                t(MASS::mvrnorm(50, mu = rep(0, p), A2)), 
                t(MASS::mvrnorm(50, mu = rep(0, p), A3)))
intervals = WBS.intervals(M = 120, lower = 1, upper = dim(X)[2])
temp = WBSIP.cov(X, X_prime, 1, dim(X)[2], intervals$Alpha, intervals$Beta, delta = 5)
tau = sqrt(p*log(ncol(X)))*1.5
sort(thresholdBS(temp, tau)$cpt_hat[,1])
}
\references{
Wang, Yu and Rinaldo (2021) <doi:10.3150/20-BEJ1249>.
}
\seealso{
\code{\link{thresholdBS}} for obtain change points estimation.
}
\author{
Haotian Xu
}
