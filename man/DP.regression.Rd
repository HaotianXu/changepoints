% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HD_regression.R
\name{DP.regression}
\alias{DP.regression}
\title{Dynamic programming algorithm for regression change points localisation with \eqn{l_0} penalisation.}
\usage{
DP.regression(y, X, gamma, lambda, delta, eps = 0.001)
}
\arguments{
\item{y}{A \code{numeric} vector of response variable.}

\item{X}{A \code{numeric} matrix of covariates with vertical axis being time.}

\item{gamma}{A positive \code{numeric} scalar stands for tuning parameter associated with \eqn{l_0} penalty.}

\item{lambda}{A positive \code{numeric} scalar stands for tuning parameter associated with the lasso penalty.}

\item{delta}{A positive \code{integer} scalar stands for minimum spacing.}

\item{eps}{A \code{numeric} scalar of precision level for convergence of lasso.}
}
\value{
An object of \code{\link[base]{class}} "DP", which is a \code{list} with the following structure:
 \item{partition}{A vector of the best partition.}
 \item{cpt}{A vector of change points estimation.}
}
\description{
Perform dynamic programming algorithm for regression change points localisation.
}
\examples{
d0 = 10
p = 20
n = 100
cpt_true = c(30, 70)
data = simu.change.regression(d0, cpt_true, p, n, sigma = 1, kappa = 9)
temp = DP.regression(y = data$y, X = data$X, gamma = 2, lambda = 1, delta = 5)
cpt_hat = temp$cpt
}
\references{
Rinaldo, Wang, Wen, Willett and Yu (2020) <arxiv:2010.10410>
}
\author{
Daren Wang & Haotian Xu
}
