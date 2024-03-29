% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/univariate_polynomial.R
\name{CV.search.DP.poly}
\alias{CV.search.DP.poly}
\title{Grid search for dynamic programming to select the tuning parameter through Cross-Validation.}
\usage{
CV.search.DP.poly(y, r, gamma_set, delta)
}
\arguments{
\item{y}{A \code{numeric} vector of observations.}

\item{r}{An \code{integer} scalar order of polynomials.}

\item{gamma_set}{A \code{numeric} vector of candidate tuning parameter associated with the l0 penalty.}

\item{delta}{A positive \code{integer} scalar of minimum spacing.}
}
\value{
A \code{list} with the following structure:
 \item{cpt_hat}{A list of vector of estimated change points locations (sorted in strictly increasing order)}
 \item{K_hat}{A list of scalar of number of estimated change points}
 \item{test_error}{A list of vector of testing errors}
 \item{train_error}{A list of vector of training errors}
}
\description{
Perform grid search for dynamic programming to select the tuning parameter through Cross-Validation.
}
\examples{
set.seed(0)
cpt_true = c(20, 50, 170)
y = rnorm(300) + c(rep(0,20),rep(2,30),rep(0,120),rep(2,130))
plot.ts(y)
gamma_set = 3:9
DP_result = CV.search.DP.poly(y, r = 2, gamma_set, delta = 5)
min_idx = which.min(DP_result$test_error)
cpt_init = unlist(DP_result$cpt_hat[min_idx])
local.refine.poly(cpt_init, y, r = 2, delta_lr = 5)
}
\references{
Yu and Chatterjee (2020) <arXiv:2007.09910>
}
\author{
Haotian Xu
}
