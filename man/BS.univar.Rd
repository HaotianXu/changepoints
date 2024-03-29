% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/univariate_mean.R
\name{BS.univar}
\alias{BS.univar}
\title{Standard binary segmentation for univariate mean change points detection.}
\usage{
BS.univar(y, s, e, delta = 2, level = 0)
}
\arguments{
\item{y}{A \code{numeric} vector of observations.}

\item{s}{A \code{integer} scalar of starting index.}

\item{e}{A \code{integer} scalar of ending index.}

\item{delta}{A positive \code{numeric} scalar of minimum spacing.}

\item{level}{Should be fixed as 0.}
}
\value{
An object of \code{\link[base]{class}} "BS", which is a \code{list} with the following structure:
 \item{S}{A vector of estimated change point locations (sorted in strictly increasing order).}
 \item{Dval}{A vector of values of CUSUM statistic.}
 \item{Level}{A vector representing the levels at which each change point is detected.}
 \item{Parent}{A matrix with the starting indices on the first row and the ending indices on the second row.}
}
\description{
Perform standard binary segmentation for univariate mean change points detection.
}
\examples{
set.seed(0)
cpt_true = c(20, 50, 170)
y = rnorm(300) + c(rep(0,20),rep(2,30),rep(0,120),rep(2,130))
temp = BS.univar(y, 1, length(y), delta = 5)
plot.ts(y)
points(x = tail(temp$S[order(temp$Dval)],4),
       y = y[tail(temp$S[order(temp$Dval)],4)], col = "red")
BS_result = thresholdBS(temp, tau = 4)
BS_result
print(BS_result$BS_tree, "value")
print(BS_result$BS_tree_trimmed, "value")
cpt_hat = sort(BS_result$cpt_hat[,1]) # the threshold tau is specified to be 4
Hausdorff.dist(cpt_hat, cpt_true)
cpt_LR = local.refine.univar(cpt_hat, y)
Hausdorff.dist(cpt_LR, cpt_true)
}
\references{
Wang, Yu and Rinaldo (2020) <doi:10.1214/20-EJS1710>.
}
\seealso{
\code{\link{thresholdBS}} for obtaining change points estimation, \code{\link{tuneBSunivar}} for a tuning version.
}
\author{
Haotian Xu
}
