% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RDPG.R
\name{simu.RDPG}
\alias{simu.RDPG}
\title{Simulate a dot product graph (without change point).}
\usage{
simu.RDPG(x_mat, n, symm = TRUE, self = FALSE)
}
\arguments{
\item{x_mat}{A \code{numeric} matrix representing the latent positions with horizontal axis being latent dimensions and vertical axis being nodes (each entry takes value in \eqn{[0,1]}).}

\item{n}{A \code{integer} scalar representing the number of observations.}

\item{symm}{A \code{logic} scalar indicating if adjacency matrices are required to be symmetric.}

\item{self}{A \code{logic} scalar indicating if adjacency matrices are required to have self-loop.}
}
\value{
A \code{list} with the following structure:
 \item{obs_mat}{A matrix, with each column be the vectorized adjacency (sub)matrix. For example, if "symm = TRUE" and "self = FALSE", only the strictly lower triangular matrix is considered.}
 \item{graphon_mat}{Underlying graphon matrix.}
}
\description{
Simulate a dot product graph (without change point). The generated data is a matrix with each column corresponding to the vectorized adjacency (sub)matrix at a time point. For example, if the network matrix is required to be symmetric and without self-loop, only the strictly lower diagonal entries are considered.
}
\examples{
p = 20 # number of nodes
n = 50 # sample size for each segment
lat_dim_num = 5 # number of latent dimensions
set.seed(1)
x_mat = matrix(runif(p*lat_dim_num), nrow = p, ncol = lat_dim_num)
x_tilde_mat = matrix(runif(p*lat_dim_num), nrow = p, ncol = lat_dim_num)
y_mat = rbind(x_tilde_mat[1:floor(p/4),], x_mat[(floor(p/4)+1):p,])
rdpg1 = simu.RDPG(x_mat, n, symm = TRUE, self = FALSE)
rdpg2 = simu.RDPG(y_mat, n, symm = TRUE, self = FALSE)
data1_mat = rdpg1$obs_mat
data2_mat = rdpg2$obs_mat
data_mat = cbind(data1_mat, data2_mat)
}
\author{
Haotian Xu
}
