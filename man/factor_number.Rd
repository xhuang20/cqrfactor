% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/factor_number.R
\name{factor_number}
\alias{factor_number}
\title{Factor number selection for composite quantile factor model}
\usage{
factor_number(y,q = 5,tol = 1e-6,max_number = 10)
}
\arguments{
\item{y}{A t by n matrix with n variables and t observations for each variable.}

\item{tau}{A vector of quantile positions, and \code{tau = (1:q)/(q+1)},
where q is the numbrer of quantiles.}

\item{q}{The number of quantiles.}

\item{tol}{Tolerance for convergence of the MM algorithm. Defaults to 1e-6
when the minconver option is 1. When the minconver option is 0, the
algorithm's convergence depends on the average parameter value change between
iteration steps, and it can take long time to converge. We recommend to use
a bigger number such as 1e-3 when minconver = 1.}

\item{maxit}{The maximum iteration number between factor estimation and
factor loading estimation. Defaults to 500.}

\item{maxit_factor}{The maximum iteration number in the MM algorithm in
factor estimation. Defaults to 500.}

\item{maxit_loading}{The maximum iteration number in the MM algorithm in
factor loading estimation. Defaults to 500.}

\item{convergence}{If the value is 0, the algorithm converges when the minimum
change among all estimated parameters between two iteration steps is less
than or equal to tol. If the value is 1, the algorithm converges when the
average change of all estimated parameters is less than or equal to tol.
This can take a long computation time. If minconver = 1, set tol to a bigger
number such as 1e-3. Defaults to 0.}

\item{seed}{A non-negative integer for random seed used in the C++ function
\code{srand} to initialize factor values. If seed is -1, the initial value
of the factors will be set to the eigenvectors of the outer product of the
data matrix y. Defaults to 1.}

\item{standardize}{A logical variable to indicate whether to standardize the
input data \code{y} during factor estimation. If \code{TRUE}, each column of
\code{y} is standardized to have a mean of zeor and a standard deviation of
one. Defaults to \code{FALSE}.}

\item{max_number}{The maximum number of factor to try in estimating the
number of factors. Defaults to \code{10}.}
}
\value{
\code{factor_number} returns a list with the following components:
\item{IC1}{A vector of values for information criterion 1 for every number
of factors.} \item{IC2}{A vector of values for information criterion 1 for
every number of factors.} \item{IC1_factor_number}{The factor number that
gives the smallest IC1 value.} \item{IC2_factor_number}{The factor number
that gives the smallest IC2 value.}
}
\description{
This function computes the value of the two information criteria
in Huang (2023) for factor number selection in composite quantile factor
model.
}
\details{
This function takes the input data and outputs the value of the two
information criteria for every factor number. The factor number that minimizes the
}
\examples{
\dontrun{
library(cqrfactor)
tlen = 20
nlen = 30
r  = 3
q  = 5
tau = (1:q) / (q + 1) 

set.seed(123)
lmat = matrix(rnorm(nlen*r), nlen, r)
fmat = matrix(rnorm(tlen*r), tlen, r)
u    = matrix(rnorm(tlen*nlen), tlen, nlen)
dat  = fmat \%*\% t(lmat) + u 

factor_number(y = dat, q = q, tol = 1e-6, maxit = 500, maxit_factor = 500, 
              maxit_loading = 500, convergence = 1, seed = 1, 
              standardize = FALSE, max_number = 10) 
}
}
