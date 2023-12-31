#' @title Composite quantile factor model
#'
#' @description This function computes the composite quantile-based factors and 
#'  factor loadings for a given number of factors.
#'
#' @param y A t by n matrix with n variables and t observations for each variable.
#' @param tau A vector of quantile positions, and \code{tau = (1:q)/(q+1)}, 
#'  where q is the numbrer of quantiles.
#' @param q The number of quantiles.
#' @param r The number of factors.
#' @param tol Tolerance for convergence of the MM algorithm. Defaults to 1e-6
#'  when the minconver option is 1. When the minconver option is 0, the 
#'  algorithm's convergence depends on the average parameter value change between
#'  iteration steps, and it can take long time to converge. We recommend to use 
#'  a bigger number such as 1e-3 when minconver = 1.
#' @param maxit The maximum iteration number between factor estimation and factor 
#'  loading estimation. Defaults to 500.
#' @param maxit_factor The maximum iteration number in the MM algorithm in factor
#'  estimation. Defaults to 500.
#' @param maxit_loading The maximum iteration number in the MM algorithm in factor
#'  loading estimation. Defaults to 500.
#' @param convergence If the value is 0, the algorithm converges when the minimum
#'  change among all estimated parameters between two iteration steps is less
#'  than or equal to tol. If the value is 1, the algorithm converges when the 
#'  average change of all estimated parameters is less than or equal to tol.
#'  This can take a long computation time. If minconver = 1, set tol to a bigger
#'  number such as 1e-3. Defaults to 0. 
#' @param seed A non-negative integer for random seed used in the C++ function 
#'  \code{srand} to initialize factor values. If seed is -1, the initial value 
#'  of the factors will be set to the eigenvectors of the outer product of the 
#'  data matrix y. Defaults to 1.
#' @param standardize A logical variable to indicate whether to standardize the
#'  input data \code{y} during factor estimation. If \code{TRUE}, each column of
#'  \code{y} is standardized to have a mean of zeor and a standard deviation of
#'  one. Defaults to \code{FALSE}.
#'
#' @details This function takes the input data and outputs the estimated factors
#'  and factor loadings. The input data is organized such that rows represent
#'  observations and columns represent variables. No missing data is allowed and
#'  the input panel data must have a balanced structure. For fast computation,
#'  we recommend to use the option \code{minconver = 1}, and one can adjust the
#'  value of the convergence tolerance option \code{tol} to check the stability
#'  of the solution. If \code{minconver = 0}, the program takes long time to 
#'  finish, and it is recommended to set \code{tol} to a relative large number
#'  such as \code{1e-3} to speed up the convergence and run the code on a 
#'  cluster if possible.
#'
#' @return \code{mmfactor} returns a list with the following components:
#'  \item{fmat}{The estimated factors, a t by r matrix, where r is the number
#'  of factors.} \item{lmat}{The estimated factor loading matrix, an n by r 
#'  matrix.} \item{quantiles}{The q estimated quantiles of the error term.} 
#'  \item{iteration_mat}{A \code{maxit} by 3 matrix. Each row corresponds to an
#'  iteration between factor estimation and factor loading estimation. The 
#'  second element is the number of iterations the MM uses for estimating the 
#'  factors while the third element is the number of iteration for the factor
#'  loadings.}
#'
#' @export
#'
#' @usage cqrfactor(y,q = 5,tol = 1e-5, maxit = 1000)
#'
#' @examples
#' \dontrun{
#' library(cqrfactor)
#' tlen = 20
#' nlen = 30
#' r  = 3
#' q  = 5
#' tau = (1:q) / (q + 1) 
#' 
#' set.seed(123)
#' lmat = matrix(rnorm(nlen*r), nlen, r)
#' fmat = matrix(rnorm(tlen*r), tlen, r)
#' u    = matrix(rnorm(tlen*nlen), tlen, nlen)
#' dat  = fmat %*% t(lmat) + u 
#' 
#' cqrfactor(y = dat, q = q, r = r, tol = 1e-6, maxit = 100, maxit_factor = 100, 
#'           maxit_loading = 100, seed = 1) 
#' }

cqrfactor <- function(y,
                      tau,
                      q = 5,
                      r = 3,
                      tol = 1e-6,
                      maxit = 500,
                      maxit_factor = 500,
                      maxit_loading = 500,
                      convergence = 0,
                      seed = 1,
                      standardize = FALSE){
  
  if (missing(tau)) {
    tau = (1:q) / (q + 1) 
  } else if (sum(tau <= 0 | tau >= 1) >= 1) {
    stop("The quantile positions in tau must be between 0 and 1.")
  } else {
    q = length(tau)
  }
  
  if (standardize) y = scale(y)
  
  result = mmfactor(y = y,
                    tau = tau,
                    q = q,
                    r = r,
                    tol = tol,
                    maxit = maxit,
                    maxit_factor = maxit_factor,
                    maxit_loading = maxit_loading,
                    convergence = convergence,
                    seed = seed)
  
  return(result)
  
}