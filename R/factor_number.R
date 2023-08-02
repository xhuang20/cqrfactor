#' @title Factor number selection for composite quantile factor model
#'
#' @description This function computes the value of the two information criteria
#'  in Huang (2023) for factor number selection in composite quantile factor 
#'  model.  
#'
#' @param y A t by n matrix with n variables and t observations for each variable.
#' @param tau A vector of quantile positions, and \code{tau = (1:q)/(q+1)}, 
#'  where q is the numbrer of quantiles.
#' @param q The number of quantiles.
#' @param tol Tolerance for convergence of the MM algorithm. Defaults to 1e-6
#'  when the minconver option is 1. When the minconver option is 0, the 
#'  algorithm's convergence depends on the average parameter value change between
#'  iteration steps, and it can take long time to converge. We recommend to use 
#'  a bigger number such as 1e-3 when minconver = 1.
#' @param maxit The maximum iteration number between factor estimation and 
#'  factor loading estimation. Defaults to 500.
#' @param maxit_factor The maximum iteration number in the MM algorithm in 
#'  factor estimation. Defaults to 500.
#' @param maxit_loading The maximum iteration number in the MM algorithm in 
#'  factor loading estimation. Defaults to 500.
#' @param minconver If the value is 1, the algorithm converges when the minimum
#'  change among all estimated parameters between two iteration steps is less
#'  than or equal to tol. If the value is 0, the algorithm converges when the 
#'  average change of all estimated parameters is less than or equal to tol.
#'  This can take a long computation time. If \code{minconver = 0}, set tol to 
#'  a bigger number such as 1e-3. Defaults to 1. 
#' @param seed A non-negative integer for random seed used in the C++ function 
#'  \code{srand} to initialize factor values. If seed is -1, the initial value 
#'  of the factors will be set to the eigenvectors of the outer product of the 
#'  data matrix y. Defaults to 1.
#' @param standardize A logical variable to indicate whether to standardize the
#'  input data \code{y} during factor estimation. If \code{TRUE}, each column of
#'  \code{y} is standardized to have a mean of zeor and a standard deviation of
#'  one. Defaults to \code{FALSE}.
#' @param max_number The maximum number of factor to try in estimating the 
#'  number of factors. Defaults to \code{10}.
#'
#' @details This function takes the input data and outputs the value of the two
#'  information criteria for every factor number. The factor number that minimizes the 
#'
#' @return \code{factor_number} returns a list with the following components:
#'   \item{IC1}{A vector of values for information criterion 1 for every number
#'   of factors.} \item{IC2}{A vector of values for information criterion 1 for 
#'   every number of factors.} \item{IC1_factor_number}{The factor number that
#'   gives the smallest IC1 value.} \item{IC2_factor_number}{The factor number
#'   that gives the smallest IC2 value.}
#'
#' @export
#'
#' @usage factor_number(y,q = 5,tol = 1e-6,max_number = 10)
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
#' factor_number(y = dat, q = q, tol = 1e-6, maxit = 500, maxit_factor = 500, 
#'               maxit_loading = 500, minconver = 1, seed = 1, 
#'               standardize = FALSE, max_number = 10) 
#' }

factor_number <- function(y,
                          tau,
                          q = 5,
                          tol = 1e-6,
                          maxit = 500,
                          maxit_factor = 500,
                          maxit_loading = 500,
                          minconver = 1,
                          seed = 1,
                          standardize = FALSE,
                          max_number = 10){
  
  if (missing(tau)) {
    tau = (1:q) / (q + 1) 
  } else if (sum(tau <= 0 | tau >= 1) >= 1) {
    stop("The quantile positions in tau must be between 0 and 1.")
  } else {
    q = length(tau)
  }
  
  if (standardize) y = scale(y)
  
  tlen = dim(y)[1]
  nlen = dim(y)[2]
  ic   = matrix(0,max_number,2)
  result = cqrfactor(y = y, 
                     tau = tau,
                     q = q, 
                     r = max_number,     
                     tol = tol,
                     maxit = maxit,
                     maxit_factor = maxit_factor,
                     maxit_loading = maxit_loading,
                     minconver = minconver,
                     seed = seed)
  
  for (rnum in max_number:1) {
    fmat = result$fmat[,1:rnum]
    lmat = result$lmat[,1:rnum]
    quantiles = result$quantiles
    V = 0
    for (k in 1:q) {
      res = y - fmat %*% t(lmat) - quantiles[k]
      V = V + ifelse(res >= 0, tau[k] * res, (tau[k] - 1) * res)
    }
    V = sum(V) / (nlen * tlen)
    ic[rnum,1] = log(V) + rnum * (nlen + tlen) / (nlen*tlen) *  # IC 1
      log(nlen * tlen / (nlen + tlen))                    
    ic[rnum,2] = log(V) + rnum * (nlen + tlen) / (nlen*tlen) *  # IC 2
      log(log(nlen * tlen / (nlen + tlen)))                    
  }
  
  return(list(IC1 = ic[,1],
              IC2 = ic[,2],
              IC1_factor_number = which.min(ic[,1]),
              IC2_factor_number = which.min(ic[,2])))
}