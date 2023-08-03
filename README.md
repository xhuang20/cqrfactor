
# cqrfactor: An R package for composite quantile factor analysis

This R package computes the factors and factor loadings for the
composite quantile factor model (CQFM) in Huang (2023). We also provide
a function to implement the two information criteria for factor number
selection.

## Installation

You can install this package from
[GitHub](https://github.com/xhuang20/cqrfactor.git) using

``` r
devtools::install_github("xhuang20/cqrfactor")
```

## An example

Consider the following example of estimating CQFM.

``` r
library(cqrfactor)
tlen = 20                                 # number of observations
nlen = 30                                 # number of variables (cross-section units)

q  = 5                                    # number of quanile positions in CQFM
r  = 3                                    # number of factors in simulation
tau = (1:q) / (q + 1)                     # compute q quantile positions
tol = 1e-6                                # convergence tolerance 

set.seed(123)                             # random seed for factor initialization
lmat = matrix(rnorm(nlen*r), nlen, r)
fmat = matrix(rnorm(tlen*r), tlen, r)
u    = matrix(rnorm(tlen*nlen), tlen, nlen)
dat  = fmat %*% t(lmat) + u               # The dimension of the data is tlen*nlen.

result1 = cqrfactor(y = dat, q = q, r = r, tol = tol, maxit = 500, maxit_factor = 500, 
                   maxit_loading = 500, convergence = 0, seed = 1)

# Print the first few rows of the estimated factors.
print(head(result1$fmat))
```

    ##            [,1]       [,2]       [,3]
    ## [1,]  1.0893460 -1.4041511  0.4153459
    ## [2,]  0.2512260 -0.4595195 -1.4545558
    ## [3,]  0.5750143 -0.2423259  1.7023292
    ## [4,] -0.2628112  0.4494284  0.5707835
    ## [5,] -1.3262489 -1.2209089 -0.1795293
    ## [6,]  0.4016536  0.4946898 -1.7261063

``` r
# Print the first few rows of the esimated factor loadings.
print(head(result1$lmat))
```

    ##            [,1]        [,2]       [,3]
    ## [1,]  0.5315626  0.41512813 -0.6055471
    ## [2,] -0.2711427  0.07450419  0.4065602
    ## [3,] -0.4411295 -1.49518582 -0.6655097
    ## [4,] -1.4279021 -0.11505804 -0.3397502
    ## [5,] -1.6456471 -0.55344470  0.2225508
    ## [6,] -0.5503292 -1.97913038 -0.7027694

Next, consider the example of factor number estimation.

``` r
result2 = factor_number(y = dat, q = q, tol = tol, maxit = 500, maxit_factor = 500, 
                        maxit_loading = 500, convergence = 0, seed = 1, 
                        max_number = 10) 

# Print the estimated number of factors.
print(result2$IC1_factor_number)
```

    ## [1] 2

## Reference

Huang, X. (2023) Composite quantile factor models, working paper.
