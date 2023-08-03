
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
                   maxit_loading = 100, convergence = 0, seed = 1)

# Print the first few rows of the estimated factors.
print(head(result1$fmat))
```

    ##            [,1]       [,2]       [,3]
    ## [1,]  1.2440284 -1.3870907  0.3169970
    ## [2,]  0.2650284 -0.3000111 -1.6193824
    ## [3,]  0.5060248 -0.1665415  1.7133201
    ## [4,] -0.2872699  0.4354499  0.6132991
    ## [5,] -1.2986753 -1.2932778 -0.1385656
    ## [6,]  0.4276419  0.4375698 -1.5541025

``` r
# Print the first few rows of the esimated factor loadings.
print(head(result1$lmat))
```

    ##            [,1]        [,2]        [,3]
    ## [1,]  0.4494047  0.47782197 -0.56562040
    ## [2,] -0.2857761  0.09051411  0.41609015
    ## [3,] -0.3115365 -1.49887024 -0.61275507
    ## [4,] -1.4140430 -0.28577138 -0.25125264
    ## [5,] -1.6580694 -0.62859565  0.03986152
    ## [6,] -0.4189660 -2.04412429 -0.63587697

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
