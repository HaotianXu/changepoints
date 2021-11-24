
<!-- README.md is generated from README.Rmd. Please edit that file -->

# changepoints

<!-- badges: start -->

[![R-CMD-check](https://github.com/HaotianXu/changepoints/workflows/R-CMD-check/badge.svg)](https://github.com/HaotianXu/changepoints/actions)
<!-- badges: end -->

The goal of changepoints is to …

## Installation

You can install the development version of changepoints from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("HaotianXu/changepoints")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(changepoints)
## basic example code
delta = 5
sigma2 = 1

set.seed(0)
y = c(rep(0, 50), rep(2, 50), rep(0, 50), rep(-2, 50)) + rnorm(200, mean = 0, sd = sqrt(sigma2))
n = length(y)
gamma_set = c(0.01, 0.5, 1, 5, 10, 50)
DP_result = CV.search.DP.univar(y, gamma_set, delta = 5)
#>             [,1]      [,2]      [,3]      [,4]      [,5]      [,6]    
#> cpt_hat     numeric,7 numeric,6 numeric,4 numeric,3 numeric,3 101     
#> K_hat       7         6         4         3         3         1       
#> test_error  89.51068  89.50708  106.2182  85.13471  85.13471  175.6256
#> train_error 83.85808  83.87686  85.52212  87.96934  87.96934  161.2432
min_idx = which.min(DP_result$test_error)
cpt_DP_hat = unlist(DP_result$cpt_hat[[min_idx]])
cpt_DP_LR = local.refine.univar(cpt_DP_hat, y)
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/v1/examples>.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
