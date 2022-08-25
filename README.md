
<!-- README.md is generated from README.Rmd. Please edit that file -->

# A collection of change-point localisation methods.

<!-- badges: start -->

[![R-CMD-check](https://github.com/HaotianXu/changepoints/workflows/R-CMD-check/badge.svg)](https://github.com/HaotianXu/changepoints/actions)
[![Last-changedate](https://img.shields.io/badge/last%20change-2022--08--25-green.svg)](https://github.com/HaotianXu/changepoints)
[![license](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)
<!-- badges: end -->

Performs a series of offline and/or online change-point localisation
algorithms for

1.  univariate mean
    -   [Wang, Yu and Rinaldo
        (2020)](https://doi.org/10.1214/20-EJS1710)
    -   [Yu, Padilla, Wang and Rinaldo
        (2020)](https://arxiv.org/abs/2006.03283)
2.  univariate polynomials
    -   [Yu, Chatterjee and Xu
        (2021)](https://doi.org/10.1214/21-EJS1963)
3.  univariate and multivariate nonparametric settings
    -   [Padilla, Yu, Wang and Rinaldo
        (2021)](https://doi.org/10.1214/21-EJS1809)
    -   [Padilla, Yu, Wang and Rinaldo
        (2021)](https://doi.org/10.1109/TIT.2021.3130330)
4.  high-dimensional covariances
    -   [Wang, Yu and Rinaldo
        (2021)](https://doi.org/10.3150/20-BEJ1249)
5.  high-dimensional networks with and without missing values
    -   [Wang, Yu and Rinaldo
        (2021)](https://doi.org/10.1214/20-AOS1953)
    -   [Yu, Padilla, Wang and Rinaldo
        (2021)](https://arxiv.org/abs/2101.05477)
    -   [Dubey, Xu and Yu (2021)](https://arxiv.org/abs/2110.06450)
6.  high-dimensional linear regression models
    -   [Rinaldo, Wang, Wen, Willett and Yu
        (2021)](https://proceedings.mlr.press/v130/rinaldo21a.html)
    -   [Xu, Wang, Zhao, and Yu
        (2022)](https://arxiv.org/abs/2207.12453)
7.  high-dimensional vector autoregressive models
    -   [Wang, Yu, Rinaldo and Willett
        (2019)](https://arxiv.org/abs/1909.06359)
8.  high-dimensional self exciting point processes
    -   [Wang, Yu and Willett (2020)](https://arxiv.org/abs/2006.03572)
9.  dependent dynamic nonparametric random dot product graphs
    -   [Padilla, Yu and Priebe
        (2019)](https://arxiv.org/abs/1911.07494)
10. robust univariate mean against adversarial attacks
    -   [Li and Yu
        (2021)](https://proceedings.neurips.cc/paper/2021/hash/c1e39d912d21c91dce811d6da9929ae8-Abstract.html)

## Installation

Users must have a (C++) compiler installed on their machine that is
compatible with R (e.g. Clang). The development version of
`changepoints` from [GitHub](https://github.com/) can be installed with:

``` r
## if not installed
## Install dependencies
install.packages(c("devtools","glmnet","gglasso","ks","data.tree"))

## install.packages("devtools")
devtools::install_github("HaotianXu/changepoints")
```

## Example

This is an example for offline univariate mean change point detection by
$l_0$ penalization:

``` r
library(changepoints)
## simulate data with true change points being 50, 100 and 150
set.seed(0)
y = c(rep(0, 50), rep(2, 50), rep(0, 50), rep(-2, 50)) + rnorm(200, mean = 0, sd = 1)
## estimate change points by l_0 penalization
gamma_set = c(0.01, 0.5, 1, 5, 10, 50) # possible value of tuning parameter
## perform cross-validation
DP_result = CV.search.DP.univar(y, gamma_set, delta = 5)
## estimate change points and perform local refinement
min_idx = which.min(DP_result$test_error)
cpt_DP_hat = unlist(DP_result$cpt_hat[[min_idx]])
cpt_DP_LR = local.refine.univar(cpt_DP_hat, y)
```

Alternatively, `wild binary segmentation` can also be performed:

``` r
## generate random intervals for WBS
intervals = WBS.intervals(M = 100, lower = 1, upper = 200)
## perform WBS
WBS_result = WBS.univar(y, 1, 200, intervals$Alpha, intervals$Beta, delta = 5)
WBS_result
## trim binary tree with threshold being 3
WBS_trimmed = thresholdBS(WBS_result, tau = 3)
## print the trimmed binary tree
print(WBS_trimmed$BS_tree_trimmed, "value")
## estimate change points and perform local refinement
cpt_WBS_hat = sort(WBS_trimmed$cpt_hat[,1])
cpt_BS_LR = local.refine.univar(cpt_WBS_hat, y)
```

`wild binary segmentation` with tuning parameter selected by information
criteria :

``` r
WBS_CPD_result = tuneBSunivar(WBS_result, y)
WBS_CPD_LR = local.refine.univar(WBS_CPD_result$cpt, y)
```
