# clustered-multistate
R code for nonparametric analysis of clustered multi-state processes

## Description

This repository contains R functions for nonparametric analysis of (continuous-time) multi-state processes with cluster-correlated observations. These functions currently support nonparametric estimation of population-averaged transition probabilities, calculation of 95% pointwise confidence intervals and simultaneous confidence bands, and two-sample Kolmogorov-Smirnov-type tests. The functions do not impose assumptions regarding the within-cluster dependence and can be used for both Markov and non-Markov processes. Right censoring, left truncation, and informative cluster size are allowed.

## Main Functions

The main functions are `patp()` and `patp_test()`. Both functions are beta version.

### Dependencies
The functions require the R packages `mstate` and `Rfast` to be installed.

### Input data
The input data need to be a dataframe in the long format required by the `mstate` package. The dataframe should contain the variables

* `from`
* `to`
* `trans`
* `Tstart`
* `Tstop`
* `status`

For more details see [mstate](https://www.jstatsoft.org/article/view/v038i07). 

### Description of `patp()`

The function `patp()` calculates the working independence Aalen-Johansen estimator of the population averaged transition probabilities. These probabilities have the form Pr(X(t) = j| X(s)=h), where X(t) is the process of interest at time t, and h,j=1,...,k are possible states of the process X(t).


* `data`
* `tmat`
* `cid`
* `id`
* `h` 
* `j`
* `s`
* `weighted`
* `LMAJ`
* `B`
* `cband`


### Description of `patp_test()`

* `data`
* `tmat`
* `cid`
* `id`
* `group`
* `h` 
* `j`
* `s`
* `weighted`
* `LMAJ`
* `B`


## Examples

