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

* `from`:
* `to`:
* `trans`:
* `Tstart`:
* `Tstop`:
* `status`:

For more details see <https://www.jstatsoft.org/article/view/v038i07>. 

### Function `patp()`

The function `patp()` calculates the working independence Aalen-Johansen estimator of the population averaged transition probabilities. These probabilities have the form Pr(X(t) = j| X(s) = h), where X(t) is the process of interest at time t, and h,j=1,...,k are possible states of the process X(t). The function has the following arguments:


* `data`: a data.frame in the long format required by the `mstate` package
* `tmat`: a matrix of possible transitions between states of the process where different transitions are identified by a different integer. If a direct transition between two states is not possible it is indicated as NA. This matrix can be obtained via the `mstate` function `transMat()`. For the special case of the illness-death model without recovery this can be achieved as follows
```
> tmat <- transMat(x = list(c(2, 3), c(3), c()), 
+                  names = c("Health", "Illness", "Death"))
```
The resulting matrix `tmat` is
```
> tmat
         to
from      Health Illness Death
  Health      NA       1     2
  Illness     NA      NA     3
  Death       NA      NA    NA
```

* `cid`: variable name that identifies the clusters.
* `id`: variable name that identifies the individual observations.
* `h`:
* `j`:
* `s`:
* `weighted`:
* `LMAJ`:
* `B`:
* `cband`:


### Function `patp_test()`

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

