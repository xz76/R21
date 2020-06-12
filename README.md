# clustered-multistate
R code for nonparametric analysis of clustered multi-state processes

## Description

This repository contains R functions for nonparametric analysis of (continuous-time) multi-state processes with cluster-correlated observations. These functions currently support nonparametric estimation of population-averaged transition probabilities, calculation of 95% pointwise confidence intervals and simultaneous confidence bands, and two-sample Kolmogorov-Smirnov-type tests. The functions do not impose assumptions regarding the within-cluster dependence and can be used for both Markov and non-Markov processes. Right censoring, left truncation, and association between cluster size and the multi-state process (informative cluster size) are allowed.

## Main Functions

The main functions are `patp()` and `patp_test()`. Both functions are beta version.

### Dependencies
The functions require the R packages `mstate` and `Rfast` to be installed and loaded.

### Input data
The input data need to be a data frame in the long format required by the `mstate` package. The data frame should contain the variables

* `Tstart`: starting time of the interval in the record.
* `Tstop`: ending time of the interval in record.
* `from`: the state of the process at `Tstart`. The possible values are 1,...,k. 
* `to`: the state of the process at `Tstop`. The possible values are 1,...,k.
* `trans`: an integer that uniquely identifies the transition.
* `status`: indicator variable. If `status=1`, the correspoding transition has been observed.

The `mstate` function `msprep()` can be used to reshape a dataset in wide format into the required long format. For more details see <https://www.jstatsoft.org/article/view/v038i07>. 


### Function `patp()`

The function `patp()` calculates the working independence Aalen-Johansen estimator of the population-averaged transition probabilities. These probabilities have the form Pr(X(t) = j| X(s) = h), where X(t) is the process of interest at time t, and h,j=1,...,k are possible states of the process X(t). The function has the following arguments:

* `data`: a data.frame in the long format required by the `mstate` package.
* `tmat`: a matrix of possible transitions between states of the process where different transitions are identified by a different integer. If a direct transition between two states is not possible it is indicated as NA. This matrix can be obtained via the `mstate` function `transMat()`.
* `cid`: variable name that identifies the clusters.
* `id`: variable name that identifies the individual observations.
* `h`: the state h in Pr(X(t) = j| X(s) = h).
* `j`: the state j in Pr(X(t) = j| X(s) = h).
* `s`: the time s in Pr(X(t) = j| X(s) = h). The default value is `0`.
* `weighted`: logical value. If `TRUE`, the estimator is weighted by the inverse of the cluster sizes. This is useful when cluster size is random and expected to be informative. The defaul value is `FALSE`.
* `LMAJ`: logical value. If `TRUE`, the landmark version of the estimator is returned. This is useful when `s>0` and the Markov assumption is not plausible. The defaul value is `FALSE`.
* `B`: number of nonparametric cluster bootstrap replications. If `B=0`, no standard errors or confidence intervals/bands are returned. The default value is `100`.
* `cband`: logical value. If `TRUE`, the limits of the 95% simultaneous confidence band are returned. The defaul value is `FALSE`.


### Function `patp_test()`

The function `patp_test()` calculates the p-value for the comparison of the population-averaged transition probability Pr(X(t) = j| X(s) = h) between two groups, using a two-sample Kolmogorov-Smirnov-type test. The function performes has following arguments:

* `data`: a data.frame in the long format required by the `mstate` package.
* `tmat`: a matrix of possible transitions between states of the process where different transitions are identified by a different integer. If a direct transition between two states is not possible it is indicated as NA. This matrix can be obtained via the `mstate` function `transMat()`.
* `cid`: variable name that identifies the clusters.
* `id`: variable name that identifies the individual observations.
* `group`: variable name of the binary grouping variable.
* `h`: the state h in Pr(X(t) = j| X(s) = h).
* `j`: the state j in Pr(X(t) = j| X(s) = h).
* `s`: the time s in Pr(X(t) = j| X(s) = h). The default value is `0`.
* `weighted`: logical value. If `TRUE`, the estimators are weighted by the inverse of cluster sizes. This is useful when cluster size is random and expected to be informative. The defaul value is `FALSE`.
* `LMAJ`: logical value. If `TRUE`, the landmark version of the estimator is used in the test. This is useful when `s>0` and the Markov assumption is not plausible. The defaul value is `FALSE`.
* `B`: number of nonparametric cluster bootstrap replications. The default value is `1000`.


## Example

The artificial dataset `example_data.csv` (included in this repository) contains clustered observations from an illness-death process without recovery. The matrix `tmat` of possible transitions for this process can be created as follows
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
The dataset can be obtained as follows
```
> library(foreign)
> data <- read.csv("example_data.csv")
> head(data)
  id cid       ill ill.s       dth dth.s group
1  1   1 1.9184301     0 1.9184301     0     1
2  2   1 1.9391350     0 1.9391350     0     1
3  3   1 2.6312586     0 2.6312586     0     1
4  4   1 0.3779283     0 0.3779283     0     1
5  5   1 2.1919740     0 2.1919740     1     1
6  6   1 1.5983068     1 2.7530976     0     1
```
The variables `cid` and `id` correspond to the cluster identification number and the individual identification number, respectively. The variable `group` is the binary group indicator, `ill` is the time of arrival at the illness state, and `ill.s` is the indicator of illness. `dth` and `dth.s` are the death time and death indicator, respectively. The data frame contains one record per individual. Under the illness-death model without recovery, there are four possible scenarios; i) no "illness" or "death" are observed (i.e. right censoring while in the "healhty" state), ii) only "illness"
 is observed (i.e. right censoring while in the "illness" state), iii) only "death" is observed (i.e. there was a direct transition from the "healhty" state to the "death" state), and iv) both "illness" and "death" are observed. An example of these four cases from the data frame `data` is presented below:
```
> data[data$id %in% c(1,5,6,16),]
   id cid       ill ill.s       dth dth.s group
1   1   1 1.9184301     0 1.9184301     0     1
5   5   1 2.1919740     0 2.1919740     1     1
6   6   1 1.5983068     1 2.7530976     0     1
16 16   1 0.3580249     1 0.7174463     1     0
```
The data frome can be reshaped in the approriate long format using the `mstate` function `msprep()` as follows
```
> data <- msprep(data = data, trans = tmat, time = c(NA, "ill", "dth"),
+                status = c(NA, "ill.s", "dth.s"),
+                keep = c("cid", "group"))
```
The data of the four records listed above are now as follows:
```
> data[data$id %in% c(1,5,6,16),]
An object of class 'msdata'

Data:
   id from to trans    Tstart     Tstop      time status cid group
1   1    1  2     1 0.0000000 1.9184301 1.9184301      0   1     1
2   1    1  3     2 0.0000000 1.9184301 1.9184301      0   1     1
9   5    1  2     1 0.0000000 2.1919740 2.1919740      0   1     1
10  5    1  3     2 0.0000000 2.1919740 2.1919740      1   1     1
11  6    1  2     1 0.0000000 1.5983068 1.5983068      1   1     1
12  6    1  3     2 0.0000000 1.5983068 1.5983068      0   1     1
13  6    2  3     3 1.5983068 2.7530976 1.1547908      0   1     1
33 16    1  2     1 0.0000000 0.3580249 0.3580249      1   1     0
34 16    1  3     2 0.0000000 0.3580249 0.3580249      0   1     0
35 16    2  3     3 0.3580249 0.7174463 0.3594214      1   1     0
```
Individuals who did not visit the "illness" state (individuals with `id` 1 and 5 in the example) have two records only. This is because such individuals do not contain information about the transition 3 (i.e. from the "illness" state to the "death" state). Individuals who visited the "illness" state (individuals with `id` 6 and 16) have three records since they contain information about all three states.

Estimating the population-averaged transition probability P(X(t) = 2| X(0) = 1) and calculating standard errors and 95% pointwise confidence intervals based on 100 cluster bootstrap replications can be achieved as follows
```
> set.seed(1234)
> P12 <- patp(data=data, tmat=tmat, cid="cid", id="id", 
+             h=1, j=2, s=0, B=100)
```
For the illness-death model without recovery, the transition probability P(X(t) = 2| X(0) = 1) is equal to the state occupation probability P(X(t) = 2). To also calculate 95% simultaneous confidence bands requires the code
```
> P12 <- patp(data=data, tmat=tmat, cid="cid", id="id", 
+             h=1, j=2, s=0, B=1000, cband=TRUE)
```
It is recommended to use at least 1000 cluster bootstrap replications when calculating 95% confidence bands. Two-sample comparison of the transition probability P(X(t) = 2| X(0) = 1) between the groups defined by the variable `group` can be performed as follows
```
> patp_test(data=data, tmat=tmat, cid="cid", id="id",
+           group="group", h=1, j=2, s=0, B=1000)
```
It is recommended to use at least 1000 cluster bootstrap replications when performing two-sample hypothesis testing.
