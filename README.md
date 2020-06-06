# clustered-multistate
R code for nonparametric analysis of clustered multistate processes

## Description

This repository contains R functions for nonparametric analysis of (continuous-time) multistate processes with cluster-correlated observations. These functions currently support nonparametric estimation of population-averaged transition probabilities, calculation of 95% pointwise confidence intervals and simultaneous confidence bands, and two-sample Kolmogorov-Smirnov-type tests. No assumptions regarding the within-cluster dependence are imposed and the functions can be for both Markov and non-Markov processes. Right censoring, left truncation, and informative cluster size are allowed.

## Main Functions

The main functions are `wiaj()` and `wiaj_test()`. Both functions are beta versions.

### Dependencies
The functions utilize the R packages `mstate` and `Rfast`.

### Input data

### Description of `wiaj()`

### Description of `wiaj_test()`


## Examples

