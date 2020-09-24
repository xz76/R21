#!/bin/bash

parallel -j 6 Rscript one_sim.R ::: {1..1000} &
