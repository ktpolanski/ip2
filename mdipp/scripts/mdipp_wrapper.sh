#!/bin/bash
set -e

#start by running MDI++
/mdipp/mdi++ $@ > mcmc_chains.csv

#follow it up by running the R post-processer
Rscript /scripts/mdipp.R $@

#and now generate some functional analysis inputs
mkdir functional_analysis_inputs
python3 /scripts/bingomeme.py