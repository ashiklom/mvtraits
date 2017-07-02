#!/bin/bash
set -e
mkdir -p 'prof'
Rscript -e 'devtools::install(".")'

R -d "valgrind --tool=callgrind" -f tests/testthat/test.fit_mvnorm_hier.R
