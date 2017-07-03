#!/bin/bash -l

module load gcc/6.2.0 armadillo/7.400.2
Rscript scripts/load_results.R results/mass_hier.rds
