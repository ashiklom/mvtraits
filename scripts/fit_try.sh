#!/bin/bash -l
#$ -q "geo*"
#$ -pe omp 8
#$ -j y
#$ -o logs/
#$ -N hier_area

module load gcc/6.2.0 armadillo/7.400.2
Rscript tests/fit2try.R
