#!/bin/bash

#SBATCH --time=02:00:00   # walltime
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH -C 'rhel7'   # features syntax (use quotes): -C 'a&b&c&d'
#SBATCH --mem=120G   # memory 

Rscript batch-sequence-effects.R
Rscript batch-cohort-effects.R
Rscript batch-method-effects.R
