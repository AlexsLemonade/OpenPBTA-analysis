#!/bin/bash

cd "$(dirname "${BASH_SOURCE[0]}")"


Rscript histologies-sequence.R
Rscript histologies-cohort.R
Rscript histologies-method.R