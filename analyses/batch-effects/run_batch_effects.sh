#!/bin/bash

cd "$(dirname "${BASH_SOURCE[0]}")"

Rscript batch-sequence-effects.R
Rscript batch-cohort-effects.R
Rscript batch-method-effects.R
