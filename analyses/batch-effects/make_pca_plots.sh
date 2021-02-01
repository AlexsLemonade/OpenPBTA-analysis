#!/bin/bash

cd "$(dirname "${BASH_SOURCE[0]}")"


rscript cohort-confounding.R
rscript method-confounding.R
rscript seq_center-confounding.R

