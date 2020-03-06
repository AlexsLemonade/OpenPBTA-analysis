#!/bin/bash
# 
# Run all figure making scripts. 

# Find current directory based on this script
cd "$(dirname "${BASH_SOURCE[0]}")" 
 
# Make an output folder for all directories
mkdir figures

## Sample distribution

## Mutational landscape figure
Rscript figure-scripts/fig2-mutational-landscape.R

## Interaction plots 

## Oncoprint plot(s)

## Copy number status heatmap

## Transcriptomic overview
