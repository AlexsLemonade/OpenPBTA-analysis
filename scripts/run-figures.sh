#!/bin/bash
# 
# Run all figure making scripts. 


cd "$(dirname "${BASH_SOURCE[0]}")" 
 
mkdir figures

## Sample distribution

## Tumor mutation burden 
Rscript fig2-mutational-landscape.R

## Interaction plots 

## Oncoprint plot(s)

## Copy number status heatmap

## Transcriptomic overview
