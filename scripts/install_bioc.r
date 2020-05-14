#!/usr/bin/env r
#
# Installs bioconductor packages, but without including
#
# Josh Shapiro for CCDL
#
# Inspired by install2.r from littler by Carl Boettiger, Dirk Eddelbuettel, and Brandon Bertelsen
# and this tweet from Carl Boettiger https://twitter.com/cboettig/status/1260766721515782145?s=20

library(docopt)

doc <- "Usage: install_bioc.r [-h] [PACKAGES ...]

-h --help           show this help text"

opts <- docopt(doc)

bioc_repo <- BiocManager::repositories()[1]

install.packages(opts$PACKAGES, repos = c(bioc_repo, getOption('repos')))