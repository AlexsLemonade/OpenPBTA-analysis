#!/usr/bin/env Rscript
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

bioc_repos <- BiocManager::repositories()
bioc_repos <- bioc_repos[names(bioc_repos) != "CRAN"]

# We want errors not just warnings
options(warn = 2)
install.packages(opts$PACKAGES, repos = c(bioc_repos, getOption('repos')))