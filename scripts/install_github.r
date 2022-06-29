#!/usr/bin/env Rscript
#
# Installs R packages from github, optionally using a PAT stored in a file
#
# Josh Shapiro for CCDL
#
# Install packages from github using remotes::install_github()
# Allows specification of a file containing a GitHub PAT to avoid rate limiting.


library(docopt)

doc <- "Usage: 
  install_github.r <repository> [options]

Options:
  -h --help           show this help text
  --ref <reference>   optional commit reference
  --pat_file <file>   a file containing a GitHub Personal Access Token
  --no_deps           don't install dependencies
"

opts <- docopt(doc)

if(is.null(opts$pat_file)){
  pat <- NULL
}else{
  pat <- scan(opts$pat_file, what = "character", n = 1)
}

remotes::install_github(repo = opts$repository,
                        ref = opts$ref,
                        dependencies = !opts$no_deps,
                        auth_token = pat,
                        upgrade = FALSE)
