# Author: Komal S. Rathi
# Date: 07/22/2020
# Function:
# Script to compare subtype classification with known findings

# load libraries
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))

option_list <- list(
  make_option(c("--observed_class"), type = "character",
              help = "Observed class types (.rds)"),
  make_option(c("--expected_class"), type = "character",
              help = "Expected class types (.rds)"),
  make_option(c("--outputfile"), type = "character",
              help = "Comparison Output (.rds)")
)

# parse parameters
opt <- parse_args(OptionParser(option_list = option_list))
obs.class <- opt$observed_class
exp.class <- opt$expected_class
outputfile <- opt$outputfile

# read data
obs.class <- readRDS(obs.class)
exp.class <- readRDS(exp.class)

# merge expected and observed data
output <- obs.class %>%
  left_join(exp.class, by = c('sample_id'))

saveRDS(output, file = outputfile)
