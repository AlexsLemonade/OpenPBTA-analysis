# Author: Komal S. Rathi
# Function: Function to classify MB subtypes

# load libraries
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(medulloPackage))

# function to run molecular subtyping
classify_mb <- function(exprs_input, method){

  # get object
  exprs_mb <- get(exprs_input)

  # run medulloPackage
  res <- medulloPackage::classify(exprs_mb)
  res <- res %>%
    mutate(classifier = method,
           dataset = exprs_input) %>%
    arrange(sample) 
}
