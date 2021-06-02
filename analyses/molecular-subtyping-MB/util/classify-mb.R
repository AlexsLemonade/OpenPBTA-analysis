# Author: Komal S. Rathi
# Date: 07/22/2020
# Function: Function to classify MB subtypes

# load libraries
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(medulloPackage))
suppressPackageStartupMessages(library(MM2S))
suppressPackageStartupMessages(library(org.Hs.eg.db))

# function to run molecular subtyping
classify.mb <- function(expr.input, method){
  # get object
  expr.mb <- get(expr.input)
  
  if(method == "medullo-classifier"){
    # run medulloPackage
    res <- medulloPackage::classify(expr.mb)
    res <- res %>%
      mutate(classifier = method,
             dataset = expr.input) %>%
      arrange(sample)
  } else {
    # MM2S package uses ENTREZ gene ids whereas expression data contains gene symbols
    hs <- org.Hs.eg.db
    mapping <- select(hs, 
                      keys = rownames(expr.mb),
                      columns = c("ENTREZID", "SYMBOL"),
                      keytype = "SYMBOL")
    mapping <- mapping %>%
      filter(!is.na(ENTREZID)) 
    
    # replace gene symbols with ENTREZ gene id as rownames
    expr.mb <- expr.mb %>%
      as.data.frame() %>%
      rownames_to_column('SYMBOL') %>%
      inner_join(mapping, by = 'SYMBOL') %>%
      dplyr::select(-c(SYMBOL)) %>%
      column_to_rownames('ENTREZID')
    
    # run MM2S
    res <- MM2S.human(InputMatrix = expr.mb,
                      parallelize = 4,
                      seed = 12345, tempdir())
    
    # get scores
    scores <- res$Predictions %>% 
      as.data.frame() %>%
      rownames_to_column('sample') %>% 
      gather(best.fit, score, -sample)
    
    # get predicted subtype
    subtype <- res$MM2S_Subtype %>%
      as.data.frame() %>%
      mutate(MM2S_Prediction = stringr::str_replace(MM2S_Prediction, "NORMAL", "Normal"), 
             SampleName = as.character(SampleName))
    
    # merge predicted subtype and associated scores
    res <- scores %>%
      inner_join(subtype, by = c("sample" = "SampleName", 
                                 "best.fit" = "MM2S_Prediction")) %>%
      mutate(classifier = method,
             dataset = expr.input) %>%
      arrange(sample)
  }
}
