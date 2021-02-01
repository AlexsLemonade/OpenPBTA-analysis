library(tidyverse)
library(BatchQC)
library(sva)

#dat_rsem_stranded = shorten(dat_rsem_stranded)
getCurrentFileLocation <-  function()
{
  this_file <- commandArgs() %>% 
    tibble::enframe(name = NULL) %>%
    tidyr::separate(col=value, into=c("key", "value"), sep="=", fill='right') %>%
    dplyr::filter(key == "--file") %>%
    dplyr::pull(value)
  if (length(this_file)==0)
  {
    this_file <- rstudioapi::getSourceEditorContext()$path
  }
  return(dirname(this_file))
}

# function used on kallisto data for summarizing transcripts at the gene level
grouper = function(df){
  df = df %>%
    group_by(gene_id) %>%
    summarise_at(colnames(df)[3:ncol(df)], mean, na.rm = TRUE)
  return(df)
}




trans = function(data, Kids_First_Biospecimen_ID, val, starting_col){
  # transpose data through gather and spread functions
  data = data %>%
    gather(Kids_First_Biospecimen_ID, val, starting_col:(ncol(data))) %>%
    spread("gene_id", val)
  return(data)
}

clean = function(data, gene_id){
  # transpose the data and rename the columns
  data_t = trans(data, "id", "gene_expression", 2)
  return(data_t)
}

run_batchQC = function(df, id_batch, gene_id, report_name, file_name, run_combat){

  # create output path
  output_dir = here('analyses', 'batch-effects','results')
 
  # clean the data (see clean comments)
  data_t = clean(df, gene_id)
 
  # join batch information with gene expression data  
  data_t = inner_join(id_batch, data_t, by = "Kids_First_Biospecimen_ID")
 
  # filter out NAs
  data_t = filter(data_t, !is.na(data_t$batch))

  # only include most common histologies
  short_histo = names(sort(table((data_t$histology)), decreasing = TRUE)[1:3])
  data_t = filter(data_t, histology %in% short_histo)

  # verify that batch is indeed a factor and numeric
  batch = as_factor(data_t$batch)
  batch = as.numeric(batch)

  # make so that batch values start at 1
  levels(batch) <- 1:length(levels(batch))

  # make a matrix of the purly the gene expression data
  # this is needed as input for combat and batchQC
  poly1matrix = t(data_t[,4:ncol(data_t)])

  # pull out genes that have no variance (these will be added back later)
  var_zero_genes = poly1matrix[which((apply(poly1matrix, 1, var)==0)),]
  poly1matrix = poly1matrix[which((apply(poly1matrix, 1, var)>0)),]


  print("DATA IS READY") 
  poly1matrix = t(data_t[,4:ncol(data_t)])
  
  
  # km = kmeans(expression)
  
  pca = prcomp(log2(data_t[,4:ncol(data_t)]+ 1))
  
  tb = as_tibble(pca$x)
  
  # df = cbind(km$cluster, data_t)
  
  shapes = c(2,3,11,12)
  colors = c(
    "#a6611a",
    "#018571",
    "#c2a5cf",
    "#d01c8b",
    "#e66101"
  )
  
  ggplot(tb, aes(x = PC1, y=PC2, color=data_t$histology, shape=data_t$batch)) +
    geom_point() +
    scale_color_manual(values = colors) +
    scale_shape_manual(values = shapes) +
    labs(color = "Histology", shape = "batch") +
    theme_bw() 
  
  ggsave(paste(report_name, "histology.png", sep = "_"))
  
  
  var_zero_genes = poly1matrix[which((apply(poly1matrix, 1, var)<0.1)),]
  poly1matrix = poly1matrix[which((apply(poly1matrix, 1, var)>0)),]
  print(length(poly1matrix[which((apply(poly1matrix, 1, var)<0.1)),]))
  
  
  poly2matrix <- log2(poly1matrix + 1)
  

  combat_matrix = ComBat(dat=poly2matrix, batch=batch)
  
  pca = prcomp(t(combat_matrix))
  
  tb = as_tibble(pca$x)
  
  
  
  ggplot(tb, aes(x = PC1, y=PC2, color=data_t$histology, shape=data_t$batch)) +
    geom_point() +
    scale_color_manual(values = colors) +
    scale_shape_manual(values = shapes) +
    labs(color = "Histology", shape = "Batch") +
    theme_bw() 
  
  ggsave(paste(report_name, "batch-adjusted-histology.png", sep = "_"))  

}

# START HERE WHILE READING CODE COMMENTS


library(here)
library(rprojroot)
path = here("data")
print(path)
setwd(path)

# download all covariate data which will be used to identify batches
covariate = read_tsv("pbta-histologies.tsv",col_types = cols(molecular_subtype = "c"))

# download gene expression data
dat_rsem_polya = readRDS("pbta-gene-expression-rsem-tpm.polya.rds")
dat_rsem_stranded = readRDS("pbta-gene-expression-rsem-tpm.stranded.rds")
dat_kallisto_stranded <- readRDS("pbta-gene-expression-kallisto.stranded.rds")
dat_kallisto_stranded = dat_kallisto_stranded[,2:ncol(dat_kallisto_stranded)]
dat_kallisto_polya <- readRDS("pbta-gene-expression-kallisto.polya.rds")


# Summarize transcripts at the gene level by averageing all transcipts values per gene
# This step only needs to be done with kallisto gene expression files
dat_kallisto_polya = grouper(dat_kallisto_polya)
dat_kallisto_stranded = grouper(dat_kallisto_stranded)

print("DATA IS GROUPED")

# Only grab RNA-seq data
covariate = filter(covariate, experimental_strategy == "RNA-Seq")

# Add a new column called batch which will be a duplicate of the cohort column
covariate = mutate(covariate, "batch" = covariate$cohort, na.rm = TRUE)

# convert 'batch' column into numbers through first converting it into a factor
covariate$batch = as_factor(covariate$batch)

# move the batch column to the 4th location (optional step which is useful for debugging)
covariate = covariate %>% select(1:4, batch, everything())

# make a new matrix which will contain patient ID's and the corresponding batch and histology
id_batch = cbind("Kids_First_Biospecimen_ID" = covariate$Kids_First_Biospecimen_ID, "batch" = covariate$batch, "histology" = covariate$short_histology)

#rename columns
colnames(id_batch) <- c("Kids_First_Biospecimen_ID", "batch", "histology")
id_batch = as_tibble(id_batch)

# convert id_batch from a matrix to a tibble
id_batch = as_tibble(id_batch)

# Run batchQC and combat. This will run batchQC before and after running combat. It will also perform join operations between id_batch and our gene expression files 
run_batchQC(dat_rsem_polya, id_batch, gene_id = dat_rsem_polya$gene_id, report_name = "rsem_poly_cohort","pbta-gene-expression-rsem-tpm-combat-seq-center.polya.rds")
#run_batchQC(dat_rsem_stranded, id_batch, gene_id = dat_rsem_stranded$gene_id, report_name = "rsem_stranded_sequence", "pbta-gene-expression-rsem-tpm-combat-seq-center.stranded.rds")
run_batchQC(dat_kallisto_polya, id_batch, gene_id = dat_kallisto_polya$gene_id, report_name = "kallisto_polya_cohort", "pbta-gene-expression-kallisto-combat-seq-center.polya.rds")
#run_batchQC(dat_kallisto_stranded, id_batch, gene_id = dat_kallisto_stranded$gene_id, report_name = "kallisto_stranded_sequence", "pbta-gene-expression-kallisto-combat-seq-center.stranded.rds")

