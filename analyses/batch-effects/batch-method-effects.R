library(tidyverse)
library(BatchQC)
library(sva)
library(stringr)


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



grouper = function(df){
  df = df %>%
    group_by(gene_id) %>%
    summarise_at(colnames(df)[3:ncol(df)], mean, na.rm = TRUE)
  return(df)
}




run_batchQC = function(df_polya, df_stranded, report_name, file_name){
  
 output_dir = here('analyses', 'batch-effects','results')
 
 #transpose first tibble
  df_t = df_polya %>%
    gather(Kids_First_Biospecimen_ID, gene_expression, 2:ncol(df_polya)) %>%
    spread(gene_id, gene_expression)


 # transpose second tibble
  df_t2 = df_stranded %>%
    gather(Kids_First_Biospecimen_ID, gene_expression, 2:ncol(df_stranded)) %>%
    spread(gene_id, gene_expression)

 # add a new column called batch and add one because first polyA tibble represents first batch
  df_t = df_t %>%
    mutate("batch" = rep(1, times = nrow(df_t)))

 # add a new batch column and insert values 2
  df_t2 = df_t2 %>%
    mutate("batch" = rep(2, times = nrow(df_t2)))
  
 # join the two tibbles together
  data_t = full_join(df_t, df_t2)

 # move the batches to the front
  data_t = data_t %>% select(1, batch, everything())
  
 # make batch a factor
  batch = as_factor(data_t$batch)

 # verify that batches start at 1
  levels(batch) <- 1:length(levels(batch))

 # make a transposed matrix of just gene expression data 
  poly1matrix = t(data_t[3:ncol(data_t)])

 # remove genes with zero varation (but will be added back in the end)
  var_zero_genes = poly1matrix[which((apply(poly1matrix, 1, var)==0)),]
  poly1matrix = poly1matrix[which((apply(poly1matrix, 1, var)>0)),]

 # log transform data so that batchQC and combat work correctly
  poly1matrix = log2(poly1matrix + 1)

  print("DATA IS READY") 
  batchQC(poly1matrix, batch=batch,
          report_file=report_name, report_dir=output_dir,
          report_option_binary="100001000",
          view_report=FALSE, interactive=FALSE, batchqc_output=FALSE)
  
  
  print("RUNNING COMBAT")

  combat_matrix = ComBat(dat=poly1matrix, batch=batch)

  batchQC(combat_matrix, batch=batch,
          report_file=paste("combat", report_name, sep="_"), report_dir=output_dir,
          report_option_binary="100001000",
          view_report=FALSE, interactive=FALSE, batchqc_output=FALSE)

  combat_matrix = 2^combat_matrix - 1

  combat_matrix = rbind(combat_matrix,var_zero_genes)

  saveRDS(combat_matrix, file = file.path(output_dir,file_name))

}

# START READING CODE COMMENTS HERE
library(here) 
path = here("data", "release-v13-20200116")
setwd(path)

# download covariate data
covariate = read_tsv("pbta-histologies.tsv",col_types = cols(molecular_subtype = "c"))

# download gene expression data
dat_rsem_polya = readRDS("pbta-gene-expression-rsem-tpm.polya.rds")
dat_rsem_stranded = readRDS("pbta-gene-expression-rsem-tpm.stranded.rds")
dat_kallisto_stranded <- readRDS("pbta-gene-expression-kallisto.stranded.rds")
dat_kallisto_stranded = dat_kallisto_stranded[,2:ncol(dat_kallisto_stranded)]
dat_kallisto_polya <- readRDS("pbta-gene-expression-kallisto.polya.rds")

print("DOWNLOADED DATA")

# summarize transcript data at the gene level
dat_kallisto_polya = grouper(dat_kallisto_polya)
dat_kallisto_stranded = grouper(dat_kallisto_stranded)

print("GROUPED TRANSCRIPTS")


("PREPARING TO RUN ANALYSIS")

run_batchQC(dat_rsem_polya, dat_rsem_stranded, report_name = "rsem_method_report.html","pbta-gene-expression-rsem-tpm.rds")
run_batchQC(dat_kallisto_polya, dat_kallisto_stranded, report_name = "kallisto_method_report.html","pbta-gene-expression-kallisto.rds")

