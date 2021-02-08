library(tidyverse)
library(BatchQC)
library(sva)
library(stringr)




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
# Get correct file paths to data
library(rprojroot)
root_dir = find_root(has_file("OpenPBTA-analysis.Rproj"))
analysis_dir = file.path(root_dir, "analyses", "batch-effects")
data_dir = file.path(root_dir, "data")
functions = file.path(analysis_dir, "util", "functions.R")
source(functions)

# download all covariate data which will be used to identify batches
covariate_file = file.path(data_dir, "pbta-histologies.tsv")
covariate = read_tsv(covariate_file, col_types = cols(molecular_subtype = "c"))

# download gene expression data
dat_rsem_polya_file = file.path(data_dir, "pbta-gene-expression-rsem-tpm.polya.rds")
dat_rsem_polya = readRDS(dat_rsem_polya_file)

dat_rsem_stranded_file = file.path(data_dir, "pbta-gene-expression-rsem-tpm.stranded.rds")
dat_rsem_stranded = readRDS(dat_rsem_stranded_file)

dat_kallisto_stranded_file = file.path(data_dir, "pbta-gene-expression-kallisto.stranded.rds")
dat_kallisto_stranded <- readRDS(dat_kallisto_stranded_file)
dat_kallisto_stranded = dat_kallisto_stranded[,2:ncol(dat_kallisto_stranded)]

dat_kallisto_polya_file = file.path(data_dir, "pbta-gene-expression-kallisto.polya.rds")
dat_kallisto_polya <- readRDS(dat_kallisto_polya_file)

print("DOWNLOADED DATA")

### coment out section below when running the full analysis ###

dat_rsem_polya = shorten(dat_rsem_polya)
dat_rsem_stranded = shorten(dat_rsem_stranded)
dat_kallisto_stranded = shorten(dat_kallisto_stranded)
dat_kallisto_polya = shorten(dat_kallisto_polya)
print("DATA HAS BEEN SHORTENED")

### coment out section above when running the full analysis ###


# summarize transcript data at the gene level
dat_kallisto_polya = grouper(dat_kallisto_polya)
dat_kallisto_stranded = grouper(dat_kallisto_stranded)

print("GROUPED TRANSCRIPTS")


("PREPARING TO RUN ANALYSIS")

run_batchQC(dat_rsem_polya, dat_rsem_stranded, report_name = "rsem_method_report.html","pbta-gene-expression-rsem-tpm.rds")
run_batchQC(dat_kallisto_polya, dat_kallisto_stranded, report_name = "kallisto_method_report.html","pbta-gene-expression-kallisto.rds")

