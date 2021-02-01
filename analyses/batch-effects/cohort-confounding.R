library(tidyverse)
library(BatchQC)
library(sva)
library(here)

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
  # df = df[1:200,]
  df = df %>%
    group_by(gene_id) %>%
    summarise_at(colnames(df)[3:ncol(df)], mean, na.rm = TRUE)
 
  return(df)
}


remove_char = function(a){
  b = NULL
  for (x in a){
    if (!is.na(as.numeric(x))){
      b = c(b, x)
    }
  }
  return(b)
}



trans = function(data, Kids_First_Biospecimen_ID, val, starting_col){
  data = data %>%
    gather(Kids_First_Biospecimen_ID, val, starting_col:(ncol(data))) %>%
    spread("gene_id", val)
  return(data)
}

clean = function(data, gene_id){
  data_t = trans(data, "id", "gene_expression", 2)
  return(data_t)
}

shorten = function(df){
  return(df[1:100,1:20])
}

run_batchQC = function(df, id_batch, gene_id, report_name){
  
  # data_t = clean(dat_rsem_polya, "gene_name" = rownames(dat_rsem_polya))
  
  
  
  data_t = clean(df, gene_id)
  
  
  data_t = inner_join(id_batch, data_t, by = "Kids_First_Biospecimen_ID")
  
  data_t = filter(data_t, !is.na(data_t$batch))
  
  
  batch = as_factor(data_t$batch)
  batch = as.numeric(batch)
  levels(batch) <- 1:length(levels(batch))
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
    labs(color = "Histology", shape = "batch") +
    theme_bw() 
    
  ggsave(paste(report_name, "histology.png", sep = "_"))

 
  # 
  # ggplot(tb, aes(x = PC1, y = PC2, color = as.factor(km$cluster))) +
  #   geom_point() +
  #   labs(color = "Cluster") +
  #   theme_bw()
  # 
  # ggsave("kmeans.png")
  
  # 
  var_zero_genes = poly1matrix[which((apply(poly1matrix, 1, var)<0.1)),]
  poly1matrix = poly1matrix[which((apply(poly1matrix, 1, var)>0)),]
  print(length(poly1matrix[which((apply(poly1matrix, 1, var)<0.1)),]))


  poly2matrix <- log2(poly1matrix + 1)

  # 
  # 
  # print("DATA IS READY") 
  # batchQC(poly1matrix, batch=batch,
  #         report_file=report_name, report_dir=".",
  #         report_option_binary="100001000",
  #         view_report=FALSE, interactive=FALSE, batchqc_output = FALSE)
  # 
  # 
  # 
  combat_matrix = ComBat(dat=poly2matrix, batch=batch)
  
  pca = prcomp(t(combat_matrix))
  
  tb = as_tibble(pca$x)
  
  # df = cbind(km$cluster, data_t)
  
  
  ggplot(tb, aes(x = PC1, y=PC2, color=data_t$histology, shape=data_t$batch)) +
    geom_point() +
    scale_color_manual(values = colors) +
    labs(color = "Histology", shape = "Batch") +
    theme_bw() 
    
  ggsave(paste(report_name, "batch-adjusted-histology.png", sep = "_"))
  # 
  # 
  # batchQC(combat_matrix, batch=batch,
  #         report_file=paste("combat", report_name, sep="_"), report_dir=".",
  #         report_option_binary="100001000",
  #         view_report=FALSE, interactive=FALSE, batchqc_output=TRUE)
  # 
  # print(combat_matrix)
  # print("MIN EQUALS")
  # print(min(poly1matrix))
  # print("Max Equals")
  # print(max(poly1matrix))
}

library(here)
library(rprojroot)
path = here("data")
print(path)
setwd(path)

covariate = read_tsv("pbta-histologies.tsv",col_types = cols(molecular_subtype = "c"))

dat_rsem_polya = readRDS("pbta-gene-expression-rsem-tpm.polya.rds")
dat_rsem_stranded = readRDS("pbta-gene-expression-rsem-tpm.stranded.rds")
dat_kallisto_stranded <- readRDS("pbta-gene-expression-kallisto.stranded.rds")
dat_kallisto_stranded = dat_kallisto_stranded[,2:ncol(dat_kallisto_stranded)]
dat_kallisto_polya <- readRDS("pbta-gene-expression-kallisto.polya.rds")


dat_kallisto_polya = grouper(dat_kallisto_polya)
dat_kallisto_stranded = grouper(dat_kallisto_stranded)

print("DATA IS GROUPED")

covariate = filter(covariate, experimental_strategy == "RNA-Seq")
covariate = mutate(covariate, "batch" = covariate$cohort)
covariate$batch = as_factor(covariate$batch)
covariate = covariate %>% select(1:4, batch, everything())

short_histo = names(sort(table((covariate$short_histology)), decreasing = TRUE)[1:3])
# long_histo = (names(sort(table((covariate$broad_histology)), decreasing = TRUE)[1:3]))

covariate = filter(covariate, short_histology %in% short_histo)

id_batch = cbind("Kids_First_Biospecimen_ID" = covariate$Kids_First_Biospecimen_ID, "batch" = covariate$batch, "histology" = covariate$short_histology)
colnames(id_batch) <- c("Kids_First_Biospecimen_ID", "batch", "histology")
id_batch = as_tibble(id_batch)



run_batchQC(dat_rsem_polya, id_batch, gene_id = dat_rsem_polya$gene_id, report_name = "rsem_polya_cohort")
run_batchQC(dat_rsem_stranded, id_batch, gene_id = dat_rsem_stranded$gene_id, report_name = "rsem_stranded_cohort")
run_batchQC(dat_kallisto_polya, id_batch, gene_id = dat_kallisto_polya$gene_id, report_name = "kallisto_polya_cohort")
run_batchQC(dat_kallisto_stranded, id_batch, gene_id = dat_kallisto_stranded$gene_id, report_name = "kallisto_stranded_cohort")




