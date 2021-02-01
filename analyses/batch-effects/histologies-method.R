library(tidyverse)
library(BatchQC)
library(sva)
grouper = function(df){
  
  # df = df[1:100,]
  
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
  return(df[1:100,])
}

run_batchQC = function(df_polya, df_stranded, id_batch, report_name){
  
  
  
  
  
  df_t = df_polya %>%
    gather(Kids_First_Biospecimen_ID, gene_expression, 2:ncol(df_polya)) %>%
    spread(gene_id, gene_expression)
  
  
  
  df_t2 = df_stranded %>%
    gather(Kids_First_Biospecimen_ID, gene_expression, 2:ncol(df_stranded)) %>%
    spread(gene_id, gene_expression)
  
  df_t = df_t %>%
    mutate("batch" = rep(1, times = nrow(df_t)))
  
  df_t2 = df_t2 %>%
    mutate("batch" = rep(2, times = nrow(df_t2)))
  
  
  data_t = full_join(df_t, df_t2)
  
  data_t = data_t %>% select(1, batch, everything())
  
  data_t = left_join(id_batch, data_t)

  # only include three most common histologies
  short_histo = names(sort(table((data_t$histology)), decreasing = TRUE)[1:3])
  data_t = filter(data_t, histology %in% short_histo)
  
  batch = as_factor(data_t$batch)
  levels(batch) <- 1:length(levels(batch))
  poly1matrix = t(data_t[4:ncol(data_t)])
  
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
  data_t$batch = as.factor(data_t$batch)
  
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
  
  # df = cbind(km$cluster, data_t)
  
  
  ggplot(tb, aes(x = PC1, y=PC2, color=data_t$histology, shape=data_t$batch)) +
    geom_point() +
    scale_color_manual(values = colors) +
    scale_shape_manual(values = shapes) +
    labs(color = "Histology", shape = "Batch") +
    theme_bw() 
  
  ggsave(paste(report_name, "batch-adjusted-histology.png", sep = "_"))
  
}

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


dat_rsem_polya = shorten(dat_rsem_polya)
dat_rsem_stranded = shorten(dat_rsem_stranded)

dat_kallisto_polya = grouper(dat_kallisto_polya)
dat_kallisto_stranded = grouper(dat_kallisto_stranded)


covariate = filter(covariate, experimental_strategy == "RNA-Seq")

id_batch = cbind("Kids_First_Biospecimen_ID" = covariate$Kids_First_Biospecimen_ID, "histology" = covariate$short_histology)
colnames(id_batch) <- c("Kids_First_Biospecimen_ID", "histology")
id_batch = as_tibble(id_batch)

run_batchQC(dat_rsem_polya, dat_rsem_stranded, id_batch, report_name = "rsem_method")
run_batchQC(dat_kallisto_polya, dat_kallisto_stranded, id_batch, report_name = "kallisto_method")

