library(tidyverse)
library(BatchQC)
library(sva)
grouper = function(df){
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

run_batchQC = function(df_polya, df_stranded, id_histology, report_name){
  
  # create output directory (if necessary)
  dir.create(file.path(analysis_dir, "results"), showWarnings = FALSE)
  
  # create full base name of output file
  report_name = file.path(analysis_dir, "results", report_name)
  
  # tranpose data and change the index column name to Kids_First_Biospecimen_ID
  df_t = df_polya %>%
    gather(Kids_First_Biospecimen_ID, gene_expression, 2:ncol(df_polya)) %>%
    spread(gene_id, gene_expression)
  
  # tranpose data and change the index column name to Kids_First_Biospecimen_ID  
  df_t2 = df_stranded %>%
    gather(Kids_First_Biospecimen_ID, gene_expression, 2:ncol(df_stranded)) %>%
    spread(gene_id, gene_expression)
  
  # add a new column called batch and set all of the values to 1
  df_t = df_t %>%
    mutate("batch" = rep(1, times = nrow(df_t)))
  
  # add a new column called batch and set all of the values to 2
  df_t2 = df_t2 %>%
    mutate("batch" = rep(2, times = nrow(df_t2)))
  
  # join the two data sets together
  data_t = full_join(df_t, df_t2)
  
  # move the batch column to the front
  data_t = data_t %>% select(1, batch, everything())
  
  # add in the histology column
  data_t = left_join(id_histology, data_t)

  #### only include three most common histologies
  # create a list of the names of the three most common histologies
  short_histo = names(sort(table((data_t$histology)), decreasing = TRUE)[1:3])
  # filter the data to only include the top 3
  data_t = filter(data_t, histology %in% short_histo)
  
  # remove NA rows (index 274 had an NA with Kallisto)
  data_t = data_t[rowSums(is.na(data_t[4:ncol(data_t)])) == 0, ]
  
  #turn the batch column into a factor
  batch = as_factor(data_t$batch)
  batch = as.numeric(batch)
  
  #rename the numbers to be levels
  levels(batch) <- 1:length(levels(batch))
  
  print(dim(data_t))
  

  
  
  print(dim(data_t))
  
  #make a matrix representing just the gene expression data (this gets plugged into combat)
  poly1matrix = t(data_t[4:ncol(data_t)])
  
  pca = prcomp(log2(data_t[,4:ncol(data_t)]+ 1))
  
  tb = as_tibble(pca$x)
  
  # df = cbind(km$cluster, data_t)
  
  shapes = c(2,3,11,12,13,14)
  colors = c(
    "#a6611a",
    "#018571",
    "#c2a5cf",
    "#d01c8b",
    "#e66101",
    "#00192B",
    "#000000"
  )
  data_t$batch = as.factor(data_t$batch)
  
  ggplot(tb, aes(x = PC1, y=PC2, color=data_t$histology, shape=data_t$batch)) +
    geom_point() +
    scale_color_manual(values = colors) +
    scale_shape_manual(values = shapes) +
    labs(color = "Histology", shape = "batch") +
    theme_bw() 
  
  ggsave(paste(report_name, "histology.png", sep = "_"))
  
  
  var_zero_genes = poly1matrix[which((apply(poly1matrix, 1, var)==0)),]
  poly1matrix = poly1matrix[which((apply(poly1matrix, 1, var)>0)),]
  print(sum(var_zero_genes))
  
  
  
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


print("DOWNLOADED DATA")

dat_kallisto_polya = grouper(dat_kallisto_polya)
dat_kallisto_stranded = grouper(dat_kallisto_stranded)


covariate = filter(covariate, experimental_strategy == "RNA-Seq")

id_histology = cbind("Kids_First_Biospecimen_ID" = covariate$Kids_First_Biospecimen_ID, "histology" = covariate$short_histology)
colnames(id_histology) <- c("Kids_First_Biospecimen_ID", "histology")
id_histology = as_tibble(id_histology)

# run_batchQC(dat_rsem_polya, dat_rsem_stranded, id_histology, report_name = "rsem_method")
run_batchQC(dat_kallisto_polya, dat_kallisto_stranded, id_histology, report_name = "kallisto_method")

