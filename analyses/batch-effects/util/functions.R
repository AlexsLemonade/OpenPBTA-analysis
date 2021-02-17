# function used to cut files
shorten = function(df) {
  return(df[1:100,])
}

# function used on kallisto data for summarizing transcripts at the gene level
grouper = function(df){
  df = df %>%
    group_by(gene_id) %>%
    summarise_at(colnames(df)[3:ncol(df)], mean, na.rm = TRUE)
  return(df)
}

# transpose data through gather and spread functions, rename the index column
trans = function(data, Kids_First_Biospecimen_ID, val, starting_col){
  data = data %>%
    gather(Kids_First_Biospecimen_ID, val, starting_col:(ncol(data))) %>%
    spread("gene_id", val)
  return(data)
}

# transpose the data and rename columns
clean = function(data, gene_id){
  data_t = trans(data, "id", "gene_expression", 2)
  return(data_t)
}

# create directories and make sure paths will work on any computer
handle_file_outputs = function(report_name){
  # create output directory (if necessary)
  dir.create(file.path(analysis_dir, "results"), showWarnings = FALSE)
  
  # create full base name of output file
  report_name = file.path(analysis_dir, "results", report_name)
  
  return(report_name)
}

# make a PCA plot. Context name refers to the name that will be contactenated with the report name
make_plot = function(tb, data_t, report_name, context_name){
  data_t$batch = as_factor(data_t$batch)
  shapes = c(2,3,11,12)
  colors = c(
    "#a6611a",
    "#018571",
    "#c2a5cf",
    "#d01c8b",
    "#e66101"
  )
  
  p = ggplot(tb, aes(x = PC1, y=PC2, color=data_t$histology, shape=data_t$batch)) +
    geom_point() +
    scale_color_manual(values = colors) +
    scale_shape_manual(values = shapes) +
    labs(color = "Histology", shape = "batch") +
    theme_bw() 
  
  ggsave(paste(report_name, context_name, sep = "_"))
  
  return(p)
  
  
}

# join batch information with gene expression and filter NA's
tibble_combat_prep = function(df, gene_id, id_batch_histology){
  # clean the data (see clean comments)
  data_t = clean(df, gene_id)
  
  # join batch information with gene expression data  
  data_t = inner_join(id_batch_histology, data_t, by = "Kids_First_Biospecimen_ID")
  
  # filter out NAs
  data_t = filter(data_t, !is.na(data_t$batch))
  
  return(data_t)
}

# for batchQC to work, the batches must be numeric values
batch_prep = function(data_t){
  # verify that batch is indeed a factor and numeric
  batch = as_factor(data_t$batch)
  batch = as.numeric(batch)
  
  # make so that batch values start at 1
  levels(batch) <- 1:length(levels(batch))
  return(batch)
}

covariate_filtering = function(covariate, batch_column){
  
  # Only grab RNA-Seq data
  covariate = filter(covariate, experimental_strategy == "RNA-Seq")
  
  # Add a new column called batch which will be a duplicate of the whichever column name is passed into the function
  covariate$batch = unlist(covariate[which(names(covariate) == batch_column)])
  
  # convert 'batch' column into numbers through first converting it into a factor
  covariate$batch = as_factor(covariate$batch)
  
  # move the batch column to the 4th location (optional step which is useful for debugging)
  covariate = covariate %>% select(1:4, batch, everything())
  
  return(covariate)
}

# read in the covariate data and create a tibble with two columns: the patient ID and the Batch
make_id_batch = function(covariate, batch_column){
  
  # move the batch column to the front and only include RNA-Seq data
  covariate = covariate_filtering(covariate, batch_column) 

  # make a new matrix which will contain patient ID's and the corresponding batch
  id_batch = cbind("Kids_First_Biospecimen_ID" = covariate$Kids_First_Biospecimen_ID, "batch" = covariate$batch)

  # rename columns
  colnames(id_batch) <- c("Kids_First_Biospecimen_ID", "batch")

  # convert id_batch from a matrix to a tibble
  id_batch = as_tibble(id_batch)
  
  return(id_batch)
}

make_id_batch_histology = function(covariate, batch_column){
  # move the batch column to the front and only include RNA-Seq data
  covariate = covariate_filtering(covariate, batch_column) 
  
  # make a new matrix which will contain patient ID's and the corresponding batch and histology
  id_batch_histology = cbind("Kids_First_Biospecimen_ID" = covariate$Kids_First_Biospecimen_ID, "batch" = covariate$batch, "histology" = covariate$short_histology)
  
  #rename columns
  colnames(id_batch_histology) <- c("Kids_First_Biospecimen_ID", "batch", "histology")
  id_batch_histology = as_tibble(id_batch_histology)
  
  # convert id_batch_histology from a matrix to a tibble
  id_batch_histology = as_tibble(id_batch_histology)
  
  return(id_batch_histology)
}

# make two PCA plots before and after applying combat
make_before_and_after_combat_plots = function(data_t, batch, report_name, run_combat) {
  
  data_t$batch = batch_prep(data_t)
  
  #make a matrix representing just the gene expression data (this gets plugged into combat)
  poly1matrix = t(data_t[,4:ncol(data_t)])
  
  pca = prcomp(log2(data_t[,4:ncol(data_t)]+ 1))
  
  tb = as_tibble(pca$x)
  
  p1 = make_plot(tb, data_t, report_name, "histology.png")
  
  if (run_combat == FALSE){
    print("Combat set to not run on this data set")
    return()
  }
  
  # filter out genes without variance
  var_zero_genes = poly1matrix[which((apply(poly1matrix, 1, var)==0)),]
  poly1matrix = poly1matrix[which((apply(poly1matrix, 1, var)>0)),]
  print(sum(var_zero_genes))
  
  poly2matrix <- log2(poly1matrix + 1)
  
  
  combat_matrix = ComBat(dat=poly2matrix, batch=batch)
  
  pca = prcomp(t(combat_matrix))
  
  tb = as_tibble(pca$x)
  
  p2 = make_plot(tb, data_t, report_name, "batch-adjusted-histology.png")
  
  return(list(p1, p2))
  
}

# only include three most common histologies
histology_filtering = function(data_t){
  # create a list of the names of the three most common histologies
  short_histo = names(sort(table((data_t$histology)), decreasing = TRUE)[1:3])
  
  # filter the data to only include the top 3
  data_t = filter(data_t, histology %in% short_histo)
  
  return(data_t)
}

# root functions that calls sub functions to make all of the graphs
make_histology_pca_plots = function(df, id_batch_histology, gene_id, report_name, file_name, run_combat=TRUE){
  
  #handle file outputs
  report_name = handle_file_outputs(report_name)
 
  # clean the data see tible_combat_prep comments
  data_t = tibble_combat_prep(df, gene_id, id_batch_histology)
  
  # only include three most common histologies
  data_t = histology_filtering(data_t)
  
  # see batch_prep comments
  batch = batch_prep(data_t)
  
  print("DATA IS READY") 
  
  #Make before and after PCA plots
  return(make_before_and_after_combat_plots(data_t, batch, report_name, run_combat))

}

# root functions that calls sub functions to combat and batchQC looking for batch affects
run_batchQC = function(df, id_batch, gene_id, report_name, file_name, run_combat){
  
  # create output directory (if necessary)
  dir.create(file.path(analysis_dir, "results"), showWarnings = FALSE)
  output_dir = file.path(analysis_dir, "results")
  
  # clean the data see tible_combat_prep comments
  data_t = tibble_combat_prep(df, gene_id, id_batch)
  
  # run batchQC, save results, then run combat and save results
  run_batchQC_and_combat(data_t, output_dir, report_name, file_name, run_combat)
  
}

# join polyA and stranded data frames & create batch column
polyA_stranded_joining = function(df_polya, df_stranded){
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
  
  return(data_t)
}

# root functions that calls sub functions to combat and batchQC looking for batch affects || modified when polyA vs stranded is the batch
run_batchQC_polyA_vs_stranded = function(df_polya, df_stranded, report_name, file_name){
  
  # create output directory (if necessary)
  dir.create(file.path(analysis_dir, "results"), showWarnings = FALSE)
  output_dir = file.path(analysis_dir, "results")
  
  # join polyA and stranded data frames & create batch column
  data_t = polyA_stranded_joining(df_polya, df_stranded)
  
  # run batchQC, save results, then run combat and save results
  run_batchQC_and_combat(data_t, output_dir, report_name, file_name)
  
}

# root functions that calls sub functions to combat and batchQC looking for batch affects
run_batchQC_and_combat = function(data_t, output_dir, report_name, file_name, run_combat = TRUE) {
  
  # set up batcb vector
   batch = batch_prep(data_t)
  
  # make a matrix of the purly the gene expression data
  # this is needed as input for combat and batchQC
  poly1matrix = t(data_t[,3:ncol(data_t)])
  
  # pull out genes that have no variance (these will be added back later)
  var_zero_genes = poly1matrix[which((apply(poly1matrix, 1, var)==0)),]
  poly1matrix = poly1matrix[which((apply(poly1matrix, 1, var)>0)),]
  
  
  # log tranform data so batchQC and combat run properly
  poly1matrix = log2(poly1matrix + 1)
  
  print("DATA IS READY") 
  
  # run batchQC and only generate the P-value distribution and PCA plot
  batchQC(poly1matrix, batch=batch,
          report_file=report_name, report_dir=output_dir,
          report_option_binary="100001000",
          view_report=FALSE, interactive=FALSE, batchqc_output=FALSE)
  
  # certain gene expression files didn't have batch effects. Therefore we 
  # won't run combat for these.
  if (run_combat == TRUE){
    print("RUNNING COMBAT")
    
    combat_matrix = ComBat(dat=poly1matrix, batch=batch)
    
    # run batchQC to see if combat successfully removed batch effects
    batchQC(combat_matrix, batch=batch,
            report_file=paste("combat", report_name, sep="_"), report_dir=output_dir,
            report_option_binary="100001000",
            view_report=FALSE, interactive=FALSE, batchqc_output=FALSE)
    
    # undo log transformation of the data
    combat_matrix = 2^combat_matrix - 1
    
    # add back in genes with zero variance
    combat_matrix = rbind(combat_matrix,var_zero_genes)
    
    # save batch corrected gene expression data
    saveRDS(combat_matrix, file = file.path(output_dir,file_name))
  }
  
}

# root functions that calls sub functions to make all of the graphs
make_histology_pca_plots_methods = function(df_polya, df_stranded, id_histology, report_name, run_combat = TRUE){
  
  #handle file outputs
  report_name = handle_file_outputs(report_name)
  
  # join polyA and stranded data frames & create batch column
  data_t = polyA_stranded_joining(df_polya, df_stranded)
  
  # add in the histology column
  data_t = left_join(id_histology, data_t)
  
  # only include three most common histologies
  data_t = histology_filtering(data_t)
  
  # remove NA rows (index 274 had an NA with Kallisto)
  print(dim(data_t))
  data_t = data_t[rowSums(is.na(data_t[4:ncol(data_t)])) == 0, ]
  print(dim(data_t))
  
  #turn the batch column into a factor
  batch = batch_prep(data_t)
  
  print("DATA IS READY") 
  
  #Make before and after PCA plots
  make_before_and_after_combat_plots(data_t, batch, report_name, run_combat)
  
}

