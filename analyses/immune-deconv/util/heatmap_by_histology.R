# function to create heatmap of average immune scores per cell type per histology 
heatmap_by_histology <- function(deconv_output, output_file) {
  
  # create labels: count of samples per histology
  deconv_output <- deconv_output %>%
    group_by(broad_histology, cell_type) %>%
    unique() %>%
    mutate(label = n()) %>%
    mutate(label = paste0(broad_histology,' (',label,')'))
  
  # calculate average scores per cell type per histology
  deconv_output <- deconv_output %>%
    filter(!cell_type %in% c("microenvironment score", "stroma score", "immune score")) %>%
    group_by(cell_type, label) %>%
    dplyr::summarise(mean = mean(fraction)) %>%
    # convert into matrix of cell type vs histology
    spread(key = label, value = mean) %>%
    column_to_rownames('cell_type')
  
  # plot non-brain and brain tumors separately
  pdf(file = output_file, width = 15, height = 8)
  
  # non brain tumors
  non_brain_tumors <- c("Histiocytic tumor", "Lymphoma")
  mat <- deconv_output %>%
    dplyr::select(grep(paste0(non_brain_tumors, collapse="|"), colnames(deconv_output), value = TRUE))
  if(ncol(mat) > 1){
    mat <- mat %>%
      rownames_to_column('celltype') %>%
      filter_at(vars(-celltype), any_vars(. != 0)) %>%
      column_to_rownames('celltype') %>%
      t() %>%
      pheatmap(fontsize = 10,
               scale = "column", angle_col = 45,
               main = "Average immune scores normalized by rows\nNon-Brain Tumors",
               annotation_legend = T, cellwidth = 15, cellheight = 15)
  }
  
  # brain tumors
  mat <- deconv_output %>%
    dplyr::select(grep(paste0(non_brain_tumors, collapse="|"), colnames(deconv_output), invert = TRUE, value = TRUE))
  if(ncol(mat) > 1){
    mat <- mat %>%
      rownames_to_column('celltype') %>%
      filter_at(vars(-celltype), any_vars(. != 0)) %>%
      column_to_rownames('celltype') %>%
      t() %>%
      pheatmap(fontsize = 10,
               scale = "column", angle_col = 45,
               main = "Average immune scores normalized by rows\nBrain Tumors",
               annotation_legend = T, cellwidth = 15, cellheight = 15)
  }
  
  dev.off()
}