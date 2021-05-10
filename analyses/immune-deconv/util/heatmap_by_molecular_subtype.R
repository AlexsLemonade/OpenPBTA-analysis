# function to create heatmap of average immune scores per cell type per molecular subtype per histology 
heatmap_by_molecular_subtype <- function(deconv_output){

  # plot title
  broad_histology <- unique(deconv_output$broad_histology)
  plot_title <- paste0("Average immune scores normalized by rows: ", broad_histology)
  
  # remove NA and to be classified
  deconv_output <- deconv_output %>%
    filter(!is.na(molecular_subtype),
           !grepl("To be classified", molecular_subtype))
  
  # create labels: count of samples per molecular subtype
  deconv_output <- deconv_output %>%
    group_by(broad_histology, molecular_subtype, cell_type) %>%
    unique() %>%
    mutate(label = n()) %>%
    mutate(label = paste0(molecular_subtype,' (',label,')'))
  
  # calculate average scores per cell type per molecular subtype
  deconv_output <- deconv_output %>%
    filter(!cell_type %in% c("microenvironment score", "stroma score", "immune score")) %>%
    group_by(cell_type, label) %>%
    dplyr::summarise(mean = mean(fraction)) %>%
    # convert into matrix of cell type vs molecular subtype
    spread(key = label, value = mean) %>%
    column_to_rownames('cell_type')
  
  # create heatmap per molecular subtype
  if(ncol(deconv_output) > 1){
    deconv_output <- deconv_output %>%
      rownames_to_column('celltype') %>%
      filter_at(vars(-celltype), any_vars(. != 0)) %>%
      column_to_rownames('celltype') %>%
      t() %>%
      pheatmap(fontsize = 10,
               scale = "column", angle_col = 45,
               main = plot_title,
               annotation_legend = T, cellwidth = 15, cellheight = 15)
  }
}