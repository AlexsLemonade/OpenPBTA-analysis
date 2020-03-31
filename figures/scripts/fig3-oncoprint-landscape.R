# Oncoprint Landscape Figure
#
# 2020
# Chante Bethell for ALSF - CCDL
#
# This script is intended to run steps needed to create Figure 3.

#### Set Up --------------------------------------------------------------------

# Install maftools
if (!("maftools" %in% installed.packages())) {
  install.packages("BiocManager")
  BiocManager::install("maftools")
}
library(maftools)

# Load in patchwork for assembling the final multipanel figure
library(patchwork)

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

#### Directories and Files -----------------------------------------------------

# Detect the ".git" folder -- this will in the project root directory.
# Use this as the root directory to ensure proper execution, no matter where
# it is called from.
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Path to the data obtained via `bash download-data.sh`.
data_dir <- file.path(root_dir, "data")

# Declare output directory
output_dir <- file.path(root_dir, "figures", "pngs")

# Source the color palette for plots
source(
  file.path(
    root_dir,
    "analyses",
    "oncoprint-landscape",
    "util",
    "oncoplot-palette.R"
  )
)

#### Read in data --------------------------------------------------------------

# Read in metadata
metadata <-
  readr::read_tsv(file.path(root_dir, "data", "pbta-histologies.tsv")) %>%
  dplyr::rename(Tumor_Sample_Barcode = sample_id)

# Read in MAF files
maf_df_primary_plus <-
  readr::read_tsv(
    file.path(
      root_dir,
      "scratch",
      "oncoprint_files",
      "all_participants_primary-plus_maf.tsv"
    )
  )

maf_df_primary_only <-
  readr::read_tsv(
    file.path(
      root_dir,
      "scratch",
      "oncoprint_files",
      "all_participants_primary-only_maf.tsv"
    )
  )

# Read in cnv files
cnv_file_primary_plus <-
  readr::read_tsv(
    file.path(
      root_dir,
      "scratch",
      "oncoprint_files",
      "all_participants_primary-plus_cnv.tsv"
    )
  )

cnv_file_primary_only <-
  readr::read_tsv(
    file.path(
      root_dir,
      "scratch",
      "oncoprint_files",
      "all_participants_primary-only_cnv.tsv"
    )
  )

# Read in fusion files
fusion_file_primary_plus <- 
  readr::read_tsv(
    file.path(
      root_dir,
      "scratch",
      "oncoprint_files",
      "all_participants_primary-plus_fusions.tsv"
    )
  )

fusion_file_primary_only <- 
  readr::read_tsv(
    file.path(
      root_dir,
      "scratch",
      "oncoprint_files",
      "all_participants_primary-only_fusions.tsv"
    )
  )

# Read in gene list
gene_list <-
  readr::read_tsv(
    file.path(
      root_dir,
      "analyses",
      "interaction-plots",
      "results",
      "gene_disease_top50.tsv"
    )
  )

#### Custom function ----------------------------------------------------------

prepare_and_plot_oncoprint <- function(maf_df,
                                       cnv_file,
                                       fusion_file,
                                       gene_list) {
  # Given the maf, cnv, and fusion data.frames prepared in `00-map-to-sample_id.R`,
  # along with a list of genes prepared in the `interaction-plots` directory,
  # plot an oncoprint landscape figure.
  
  # Bind rows of maf and fusion data frames
  maf_df <- dplyr::bind_rows(maf_df, fusion_file)
  
  # Convert into MAF object
  maf_object <-
    read.maf(
      maf = maf_df,
      clinicalData = metadata,
      cnTable = cnv_file,
      removeDuplicatedVariants = FALSE,
      vc_nonSyn = c(
        "Frame_Shift_Del",
        "Frame_Shift_Ins",
        "Splice_Site",
        "Nonsense_Mutation",
        "Nonstop_Mutation",
        "In_Frame_Del",
        "In_Frame_Ins",
        "Missense_Mutation",
        "Fusion",
        "Multi_Hit",
        "Multi_Hit_Fusion",
        "Hom_Deletion",
        "Hem_Deletion",
        "amplification",
        "gain",
        "loss"
      )
    )
  
  #### Specify genes
  # Get top mutated this data and goi list
  gene_sum <- mafSummary(maf_object)$gene.summary
  
  # Subset for genes in the histology-specific list
  subset_gene_sum <-
    subset(gene_sum, Hugo_Symbol %in% gene_list$gene)
  
  # Get top altered genes
  goi_ordered <-
    subset_gene_sum[order(subset_gene_sum$AlteredSamples, decreasing = TRUE),]
  
  # Select n top genes
  num_genes <- ifelse(nrow(goi_ordered) > 20, 20, nrow(goi_ordered))
  goi_ordered_num <- goi_ordered[1:num_genes,]
  genes <- goi_ordered_num$Hugo_Symbol
  
  #### Plot Oncoprint
  
  # Given a maf file, plot an oncoprint of the variants in the
  # dataset and save as a png file.
  oncoplot(
    maf_object,
    clinicalFeatures = c(
      "broad_histology",
      "short_histology",
      "reported_gender",
      "tumor_descriptor",
      "molecular_subtype"
    ),
    genes = genes,
    logColBar = TRUE,
    sortByAnnotation = TRUE,
    showTumorSampleBarcodes = TRUE,
    removeNonMutated = TRUE,
    annotationFontSize = 0.7,
    SampleNamefontSize = 0.5,
    fontSize = 0.7,
    colors = color_palette
    ## TODO: Implement use of color palettes in `figures/palettes` directory
  )

}

#### Generate Oncoprints ------------------------------------------------------

primary_plus_oncoprint <- prepare_and_plot_oncoprint(maf_df_primary_plus,
                                                     cnv_file_primary_plus,
                                                     fusion_file_primary_plus,
                                                     gene_list = gene_list)

primary_only_oncoprint <- prepare_and_plot_oncoprint(maf_df_primary_only,
                                                     cnv_file_primary_only,
                                                     fusion_file_primary_only,
                                                     gene_list = gene_list)

#### Assemble multipanel plot -------------------------------------------------

# Combine plots with patchwork
# Layout of the two plots will be one over the other (1 column and 2 rows)
combined_plot <- primary_plus_oncoprint + primary_only_oncoprint +
  plot_layout(ncol = 1, nrow = 2) +
  plot_annotation(tag_levels = 'A') &
  theme(# add uniform labels
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 9))

# Save to PNG
ggplot2::ggsave(file.path(output_dir, "fig3-oncoprint-lanscape.png"),
                width = 12, height = 8,
                units = "in"
)
