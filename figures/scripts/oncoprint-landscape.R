# Chante Bethell for ALSF CCDL 2021
#
# Makes a multipanel plot for broad histology specific oncoprints.
#
# Code for this script was adapted from 
# `figures/scripts/transcriptomic-overview.R`

#### Set Up --------------------------------------------------------------------

# Load maftools
library(maftools)

# Load multipanelfigure for figure assembly later
library(multipanelfigure)

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

#### Directories and Files -----------------------------------------------------

# Detect the ".git" folder -- this will in the project root directory.
# Use this as the root directory to ensure proper execution, no matter where
# it is called from.
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Declare output directory
output_dir <- file.path(root_dir, "figures", "pngs")

# Final figure filename
output_png <- file.path(output_dir, "oncoprint-landscape.png")

# Data directory
data_dir <- file.path(root_dir, "data")

# Analysis directory
analyses_dir <- file.path(root_dir, "analyses")

# Oncoprint scratch files
scratch_dir <- file.path(root_dir, "scratch", "oncoprint_files")

# Source the custom functions script
source(
  file.path(
    analyses_dir,
    "oncoprint-landscape",
    "util",
    "oncoplot-functions.R"
  )
)

# Define cnv_file, fusion_file, and genes object here as they still need to
# be defined for the `prepare_and_plot_oncoprint` custom function (the
# cnv_file specifically for the `read.maf` function within the custom function),
# even if they are NULL
cnv_df <- file.path(scratch_dir, "primary_only_cnv.tsv")
fusion_df <- file.path(scratch_dir, "primary_only_fusions.tsv")
maf_file <- file.path(scratch_dir, "primary_only_maf.tsv")
metadata_file <- file.path(data_dir, "pbta-histologies.tsv")

#### Temporary PNG names -------------------------------------------------------
# Because we are mixing Heatmaps and ggplots in this figure, one of the best
# ways to control sizes is to create PNGs and then use multipanelfigure to
# put the final figure together.
#
# Here we'll set up the file names used throughout
lgat_png <- file.path(output_dir, "temp_lgat.png")
embryonal_png <- file.path(output_dir, "temp_embryonal.png")
hgat_png <- file.path(output_dir, "temp_hgat.png")
ependymal_png <- file.path(output_dir, "temp_ependymal.png")
other_cns_png <- file.path(output_dir, "temp_other_cns.png")

#### Functions ----------------------------------------------------------------

plot_oncoprint <- function (broad_histology,
                            metadata,
                            cnv_df,
                            fusion_df,
                            maf_df,
                            png_name) {
  # This function takes in the approprate data for plotting oncoprints, as well
  # as a specific broad histology string, to prepare and plot an oncoprint
  
  # Read in the histology specific genes of interest list
  goi_list <-
    readr::read_tsv(tolower(gsub(
      " ",
      "-",
      file.path(
        analyses_dir,
        "oncoprint-landscape",
        "data",
        paste0(broad_histology, "_goi_list.tsv")
      )
    ))) %>%
    as.matrix()
  
  # Filter to the metadata associated with the broad histology value
  if (broad_histology != "Other CNS") {
    metadata <- metadata %>%
      dplyr::filter(broad_histology == as.vector(broad_histology))
  } else {
    metadata <- metadata %>%
      dplyr::filter(
        broad_histology %in% c(
          "Tumors of sellar region",
          "Neuronal and mixed neuronal-glial tumor",
          "Tumor of cranial and paraspinal nerves",
          "Meningioma",
          "Mesenchymal non-meningothelial tumor",
          "Germ cell tumor",
          "Choroid plexus tumor",
          "Histiocytic tumor",
          "Tumor of pineal region",
          "Metastatic tumors",
          "Other astrocytic tumor",
          "Lymphoma",
          "Melanocytic tumor",
          "Other tumor"
        )
      )
  }
  
  # Now filter the remaining data files
  maf_df <- maf_df %>%
    dplyr::filter(Tumor_Sample_Barcode %in% metadata$Tumor_Sample_Barcode)
  
  cnv_df <- cnv_df %>%
    # as.data.frame() %>%
    dplyr::filter(Tumor_Sample_Barcode %in% metadata$Tumor_Sample_Barcode)
  
  fusion_df <- fusion_df %>%
    dplyr::filter(Tumor_Sample_Barcode %in% metadata$Tumor_Sample_Barcode)
  
  # Read in the oncoprint color palette
  oncoprint_col_palette <- readr::read_tsv(file.path(
    root_dir,
    "figures",
    "palettes",
    "oncoprint_color_palette.tsv"
  )) %>%
    # Use deframe so we can use it as a recoding list
    tibble::deframe()
  
  # Color coding for `display_group` classification
  # Get unique tumor descriptor categories
  histologies_color_key_df <- metadata %>%
    dplyr::arrange(display_order) %>%
    dplyr::select(display_group, hex_codes) %>%
    dplyr::distinct()
  
  # Make color key specific to these samples
  histologies_color_key <- histologies_color_key_df$hex_codes
  names(histologies_color_key) <-
    histologies_color_key_df$display_group
  
  # Now format the color key objet into a list
  annotation_colors <- list(display_group = histologies_color_key)
  
  #### Prepare MAF object for plotting ------------------------------------------
  
  maf_object <- prepare_maf_object(
    maf_df = maf_df,
    cnv_df = cnv_df,
    metadata = metadata,
    fusion_df = fusion_df
  )
  
  #### Plot and Save Oncoprint --------------------------------------------------
  
  # Given a maf object, plot an oncoprint of the variants in the
  # dataset and save as a png file.
  png(png_name,
      width = 65,
      height = 30,
      units = "cm",
      res = 300
  )
  oncoplot(
    maf_object,
    clinicalFeatures = "display_group",
    genes = goi_list,
    logColBar = TRUE,
    sortByAnnotation = TRUE,
    showTumorSampleBarcodes = TRUE,
    removeNonMutated = TRUE,
    annotationFontSize = 1.0,
    SampleNamefontSize = 0.5,
    fontSize = 0.7,
    colors = oncoprint_col_palette,
    annotationColor = annotation_colors,
    bgCol = "#F5F5F5",
    top = 25
  )
  dev.off()
  
}

#### Read in data --------------------------------------------------------------

# Read in metadata
metadata <- readr::read_tsv(metadata_file, guess_max = 10000) %>%
  dplyr::rename(Tumor_Sample_Barcode = sample_id)

# Read in MAF file
maf_df <- data.table::fread(maf_file,
                            stringsAsFactors = FALSE,
                            data.table = FALSE)

# Read in cnv file
cnv_df <- readr::read_tsv(cnv_df) %>%
  dplyr::mutate(
    Variant_Classification = dplyr::case_when(
      Variant_Classification == "loss" ~ "Del",
      Variant_Classification == "gain" ~ "Amp",
      TRUE ~ as.character(Variant_Classification)
    )
  )

# Read in fusion file and join
fusion_df <- readr::read_tsv(fusion_df)

#### Set up oncoprint annotation objects --------------------------------------
# Read in histology standard color palette for project
histology_label_mapping <- readr::read_tsv(
  file.path(root_dir,
            "figures",
            "palettes", 
            "histology_label_color_table.tsv")) %>%
  # Select just the columns we will need for plotting
  dplyr::select(Kids_First_Biospecimen_ID, display_group, display_order, hex_codes)

# Join on these columns to the metadata
metadata <- metadata %>% 
  dplyr::inner_join(histology_label_mapping, by = "Kids_First_Biospecimen_ID") %>% 
  # Reorder display_group based on display_order
  dplyr::mutate(display_group = forcats::fct_reorder(display_group, display_order))

#### Run plot oncoprint function ----------------------------------------------

lgat_onco <- plot_oncoprint("Low-grade astrocytic tumor",
                            metadata,
                            cnv_df,
                            fusion_df,
                            maf_df,
                            lgat_png)

embryonal_onco <- plot_oncoprint("Embryonal tumor",
                            metadata,
                            cnv_df,
                            fusion_df,
                            maf_df,
                            embryonal_png)

hgat_onco <- plot_oncoprint("Diffuse astrocytic and oligodendroglial tumor",
                            metadata,
                            cnv_df,
                            fusion_df,
                            maf_df,
                            hgat_png)

ependymal_onco <- plot_oncoprint("Ependymal tumor",
                            metadata,
                            cnv_df,
                            fusion_df,
                            maf_df,
                            ependymal_png)

other_cns_onco <- plot_oncoprint("Other CNS",
                            metadata,
                            cnv_df,
                            fusion_df,
                            maf_df,
                            other_cns_png)


#### Assemble multipanel figure ------------------------------------------------

transcriptomic_figure <- multi_panel_figure(columns = 6,
                                            rows = 2,
                                            width = 1200,
                                            height = 300,
                                            panel_label_type = "none")

transcriptomic_figure <- fill_panel(transcriptomic_figure,
                                    lgat_png,
                                    col = 1:2,
                                    row = 1,
                                    scaling = "fit")

transcriptomic_figure <- fill_panel(transcriptomic_figure,
                                    embryonal_png,
                                    col = 3:4,
                                    row =1,
                                    scaling = "fit")

transcriptomic_figure <- fill_panel(transcriptomic_figure,
                                    hgat_png,
                                    col = 5:6,
                                    row = 1,
                                    scaling = "fit")

transcriptomic_figure <- fill_panel(transcriptomic_figure,
                                    ependymal_png,
                                    col = 1:2,
                                    row = 2,
                                    scaling = "fit")

transcriptomic_figure <- fill_panel(transcriptomic_figure,
                                    other_cns_png,
                                    col = 3:4,
                                    row = 2,
                                    scaling = "fit")

save_multi_panel_figure(transcriptomic_figure, output_png)

#### Remove temporary PNGs -----------------------------------------------------

# Now that we've constructed and saved the final output we can remove the PNG
# files for the individual panels and the legend
file.remove(c(lgat_png,
              embryonal_png,
              hgat_png,
              ependymal_png,
              other_cns_png))

