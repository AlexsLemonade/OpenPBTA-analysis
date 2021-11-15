# JN Taroni for ALSF CCDL 2021
#
# Makes publication ready OncoPrint panels (PDFs)

#### Directories ---------------------------------------------------------------

# Detect the ".git" folder -- this will in the project root directory.
# Use this as the root directory to ensure proper execution, no matter where
# it is called from.
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Declare output directory
figure_dir <- file.path(root_dir, "figures", "pdfs")

# The panels for the main figure
main_output_dir <- file.path(figure_dir, "fig2", "panels")
dir.create(main_output_dir, recursive = TRUE, showWarnings = FALSE)

# The supplemental figure (this will be the Other CNS cancer groups without
# their own )
supp_output_dir <- file.path(figure_dir, "supp")
dir.create(supp_output_dir, recursive = TRUE, showWarnings = FALSE)

# Oncoprint directory
module_dir <- file.path(root_dir, "analyses", "oncoprint-landscape")

# Intermediate files from running the module itself
data_input_dir <- file.path(root_dir, "scratch", "oncoprint_files")

# We need two of our standardized palettes
palette_dir <- file.path(root_dir, "figures", "palettes")

#### Libraries & custom functions ----------------------------------------------

library(maftools)
library(magrittr)

# Source the custom functions script from the oncoprint module
source(
  file.path(
    module_dir,
    "util",
    "oncoplot-functions.R"
  )
)

# Not to be used outside of this context!!
# Given a file that contains genes of interest, filter to the top 20 genes
# We need to do this once per oncoprint for a total of 4x
get_histology_goi <- function(filename) {

  # Read in genes of interest information using the `read_tsv()` function
  goi_list <- readr::read_tsv(filename) %>%
    as.matrix()

  filtered_maf_object <- subsetMaf(
    maf = histology_maf_object,
    tsb = histologies_df$Tumor_Sample_Barcode,
    genes = goi_list,
    mafObj = TRUE
  )

  # Get top mutated genes per this subset object
  gene_sum <- mafSummary(filtered_maf_object)$gene.summary

  # Sort to get top altered genes rather than mutated only genes
  goi_list <- gene_sum %>%
    dplyr::arrange(dplyr::desc(AlteredSamples)) %>%
    # Filter to genes where multiple samples have an alteration
    dplyr::filter(AlteredSamples > 1) %>%
    dplyr::pull(Hugo_Symbol)

  # Now let's filter to the top 20 genes -- this matches what's done in
  # the module
  goi_list <- goi_list[1:20]

}

# Also not to be used outside of this context!
# Given a vector of cancer groups, prep a MAF object to be used with
# oncoplot() -- filters & combines the SNV, CNV, and fusion data
prep_histology_maf <- function(included_cancer_groups) {

  # Get the sample ids to be included
  included_sample_ids <- histologies_df %>%
    dplyr::filter(cancer_group %in% included_cancer_groups) %>%
    dplyr::pull(Tumor_Sample_Barcode)

  # MAF
  histology_maf_df <- maf_df %>%
    dplyr::filter(Tumor_Sample_Barcode %in% included_sample_ids)

  # CNV
  histology_cnv_df <- cnv_df %>%
    dplyr::filter(Tumor_Sample_Barcode %in% included_sample_ids)

  # Fusions
  histology_fusion_df <- fusion_df %>%
    dplyr::filter(Tumor_Sample_Barcode %in% included_sample_ids)

  # MAF object
  histology_maf_object <- prepare_maf_object(
    maf_df = histology_maf_df,
    cnv_df = histology_cnv_df,
    metadata = histologies_df,
    fusion_df = histology_fusion_df
  )
}

# Save PDF versions of legends when given a vector of labels (group_vector),
# a vector of hexcodes and the output file name
create_legend <- function(group_vector,
                          hexcode_vector,
                          output_pdf,
                          label = "Cancer Group") {
  pdf(output_pdf)
  plot(NULL, xaxt = "n", yaxt = "n", bty = "n", ylab = "", xlab = "",
       xlim = 0:1, ylim = 0:1)
  legend("topleft",
         legend = group_vector,
         col = hexcode_vector,
         pch = 15, pt.cex = 1.5, cex = 0.75, bty = "n")
  mtext(label, at = 0.0625, cex = 1)
  dev.off()
}

#### Required data files -------------------------------------------------------

# File suffixes for different data types, from the module
maf_suffix <- "_maf.tsv"
cnv_suffix <- "_cnv.tsv"
fusion_suffix <- "_fusions.tsv"

# Tiny function to make both lists of files
create_data_input_list <- function(prefix) {
  list(
    maf = file.path(data_input_dir, paste0(prefix, maf_suffix)),
    cnv = file.path(data_input_dir, paste0(prefix, cnv_suffix)),
    fusion = file.path(data_input_dir, paste0(prefix, fusion_suffix))
  )
}

# Create the data input file lists for both primary and primary-plus
data_input_list <- list(
  primary_only = create_data_input_list("primary_only"),
  primary_plus = create_data_input_list("primary-plus")
)

#### Metadata ------------------------------------------------------------------

histologies_df <- readr::read_tsv(
  file.path(root_dir, "data", "pbta-histologies.tsv")
) %>%
  dplyr::rename(Tumor_Sample_Barcode = sample_id)

#### Palette files -------------------------------------------------------------

# For coloring the cells of the oncoprint -- the data frame is for the legend
oncoprint_palette_df <- readr::read_tsv(
  file.path(palette_dir, "oncoprint_color_palette.tsv")
)

# For use with oncoplot()
oncoprint_palette <- oncoprint_palette_df %>%
  # Use deframe so we can use it as a recoding list
  tibble::deframe()

# Generate palette for disease labels, etc.
group_palette_df <-  readr::read_tsv(
  file.path(palette_dir, "broad_histology_cancer_group_palette.tsv")
)

# Get palette for cancer group that is *specifically* for the oncoprint
cancer_group_palette <- group_palette_df %>%
  dplyr::select(cancer_group, oncoprint_hex) %>%
  # Remove NA values
  dplyr::filter(complete.cases(.))

# Make color palette suitable for use later
cancer_group_colors <- cancer_group_palette$oncoprint_hex
names(cancer_group_colors) <- cancer_group_palette$cancer_group

# Now format the color key object into a list
annotation_colors <- list(cancer_group = cancer_group_colors,
                          germline_sex_estimate = c("Male" = "#2166ac",
                                                    "Female" = "#b2182b"))

#### Hardcoding legend ordering ------------------------------------------------
# These reflect the cancer groups that are included, and their ordering, for
# the primary plus oncoprints

legend_ordering <- list(
  lgat = c(
    "Low-grade glioma astrocytoma",
    "Ganglioglioma",
    "Pleomorphic xanthoastrocytoma"
  ),
  hgat = c(
    "Diffuse midline glioma",
    "High-grade glioma astrocytoma",
    "Oligodendroglioma"
  ),
  embryonal = c(
    "Medulloblastoma",
    "Atypical Teratoid Rhabdoid Tumor",
    "Embryonal tumor with multilayer rosettes",
    "CNS Embryonal tumor",
    "Ganglioneuroblastoma"
  ),
  other = c(
    "Ependymoma",
    "Craniopharyngioma",
    "Meningioma",
    "Dysembryoplastic neuroepithelial tumor",
    "Ewing sarcoma",
    "Schwannoma",
    "Neurofibroma Plexiform"
  )

)

#### Genes of interest lists ---------------------------------------------------

goi_dir <- file.path(module_dir, "data")
goi_files_list <- list(
  "Low-grade astrocytic tumor" = list(
    file = file.path(goi_dir, "lgat_goi_list.tsv"),
    shorthand = "lgat"
  ),
  "Diffuse astrocytic and oligodendroglial tumor" = list(
    file = file.path(goi_dir, "hgat_goi_list.tsv"),
    shorthand = "hgat"
  ),
  "Embryonal tumor" = list(
    file = file.path(goi_dir, "embryonal-tumor_goi_list.tsv"),
    shorthand = "embryonal"
  ),
  "Other CNS" = list(
    file = file.path(goi_dir, "other_goi_list.tsv"),
    shorthand = "other"
  )
)

#### OncoPrints and cancer group legends ---------------------------------------

for (type_iter in seq_along(data_input_list)) {

  # primary or primary plus?
  specimen_type <- names(data_input_list)[type_iter]

  # Read in the input data for the relevant set of specimens
  maf_df <- data.table::fread(data_input_list[[specimen_type]]$maf,
                              stringsAsFactors = FALSE,
                              data.table = FALSE) %>%
    dplyr::filter(Variant_Classification != "Intron")
  cnv_df <- readr::read_tsv(data_input_list[[specimen_type]]$cnv)
  fusion_df <- readr::read_tsv(data_input_list[[specimen_type]]$fusion)

  for (histology in unique(group_palette_df$oncoprint_group)) {

    # There are NA values in this column, so if we encounter that skip to the
    # next one
    if (is.na(histology)) next

    # Get vector of cancer groups to include
    included_cancer_groups <- group_palette_df %>%
      dplyr::filter(oncoprint_include,
                    oncoprint_group == histology) %>%
      dplyr::pull(cancer_group)

    # For convenience, save the shorthand (e.g., "lgat") for filenames, etc.
    histology_shorthand <- goi_files_list[[histology]]$shorthand

    # Prep MAF object for plot
    histology_maf_object <- prep_histology_maf(included_cancer_groups)

    # Prepare the genes of interest list for this histology
    histology_goi <- get_histology_goi(goi_files_list[[histology]]$file)

    # Construct the output PDF name
    output_pdf <- paste(specimen_type,
                        histology_shorthand,
                        "oncoprint.pdf",
                        sep = "_")

    pdf(
      file.path(main_output_dir, output_pdf),
      width = 13.75,
      height = 7.9
    )

    oncoplot(
      histology_maf_object,
      clinicalFeatures = c("cancer_group", "germline_sex_estimate"),
      genes = histology_goi,
      logColBar = TRUE,
      sortByAnnotation = TRUE,
      showTumorSampleBarcodes = FALSE,
      removeNonMutated = TRUE,
      annotationFontSize = 1.25,
      SampleNamefontSize = 1,
      legend_height = 0.1,
      fontSize = 1,
      colors = oncoprint_palette,
      annotationColor = annotation_colors,
      bgCol = "#F5F5F5",
      drawRowBar = FALSE,
      titleText = histology,
      titleFontSize = 1.3,
      gene_mar = 10
    )

    dev.off()

    #### Save the cancer group legend ####

    # Construct the legend output PDF filename
    legend_output_pdf <- file.path(
      main_output_dir,
      paste0(histology_shorthand,
             "_oncoprint_cancer_group_legend.pdf")
    )

    # This legend should only be for the relevant histology (i.e., panel) and
    # only include cancer groups that have mutated samples -- this is captured
    # in the legend ordering list
    legend_df <- group_palette_df %>%
      dplyr::filter(oncoprint_include,
                    oncoprint_group == histology,
                    cancer_group %in% legend_ordering[[histology_shorthand]])

    # Order according to the legend ordering list such that the cancer groups
    # are in the same order in the legend and the oncoprint
    legend_df <- legend_df[match(legend_ordering[[histology_shorthand]],
                                 legend_df$cancer_group), ]

    # Save the legend as PDF
    create_legend(legend_df$cancer_group,
                  legend_df$oncoprint_hex,
                  legend_output_pdf)

    #### Supplemental display item ####
    #### Other CNS only ###############

    if (histology == "Other CNS") {

      # Get vector of cancer groups to include
      # These will be Other CNS oncoprint groups with a FALSE oncoprint_include
      included_cancer_groups <- group_palette_df %>%
        dplyr::filter(!oncoprint_include,
                      oncoprint_group == histology) %>%
        dplyr::pull(cancer_group)

      # Prep MAF object for plot
      histology_maf_object <- prep_histology_maf(included_cancer_groups)

      # We need a new palette for the cancer groups but we only show samples
      # with alterations in our genes of interest list
      mutated_samples <- histology_maf_object@data %>%
        dplyr::filter(Hugo_Symbol %in% histology_goi) %>%
        dplyr::pull(Tumor_Sample_Barcode) %>%
        unique()

      # Only the cancer groups that will be represented in the oncoprint
      mutated_cancer_groups <- histologies_df %>%
        dplyr::filter(Tumor_Sample_Barcode %in% mutated_samples,
                      !is.na(cancer_group)) %>%
        dplyr::pull(cancer_group) %>%
        unique()

      # Use the paired color brewer palette for the cancer groups that will
      # be included
      supp_colors <- RColorBrewer::brewer.pal(length(mutated_cancer_groups),
                                              "Paired")
      names(supp_colors) <- mutated_cancer_groups

      # Construct the output PDF name
      output_pdf <- paste(specimen_type,
                          goi_files_list[[histology]]$shorthand,
                          "supplemental",
                          "oncoprint.pdf",
                          sep = "_")

      pdf(
        file.path(supp_output_dir, output_pdf),
        width = 13.75,
        height = 7.9
      )

      oncoplot(
        histology_maf_object,
        clinicalFeatures = c("cancer_group", "germline_sex_estimate"),
        genes = histology_goi,
        logColBar = TRUE,
        sortByAnnotation = TRUE,
        showTumorSampleBarcodes = FALSE,
        removeNonMutated = TRUE,
        annotationFontSize = 1.25,
        SampleNamefontSize = 1,
        annotationColor = list(cancer_group = supp_colors,
                               germline_sex_estimate = c("Male" = "#2166ac",
                                                         "Female" = "#b2182b")),
        fontSize = 1,
        colors = oncoprint_palette,
        bgCol = "#F5F5F5",
        drawRowBar = FALSE,
        titleText = histology,
        titleFontSize = 1.3,
        gene_mar = 10
      )

      dev.off()

    }

  }

}

#### Other legends -------------------------------------------------------------

# Germline sex estimate
create_legend(names(annotation_colors$germline_sex_estimate),
              annotation_colors$germline_sex_estimate,
              file.path(main_output_dir, "germline_sex_legend.pdf"),
              "Germline Sex Estimate")

# Oncoprint alterations
oncoprint_palette_df <- oncoprint_palette_df %>%
  dplyr::mutate(color_names = stringr::str_replace_all(color_names,
                                                       "_", " ") %>%
                  stringr::str_replace_all("'", "' ") %>%
                  stringr::str_replace(" Hit", "-Hit"))

create_legend(oncoprint_palette_df$color_names,
              oncoprint_palette_df$hex_codes,
              file.path(main_output_dir, "oncoprint_alteration_legend.pdf"),
              "Alteration")
