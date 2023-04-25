# JN Taroni and SJ Spielman for ALSF CCDL 2021-2022
#
# Makes publication ready OncoPrint panels, specifically for primary-only (PDFs)

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
supp_output_dir <- file.path(figure_dir, "supp", "figs3", "panels")
dir.create(supp_output_dir, recursive = TRUE, showWarnings = FALSE)

# Oncoprint directory
module_dir <- file.path(root_dir, "analyses", "oncoprint-landscape")

# Intermediate files from running the module itself
data_input_dir <- file.path(root_dir, "scratch", "oncoprint_files")

# We need two of our standardized palettes
palette_dir <- file.path(root_dir, "figures", "palettes")

# Directory where CSV files for Zenodo upload will be saved
# Note that only the S3b filename is defined here, as fig2 file names
#  are defined making use of code later in this script
zenodo_upload_dir <- file.path(root_dir, "tables", "zenodo-upload")
figS3b_csv <- file.path(zenodo_upload_dir, "figure-S3b-data.csv")


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
# oncoplot() -- filters & combines the SNV, CNV, and fusion data.
# This function also prepares and returns a data frame that can be exported containing
# all oncoplot data.
prep_histology_maf <- function(included_cancer_groups,
                               main = TRUE) {
  # The `main` argument indicates whether this is main text figure (TRUE),
  #  or a supplemental figure (FALSE). Main and supp use different cancer group
  #  columns.

  # Get the sample ids to be included
  if (main) {
    included_sample_ids <- histologies_df %>%
      dplyr::filter(cancer_group_display %in% included_cancer_groups) %>%
      dplyr::pull(Tumor_Sample_Barcode)

    # Update histologies_df for palette compatibility
    histologies_df_temp <- histologies_df %>%
      dplyr::select(-cancer_group) %>%
      dplyr::rename(cancer_group = cancer_group_display)

  } else {
    included_sample_ids <- histologies_df %>%
      dplyr::filter(cancer_group %in% included_cancer_groups) %>%
      dplyr::pull(Tumor_Sample_Barcode)

    # Create the temp variable as is
    histologies_df_temp <- histologies_df
  }

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
    metadata = histologies_df_temp,
    fusion_df = histology_fusion_df
  )

  # Create a data frame that can be used to export this information
  select_cols <- c("Hugo_Symbol", "Tumor_Sample_Barcode", "Variant_Classification")

  maf_export_df <- histology_maf_df %>%
    dplyr::select(select_cols) %>%
    dplyr::bind_rows(
      dplyr::select(histology_fusion_df, select_cols)
    ) %>%
    dplyr::bind_rows(
      dplyr::mutate(histology_cnv_df, Variant_Classification = NA_character_)
    )  %>%
    # Join with metadata information that is shown in the plot
    dplyr::inner_join(
      dplyr::select(histologies_df,
                    Tumor_Sample_Barcode,
                    cancer_group_display,
                    broad_histology_display,
                    germline_sex_estimate)
    ) %>%
    # reorder columns and rename for consistency with plot and data release histologies file
    dplyr::select(sample_id = Tumor_Sample_Barcode,
                  Hugo_Symbol,
                  alteration = Variant_Classification,
                  dplyr::everything()) %>%
    # arrange on sample_id, within groups of diagnoses
    dplyr::group_by(cancer_group_display) %>%
    dplyr::arrange(sample_id, .by_group = TRUE)

  # Return the maf object and the df for export
  return(
    list(
      maf_object = histology_maf_object,
      maf_df  = maf_export_df
    )
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

# Ensure data_input_dir exists
if (!(dir.exists(data_input_dir))) {
  stop("Error: The oncoprint-landscape module must be run first to create and populate the directory `scratch/oncoprint_files/` in order to create these figures.")
}

# File suffixes for different data types, from the module
maf_suffix <- "_maf.tsv"
cnv_suffix <- "_cnv.tsv"
fusion_suffix <- "_fusions.tsv"

# Tiny function to make lists of files
create_data_input_list <- function(prefix) {
  list(
    maf = file.path(data_input_dir, paste0(prefix, maf_suffix)),
    cnv = file.path(data_input_dir, paste0(prefix, cnv_suffix)),
    fusion = file.path(data_input_dir, paste0(prefix, fusion_suffix))
  )
}

# Create the data input file lists for both primary and primary-plus
# Primary-plus is not currently used in main text, so commented out
data_input_list <- list(
  primary_only = create_data_input_list("primary_only") #,
  #primary_plus = create_data_input_list("primary-plus")
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
  ) %>%
  # Remove NAs from oncoprint_group
  tidyr::drop_na(oncoprint_group)

# Add cancer_group_display into histologies_df
histologies_df <- dplyr::inner_join(
  histologies_df,
  dplyr::select(group_palette_df, broad_histology, broad_histology_display, cancer_group, cancer_group_display)
)

# Get palette for cancer group that is *specifically* for the oncoprint
cancer_group_palette <- group_palette_df %>%
  dplyr::select(cancer_group_display, cancer_group_hex) %>%
  # Remove NA values
  dplyr::filter(complete.cases(.)) %>%
  dplyr::distinct()

# Make color palette suitable for use later
cancer_group_colors <- cancer_group_palette$cancer_group_hex
names(cancer_group_colors) <- cancer_group_palette$cancer_group_display

# Define array of colors for sex estimates
germline_sex_estimate_colors <- c("Male"   = "#2166ac",
                                  "Female" = "#b2182b")

# Now format the color key object into a list
annotation_colors <- list(cancer_group = cancer_group_colors,
                          germline_sex_estimate = germline_sex_estimate_colors)

#### Hardcoding legend ordering ------------------------------------------------
# These reflect the cancer groups that are included, and their ordering, for
# the primary only oncoprints

legend_ordering <- list(
  lgat = c(
    "Pilocytic astrocytoma",
    "Other low-grade glioma",
    "Ganglioglioma",
    "Pleomorphic xanthoastrocytoma",
    "Subependymal Giant Cell Astrocytoma"
  ),
  hgat = c(
    "Diffuse midline glioma",
    "Other high-grade glioma"
  ),
  embryonal = c(
    "Medulloblastoma",
    "Other embryonal tumor",
    "Atypical Teratoid Rhabdoid Tumor"
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

#### Define file names for Zenodo CSV export.
# These files are being defined here so code can make use of `goi_files_list` names,
#  which is looped over to create the oncoplots
zenodo_csv_filenames <- file.path(
  zenodo_upload_dir,
  glue::glue("figure-2{letters[1:4]}-data.csv")
) %>%
  # Specify order of manuscript figure 2:
  # https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/ae3eb012df4a5df26ee81fbd9dcc0a9ffbe12446/figures/pngs/figure2.png
  set_names( names(goi_files_list)[c(1, 3, 2, 4)] )


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

  for (histology in names(goi_files_list)) {

    # For convenience, save the shorthand (e.g., "lgat") for filenames, etc.
    histology_shorthand <- goi_files_list[[histology]]$shorthand

    # Get vector of cancer groups to include from hard-coded legend order
    included_cancer_groups <- legend_ordering[[histology_shorthand]]

    # Prep MAF object for plot and data create CSV for export
    histology_maf_pair <- prep_histology_maf(included_cancer_groups)
    histology_maf_object <- histology_maf_pair$maf_object # MUST define since used by `get_histology_goi()`

    # Prepare the genes of interest list for this histology
    histology_goi <- get_histology_goi(goi_files_list[[histology]]$file)
    
    # Construct the output PDF name
    output_pdf <- paste(specimen_type,
                        histology_shorthand,
                        "oncoprint.pdf",
                        sep = "_")

    pdf(
      file.path(main_output_dir, output_pdf),
      width = 9,
      height = 5.17
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
      # largest size before title is outside of margins
      titleFontSize = 1.5, 
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
      dplyr::filter(oncoprint_group == histology,
                    cancer_group_display %in% legend_ordering[[histology_shorthand]])

    # Order according to the legend ordering list such that the cancer groups
    # are in the same order in the legend and the oncoprint
    legend_df <- legend_df[match(legend_ordering[[histology_shorthand]],
                                 legend_df$cancer_group_display), ]

    # Save the legend as PDF
    create_legend(legend_df$cancer_group_display,
                  legend_df$cancer_group_hex,
                  legend_output_pdf)

    # Export CSV for Zenodo upload, after filtering to histology_goi genes
    maf_export <- histology_maf_pair$maf_df %>%
      dplyr::filter(Hugo_Symbol %in% histology_goi)
    readr::write_csv(
      maf_export,
      zenodo_csv_filenames[[histology]]
    )

    #### Supplemental display item ####
    #### Other CNS only ###############

    if (histology == "Other CNS") {

      # Get vector of cancer groups to include
      # These will be Other CNS oncoprint groups with a FALSE oncoprint_main
      included_cancer_groups <- group_palette_df %>%
        dplyr::filter(oncoprint_group == histology,
                      oncoprint_main == FALSE) %>%
        dplyr::pull(cancer_group)

      # Prep MAF object for plot
      histology_maf_pair <- prep_histology_maf(included_cancer_groups,
                                                 main = FALSE)
      histology_maf_object <- histology_maf_pair$maf_object

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
      names(supp_colors) <- sort(mutated_cancer_groups) # oncoprint will only do alphabetical

      # Have samples go in order of these mutated_cancer_groups
      sample_order <- histologies_df %>%
        dplyr::filter(cancer_group %in% mutated_cancer_groups) %>%
        dplyr::distinct() %>%
        dplyr::arrange(cancer_group) %>%
        dplyr::pull(Tumor_Sample_Barcode)

      # Construct the output PDF name
      output_pdf <- paste(specimen_type,
                          goi_files_list[[histology]]$shorthand,
                          "supplemental",
                          "oncoprint.pdf",
                          sep = "_")

      pdf(
        file.path(supp_output_dir, output_pdf),
        width = 18, # wider to fit the full legend
        height = 10.35 # height to match increased width so this has the same aspect ratio as other oncoplots
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
                               germline_sex_estimate = germline_sex_estimate_colors),
        fontSize = 1,
        colors = oncoprint_palette,
        bgCol = "#F5F5F5",
        drawRowBar = FALSE,
        titleText = histology,
        titleFontSize = 1.3,
        gene_mar = 10,
        sampleOrder = sample_order
      )

      dev.off()

      # Export CSV for Zenodo upload
      # First, filter to GOI:
      maf_export <- histology_maf_pair$maf_df %>%
        dplyr::filter(Hugo_Symbol %in% histology_goi)

      readr::write_csv(
        maf_export,
        figS3b_csv
      )

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
