# Plot and compare detected CNV aberrations given CNVkit and Control-FREEC 
# output
#
# Chante Bethell for CCDL 2019
#
# Usage:
# This script is intended to be run via the command line from the top directory 
# of the repository as follows:
#
# Rscript analyses/cnv-comparison/01-cnv-comparison-plotting.R
#

#### Install packages ----------------------------------------------------------
if (!("GenVisR" %in% installed.packages())) {
  install.packages("BiocManager")
  BiocManager::install("GenVisR")
}

if (!("cowplot" %in% installed.packages())) {
  install.packages("BiocManager")
  BiocManager::install("cowplot")
}

##### Set up functions ---------------------------------------------------------

# Magrittr pipe
`%>%` <- dplyr::`%>%`

# Detect the ".git" folder -- this will in the project root directory.
# Use this as the root directory to ensure proper sourcing of functions no
# matter where this is called from
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Source custom functions script
source(file.path(root_dir, "analyses", "cnv-comparison", "util",
                 "cnv-comparison-functions.R"))

#### Set up file paths ---------------------------------------------------------

input_directory <- file.path(root_dir, "data")

# Create the output directory if it does not exist
output_directory <- file.path(root_dir, "analyses", "cnv-comparison", "plots")

if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)
}

#### Read in data --------------------- ----------------------------------------

# Create list of file paths to the CNV data 
cnv_list <-
  list(
    file.path(input_directory, "pbta-cnv-cnvkit.seg.gz"),
    file.path(input_directory, "pbta-cnv-controlfreec.seg.gz")
  )

# Read in list of CNV data
cnv_data <- lapply(cnv_list, read_in_cnv)

# Read in metadata
metadata <-
  readr::read_tsv(file.path(input_directory, "pbta-histologies.tsv"))

#### Filter data ---------------------------------------------------------------

# Filter CNV data by cutoff segmean using `filter_segmean`
cnv_filtered <-
  lapply(cnv_data, filter_segmean, segmean_cutoff = 0.5)

#### Plot frequency and proportion using GenVisR -------------------------------

# Run `plot_cnFreq` 
cnv_proportion_plot <- lapply(cnv_filtered, plot_cn_freq, 0, .2, "proportion")
cnv_frequency_plot <- lapply(cnv_filtered, plot_cn_freq, 0, .2, "frequency")

# Plot cowplot of frequency plots and save
plot_cowplot(
  cnv_proportion_plot,
  output_directory,
  "compare_cnv_output_proportion.pdf"
)

# Plot cowplot of proportion plots and save
plot_cowplot(
  cnv_frequency_plot,
  output_directory,
  "compare_cnv_output_frequency.pdf"
)

#### Plot boxplots using ggplot2 -----------------------------------------------

# Run `plot_violin` on CNV data
cnv_violin_plots <- lapply(cnv_filtered, plot_violin)

# Save the plot combining the cnvkit and controlfreec boxplots
plot_cowplot(
  cnv_violin_plots,
  output_directory,
  "compare_cnv_output_violin_plot.pdf"
)

##### Plot barplots using ggplot2 ----------------------------------------------

# Run `plot_histology_barplot`
cnv_histology_barplots <-
  lapply(cnv_filtered, plot_histology_barplot, metadata)

# Run `plot_aberration_barplot`
cnv_aberration_barplots <-
  lapply(cnv_filtered, plot_aberration_barplot)

# Save the plot combining the cnvkit and controlfreec barplots
plot_cowplot(
  cnv_histology_barplots,
  output_directory,
  "compare_cnv_output_barplot_histology.pdf"
)
plot_cowplot(
  cnv_aberration_barplots,
  output_directory,
  "compare_cnv_output_barplot_aberration.pdf"
)