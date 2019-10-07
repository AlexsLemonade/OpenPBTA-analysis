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

# Read in metadata
metadata <-
  readr::read_tsv(file.path(input_directory, "pbta-histologies.tsv"))

# Read in cnvkit data
cnvkit <- read_in_cnv(input_directory, "pbta-cnv-cnvkit.seg.gz")

# Read in cnvkit subset data
cnvkit_subset <-
  read_in_cnv(input_directory, "testing/pbta-cnv-cnvkit.seg.gz")

# Read in controlfreec data
controlfreec <-
  read_in_cnv(input_directory, "pbta-cnv-controlfreec.seg.gz")

# Read in controlfreec subset data
controlfreec_subset <-
  read_in_cnv(input_directory, "testing/pbta-cnv-controlfreec.seg.gz")

#### Filter data ---------------------------------------------------------------

# Filter the cnvkit data by cutoff segmean in preparation for plotting
cnvkit_format <- filter_segmean(cnvkit)

# Filter the cnvkit subset data by cutoff segmean in preparation for plotting
cnvkit_subset <- filter_segmean(cnvkit_subset)

# Filter the controlfreec data by cutoff segmean in preparation for plotting
controlfreec_format <- filter_segmean(controlfreec)

# Filter the controlfreec subset data for cutoff in preparation for plotting
controlfreec_subset <- filter_segmean(controlfreec_subset)

#### Plot frequency and proportion using GenVisR -------------------------------

# Run `plot_cnFreq` for cnvkit (has to be run with subset because
# original dataset is too large for function)
cnvkit_prop_plot <-
  plot_cnFreq(cnvkit_format, "proportion", "CNVkit proportion")
cnvkit_freq_plot <-
  plot_cnFreq(cnvkit_format, "frequency", "CNVkit frequency")

# Run `plot_cnFreq` for controlfreec
controlfreec_prop_plot <-
  plot_cnFreq(controlfreec_format,
              "proportion",
              "Control-FREEC proportion")
controlfreec_freq_plot <-
  plot_cnFreq(controlfreec_format,
              "frequency",
              "Control-FREEC frequency")

# Plot cowplot of frequency plots and save
plot_cowplot(
  cnvkit_freq_plot,
  controlfreec_freq_plot,
  output_directory,
  "compare_cnv_output_frequency.pdf"
)

# Plot cowplot of proportion plots and save
plot_cowplot(
  cnvkit_prop_plot,
  controlfreec_prop_plot,
  output_directory,
  "compare_cnv_output_proportion.pdf"
)

#### Plot boxplots using ggplot2 -----------------------------------------------

# Run `plot_boxplot` on cnvkit
cnvkit_boxplot <- plot_boxplot(cnvkit_format, "CNVkit boxplot")

# Run `plot_boxplot` on controlfreec
controlfreec_boxplot <-
  plot_boxplot(controlfreec_format, "Control-FREEC boxplot")

# Save the plot combining the cnvkit and controlfreec boxplots
plot_cowplot(
  cnvkit_boxplot,
  controlfreec_boxplot,
  output_directory,
  "compare_cnv_output_boxplot.pdf"
)

##### Plot barplots using ggplot2 ----------------------------------------------

# Run `plot_histology_barplot` on cnvkit
cnvkit_annotated_barplot_histology <-
  plot_histology_barplot(cnvkit_format, metadata, "CNVkit")

# Run `plot_histology_barplot` on controlfreec
controlfreec_annotated_barplot_histology <-
  plot_histology_barplot(controlfreec_format,
               metadata,
               "Control-FREEC")

# Run `plot_aberration_barplot` on cnvkit
cnvkit_annotated_barplot_aberration <-
  plot_aberration_barplot(cnvkit_format, "CNVkit")

# Run `plot_abberration_barplot` on controlfreec
controlfreec_annotated_barplot_aberration <-
  plot_aberration_barplot(controlfreec_format,
                         "Control-FREEC")

# Save the plot combining the cnvkit and controlfreec barplots
plot_cowplot(
  cnvkit_annotated_barplot_histology,
  controlfreec_annotated_barplot_histology,
  output_directory,
  "compare_cnv_output_barplot_histology.pdf"
)
plot_cowplot(
  cnvkit_annotated_barplot_aberration,
  controlfreec_annotated_barplot_aberration,
  output_directory,
  "compare_cnv_output_barplot_aberration.pdf"
)
