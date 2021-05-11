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

# Oncoprint directory
onco_dir <- file.path(root_dir, "analyses", "oncoprint-landscape", "plots")

#### Command line options ------------------------------------------------------

# Declare command line options
option_list <- list(
  optparse::make_option(
    c("-l", "--lead_filename"),
    type = "character",
    default = NULL,
    help = "the leading filename for the oncoprints -- can be primary_only or primary-plus"
  ),
  optparse::make_option(
    c("-p", "--png_name"),
    type = "character",
    default = NULL,
    help = "oncoprint output png file name"
  )
)

# Read the arguments passed
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

#### PNG names ----------------------------------------------------------------

# Here we'll set up the file names for the broad histology specific oncoprints
# stored in the `oncoprint-landscape` module
lgat_png <- file.path(onco_dir, paste0(opt$lead_filename, "_lgat_goi_oncoprint.png"))
embryonal_png <- file.path(onco_dir, paste0(opt$lead_filename, "_embryonal_goi_oncoprint.png"))
hgat_png <- file.path(onco_dir, paste0(opt$lead_filename, "_hgat_goi_oncoprint.png"))
ependymal_png <- file.path(onco_dir, paste0(opt$lead_filename, "_ependymal_goi_oncoprint.png"))
other_cns_png <- file.path(onco_dir, paste0(opt$lead_filename, "_other_goi_oncoprint.png"))

#### Assemble multipanel figure ------------------------------------------------

oncoprint_figure <- multi_panel_figure(columns = 8,
                                            rows = 2,
                                            width = 1200,
                                            height = 300)

oncoprint_figure <- fill_panel(oncoprint_figure,
                                    lgat_png,
                                    col = 2:3,
                                    row = 1,
                                    scaling = "stretch",
                                    label = "Low-grade astrocytic tumor")

oncoprint_figure <- fill_panel(oncoprint_figure,
                                    embryonal_png,
                                    col = 4:5,
                                    row = 1,
                                    scaling = "stretch",
                                    label = "Embryonal tumor")

oncoprint_figure <- fill_panel(oncoprint_figure,
                                    hgat_png,
                                    col = 6:7,
                                    row = 1,
                                    scaling = "stretch",
                                    label = "Diffuse astrocytic and oligodendroglial tumor")

oncoprint_figure <- fill_panel(oncoprint_figure,
                                    ependymal_png,
                                    col = 2:3,
                                    row = 2,
                                    scaling = "stretch",
                                    label = "Ependymal tumor")

oncoprint_figure <- fill_panel(oncoprint_figure,
                                    other_cns_png,
                                    col = 4:5,
                                    row = 2,
                                    scaling = "stretch",
                                    label = "Other CNS")

save_multi_panel_figure(oncoprint_figure, opt$png_name)
