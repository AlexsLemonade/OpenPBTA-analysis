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
other_cns_png <- file.path(onco_dir, paste0(opt$lead_filename, "_other_goi_oncoprint.png"))

#### Assemble multipanel figure ------------------------------------------------

oncoprint_figure <- multi_panel_figure(columns = 4,
                                       rows = 2,
                                       width = 450,
                                       height = 265)

oncoprint_figure <- fill_panel(oncoprint_figure,
                                    lgat_png,
                                    col = 1:2,
                                    row = 1,
                                    scaling = "stretch",
                                    label = "A")

oncoprint_figure <- fill_panel(oncoprint_figure,
                                    embryonal_png,
                                    col = 3:4,
                                    row = 1,
                                    scaling = "stretch",
                                    label = "B")

oncoprint_figure <- fill_panel(oncoprint_figure,
                                    hgat_png,
                                    col = 1:2,
                                    row = 2,
                                    scaling = "stretch",
                                    label = "C")

oncoprint_figure <- fill_panel(oncoprint_figure,
                                    other_cns_png,
                                    col = 3:4,
                                    row = 2,
                                    scaling = "stretch",
                                    label = "D")

save_multi_panel_figure(oncoprint_figure, opt$png_name)
