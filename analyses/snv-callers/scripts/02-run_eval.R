# Plot the variant caller data and print out a report html file.
#
# 2019
#
# C. Savonen for ALSF - CCDL
#
# Option descriptions
# --label : Label to be used for folder and all output. eg. 'strelka2'. Optional.
#      Default is 'maf'
# --plot_type : Specify what kind of plots you want printed out. Must be
#               compatible with ggsave. eg pdf. Default is png
# --output : Where you would like the output from this script to be stored.
# --sql_file : File path to where the SQL file was saved in 00-set_up.R.
# --strategy : Specify whether you would like WXS and WGS separated for the plots.
#              Analysis is still done on all data in the MAF file regardless.
#              Acceptable options are 'wgs', 'wxs' or 'both', both for if you
#              don't want to separate them. Default is both.
# --cosmic : Relative file path to COSMIC file to be analyzed. Assumes file path
#            is given from top directory of 'OpenPBTA-analysis'.
# --overwrite : If TRUE, will overwrite any reports of the same name. Default is
#              FALSE
# --no_region : If used, regional analysis will not be done.

#
# Command line example:
#
# Rscript 02-run_eval.R \
# -l strelka2 \
# -p png \
# -o strelka2 \
# -s wxs 
#
# Establish base dir
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Magrittr pipe
`%>%` <- dplyr::`%>%`

# Import special functions
source(file.path(root_dir, "analyses", "snv-callers", "util", "wrangle_functions.R"))
source(file.path(root_dir, "analyses", "snv-callers", "util", "plot_functions.R"))

# Load library:
library(optparse)

#--------------------------------Set up options--------------------------------#
# Set up optparse options
option_list <- list(
  make_option(
    opt_str = c("-l", "--label"), type = "character",
    default = "maf", help = "Label to be used for folder and all
                output. eg. 'strelka2'. Optional. Default is 'maf'",
    metavar = "character"
  ),
  make_option(
    opt_str = c("-p", "--plot_type"), type = "character",
    default = "png", help = "Specify what kind of plots you want
                printed out. Must be compatible with ggsave. eg pdf.
                Default is png.",
    metavar = "character"
  ),
  make_option(
    opt_str = "--sql_file", type = "character", default = "none",
    help = "File path to where the SQL file was saved in 00-set_up.R",
    metavar = "character"
  ),
  make_option(
    opt_str = c("-o", "--output"), type = "character",
    default = NULL, help = "Path to folder where you would like the
              output from this script to be stored.",
    metavar = "character"
  ),
  make_option(
    opt_str = c("-s", "--strategy"), type = "character",
    default = "both", help = "Specify whether you would like WXS and
                WGS separated for the plots. Can state all three with commas in
                between such as 'wgs,wxs,both'. Acceptable options are 'wgs',
                'wxs' or 'both', both for if you don't want to separate them.
                Default is both.",
    metavar = "character"
  ),
  make_option(
    opt_str = c("-c", "--cosmic"), type = "character", default = "none",
    help = "Relative file path (assuming from top directory of
              'OpenPBTA-analysis') to COSMIC file to be analyzed.",
    metavar = "character"
  ),
  make_option(
    opt_str = "--overwrite", action = "store_true",
    default = FALSE, help = "If TRUE, will overwrite any reports of
              the same name. Default is FALSE",
    metavar = "character"
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

########################### Check options specified ############################
# Normalize this file path
opt$sql_file <- file.path(root_dir, opt$sql_file)
opt$cosmic <- file.path(root_dir, opt$cosmic)

# Check the output directory exists
if (!dir.exists(opt$sql_file)) {
  stop(paste("Error:", opt$sql_file, "does not exist"))
}

# Check for COSMIC file
if (!file.exists(opt$cosmic)) {
  stop(paste("Error:", opt$cosmic, "does not exist"))
}

# The list of needed file suffixes
needed_files <- c(opt$sql_file, opt$cosmic)

# Get list of which files were found
files_found <- file.exists(needed_files)

# Report error if any of them aren't found
if (!all(files_found)) {
  stop(paste("\n Could not find needed file(s):",
             needed_files[which(!files_found)],
             "Check your options and set up.",
             sep = "\n"
  ))
}

# List plot types we can take. This is base on ggsave's documentation
acceptable_plot_types <- c(
  "eps", "ps", "tex", "pdf", "jpeg", "tiff", "png",
  "bmp", "svg"
)

# Check the plot type option
if (!(opt$plot_type %in% acceptable_plot_types)) {
  stop("Error: unrecognized plot type specified. Only plot types accepted by
       ggplot2::ggsave may be used.")
}
# Add the period
opt$plot_type <- paste0(".", opt$plot_type)

################################### Set Up #####################################
# Set and make the plots directory
opt$output <- file.path(root_dir, opt$output)

# Make caller specific plots folder
if (!dir.exists(opt$output)) {
  dir.create(opt$output, recursive = TRUE)
}

# Make a list of the plot suffixes
plot_suffixes <- c("_base_change", "_depth_vs_vaf", "_cosmic_plot", "_tmb_plot")

if (opt$no_region) {
  plot_suffixes <- c(plot_suffixes, "_snv_region")
}

# Make the plot names with specified prefix
plot_names <- paste0(plot_suffixes, opt$plot_type)

############################ Check strategy option #############################
# Reformat the strategy option into lower case and vector
opt$strategy <- tolower(unlist(strsplit(opt$strategy, ",")))

# Check strategy options
if (!all(opt$strategy %in% c("wgs", "wxs", "both"))) {
  stop("Error: unrecognized --strategy option. Acceptable options are 'wgs',
       'wxs' or 'both'. Multiple can be specified at once.")
}

#################### Run this for each experimental strategy ###################
for (strategy in opt$strategy) {
  # File paths plots we will create
  plot_paths <- file.path(
    opt$output,
    paste0(opt$label, "_", strategy, plot_names)
  )
  # Bring along the plot names
  names(plot_paths) <- plot_names

  ################## Plot the data using special functions #####################
  # Base call barplot
  base_change_plot(vaf_df, exp_strategy = strategy)
  ggplot2::ggsave(filename = plot_paths["_base_change.png"], plot = ggplot2::last_plot())

  # Read depth and VAF
  depth_vs_vaf_plot(vaf_df, exp_strategy = strategy)
  ggplot2::ggsave(filename = plot_paths["_depth_vs_vaf.png"], plot = ggplot2::last_plot())

  # Percent variants in COSMIC
  cosmic_plot(vaf_df, exp_strategy = strategy, opt$cosmic)
  ggplot2::ggsave(filename = plot_paths["_cosmic_plot.png"], plot = ggplot2::last_plot())

  # TMB by histology
  tmb_plot(tmb_df, x_axis = "short_histology", exp_strategy = strategy)
  ggplot2::ggsave(filename = plot_paths["_tmb_plot.png"], plot = ggplot2::last_plot())

  if (opt$no_region) {
    # Genomic region breakdown
    snv_region_plot(maf_annot, exp_strategy = strategy)
    ggplot2::ggsave(filename = plot_paths["_snv_region.png"], plot = ggplot2::last_plot())
  }

  ######################## Make plots into a report ############################
  # Make a summary report about the variant caller and strategy
  output_file <- file.path(
    opt$output,
    paste0(opt$label, "_", strategy, "_report.Rmd")
  )

  # Path to the template file
  template_folder <- file.path(
    root_dir, "analyses", "snv-callers", "template"
  )

  # Designate which template file name
  if (opt$no_region) {
    template_file <- file.path(template_folder, "variant_caller_report_template.Rmd")
  } else {
    template_file <- file.path(template_folder, "variant_caller_report_no_region_template.Rmd")
  }

  # Make copy of template
  if (file.exists(template_file)) {
    file.copy(from = template_file, to = output_file, overwrite = opt$overwrite)
  } else {
    stop(cat("The Rmd template file ", template_file, " does not exist."))
  }

  # Run this notebook
  rmarkdown::render(output_file, "html_document")
}
