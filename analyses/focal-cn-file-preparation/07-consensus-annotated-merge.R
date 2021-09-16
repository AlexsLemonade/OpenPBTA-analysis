# This script converts merges consensus seg files with cnvkit annotated files 
# for both autosomes and x_and_y. The autosomes and x_and_y files are then merged
# to generate one single file

# #### Example Usage
#
# This script is intended to be run via the command line.
# This example assumes it is being run from the root of the repository.
#
# Rscript --vanilla 07-consensus-annotated-merge.R \
#   --cnvkit_auto results/cnvkit_annotated_cn_wxs_autosomes.tsv.gz \
#   --cnvkit_x_and_y results/cnvkit_annotated_cn_wxs_x_and_y.tsv.gz \
#   --consensus_auto results/consensus_seg_annotated_cn_autosomes.tsv.gz \
#   --consensus_x_and_y results/consensus_seg_annotated_cn_x_and_y.tsv.gz \
#   --outdir results

#### Set Up --------------------------------------------------------------------

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

#### Start --------------------------------------------------------------------

# Declare command line options
option_list <- list(
  optparse::make_option(
    c("--cnvkit_auto"),
    type = "character",
    default = "results/cnvkit_annotated_cn_wxs_autosomes.tsv.gz",
    help = "annoatated cnvkit cnv calls on autosomes for WXS samples"
  ),
  optparse::make_option(
    c("--cnvkit_x_and_y"),
    type = "character",
    default = "results/cnvkit_annotated_cn_wxs_x_and_y.tsv.gz",
    help = "annoatated cnvkit cnv calls on x and y for WXS samples"
  ),
  optparse::make_option(
    c("--consensus_auto"),
    type = "character",
    default = "results/consensus_seg_annotated_cn_autosomes.tsv.gz",
    help = "annoatated consensus cnv calls on autosomes for WGS samples"
  ),
  optparse::make_option(
    c("--consensus_x_and_y"),
    type = "character",
    default = "results/consensus_seg_annotated_cn_x_and_y.tsv.gz",
    help = "annoatated consensus cnv calls on x and y for WGS samples"
  ),
  optparse::make_option(
    c("--outdir"),
    type = "character",
    default = "results",
    help = "path to save the final merged files"
  )
)

# Read the arguments passed
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

# Read in the files 
cnvkit_auto <- data.table::fread(opt$cnvkit_auto)
cnvkit_x_and_y <- data.table::fread(opt$cnvkit_x_and_y)
consensus_auto <- data.table::fread(opt$consensus_auto)
consensus_x_and_y <- data.table::fread(opt$consensus_x_and_y)

# merge WGS (consensus) and WXS (cnvkit) autosomes and x_and_y respectively
merged_auto <- rbind(cnvkit_auto, consensus_auto)
merged_x_and_y <- rbind(cnvkit_x_and_y, consensus_x_and_y)

# write out the files
readr::write_tsv(merged_auto, 
          file.path(opt$outdir, "consensus_wgs_plus_cnvkit_wxs_autosomes.tsv.gz"))

readr::write_tsv(merged_x_and_y, 
          file.path(opt$outdir, "consensus_wgs_plus_cnvkit_wxs_x_and_y.tsv.gz"))

# Merge autosomes with X and Y and output a combined file 
# For the combined file, we do not need germline_sex_estimate from x_and_y
merged_x_and_y <- merged_x_and_y %>% dplyr::select(-germline_sex_estimate)
combined <- rbind(merged_auto, merged_x_and_y)

readr::write_tsv(combined, 
          file.path(opt$outdir, "consensus_wgs_plus_cnvkit_wxs.tsv.gz"))




