# Consensus mutation file creation
# 2019
# C. Savonen for ALSF - CCDL
#
# Purpose: Save consensus mutation calls
# Because of VarDict's oversensitivity, we will ignore VarDict's calls and
# keep mutations that are agreed upon by Lancet, Mutect2, and Strelka2.

#################################### Set Up ####################################
# Establish base dir
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Import this for the recalculation of TMB for the consensus mutation calls
source(file.path(root_dir, "analyses", "snv-callers", "util", "wrangle_functions.R"))

# Magrittr pipe
`%>%` <- dplyr::`%>%`

# Let's get the mutation ids for these.
consen_mutations <- callers_per_mutation %>%
  dplyr::filter(caller_combo == "lancet-mutect2-strelka2") %>%
  dplyr::select(-caller_combo) %>%
  dplyr::pull(mutation_id)

# Isolate the mutations to only these mutations, use the Strelka2 stats.
consen_mutation <- vaf_df %>%
  dplyr::filter(
    caller == "strelka2",
    mutation_id %in% consen_mutations
  ) %>%
  readr::write_tsv(file.path(consensus_results_dir, "consensus_mutation.maf.tsv"))

# Zip this file up.
zip(
  file.path(consensus_results_dir, "consensus_mutation.maf.tsv.zip"),
  file.path(consensus_results_dir, "consensus_mutation.maf.tsv")
)

# Re-calculate TMB 
# Set up BED region files for TMB calculations
wgs_bed <- readr::read_tsv(file.path("..", "..", "data", "WGS.hg38.strelka2.unpadded.bed"),
  col_names = FALSE
)
wxs_bed <- readr::read_tsv(file.path("..", "..", "data", "WXS.hg38.100bp_padded.bed"),
  col_names = FALSE
)

# Calculate size of genome surveyed
wgs_genome_size <- sum(wgs_bed[, 3] - wgs_bed[, 2])
wxs_exome_size <- sum(wxs_bed[, 3] - wxs_bed[, 2])

# Calculate TMBs and write to TMB file
tmb_df <- calculate_tmb(vaf_df,
  wgs_size = wgs_genome_size,
  wxs_size = wxs_exome_size
) %>%
  readr::write_tsv(file.path(consensus_results_dir, "consensus_mutation_tmb.tsv"))

# Give message
message("Consensus mutations and TMB re-calculations saved in: results/consensus")
