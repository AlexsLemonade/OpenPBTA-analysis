suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(rtracklayer))

option_list <- list(
  make_option(c("-i", "--freq_file"), type = "character",
              help = "Input filtered fusion from `fusion_filtering` module "),
  make_option(c("-a", "--efo_mondo_file"), type = "character",
              help = "efo-mondo-map.tsv from data release"),
  make_option(c("-c", "--ensg_hugo_rmtl_file"), type  = "character",
              help = "ensg-hugo-rmtl-v1-mapping.tsv from data release "),
  make_option(c("-o", "--output_filename"), type = "character", 
              default = "putative-oncogene-fusion-freq",
              help = "Output filename suffix ")
)

# parse the parameters
opt <- parse_args(OptionParser(option_list = option_list))
freq_file <- opt$freq_file
efo_mondo_file <- opt$efo_mondo
ensg_hugo_rmtl_file <- opt$ensg_hugo_rmtl
output_filename <- opt$output_filename

# Detect the ".git" folder -- this will in the project root directory.
# Use this as the root directory to ensure proper sourcing of functions no
# matter where this is called from
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Set path to results, plots, and subset files directories
data_dir <- file.path(root_dir, "data")
module_dir <- file.path(root_dir, "analyses", "fusion-frequencies")
results_dir <- file.path(module_dir, "results")

# Source function to gather counts and frequency calculation per cohort and cancer_group
source(file.path(module_dir,"utils","freq_counts.R"))

# Read in output file from 01-fusion-frequencies.R
m_fus_freq_tbl <- read.delim(freq_file,stringsAsFactors = FALSE) %>%
  # replace NA to "" in columns that have NA
replace_na(list(
                "Gene1A_anno"="",
                "Gene1B_anno"="",
                "Gene2A_anno"="",
                "Gene2B_anno"="",
                "Fusion_anno"="")
)
#efo cancer_group mappings
efo_mondo_cg_df <- read_tsv(efo_mondo_file,
                            col_types = cols(.default = col_guess())) %>%
  distinct()
# efo_mondo_cg_df$cancer_group[
#   !efo_mondo_cg_df$cancer_group %in% htl_df$cancer_group]

# assert all cancer_groups are not NA
stopifnot(identical(sum(is.na(efo_mondo_cg_df$cancer_group)), as.integer(0)))
# assert all cancer_groups are unique.
# result SNV table is left joined by cancer_groups
stopifnot(identical(length(unique(efo_mondo_cg_df$cancer_group)),
                    nrow(efo_mondo_cg_df)))

# ensg hugo rmtl mappings
ensg_hugo_rmtl_df <- read_tsv(ensg_hugo_rmtl_file,
                              col_types = cols(.default = col_guess())) %>%
  distinct()
# assert all ensg_ids and gene_symbols are not NA
stopifnot(identical(sum(is.na(ensg_hugo_rmtl_df$ensg_id)), as.integer(0)))
stopifnot(identical(sum(is.na(ensg_hugo_rmtl_df$gene_symbol)), as.integer(0)))
# assert all ensg_id are unique
# result SNV table is left joined by ensg_id
stopifnot(identical(length(unique(ensg_hugo_rmtl_df$ensg_id)),
                    nrow(ensg_hugo_rmtl_df)))


# Add EFO, MONDO, and RMTL to the output table ---------------------------------
ann_efo_mondo_cg_df <- efo_mondo_cg_df %>%
  dplyr::rename(Disease = cancer_group, EFO = efo_code, MONDO = mondo_code)

m_fus_freq_tbl <- m_fus_freq_tbl %>%
  left_join(ann_efo_mondo_cg_df, by = 'Disease') %>%
  replace_na(list(EFO = '', MONDO = ''))

# asert all rmtl NAs have version NAs, vice versa
stopifnot(identical(is.na(ensg_hugo_rmtl_df$rmtl),
                    is.na(ensg_hugo_rmtl_df$version)))
ann_ensg_hugo_rmtl_df <- ensg_hugo_rmtl_df %>%
  mutate(RMTL = if_else((!is.na(rmtl) & !is.na(version)),
                        paste0(rmtl, ' (', version, ')'),
                        '')) %>%
  dplyr::select(ensg_id, RMTL,gene_symbol) %>%
  dplyr::rename(Gene_Ensembl_ID = ensg_id,
                Gene_Symbol = gene_symbol) 

m_fus_freq_tbl <- m_fus_freq_tbl %>%
  left_join(ann_ensg_hugo_rmtl_df, by = "Gene_Symbol") %>%
  replace_na(list(RMTL = '',
                  Gene_Ensembl_ID = '')) 

stopifnot(identical(sum(is.na(m_fus_freq_tbl)), as.integer(0)))


# Add additional gene annotations ---------------------------------------------------
message('Retrieve Gene_full_name and Protein_RefSeq_ID from mygene.info...')
ens_gids <- ensg_hugo_rmtl_df$ensg_id[which(ensg_hugo_rmtl_df$gene_symbol %in%
                                              m_fus_freq_tbl$Gene_Symbol)]
ens_gids <- ens_gids[!is.na(ens_gids)]

mg_qres_list <- mygene::queryMany(
  ens_gids, scopes = 'ensembl.gene', fields = c('refseq', 'name'),
  species = 'human', returnall = TRUE, return.as = 'DataFrame')

mg_qres_df <- as_tibble(
  mg_qres_list$response[, c('query', 'notfound', 'name', 'refseq.protein')]) %>%
  replace_na(list(notfound = FALSE)) %>%
  dplyr::filter(!notfound) %>%
  group_by(query) %>%
  summarise(name = paste(unique(name), collapse = ', '),
            refseq_protein = collapse_rp_lists(refseq.protein)) %>%
  dplyr::rename(Gene = query, Gene_full_name = name,
         Protein_RefSeq_ID = refseq_protein)

# add additional fields to alteration df
m_fus_freq_tbl <- m_fus_freq_tbl %>%
  left_join(mg_qres_df, by = c('Gene_Ensembl_ID'='Gene')) %>%
  replace_na(list(Gene_Ensembl_ID='',
                  Gene_full_name=''))

suppressMessages(gc(reset = TRUE))


jsonlite::write_json(
  m_fus_freq_tbl,
  file.path(results_dir, paste0(output_filename,'.json'))
  )
