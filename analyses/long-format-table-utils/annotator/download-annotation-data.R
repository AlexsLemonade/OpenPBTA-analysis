# Download gene Ensembl ENSG ID to gene full name and protein RefSeq IDs mapping
# file `annotation-data/ensg-gene-full-name-refseq-protein.tsv` from
# https://mygene.info/

# Import functions -------------------------------------------------------------
# Get %>% without loading the whole library
`%>%` <- dplyr::`%>%`



# Define functions -------------------------------------------------------------
# Collapse a list of refseq.protein character vectors from mygene.info query
# results
#
# Args:
# - xl: list of refseq.protein character vectors from mygene.info query results
collapse_rp_lists <- function(xl) {
  fxl <- purrr::discard(xl, is.null)
  # assert all characters
  sapply(fxl, function(x) {
    stopifnot(is.character(x))
  })
  fxv <- unique(purrr::reduce(fxl, c, .init = character(0)))
  np_fxv <- fxv[stringr::str_detect(fxv, "^NP_")]
  c_np_fxv <- paste(np_fxv, collapse = ",")
  return(c_np_fxv)
}
# # test cases
# collapse_rp_lists(list(NULL, NULL, NULL))
# collapse_rp_lists(list())
# collapse_rp_lists(list(c("NP_1", "NP_2"), c("NP_1", "NP_3")))
# collapse_rp_lists(list(c("NP_1", "NP_2"), NULL, c("NP_1", "NP_3")))
# collapse_rp_lists(list(c("NP_1", "NP_2")))



# Set up directory paths -------------------------------------------------------
# Adapted from the oncoprint-landscape module.
#
# Detect the ".git" folder -- this will in the project root directory. Use this
# as the root directory to ensure proper execution, no matter where it is called
# from.
#
# This only works if the working directory is OpenPedCan-analysis or a
# subdirectory of OpenPedCan-analysis
#
# root_dir is the absolute path of OpenPedCan-analysis
tryCatch(
  {
    root_dir <- rprojroot::find_root(rprojroot::has_file(".git/index"))
  },
  error = function(err_cond) {
    # adapted from http://adv-r.had.co.nz/Exceptions-Debugging.html
    err_cond$message <- paste0(
      err_cond$message,
      "\nTry re-running this script with working directory as ",
      "OpenPedCan-analysis or a subdirectory of OpenPedCan-analysis.\n")
    stop(err_cond)
  }
)

input_data_dir <- file.path(root_dir, "data")

output_data_dir <- file.path(
  root_dir, "analyses", "long-format-table-utils", "annotator",
  "annotation-data")

if (!dir.exists(output_data_dir)) {
  dir.create(output_data_dir)
}



# Read input data --------------------------------------------------------------
# ensg hugo rmtl mappings
ensg_hugo_rmtl_df <- dplyr::distinct(
  readr::read_tsv(file.path(input_data_dir, "ensg-hugo-rmtl-v1-mapping.tsv"),
                  col_types = readr::cols(.default = readr::col_guess())))

# assert all ensg_ids and gene_symbols are not NA
if (!identical(sum(is.na(ensg_hugo_rmtl_df$ensg_id)), as.integer(0))) {
  stop(paste0("Found NA in ensg-hugo-rmtl-v1-mapping.tsv ensg_id.\n",
              "Check if PedOT release data are downloaded properly.\n",
              "If data is downloaded properly, submit a GitHub data question."))
}

if (!identical(sum(is.na(ensg_hugo_rmtl_df$gene_symbol)), as.integer(0))) {
  stop(paste0("Found NA in ensg-hugo-rmtl-v1-mapping.tsv gene_symbol.\n",
              "Check if PedOT release data are downloaded properly.\n",
              "If data is downloaded properly, submit a GitHub data question."))
}

# assert all ensg_id are unique
if (!identical(length(unique(ensg_hugo_rmtl_df$ensg_id)),
               nrow(ensg_hugo_rmtl_df))) {
  stop(paste0("Found duplicated ensg_id in ensg-hugo-rmtl-v1-mapping.tsv.\n",
              "Check if PedOT release data are downloaded properly.\n",
              "If data is downloaded properly, submit a GitHub data question."))
}

# Download data from https://mygene.info/ --------------------------------------
message("Retrieve Gene_full_name and Protein_RefSeq_ID from mygene.info...")
ens_gids <- ensg_hugo_rmtl_df$ensg_id

mg_qres_list <- mygene::queryMany(
  ens_gids, scopes = "ensembl.gene", fields = c("refseq", "name"),
  species = "human", returnall = TRUE, return.as = "DataFrame")

mg_qres_df <- tibble::as_tibble(
  mg_qres_list$response[, c("query", "notfound", "name", "refseq.protein")]) %>%
  tidyr::replace_na(list(notfound = FALSE)) %>%
  dplyr::filter(!notfound) %>%
  dplyr::group_by(query) %>%
  dplyr::summarise(name = paste(unique(name), collapse = ", "),
                   refseq_protein = collapse_rp_lists(refseq.protein)) %>%
  dplyr::rename(Gene_Ensembl_ID = query, Gene_full_name = name,
                Protein_RefSeq_ID = refseq_protein)



# Output TSV file --------------------------------------------------------------
readr::write_tsv(
  mg_qres_df,
  file.path(output_data_dir, "ensg-gene-full-name-refseq-protein.tsv"))
