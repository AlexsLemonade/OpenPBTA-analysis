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
# - rp_vec_list: list of refseq.protein character vectors from mygene.info query
#   results
#
# Returns a single character value of a comma-separated refseq.protein value or
# NA
collapse_rp_lists <- function(rp_vec_list) {
  # remove list elements that are NULL
  rm_null_rp_vec_list <- purrr::discard(rp_vec_list, is.null)
  # assert non-null list elements are all characters
  lapply(rm_null_rp_vec_list, function(x) {
    if (!is.character(x)) {
      stop(paste0("mygene query returns non-character refseq.protein.\n",
                  "Check query results. Revise download-annotation-data.R ",
                  "to handle non-character values."))
    }
  })

  # combined vector of unique refseq.protein values
  uniq_c_rm_null_rp_vec <- sort(unique(
    purrr::reduce(rm_null_rp_vec_list, c, .init = character(0))))
  # remove NA
  rm_na_uniq_c_rm_null_rp_vec <- purrr::discard(
    uniq_c_rm_null_rp_vec, is.na)
  # keep only NP_### refseq.protein values
  np_rp_vec <- purrr::keep(
    rm_na_uniq_c_rm_null_rp_vec, function(x) stringr::str_detect(x, "^NP_"))

  clp_np_rp_str <- paste(np_rp_vec, collapse = ",")

  if (identical(clp_np_rp_str, "")) {
    # no NP_### value for the gene query
    #
    # dplyr::summarise needs function return value to be character for a
    # character column. NA_character_ is still NA in character class rather than
    # "NA".
    return(NA_character_)
  } else {
    return(clp_np_rp_str)
  }
}



# Collapse a name character vector from mygene.info query results
#
# Args:
# - gn_vec: character vector of name values from mygene.info query results
#
# Returns a single character value of a comma-separated name value or NA
collapse_name_vec <- function(gn_vec) {
  # assert non-null list elements are all characters
  if (!is.character(gn_vec)) {
    stop(paste0("mygene query returns non-character refseq.protein.\n",
                "Check query results. Revise download-annotation-data.R ",
                "to handle non-character values."))
  }
  # remove NAs
  uniq_rm_na_gn_vec <- sort(unique(purrr::discard(gn_vec, is.na)))

  clp_gn_str <- paste(uniq_rm_na_gn_vec, collapse = ",")
  if (identical(clp_gn_str, "")) {
    # no non-NA name for the gene query
    #
    # dplyr::summarise needs function return value to be character for a
    # character column. NA_character_ is still NA in character class rather than
    # "NA".
    return(NA_character_)
  } else {
    return(clp_gn_str)
  }
}



# Set up directory paths -------------------------------------------------------
# Detect the ".git" folder -- this will be in the project root directory. Use
# this as the root directory to ensure proper execution, no matter where it is
# called from.
#
# This only works if the working directory is OpenPedCan-analysis or a
# subdirectory of OpenPedCan-analysis
#
# root_dir is the absolute path of OpenPedCan-analysis
#
# Adapted from the oncoprint-landscape module
#
# rprojroot::has_file(".git/index") returns a rprojroot::root_criterion, and
# main git working tree, created by git clone and git init, has the .git/index
# file
#
# rprojroot::has_file(".git") returns a rprojroot::root_criterion, and linked
# git working tree, created by git worktree add, has the .git file
#
# "Root criteria can be combined with the | operator. The result is a
# composite root criterion that requires either of the original criteria to
# match." -- help("root_criterion", "rprojroot") rprojroot_1.3-2
tryCatch(
  {
    root_dir <- rprojroot::find_root(
      rprojroot::has_file(".git/index") | rprojroot::has_file(".git"))
  },
  error = function(err_cond) {
    # adapted from http://adv-r.had.co.nz/Exceptions-Debugging.html
    err_cond$message <- paste0(
      err_cond$message,
      "\nTry re-running this function with working directory as ",
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
  readr::read_tsv(file.path(input_data_dir, "ensg-hugo-rmtl-mapping.tsv"),
                  col_types = readr::cols(.default = readr::col_guess())))

# assert all ensg_ids and gene_symbols are not NA
if (!identical(sum(is.na(ensg_hugo_rmtl_df$ensg_id)), as.integer(0))) {
  stop(paste0("Found NA in ensg-hugo-rmtl-mapping.tsv ensg_id.\n",
              "Check if PedOT release data are downloaded properly.\n",
              "If data is downloaded properly, submit a GitHub data issue."))
}

if (!identical(sum(is.na(ensg_hugo_rmtl_df$gene_symbol)), as.integer(0))) {
  stop(paste0("Found NA in ensg-hugo-rmtl-mapping.tsv gene_symbol.\n",
              "Check if PedOT release data are downloaded properly.\n",
              "If data is downloaded properly, submit a GitHub data issue."))
}

# assert all ensg_id are unique
if (!identical(length(unique(ensg_hugo_rmtl_df$ensg_id)),
               nrow(ensg_hugo_rmtl_df))) {
  stop(paste0("Found duplicated ensg_id in ensg-hugo-rmtl-mapping.tsv.\n",
              "Check if PedOT release data are downloaded properly.\n",
              "If data is downloaded properly, submit a GitHub data issue."))
}

# Download data from https://mygene.info/ --------------------------------------
message("Retrieve Gene_full_name and Protein_RefSeq_ID from mygene.info...")
ens_gids <- ensg_hugo_rmtl_df$ensg_id

mg_qres_list <- mygene::queryMany(
  ens_gids, scopes = "ensembl.gene", fields = c("refseq", "name"),
  species = "human", returnall = TRUE, return.as = "DataFrame")

found_mg_qres_df <- tibble::as_tibble(
  mg_qres_list$response[, c("query", "notfound", "name", "refseq.protein")]) %>%
  tidyr::replace_na(list(notfound = FALSE)) %>%
  dplyr::filter(!notfound)

# remove rows that have both NA name and NULL refseq.protein
rm_bnn_found_mg_qres_df <- found_mg_qres_df %>%
  dplyr::filter(!(is.na(name) & purrr::map_lgl(refseq.protein, is.null)))

# collapses name and refseq.protein for output
out_rm_bnn_found_mg_qres_df <- rm_bnn_found_mg_qres_df %>%
  dplyr::group_by(query) %>%
  dplyr::summarise(name = collapse_name_vec(name),
                   refseq_protein = collapse_rp_lists(refseq.protein)) %>%
  dplyr::rename(Gene_Ensembl_ID = query, Gene_full_name = name,
                Protein_RefSeq_ID = refseq_protein)



# Output TSV file --------------------------------------------------------------
readr::write_tsv(
  out_rm_bnn_found_mg_qres_df,
  file.path(output_data_dir, "ensg-gene-full-name-refseq-protein.tsv"))

message("Done running download-annotation-data.R")
