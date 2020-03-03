# cooccur_function.R
# Functions for calculating co-occurence between mutations

#' Calculate Fisher's exact test for row-wise data
#'
#' Data order follows a two by two matrix, filled by rows or columns
#' (as these are equivalent). Not vectorized!
#'
#' @param w row 1 column 1 value
#' @param x row 2 column 1 value
#' @param y row 1 column 2 value
#' @param z row 2 column 2 value
#'
#' @return p-value from a Fisher's exact test
row_fisher <- function(w, x, y, z) {
  # function to calculate fisher test from row elements
  mat <- matrix(c(w, x, y, z), ncol = 2)
  fisher <- fisher.test(mat)
  return(fisher$p.value)
}

#' Filter mutations based on VAF and effect
#'
#' @param maf_df a data frame of a maf file or subset of one. Minimally includes
#'   `Consequence` column and either `vaf` or both `t_ref_count` and `t_alt_count`.
#' @param min_vaf The minimum VAF to include.
#' @param min_depth The minimum sequencing depth to include.
#' @param exclude_consequence A vector of consequences (as strings) that should
#'   be excluded in the filtered results. Note that if a mutation has multiple
#'   consequences, it will not be excluded unless ALL consequences are excluded.
#' @param include_consequence A vector of consequences (as strings) that should
#'   be included in the filtered results.
#'
#' @return A filtered data frame, with the same columns as the input `maf_df`,
#'   with VAF added if not already present.
filter_mutations <- function(maf_df,
                             min_vaf = 0,
                             min_depth = 0,
                             exclude_var_class = c(),
                             include_var_class = c()) {
  if (length(exclude_var_class) > 0 && length(include_var_class) > 0) {
    warning(
      "Both excluding and including consequences at the same time may have unexpected results."
    )
  }
  if (!("vaf" %in% colnames(maf_df))) {
    maf_df <- maf_df %>%
      dplyr::mutate(vaf = t_alt_count / (t_ref_count + t_alt_count))
  }

  # filter on vaf, depth, and variant class
  maf_df <- maf_df %>%
    dplyr::filter(
      vaf >= min_vaf,
      (t_ref_count + t_alt_count) > min_depth
    ) %>%
    dplyr::filter(!(Variant_Classification %in% exclude_var_class))

  if (length(include_var_class) > 0) {
    maf_df <- maf_df %>%
      dplyr::filter(Variant_Classification %in% include_var_class)
  }
}




#' Calculate co-occurence relationships for a list of genes
#'
#' @param gene_sample_df a data frame with columns `gene`, `sample`, and
#'   `mutations` where each row represents a gene mutated in the named
#'   sample, and `mutations` is the number of mutations in that gene-sample
#'   combination.
#' @param genes a vector of genes for which co-occurence data is to be calculated.
#'   Default is to pick a random subset, though this is almost never what
#'   is wanted! (Note that all-by-all comparisons is likely to be extremely
#'   slow).
#' @param samples Which samples should the comparison include. Defaults to all
#'   samples present in the data frame.
#'
#' @return A data frame summarizing the co-occurence of pairs of genes in the
#'   gene list with columns `gene1`; `gene2`; counts of each mutations in
#'   each category of sharing (`mut11`: both mutated; `mut10`: mutated in
#'   the first but not second gene, etc.); `odds_ratio` for co-occurence,
#'   `cooccur_sign`: 1 if co-occurence greater than by chance, -1 if less
#'   frequent than expected; `p` the fisher's exact test p value; and a
#'   `cooccur_score` calucated as `cooccur_sign * -log10(p)`.
coocurrence <- function(gene_sample_df,
                        genes = sample(unique(gene_sample_df$gene), 25),
                        samples = unique(gene_sample_df$sample)) {
  # gene_list in order of interest
  #
  # get all pairs of genes
  gene_pairs <- t(combn(genes, m = 2))
  colnames(gene_pairs) <- c("gene1", "gene2")

  # get mutation counts for all genes/sample pairs in gene list
  # fills in any missing values with 0
  all_sample_counts <- expand.grid(
    gene = genes,
    sample = samples,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::left_join(gene_sample_df, by = c("gene", "sample")) %>%
    tidyr::replace_na(list(mutations = 0))

  gene_pair_counts <- gene_pairs %>%
    tibble::as_tibble() %>%
    dplyr::left_join(all_sample_counts,
      by = c("gene1" = "gene")
    ) %>%
    dplyr::rename(muts1 = mutations) %>%
    dplyr::left_join(all_sample_counts,
      by = c(
        "gene2" = "gene",
        "sample" = "sample"
      )
    ) %>%
    dplyr::rename(muts2 = mutations)

  gene_pair_summary <- gene_pair_counts %>%
    dplyr::group_by(gene1, gene2) %>%
    dplyr::summarize(
      mut11 = sum(muts1 > 0 & muts2 > 0),
      mut10 = sum(muts1 > 0 & muts2 == 0),
      mut01 = sum(muts1 == 0 & muts2 > 0),
      mut00 = sum(muts1 == 0 & muts2 == 0),
      odds_ratio = (mut11 * mut00) / (mut10 * mut01),
      cooccur_sign = ifelse(odds_ratio > 1, 1, -1)
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(p = row_fisher(mut11, mut10, mut01, mut00)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      q = p.adjust(p, method = "BH"),
      cooccur_score = cooccur_sign * -log10(q)
    )

  return(gene_pair_summary)
}
