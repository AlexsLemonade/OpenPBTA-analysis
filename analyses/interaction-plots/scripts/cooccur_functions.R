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

#' Calculate confidence interval for odds ratio
#' 
#' Uses Haldane correction if any counts are zero
#' @return CI of odds ratio 
calc_ci <- function(mut11, mut10, mut01, mut00, odds_ratio){
  if(any(mut11==0, mut10==0, mut01==0, mut00==0)){
    adjusted_or <- (mut11 + 0.5)*(mut00 + 0.5)/((mut01 + 0.5)*(mut10 + 0.5))
    standard_error_or <- sqrt(1/(mut11+0.5) + 1/(mut10+0.5) + 1/(mut01+0.5) + 1/(mut00+0.5))
    or_ci_lower_bound <- min(exp(log(adjusted_or) - 1.96 * standard_error_or), odds_ratio)
    or_ci_upper_bound <- max(exp(log(adjusted_or) + 1.96 * standard_error_or), odds_ratio)
  } else {
    standard_error_or <- sqrt(1/mut11 + 1/mut10 + 1/mut01 + 1/mut00)
    or_ci_lower_bound <- exp(log(odds_ratio) - 1.96 * standard_error_or)
    or_ci_upper_bound <- exp(log(odds_ratio) + 1.96 * standard_error_or)
  }
  ci <- paste0(or_ci_lower_bound, ",", or_ci_upper_bound)
  return(ci)
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
#'   gene list with columns:
#'  `gene1`; `gene2` (gene names);
#'  `mut11`; `mut10`; `mut01` (counts of mutations in each category of sharing:
#'     `mut11`: both mutated; `mut10`: mutated in the first but not second gene, etc.);
#'  `odds_ratio` (odds ratio for for co-occurence);
#'  `cooccur_sign` (1 if co-occurence greater than by chance, -1 if less frequent than expected)
#'  `p` (the fisher's exact test p value);
#'  `q` (Benjamini-Hochberg adjusted p value);
#'  `cooccur_score` (calculated as `cooccur_sign * -log10(p)`);
#'  `n_mutated_gene1`; `n_mutated_gene2` (count of samples with gene1 or gene2 mutations); 
#'  `perc_mutated_gene1`; `perc_mutated_gene2` (percent of samples with gene 1 or gene 2 mutations);
#'  `perc_cooccur`;
#'  `perc_mutexl`;#' @return A data frame summarizing the co-occurence of pairs of genes in the
#'   gene list with columns:
#'  `gene1`; `gene2` (gene names);
#'  `mut11`; `mut10`; `mut01` (counts of mutations in each category of sharing:
#'     `mut11`: both mutated; `mut10`: mutated in the first but not second gene, etc.);
#'  `odds_ratio` (odds ratio for co-occurence);
#'  `cooccur_sign` (1 if co-occurence greater than by chance, -1 if less frequent than expected)
#'  `odds_ratio_ci` (confidence interval for odds_ratio);
#'  `p` (the fisher's exact test p value);
#'  `q` (Benjamini-Hochberg adjusted p value);
#'  `cooccur_score` (calculated as `cooccur_sign * -log10(p)`);
#'  `n_mutated_gene1`; `n_mutated_gene2` (count of samples with gene1 or gene2 mutations); 
#'  `perc_mutated_gene1`; `perc_mutated_gene2` (percent of samples with gene 1 or gene 2 mutations);
#'  `perc_cooccur`; (Percent of samples with mutations in gene 1 and/or gene 2 with co-occurring gene1/2 mutations);
#'  `perc_mutexcl`; (Percent of samples with mutations in gene 1 and/or gene 2 where gene1/2 mutations are mutually exclusive);
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
    dplyr::rename(muts2 = mutations) %>%
    dplyr::mutate(
      gene1 = factor(gene1, levels = genes),
      gene2 = factor(gene2, levels = genes)
    )

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
    dplyr::mutate(odds_ratio_ci = calc_ci(mut11, mut10, mut01, mut00, odds_ratio),
                  p = row_fisher(mut11, mut10, mut01, mut00)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      q = p.adjust(p, method = "BH"),
      cooccur_score = cooccur_sign * -log10(p)) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(    
      n_mutated_gene1 = sum(mut11, mut10),
      n_mutated_gene2 = sum(mut11, mut01),
      perc_mutated_gene1 = sum(mut11, mut10)*100/sum(mut11, mut10, mut01, mut00),
      perc_mutated_gene2 = sum(mut11, mut01)*100/sum(mut11, mut10, mut01, mut00),
      perc_cooccur = mut11*100/sum(mut11, mut10, mut01),
      perc_mutexcl = sum(mut10,mut01)*100/sum(mut11, mut10, mut01)
)
    
  return(gene_pair_summary)
}
