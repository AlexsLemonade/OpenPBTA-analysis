### calculate fishers exact test

row_fisher <- function(w, x, y, z){
  # function to calculate fisher test from row elements
  mat <- matrix(c(w, x, y, z), ncol = 2)
  fisher <- fisher.test(mat)
  return(fisher$p.value)
}


#### Calculate co-occurence for top genes

coocurrence <- function(gene_sample_df, 
                        genes= sample(unique(gene_sample_df$gene), 25), 
                        samples = unique(gene_sample_df$sample)){
  # gene_list in order of interest
  #
  # get all pairs of genes
  gene_pairs <- t(combn(genes, m = 2))
  colnames(gene_pairs) <- c("gene1", "gene2")
  
  # get mutation counts for all genes/sample pairs in gene list
  # fills in any missing values with 0
  all_sample_counts <- expand.grid(gene = genes,
                                   sample = samples, 
                                   stringsAsFactors = FALSE) %>%
    dplyr::left_join(gene_sample_df, by = c("gene", "sample")) %>%
    tidyr::replace_na(list(mutations = 0))
  
  gene_pair_counts <- gene_pairs %>% 
    tibble::as_tibble() %>%
    dplyr::left_join(all_sample_counts,
                     by = c("gene1" = "gene")) %>%
    dplyr::rename(muts1 = mutations) %>%
    dplyr::left_join(all_sample_counts,
                     by = c("gene2" = "gene", 
                            "sample" = "sample")) %>%
    dplyr::rename(muts2 = mutations) %>%
    dplyr::mutate(gene1 = factor(gene1, levels = genes),
                  gene2 = factor(gene2, levels = genes))
  
  gene_pair_summary <- gene_pair_counts %>%
    dplyr::group_by(gene1, gene2) %>%
    dplyr::summarize(mut11 = sum(muts1 > 0  & muts2 > 0), 
                     mut10 = sum(muts1 > 0  & muts2 == 0), 
                     mut01 = sum(muts1 == 0 & muts2 > 0),
                     mut00 = sum(muts1 == 0 & muts2 ==0),
                     odds_ratio = (mut11 * mut00) / (mut10 * mut01),
                     cooccur_sign = ifelse(odds_ratio > 1, 1, -1)) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(p = row_fisher(mut11, mut10, mut01, mut00), 
                  cooccur_score = cooccur_sign * -log10(p))
  
  return(gene_pair_summary)
}