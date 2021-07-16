# Add gene and cancer_group annotations to long-format tables
#
# Args:
# - long_format_table: A data.frame or tibble that contains the following
#   columns as character:
#   - Gene_symbol: HUGO symbols, e.g. PHLPP1, TM6SF1, and DNAH5
#   - Gene_Ensembl_ID: Ensembl ENSG IDs without `.#` versions, e.g.
#     ENSG00000039139, ENSG00000111261, and ENSG00000169710
#   - Disease: The `cancer_group` in the `histologies.tsv`, e.g.
#     Adamantinomatous Craniopharyngioma, Atypical Teratoid Rhabdoid Tumor, and
#     Low-grade glioma/astrocytoma
# - is_gene_level_table: TRUE or FALSE on whether the table is gene-level. If is
#   gene-level, Gene_type column will be added, otherwise will not. Default
#   value is FALSE.
# - add_Protein_RefSeq_ID: TRUE or FALSE on whether to add Protein_RefSeq_ID
#   column. Default value is FALSE.
#
# Returns a data.frame or tibble, based on input table type, with additonal
# annotation columns
annotate_long_format_table <- function(long_format_table,
                                       is_gene_level_table = FALSE,
                                       add_Protein_RefSeq_ID = FALSE) {

}
