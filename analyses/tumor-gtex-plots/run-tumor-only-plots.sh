# For this module, we will be combining CBTN and PNOC to a single cohort PBTA

# 1. cohort + cancer_group level plots
# only tumors: GMKF + PBTA (CBTN, PNOC)
Rscript 01-tumor-gtex-plots.R \
--expr_mat '../../data/v5/gene-expression-rsem-tpm-collapsed.rds' \
--hist_file '../../data/v5/histologies.tsv' \
--cohort_list 'GMKF, CBTN, PNOC' \
--tumor_vs_normal FALSE \
--analysis_type 'cohort_cancer_group_level' \
--plot_width 10 \
--plot_height 9 \
--mapping_file 'metadata.tsv'

# 2. cancer_group level plots
# only tumors: GMKF + PBTA (CBTN, PNOC) 
Rscript 01-tumor-gtex-plots.R \
--expr_mat '../../data/v5/gene-expression-rsem-tpm-collapsed.rds' \
--hist_file '../../data/v5/histologies.tsv' \
--cohort_list 'GMKF, CBTN, PNOC' \
--tumor_vs_normal FALSE \
--analysis_type 'cancer_group_level' \
--plot_width 10 \
--plot_height 9 \
--mapping_file 'metadata.tsv'

