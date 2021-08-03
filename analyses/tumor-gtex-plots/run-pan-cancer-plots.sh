# 1. cohort + cancer_group level plots
# only tumors: GMKF + PBTA + TARGET
Rscript 01-tumor-gtex-plots.R \
--expr_mat '../../data/gene-expression-rsem-tpm-collapsed.rds' \
--hist_file '../../data/histologies.tsv' \
--map_file '../../data/ensg-hugo-rmtl-mapping.tsv' \
--cohort_list 'GMKF, PBTA, TARGET' \
--tumor_vs_normal FALSE \
--analysis_type 'cohort_cancer_group_level' \
--plot_width 10 \
--plot_height 9 \
--meta_file 'metadata.tsv'

# 2. cancer_group level plots
# only tumors: GMKF + PBTA + TARGET  
Rscript 01-tumor-gtex-plots.R \
--expr_mat '../../data/gene-expression-rsem-tpm-collapsed.rds' \
--hist_file '../../data/histologies.tsv' \
--map_file '../../data/ensg-hugo-rmtl-mapping.tsv' \
--cohort_list 'GMKF, PBTA, TARGET' \
--tumor_vs_normal FALSE \
--analysis_type 'cancer_group_level' \
--plot_width 10 \
--plot_height 9 \
--meta_file 'metadata.tsv'

