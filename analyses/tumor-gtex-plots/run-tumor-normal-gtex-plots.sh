# 1. cohort + cancer_group level plots
# tumor vs normal: GMKF + PBTA + TARGET + GTEx
Rscript 01-tumor-gtex-plots.R \
--expr_mat '../../data/gene-expression-rsem-tpm-collapsed.rds' \
--hist_file '../../data/histologies.tsv' \
--map_file '../../data/ensg-hugo-rmtl-mapping.tsv' \
--cohort_list 'GMKF, PBTA, TARGET, GTEx' \
--tumor_vs_normal TRUE \
--analysis_type 'cohort_cancer_group_level' \
--plot_width 13 \
--plot_height 9 \
--meta_file 'metadata.tsv'

# 2. cancer_group level plots
# tumor vs normal: GMKF + PBTA + TARGET + GTEx 
Rscript 01-tumor-gtex-plots.R \
--expr_mat '../../data/gene-expression-rsem-tpm-collapsed.rds' \
--hist_file '../../data/histologies.tsv' \
--map_file '../../data/ensg-hugo-rmtl-mapping.tsv' \
--cohort_list 'GMKF, PBTA, TARGET, GTEx' \
--tumor_vs_normal TRUE \
--analysis_type 'cancer_group_level' \
--plot_width 13 \
--plot_height 9 \
--meta_file 'metadata.tsv'