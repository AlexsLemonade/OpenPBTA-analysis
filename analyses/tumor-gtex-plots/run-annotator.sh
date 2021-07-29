# run annotator on output files

# pan-cancer cancer-group level
Rscript --vanilla ../../analyses/long-format-table-utils/annotator/annotator-cli.R \
-r -v -c MONDO,RMTL,EFO \
-i results/pan_cancer_plots_cancer_group_level.tsv \
-o results/pan_cancer_plots_cancer_group_level.tsv

# pan-cancer cohort-cancer-group level
Rscript --vanilla ../../analyses/long-format-table-utils/annotator/annotator-cli.R \
-r -v -c MONDO,RMTL,EFO \
-i results/pan_cancer_plots_cohort_cancer_group_level.tsv \
-o results/pan_cancer_plots_cohort_cancer_group_level.tsv

# tumor-normal cancer-group level
Rscript --vanilla ../../analyses/long-format-table-utils/annotator/annotator-cli.R \
-r -v -c MONDO,RMTL,EFO \
-i results/tumor_normal_gtex_plots_cancer_group_level.tsv \
-o results/tumor_normal_gtex_plots_cancer_group_level.tsv

# tumor-normal cohort-cancer-group level
Rscript --vanilla ../../analyses/long-format-table-utils/annotator/annotator-cli.R \
-r -v -c MONDO,RMTL,EFO \
-i results/tumor_normal_gtex_plots_cohort_cancer_group_level.tsv \
-o results/tumor_normal_gtex_plots_cohort_cancer_group_level.tsv