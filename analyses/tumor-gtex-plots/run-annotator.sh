# run annotator on tsv files, convert to jsonl and gzip

# pan-cancer cancer-group level
Rscript --vanilla 02-chunkwise-annotate-and-zip.R \
--input_file results/pan_cancer_plots_cancer_group_level.tsv

# pan-cancer cohort-cancer-group level
Rscript --vanilla 02-chunkwise-annotate-and-zip.R \
--input_file results/pan_cancer_plots_cohort_cancer_group_level.tsv

# tumor-normal cancer-group level
Rscript --vanilla 02-chunkwise-annotate-and-zip.R \
--input_file results/tumor_normal_gtex_plots_cancer_group_level.tsv

# tumor-normal cohort-cancer-group level
Rscript --vanilla 02-chunkwise-annotate-and-zip.R \
--input_file results/tumor_normal_gtex_plots_cohort_cancer_group_level.tsv
