library(tidyverse)
library(BatchQC)
library(sva)

# Get correct file paths to data
library(rprojroot)
root_dir = find_root(has_file("OpenPBTA-analysis.Rproj"))
analysis_dir = file.path(root_dir, "analyses", "batch-effects")
data_dir = file.path(root_dir, "data")
functions = file.path(analysis_dir, "util", "functions.R")
source(functions)

# download all covariate data which will be used to identify batches
covariate_file = file.path(data_dir, "pbta-histologies.tsv")
covariate = read_tsv(covariate_file, col_types = cols(molecular_subtype = "c"))

# download gene expression data
dat_rsem_polya_file = file.path(data_dir, "pbta-gene-expression-rsem-tpm.polya.rds")
dat_rsem_polya = readRDS(dat_rsem_polya_file)

dat_rsem_stranded_file = file.path(data_dir, "pbta-gene-expression-rsem-tpm.stranded.rds")
dat_rsem_stranded = readRDS(dat_rsem_stranded_file)

dat_kallisto_stranded_file = file.path(data_dir, "pbta-gene-expression-kallisto.stranded.rds")
dat_kallisto_stranded <- readRDS(dat_kallisto_stranded_file)
dat_kallisto_stranded = dat_kallisto_stranded[,2:ncol(dat_kallisto_stranded)]

dat_kallisto_polya_file = file.path(data_dir, "pbta-gene-expression-kallisto.polya.rds")
dat_kallisto_polya <- readRDS(dat_kallisto_polya_file)

print("DOWNLOADED DATA")

### coment out section below when running the full analysis ###

# dat_rsem_polya = shorten(dat_rsem_polya)
# dat_rsem_stranded = shorten(dat_rsem_stranded)
# dat_kallisto_stranded = shorten(dat_kallisto_stranded)
# dat_kallisto_polya = shorten(dat_kallisto_polya)
# print("DATA HAS BEEN SHORTENED")

### coment out section above when running the full analysis ###

# Summarize transcripts at the gene level by averageing all transcipts values per gene
# This step only needs to be done with kallisto gene expression files
dat_kallisto_polya = grouper(dat_kallisto_polya)
dat_kallisto_stranded = grouper(dat_kallisto_stranded)

print("DATA IS GROUPED")

# See comments on make_id_batch_histology
id_batch_histology = make_id_batch_histology(covariate, "cohort")

# Run batchQC and combat. This will run batchQC before and after running combat. It will also perform join operations between id_batch_histology and our gene expression files 
make_histology_pca_plots(dat_rsem_polya, id_batch_histology, gene_id = dat_rsem_polya$gene_id, report_name = "rsem_poly_cohort","pbta-gene-expression-rsem-tpm-combat-seq-center.polya.rds")

# commenting out because there's only 1 cohort in stranded data sets
make_histology_pca_plots(dat_rsem_stranded, id_batch_histology, gene_id = dat_rsem_stranded$gene_id, report_name = "rsem_stranded_sequence", "pbta-gene-expression-rsem-tpm-combat-seq-center.stranded.rds", FALSE)
make_histology_pca_plots(dat_kallisto_polya, id_batch_histology, gene_id = dat_kallisto_polya$gene_id, report_name = "kallisto_polya_cohort", "pbta-gene-expression-kallisto-combat-seq-center.polya.rds")

# commenting out because there's only 1 cohort in stranded data sets
make_histology_pca_plots(dat_kallisto_stranded, id_batch_histology, gene_id = dat_kallisto_stranded$gene_id, report_name = "kallisto_stranded_sequence", "pbta-gene-expression-kallisto-combat-seq-center.stranded.rds", FALSE)

