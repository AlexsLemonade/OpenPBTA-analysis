library(tidyverse)
library(BatchQC)
library(sva)



# START HERE WHILE READING CODE COMMENTS

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

print("GROUPED DATA")

# See comments on make_id_batch
id_batch = make_id_batch(covariate, "cohort")

# Run batchQC and combat. This will run batchQC before and after running combat. It will also perform join operations between id_batch and our gene expression files 
run_batchQC(dat_rsem_polya, id_batch, gene_id = dat_rsem_polya$gene_id, report_name = "rsem_polya_cohort_report.html","pbta-gene-expression-rsem-tpm-combat-cohort.polya.rds", TRUE)
run_batchQC(dat_rsem_stranded, id_batch, gene_id = dat_rsem_stranded$gene_id, report_name = "rsem_stranded_cohort_report.html", "pbta-gene-expression-rsem-tpm-combat-cohort.stranded.rds", run_combat =  FALSE)
run_batchQC(dat_kallisto_polya, id_batch, gene_id = dat_kallisto_polya$gene_id, report_name = "kallisto_polya_cohort_report.html", "pbta-gene-expression-kallisto-combat-cohort.polya.rds", TRUE)
run_batchQC(dat_kallisto_stranded, id_batch, gene_id = dat_kallisto_stranded$gene_id, report_name = "kallisto_stranded_cohort_report.html", "pbta-gene-expression-kallisto-combat-cohort.stranded.rds", run_combat = FALSE)
