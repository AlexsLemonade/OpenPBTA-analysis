library(tidyverse)
library(BatchQC)
library(sva)
library(stringr)

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

dat_rsem_polya = shorten(dat_rsem_polya)
dat_rsem_stranded = shorten(dat_rsem_stranded)
dat_kallisto_stranded = shorten(dat_kallisto_stranded)
dat_kallisto_polya = shorten(dat_kallisto_polya)
print("DATA HAS BEEN SHORTENED")

### coment out section above when running the full analysis ###


# summarize transcript data at the gene level
dat_kallisto_polya = grouper(dat_kallisto_polya)
dat_kallisto_stranded = grouper(dat_kallisto_stranded)

print("GROUPED TRANSCRIPTS")


("PREPARING TO RUN ANALYSIS")

run_batchQC_polyA_vs_stranded(dat_rsem_polya, dat_rsem_stranded, report_name = "rsem_method_report.html","pbta-gene-expression-rsem-tpm.rds")
run_batchQC_polyA_vs_stranded(dat_kallisto_polya, dat_kallisto_stranded, report_name = "kallisto_method_report.html","pbta-gene-expression-kallisto.rds")

