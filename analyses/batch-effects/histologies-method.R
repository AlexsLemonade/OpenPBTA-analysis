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

dat_rsem_polya = shorten(dat_rsem_polya)
dat_rsem_stranded = shorten(dat_rsem_stranded)
dat_kallisto_stranded = shorten(dat_kallisto_stranded)
dat_kallisto_polya = shorten(dat_kallisto_polya)
print("DATA HAS BEEN SHORTENED")

### coment out section above when running the full analysis ###


print("DOWNLOADED DATA")

dat_kallisto_polya = grouper(dat_kallisto_polya)
dat_kallisto_stranded = grouper(dat_kallisto_stranded)


covariate = filter(covariate, experimental_strategy == "RNA-Seq")

id_histology = cbind("Kids_First_Biospecimen_ID" = covariate$Kids_First_Biospecimen_ID, "histology" = covariate$short_histology)
colnames(id_histology) <- c("Kids_First_Biospecimen_ID", "histology")
id_histology = as_tibble(id_histology)

make_histology_pca_plots_methods(dat_rsem_polya, dat_rsem_stranded, id_histology, report_name = "rsem_method")
make_histology_pca_plots_methods(dat_kallisto_polya, dat_kallisto_stranded, id_histology, report_name = "kallisto_method")

