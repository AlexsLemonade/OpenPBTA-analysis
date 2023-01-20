# S Spielman for CCDL, 2023
# See https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/1651
# This script assesses whether there is a relationship between EXTEND scores and the presence of TERT promoter mutations
# From CIVIC we can get their SNP ids:
# C228T: https://civicdb.org/variants/248/summary
# C250T: https://civicdb.org/variants/4004/summary


# Libraries and paths -----------
library(tidyverse)
set.seed(2023) # set for jitter plot

# Input:
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, 
                          "analyses", 
                          "telomerase-activity-prediction")

extend_file <- file.path(analysis_dir, 
                         "results", 
                         "TelomeraseScores_PTBAStranded_FPKM.txt")

metadata_file <- file.path(data_dir, 
                           "pbta-histologies.tsv")

maf_file <- file.path(data_dir, 
                      "pbta-snv-consensus-mutation.maf.tsv.gz")

hotspot_file <- file.path(data_dir, 
                          "pbta-snv-scavenged-hotspots.maf.tsv.gz")


# Output:
tertp_plot_file <- file.path(analysis_dir, 
                             "plots", 
                             "TERTp_mutations.pdf")


# Read in data --------------------

extend_df <- read_tsv(extend_file)
metadata_df <- read_tsv(metadata_file, guess_max = 10000)
maf_df <- read_tsv(maf_file)
hotspot_df <- read_tsv(hotspot_file)


# Let's go find some variants ----------------------
snps <- c("rs1561215364", "rs1242535815")
maf_snp_bio_ids <- maf_df %>%
  filter(dbSNP_RS %in% snps) %>%
  pull(Tumor_Sample_Barcode) 

hotspot_snp_bio_ids <- hotspot_df %>%
  filter(dbSNP_RS %in% snps) %>%
  pull(Tumor_Sample_Barcode)
  
# There are 9 samples:
snp_biospecimen_ids <- unique(c(maf_snp_bio_ids, hotspot_snp_bio_ids))

# Get the corresponding sample ids
tertp_sample_ids <- metadata_df %>%
  filter(Kids_First_Biospecimen_ID %in% snp_biospecimen_ids) %>%
  pull(sample_id) 

# Now get the RNASeq stranded IDs for those sample ids
# note that not all of these may be stranded!
rna_tertp_biospecimen_ids <- metadata_df %>%
  filter(sample_id %in% tertp_sample_ids, 
         experimental_strategy == "RNA-Seq",
         RNA_library == "stranded") %>%
  pull(Kids_First_Biospecimen_ID) 
# ^ Only 6 of those samples were apparently stranded

# Compare back to EXTEND
extend_df <- extend_df %>% 
  mutate(tertp = ifelse(SampleID %in% rna_tertp_biospecimen_ids, 
                        "TERTp mutation present", 
                        "TERTp mutation not observed")) %>%
  arrange(tertp)


# Test is not significant:
# using nonparametric since we have a very small N=6 sample size for TERTp-positive
wilcox.test(NormEXTENDScores~tertp, data = extend_df)


# Visualize:
tertp_plot <- ggplot(extend_df) + 
  aes(x = "", y = NormEXTENDScores, 
      color = tertp, 
      fill = tertp, 
      size = tertp) + 
  geom_jitter(width = 0.15, shape = 21) + 
  scale_color_manual(values = c("gray70", "black")) + 
  scale_fill_manual(values = c("gray70", "red")) + 
  scale_size_manual(values = c(1, 2.5)) +
  labs(x = "", y = "Normalized EXTEND score") +
  theme_bw()

ggsave(tertp_plot_file, 
       tertp_plot,
       width = 6, height = 5)
