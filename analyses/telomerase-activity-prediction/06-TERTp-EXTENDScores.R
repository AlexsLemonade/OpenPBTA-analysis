# S Spielman for CCDL, 2023
# See https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/1651
# This script assesses whether there is a relationship between EXTEND scores and the presence of TERT SVs and/or TERT promoter (TERTp) mutations
# It additionally explores where samples with TERTp mutations lie in the TERC and TERT expression relationships with EXTEND scores
# From CIVIC we can get SNP ids of TERTp variants of interest:
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
input_dir <- file.path(analysis_dir, "input")
results_dir <- file.path(analysis_dir, "results")

extend_file_stranded <- file.path(analysis_dir, 
                         "results", 
                         "TelomeraseScores_PTBAStranded_FPKM.txt")
extend_file_polya <- file.path(analysis_dir, 
                         "results", 
                         "TelomeraseScores_PTBAPolya_FPKM.txt")

fpkm_file <- file.path(data_dir, 
                       "pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds")

metadata_file <- file.path(data_dir, 
                           "pbta-histologies.tsv")

maf_file <- file.path(data_dir, 
                      "pbta-snv-consensus-mutation.maf.tsv.gz")

hotspot_file <- file.path(data_dir, 
                          "pbta-snv-scavenged-hotspots.maf.tsv.gz")

sv_file <- file.path(data_dir, 
                  "pbta-sv-manta.tsv.gz")


# Output:
tert_alteration_file <- file.path(analysis_dir, 
                                  "results", 
                                  "extend_scores_tert_alterations.tsv")
tertp_plot_file <- file.path(analysis_dir, 
                             "plots", 
                             "TERTp_mutations.pdf")
tert_terc_plot_file <- file.path(analysis_dir, 
                                 "plots", 
                                 "TERTp_mutations_TERC_TERT_expression.pdf")


# Read in data --------------------

extend_stranded_df <- read_tsv(extend_file_stranded) %>%
  mutate(RNA_library = "stranded")
extend_polya_df <- read_tsv(extend_file_polya) %>%
  mutate(RNA_library = "polya")

extend_df <- extend_stranded_df %>%
  bind_rows(extend_polya_df)
  
fpkm_df <- read_rds(fpkm_file) %>%
  # wrangle into tidy (long) data frame
  as_tibble(rownames = "gene_symbol") %>%
  gather(Kids_First_Biospecimen_ID, 
         fpkm, 
         starts_with("BS_"))
  
metadata_df <- read_tsv(metadata_file, guess_max = 10000)
maf_df <- read_tsv(maf_file)
hotspot_df <- read_tsv(hotspot_file)

# TERT SV investigation -------------


# How many samples have TERT SVs and how many are pathogenic/likely pathogenic?
sv_df_tert <- read_tsv(sv_file) %>%
  # only keep PASS variants which are P/LP in dbVar
  filter(FILTER == "PASS",
         Gene.name == "TERT")

# three samples but none are P/LP
sv_df_tert %>%
  count(Kids.First.Biospecimen.ID.Tumor, dbVar_status)

tert_sv_bs_ids <- sv_df_tert %>%
  pull(Kids.First.Biospecimen.ID.Tumor) %>%
  unique()

# Get sample ids for TERT SV bs ids
tertsv_sample_ids <- metadata_df %>%
  filter(Kids_First_Biospecimen_ID %in% tert_sv_bs_ids) %>%
  pull(sample_id) %>%
  unique()

  
# Now get the RNASeq stranded IDs for those sample ids
# note that not all of these may be stranded!
rna_tertsv_bs_ids <- metadata_df %>%
  filter(sample_id %in% tertsv_sample_ids, 
         experimental_strategy == "RNA-Seq") %>%
  select(sample_id, Kids_First_Biospecimen_ID, RNA_library) 
# 2/3 sample ids have RNA-Seq. 1 polyA, 1 stranded.

# TERTp-positive samples and relationship to EXTEND scores ----------------------
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
         experimental_strategy == "RNA-Seq") %>%
  select(Kids_First_Biospecimen_ID, RNA_library) 
# ^ Only 6 of those samples were apparently stranded

# Compare back to EXTEND
extend_df <- extend_df %>% 
  mutate(tertp = ifelse(SampleID %in% rna_tertp_biospecimen_ids$Kids_First_Biospecimen_ID, 
                        "TERTp mutation present", 
                        "TERTp mutation not observed"),
         tertSV = ifelse(SampleID %in% rna_tertsv_bs_ids$Kids_First_Biospecimen_ID, 
                        "TERT SV present", 
                        "TERT SV not observed")) %>%
  arrange(tertp)

# Print table of samples with TERT SV or TERTp and extend scores for response to reviewers
extend_df %>%
  filter(tertp == "TERTp mutation present" | tertSV == "TERT SV present") %>%
  select(SampleID, tertp, tertSV, NormEXTENDScores, RNA_library) %>%
  write_tsv(tert_alteration_file)

# Subset to stranded only for plot
extend_df_stranded <- extend_df %>%
  filter(RNA_library == "stranded")

# Test is not significant:
# using nonparametric since we have a very small N=6 sample size for TERTp-positive
wilcox.test(NormEXTENDScores~tertp, data = extend_df_stranded)


# Visualize:
tertp_plot <- ggplot(extend_df_stranded) + 
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


# TERC and TERT ----------------------------------
fpkm_subset_df <- fpkm_df %>% 
  filter(gene_symbol %in% c("TERC", "TERT")) %>%
  inner_join(
    select(extend_df, 
           Kids_First_Biospecimen_ID = SampleID,
           NormEXTENDScores), 
    by = "Kids_First_Biospecimen_ID"
  ) %>%
  mutate(fpkm = log(fpkm + 1, 2),
         # add column for whether sample is TERTp
         tertp = ifelse(Kids_First_Biospecimen_ID %in% rna_tertp_biospecimen_ids$Kids_First_Biospecimen_ID, 
                        "TERTp mutation present", 
                        "TERTp mutation not observed"))

tert_terc_plot <- ggplot(fpkm_subset_df) + 
  aes(x = NormEXTENDScores, 
      y = fpkm) + 
  geom_point(aes(color = tertp, 
                 alpha = tertp)) + 
  geom_smooth(method = "lm") +
  scale_color_manual(values = c("grey70", "red")) +
  scale_alpha_manual(values = c(0.5, 1)) +
  facet_wrap(~gene_symbol) +
  labs(y = "log2(fpkm+1)")


ggsave(tert_terc_plot_file, 
       tert_terc_plot,
       width = 8, height = 5)
