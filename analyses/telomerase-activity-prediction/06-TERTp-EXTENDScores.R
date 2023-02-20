# S Spielman for CCDL, 2023
# See https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/1651
# This script assesses whether there is a relationship between EXTEND scores and the presence of TERT promoter (TERTp) mutations
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

extend_file <- file.path(analysis_dir, 
                         "results", 
                         "TelomeraseScores_PTBAStranded_FPKM.txt")

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

sv_annot_file <- file.path(input_dir, 
                     "AnnotSV_VJFSh10MHo.tsv")


# Output:
tertp_plot_file <- file.path(analysis_dir, 
                             "plots", 
                             "TERTp_mutations.pdf")
tert_terc_plot_file <- file.path(analysis_dir, 
                                 "plots", 
                                 "TERTp_mutations_TERC_TERT_expression.pdf")


# Read in data --------------------

extend_df <- read_tsv(extend_file)
fpkm_df <- read_rds(fpkm_file) %>%
  # wrangle into tidy (long) data frame
  as_tibble(rownames = "gene_symbol") %>%
  gather(Kids_First_Biospecimen_ID, 
         fpkm, 
         starts_with("BS_"))
  
metadata_df <- read_tsv(metadata_file, guess_max = 10000)
maf_df <- read_tsv(maf_file)
hotspot_df <- read_tsv(hotspot_file)

# Read in and export TERT SV data in BED format for knotAnnotSV to assess pathogenicity of the SV call
# This BED file was put into https://lbgi.fr/AnnotSV/ and the TSV downloaded for ACMG assessment below
sv_df <- read_tsv(sv_file) %>%
  # only keep PASS variants
  filter(FILTER == "PASS",
         # let's check all 
         grepl("TERT", Gene.name)) %>%
  select(chrom = SV.chrom, chromStart = SV.start, chromEnd = SV.end, name = AnnotSV.ID, 
                  SVtype = SV.type, Sample = Kids.First.Biospecimen.ID.Tumor) %>%
  write_tsv(file.path(results_dir, "mantaSV_tert.bed"))

# Read in knotAnnotSV results
# knotAnnotSV will create one line for the "full" SV and one line for the SV "split" into each gene if the SV covers multiple genes.
# There were cases of SVs overlapping TERT (above in `sv_df`), but when split, none of these had an ACMG LP/P call.
tert_sv_bs_ids <- read_tsv(sv_annot_file) %>%
  # filter for ACMG class 4 or 5 (LP or P)
  filter((ACMG_class == 4 | ACMG_class == 5),
         Gene_name == "TERT") %>%
  pull(Samples_ID) %>%
  unique()

# Get sample ids for TERT SV bs ids
tertsv_sample_ids <- metadata_df %>%
  filter(Kids_First_Biospecimen_ID %in% tert_sv_bs_ids) %>%
  select(sample_id)
  
# Now get the RNASeq stranded IDs for those sample ids
# note that not all of these may be stranded!
rna_tertsv_bs_ids <- metadata_df %>%
  filter(sample_id %in% tertsv_sample_ids, 
         experimental_strategy == "RNA-Seq",
         RNA_library == "stranded") %>%
  pull(Kids_First_Biospecimen_ID) 
# These are both polyA(!), so we won't add to this plot

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
         tertp = ifelse(Kids_First_Biospecimen_ID %in% rna_tertp_biospecimen_ids, 
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
