# Author: Krutika Gaonkar
#
# Read in concensus snv calls to gather alterations in TP53 and NF1
# to evaluate classifier
# @params snvConcensus multi-caller concensus snv calls
# @params clincalFile clinical file: pbta-histologies.tsv
# @params outputFolder output folder for alteration file
# @params gencode cds bed file from gencode

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("GenomicRanges"))

#### Source functions ----------------------------------------------------------
# We can use functions from the `snv-callers` module of the OpenPBTA project
# TODO: if a common util folder is established, use that instead
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "analyses", "snv-callers", "util",
                 "tmb_functions.R"))

#### Parse command line options ------------------------------------------------

option_list <- list(
  make_option(c("-s", "--snvConsensus"),type="character",
              help="Consensus snv calls (.tsv) "),
  make_option(c("-c","--clinicalFile"),type="character",
              help="clinical file for all samples (.tsv)"),
  make_option(c("-o","--outputFolder"),type="character",
              help="output folder for results "),
  make_option(c("-g","--gencode"),type="character",
              help="cds gencode bed file")
)

opt <- parse_args(OptionParser(option_list=option_list))
snvConsensusFile <- opt$snvConsensus
clinicalFile <- opt$clinicalFile
outputFolder <- opt$outputFolder
gencodeBed <- opt$gencode

#### Generate files with TP53, NF1 mutations -----------------------------------

# read in consensus SNV files
consensus_snv <- read_tsv(snvConsensusFile)
# gencode cds region BED file
gencode_cds <- read_tsv(gencodeBed, col_names = FALSE)
# clinical file
clinical <- read_tsv(clinicalFile)

# filter the MAF data.frame to only include entries that fall within the
# CDS bed file regions
coding_consensus_snv <- snv_ranges_filter(maf_df = consensus_snv,
                                          keep_ranges = gencode_cds)

# only include TP53 and NF1 and remove silent mutations and mutations in introns
tp53_nf1_coding <- coding_consensus_snv %>%
  dplyr::filter(Hugo_Symbol %in% c("TP53", "NF1")) %>%
  dplyr::filter(!(Variant_Classification %in% c("Silent", "Intron"))) %>%
  dplyr::arrange(dplyr::desc(Hugo_Symbol))

# include only the relevant columns from the MAF file
tp53_nf1_coding <- tp53_nf1_coding %>%
  dplyr::select(Chromosome, Start_Position, End_Position, Strand,
                Variant_Classification, Tumor_Sample_Barcode, Hugo_Symbol)

# biospecimen IDs for tumor or cell line DNA-seq
bs_ids <- clinical %>%
  dplyr::filter(sample_type != "Normal",
                experimental_strategy != "RNA-Seq") %>%
  dplyr::pull(Kids_First_Biospecimen_ID)

# all BS ids that are not in the data frame that contain the TP53 and NF1
# coding mutations should be labeled as not having either
bs_ids_without_mut <- setdiff(bs_ids,
                              unique(tp53_nf1_coding$Tumor_Sample_Barcode))

# create a data.frame with wildtype BS IDs for joining
without_mut_df <- data.frame(
  Chromosome = NA,
  Start_Position = NA,
  End_Position = NA,
  Strand = NA,
  Variant_Classification = NA,
  Tumor_Sample_Barcode = bs_ids_without_mut,
  Hugo_Symbol = "No_TP53_NF1_alt"
)

tp53_nf1_coding <- dplyr::bind_rows(tp53_nf1_coding,
                                    without_mut_df)

# save TP53 and NF1 SNV alterations
write_tsv(tp53_nf1_coding,
          file.path(outputFolder,"TP53_NF1_snv_alteration.tsv"))
