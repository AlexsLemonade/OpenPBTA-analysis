# In this script we will find CNV losses that overlap with TP53 domains:
#  
# TAD = trans-activating domain (essential for function)
# DBD = DNA-binding domain (residues 102–292)
# TD = tetramerization domain (residues 326–356)
#
# We want to subset CNV calls where the domain are lost which will 
# possibly	lead to loss of function to use for evaluation of TP53 
# inactivation score at a later step.
# Krutika Gaonkar for D3b
# Nov 13th 2020

### Setup
library("tidyverse")
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
# File path to subset directory
subset_dir <-
  file.path(root_dir, "analyses", "tp53_nf1_score", "tp53-subset")
# File path to results directory
results_dir <-
  file.path(root_dir, "analyses", "tp53_nf1_score", "results")
# data folder
data_dir <-
  file.path(root_dir,"data")

if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}


# Input files

# `02-add-ploidy-consensus.Rmd` in focal-cn-file-preparation adds ploidy information to the consensus SEG file and adds a status column that defines gain and loss broadly.

#```
# Rscript -e "rmarkdown::render('analyses/focal-cn-file-preparation/02-add-ploidy-consensus.Rmd', clean = TRUE)"
#```

# TO-DO consensus_seg_with_status.tsv should go into results in focal-cn-file-preparation so we can access the file in this notebook 

# consensus seg for file for the location of CNVs
consensus_seg <- read_tsv(file.path(root_dir, 
                                    "scratch",
                                    "consensus_seg_with_status.tsv"))

# Gene location and domain overlap file 
bioMartDataPfamTp53 <- 
  readRDS(system.file("extdata", "pfamDataBioMart.RDS", package = "annoFuse")) %>%
  dplyr::filter(hgnc_symbol=="TP53")

# Genomic range for cnv seg file 
cnv_gr <- consensus_seg %>%
  dplyr::rename(chr = chrom, start = loc.start, end = loc.end,
                copy_number = copy.num) %>%
  dplyr::select(-num.mark, -seg.mean) %>%
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE,
                                          starts.in.df.are.0based = FALSE)

# Genomic range for gene location and domain overlap file 
domain_gr <- bioMartDataPfamTp53 %>%
  dplyr::filter(!is.na(domain_start),!is.na(domain_end)) %>%
  dplyr::mutate(strand = if_else(strand=="-1","-","+"),
                chromosome_name = paste0("chr",chromosome_name)) %>%
  dplyr::rename(chr = chromosome_name, start = domain_start, end = domain_end) %>%
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE,
                                          starts.in.df.are.0based = FALSE)

# overlap cnv and domain
overlaps <- IRanges::mergeByOverlaps(cnv_gr, domain_gr)

# get CNV and domain overlap per BS id
annotated_cn <- data.frame(
  biospecimen_id = overlaps$Kids_First_Biospecimen_ID,
  status = overlaps$status,
  copy_number = overlaps$copy_number,
  ploidy = overlaps$tumor_ploidy,
  hgnc_symbol = overlaps$hgnc_symbol,
  pfam_id = overlaps$pfam_id,
  NAME = overlaps$NAME,
  stringsAsFactors = FALSE
) %>%
  dplyr::distinct() %>%
  # select loss that overlaps the TP53 core domains 
  dplyr::filter(status=="loss") %>%
  write_tsv(file.path(subset_dir,"loss_overlap_domains_tp53.tsv"))



