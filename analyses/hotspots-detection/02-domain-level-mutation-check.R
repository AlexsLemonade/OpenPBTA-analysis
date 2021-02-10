# hg38
library(EnsDb.Hsapiens.v86)
library(ensembldb)
library(tidyverse)
calls_base <- readRDS("OpenPBTA/analyses/hotspots-detection/results/combined_maf_subsets.RDS")


# to match geneID hg38
edb <- EnsDb.Hsapiens.v86
ids_to_match <- ensembldb::select(edb, keys=ensembldb::keys(edb, "GENEID"),
                                  column = c("GENENAME","GENEID","TXID"),
                                  keytype = "GENEID") 

txs <- getGeneRegionTrackForGviz(
  edb, filter = ~ protein_domain_id == "PF00069")

pdoms <- proteins(edb, filter = ~ tx_biotype == "protein_coding" &
                    gene_name %in% unique(calls_base$Hugo_Symbol) &
                    protein_domain_source == "pfam" &
                    protein_domain_id == "PF00069",
                  columns = c("protein_domain_id", "prot_dom_start",
                              "prot_dom_end"))
pdoms

pdoms_rng <- IRanges(start = pdoms$prot_dom_start, end = pdoms$prot_dom_end,
                     names = pdoms$protein_id)

pdoms_gnm <- proteinToGenome(pdoms_rng, edb)

# convert the list to a GRanges
pdoms_gnm_grng <- unlist(GRangesList(pdoms_gnm))

# check if cds_ok == FALSE remove these sites


