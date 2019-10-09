# Process mutations for interaction plot from  MAF file.
#
# J Shapiro for ALSF - CCDL
#
# 2019
#
# Option descriptions
#
# --maf :  Relative file path to MAF file to be analyzed. Can be .gz compressed.
#          Assumes file path is given from top directory of 'OpenPBTA-analysis'.
# --metadata : Relative file path to MAF file to be analyzed. Can be .gz compressed.
#              Assumes file path is given from top directory of 'OpenPBTA-analysis'.
#
# Command line example:
#
# Rscript analyses/interaction-plots/01-cprocess-mutations.R \
#   --maf scratch/snv_dummy_data/strelka2 \
#   --metadata data/pbta-histologies.tsv \
#   --bed_wgs data/WGS.hg38.mutect2.unpadded.bed \
#   --bed_wxs data/WXS.hg38.100bp_padded.bed \
#   --annot_rds scratch/hg38_genomic_region_annotation.rds

#### Initial Set Up
# Establish base dir
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Magrittr pipe
`%>%` <- dplyr::`%>%`

# Load libraries:
library(ggplot2)

#### Set up options 
library(optparse)

option_list <- list(
  make_option(
    opt_str = "--maf", type = "character", 
    default = file.path("data", "pbta-snv-lancet.vep.maf.gz"),
    help = "Relative file path (assuming from top directory of
              'OpenPBTA-analysis') to MAF file to be analyzed. Can be .gz compressed.",
    metavar = "character"
  ),
  make_option(
    opt_str = "--metadata", type = "character", 
    default = file.path("data", "pbta-histologies.tsv"),
    help = "Relative file path (assuming from top directory of
              'OpenPBTA-analysis') to MAF file to be analyzed. Can be .gz compressed.",
    metavar = "character"
  ),
  make_option(
    opt_str = "--vaf", type = "numeric", 
    default = 0.2,
    help = "Minimum variant allele fraction",
    metavar = "character"
  )
)

# Parse options
opts <- parse_args(OptionParser(option_list = option_list))

# File locations

maf_file <- file.path(root_dir, opts$maf)
meta_file <- file.path(root_dir, opts$metadata)


#### Read files


maf_df <- data.table::fread(maf_file, data.table = FALSE)
meta_df <- data.table::fread(meta_file, data.table = FALSE)


### Reduce MAF to a smaller set of relevant columns and add vaf

maf_df <- maf_df %>% 
  dplyr::select(Hugo_Symbol, 
                Entrez_Gene_Id, 
                Chromosome, 
                Start_Position,
                End_Position, 
                Strand, 
                Variant_Classification,
                Variant_Type,
                Reference_Allele,
                Tumor_Seq_Allele1,
                Tumor_Seq_Allele2,
                Tumor_Sample_Barcode,
                t_depth, 
                t_ref_count,
                t_alt_count,
                Consequence,
  ) %>%
  dplyr::mutate(vaf = t_alt_count / (t_ref_count + t_alt_count)) 

samples <- unique(maf_df$Tumor_Sample_Barcode)
genes <- unique(maf_df$Hugo_Symbol)

#### Create gene by sample summary table 

gene_sample_counts <- maf_df %>%
  dplyr::filter(vaf > opts$vaf,
                Consequence != "intergenic_variant", 
                Entrez_Gene_Id > 0) %>% #filter unknowns
  dplyr::group_by(gene = Hugo_Symbol, sample = Tumor_Sample_Barcode) %>%
  dplyr::summarize(mutations = dplyr::n())


gene_counts <- gene_sample_counts %>%
  dplyr::group_by(gene) %>%
  dplyr::summarize(mutant_samples = dplyr::n(), 
                   total_muts = sum(mutations),
                   muts_per_sample = mean(mutations)
            )

# get top genes
top_count_genes <- gene_counts %>%
  dplyr::arrange(desc(mutant_samples),
                 desc(muts_per_sample)) %>%
  head(20)


#### Calculate co-occurence for top genes

gene_pairs <- t(combn(top_count_genes$gene, m = 2)) %>%
  tibble::as.tibble() %>%
  dplyr::rename(gene1 = V1, gene2 = V2) 

all_sample_counts <- expand.grid(gene = top_count_genes$gene,
                                 sample = samples, 
                                 stringsAsFactors = FALSE) %>%
  dplyr::left_join(gene_sample_counts) %>%
  tidyr::replace_na(list(mutations = 0)) %>%

gene_pair_counts <- gene_pairs %>% 
  dplyr::left_join(all_sample_counts,
                   by = c("gene1" = "gene")) %>%
  dplyr::rename(muts1 = mutations) %>%
  dplyr::left_join(all_sample_counts,
                   by = c("gene2" = "gene", 
                          "sample" = "sample")) %>%
  dplyr::rename(muts2 = mutations) %>%
  # use factors to preserve order later
  dplyr::mutate(gene1 = factor(gene1, levels = top_count_genes$gene),
                gene2 = factor(gene2, levels = top_count_genes$gene))

gene_pair_summary <- gene_pair_counts %>%
  dplyr::group_by(gene1, gene2) %>%
  dplyr::summarize(mut11 = sum(muts1 > 0  & muts2 > 0), 
                   mut10 = sum(muts1 > 0  & muts2 == 0), 
                   mut01 = sum(muts1 == 0 & muts2 > 0), 
                   mut00 = sum(muts1 == 0 & muts2 ==0))

### calculate fishers exact test

row_fisher <- function(w, x, y, z){
# function to calculate fisher test from row elements
  mat <- matrix(c(w, x, y, z), ncol = 2)
  fisher <- fisher.test(mat)
  return(fisher$p.value)
}

gene_pair_summary <- gene_pair_summary %>%
  dplyr::mutate(odds_ratio = (mut11 * mut00) / (mut10 * mut01),
                cooccur_sign = ifelse(odds_ratio > 1, 1, -1)) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(p = row_fisher(mut11, mut10, mut01, mut00), 
                cooccur_score = cooccur_sign * -log10(p))


### make plot
ggplot(gene_pair_summary, aes(x = gene1, y = gene2, fill = cooccur_score))+
  geom_tile(color = "white", size = 1) +
  xlab('') + 
  ylab('') +
  scale_x_discrete(position = "top") +
  scale_fill_distiller(type = "div", palette = 5, 
                       limits = c(-20, 20), 
                       oob = scales::squish) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = -45, hjust = 1), 
        axis.line = element_blank())

