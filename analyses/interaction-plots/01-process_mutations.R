# Process mutations for interaction plot from  MAF file.
#
# J Shapiro for ALSF - CCDL
#
# 2019
#
# Option descriptions
#
# --maf :  Relative file path to MAF file to be analyzed. Can be .gz compressed.
#       File path is given from top directory of 'OpenPBTA-analysis'.
# --metadata : Relative file path to MAF file to be analyzed. Can be .gz compressed.
#       File path is given from top directory of 'OpenPBTA-analysis'.
# --specimen_list: A file of specimens to include. Ideally, this list should consist
#       of independent samples (at most one from each individual). File path is given
#       from the top directory of 'OpenPBTA-analysis'.
# --vaf: Minimum Variant allele fraction of mutations to include.
#                 
#
# Command line example:
#
# Rscript analyses/interaction-plots/01-process-mutations.R \
#   --maf data/pbta-snv-lancet.vep.maf.gz
#   --metadata data/pbta-histologies.tsv 
#   --specimen_list analysis/independent-samples/results/independent-specimens.wgs.primary.tsv

#### Initial Set Up
# Establish base dir
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
script_root <- file.path(root_dir, "analyses", "interaction-plots")

# Magrittr pipe
`%>%` <- dplyr::`%>%`

# Load libraries:
library(optparse)
library(ggplot2)


# Load functions
source(file.path(script_root, "cooccur_functions.R"))

#### Set up options 
option_list <- list(
  make_option(
    opt_str = "--maf", type = "character", 
    default = file.path("data", "pbta-snv-lancet.vep.maf.gz"),
    help = "Relative file path (from top directory of 'OpenPBTA-analysis') 
              to MAF file to be analyzed. Can be .gz compressed.",
    metavar = "character"
  ),
  make_option(
    opt_str = "--seg", type = "character", 
    default = file.path("data", "pbta-cnv-cnvkit.seg.gz"),
    help = "Relative file path (from top directory of 'OpenPBTA-analysis') 
            to CNV file to be analyzed in seg format. Can be .gz compressed.",
    metavar = "character"
  ),
  make_option(
    opt_str = "--metadata", type = "character", 
    default = file.path("data", "pbta-histologies.tsv"),
    help = "Relative file path (from top directory of 'OpenPBTA-analysis')
            to MAF file to be analyzed. Can be .gz compressed.",
    metavar = "character"
  ),
  make_option(
    opt_str = "--specimen_list", type = "character", 
    default = NA,
    help = "Relative file path (from top directory of 'OpenPBTA-analysis')
            to MAF file to be analyzed. Can be .gz compressed.",
    metavar = "character"
  ),
  make_option(
    opt_str = "--include_syn", action = "store_true", 
    default = FALSE,
    help = "Include synonymous mutations"
  ),
  make_option(
    opt_str = "--vaf", type = "numeric", 
    default = 0.2,
    help = "Minimum variant allele fraction to include",
    metavar = "numeric"
  )
)

# Parse options
opts <- parse_args(OptionParser(option_list = option_list))

# File locations

maf_file <- file.path(root_dir, opts$maf)
seg_file <- file.append(root_dir, opts$seg)
meta_file <- file.path(root_dir, opts$metadata)
if(!is.na(opts.specimen_list)){
  specimen_file <- file.path(root_dir, opts$specimen_list)
}



#### Read files


maf_df <- data.table::fread(maf_file, data.table = FALSE)
seg_df <- data.table::fread(seg_file, data.table = FALSE)
meta_df <- data.table::fread(meta_file, data.table = FALSE)
if(exists(specimen_file)){
  specimen_df <- data.table::fread(specimen_file, data.table = FALSE)
}




### Reduce MAF to a smaller set of relevant columns

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
  )  

# get sample and gene lists
if (exits(specimen_df)){
  samples <- specimen_df$Kids_First_Biospecimen_ID
} else {
  samples <- unique(maf_df$Tumor_Sample_Barcode)
}
genes <- unique(maf_df$Hugo_Symbol)

# reduce metadata to only chosen samples
sample_meta <- meta_df %>%
  dplyr::filter(Kids_First_Biospecimen_ID %in% samples)

# reduce maf to chosen samples & calculat VAF
sample_maf_df <- maf_df %>%
  dplyr::filter(Tumor_Sample_Barcode %in% samples) %>%
  dplyr::mutate(vaf = t_alt_count / (t_ref_count + t_alt_count))

# generate consequence lists and filter
exclusions <- c("intergenic_variant")
if (!opts$include_syn){exclusions <- c(exclusions, "synonymous_variant")}
if (!opts$include_intron){exclusions <- c(exclusions, "intron_variant")}
sample_maf_filtered <- sample_maf %>%
  filter_mutations(min_vaf = opts$vaf,
                   exclude_consequence = exclusions)

# count mutations by gene/sample pair and syn
gene_sample_counts <- sample_maf_df %>%
  dplyr::filter(Entrez_Gene_Id > 0) %>% #filter unknowns
  dplyr::group_by(gene = Hugo_Symbol, sample = Tumor_Sample_Barcode) %>%
  dplyr::tally(name = "mutations")

# count # of samples mutated by gene
gene_counts <- gene_sample_counts %>%
  dplyr::filter(sample %in% samples) %>%
  dplyr::group_by(gene) %>%
  dplyr::summarize(mutant_samples = dplyr::n(), 
                   total_muts = sum(mutations),
                   muts_per_sample = mean(mutations)
                   ) %>%
  dplyr::arrange(desc(mutant_samples),
                 desc(muts_per_sample))

# get most often mutated genes
top_count_genes <- head(gene_counts, 50)


gene_pair_summary <- coocurrence(gene_sample_counts, top_count_genes)

 
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
        axis.line = element_blank(), 
        axis.ticks = element_blank(),
        legend.justification=c(1,0), 
        legend.position=c(1,0),
        legend.key.size = unit(2, "char"))

