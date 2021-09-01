# Process mutations for interaction plot from  MAF file.
#
# JA Shapiro for ALSF - CCDL
#
# 2019
#
# Generates a table of gene by gene co-occurence data with p values from Fisher's exact test.
# By default, performs analysis of the top 50 most mutated genes.
#
# Option descriptions
#
# --maf :  File path to MAF file to be analyzed. Can be .gz compressed.
# --metadata : File path to metadata with sample information.
# --specimen_list: A file of specimens to include. Ideally, this list should consist
#       of independent samples (at most one from each individual). 
# --vaf: Minimum variant allele fraction of mutations to include.
# --min_depth: Minimum sequencing depth to call mutations.
# --min_mutated: Minimum number of mutated samples required to include a gene in 
#       the plot
# --max_genes: Maximum number of genes to plot interation data for 
#       (uses the most mutated n genes)
# --out: Output file location
# --disease_table: Location for disease table output
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
script_root <-
  file.path(root_dir, "analyses", "interaction-plots", "scripts")

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
    opt_str = "--maf",
    type = "character",
    default = file.path("..", "..", "data", "pbta-snv-consensus-mutation.maf.tsv.gz"),
    help = "File path of MAF file to be analyzed. Can be .gz compressed.",
    metavar = "character"
  ),
  make_option(
    opt_str = "--cnv",
    type = "character",
    default = file.path("..", "..", "data", "pbta-cnv-consensus.seg.gz"),
    help = "File path of CNV file to be analyzed in seg format. Can be .gz compressed.",
    metavar = "character"
  ),
  make_option(
    opt_str = "--metadata",
    type = "character",
    default = file.path("..", "..", "data", "pbta-histologies.tsv"),
    help = "File path of MAF file to be analyzed. Can be .gz compressed.",
    metavar = "character"
  ),
  make_option(
    opt_str = "--out",
    type = "character",
    default = file.path("results", "cooccurence.tsv"),
    help = "File path where output table will be placed.",
    metavar = "character"
  ),
  make_option(
    opt_str = "--disease_table",
    type = "character",
    default = NA,
    help = "File path where table of gene X disease mutation counts will be placed. (optional)",
    metavar = "character"
  ),
  make_option(
    opt_str = "--specimen_list",
    type = "character",
    default = file.path(
      "..", "..", "data",
      "independent-specimens.wgs.primary.tsv"
    ),
    help = "File path of specimen list file to be analyzed.
            A tsv file which must contain a column named 'Kids_First_Biospecimen_ID`",
    metavar = "character"
  ),
  make_option(
    opt_str = "--max_genes",
    type = "numeric",
    default = 50,
    help = "Maximum number of genes to include in figure. Will be filtered by number of
            mutations, so the n most mutated genes will have their co-occurence
            calculated and will appear in the resulting figure.",
    metavar = "character"
  ),
  make_option(
    opt_str = "--exclude_genes",
    type = "character",
    default = NA,
    help = "File path with a table of genes to be excluded from the figure.
            A tsv file which must contain a column named 'gene` that contains Hugo Symbols"
  ),
  make_option(
    opt_str = "--min_mutated",
    type = "numeric",
    default = 5,
    help = "Number of samples that must be mutated for a given gene to be included the 
            co-occurence calculations",
    metavar = "character"
  ),
  make_option(
    opt_str = "--include_syn",
    action = "store_true",
    default = FALSE,

    help = "Include synonymous coding mutations"
  ),
  make_option(
    opt_str = "--include_noncoding",
    action = "store_true",
    default = FALSE,
    help = "Include noncoding mutations (within transcript)"
  ),
  make_option(
    opt_str = "--include_nontranscribed",
    action = "store_true",
    default = FALSE,
    help = "Include nontranscribed (upstream & downstream) mutations"
  ),
  make_option(
    opt_str = "--vaf",
    type = "numeric",
    default = 0.05,
    help = "Minimum variant allele fraction to include",
    metavar = "numeric"
  ),
  make_option(
    opt_str = "--min_depth",
    type = "numeric",
    default = 0,
    help = "Minimum sequencing depth for called mutations",
    metavar = "numeric"
  )
)

# Parse options
opts <- parse_args(OptionParser(option_list = option_list))

# File locations

maf_file <- opts$maf
cnv_file <- opts$cnv
meta_file <- opts$metadata
out_file <- opts$out
disease_file <- opts$disease_table
if (!is.na(opts$specimen_list)) {
  specimen_file <- opts$specimen_list
}

if (!is.na(opts$exclude_genes)) {
  exclude_file <- opts$exclude_genes
}



#### Read files

maf_df <- data.table::fread(maf_file, data.table = FALSE)
cnv_df <- data.table::fread(cnv_file, data.table = FALSE)
meta_df <- data.table::fread(meta_file, data.table = FALSE)
if (exists("specimen_file")) {
  specimen_df <- data.table::fread(specimen_file, data.table = FALSE)
}
if (exists("exclude_file")) {
  exclude_df <- data.table::fread(exclude_file, data.table = FALSE)
}




### Reduce MAF to a smaller set of relevant columns

maf_df <- maf_df %>%
  dplyr::select(
    Hugo_Symbol,
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
    Consequence
  )

# get sample and gene lists
if (exists("specimen_df")) {
  samples <- specimen_df$Kids_First_Biospecimen_ID
} else {
  samples <- unique(maf_df$Tumor_Sample_Barcode)
}

genes <- unique(maf_df$Hugo_Symbol)
if (exists("exclude_df")){
  genes <- genes[!genes %in% exclude_df$gene]
}

# reduce metadata to only chosen samples
sample_meta <- meta_df %>%
  dplyr::filter(Kids_First_Biospecimen_ID %in% samples)

# reduce maf to chosen samples & calculate VAF
maf_df <- maf_df %>%
  dplyr::filter(Tumor_Sample_Barcode %in% samples) %>%
  dplyr::mutate(vaf = t_alt_count / (t_ref_count + t_alt_count))

# generate consequence lists and filter
intergenic <- c("IGR")

nontranscribed <- c(
  "3'Flank",
  "5'Flank",
  "Targeted_Region"
)

noncoding <- c(
  "RNA",
  "Intron",
  "3'UTR",
  "5'UTR",
  "Splice_Region",
  "lincRNA"
)

# Variant Classification with Low/Modifier variant consequences 
#  from maftools http://asia.ensembl.org/Help/Glossary?id=535
synonymous <- c(
  "Silent",
  "Start_Codon_Ins",
  "Start_Codon_SNP",
  "Stop_Codon_Del",
  "De_novo_Start_InFrame",
  "De_novo_Start_OutOfFrame"
)
# Variant Classification with High/Moderate variant consequences from maftools
nonsynonymous <- c(
  "Missense_Mutation",
  "Frame_Shift_Del",
  "In_Frame_Ins",
  "Frame_Shift_Ins",
  "Splice_Site",
  "Nonsense_Mutation",
  "In_Frame_Del",
  "Nonstop_Mutation",
  "Translation_Start_Site"
)


include <- nonsynonymous # always want nonsyn
if (opts$include_syn) {
  include <- c(include, synonymous)
}
if (opts$include_noncoding) {
  include <- c(include, noncoding)
}
if (opts$include_nontranscribed) {
  include <- c(include, nontranscribed)
}

maf_filtered <- maf_df %>%
  filter_mutations(
    min_vaf = opts$vaf,
    min_depth = opts$min_depth,
    include_var_class = include
  )

# count mutations by gene/sample pair
gene_sample_counts <- maf_filtered %>%
  dplyr::filter(Entrez_Gene_Id > 0, # remove unknowns
                Hugo_Symbol %in% genes) %>% # include only desired genes
  dplyr::group_by(gene = Hugo_Symbol, sample = Tumor_Sample_Barcode) %>%
  dplyr::tally(name = "mutations") %>%
  dplyr::ungroup()

# count # of samples mutated by gene
gene_counts <- gene_sample_counts %>%
  dplyr::group_by(gene) %>%
  dplyr::summarize(
    mutant_samples = dplyr::n(),
    total_muts = sum(mutations),
    mean_muts_per_sample = mean(mutations)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(
    desc(mutant_samples),
    desc(mean_muts_per_sample)
  ) %>%
  dplyr::filter(mutant_samples >= opts$min_mutated |
    dplyr::row_number() <= 2) # keep at least 2 genes


# get most often mutated genes
top_count_genes <- head(gene_counts, opts$max_genes)$gene

cooccur_summary <- coocurrence(gene_sample_counts, top_count_genes)

readr::write_tsv(cooccur_summary, out_file)

# only count genes if requested
if (is.na(disease_file)){
  quit()
}
# count mutated samples by disease types
gene_disease_counts <- gene_sample_counts %>%
  dplyr::filter(gene %in% top_count_genes) %>%
  dplyr::left_join(sample_meta, 
                   by = c("sample" = "Kids_First_Biospecimen_ID")) %>%
  dplyr::group_by(gene, disease = cancer_group) %>%
  dplyr::summarize(mutant_samples = dplyr::n(),
                   total_muts = sum(mutations),
                   mean_muts_per_sample = mean(mutations)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(
    desc(mutant_samples),
    desc(mean_muts_per_sample)
  )

readr::write_tsv(gene_disease_counts, disease_file)

