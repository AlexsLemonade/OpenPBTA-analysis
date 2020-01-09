#' @description  Generate fusion files specifically for consumption by molecular subtyping analyses
#' @author Daniel Miller <millerd15@@email.chop.edu> (D3b)
#' @note Date: January 2020

suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("optparse"))

#' **Filters**
#' 
#' *Fusions Filters*
#' 1: Exact match a list of fusions common in Ependymoma tumors
ependFuses = c(
  "C11orf95--MAML2",
  "C11orf95--RELA",
  "C11orf95--YAP1", 
  "LTBP3--RELA",
  "PTEN--TAS2R1",
  "YAP1--FAM118B",
  "YAP1--MAMLD1",
  "YAP1--MAMLD2"
)
ependGenes = c(
  "RELA"
)
#' 2: Exact match a list of fusions common in Embryonal tumors
#' as well as fusions containing a particular gene with any other gene
embryFuses = c(
  "CIC--NUTM1",
  "MN1--BEND2",
  "MN1--CXXC5"
)
embryGenes = c(
  "FOXR2",
  "MN1",
  "TTYH1"
)

#' Generate filtered fusion frame
#' @param df Unfiltered fusion data frame
#' @param bioid List of biospecimen IDs
#' @param fuses List of explicit fusion names
#' @param genes List of gene names
#' @return the filtered fusion data frame
filterFusion <- function(df, bioid, fuses, genes) {
  if (!missing(bioid)) {
    df <- filter(df, Sample %in% bioid)
  }
  if (!missing(fuses) & !missing(genes)) {
    df <- filter(df, FusionName %in% fuses | 
                   Gene1A %in% genes |
                   Gene2A %in% genes |
                   Gene1B %in% genes |
                   Gene2B %in% genes)
  } else if (!missing(fuses)) {
    df <- filter(df, FusionName %in% fuses)
  } else if (!missing(genes)) {
    df <- filter(df,     
                 Gene1A %in% genes |
                   Gene2A %in% genes |
                   Gene1B %in% genes |
                   Gene2B %in% genes)
  }
  return(df)
}

#' Creates a TSV of the filtered fusion sheet
#' @param df The filtered fusion data frame
#' @param fuses List of explicit fusion names
#' @param outputPath Path of the output file
#' @return Writes a TSV to the output path
generateOutput <- function(df, fuses, outputPath) {
  # create a list to be used to reduce the table to its minimal set
  # minimal set is the explicit fusions + any other non-specific fusions in the table
  lvls = unique(c(fuses, sort(as.character(unique(df$FusionName)))))
  df$FusionName = factor(df$FusionName, levels=lvls)
  # convert the data frame to a table
  tbl = table(df$Sample,df$FusionName)
  # convert back to a matrix
  mtx = as.data.frame.matrix(tbl)
  # convert the rownames into the first row
  mtx = setDT(mtx,keep.rownames=TRUE)
  # rename the column of the rownames
  colnames(mtx)[colnames(mtx) == 'rn'] <- 'Kids_First_Biospecimen_ID'
  write.table(mtx, file=outputPath, sep="\t", quote=FALSE, row.names=FALSE)
}

#' Set up the options
optionList <- list(
  make_option(
    opt_str = c("-d","--demographic_file"), type = "character",
    default = NULL, help = "Path to the demographic file."
  ),
  make_option(
    opt_str = c("-f","--fusions_file"), type = "character",
    default = NULL, help = "Path to the fusions file."
  ),
  make_option(
    opt_str = c("-o", "--output_dir"), type = "character",
    default = "pbta-consolidated-fusions", help = "Output directory for output files."
  )
)

#' Parse the options
opt <- parse_args(OptionParser(option_list = optionList))

#' Check the output directory
if (!dir.exists(opt$output_dir)) {
  dir.create(opt$output_dir)
}

#' Check that the files exist
if (!file.exists(opt$demographic_file)) {
  demo_missing <- paste("Error:", opt$demographic_file, "does not exist")
  if (!file.exists(opt$fusions_file)) {
    stop(paste(demo_missing, "\nError:", opt$fusions_file, "does not exist"))
  } else
    stop(demo_missing)
} else if (!file.exists(opt$fusions_file)) {
  stop(paste("Error:", opt$fusions_file, "does not exist"))
}
#' Load the files
demo <- read.csv(opt$demographic_file, sep="\t")
fuse <- read.csv(opt$fusions_file, sep="\t")

#' Filter the fusion files for your two populations
allFuseEpend <- filterFusion(df = fuse,
                          fuses = ependFuses,
                          genes = ependGenes)
allFuseEmbry <- filterFusion(df = fuse,
                          fuses = embryFuses,
                          genes = embryGenes)
#' Write the fusion frames to file
generateOutput(df  = allFuseEpend,
               fuses = ependFuses,
               outputPath = file.path(
                 opt$output_dir,"fusion_summary_ependymoma_foi.tsv"))
generateOutput(df = allFuseEmbry,
               fuses = embryFuses,
               outputPath = file.path(
                 opt$output_dir,"fusion_summary_embryonal_foi.tsv"))
