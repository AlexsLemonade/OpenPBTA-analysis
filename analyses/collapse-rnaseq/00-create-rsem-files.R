# Author: Komal S. Rathi
# Date: 11/09/2019
# Function: Merge RSEM files and split into polya and stranded

# load libraries
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(reshape2))

# read params
option_list <- list(
  make_option(c("-i", "--inputdir"), type = "character",
              help = "Input directory for RSEM files"),
  make_option(c("-c", "--clinical"), type = "character",
              help = "Histology file (.TSV)"),
  make_option(c("-m", "--manifest"), type = "character",
              help = "Manifest file (.csv)"),
  make_option(c("-o", "--outdir"), type = "character",
              help = "Path to output directory")
)

# parse the parameters
opt <- parse_args(OptionParser(option_list = option_list))
topDir <- opt$inputdir
clin <- opt$clinical
manifest <- opt$manifest
outdir <- opt$outdir

# read manifest file
manifest <- read.csv(manifest, stringsAsFactors = F)
manifest <- manifest[,c("Kids.First.Biospecimen.ID","name")]
manifest$name <- gsub('[.].*', '', manifest$name)

# read histology file and split into polyA and stranded
clin <- read.delim(clin, stringsAsFactors = F)
polya <- clin %>% 
  filter(experimental_strategy == "RNA-Seq" & RNA_library == "poly-A") %>%
  as.data.frame()

stranded <- clin %>% 
  filter(experimental_strategy == "RNA-Seq" & RNA_library == "stranded") %>%
  as.data.frame()

# read and merge RSEM genes files
lfiles <- list.files(path = topDir, pattern = "*.rsem.genes.results.gz", recursive = TRUE, full.names = T)
read.rsem <- function(x){
  print(x)
  dat <- data.table::fread(x)
  filename <- gsub('.*[/]|.rsem.genes.results.gz','',x)
  sample.id <- manifest[which(manifest$name %in% filename),'Kids.First.Biospecimen.ID']
  dat$Sample <- sample.id
  return(dat)
}
expr <- lapply(lfiles, read.rsem)
expr <- data.table::rbindlist(expr)
expr.fpkm <- dcast(expr, gene_id~Sample, value.var = 'FPKM')  # FPKM
expr.counts <- dcast(expr, gene_id~Sample, value.var = 'expected_count')  # counts

# split into polya and stranded matrices
polya.fpkm <- expr.fpkm[,polya$Kids_First_Biospecimen_ID]
stranded.fpkm <- expr.fpkm[,stranded$Kids_First_Biospecimen_ID]
polya.counts <- expr.counts[,polya$Kids_First_Biospecimen_ID]
stranded.counts <- expr.counts[,stranded$Kids_First_Biospecimen_ID]

# save output
saveRDS(polya.fpkm, file = paste0(outdir,'/pbta-gene-expression-rsem-fpkm.polya.rds'))
saveRDS(polya.counts, file = paste0(outdir, '/pbta-gene-counts-rsem-expected_count.polya.rds'))
saveRDS(stranded.fpkm, file = paste0(outdir, '/pbta-gene-expression-rsem-fpkm.stranded.rds'))
saveRDS(stranded.counts, file = paste0(outdir, '/pbta-gene-counts-rsem-expected_count.stranded.rds'))
print("Done!")

