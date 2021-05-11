#!/bin/bash

# Chante Bethell for CCDL 2020
#
# Run the HGG molecular subtyping pipeline.
# Note: A local install of BEDOPS is required and can be installed using
# conda install -c bioconda bedops
# When OPENPBTA_SUBSET=1 (default), new HGG subset files will be generated.

set -e
set -o pipefail

# This option controls whether on not the step that generates the HGG only
# files gets run -- it will be turned off in CI
SUBSET=${OPENPBTA_SUBSET:-1}

# cds gencode bed file is used by other analyses where mutation data is
# filtered to only coding regions
exon_file="../../scratch/gencode.v27.primary_assembly.annotation.bed"

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

# Gather pathology diagnosis and pathology free text diagnosis for HGG sample selection
Rscript 00-HGG-select-pathology-dx.R

# Run the first script in this module that reclassifies high-grade gliomas
Rscript -e "rmarkdown::render('01-HGG-molecular-subtyping-defining-lesions.Rmd', clean = TRUE)"

# Run the second script in this module that subset files using the samples in the output
# file generated with `01-HGG-molecular-subtyping-defining-lesions.Rmd`.
if [ "$SUBSET" -gt "0" ]; then
  Rscript --vanilla 02-HGG-molecular-subtyping-subset-files.R
fi

#### Copy number data ----------------------------------------------------------

# Run the copy number data cleaning notebook
Rscript -e "rmarkdown::render('03-HGG-molecular-subtyping-cnv.Rmd', clean = TRUE)"

#### Mutation data -------------------------------------------------------------

# if the cds gencode bed file is not available from another analysis, generate
# it here
if [ ! -f "$exon_file" ]; then
  gunzip -c "../../data/gencode.v27.primary_assembly.annotation.gtf.gz" \
    | awk '$3 ~ /CDS/' \
    | convert2bed --do-not-sort --input=gtf - \
    > $exon_file
fi

# Run notebook that cleans the mutation data
Rscript -e "rmarkdown::render('04-HGG-molecular-subtyping-mutation.Rmd', clean = TRUE)"

#### Fusion data ---------------------------------------------------------------

# Run notebook that cleans the fusion data
Rscript -e "rmarkdown::render('05-HGG-molecular-subtyping-fusion.Rmd', clean = TRUE)"

#### Gene expression data ------------------------------------------------------

# Run notebook that cleans the gene expression data
Rscript -e "rmarkdown::render('06-HGG-molecular-subtyping-gene-expression.Rmd', clean = TRUE)"

#### Combine DNA data ----------------------------------------------------------

Rscript -e "rmarkdown::render('07-HGG-molecular-subtyping-combine-table.Rmd', clean = TRUE)"

#### 1p/19q co-deleted oligodendrogliomas notebook -----------------------------

Rscript -e "rmarkdown::render('08-1p19q-codeleted-oligodendrogliomas.Rmd', clean = TRUE)"

#### HGAT with `BRAF V600E` mutations clustering ------------------------------

# Run notebook that looks at how HGAT samples with `BRAF V600E` mutations cluster
Rscript -e "rmarkdown::render('09-HGG-with-braf-clustering.Rmd', clean = TRUE)"

# Add TP53 annotation
Rscript -e "rmarkdown::render('10-HGG-TP53-annotation.Rmd',clean=TRUE)"
