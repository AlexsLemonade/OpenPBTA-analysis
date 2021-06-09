#!/bin/bash
# Module author: Komal S. Rathi
# Shell script author: Jaclyn Taroni for ALSF CCDL
# 2019

# This script runs the steps for generating collapsed RNA-seq matrices
# and analyzing the correlation levels of multi-mapped Ensembl genes

set -e
set -o pipefail

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

# create results directory if it doesn't already exist
mkdir -p results

# Usage: project acronym to use as prefix for input out files 
usage(){ echo "Usage: $0 [-h] [-x <expression or counts>]  [-q <quantification type>] " 1>&2; exit 1; }

while getopts ":hq:x:" opt; do
    case "${opt}" in
	h)
	    usage
	    ;;
        x)
            expr_count=$OPTARG
            ;;
        q)
            quantificationType=$OPTARG
            ;;
	:)
	    printf "missing argument for -%s\n" "$OPTARG" 1>&2
	    usage
	    ;;
	\?)
	    printf "illegal option: -%s\n" "$OPTARG" 1>&2
	    usage
	    ;;
    esac
done
shift $((OPTIND - 1))


if [ -z "${quantificationType}" ]; then
    usage
fi

if [ -x "${expr_count}" ]; then
    usage
fi

# generate collapsed matrices for poly-A and stranded datasets


Rscript --vanilla 01-summarize_matrices.R \
  -i ../../data/gene-${expr_count}-rsem-${quantificationType}.rds \
  -g ../../data/gencode.v27.primary_assembly.annotation.gtf.gz \
  -m results/gene-${expr_count}-rsem-${quantificationType}-collapsed.rds \
  -t results/gene-${expr_count}-rsem-${quantificationType}-collapsed_table.rds


  # run the notebook for analysis of dropped genes
Rscript -e "rmarkdown::render(input = '02-analyze-drops.Rmd', output_file = paste0('02-analyze-drops','-${quantificationType}'),params = list(annot.table = 'results/gene-${expr_count}-rsem-${quantificationType}-collapsed_table.rds'), clean = TRUE)"


