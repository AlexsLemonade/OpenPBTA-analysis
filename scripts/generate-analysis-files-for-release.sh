#!/bin/sh

set -e
set -o pipefail

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# If RUN_LOCAL is used, the time-intensive steps are skipped because they cannot
# be run on a local computer -- the idea is that setting RUN_LOCAL=1 will allow for
# local testing running/testing of all other steps
RUN_LOCAL=${RUN_LOCAL:-0}

# Get base directory of project
cd ..
BASEDIR="$(pwd)"
cd -

analyses_dir="$BASEDIR/analyses"
data_dir="$BASEDIR/data"
scratch_dir="$BASEDIR/scratch"

# Collapsed RNA-seq files
echo "Create collapse RSEM files"
bash ${analyses_dir}/collapse-rnaseq/run-collapse-rnaseq.sh

# Create the independent sample list using the *BASE* histology file
echo "Create independent sample list"
OPENPBTA_BASE_RELEASE=1 bash ${analyses_dir}/independent-samples/run-independent-samples.sh

# Fusion filtering
echo "Create fusion filtered list"
OPENPBTA_BASE_RELEASE=1 bash ${analyses_dir}/fusion_filtering/run_fusion_merged.sh

# Fusion summary
echo "Run fusion summary for subtypes"
bash ${analyses_dir}/fusion-summary/run-new-analysis.sh

# Compile all the files that need to be included in the release in one place
# in the scratch directory
compiled_dir=${scratch_dir}/analysis_files_for_release
mkdir -p ${compiled_dir}

# Copy over collapsed RNA-seq files
cp ${analyses_dir}/collapse-rnaseq/results/pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds ${compiled_dir}
cp ${analyses_dir}/collapse-rnaseq/results/pbta-gene-expression-rsem-fpkm-collapsed.polya.rds ${compiled_dir}

# Copy over independent specimen lists
cp ${analyses_dir}/independent-samples/results/independent-specimens.*  ${compiled_dir}

# Copy over fusions lists
cp ${analyses_dir}/fusion_filtering/results/pbta-fusion-putative-oncogenic.tsv ${compiled_dir}
cp ${analyses_dir}/fusion_filtering/results/pbta-fusion-recurrently-fused-genes-* ${compiled_dir}

# Copy over fusion summary
cp ${analyses_dir}/fusion-summary/results/* ${compiled_dir}

# Run modules that cannot be run locally due to memory requirements
if [ "$RUN_LOCAL" -lt "1" ]; then
  
  # Run SNV consensus & TMB step for PBTA data
  echo "Run SNV callers module for PBTA data"
  bash ${analyses_dir}/snv-callers/run_caller_consensus_analysis-pbta.sh
  
  # Run SNV consensus & TMB step for TCGA data
  echo "Run SNV callers module for TCGA data"  
  bash ${analyses_dir}/snv-callers/run_caller_consensus_analysis-tcga.sh
  
  # Copy over SNV callers
  ## PBTA
  cp ${analyses_dir}/snv-callers/results/consensus/pbta-snv-consensus-mutation.maf.tsv ${compiled_dir}
  cp ${analyses_dir}/snv-callers/results/consensus/pbta-snv-mutation-tmb-coding.tsv ${compiled_dir}
  cp ${analyses_dir}/snv-callers/results/consensus/pbta-snv-mutation-tmb-all.tsv ${compiled_dir}
  ## TCGA
  cp ${analyses_dir}/snv-callers/results/consensus/tcga-snv-consensus-snv.maf.tsv
  cp ${analyses_dir}/snv-callers/results/consensus/tcga-snv-mutation-tmb-coding.tsv ${compiled_dir}
  cp ${analyses_dir}/snv-callers/results/consensus/tcga-snv-mutation-tmb-all.tsv ${compiled_dir}
  
  # Run hotspot detection
  echo "Run hotspots detection"
  bash ${analyses_dir}/hotspots-detection/run_overlaps_hotspots.sh
  
  # Copy over hotspots detection
  cp ${analyses_dir}/hotspots-detection/results/pbta-snv-scavenged-hotspot.maf.tsv.gz ${compiled_dir}
  
  # Run consensus CN caller step
  echo "Run CNV consensus"
  bash ${analyses_dir}/copy_number_consensus_call/run_consensus_call.sh
  
  # Copy over CNV consensus
  cp ${analyses_dir}/copy_number_consensus_call/results/pbta-cnv-consensus.seg.gz ${compiled_dir}
  
  # Run GISTIC step -- only the part that generates ZIP file
  echo "Run GISTIC"
  # Run a step that subs ploidy for NA to allow GISTIC to run
  Rscript ${analyses_dir}/run-gistic/scripts/prepare_seg_for_gistic.R \
    --in_consensus ${analyses_dir}/copy_number_consensus_call/results/pbta-cnv-consensus.seg.gz \
    --out_consensus ${analyses_dir}/run-gistic/results/pbta-cnv-consensus-gistic-only.seg.gz \
    --histology ${data_dir}/pbta-histologies-base.tsv
    
  # This will use the file that just got generated above
  bash ${analyses_dir}/run-gistic/scripts/run-gistic-openpbta.sh
  
  # Copy over GISTIC
  cp ${analyses_dir}/run-gistic/results/pbta-cnv-consensus-gistic.zep ${compiled_dir}
  
  # Run step that generates "most focal CN" files (annotation) using the *BASE* histology file
  echo "Run focal CN file preparation"
  OPENPBTA_BASE_RELEASE=1 bash ${analyses_dir}/focal-cn-file-preparation/run-prepare-cn.sh

  # Copy over focal CN
  cp ${analyses_dir}/focal-cn-file-preparation/results/consensus_seg_annotated_cn_autosomes.tsv.gz ${compiled_dir}
  cp ${analyses_dir}/focal-cn-file-preparation/results/consensus_seg_annotated_cn_x_and_y.tsv.gz ${compiled_dir}
  
fi


