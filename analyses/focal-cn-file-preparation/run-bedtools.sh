# C. Bethell for CCDL 2020
# Run focal-cn-file-preparation module
#
# Usage: bash run-bedtools.sh

scratch_dir=../../scratch
callable_intersect_with_cytoband_file=${scratch_dir}/cytoband_status/intersect_with_cytoband_callable.bed
loss_intersect_with_cytoband_file=${scratch_dir}/cytoband_status/intersect_with_cytoband_losses.bed
gain_intersect_with_cytoband_file=${scratch_dir}/cytoband_status/intersect_with_cytoband_gains.bed

# Download and save UCSC cytoband file as bed file
wget -O ${scratch_dir}/ucsc_cytoband.bed http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBand.txt.gz

# Use bedtools coverage to find the intersection between the UCSC file with
# cytoband data and the `consensus_status` bed files prepared in
# `02-add-ploidy-consensus.Rmd`
for bedfile in ${scratch_dir}/cytoband_status/consensus_status.*.bed; do
  file_id=$(basename -s .bed $bedfile)
  sample_id=${file_id#consensus_status.}
  bedtools coverage \
    -a ${scratch_dir}/ucsc_cytoband.bed \
    -b ${bedfile} \
    -sorted \
    | sed "s/$/\t${sample_id}/" \
    > ${scratch_dir}/cytoband_status/${file_id}.coverage.bed
done

cat ${scratch_dir}/cytoband_status/consensus_status.*.coverage.bed \
 > $callable_intersect_with_cytoband_file

# Use bedtools coverage to find the intersection between the UCSC file with
# cytoband data and the `consensus_loss_status` bed files prepared in
# `02-add-ploidy-consensus.Rmd`
for bedfile in ${scratch_dir}/cytoband_status/consensus_loss_status.*.bed; do
  file_id=$(basename -s .bed $bedfile)
  sample_id=${file_id#consensus_loss_status.}
  bedtools coverage \
    -a ${scratch_dir}/ucsc_cytoband.bed \
    -b ${bedfile} \
    -sorted \
    | sed "s/$/\t${sample_id}/" \
    > ${scratch_dir}/cytoband_status/${file_id}.coverage.bed
done

cat ${scratch_dir}/cytoband_status/consensus_loss_status.*.coverage.bed \
 > $loss_intersect_with_cytoband_file
 
 # Use bedtools coverage to find the intersection between the UCSC file with
# cytoband data and the `consensus_gain_status` bed files prepared in
# `02-add-ploidy-consensus.Rmd`
for bedfile in ${scratch_dir}/cytoband_status/consensus_gain_status.*.bed; do
  file_id=$(basename -s .bed $bedfile)
  sample_id=${file_id#consensus_gain_status.}
  bedtools coverage \
    -a ${scratch_dir}/ucsc_cytoband.bed \
    -b ${bedfile} \
    -sorted \
    | sed "s/$/\t${sample_id}/" \
    > ${scratch_dir}/cytoband_status/${file_id}.coverage.bed
done

cat ${scratch_dir}/cytoband_status/consensus_gain_status.*.coverage.bed \
 > $gain_intersect_with_cytoband_file
 