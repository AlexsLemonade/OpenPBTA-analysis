
set -e
set -o pipefail

URL=${BASEURL:-https://s3.amazonaws.com/kf-openaccess-us-east-1-prd-pbta/data}
RELEASE=${REL:-release-v2-20190809}

FILES=(pbta_histologies.tsv pbta-cnv-cnvkit.seg.gz pbta-cnv-controlfreec.seg.gz pbta-fusion-arriba.tsv.gz pbta-fusion-starfusion.tsv.gz pbta-gene-expression-kallisto.rds pbta-snv-mutect2.vep.maf.gz pbta-snv-strelka2.vep.maf.gz pbta-sv-manta.tsv.gz md5sum.txt)

for file in "${FILES[@]}"
do
  curl --create-dirs $URL/$RELEASE/$file -o data/$RELEASE/$file -z data/$RELEASE/$file
done

cd data/$RELEASE
md5sum -c md5sum.txt
cd ../../

for file in "${FILES[@]}"
do
  ln -s data/$RELEASE/$file data/$file
done
