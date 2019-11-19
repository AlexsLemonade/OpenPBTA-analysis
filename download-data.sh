
set -e
set -o pipefail

# Use the OpenPBTA bucket as the default.
URL=${OPENPBTA_URL:-https://s3.amazonaws.com/kf-openaccess-us-east-1-prd-pbta/data}
RELEASE=${OPENPBTA_RELEASE:-release-v10-20191115}

# Remove symlinks in data
find data -type l -delete

# The md5sum file provides our single point of truth for which files are in a release.
curl --create-dirs $URL/$RELEASE/md5sum.txt -o data/$RELEASE/md5sum.txt -z data/$RELEASE/md5sum.txt

# Consider the filenames in the md5sum file + CHANGELOG.md
FILES=(`tr -s ' ' < data/$RELEASE/md5sum.txt | cut -d ' ' -f 2` CHANGELOG.md)

# Download the items in FILES if newer than what's on server
for file in "${FILES[@]}"
do
  curl --create-dirs $URL/$RELEASE/$file -o data/$RELEASE/$file -z data/$RELEASE/$file
done

# Download reference and gencode file from public ftp
GENCODE="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.primary_assembly.annotation.gtf.gz"
cd data
curl -JO $GENCODE

# if in CI, then we want to generate the reference FASTA from the BSgenome.Hsapiens.UCSC.hg38 R package
# because it is considerably faster to do so

if [ "$RELEASE" == "testing" ]; then
  Rscript -e "BSgenome::export(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38, 'GRCh38.primary_assembly.genome.fa.gz', compress = 'gzip')"
else
  REFERENCE="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/GRCh38.primary_assembly.genome.fa.gz"
  curl -JO $REFERENCE
fi
cd -

# Check the md5s for everything we downloaded except CHANGELOG.md
cd data/$RELEASE
md5sum -c md5sum.txt
cd ../../

# Make symlinks in data/ to the files in the just downloaded release folder.
for file in "${FILES[@]}"
do
  ln -sfn $RELEASE/$file data/$file
done

# Unzip any zip files in the data directory using the update flag
unzip -u -d data data/*.zip 
