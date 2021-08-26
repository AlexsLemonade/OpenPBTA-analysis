
set -e
set -o pipefail

# Use the bucket that contains the CI data
URL=https://s3.amazonaws.com/d3b-openaccess-us-east-1-prd-pbta/data 
RELEASE=testing

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
REFERENCE="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/GRCh38.primary_assembly.genome.fa.gz"
cd data
curl -JO $GENCODE -z gencode.v27.primary_assembly.annotation.gtf.gz
curl -JO $REFERENCE -z GRCh38.primary_assembly.genome.fa.gz
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
