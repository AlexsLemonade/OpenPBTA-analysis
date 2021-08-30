
set -e
set -o pipefail

# Use the OpenPBTA bucket as the default.
URL=${OPENPBTA_URL:-https://s3.amazonaws.com/d3b-openaccess-us-east-1-prd-pbta/data}
RELEASE=${OPENPBTA_RELEASE:-release-v21-20210820}
PREVIOUS=${OPENPBTA_RELEASE:-release-v20-20210726}

# Remove old symlinks in data
find data -type l -delete

# The md5sum file provides our single point of truth for which files are in a release.
curl --create-dirs $URL/$RELEASE/md5sum.txt -o data/$RELEASE/md5sum.txt -z data/$RELEASE/md5sum.txt

# Consider the filenames in the md5sum file and the release notes
FILES=(`tr -s ' ' < data/$RELEASE/md5sum.txt | cut -d ' ' -f 2` release-notes.md)

if [ -d "data/$PREVIOUS" ]
then
  # Find unchanged files
  echo "Checking for unchanged files..."
  cd data/$PREVIOUS
  UNCHANGED=(`md5sum -c ../$RELEASE/md5sum.txt 2>/dev/null | grep OK |cut -d ':' -f 1  || true`)
  echo $UNCHANGED
  cd ../../

  # Hard link unchanged files
  for oldfile in "${UNCHANGED[@]}"
  do
    if [ ! -e "data/$RELEASE/$oldfile" ]
    then
      echo "Hard linking $oldfile"
      ln data/$PREVIOUS/$oldfile data/$RELEASE/$oldfile
    fi
  done
fi

# Download the items in FILES if not already present
for file in "${FILES[@]}"
do
  if [ ! -e "data/$RELEASE/$file" ]
  then
    echo "Downloading $file"
    curl $URL/$RELEASE/$file -o data/$RELEASE/$file
  fi
done

# Download reference and gencode file from public ftp if it does not already exist
GENCODE="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.primary_assembly.annotation.gtf.gz"
cd data
if [ ! -e ${GENCODE##*/} ]
then
  echo "Downloading ${GENCODE##*/}"
  curl -O $GENCODE
fi

# if in CI, then we want to generate the reference FASTA from the BSgenome.Hsapiens.UCSC.hg38 R package
# because it is considerably faster to do so

if [ "$RELEASE" == "testing" ]; then
  Rscript -e "Biostrings::writeXStringSet(Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38, c('chr21', 'chr22', 'chrX', 'chrY')), 'GRCh38.primary_assembly.genome.fa.gz', format = 'fasta', compress = 'gzip')"
else
  REFERENCE="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/GRCh38.primary_assembly.genome.fa.gz"
  if [ ! -e ${REFERENCE##*/} ]
  then
    echo "Downloading ${REFERENCE##*/}"
    curl -O $REFERENCE
  fi
fi
cd -

# Check the md5s for everything we downloaded except CHANGELOG.md
cd data/$RELEASE
echo "Checking MD5 hashes..."
md5sum -c md5sum.txt
cd ../../

# Make symlinks in data/ to the files in the just downloaded release folder.
for file in "${FILES[@]}"
do
  ln -sfn $RELEASE/$file data/$file
done

# make data directory unwritable in CI
if [ "$RELEASE" == "testing" ]; then
  chmod u-w data
fi
