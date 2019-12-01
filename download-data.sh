
set -e
set -o pipefail

# Use the OpenPBTA bucket as the default.
URL=${OPENPBTA_URL:-https://s3.amazonaws.com/kf-openaccess-us-east-1-prd-pbta/data}
RELEASE=${OPENPBTA_RELEASE:-release-v11-20191126}
PREVIOUS=${OPENPBTA_RELEASE:-release-v10-20191115}


# The md5sum file provides our single point of truth for which files are in a release.
curl --create-dirs $URL/$RELEASE/md5sum.txt -o data/$RELEASE/md5sum.txt -z data/$RELEASE/md5sum.txt

# Consider the filenames in the md5sum file + CHANGELOG.md
FILES=(`tr -s ' ' < data/$RELEASE/md5sum.txt | cut -d ' ' -f 2` CHANGELOG.md)

if [ -d "data/$PREVIOUS" ]
then
  # Find unchanged files
  echo "Checking for unchanged files..."
  cd data/$PREVIOUS
  UNCHANGED=(`md5sum -c ../$RELEASE/md5sum.txt --ignore-missing| grep OK |cut -d ':' -f 1  || true`)
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
  Rscript -e "BSgenome::export(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38, 'GRCh38.primary_assembly.genome.fa.gz', compress = 'gzip', format = 'fasta')"
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
md5sum -c md5sum.txt
cd ../../

# Remove old symlinks in data
find data -type l -delete

# Make symlinks in data/ to the files in the just downloaded release folder.
for file in "${FILES[@]}"
do
  ln -sfn $RELEASE/$file data/$file
done
