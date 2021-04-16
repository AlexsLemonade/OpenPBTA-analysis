#!/bin/bash

set -e
set -o pipefail

# This script checks that all python files in the docker image match 
# requirements.txt

# run from this file location but move up to root
cd "$(dirname "${BASH_SOURCE[0]}")"
cd ..

req_diff=/tmp/package_diffs.txt

## diff will exit code 1 with differences, so we need to pass true
pip3 freeze | diff requirements.txt - > $req_diff || true

# check if there are any differences in the file
if [ -s $req_diff ] 
then
  cat $req_diff && rm $req_diff
  echo "Python packages do not match requirements.txt, please check."
  exit 1
fi 

# if the diffs file was not produced for some reason, we should be sure to fail the same way
if [ ! -e $req_diff ]
then 
  pip3 freeze | diff requirements.txt -
fi

# clean up 
rm $req_diff
