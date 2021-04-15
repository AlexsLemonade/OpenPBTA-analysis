#!/bin/bash

set -e
set -o pipefail

# This script checks that all python files in the docker image match 
# requirements.txt

# run from this file location but move up to root
cd "$(dirname "${BASH_SOURCE[0]}")"
cd ..

diffs=$(pip3 freeze | diff requirements.txt - | wc -c)

if [ ! $diffs -eq 0 ] 
then
  echo "Python packages do not match requirements.txt, please check:" 
  pip3 freeze | diff requirements.txt  - 
  exit 1
fi 
  

