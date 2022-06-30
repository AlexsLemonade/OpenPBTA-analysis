#!/bin/sh

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
 use Cwd "abs_path";
 print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

# However in order to give Docker access to all the code we have to
# move up a level
cd ..

echo "Rebuilding the Docker image."
finished=1
attempts=0

# Use BuildKit
export DOCKER_BUILDKIT=1
# Simpler output for progress tracking
export BUILDKIT_PROGRESS=plain

while [ $finished != 0 ] && [ $attempts -lt 3 ]; do
    if [ $attempts -gt 0 ]; then
        echo "Failed to build Docker image, trying again."
    fi
    
    docker build \
           --secret id=gh_pat,env=GH_PAT \
           --tag "open-pbta" \
           --file "Dockerfile" .
    finished=$?
    attempts=$((attempts+1))
done

if [ $finished != 0 ] && [ $attempts -ge 3 ]; then
    echo "Could not build the Docker image after three attempts."
    exit 1
fi

env | grep "OPENPBTA_.*" > open_pbta_envs.txt

docker run \
       --env-file=open_pbta_envs.txt \
       --volume "$(pwd)":/rocker-build/ \
       -it "open-pbta" "$@"
