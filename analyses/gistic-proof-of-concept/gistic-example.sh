#!/bin/bash

set -e 
set -o pipefail

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/mcr/v83/runtime/glnxa64:/opt/mcr/v83/bin/glnxa64:/opt/mcr/v83/sys/os/glnxa64:
export XAPPLRESDIR=/opt/mcr/v83/X11/app-defaults

cd /home/rstudio/gistic_install && ./run_gistic_example
