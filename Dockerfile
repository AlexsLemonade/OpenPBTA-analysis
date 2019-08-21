FROM rocker/tidyverse:3.6.0
MAINTAINER ccdl@alexslemonade.org
WORKDIR /rocker-build/

RUN apt-get update && apt-get install -y --no-install-recommends apt-utils

RUN apt-get install dialog apt-utils -y

# This is needed for coloblindr to install
RUN apt-get update -qq && apt-get -y --no-install-recommends install \
    libmagick++-dev 

# Need hexbin for making a hex plot 
RUN apt-get update -qq && apt-get -y --no-install-recommends install \
    && install2.r --error \
    --deps TRUE \
    hexbin \
    R.utils

# Use maftools for reading MAF files
RUN R -e "BiocManager::install(c('maftools'), update = FALSE)"

# Need this package to make plots colorblind friendly
RUN R -e "devtools::install_github('clauswilke/colorblindr', ref = '1ac3d4d62dad047b68bb66c06cee927a4517d678', dependencies = TRUE)"

#### Please install your dependencies here
#### Add a comment to indicate what analysis it is required for

