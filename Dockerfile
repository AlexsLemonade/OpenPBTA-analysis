FROM rocker/tidyverse:3.6.0
MAINTAINER ccdl@alexslemonade.org
WORKDIR /rocker-build/

RUN apt-get update && apt-get install -y --no-install-recommends apt-utils

RUN apt-get install dialog apt-utils -y

# Required for installing mapview for interactive sample distribution plots
# libmagick++-dev is needed for coloblindr to install
RUN apt-get update -qq && apt-get -y --no-install-recommends install \
    libgdal-dev \
    libudunits2-dev \
    libmagick++-dev

# Required forinteractive sample distribution plots
# map view is needed to create HTML outputs of the interactive plots
RUN apt-get update -qq && apt-get -y --no-install-recommends install \    
    && install2.r --error \
    --deps TRUE \
    gdalUtils \
    leafem \
    lwgeom \
    stars \
    leafpop \
    plainview \
    sf \
    mapview

# Installs packages needed for still treemap, interactive plots, and hex plots
RUN apt-get update -qq && apt-get -y --no-install-recommends install \    
    && install2.r --error \
    --deps TRUE \
    R.utils \
    treemap \
    d3r \
    hexbin \ 
    VennDiagram \
    rprojroot

# Use maftools for reading MAF files
RUN R -e "BiocManager::install(c('maftools'), update = FALSE)"

# This is needed to create the interactive pie chart
RUN R -e "devtools::install_github('timelyportfolio/sunburstR', ref = 'd40d7ed71ee87ca4fbb9cb8b7cf1e198a23605a9', dependencies = TRUE)"

# This is needed to create the interactive treemap
RUN R -e "devtools::install_github('timelyportfolio/d3treeR', ref = '0eaba7f1c6438e977f8a5c082f1474408ac1fd80', dependencies = TRUE)"

# Need this package to make plots colorblind friendly
RUN R -e "devtools::install_github('clauswilke/colorblindr', ref = '1ac3d4d62dad047b68bb66c06cee927a4517d678', dependencies = TRUE)"

#Need this to read command args
RUN R -e "install.packages('optparse')" 

#### Please install your dependencies here
#### Add a comment to indicate what analysis it is required for
