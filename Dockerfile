FROM rocker/tidyverse:3.6.0
MAINTAINER ccdl@alexslemonade.org
WORKDIR /rocker-build/

RUN apt-get update && apt-get install -y --no-install-recommends apt-utils

RUN apt-get install dialog apt-utils -y

# Required for installing mapview for interactive sample distribution plots
RUN apt-get update -qq && apt-get -y --no-install-recommends install \
    libgdal-dev \
    libudunits2-dev \
    # This is needed for coloblindr to install
    libmagick++-dev 

# Required forinteractive sample distribution plots
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
    # This is needed to create HTML outputs of the interactive plots
    mapview

RUN apt-get update -qq && apt-get -y --no-install-recommends install \    
    && install2.r --error \
    --deps TRUE \
    R.utils \
    # This is needed to create a still treemap
    treemap \
    # This is needed to convert a data.frame into a d3.js hierarchy object
    d3r \
    # Need hexbin for making a hex plot 
    hexbin \ 
    VennDiagram 

# Use maftools for reading MAF files
RUN R -e "BiocManager::install(c('maftools'), update = FALSE)"

# This is needed to create the interactive pie chart
RUN R -e "devtools::install_github('timelyportfolio/sunburstR', ref = 'd40d7ed71ee87ca4fbb9cb8b7cf1e198a23605a9', dependencies = TRUE)"

# This is needed to create the interactive treemap
RUN R -e "devtools::install_github('timelyportfolio/d3treeR', ref = '0eaba7f1c6438e977f8a5c082f1474408ac1fd80', dependencies = TRUE)"

# Need this package to make plots colorblind friendly
RUN R -e "devtools::install_github('clauswilke/colorblindr', ref = '1ac3d4d62dad047b68bb66c06cee927a4517d678', dependencies = TRUE)"

#### Please install your dependencies here
#### Add a comment to indicate what analysis it is required for
