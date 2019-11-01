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
# Rtsne and umap are required for dimension reduction analyses 
# optparse is needed for passing arguments from the command line to R script
RUN apt-get update -qq && apt-get -y --no-install-recommends install \    
    && install2.r --error \
    --deps TRUE \
    R.utils \
    treemap \
    d3r \
    hexbin \
    VennDiagram \
    Rtsne \
    umap \
    rprojroot \
    optparse \
    pheatmap \
    RColorBrewer \
    viridis \
    data.table 

# maftools for proof of concept in create-subset-files
RUN R -e "BiocManager::install(c('maftools'), update = FALSE)"

# This is needed for the CNV frequency and proportion aberration plots
RUN R -e "BiocManager::install(c('GenVisR'), update = FALSE)"

# These packages are for the genomic region analysis for snv-callers
RUN R -e "BiocManager::install(c('annotatr', 'TxDb.Hsapiens.UCSC.hg38.knownGene', 'org.Hs.eg.db'), update = FALSE)"

# Packages for expression normalization and batch correction
RUN R -e "BiocManager::install(c('preprocessCore', 'sva'), update = FALSE)"

# These packages are for single-sample GSEA analysis
RUN R -e "BiocManager::install(c('GSEABase', 'GSVA'), update = FALSE)"

# This is needed to create the interactive pie chart
RUN R -e "devtools::install_github('timelyportfolio/sunburstR', ref = 'd40d7ed71ee87ca4fbb9cb8b7cf1e198a23605a9', dependencies = TRUE)"

# This is needed to create the interactive treemap
RUN R -e "devtools::install_github('timelyportfolio/d3treeR', ref = '0eaba7f1c6438e977f8a5c082f1474408ac1fd80', dependencies = TRUE)"

# Need this package to make plots colorblind friendly
RUN R -e "devtools::install_github('clauswilke/colorblindr', ref = '1ac3d4d62dad047b68bb66c06cee927a4517d678', dependencies = TRUE)"


# Required for sex prediction from RNA-seq data
RUN apt-get update -qq && apt-get -y --no-install-recommends install \    
    && install2.r --error \
    --deps TRUE \
    glmnet \
    glmnetUtils


# Install java and rJava for some of the snv plotting comparison packages
RUN apt-get -y update && apt-get install -y \
   default-jdk \
   r-cran-rjava \
   && apt-get clean \
   && rm -rf /var/lib/apt/lists/

# Install for SNV comparison plots
RUN apt-get update -qq && apt-get -y --no-install-recommends install \    
    && install2.r --error \
    --deps TRUE \
    UpSetR

RUN R -e "devtools::install_github('const-ae/ggupset', ref = '7a33263cc5fafdd72a5bfcbebe5185fafe050c73', dependencies = TRUE)"

# GGally and its required packages
RUN apt-get update -qq && apt-get -y --no-install-recommends install \    
    && install2.r --error \
    lattice \
    rpart \
    class \
    MASS \
    GGally


#### Please install your dependencies here
#### Add a comment to indicate what analysis it is required for
