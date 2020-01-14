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


## This is deprecated 
#  # These packages are for single-sample GSEA analysis
#  RUN R -e "BiocManager::install(c('GSEABase', 'GSVA'), update = FALSE)"


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
    glmnetUtils \
    caret \
    e1071

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
    GGally \
    Matrix

# Help display tables in R Notebooks
RUN apt-get update -qq && apt-get -y --no-install-recommends install \
    && install2.r --error \
    flextable

# Required for mapping segments to genes
# Add bedtools
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.28.0/bedtools-2.28.0.tar.gz
RUN tar -zxvf bedtools-2.28.0.tar.gz
RUN cd bedtools2 && \
    make && \
    mv bin/* /usr/local/bin

# Required for installing htslib
RUN apt-get update -qq && apt-get -y --no-install-recommends install \
    zlib1g \
    libbz2-dev \
    liblzma-dev

# Add bedops per the BEDOPS documentation
RUN wget https://github.com/bedops/bedops/releases/download/v2.4.37/bedops_linux_x86_64-v2.4.37.tar.bz2
RUN tar -jxvf bedops_linux_x86_64-v2.4.37.tar.bz2
RUN cp bin/* /usr/local/bin

# HTSlib
RUN wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2
RUN tar -jxvf htslib-1.9.tar.bz2
RUN cd htslib-1.9 && \
    ./configure && \
    make && \
    make install
RUN mv bin/* /usr/local/bin

# bedr package
RUN apt-get update -qq && apt-get -y --no-install-recommends install \
    && install2.r --error \
    bedr

# Check to make sure the binaries are available by loading the bedr library
RUN Rscript -e "library(bedr)"

# Install for mutation signature analysis
RUN R -e "BiocManager::install(c('BSgenome.Hsapiens.UCSC.hg19', 'BSgenome.Hsapiens.UCSC.hg38'))"

# Also install for mutation signature analysis
# qdapRegex is for the fusion analysis
RUN apt-get update -qq && apt-get -y --no-install-recommends install \
    && install2.r --error \
    --deps TRUE \
    deconstructSigs \
    qdapRegex

# packages required for collapsing RNA-seq data by removing duplicated gene symbols
RUN R -e "install.packages('DT', dependencies = TRUE)"
RUN R -e "BiocManager::install(c('rtracklayer'), update = FALSE)"

# Needed to install TCGAbiolinks
RUN apt-get update -qq && apt-get -y --no-install-recommends install \
    && install2.r --error \
    --deps TRUE \
    survival \
    nlme \
    cluster \
    foreign \
    nnet \
    mgcv

# TCGAbiolinks for TMB compare analysis
RUN R -e "BiocManager::install(c('TCGAbiolinks'), update = FALSE)"

# Install python3 data science basics (pandas)
# using pip to get more current versions
RUN apt-get update -qq && apt-get -y --no-install-recommends install python3-pip  python3-dev 
RUN pip3 install "numpy==1.17.3" && \
   pip3 install "six==1.13.0" "setuptools==41.6.0" && \
   pip3 install "cycler==0.10.0" "kiwisolver==1.1.0" "pyparsing==2.4.5" "python-dateutil==2.8.1" "pytz==2019.3" && \
   pip3 install "matplotlib==3.0.3" && \
   pip3 install "scipy==1.3.2" && \
   pip3 install "pandas==0.25.3" && \
   pip3 install "snakemake==5.8.1"


# pip install for modules Ras, NF1, and TP53 Classifiers
RUN pip3 install "statsmodels==0.10.2" && \
   pip3 install "plotnine==0.3.0" && \
   pip3 install "scikit-learn==0.19.1" &&\
   pip3 install "rpy2==2.9.3" && \
   pip3 install "seaborn==0.8.1" && \
   pip3 install "jupyter==1.0.0" && \
   pip3 install "ipykernel==4.8.1" && \
   pip3 install "widgetsnbextension==2.0.0" && \
   pip3 install "tzlocal"


# Add curl
RUN apt-get update && apt-get install -y --no-install-recommends curl

# Need for survminer for doing survival analysis
RUN apt-get update -qq && apt-get -y --no-install-recommends install \
    && install2.r --error \
    --deps TRUE \
    survival \
    cmprsk \
    survMisc \
    survminer 

# pyreadr for comparative-RNASeq-analysis
RUN pip3 install "pyreadr==0.2.1"

# ggfortify for plotting
RUN apt-get update -qq && apt-get -y --no-install-recommends install \
    && install2.r --error \
    --deps TRUE \
    spatial \
    ggfortify

# package required for immune deconvolution
RUN R -e "install.packages('remotes', dependencies = TRUE)"
RUN R -e "remotes::install_github('icbi-lab/immunedeconv', ref = '493bcaa9e1f73554ac2d25aff6e6a7925b0ea7a6', dependencies = TRUE)"
RUN R -e "install.packages('corrplot', dependencies = TRUE)"

# Install for mutation signature analysis
RUN R -e "BiocManager::install('ggbio')"

#### Please install your dependencies here
#### Add a comment to indicate what analysis it is required for

# CRAN package msigdbr needed for gene-set-enrichment-analysis
RUN apt-get update -qq && apt-get -y --no-install-recommends install \
    && install2.r --error \
    --deps TRUE \
    msigdbr

# Bioconductor package GSVA needed for gene-set-enrichment-analysis
RUN R -e "BiocManager::install(c('GSVA'), update = FALSE)"
