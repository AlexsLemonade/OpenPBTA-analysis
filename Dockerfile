FROM rocker/tidyverse:3.6.2
MAINTAINER ccdl@alexslemonade.org
WORKDIR /rocker-build/

### Install apt-getable packages to start
#########################################
RUN apt-get update && apt-get install -y --no-install-recommends apt-utils dialog

# Add curl & bzip2
RUN apt-get update -qq && apt-get -y --no-install-recommends install \
    curl \
    bzip2

# Required for installing htslib
RUN apt-get -y --no-install-recommends install \
    zlib1g \
    libbz2-dev \
    liblzma-dev

# libmagick++-dev is needed for coloblindr to install
RUN apt-get -y --no-install-recommends install \
    libgdal-dev \
    libudunits2-dev \
    libmagick++-dev

# Required for installing pdftools, which is a dependency of gridGraphics
RUN apt-get -y --no-install-recommends install \
    libpoppler-cpp-dev

# Install pip3 and instalation tools
RUN apt-get -y --no-install-recommends install \
    python3-pip  python3-dev
RUN pip3 install "setuptools==46.2.0" "wheel==0.34.2" "six==1.14.0"

# Install java and rJava for some of the snv plotting comparison packages
RUN apt-get -y --no-install-recommends install \
   default-jdk \
   r-cran-rjava # \
#    && apt-get clean \
#    && rm -rf /var/lib/apt/lists/

#### R packages and python below
################################

# Commonly used R packages
RUN install2.r --error \
    --deps TRUE \
    rprojroot \
    optparse \
    data.table \
    RColorBrewer \
    viridis \
    R.utils \
    lattice \
    rpart \
    class \
    MASS \
    GGally \
    Matrix \
    survival \
    nlme \
    cluster \
    foreign \
    nnet \
    mgcv \
    flextable


# Required for interactive sample distribution plots
# map view is needed to create HTML outputs of the interactive plots
RUN install2.r --error \
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
RUN install2.r --error \
    --deps TRUE \
    treemap \
    hexbin \
    VennDiagram \
    Rtsne \
    umap  \
    d3r


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
RUN R -e "remotes::install_github('timelyportfolio/sunburstR', ref = 'd40d7ed71ee87ca4fbb9cb8b7cf1e198a23605a9', dependencies = TRUE)"

# This is needed to create the interactive treemap
RUN R -e "remotes::install_github('timelyportfolio/d3treeR', ref = '0eaba7f1c6438e977f8a5c082f1474408ac1fd80', dependencies = TRUE)"

# Need this package to make plots colorblind friendly
RUN R -e "remotes::install_github('clauswilke/colorblindr', ref = '1ac3d4d62dad047b68bb66c06cee927a4517d678', dependencies = TRUE)"

# Required for sex prediction from RNA-seq data
RUN install2.r --error \
    --deps TRUE \
    glmnet \
    glmnetUtils \
    caret \
    e1071


# Install for SNV comparison plots
RUN install2.r --error \
    --deps TRUE \
    UpSetR

RUN R -e "remotes::install_github('const-ae/ggupset', ref = '7a33263cc5fafdd72a5bfcbebe5185fafe050c73', dependencies = TRUE)"

# Required for mapping segments to genes
# Add bedtools
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.28.0/bedtools-2.28.0.tar.gz && \
    tar -zxvf bedtools-2.28.0.tar.gz && \
    cd bedtools2 && \
    make && \
    mv bin/* /usr/local/bin

# Add bedops per the BEDOPS documentation
RUN wget https://github.com/bedops/bedops/releases/download/v2.4.37/bedops_linux_x86_64-v2.4.37.tar.bz2 && \
    tar -jxvf bedops_linux_x86_64-v2.4.37.tar.bz2 && \
    rm -f bedops_linux_x86_64-v2.4.37.tar.bz2 && \
    cp bin/* /usr/local/bin

# HTSlib
RUN wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2
RUN tar -jxvf htslib-1.9.tar.bz2 && rm -f htslib-1.9.tar.bz2
RUN cd htslib-1.9 && \
    ./configure && \
    make && \
    make install && \
    mv bin/* /usr/local/bin

# bedr package
RUN install2.r --error \
    bedr

# Check to make sure the binaries are available by loading the bedr library
RUN Rscript -e "library(bedr)"

# Install for mutation signature analysis
RUN R -e "BiocManager::install(c('BSgenome.Hsapiens.UCSC.hg19', 'BSgenome.Hsapiens.UCSC.hg38'))"

# Also install for mutation signature analysis
# qdapRegex is for the fusion analysis
RUN install2.r --error \
    --deps TRUE \
    deconstructSigs \
    qdapRegex

# packages required for collapsing RNA-seq data by removing duplicated gene symbols
RUN R -e "install.packages('DT', dependencies = TRUE)"
RUN R -e "BiocManager::install(c('rtracklayer'), update = FALSE)"


# TCGAbiolinks for TMB compare analysis
RUN R -e "BiocManager::install(c('TCGAbiolinks'), update = FALSE)"

# Install python3 data science tools
RUN pip3 install "numpy==1.18.4" && \
    pip3 install "cycler==0.10.0" "kiwisolver==1.1.0" "pyparsing==2.4.5" "python-dateutil==2.8.1" "pytz==2019.3" && \
    pip3 install "matplotlib==3.2.1" && \
    pip3 install "scipy==1.4.1" && \
    pip3 install "pandas==1.0.3" && \
    pip3 install "scikit-learn==0.23.0" &&\
    pip3 install "jupyter==1.0.0" && \
    pip3 install "ipykernel==5.2.1" && \
    pip3 install "widgetsnbextension==3.5.1" && \
    pip3 install "snakemake==5.17.0" && \
    pip3 install "statsmodels==0.11.1" && \
    pip3 install "plotnine==0.6.0" && \
    pip3 install "rpy2==3.3.2" && \
    pip3 install "seaborn==0.10.1" && \
    pip3 install "tzlocal==2.1"


# pyreadr for comparative-RNASeq-analysis
RUN pip3 install "pyreadr==0.2.1" && \
    pip3 install "pyarrow==0.16.0"

# Install CrossMap for liftover
RUN pip3 install "cython==0.29.15" && \
    pip3 install "bx-python==0.8.8" && \
    pip3 install "pybigwig==0.3.17" && \
    pip3 install "pysam==0.15.4" && \
    pip3 install "CrossMap==0.3.9"

# Need for survminer for doing survival analysis
RUN install2.r --error \
    --deps TRUE \
    cmprsk \
    survMisc \
    survminer

# ggfortify for plotting
RUN install2.r --error \
    --deps TRUE \
    spatial \
    ggfortify

# package required for immune deconvolution
RUN R -e "remotes::install_github('icbi-lab/immunedeconv', ref = '493bcaa9e1f73554ac2d25aff6e6a7925b0ea7a6', dependencies = TRUE)"
RUN install2.r --error \
    --deps TRUE \
    corrplot

# Install for mutation signature analysis
RUN R -e "BiocManager::install('ggbio')"

# CRAN package msigdbr needed for gene-set-enrichment-analysis
RUN install2.r --error \
    --deps TRUE \
    msigdbr

# Bioconductor package GSVA needed for gene-set-enrichment-analysis
RUN R -e "BiocManager::install(c('GSVA'), update = FALSE)"

# remote package EXTEND needed for telomerase-activity-prediciton analysis
RUN R -e "remotes::install_github('NNoureen/EXTEND', ref = '467c2724e1324ef05ad9260c3079e5b0b0366420', dependencies = TRUE)"



# CRAN package gridGraphics needed for telomerase-activity-prediction
RUN install2.r --error \
    --deps TRUE \
    gridGraphics

# package required for shatterseek
RUN R -e "withr::with_envvar(c(R_REMOTES_NO_ERRORS_FROM_WARNINGS='true'), remotes::install_github('parklab/ShatterSeek', ref = '83ab3effaf9589cc391ecc2ac45a6eaf578b5046', dependencies = TRUE))"


# ComplexHeatmap and circlize were apparently not explicitly installed anywhere, but is used throughout
RUN R -e "BiocManager::install(c('ComplexHeatmap', 'circlize'), update = FALSE)"

# MATLAB Compiler Runtime is required for GISTIC, MutSigCV
# Install steps are adapted from usuresearch/matlab-runtime
# https://hub.docker.com/r/usuresearch/matlab-runtime/dockerfile

ENV DEBIAN_FRONTEND noninteractive
RUN apt-get -q update && \
    apt-get install -q -y --no-install-recommends \
    xorg

# This is the version of MCR required to run the precompiled version of GISTIC
RUN mkdir /mcr-install-v83 && \
    mkdir /opt/mcr && \
    cd /mcr-install-v83 && \
    wget https://www.mathworks.com/supportfiles/downloads/R2014a/deployment_files/R2014a/installers/glnxa64/MCR_R2014a_glnxa64_installer.zip && \
    unzip -q MCR_R2014a_glnxa64_installer.zip && \
    rm -f MCR_R2014a_glnxa64_installer.zip && \
    ./install -destinationFolder /opt/mcr -agreeToLicense yes -mode silent && \
    cd / && \
    rm -rf mcr-install-v83

WORKDIR /home/rstudio/

# GISTIC installation
RUN mkdir -p gistic_install && \
    cd gistic_install && \
    wget -q ftp://ftp.broadinstitute.org/pub/GISTIC2.0/GISTIC_2_0_23.tar.gz && \
    tar zxf GISTIC_2_0_23.tar.gz && \
    rm -f GISTIC_2_0_23.tar.gz && \
    rm -rf MCR_Installer

RUN chown -R rstudio:rstudio /home/rstudio/gistic_install
RUN chmod 755 /home/rstudio/gistic_install


# Packages required for rna-seq-composition
RUN install2.r --error \
    --deps TRUE \
    EnvStats \
    janitor

# Patchwork for plot compositions
RUN R -e "remotes::install_github('thomasp85/patchwork', ref = 'c67c6603ba59dd46899f17197f9858bc5672e9f4')"

# This is required for creating a treemap of the broad histology and integrated diagnoses
RUN R -e "remotes::install_github('wilkox/treemapify', ref = 'e70adf727f4d13223de8146458db9bef97f872cb', dependencies = TRUE)"

#### Please install your dependencies immediately above this comment.
#### Add a comment to indicate what analysis it is required for


WORKDIR /rocker-build/
