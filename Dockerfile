FROM rocker/tidyverse:3.6.2
MAINTAINER ccdl@alexslemonade.org
WORKDIR /rocker-build/

COPY scripts/install_bioc.r .

### Install apt-getable packages to start
#########################################
RUN apt-get update && apt-get install -y --no-install-recommends apt-utils dialog

# Add curl, bzip2 and some dev libs
RUN apt-get update -qq && apt-get -y --no-install-recommends install \
    curl \
    bzip2 \
    zlib1g \
    libbz2-dev \
    liblzma-dev \
    libreadline-dev

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
RUN pip3 install "setuptools==46.3.0" "six==1.14.0" "wheel==0.34.2"

# Install java
RUN apt-get -y --no-install-recommends install \
   default-jdk

# Standalone tools and libraries
################################

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
    mv bin/* /usr/local/bin

# HTSlib
RUN wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 && \
    tar -jxvf htslib-1.9.tar.bz2 && rm -f htslib-1.9.tar.bz2 && \
    cd htslib-1.9 && \
    ./configure && \
    make && \
    make install


#### R packages
###############

# Commonly used R packages
RUN install2.r --error --deps TRUE \
    class \
    cluster \
    data.table \
    DT \
    flextable \
    foreign \
    GGally \
    lattice \
    MASS \
    Matrix \
    mgcv \
    nlme \
    nnet \
    optparse \
    R.utils \
    RColorBrewer \
    rpart \
    rprojroot \
    survival \
    viridis


# Required for interactive sample distribution plots
# map view is needed to create HTML outputs of the interactive plots
RUN install2.r --error --deps TRUE \
    gdalUtils \
    leafem \
    leafpop \
    lwgeom \
    mapview \
    plainview \
    sf \
    stars

# Installs packages needed for plottings
# treemap, interactive plots, and hex plots
# Rtsne and umap are required for dimension reduction analyses
RUN install2.r --error --deps TRUE \
    corrplot \
    d3r \
    ggfortify \
    ggpubr \
    ggrepel \
    ggsci \
    ggsignif \
    gridGraphics \
    hexbin \
    pheatmap \
    Rtsne \
    spatial \
    treemap \
    umap  \
    UpSetR \
    VennDiagram

# Install rjava
RUN install2.r --error --deps TRUE \
    rJava

# Need for survminer for doing survival analysis
RUN install2.r --error --deps TRUE \
    cmprsk \
    survMisc \
    survminer

# maftools for proof of concept in create-subset-files
RUN ./install_bioc.r \
    maftools

# ComplexHeatmap & circlize
RUN ./install_bioc.r \
    ComplexHeatmap \
    circlize

# This is needed for the CNV frequency and proportion aberration plots
RUN ./install_bioc.r \
    GenVisR

# These packages are for the genomic region analysis for snv-callers
RUN ./install_bioc.r \
    annotatr \
    TxDb.Hsapiens.UCSC.hg38.knownGene \
    org.Hs.eg.db \
    BSgenome.Hsapiens.UCSC.hg19 \
    BSgenome.Hsapiens.UCSC.hg38

# Packages for expression normalization and batch correction
RUN ./install_bioc.r \
    preprocessCore \
    sva


## This is deprecated
#  # These packages are for single-sample GSEA analysis
#  RUN ./install_bioc.r 'GSEABase', 'GSVA'

# Required for sex prediction from RNA-seq data
RUN install2.r --error --deps TRUE \
    glmnet \
    glmnetUtils \
    caret \
    e1071


# bedr package
RUN install2.r --error --deps TRUE \
    bedr
# Check to make sure the binaries are available by loading the bedr library
RUN Rscript -e "library(bedr)"

# Also install for mutation signature analysis
# qdapRegex is for the fusion analysis
RUN install2.r --error --deps TRUE \
    deconstructSigs \
    qdapRegex

# packages required for collapsing RNA-seq data by removing duplicated gene symbols
RUN ./install_bioc.r \
    rtracklayer

# TCGAbiolinks for TMB compare analysis
RUN ./install_bioc.r \
    TCGAbiolinks

# Install for mutation signature analysis
RUN ./install_bioc.r \
    ggbio

# CRAN package msigdbr needed for gene-set-enrichment-analysis
RUN install2.r --error --deps TRUE \
    msigdbr
# Bioconductor package GSVA needed for gene-set-enrichment-analysis
RUN ./install_bioc.r \
    GSVA


# package required for immune deconvolution
RUN R -e "remotes::install_github('icbi-lab/immunedeconv', ref = '493bcaa9e1f73554ac2d25aff6e6a7925b0ea7a6', dependencies = TRUE)"

RUN R -e "remotes::install_github('const-ae/ggupset', ref = '7a33263cc5fafdd72a5bfcbebe5185fafe050c73', dependencies = TRUE)"

# This is needed to create the interactive pie chart
RUN R -e "remotes::install_github('timelyportfolio/sunburstR', ref = 'd40d7ed71ee87ca4fbb9cb8b7cf1e198a23605a9', dependencies = TRUE)"

# This is needed to create the interactive treemap
RUN R -e "remotes::install_github('timelyportfolio/d3treeR', ref = '0eaba7f1c6438e977f8a5c082f1474408ac1fd80', dependencies = TRUE)"

# Need this package to make plots colorblind friendly
RUN R -e "remotes::install_github('clauswilke/colorblindr', ref = '1ac3d4d62dad047b68bb66c06cee927a4517d678', dependencies = TRUE)"

# remote package EXTEND needed for telomerase-activity-prediciton analysis
RUN R -e "remotes::install_github('NNoureen/EXTEND', ref = '467c2724e1324ef05ad9260c3079e5b0b0366420', dependencies = TRUE)"

# package required for shatterseek
RUN R -e "withr::with_envvar(c(R_REMOTES_NO_ERRORS_FROM_WARNINGS='true'), remotes::install_github('parklab/ShatterSeek', ref = '83ab3effaf9589cc391ecc2ac45a6eaf578b5046', dependencies = TRUE))"

# Packages required for rna-seq-composition
RUN install2.r --error --deps TRUE \
    EnvStats \
    janitor

# Patchwork for plot compositions
RUN R -e "remotes::install_github('thomasp85/patchwork', ref = 'c67c6603ba59dd46899f17197f9858bc5672e9f4')"

# This is required for creating a treemap of the broad histology and integrated diagnoses
RUN R -e "remotes::install_github('wilkox/treemapify', ref = 'e70adf727f4d13223de8146458db9bef97f872cb', dependencies = TRUE)"

# Install python libraries
##########################

# Install python3 data science tools and all dependencies (versioned)
RUN pip3 install \
    "MarkupSafe==1.1.1" \
    "Send2Trash==1.5.0" \
    "appdirs==1.4.4" \
    "attrs==19.3.0" \
    "backcall==0.1.0" \
    "bleach==3.1.5" \
    "certifi==2020.4.5.1" \
    "chardet==3.0.4" \
    "configargparse==1.2.3" \
    "cycler==0.10.0" \
    "datrie==0.8.2" \
    "decorator==4.4.2" \
    "defusedxml==0.6.0" \
    "docutils==0.16" \
    "entrypoints==0.3" \
    "gitdb==4.0.5" \
    "gitpython==3.1.2" \
    "idna==2.9" \
    "importlib-metadata==1.6.0" \
    "ipykernel==5.2.1" \
    "ipython==7.14.0" \
    "ipython-genutils==0.2.0" \
    "ipywidgets==7.5.1" \
    "jedi==0.17.0" \
    "jinja2==2.11.2" \
    "joblib==0.15.1" \
    "jsonschema==3.2.0" \
    "jupyter==1.0.0" \
    "jupyter-client==6.1.3" \
    "jupyter-console==6.1.0" \
    "jupyter-core==4.6.3" \
    "kiwisolver==1.1.0" \
    "matplotlib==3.2.1" \
    "mistune==0.8.4" \
    "mizani==0.6.0" \
    "nbconvert==5.6.1" \
    "nbformat==5.0.6" \
    "notebook==6.0.3" \
    "numpy==1.18.4" \
    "packaging==20.3" \
    "palettable==3.3.0" \
    "pandas==0.25.3" \
    "pandocfilters==1.4.2" \
    "parso==0.7.0" \
    "patsy==0.5.1" \
    "pexpect==4.8.0" \
    "pickleshare==0.7.5" \
    "plotnine==0.3.0" \
    "prometheus-client=0.7.1" \
    "prompt-toolkit=3.0.5" \
    "psutil==5.7.0" \
    "ptyprocess==0.6.0" \
    "pyarrow==0.16.0" \
    "pygments==2.6.1" \
    "pyparsing==2.4.5" \
    "pyreadr==0.2.6" \
    "pyrsistent==0.16.0" \
    "python-dateutil=2.8.1" \
    "pytz==2019.3" \
    "pyyaml==5.3.1" \
    "pyzmq==19.0.1" \
    "qtconsole==4.7.4" \
    "qtpy==1.9.0" \
    "ratelimiter==1.2.0.post0" \
    "requests==2.23.0" \
    "scikit-learn==0.22.2" \
    "scipy==1.4.1" \
    "seaborn==0.9.0" \
    "smmap==3.0.4" \
    "snakemake==5.8.1" \
    "statsmodels==0.10.2" \
    "terminado==0.8.3" \
    "testpath==0.4.4" \
    "tornado==6.0.4" \
    "traitlets==4.3.3" \
    "tzlocal==2.0.0" \
    "urllib3==1.25.9" \
    "wcwidth==0.1.9" \
    "webencodings==0.5.1" \
    "widgetsnbextension==3.5.1" \
    "wrapt==1.12.1" \
    "zipp==3.1.0"

# Install Rpy2
RUN pip3 install "rpy2==3.2.7"

# Install CrossMap for liftover
RUN pip3 install "cython==0.29.15" && \
    pip3 install "bx-python==0.8.8" && \
    pip3 install "pybigwig==0.3.17" && \
    pip3 install "pysam==0.15.4" && \
    pip3 install "CrossMap==0.3.9"


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
    wget -q https://www.mathworks.com/supportfiles/downloads/R2014a/deployment_files/R2014a/installers/glnxa64/MCR_R2014a_glnxa64_installer.zip && \
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
RUN chown -R rstudio:rstudio /home/rstudio/gistic_install && \
    chmod 755 /home/rstudio/gistic_install
WORKDIR /rocker-build/



#### Please install your dependencies immediately above this comment.
#### Add a comment to indicate what analysis it is required for


WORKDIR /rocker-build/
