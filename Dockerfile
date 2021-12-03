FROM rocker/tidyverse:3.6.0
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

# Install pip3 and low-level python installation reqs
RUN apt-get -y --no-install-recommends install \
    python3-pip  python3-dev
RUN pip3 install \
  "Cython==0.29.15" \
  "setuptools==46.3.0" \
  "six==1.14.0" \
  "wheel==0.34.2"

# Install java
RUN apt-get -y --no-install-recommends install \
   default-jdk


# Standalone tools and libraries
################################

# Required for mapping segments to genes
# Add bedtools
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.28.0/bedtools-2.28.0.tar.gz && \
    tar -zxvf bedtools-2.28.0.tar.gz && rm -f bedtools-2.28.0.tar.gz && \
    cd bedtools2 && \
    make && \
    mv bin/* /usr/local/bin && \
    cd .. && rm -rf bedtools2

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
    make install && \
    cd .. && rm -rf htslib-1.9


#### R packages
###############

# Commonly used R packages
RUN ./install_bioc.r \
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
    viridis \
    openxlsx


# Required for interactive sample distribution plots
# map view is needed to create HTML outputs of the interactive plots
RUN ./install_bioc.r \
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
RUN ./install_bioc.r \
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
RUN ./install_bioc.r \
    rJava

# Need for survminer for doing survival analysis
RUN ./install_bioc.r \
    cmprsk \
    survMisc \
    survminer

# maftools for proof of concept in create-subset-files
RUN R -e "remotes::install_github('PoisonAlien/maftools', ref = '9719868262f946e0b8eb2e7ec2510ee18c6cafa3')"

# ComplexHeatmap
RUN ./install_bioc.r \
    ComplexHeatmap

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
RUN ./install_bioc.r \
    glmnet \
    glmnetUtils \
    caret \
    e1071


# bedr package & check to make sure binaries are available by loading
RUN ./install_bioc.r \
    bedr \
    && Rscript -e "library(bedr)"

# Also install for mutation signature analysis
# qdapRegex is for the fusion analysis
RUN ./install_bioc.r \
    deconstructSigs \
    qdapRegex

# packages required for collapsing RNA-seq data by removing duplicated gene symbols
RUN ./install_bioc.r \
    rtracklayer

# TCGAbiolinks for TMB compare analysis
RUN R -e "remotes::install_github('RDocTaskForce/parsetools', ref = '1e682a9f4c5c7192d22e8985ce7723c09e98d62b', dependencies = TRUE)" \
    && R -e "remotes::install_github('RDocTaskForce/testextra', ref = '4e5dfac8853c08d5c2a8790a0a1f8165f293b4be', dependencies = TRUE)" \
    && R -e "remotes::install_github('halpo/purrrogress', ref = '54f2130477f161896e7b271ed3ea828c7e4ccb1c', dependencies = TRUE)" \
    && ./install_bioc.r TCGAbiolinks

# Install for mutation signature analysis
RUN ./install_bioc.r \
    ggbio

# CRAN package msigdbr and GSVA for gene-set-enrichment-analysis
RUN ./install_bioc.r \
    msigdbr \
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
RUN ./install_bioc.r \
    EnvStats \
    janitor

# Patchwork for plot compositions
RUN R -e "remotes::install_github('thomasp85/patchwork', ref = 'c67c6603ba59dd46899f17197f9858bc5672e9f4')"

# This is required for creating a treemap of the broad histology and integrated diagnoses
RUN R -e "remotes::install_github('wilkox/treemapify', ref = 'e70adf727f4d13223de8146458db9bef97f872cb', dependencies = TRUE)"

# Need this specific version of circlize so it has hg38
RUN R -e "remotes::install_github('jokergoo/circlize', ref = 'b7d86409d7f893e881980b705ba1dbc758df847d', dependencies = TRUE)"

# Install python packages
##########################

# Install python3 tools and ALL dependencies
RUN pip3 install \
    "appdirs==1.4.4" \
    "attrs==20.3.0" \
    "backcall==0.2.0" \
    "bleach==3.3.0" \
    "bx-python==0.8.8" \
    "certifi==2020.12.5" \
    "chardet==4.0.0" \
    "ConfigArgParse==1.4" \
    "CrossMap==0.3.9" \
    "cycler==0.10.0" \
    "datrie==0.8.2" \
    "decorator==4.4.2" \
    "defusedxml==0.7.1" \
    "docutils==0.16" \
    "entrypoints==0.3" \
    "gitdb==4.0.7" \
    "GitPython==3.1.14" \
    "idna==2.10" \
    "importlib-metadata==2.1.1" \
    "ipykernel==4.8.1" \
    "ipython==7.9.0" \
    "ipython-genutils==0.2.0" \
    "jedi==0.17.2" \
    "Jinja2==2.11.3" \
    "jsonschema==3.2.0" \
    "jupyter-client==6.1.12" \
    "jupyter-core==4.6.3" \
    "kiwisolver==1.1.0" \
    "MarkupSafe==1.1.1" \
    "matplotlib==3.0.3" \
    "mistune==0.8.4" \
    "mizani==0.5.4" \
    "nbconvert==5.6.1" \
    "nbformat==5.1.2" \
    "notebook==6.0.0" \
    "numpy==1.17.3" \
    "packaging==20.9" \
    "palettable==3.3.0" \
    "pandas==0.25.3" \
    "pandocfilters==1.4.3" \
    "parso==0.7.1" \
    "patsy==0.5.1" \
    "pexpect==4.8.0" \
    "pickleshare==0.7.5" \
    "plotnine==0.3.0" \
    "prometheus-client==0.9.0" \
    "prompt-toolkit==2.0.10" \
    "psutil==5.8.0" \
    "ptyprocess==0.7.0" \
    "pyarrow==0.16.0" \
    "pybedtools==0.8.1" \
    "pyBigWig==0.3.17" \
    "Pygments==2.8.1" \
    "pyparsing==2.4.5" \
    "pyreadr==0.2.1" \
    "pyrsistent==0.17.3" \
    "pysam==0.15.4" \
    "python-dateutil==2.8.1" \
    "pytz==2019.3" \
    "PyYAML==5.3.1" \
    "pyzmq==20.0.0" \
    "ratelimiter==1.2.0.post0" \
    "requests==2.25.1" \
    "rpy2==2.9.3" \
    "scikit-learn==0.19.1" \
    "scipy==1.3.2" \
    "seaborn==0.8.1" \
    "Send2Trash==1.5.0" \
    "six==1.14.0" \
    "smmap==4.0.0" \
    "snakemake==5.8.1" \
    "statsmodels==0.10.2" \
    "terminado==0.8.3" \
    "testpath==0.4.4" \
    "tornado==6.1" \
    "traitlets==4.3.3" \
    "tzlocal==2.0.0" \
    "urllib3==1.26.4" \
    "wcwidth==0.2.5" \
    "webencodings==0.5.1" \
    "widgetsnbextension==2.0.0" \
    "wrapt==1.12.1" \
    "zipp==1.2.0" \
    && rm -rf /root/.cache/pip/wheels


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
    rm -rf MCR_Installer && \
    chown -R rstudio:rstudio /home/rstudio/gistic_install && \
    chmod 755 /home/rstudio/gistic_install
WORKDIR /rocker-build/

# Install multipanelfigure, required for transcriptomic overview figure
# gplots for gistic comparison
RUN ./install_bioc.r \
    multipanelfigure \
    gplots

# Molecular subtyping MB
RUN R -e "remotes::install_github('d3b-center/medullo-classifier-package', ref = 'e3d12f64e2e4e00f5ea884f3353eb8c4b612abe8', dependencies = TRUE, upgrade = FALSE)" \
    && ./install_bioc.r MM2S
# More recent version of sva required for molecular subtyping MB
RUN R -e "remotes::install_github('jtleek/sva-devel@123be9b2b9fd7c7cd495fab7d7d901767964ce9e', dependencies = FALSE, upgrade = FALSE)"

# Packages required for de novo mutational signatures
RUN install2.r --error --deps TRUE \
  lsa

# To install sigfit, we need a more recent version of rstantools than we can obtain via the MRAN snapshot route
# We're using the ref for the most recent release on GitHub (2.0.0)
RUN R -e "remotes::install_github('stan-dev/rstantools', ref = 'd43bf9fb6120d40a60e708853e4b80cdb4689d19', dependencies = TRUE)"

# Build arguments are according to the sigfit instructions
RUN R -e "remotes::install_github('kgori/sigfit', ref = '209776ee1d2193ad4b682b2e2472f848bd7c67a6', build_vignettes = TRUE, build_opts = c('--no-resave-data', '--no-manual'), dependencies = TRUE)"

# Package for kinase domain retention for fusions
RUN ./install_bioc.r \
     EnsDb.Hsapiens.v86 \
     ensembldb

RUN R -e "remotes::install_github('d3b-center/annoFuse',ref = 'c6a2111b5949ca2aae3853f7f34de3d0db4ffa33', dependencies = TRUE)"


# gmp, dependency for signature.tools.lib
RUN apt-get -y --no-install-recommends install \
    libgmp-dev

# CNS signatures can be obtained from signature.tools.lib
RUN R -e "remotes::install_github('Nik-Zainal-Group/signature.tools.lib', ref = '73e899c9090a215a76a307480bda76c241a4a489')"


#### Please install your dependencies immediately above this comment.
#### Add a comment to indicate what analysis it is required for


WORKDIR /rocker-build/
