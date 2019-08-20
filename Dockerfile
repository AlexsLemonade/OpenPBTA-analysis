FROM rocker/tidyverse:3.6.0
MAINTAINER ccdl@alexslemonade.org

RUN apt-get update && apt-get install -y --no-install-recommends apt-utils

RUN apt-get install dialog apt-utils -y

RUN groupadd user && useradd --create-home --home-dir /home/user -g user user
WORKDIR /home/user

USER user

ENTRYPOINT []

RUN apt-get update -qq && apt-get -y --no-install-recommends install \
    libmagick++-dev 

RUN apt-get update -qq && apt-get -y --no-install-recommends install \
    && install2.r --error \
    --deps TRUE \
    hexbin

RUN R -e "BiocManager::install(c('maftools'), update = FALSE)"

RUN R -e "devtools::install_github('clauswilke/colorblindr', ref = '1ac3d4d62dad047b68bb66c06cee927a4517d678', dependencies = TRUE)"
