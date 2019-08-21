FROM rocker/tidyverse:3.6.0
MAINTAINER ccdl@alexslemonade.org
WORKDIR /rocker-build/

RUN apt-get update && apt-get install -y --no-install-recommends apt-utils
RUN apt-get install dialog apt-utils -y

#### Please install your dependencies here

#### Add a comment to indicate what analysis it is required for

RUN apt-get update -qq && apt-get -y --no-install-recommends install \
    && install2.r --error \
    --deps TRUE \
    readr \    # This is needed to read in the dataset files
    d3r \      # This is needed to convert a data.frame into a d3.js hierarchy object
    treemap \  # This is needed to create a still treemap
    mapview \  # This is needed to create HTML outputs of the interactive plots
    gridExta \ # This is needed to arrange multiple plots into a grid 
    R.utils

# This is needed to create the interactive pie chart
RUN R -e "devtools::install_github('timelyportfolio/sunburstR', dependencies = TRUE)"

# This is needed to create the interactive treemap
RUN R -e "devtools::install_github('timelyportfolio/d3treeR', dependencies = TRUE)"
