FROM rocker/tidyverse:3.6.0
MAINTAINER ccdl@alexslemonade.org
WORKDIR /rocker-build/

RUN apt-get update && apt-get install -y --no-install-recommends apt-utils
RUN apt-get install dialog apt-utils -y

#### Please install your dependencies here

#### Add a comment to indicate what analysis it is required for

RUN apt-get update -qq && apt-get -y --no-install-recommends install \
    && install2.r --error \
    treemap \  # This is needed to create a still treemap

RUN apt-get update -qq && apt-get -y --no-install-recommends install \
    && install2.r --error \
    --deps TRUE \
    d3r \      # This is needed to convert a data.frame into a d3.js hierarchy object
    mapview \  # This is needed to create HTML outputs of the interactive plots
    gridExta \ # This is needed to arrange multiple plots into a grid 
    R.utils

# This is needed to create the interactive pie chart
RUN R -e "devtools::install_github('timelyportfolio/sunburstR', ref = 'd40d7ed71ee87ca4fbb9cb8b7cf1e198a23605a9', dependencies = TRUE)"

# This is needed to create the interactive treemap
RUN R -e "devtools::install_github('timelyportfolio/d3treeR', ref = '0eaba7f1c6438e977f8a5c082f1474408ac1fd80', dependencies = TRUE)"
