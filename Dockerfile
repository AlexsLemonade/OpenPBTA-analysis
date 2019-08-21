FROM rocker/tidyverse:3.6.0
MAINTAINER ccdl@alexslemonade.org
WORKDIR /rocker-build/

RUN apt-get update && apt-get install -y --no-install-recommends apt-utils
RUN apt-get install dialog apt-utils -y

#### Please install your dependencies here
#### Add a comment to indicate what analysis it is required for

