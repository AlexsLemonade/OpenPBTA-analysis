FROM rocker/tidyverse:3.6.0
MAINTAINER ccdl@alexslemonade.org

RUN apt-get update && apt-get install -y --no-install-recommends apt-utils
RUN apt-get install dialog apt-utils -y

#### Please install your dependencies here
#### Add a comment to indicate what analysis it is required for

RUN groupadd user && useradd --create-home --home-dir /home/user -g user user
WORKDIR /home/user

