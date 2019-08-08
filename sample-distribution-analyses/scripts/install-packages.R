####install required packages
if(!require(tidyr)){
  install.packages("tidyr", repos='http://cran.us.r-project.org')
  library(tidyr)
}
if(!require(dplyr)){
  install.packages("dplyr", repos='http://cran.us.r-project.org')
  library(dplyr)
}
if(!require(readr)){
  install.packages("readr", repos='http://cran.us.r-project.org')
  library(readr)
}
if(!require(ggplot2)){
  install.packages("ggplot2", repos='http://cran.us.r-project.org')
  library(ggplot2)
}
if(!require(stringr)){
  install.packages("stringr", repos='http://cran.us.r-project.org')
  library(stringr)
}
if(!require(knitr)){
  install.packages("knitr", repos='http://cran.us.r-project.org')
  library(knitr)
}
if(!require(devtools)){
  install.packages("devtools", repos='http://cran.us.r-project.org')
  library(devtools)
}
if(!require(sunburstR)){
  devtools::install_github("timelyportfolio/sunburstR")
  library(sunburstR)
}
if(!require(formattable)){
  install.packages("formattable", repos='http://cran.us.r-project.org')
  library(formattable)
}
if(!require(kableExtra)){
  install.packages("kableExtra", repos='http://cran.us.r-project.org')
  library(kableExtra)
}
if(!require(d3r)){
  install.packages("d3r", repos='http://cran.us.r-project.org')
  library(d3r)
}
if(!require(treemap)){
  install.packages("treemap", repos='http://cran.us.r-project.org')
  library(treemap)
}
if(!require(htmltools)){
  install.packages("htmltools", repos='http://cran.us.r-project.org')
  library(htmltools)
}
if(!require(colorblindr)){
  install.packages("colorspace", repos = "http://R-Forge.R-project.org")
  library(colorspace)
}
if(!require(colorblindr)){
  devtools::install_github("clauswilke/colorblindr")
  library(colorblindr)
}
if(!require(mapview)){
  install.packages("mapview", repos = 'http://cran.us.r-project.org')
  library(mapview)
}

