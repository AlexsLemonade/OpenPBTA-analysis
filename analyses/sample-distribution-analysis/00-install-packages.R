# Install required packages
if (!("dplyr" %in% installed.packages())) {
    install.packages("dplyr", repos = "http://cran.us.r-project.org")
}
if (!("readr" %in% installed.packages())) {
    install.packages("readr", repos = "http://cran.us.r-project.org")
}
if (!("ggplot2" %in% installed.packages())) {
    install.packages("ggplot2", repos = "http://cran.us.r-project.org")
}
if (!("stringr" %in% installed.packages())) {
    install.packages("stringr", repos = "http://cran.us.r-project.org")
}
if (!("knitr" %in% installed.packages())) {
    install.packages("knitr", repos = "http://cran.us.r-project.org")
}
if (!("devtools" %in% installed.packages())) {
    install.packages("devtools", repos = "http://cran.us.r-project.org")
}
if (!("sunburstR" %in% installed.packages())) {
    devtools::install_github("timelyportfolio/sunburstR")
}
if (!("d3r" %in% installed.packages())) {
    install.packages("d3r", repos = "http://cran.us.r-project.org")
}
if (!("treemap" %in% installed.packages())) {
    install.packages("treemap", repos = "http://cran.us.r-project.org")
}
if (!("colorspace" %in% installed.packages())) {
    install.packages("colorspace", repos = "http://R-Forge.R-project.org")
}
if (!("colorblindr" %in% installed.packages())) {
    devtools::install_github("clauswilke/colorblindr")
}
if (!("mapview" %in% installed.packages())) {
    install.packages("mapview", repos = "http://cran.us.r-project.org")
}
if (!("d3treeR" %in% installed.packages())) {
    devtools::install_github("timelyportfolio/d3treeR")
}


