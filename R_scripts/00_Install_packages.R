# Run the script below to install the required packages
# Rtools is necessary for GOTMr and gotmtools: https://cran.r-project.org/bin/windows/Rtools/ 

install.packages('glmtools', repos = "https://owi.usgs.gov/R")
install.packages("devtools")
devtools::install_github("tadhg-moore/GOTMr")
devtools::install_github("tadhg-moore/gotmtools")
install.packages('ggplot2')
install.packages('ggpubr')
install.packages('reshape')

