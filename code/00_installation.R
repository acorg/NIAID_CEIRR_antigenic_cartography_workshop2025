# Run this file to install all necessary packages

# Package names
packages <- c(
  "Racmacs",
  "tidyverse",
  "remotes",
  "knitr",
  "here",
  "seqinr",
  "r3js",
  "chromote",
  "htmlwidgets",
  "plyr"
)

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# The ablandscapes package needs to be installed from github
# You need to have developer tools installed as described in the README
remotes::install_github("acorg/ablandscapes")