# Run this file to install all necessary packages

# Package names
packages <- c(
  "Racmacs",
  "tidyverse",
  "remotes",
  "knitr",
  "here"
)

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# The ablandscapes package needs to be installed from github
remotes::install_github("acorg/ablandscapes")