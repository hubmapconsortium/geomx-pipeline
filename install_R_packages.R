install.packages(
c(
  'optparse'
  ),
 Ncpus=6
)

# ArchR is designed to be run on Unix-based operating systems such as macOS and linux.
# ArchR is NOT supported on Windows or other operating systems.

# ArchR installation currently requires devtools and BiocManager for installation
# of GitHub and Bioconductor packages. Run the following commands to install the
# various dependencies used by ArchR:

# Docker image rocker/tidyverse includes devtools so no need to install it
# First, install devtools (for installing GitHub packages) if it isn’t already installed:
#if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")

# Then, install BiocManager (for installing bioconductor packages) if it isn’t already installed:
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(version = "3.13", ask=FALSE)

if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")

BiocManager::install("GeomxTools")

install.packages("remotes")
remotes::install_github("Nanostring-Biostats/GeomxTools")
