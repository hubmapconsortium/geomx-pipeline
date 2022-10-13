## ----main,echo=TRUE,include=FALSE-------------------------------------------------------------------------
library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)
library(ggplot2)
library(DESeq2)
library(pheatmap)
library(knitr)
library(dplyr)
library(ggforce)
library(optparse)

option_list <- list(
  make_option(
    c("-d", "--data_directory"),
    type = "character",
    default = NULL,
    help = "Data directory to use."
  ),
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$data_directory)) {
  print_help(opt_parser)
  stop("--data_directory argument must be supplied (input data directory).", call. = FALSE)
}

DCCFiles <- dir(data_directory, pattern=".dcc$", full.names=TRUE)
PKCFiles <- dir(data_directory, pattern=".pkc$", full.names=TRUE)

data <-readNanoStringGeoMxSet(dccFiles = DCCFiles,
                                      pkcFiles = PKCFiles,)

target_data <- aggregateCounts(data)

filename = sprintf("%s.rds", data_directory)
save_rds(target_data, )