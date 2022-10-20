## ----main,echo=TRUE,include=FALSE-------------------------------------------------------------------------
library(NanoStringNCTools)
library(GeomxTools)
library(optparse)
library(Seurat)

option_list <- list(
  make_option(
    c("-d", "--data_directory"),
    type = "character",
    default = NULL,
    help = "Data directory to use."
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$data_directory)) {
  print_help(opt_parser)
  stop("--data_directory argument must be supplied (input data directory).", call. = FALSE)
}

data_directory = "/var/lib/cwl/"

DCCFiles <- dir(data_directory, pattern=".dcc$", full.names=TRUE, recursive=TRUE)
PKCFiles <- dir(data_directory, pattern=".pkc$", full.names=TRUE, recursive=TRUE)
annotations_file = dir(data_directory, pattern=".xlsx$", full.names=TRUE, recursive=TRUE)

data <-readNanoStringGeoMxSet(dccFiles = DCCFiles,
                                      pkcFiles = PKCFiles, phenoDataFile=annotations_file,
                                      phenoDataSheet = "Annotations")

target_data <- aggregateCounts(data)

saveRDS(target_data, file="sample_by_gene.rds")

#@TODO:Find out what normalization if any is appropriate
#seurat_target_data <- as.Seurat(
#  target_data,
#)

#saveRDS(seurat_target_data, file="seurat_sample_by_gene.rds")
