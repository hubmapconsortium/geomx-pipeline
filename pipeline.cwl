#!/usr/bin/env cwl-runner

class: Workflow
cwlVersion: v1.0
label: Pipeline for converting dcc files output by GeoMX into sample by gene matrices

inputs:

  enable_manhole:
    label: "Whether to enable remote debugging via 'manhole'"
    type: boolean?

  data_directory:
    label: "Path to directory containing dcc files and pkc file"
    type: Directory


outputs:

  h5ad_file:
    outputSource: make-sample-by-gene/h5ad_file
    type: File

  rds_files:
    outputSource: make-sample-by-gene-r/rds_files
    type: File[]

steps:

  - id: make-sample-by-gene
    in:
      - id: data_directory
        source: data_directory
      - id: enable_manhole
        source: enable_manhole

    out:
      - h5ad_file

    run: steps/make-sample-by-gene.cwl
    label: "Converts several dcc files into an annotated sample by gene matrix"

  - id: make-sample-by-gene-r
    in:
      - id: data_directory
        source: data_directory
      - id: h5ad_file
        source: h5ad_file
      - id: enable_manhole
        source: enable_manhole

    out:
      - rds_files

    run: steps/make-sample-by-gene-r.cwl
    label: "Converts several dcc files into an annotated sample by gene matrix"
