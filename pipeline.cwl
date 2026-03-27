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
  ome_tiff:
    outputSource: adjust-ometiff-channels/ome_tiff
    type: File

steps:
  make-sample-by-gene:
    in:
      data_directory:
        source: data_directory
      enable_manhole:
        source: enable_manhole
    out:
      - h5ad_file
    run: steps/make-sample-by-gene.cwl
    label: "Converts several dcc files into an annotated sample by gene matrix"
  adjust-ometiff-channels:
    in:
      data_directory:
        source: data_directory
    out:
      - ome_tiff
    run: steps/adjust-ometiff-channels.cwl
    label: "Write new OME-TIFF with adjusted channel names"
