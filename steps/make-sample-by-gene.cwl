cwlVersion: v1.1
class: CommandLineTool
label: Converts vendor output to H5AD

requirements:
  DockerRequirement:
    dockerPull: hubmap/geomx-pipeline:0.3

baseCommand: /opt/make_sample_by_gene.py

inputs:
  enable_manhole:
    label: "Whether to enable remote debugging via 'manhole'"
    type: boolean?
    inputBinding:
      position: 0

  data_directory:
    type: Directory
    doc: Path to processed dataset directory
    inputBinding:
      position: 1


outputs:
  h5ad_file:
    type: File
    outputBinding:
      glob: "sample_by_gene.h5ad"
    doc: sample by gene matrix in h5ad format
