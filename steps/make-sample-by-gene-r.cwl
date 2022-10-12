cwlVersion: v1.1
class: CommandLineTool
label: segments each image in the directory for FTUs

requirements:
  DockerRequirement:
    dockerPull: hubmap/geomx-pipeline-r
  DockerGpuRequirement: {}

baseCommand: /opt/make_sample_by_gene.R

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

  h5ad_file:
    type: File
    doc: Path to h5ad file
    inputBinding:
      position: 2

outputs:
  rds_files:
    type: File[]
    outputBinding:
      glob: "*.rds"
    doc: rds_files containing serialized Seurat object and serialized geoMX experiment objects