cwlVersion: v1.1
class: CommandLineTool
label: segments each image in the directory for FTUs

requirements:
  DockerRequirement:
    dockerPull: hubmap/geomx-pipeline-r

baseCommand: [Rscript, /opt/make_sample_by_gene.R]

inputs:

  data_directory:
    type: Directory
    doc: Path to processed dataset directory
    inputBinding:
      position: 1
      prefix: --data_directory
      valueFrom: $(self.basename)

outputs:
  rds_files:
    type: File[]
    outputBinding:
      glob: "*.rds"
    doc: rds_files containing serialized Seurat object and serialized geoMX experiment objects