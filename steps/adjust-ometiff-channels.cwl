cwlVersion: v1.1
class: CommandLineTool
label: Replaces channels in OME-TIFF metadata

requirements:
  DockerRequirement:
    dockerPull: hubmap/geomx-adjust-ometiff-channels

baseCommand: /opt/adjust_ometiff_channels.py

inputs:
  data_directory:
    type: Directory
    doc: Source dataset directory
    inputBinding:
      position: 1

outputs:
  ome_tiff_directory:
    type: Directory
    outputBinding:
      glob: "lab_processed"
    doc: OME-TIFF with updated channel names
