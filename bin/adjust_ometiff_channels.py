#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path
from shutil import copy

import bioio
import pandas as pd
from bioio.writers import OmeTiffWriter
from bioio_base.dimensions import DimensionNames


def find_ome_tiff(data_dir: Path) -> Path:
    images_dir = data_dir / "lab_processed/images"
    ome_tiffs = list(images_dir.glob("*.ome.tiff"))
    if (l := len(ome_tiffs)) != 1:
        raise ValueError(f"Need 1 OME-TIFF in {images_dir}; found {l}")
    return ome_tiffs[0]


def main(data_dir: Path):
    ome_tiff = find_ome_tiff(data_dir)
    print("Found OME-TIFF:", ome_tiff)
    output_ometiff_path = Path(ome_tiff.name)
    maybe_markers_csv = data_dir / "raw/markers.csv"
    if maybe_markers_csv.is_file():
        print("Found markers CSV:", maybe_markers_csv)
        markers = pd.read_csv(maybe_markers_csv)
        image = bioio.BioImage(ome_tiff)
        print("Old channel names:", image.channel_names)
        # handle CSVs with junk rows at bottom
        new_markers = list(markers["marker"][: len(image.channel_names)])
        print("Instantiating new image with channel names:", new_markers)
        new_image = bioio.BioImage(
            image.xarray_dask_data.assign_coords({DimensionNames.Channel: new_markers}),
            physical_pixel_sizes=image.physical_pixel_sizes,
        )
        print("Saving image to", output_ometiff_path)
        new_image.save(output_ometiff_path)
    else:
        print("No markers CSV; copying OME-TIFF to", output_ometiff_path)
        copy(ome_tiff, output_ometiff_path)


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("data_dir", type=Path)
    args = p.parse_args()

    main(args.data_dir)
