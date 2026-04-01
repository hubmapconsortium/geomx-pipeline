#!/usr/bin/env python3
import shlex
from argparse import ArgumentParser
from os import fspath
from pathlib import Path
from shutil import copy
from subprocess import run

import pandas as pd
import tifffile
from lxml import etree


def find_ome_tiff(data_dir: Path) -> Path:
    images_dir = data_dir / "lab_processed/images"
    ome_tiffs = list(images_dir.glob("*.ome.tiff"))
    if (l := len(ome_tiffs)) != 1:
        raise ValueError(f"Need 1 OME-TIFF in {images_dir}; found {l}")
    return ome_tiffs[0]


def main(data_dir: Path):
    ome_tiff = find_ome_tiff(data_dir)
    print("Found OME-TIFF:", ome_tiff)
    output_ometiff_path = Path(ome_tiff.relative_to(data_dir))
    output_ometiff_path.parent.mkdir(exist_ok=True, parents=True)
    maybe_markers_csv = data_dir / "raw/markers.csv"
    print("Copying", ome_tiff, "to", output_ometiff_path)
    copy(ome_tiff, output_ometiff_path)
    if maybe_markers_csv.is_file():
        print("Found markers CSV:", maybe_markers_csv)
        markers = pd.read_csv(maybe_markers_csv)
        image = tifffile.TiffFile(ome_tiff)
        tree = etree.fromstring(image.ome_metadata)
        namespaces = tree.nsmap.copy()
        namespaces["OME"] = namespaces[None]
        del namespaces[None]
        channel_nodes = tree.xpath("//OME:Pixels/OME:Channel", namespaces=namespaces)
        print("Old channel names:", [c.attrib["Name"] for c in channel_nodes])
        # handle CSVs with junk rows at bottom
        new_markers = list(markers["marker"][: len(channel_nodes)])
        print("Adjusting OME-XML with new channel names:", new_markers)
        for channel_node, new_marker in zip(channel_nodes, new_markers):
            # TODO: figure out a better attribute name
            channel_node.attrib["OrigName"] = channel_node.attrib["Name"]
            channel_node.attrib["Name"] = new_marker
        new_ome_xml = etree.tostring(tree, encoding="utf-8")
        with open("temp.xml", "wb") as f:
            f.write(new_ome_xml)
        # Hack: adjust image file as little as possible; the Python tifffile library
        # could adjust image plane order, and the data could get out of sync with
        # the original OME-XML metadata (which we also modified as little as possible).
        # So, call out to `tiffcomment` to do this adjustment with the new OME-XML.
        command = [
            "tiffcomment",
            "--set-file",
            "temp.xml",
            fspath(output_ometiff_path),
        ]
        print("Running:", shlex.join(command))
        run(command)
    else:
        print("No markers CSV; leaving OME-TIFF unchanged")


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("data_dir", type=Path)
    args = p.parse_args()

    main(args.data_dir)
