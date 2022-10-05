#!/usr/bin/env python3

from pathlib import Path
from typing import Iterable
from argparse import ArgumentParser
import os
import anndata

def find_files(directory: Path, pattern: str) -> Iterable[Path]:
    for dirpath_str, dirnames, filenames in os.walk(directory):
        dirpath = Path(dirpath_str)
        for filename in filenames:
            filepath = dirpath / filename
            if filepath.match(pattern):
                yield filepath

def main(data_directory: Path):
    dcc_files = find_files(data_directory, "*.dcc")
    pkc_files = find_files(data_directory, "*.pkc")
    #get obs
    #get vars
    #initialize adata
    #fill adata
    #write adata

if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('data_directory', type=Path)
    p.add_argument("--enable-manhole", action="store_true")
    args = p.parse_args()

    if args.enable_manhole:
        import manhole

        manhole.install(activate_on="USR1")

    main(args.data_directorye)