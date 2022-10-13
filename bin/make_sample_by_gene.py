#!/usr/bin/env python3

from pathlib import Path
from typing import Iterable, Dict
from argparse import ArgumentParser
import os
import anndata
import json
import pandas as pd
import numpy as np

GENE_MAPPING_DIRECTORIES = [
    Path(__file__).parent.parent / 'data',
    Path('/opt/data'),
    Path('/opt/'),

]

def find_files(directory: Path, pattern: str) -> Iterable[Path]:
    for dirpath_str, dirnames, filenames in os.walk(directory):
        dirpath = Path(dirpath_str)
        for filename in filenames:
            filepath = dirpath / filename
            if filepath.match(pattern):
                yield filepath


def read_gene_mapping(reverse=False) -> Dict[str, str]:
    """
    Try to find the Ensembl to HUGO symbol mapping, with paths suitable
    for running this script inside and outside a Docker container.
    :return:
    """

    filename = 'ensembl_to_symbol.json' if reverse else 'symbol_to_ensembl.json'

    for directory in GENE_MAPPING_DIRECTORIES:
        mapping_file = directory / filename
        if mapping_file.is_file():
            with open(mapping_file) as f:
                return json.load(f)
    message_pieces = ["Couldn't find Ensembl → HUGO mapping file. Tried:"]
    message_pieces.extend(f'\t{path}' for path in GENE_MAPPING_DIRECTORIES)
    raise ValueError('\n'.join(message_pieces))

def map_gene_ids(adata):
    obsm = adata.obsm
    gene_mapping = read_gene_mapping()
    keep_vars = [gene in gene_mapping for gene in adata.var.index]
    adata = adata[:, keep_vars]
    temp_df = pd.DataFrame(adata.X, index=adata.obs.index, columns=adata.var.index)
    aggregated = temp_df.groupby(level=0, axis=1).sum()
    adata = anndata.AnnData(aggregated, obs=adata.obs)
    adata.var.index = [gene_mapping[var] for var in adata.var.index]
    reverse_dict = read_gene_mapping(reverse=True)
    adata.var["gene_symbol"] = pd.Series([reverse_dict[ensembl_id] for ensembl_id in adata.var.index], index=adata.var.index)
    adata.obsm = obsm
    adata.var_names_make_unique()
    return adata

def get_probe_to_gene_symbol_mapping(pkc_file):
    with open(pkc_file) as f:
        target_dict = json.load(f)
    probe_to_gene_symbol_mapping = {entry['Probes'][0]["RTS_ID"]: entry['DisplayName'] for entry in target_dict["Targets"]}
    for entry in target_dict["Targets"]:
        if len(entry["Probes"]) != 1:
            for item in entry["Probes"]:
                probe_to_gene_symbol_mapping[item["RTS_ID"]] = item["DisplayName"]

    return probe_to_gene_symbol_mapping

def parse_dcc(dcc_file, section_tag):
    """Takes a file path and a tag, returns a list of strings, one for each line in the section marked with that tag"""
    with open(dcc_file) as f:
        text = f.read()
    text_split = text.split(f'{section_tag}>')

    assert len(text_split) == 3
    section = text_split[1]
    section = section.strip('</')
    section_lines = section.split('\n')
    return section_lines

def lines_to_dict(section_lines):
    line_splits = [line.split(',') for line in section_lines if len(line.split(',')) == 2]
    line_dict = {line_split[0]:line_split[1] for line_split in line_splits}
    for attribute in line_dict:
        if line_dict[attribute].isnumeric():
            line_dict[attribute] = int(line_dict[attribute])
        elif attribute in ["UMI merge factor", "umiQ30", "rtsQ30"]:
            line_dict[attribute] = float(line_dict[attribute])

    return line_dict

def get_annotations(dcc_file):
    annotations = {}
    scan_attributes_lines = parse_dcc(dcc_file, "Scan_Attributes")
    scan_attribute_dict = lines_to_dict(scan_attributes_lines)
    ngs_processing_lines = parse_dcc(dcc_file, "NGS_Processing_Attributes")
    ngs_processing_dict = lines_to_dict(ngs_processing_lines)
    annotations.update(scan_attribute_dict)
    annotations.update(ngs_processing_dict)
    return annotations


def get_aoi_id(dcc_file):
    return dcc_file.stem

def main(data_directory: Path):
    pkc_files = list(find_files(data_directory, "*.pkc"))
    assert len(pkc_files) == 1
    pkc_file = pkc_files[0]
    probe_to_gene_mapping = get_probe_to_gene_symbol_mapping(pkc_file)
    dcc_files = list(find_files(data_directory, "*.dcc"))
    code_summaries = {get_aoi_id(dcc_file): lines_to_dict(parse_dcc(dcc_file, "Code_Summary")) for dcc_file in dcc_files}
    annotations = [get_annotations(dcc_file) for dcc_file in dcc_files]
    obs = pd.DataFrame(annotations)
    obs = obs.set_index("ID", drop=True, inplace=False)
    vars = list(probe_to_gene_mapping.values())
    vars = pd.DataFrame(index=vars)
    adata = anndata.AnnData(obs=obs, var=vars, X=np.zeros((len(obs.index), len(vars.index))))
    for aoi_id in code_summaries:
        for probe_id in code_summaries[aoi_id]:
            if aoi_id not in adata.obs.index:
                print(f"{aoi_id} not in obs_index")
            if probe_id not in probe_to_gene_mapping:
                print(f"{probe_id} not found in probe set")
            else:
                gene_symbol = probe_to_gene_mapping[probe_id]
                if gene_symbol not in adata.var.index:
                    print(f"{gene_symbol} not in var_index")
                adata[aoi_id, gene_symbol].X = code_summaries[aoi_id][probe_id]
    bool_mask = ["NegProbe" not in var for var in adata.var.index]
    adata = adata[:, bool_mask]

    #Drop negative probes
    adata = map_gene_ids(adata)
    adata.write_h5ad('sample_by_gene.h5ad')

if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('data_directory', type=Path)
    p.add_argument("--enable-manhole", action="store_true")
    args = p.parse_args()

    if args.enable_manhole:
        import manhole

        manhole.install(activate_on="USR1")

    main(args.data_directory)