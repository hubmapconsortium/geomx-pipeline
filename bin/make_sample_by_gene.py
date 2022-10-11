#!/usr/bin/env python3

from pathlib import Path
from typing import Iterable, Dict
from argparse import ArgumentParser
import os
import anndata
import json
import pandas as pd

GENE_MAPPING_DIRECTORIES = [
    Path(__file__).parent.parent / 'data',
    Path('/opt/data'),
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
    temp_df = pd.DataFrame(adata.X.todense(), index=adata.obs.index, columns=adata.var.index)
    aggregated = temp_df.groupby(level=0, axis=1).sum()
    adata = anndata.AnnData(aggregated, obs=adata.obs)
    adata.var.index = [gene_mapping[var] for var in adata.var.index]
    reverse_dict = read_gene_mapping(reverse=True)
    adata.var["gene_symbol"] = pd.Series([reverse_dict[ensembl_id] for ensembl_id in adata.var.index])
    adata.obsm = obsm
    adata.var_names_make_unique()
    return adata

def get_probe_to_gene_symbol_mapping(pkc_file):
    with open(pkc_file) as f:
        target_dict = json.load(f)
    probe_to_gene_symbol_mapping = {entry["RTS_ID"]: entry["SystematicName"] for entry in target_dict["targets"]}
    return probe_to_gene_symbol_mapping

def get_code_summary(dcc_file):
    with open(dcc_file) as f:
        text = f.read()
    text_split = text.split('<Code_Summary')
    assert len(text_split) == 3
    code_summary = text_split[1]
    code_summary = code_summary.strip('>')
    code_summary_lines = code_summary.split('\n')
    return {code_summary_line.split(',')[0]:int(code_summary_line.split(',')[1].strip()) for code_summary_line in code_summary_lines}

def get_aoi_id(dcc_file):
    return dcc_file.stem

def main(data_directory: Path):
    pkc_file = find_files(data_directory, "*.pkc")
    probe_to_gene_mapping = get_probe_to_gene_symbol_mapping(pkc_file)
    dcc_files = find_files(data_directory, "*.dcc")
    code_summaries = {get_aoi_id(dcc_file): get_code_summary(dcc_file) for dcc_file in dcc_files}
    obs = list(code_summaries.keys())
    vars = list(probe_to_gene_mapping.values())
    adata = anndata.AnnData(obs=obs, vars=vars)
    for aoi_id in code_summaries:
        for probe_id in code_summaries[aoi_id]:
            gene_symbol = probe_to_gene_mapping[probe_id]
            adata[aoi_id, gene_symbol] = code_summaries[aoi_id][probe_id]

    adata = map_gene_ids(adata)
    adata.write_h5ad('cell_by_gene.h5ad')

if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('data_directory', type=Path)
    p.add_argument("--enable-manhole", action="store_true")
    args = p.parse_args()

    if args.enable_manhole:
        import manhole

        manhole.install(activate_on="USR1")

    main(args.data_directorye)