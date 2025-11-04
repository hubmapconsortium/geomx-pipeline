#!/usr/bin/env python3

from pathlib import Path
from typing import Iterable, Dict
from argparse import ArgumentParser
import os
import anndata
import muon
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

def get_probe_mappings(pkc_files):
    probe_to_gene_symbol_mapping = {}
    probe_to_protein_mapping = {}
    for pkc_file in pkc_files:
        with open(pkc_file) as f:
            target_dict = json.load(f)

        probe_mapping = probe_to_gene_symbol_mapping if target_dict['AnalyteType'] == 'RNA' else probe_to_protein_mapping
        pkc_mapping = {entry['Probes'][0]["RTS_ID"]: entry['DisplayName'] for entry in target_dict["Targets"]}
        for entry in target_dict["Targets"]:
            if len(entry["Probes"]) != 1:
                for item in entry["Probes"]:
                    pkc_mapping[item["RTS_ID"]] = item["DisplayName"]
        probe_mapping.update(pkc_mapping)

    return probe_to_gene_symbol_mapping, probe_to_protein_mapping

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

def drop_negative_probes(adata):
    bool_mask = ["NegProbe" not in var for var in adata.var.index]
    neg_probes = [not val for val in bool_mask]
    neg_probe_counts_df = adata[:,neg_probes].to_df()
    adata = adata[:, bool_mask]
    for column in neg_probe_counts_df.columns:
        adata.obsm[column] = neg_probe_counts_df[column].to_numpy()
    return adata

def main(data_directory: Path):
    pkc_files = list(find_files(data_directory, "*.pkc"))
    probe_to_gene_mapping, probe_to_protein_mapping = get_probe_mappings(pkc_files)
    dcc_files = list(find_files(data_directory, "*.dcc"))
    code_summaries = {get_aoi_id(dcc_file): lines_to_dict(parse_dcc(dcc_file, "Code_Summary")) for dcc_file in dcc_files}
    annotations = [get_annotations(dcc_file) for dcc_file in dcc_files]
    obs = pd.DataFrame(annotations)
    obs = obs.set_index("ID", drop=True, inplace=False)

    rna_vars = list(probe_to_gene_mapping.values())
    rna_vars = pd.DataFrame(index=rna_vars)
    rna_adata = anndata.AnnData(obs=obs, var=rna_vars, X=np.zeros((len(obs.index), len(rna_vars.index))))

    protein_vars = list(probe_to_protein_mapping.values())
    protein_vars = pd.DataFrame(index=protein_vars)
    protein_adata = anndata.AnnData(obs=obs, var=protein_vars, X=np.zeros((len(obs.index), len(protein_vars.index))))

    for aoi_id in code_summaries:
        for probe_id in code_summaries[aoi_id]:
            if aoi_id in rna_adata.obs.index:
                if probe_id in probe_to_gene_mapping:
                    gene_symbol = probe_to_gene_mapping[probe_id]
                    rna_adata[aoi_id, gene_symbol].X = code_summaries[aoi_id][probe_id]
                elif probe_id in probe_to_protein_mapping:
                    protein_id = probe_to_protein_mapping[probe_id]
                    protein_adata[aoi_id, protein_id].X = code_summaries[aoi_id][probe_id]
                else:
                    print(f"{probe_id} not found in probe set")
            else:
                print(f"{aoi_id} not in obs_index")

    rna_adata = drop_negative_probes(rna_adata)
    protein_adata = drop_negative_probes(protein_adata)
    #Drop negative probes

    adata_dict = {}

    rna_adata = map_gene_ids(rna_adata)
    if len(rna_adata.var.index) != 0:
        adata_dict['rna'] = rna_adata
    if len(protein_adata.var.index) != 0:
        adata_dict['protein'] = protein_adata
    mudata = muon.MuData(adata_dict)
    mudata.write_h5mu('sample_by_gene.h5mu')

if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('data_directory', type=Path)
    p.add_argument("--enable-manhole", action="store_true")
    args = p.parse_args()

    if args.enable_manhole:
        import manhole

        manhole.install(activate_on="USR1")

    main(args.data_directory)