from argparse import ArgumentParser
from pathlib import Path

import anndata as ad
import loompy
import numpy as np
from scipy import sparse

parser = ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--output", required=True)
parser.add_argument("--drug", required=True)
parser.add_argument("--timepoint", required=True)
parser.add_argument("--drug-column", required=True)
parser.add_argument("--timepoint-column", required=True)
parser.add_argument("--celltype-column", required=True)
parser.add_argument("--counts-layer", required=True)
args = parser.parse_args()

adata = ad.read_h5ad(args.input, backed="r")

for column in [
    args.drug_column,
    args.timepoint_column,
    args.celltype_column,
]:
    if column not in adata.obs:
        raise KeyError(f"Missing obs column: {column}")

if args.counts_layer not in adata.layers:
    raise KeyError(f"Missing layer: {args.counts_layer}")

drug_values = (
    adata.obs[args.drug_column]
    .astype(str)
    .str.lower()
)

timepoint_values = (
    adata.obs[args.timepoint_column]
    .astype(str)
    .str.lower()
)

mask = (
    drug_values.eq(args.drug.lower())
    & timepoint_values.eq(args.timepoint.lower())
)

print(f"Selected cells: {int(mask.sum())}")

if mask.sum() == 0:
    raise ValueError(
        f"No cells found for {args.drug}, {args.timepoint}"
    )

subset = adata[mask].to_memory()
adata.file.close()

if not subset.obs_names.is_unique:
    raise ValueError("Cell IDs are not unique")

if not subset.var_names.is_unique:
    raise ValueError("Gene names are not unique")

counts = subset.layers[args.counts_layer]

if sparse.issparse(counts):
    counts = counts.T.tocsc()
else:
    counts = np.asarray(counts).T

row_attrs = {
    "Gene": subset.var_names.astype(str).to_numpy()
}

col_attrs = {
    "CellID": subset.obs_names.astype(str).to_numpy(),
    "Drug": subset.obs[args.drug_column].astype(str).to_numpy(),
    "Timepoint": subset.obs[args.timepoint_column].astype(str).to_numpy(),
    "CellType": subset.obs[args.celltype_column].astype(str).to_numpy(),
}

Path(args.output).parent.mkdir(parents=True, exist_ok=True)

loompy.create(
    args.output,
    counts,
    row_attrs,
    col_attrs,
)

print(f"Created: {args.output}")
print(f"Shape: {subset.n_obs} cells x {subset.n_vars} genes")