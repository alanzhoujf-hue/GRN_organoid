from argparse import ArgumentParser
import os
from pathlib import Path

import loompy
import pandas as pd
from pyscenic.rss import regulon_specificity_scores

parser = ArgumentParser()
parser.add_argument("--auc", required=True)
parser.add_argument("--loom", required=True)
parser.add_argument("--output", required=True)

args = parser.parse_args()

auc = pd.read_csv(args.auc, index_col=0)
auc.index = auc.index.astype(str)
auc = auc.apply(pd.to_numeric, errors="raise")


with loompy.connect(args.loom, mode="r", validate=False) as ds:
    cell_ids = pd.Index(ds.ca["CellID"].astype(str))
    celltypes = pd.Series(
        ds.ca["CellType"].astype(str),
        index=cell_ids,
        name="CellType",
    )

if not auc.index.is_unique:
    raise ValueError("AUC cell IDs are not unique")

if not celltypes.index.is_unique:
    raise ValueError("Loom cell IDs are not unique")

common = auc.index.intersection(celltypes.index)

if len(common) == 0:
    raise ValueError("No matching cell IDs between AUC and Loom")

auc = auc.loc[common]
celltypes = celltypes.loc[common]

valid = (
    celltypes.notna()
    & celltypes.str.strip().ne("")
    & celltypes.str.lower().ne("nan")
)

auc = auc.loc[valid]
celltypes = celltypes.loc[valid]

rss = regulon_specificity_scores(
    auc,
    celltypes,
)

Path(args.output).parent.mkdir(parents=True, exist_ok=True)
rss.to_csv(args.output)

print(f"Cells used: {len(auc)}")
print(f"Cell types: {celltypes.nunique()}")
print(f"Regulons: {auc.shape[1]}")
print(f"Created: {args.output}")

