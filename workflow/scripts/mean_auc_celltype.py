"""Average cell-level AUCell scores by CellType via a Snakemake script rule."""

import os
from pathlib import Path

import loompy
import numpy as np
import pandas as pd


def as_strings(values):
    return np.asarray([
        value.decode("utf-8") if isinstance(value, bytes) else str(value)
        for value in values
    ])


def mean_auc_celltype(auc_file, loom_file, output_file):
    # Read IDs as strings to preserve leading zeros and other literal IDs.
    auc = pd.read_csv(auc_file, index_col=0, dtype=str)
    auc.index = auc.index.astype(str)
    auc = auc.apply(pd.to_numeric, errors="raise")
    if auc.empty:
        raise ValueError("AUC matrix is empty")
    if not auc.index.is_unique or not auc.columns.is_unique:
        raise ValueError("AUC cell IDs and regulon names must be unique")
    if not np.isfinite(auc.to_numpy()).all():
        raise ValueError("AUC matrix contains NaN or infinite values")

    with loompy.connect(str(loom_file), mode="r", validate=False) as ds:
        for key in ("CellID", "CellType"):
            if key not in ds.ca:
                raise KeyError(f"Missing Loom column attribute: {key}")
        celltypes = pd.Series(
            as_strings(ds.ca["CellType"]),
            index=pd.Index(as_strings(ds.ca["CellID"])),
            name="cell_type",
        )

    if not celltypes.index.is_unique:
        raise ValueError("Loom cell IDs must be unique")
    missing = auc.index.difference(celltypes.index)
    extra = celltypes.index.difference(auc.index)
    if len(missing) or len(extra):
        raise ValueError(
            "AUC and Loom cell IDs do not match: "
            f"{len(missing)} AUC IDs missing from Loom, "
            f"{len(extra)} Loom IDs missing from AUC"
        )

    celltypes = celltypes.reindex(auc.index)
    invalid = celltypes.str.strip().str.lower().isin(["", "nan", "none", "<na>"])
    celltypes.loc[invalid] = "Unknown"
    mean_auc = auc.groupby(celltypes, observed=True, sort=True).mean()
    mean_auc.index.name = "cell_type"

    output = Path(output_file)
    if output.resolve() in {Path(auc_file).resolve(), Path(loom_file).resolve()}:
        raise ValueError("Output must differ from the input files")
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(output.name + ".tmp")
    try:
        mean_auc.to_csv(temporary)
        os.replace(temporary, output)
    finally:
        if temporary.exists():
            temporary.unlink()

    print(f"Cells used: {len(auc)}", flush=True)
    print(f"Cell types: {len(mean_auc)}", flush=True)
    print(f"Regulons: {auc.shape[1]}", flush=True)
    print(f"Created: {output}", flush=True)


if __name__ == "__main__":
    mean_auc_celltype(
        snakemake.input.auc,
        snakemake.input.loom,
        snakemake.output[0],
    )
