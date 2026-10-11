"""Merge every cell from per-drug Looms through Snakemake's script directive."""

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


def prepare_pooled_loom(input_files, output_file):
    if not input_files:
        raise ValueError("No input Loom files supplied")

    reference_genes = None
    inputs = []
    seen_drugs = set()

    # Check the shared gene space before creating any output.
    for filename in input_files:
        path = Path(filename)
        drug = path.parent.name
        if drug in seen_drugs:
            raise ValueError(f"Duplicate drug input: {drug}")
        seen_drugs.add(drug)

        with loompy.connect(str(path), mode="r", validate=False) as ds:
            if "Gene" not in ds.ra:
                raise KeyError(f"{path}: missing row attribute Gene")
            for key in ("CellID", "CellType"):
                if key not in ds.ca:
                    raise KeyError(f"{path}: missing column attribute {key}")

            genes = as_strings(ds.ra["Gene"])
            cells = as_strings(ds.ca["CellID"])
            cell_types = as_strings(ds.ca["CellType"])
            if ds.shape[0] == 0 or ds.shape[1] == 0:
                raise ValueError(f"{path}: empty expression matrix")
            if not pd.Index(genes).is_unique:
                raise ValueError(f"{path}: duplicate gene names")
            if not pd.Index(cells).is_unique:
                raise ValueError(f"{path}: duplicate CellID values")

            if reference_genes is None:
                reference_genes = genes
            elif not np.array_equal(genes, reference_genes):
                raise ValueError(
                    f"{path}: Gene names/order differ from the first Loom; "
                    "all input Looms must use the same genes in the same order"
                )

            invalid = np.isin(
                np.char.lower(np.char.strip(cell_types)),
                ["", "nan", "none", "<na>"],
            )
            cell_types = cell_types.astype(object)
            cell_types[invalid] = "Unknown"
            for cell_type in sorted(np.unique(cell_types)):
                count = np.count_nonzero(cell_types == cell_type)
                print(
                    f"{drug} / {cell_type}: keeping all {count} cells",
                    flush=True,
                )

            inputs.append((path, drug, cells, cell_types))

    output = Path(output_file)
    if output.resolve() in {Path(p).resolve() for p in input_files}:
        raise ValueError("Output must differ from every input Loom")
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(output.name + ".tmp.loom")
    total_cells = 0

    try:
        with loompy.new(str(temporary)) as pooled:
            for path, drug, cells, cell_types in inputs:
                with loompy.connect(str(path), mode="r", validate=False) as ds:
                    # Bound merge memory without discarding any cells.
                    for start in range(0, len(cells), 128):
                        stop = min(start + 128, len(cells))
                        matrix = ds[:, start:stop]
                        if not np.isfinite(matrix).all() or (matrix < 0).any():
                            raise ValueError(f"{path}: nonfinite/negative expression")
                        attrs = {
                            "CellID": np.asarray([
                                f"{drug}::{cell}" for cell in cells[start:stop]
                            ]),
                            "OriginalCellID": cells[start:stop],
                            "CellType": cell_types[start:stop],
                            "Drug": np.full(stop - start, drug, dtype=object),
                        }
                        pooled.add_columns(
                            matrix,
                            col_attrs=attrs,
                            row_attrs={"Gene": reference_genes},
                        )
                        total_cells += stop - start
            pooled.attrs["PoolingMode"] = "all_cells"
        os.replace(temporary, output)
    finally:
        if temporary.exists():
            temporary.unlink()

    print(f"Created: {output}", flush=True)
    print(f"Shape: {len(reference_genes)} genes x {total_cells} cells", flush=True)


if __name__ == "__main__":
    prepare_pooled_loom(
        list(snakemake.input),
        snakemake.output[0],
    )
