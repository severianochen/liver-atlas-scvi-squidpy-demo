from pathlib import Path
import re

import numpy as np
import pandas as pd
import scanpy as sc

IN = Path("data/interim/lca_human_visium_imported.h5ad")
OUT = Path("data/interim/lca_human_visium_imported_spatial.h5ad")

GEO_DIR = Path("data/raw/liver_cell_atlas/geo_visium_per_sample")

FALLBACK_SAMPLE_SUFFIX = {
    "JBO14": "1",
    "JBO15": "2",
    "JBO18": "3",
    "JBO19": "4",
    "JBO22": "5",
}


def normalize_sample_name(value: str) -> str:
    """
    JBO014 -> JBO14
    JBO015 -> JBO15
    """
    m = re.search(r"JBO0*(\d+)", value)
    if not m:
        raise ValueError(f"Could not parse sample name from: {value}")
    return f"JBO{int(m.group(1))}"


def sample_from_filename(path: Path) -> str:
    m = re.search(r"JBO0*\d+", path.name)
    if not m:
        raise ValueError(f"Could not find JBO sample ID in filename: {path}")
    return normalize_sample_name(m.group(0))


def read_tissue_positions(path: Path) -> pd.DataFrame:
    """
    Reads old Space Ranger tissue_positions_list.csv files without headers,
    and also tolerates newer tissue_positions.csv files with headers.
    """
    raw = pd.read_csv(path, header=None, compression="infer")

    if raw.shape[1] < 6:
        raise ValueError(f"Expected at least 6 columns in {path}, found {raw.shape[1]}")

    first_cell = str(raw.iloc[0, 0]).lower()

    if "barcode" in first_cell:
        df = pd.read_csv(path, header=0, compression="infer")
        rename = {}
        for col in df.columns:
            c = str(col).lower()
            if c in {"barcode", "spot", "spot_id"}:
                rename[col] = "barcode"
            elif c == "in_tissue":
                rename[col] = "in_tissue"
            elif c == "array_row":
                rename[col] = "array_row"
            elif c == "array_col":
                rename[col] = "array_col"
            elif c == "pxl_row_in_fullres":
                rename[col] = "pxl_row_in_fullres"
            elif c == "pxl_col_in_fullres":
                rename[col] = "pxl_col_in_fullres"
        df = df.rename(columns=rename)
    else:
        df = raw.iloc[:, :6].copy()
        df.columns = [
            "barcode",
            "in_tissue",
            "array_row",
            "array_col",
            "pxl_row_in_fullres",
            "pxl_col_in_fullres",
        ]

    required = [
        "barcode",
        "in_tissue",
        "array_row",
        "array_col",
        "pxl_row_in_fullres",
        "pxl_col_in_fullres",
    ]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns in {path}: {missing}")

    df = df[required].copy()
    df["barcode"] = df["barcode"].astype(str)

    for col in required[1:]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    return df


def infer_sample_suffix_map(adata) -> dict:
    """
    Uses the LCA spot IDs such as AAACTGCTGGCTCCAA-1_1 and obs['sample']
    to infer which sample corresponds to which suffix.
    """
    if "sample" not in adata.obs.columns:
        return {}

    suffix = adata.obs_names.to_series().str.extract(r"_([0-9]+)$", expand=False)
    tmp = pd.DataFrame(
        {
            "sample": adata.obs["sample"].astype(str),
            "suffix": suffix,
        },
        index=adata.obs_names,
    )

    tmp = tmp[
        tmp["sample"].notna()
        & tmp["suffix"].notna()
        & (tmp["sample"] != "nan")
        & (tmp["sample"] != "None")
    ]

    mapping = {}
    for sample, sub in tmp.groupby("sample"):
        if len(sub) > 0:
            mapping[sample] = sub["suffix"].value_counts().index[0]

    return mapping


if not IN.exists():
    raise SystemExit(f"Missing {IN}. Run python src/05_import_lca_visium.py first.")

files = sorted(GEO_DIR.rglob("*tissue_positions*.csv*"))
if not files:
    raise SystemExit(
        f"No tissue position files found under {GEO_DIR}. "
        "Check the tar extraction and filenames."
    )

print("Found tissue position files:")
for f in files:
    print(" ", f)

adata = sc.read_h5ad(IN)
print(f"\nInput AnnData: {adata.n_obs} spots × {adata.n_vars} genes")

inferred = infer_sample_suffix_map(adata)
sample_to_suffix = dict(FALLBACK_SAMPLE_SUFFIX)
sample_to_suffix.update(inferred)

print("\nSample-to-suffix mapping:")
for k, v in sorted(sample_to_suffix.items()):
    print(f"  {k} -> _{v}")

coord_tables = []

for path in files:
    sample = sample_from_filename(path)
    if sample not in sample_to_suffix:
        raise SystemExit(f"No suffix mapping available for {sample} from {path}")

    suffix = sample_to_suffix[sample]
    df = read_tissue_positions(path)

    df["spatial_sample"] = sample
    df["spot"] = df["barcode"].astype(str) + "_" + str(suffix)

    coord_tables.append(df.set_index("spot"))

coords = pd.concat(coord_tables, axis=0)
coords = coords[~coords.index.duplicated(keep="first")]

common = adata.obs_names.intersection(coords.index)
print(f"\nMatched coordinates: {len(common)} / {adata.n_obs}")

if len(common) < adata.n_obs:
    missing = adata.obs_names.difference(coords.index)
    print("\nExample missing AnnData spot IDs:")
    print(list(missing[:10]))

    if len(common) < 0.95 * adata.n_obs:
        raise SystemExit(
            "Too many missing coordinates. Stop here and inspect the barcode/suffix mapping."
        )

    print("\nWarning: subsetting AnnData to spots with coordinates.")
    adata = adata[common].copy()

coords = coords.loc[adata.obs_names]

for col in [
    "spatial_sample",
    "barcode",
    "in_tissue",
    "array_row",
    "array_col",
    "pxl_row_in_fullres",
    "pxl_col_in_fullres",
]:
    adata.obs[col] = coords[col].values

adata.obs["spatial_sample"] = adata.obs["spatial_sample"].astype("category")

# Scanpy/Squidpy convention: x = pixel column, y = pixel row.
adata.obsm["spatial"] = coords[
    ["pxl_col_in_fullres", "pxl_row_in_fullres"]
].to_numpy(dtype=float)

print("\nFinal spatial summary:")
print(adata)
print("obsm keys:", list(adata.obsm.keys()))
print("spatial shape:", adata.obsm["spatial"].shape)
print(adata.obs["spatial_sample"].value_counts().sort_index())

# Avoid AnnData write error:
# obs index was named "barcode", and we also store the raw Space Ranger barcode.
adata.obs.index.name = "spot"

adata.write_h5ad(OUT, compression="gzip")
print(f"\nWrote updated AnnData with spatial coordinates to {OUT}")
