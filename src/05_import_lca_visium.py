from pathlib import Path
import gzip

import anndata as ad
import pandas as pd
from scipy import sparse
from scipy.io import mmread

RAW = Path("data/raw/liver_cell_atlas")

COUNT_DIR = (
    RAW
    / "unpacked"
    / "rawData_humanVisium"
    / "rawData_humanVisium"
    / "countTable_humanVisium"
)

SAMPLECOMP = (
    RAW
    / "unpacked"
    / "rawData_humanVisium"
    / "rawData_humanVisium"
    / "sampleComp_humanVisium.txt"
)

ANNOT = RAW / "annot_humanVisium.csv"
OUT = Path("data/interim/lca_human_visium_imported.h5ad")
OUT.parent.mkdir(parents=True, exist_ok=True)


def read_table_gz(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing file: {path}")
    return pd.read_csv(path, sep="\t", header=None, compression="gzip")


def strip_10x_suffix(values: pd.Index) -> pd.Index:
    return pd.Index(pd.Series(values.astype(str)).str.replace(r"-1$", "", regex=True))


def read_lca_visium_count_table(count_dir: Path) -> ad.AnnData:
    matrix_path = count_dir / "matrix.mtx.gz"
    barcodes_path = count_dir / "barcodes.tsv.gz"
    features_path = count_dir / "features.tsv.gz"

    print(f"Reading matrix:   {matrix_path}")
    print(f"Reading barcodes: {barcodes_path}")
    print(f"Reading features: {features_path}")

    barcodes_df = read_table_gz(barcodes_path)
    features_df = read_table_gz(features_path)

    barcodes = barcodes_df.iloc[:, 0].astype(str)

    if features_df.shape[1] >= 2:
        gene_ids = features_df.iloc[:, 0].astype(str)
        gene_names = features_df.iloc[:, 1].astype(str)
    else:
        gene_ids = features_df.iloc[:, 0].astype(str)
        gene_names = gene_ids.copy()

    var = pd.DataFrame(index=pd.Index(gene_names, name="gene"))
    var["gene_id"] = gene_ids.values

    if features_df.shape[1] >= 3:
        var["feature_type"] = features_df.iloc[:, 2].astype(str).values

    obs = pd.DataFrame(index=pd.Index(barcodes, name="barcode"))

    with gzip.open(matrix_path, "rb") as handle:
        matrix = mmread(handle)

    if sparse.issparse(matrix):
        matrix = matrix.tocsr()
    else:
        matrix = sparse.csr_matrix(matrix)

    print(f"Raw matrix shape from file: {matrix.shape}")
    print(f"Number of barcodes: {len(obs)}")
    print(f"Number of features: {len(var)}")

    if matrix.shape == (len(var), len(obs)):
        X = matrix.T.tocsr()
    elif matrix.shape == (len(obs), len(var)):
        X = matrix.tocsr()
    else:
        raise SystemExit(
            "Matrix dimensions do not match features/barcodes.\n"
            f"matrix shape: {matrix.shape}\n"
            f"features:     {len(var)}\n"
            f"barcodes:     {len(obs)}"
        )

    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.var_names_make_unique()
    adata.layers["counts"] = adata.X.copy()
    adata.obs["dataset"] = "Liver Cell Atlas human Visium"

    return adata


def read_metadata_table(path: Path) -> pd.DataFrame | None:
    if not path.exists():
        print(f"Metadata file not found: {path}")
        return None

    print(f"\nReading metadata: {path}")

    if path.suffix == ".csv":
        df = pd.read_csv(path)
    else:
        df = pd.read_csv(path, sep="\t")

    print("Columns:")
    print(list(df.columns))
    print("Preview:")
    print(df.head())

    return df


def merge_metadata(adata: ad.AnnData, df: pd.DataFrame | None, label: str) -> ad.AnnData:
    if df is None:
        return adata

    candidate_cols = []
    for col in df.columns:
        low = col.lower()
        if (
            "barcode" in low
            or "spot" in low
            or "cell" in low
            or "sample" in low
            or low == "id"
            or low.endswith("_id")
        ):
            candidate_cols.append(col)

    if "Unnamed: 0" in df.columns:
        candidate_cols.insert(0, "Unnamed: 0")

    candidate_cols = list(dict.fromkeys(candidate_cols + list(df.columns)))
    obs_names = pd.Index(adata.obs_names.astype(str))

    best_col = None
    best_common = pd.Index([])

    for col in candidate_cols:
        vals = pd.Index(df[col].astype(str))
        common = obs_names.intersection(vals)
        if len(common) > len(best_common):
            best_col = col
            best_common = common

    if best_col is not None and len(best_common) > 0:
        print(f"\nMerging {label} using exact key column: {best_col}")
        print(f"Exact common spots: {len(best_common)} / {adata.n_obs}")

        meta = df.copy()
        meta = meta.set_index(best_col)
        meta = meta.loc[~meta.index.duplicated(keep="first")]
        adata.obs = adata.obs.join(meta, how="left", rsuffix=f"_{label}")
        return adata

    print(f"\nNo exact match for {label}. Trying without trailing '-1' suffix...")

    obs_norm = strip_10x_suffix(obs_names)
    obs_norm_map = pd.Series(obs_norm.values, index=obs_names)

    best_col = None
    best_common_norm = pd.Index([])

    for col in candidate_cols:
        vals_original = pd.Index(df[col].astype(str))
        vals_norm = strip_10x_suffix(vals_original)
        common_norm = pd.Index(obs_norm).intersection(vals_norm)

        if len(common_norm) > len(best_common_norm):
            best_col = col
            best_common_norm = common_norm

    if best_col is not None and len(best_common_norm) > 0:
        print(f"Merging {label} using normalized key column: {best_col}")
        print(f"Normalized common spots: {len(best_common_norm)} / {adata.n_obs}")

        meta = df.copy()
        meta["_join_key_norm"] = strip_10x_suffix(pd.Index(meta[best_col].astype(str))).values
        meta = meta.set_index("_join_key_norm")
        meta = meta.loc[~meta.index.duplicated(keep="first")]

        join_keys = obs_norm_map.loc[adata.obs_names].values
        meta2 = meta.reindex(join_keys).copy()
        meta2.index = adata.obs_names

        adata.obs = adata.obs.join(meta2, how="left", rsuffix=f"_{label}")
        return adata

    print(f"Could not merge {label} automatically.")
    return adata


def add_spatial_coordinates_if_possible(adata: ad.AnnData) -> ad.AnnData:
    possible_pairs = [
        ("pxl_col_in_fullres", "pxl_row_in_fullres"),
        ("imagecol", "imagerow"),
        ("image_col", "image_row"),
        ("x", "y"),
        ("X", "Y"),
        ("array_col", "array_row"),
        ("col", "row"),
    ]

    for x_col, y_col in possible_pairs:
        if x_col in adata.obs.columns and y_col in adata.obs.columns:
            coords = adata.obs[[x_col, y_col]].apply(pd.to_numeric, errors="coerce")

            if coords.notna().all(axis=None):
                adata.obsm["spatial"] = coords.to_numpy()
                print(f"Added adata.obsm['spatial'] from columns: {x_col}, {y_col}")
                return adata

    print("No obvious spatial coordinate columns found yet.")
    print("This import still works as expression + metadata.")
    print("If sampleComp_humanVisium has coordinate columns with different names, adapt add_spatial_coordinates_if_possible().")
    return adata


def print_summary(adata: ad.AnnData, label: str) -> None:
    print("\n" + "=" * 80)
    print(label)
    print("=" * 80)
    print(adata)
    print(f"Spots/cells: {adata.n_obs}")
    print(f"Genes/features: {adata.n_vars}")
    print("obs columns:", list(adata.obs.columns))
    print("var columns:", list(adata.var.columns))
    print("obsm keys:", list(adata.obsm.keys()))
    print("layers:", list(adata.layers.keys()))
    print("=" * 80 + "\n")


if not COUNT_DIR.exists():
    raise SystemExit(f"Expected count table folder does not exist: {COUNT_DIR}")

adata = read_lca_visium_count_table(COUNT_DIR)
samplecomp = read_metadata_table(SAMPLECOMP)
annot = read_metadata_table(ANNOT)

adata = merge_metadata(adata, samplecomp, "samplecomp")
adata = merge_metadata(adata, annot, "annot")
adata = add_spatial_coordinates_if_possible(adata)

print_summary(adata, "LCA human Visium imported")
adata.write_h5ad(OUT, compression="gzip")
print(f"Wrote {OUT}")
