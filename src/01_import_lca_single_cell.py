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
    / "rawData_human"
    / "rawData_human"
    / "countTable_human"
)

ANNOT = RAW / "annot_humanAll.csv"
OUT = Path("data/interim/lca_human_all_imported.h5ad")
OUT.parent.mkdir(parents=True, exist_ok=True)


def read_table_gz(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing file: {path}")
    return pd.read_csv(path, sep="\t", header=None, compression="gzip")


def strip_10x_suffix(values: pd.Index) -> pd.Index:
    return pd.Index(pd.Series(values.astype(str)).str.replace(r"-1$", "", regex=True))


def read_lca_count_table(count_dir: Path) -> ad.AnnData:
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

    # Standard 10x-style Matrix Market is features x barcodes.
    # AnnData wants cells/barcodes x genes/features.
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
    adata.obs["dataset"] = "Liver Cell Atlas human single-cell RNA"

    return adata


def merge_annotations(adata: ad.AnnData) -> ad.AnnData:
    if not ANNOT.exists():
        print(f"No annotation file found at {ANNOT}; continuing without annotations.")
        return adata

    annot = pd.read_csv(ANNOT)
    print("Annotation columns:")
    print(list(annot.columns))
    print("Annotation preview:")
    print(annot.head())

    candidate_cols = []
    for col in annot.columns:
        low = col.lower()
        if "cell" in low or "barcode" in low or "sample" in low:
            candidate_cols.append(col)

    if "Unnamed: 0" in annot.columns:
        candidate_cols.insert(0, "Unnamed: 0")

    # As a fallback, test all columns. This is slower but safer for unknown atlas tables.
    candidate_cols = list(dict.fromkeys(candidate_cols + list(annot.columns)))

    obs_names = pd.Index(adata.obs_names.astype(str))

    best_col = None
    best_common = pd.Index([])

    for col in candidate_cols:
        vals = pd.Index(annot[col].astype(str))
        common = obs_names.intersection(vals)
        if len(common) > len(best_common):
            best_col = col
            best_common = common

    if best_col is not None and len(best_common) > 0:
        print(f"Using annotation key column: {best_col}")
        print(f"Exact common cells: {len(best_common)} / {adata.n_obs}")

        annot = annot.set_index(best_col)
        annot = annot.loc[~annot.index.duplicated(keep="first")]

        adata = adata[best_common].copy()
        adata.obs = adata.obs.join(annot, how="left")
        return adata

    print("No exact barcode/cell match found. Trying without trailing '-1' suffix...")

    obs_norm = strip_10x_suffix(obs_names)
    obs_norm_map = pd.Series(obs_norm.values, index=obs_names)

    for col in candidate_cols:
        vals_original = pd.Index(annot[col].astype(str))
        vals_norm = strip_10x_suffix(vals_original)
        common_norm = pd.Index(obs_norm).intersection(vals_norm)

        if len(common_norm) > 0:
            print(f"Using normalized annotation key column: {col}")
            print(f"Normalized common cells: {len(common_norm)} / {adata.n_obs}")

            annot = annot.copy()
            annot["_join_key_norm"] = vals_norm.values
            annot = annot.set_index("_join_key_norm")
            annot = annot.loc[~annot.index.duplicated(keep="first")]

            keep_mask = obs_norm_map.isin(common_norm)
            adata = adata[keep_mask.values].copy()

            join_keys = obs_norm_map.loc[adata.obs_names].values
            ann = annot.loc[join_keys].copy()
            ann.index = adata.obs_names

            adata.obs = adata.obs.join(ann, how="left")
            return adata

    print("Could not merge annotations automatically.")
    print("First adata barcodes:")
    print(list(adata.obs_names[:5]))
    print("Continuing without annotations.")
    return adata


def print_summary(adata: ad.AnnData, label: str) -> None:
    print("\n" + "=" * 80)
    print(label)
    print("=" * 80)
    print(adata)
    print(f"Cells/spots: {adata.n_obs}")
    print(f"Genes/features: {adata.n_vars}")
    print("obs columns:", list(adata.obs.columns))
    print("var columns:", list(adata.var.columns))
    print("=" * 80 + "\n")


if not COUNT_DIR.exists():
    raise SystemExit(f"Expected count table folder does not exist: {COUNT_DIR}")

adata = read_lca_count_table(COUNT_DIR)
adata = merge_annotations(adata)

print_summary(adata, "LCA human single-cell RNA imported")

adata.write_h5ad(OUT, compression="gzip")
print(f"Wrote {OUT}")
