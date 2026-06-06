from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
from scipy import sparse

from utils import load_yaml

config = load_yaml("config/project_config.yaml")
markers = load_yaml("config/marker_genes.yaml")

IN = Path("data/raw/holdout/holdout_liver.h5ad")
OUT = Path("data/processed/holdout_liver_pipeline_result.h5ad")
FIG = Path("figures")
FIG.mkdir(exist_ok=True)
OUT.parent.mkdir(parents=True, exist_ok=True)


def sample_values(matrix, n=10000):
    if sparse.issparse(matrix):
        vals = matrix.data[: min(n, matrix.data.size)]
    else:
        arr = np.asarray(matrix)
        vals = arr.ravel()[: min(n, arr.size)]
    return vals


def is_integer_like(matrix):
    vals = sample_values(matrix)
    vals = vals[np.isfinite(vals)]
    if len(vals) == 0:
        return False
    return bool(np.allclose(vals, np.round(vals)) and np.min(vals) >= 0)


def add_simple_marker_scores(adata, marker_dict):
    """
    Score marker programs as simple mean expression of present marker genes.
    This is safer for preprocessed holdout objects than Scanpy score_genes(),
    because we are not pretending to have raw counts.
    """
    for program, genes in marker_dict.items():
        present = [g for g in genes if g in adata.var_names]
        if len(present) < 2:
            print(f"Skipping {program}: only {len(present)} marker(s) present")
            continue

        x = adata[:, present].X
        if sparse.issparse(x):
            score = np.asarray(x.mean(axis=1)).ravel()
        else:
            score = np.asarray(x).mean(axis=1)

        score = np.nan_to_num(score, nan=0.0, posinf=0.0, neginf=0.0)
        adata.obs[f"score_{program}"] = score
        print(f"Scored {program}: {len(present)} markers")


if not IN.exists():
    raise SystemExit(f"Missing {IN}")

adata = sc.read_h5ad(IN)
adata.var_names_make_unique()

print("Loaded holdout:")
print(adata)
print("obs columns:", list(adata.obs.columns))
print("var columns:", list(adata.var.columns))
print("layers:", list(adata.layers.keys()))
print("obsm:", list(adata.obsm.keys()))
print("raw present:", adata.raw is not None)

max_cells = 30000
if adata.n_obs > max_cells:
    rng = np.random.default_rng(config["random_seed"])
    idx = rng.choice(adata.n_obs, size=max_cells, replace=False)
    adata = adata[idx].copy()
    print(f"Subsampled holdout to {adata.n_obs} cells")

has_counts = False

for layer in ["counts", "raw_counts", "count", "umi", "UMI", "X_counts"]:
    if layer in adata.layers and is_integer_like(adata.layers[layer]):
        has_counts = True
        print(f"Found integer-like count layer: {layer}")
        break

if is_integer_like(adata.X):
    has_counts = True
    print("adata.X looks integer-count-like")
else:
    print("adata.X is not integer-count-like; treating object as already processed")

if has_counts:
    print("This script currently stops here because this holdout branch was expected to be preprocessed.")
    print("Use the main LCA scripts for raw-count scVI training.")
else:
    # Keep existing processed embedding/neighborhoods.
    if "X_umap" not in adata.obsm:
        if "X_pca" in adata.obsm:
            print("No X_umap found, computing UMAP from existing X_pca")
            sc.pp.neighbors(adata, use_rep="X_pca")
            sc.tl.umap(adata, random_state=config["random_seed"])
        else:
            raise SystemExit("No X_umap or X_pca found; cannot make holdout overview safely.")

    add_simple_marker_scores(adata, markers)

    if "louvain" in adata.obs.columns:
        cluster_key = "louvain"
    elif "leiden" in adata.obs.columns:
        cluster_key = "leiden"
    else:
        cluster_key = None

    plot_cols = []
    for c in [cluster_key, "NormalvsTumor", "patientno", "ViralvsNonViral"]:
        if c is not None and c in adata.obs.columns:
            plot_cols.append(c)

    for c in [
        "score_hepatocyte",
        "score_cholangiocyte_biliary",
        "score_macrophage_kupffer_like",
        "score_t_cell",
        "score_nk_cell_cytotoxic",
        "score_hcc_or_fetal_like_program_careful",
    ]:
        if c in adata.obs.columns:
            plot_cols.append(c)

    print("Plot columns:", plot_cols)

    sc.pl.umap(
        adata,
        color=plot_cols,
        show=False,
        wspace=0.4,
        ncols=3,
    )
    plt.savefig(FIG / "10_holdout_marker_overview.png", bbox_inches="tight", dpi=200)
    plt.close()

    adata.uns["holdout_note"] = (
        "GSE156625 holdout object appears preprocessed and lacks an integer-like counts layer. "
        "scVI was therefore skipped; existing UMAP/louvain metadata and marker-program means were used."
    )

adata.write_h5ad(OUT, compression="gzip")
print(f"Wrote {OUT}")
print("Wrote figures/10_holdout_marker_overview.png")
