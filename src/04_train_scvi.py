from pathlib import Path
import scanpy as sc
import scvi
import matplotlib.pyplot as plt
from utils import load_yaml, ensure_dirs

config = load_yaml("config/project_config.yaml")
ensure_dirs(config)
scvi_cfg = config["scvi"]

IN = Path("data/processed/lca_human_marker_scored.h5ad")
OUT = Path("data/processed/lca_human_scvi.h5ad")
MODEL_DIR = Path(config["paths"]["models_dir"]) / "scvi_lca_human"
FIG = Path(config["paths"]["figures_dir"])

adata = sc.read_h5ad(IN)
if "counts" not in adata.layers:
    raise ValueError("Expected raw counts in adata.layers['counts']")

batch_key = None
for candidate in scvi_cfg["batch_key_candidates"]:
    if candidate in adata.obs.columns and adata.obs[candidate].nunique() > 1:
        batch_key = candidate
        break
print("Using batch_key:", batch_key)

scvi.settings.seed = config["random_seed"]
scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key=batch_key)

model = scvi.model.SCVI(
    adata,
    n_latent=scvi_cfg["n_latent"],
    n_layers=scvi_cfg["n_layers"],
    gene_likelihood=scvi_cfg["gene_likelihood"],
)
model.train(max_epochs=scvi_cfg["max_epochs"])

adata.obsm["X_scVI"] = model.get_latent_representation()
sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=15)
sc.tl.umap(adata, random_state=config["random_seed"])
sc.tl.leiden(adata, resolution=0.6, key_added="leiden_scvi")

color_cols = ["leiden_scvi", "leiden_scanpy"]
for c in [batch_key, "cell_type", "celltype", "annotation", "dataset"]:
    if c and c in adata.obs.columns and c not in color_cols:
        color_cols.append(c)

sc.pl.umap(adata, color=color_cols[:5], wspace=0.4, show=False)
plt.savefig(FIG / "05_scvi_umap_overview.png", bbox_inches="tight")
plt.close()

plot_scores = [
    "score_hepatocyte",
    "score_cholangiocyte_biliary",
    "score_endothelial_general",
    "score_macrophage_kupffer_like",
    "score_stellate_fibroblast",
    "score_t_cell",
]
plot_scores = [c for c in plot_scores if c in adata.obs.columns]
if plot_scores:
    sc.pl.umap(adata, color=plot_scores, cmap="viridis", show=False, wspace=0.4)
    plt.savefig(FIG / "06_scvi_umap_marker_scores.png", bbox_inches="tight")
    plt.close()

MODEL_DIR.mkdir(parents=True, exist_ok=True)
model.save(MODEL_DIR, overwrite=True)
adata.write_h5ad(OUT, compression="gzip")
print(f"Wrote {OUT}")
