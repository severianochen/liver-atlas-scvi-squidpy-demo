from pathlib import Path
import scanpy as sc
import matplotlib.pyplot as plt
from utils import load_yaml, ensure_dirs, print_adata_summary

config = load_yaml("config/project_config.yaml")
ensure_dirs(config)
sc_cfg = config["single_cell"]

IN = Path("data/interim/lca_human_all_imported.h5ad")
OUT = Path("data/processed/lca_human_qc_scanpy.h5ad")
FIG = Path(config["paths"]["figures_dir"])

sc.settings.verbosity = 2
sc.settings.set_figure_params(dpi=120, facecolor="white", figsize=(6, 5))

adata = sc.read_h5ad(IN)
adata.var_names_make_unique()
if "counts" not in adata.layers:
    adata.layers["counts"] = adata.X.copy()

adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

sc.pl.violin(adata, ["n_genes_by_counts", "total_counts", "pct_counts_mt"], jitter=0.4, multi_panel=True, show=False)
plt.savefig(FIG / "01_qc_violin_before_filtering.png", bbox_inches="tight")
plt.close()

print("Before filtering:", adata.shape)
adata = adata[adata.obs["n_genes_by_counts"] >= sc_cfg["min_genes"]].copy()
adata = adata[adata.obs["n_genes_by_counts"] <= sc_cfg["max_genes"]].copy()
adata = adata[adata.obs["pct_counts_mt"] <= sc_cfg["max_pct_mito"]].copy()
sc.pp.filter_genes(adata, min_cells=sc_cfg["min_cells_per_gene"])
print("After filtering:", adata.shape)

sc.pp.normalize_total(adata, target_sum=sc_cfg["target_sum"])
sc.pp.log1p(adata)
adata.raw = adata
sc.pp.highly_variable_genes(adata, n_top_genes=sc_cfg["n_top_genes"], subset=False)

adata_hvg = adata[:, adata.var["highly_variable"]].copy()
sc.pp.scale(adata_hvg, max_value=10)
sc.tl.pca(adata_hvg, svd_solver="arpack")
sc.pp.neighbors(adata_hvg, n_neighbors=sc_cfg["n_neighbors"], n_pcs=sc_cfg["n_pcs"])
sc.tl.umap(adata_hvg, random_state=config["random_seed"])
sc.tl.leiden(adata_hvg, resolution=sc_cfg["leiden_resolution"], key_added="leiden_scanpy")

adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"]
adata.obsm["X_umap"] = adata_hvg.obsm["X_umap"]
adata.obsp["distances"] = adata_hvg.obsp["distances"]
adata.obsp["connectivities"] = adata_hvg.obsp["connectivities"]
adata.obs["leiden_scanpy"] = adata_hvg.obs["leiden_scanpy"].astype(str)

color_cols = ["leiden_scanpy"]
for c in ["cell_type", "celltype", "annotation", "cluster", "sample", "donor", "dataset"]:
    if c in adata.obs.columns:
        color_cols.append(c)

sc.pl.umap(adata, color=color_cols[:4], wspace=0.4, show=False)
plt.savefig(FIG / "02_scanpy_umap_overview.png", bbox_inches="tight")
plt.close()

print_adata_summary(adata, "QC + Scanpy baseline")
adata.write_h5ad(OUT, compression="gzip")
print(f"Wrote {OUT}")
