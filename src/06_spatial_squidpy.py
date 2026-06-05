from pathlib import Path
import scanpy as sc
import squidpy as sq
import matplotlib.pyplot as plt
from utils import load_yaml, ensure_dirs

config = load_yaml("config/project_config.yaml")
markers = load_yaml("config/marker_genes.yaml")
ensure_dirs(config)
sp_cfg = config["spatial"]

IN = Path("data/interim/lca_human_visium_imported_spatial.h5ad")
OUT = Path("data/processed/lca_human_visium_squidpy.h5ad")
FIG = Path(config["paths"]["figures_dir"])
REPORT = Path(config["paths"]["reports_dir"])

adata = sc.read_h5ad(IN)
adata.var_names_make_unique()
if "counts" not in adata.layers:
    adata.layers["counts"] = adata.X.copy()

sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata
sc.pp.highly_variable_genes(adata, n_top_genes=sp_cfg["n_top_genes"])

adata_hvg = adata[:, adata.var["highly_variable"]].copy()
sc.pp.scale(adata_hvg, max_value=10)
sc.tl.pca(adata_hvg)
sc.pp.neighbors(adata_hvg, n_neighbors=sp_cfg["n_neighbors"], n_pcs=sp_cfg["n_pcs"])
sc.tl.umap(adata_hvg, random_state=config["random_seed"])
sc.tl.leiden(adata_hvg, resolution=sp_cfg["leiden_resolution"], key_added="leiden_spatial")
adata.obsm["X_umap"] = adata_hvg.obsm["X_umap"]
adata.obs["leiden_spatial"] = adata_hvg.obs["leiden_spatial"].astype(str)

for program, genes in markers.items():
    present = [g for g in genes if g in adata.var_names]
    if len(present) >= 2:
        sc.tl.score_genes(adata, gene_list=present, score_name=f"score_{program}", use_raw=True)

sc.pl.umap(adata, color=["leiden_spatial"], show=False)
plt.savefig(FIG / "07_spatial_expression_umap.png", bbox_inches="tight")
plt.close()

colors = ["leiden_spatial", "score_hepatocyte", "score_cholangiocyte_biliary", "score_macrophage_kupffer_like", "score_stellate_fibroblast"]
colors = [c for c in colors if c in adata.obs.columns]

if "spatial" in adata.obsm.keys():
    try:
        sq.pl.spatial_scatter(adata, color=colors, show=False)
        plt.savefig(FIG / "08_spatial_marker_scores.png", bbox_inches="tight")
        plt.close()
    except Exception as e:
        print("Image-based spatial plot failed, using generic spatial embedding:", e)
        sc.pl.embedding(adata, basis="spatial", color=colors, show=False)
        plt.savefig(FIG / "08_spatial_marker_scores_generic.png", bbox_inches="tight")
        plt.close()

    sq.gr.spatial_neighbors(adata, coord_type="generic", library_key="spatial_sample")
    try:
        sq.gr.nhood_enrichment(adata, cluster_key="leiden_spatial", n_perms=sp_cfg["n_perms"])
        sq.pl.nhood_enrichment(adata, cluster_key="leiden_spatial", show=False)
        plt.savefig(FIG / "09_spatial_neighborhood_enrichment.png", bbox_inches="tight")
        plt.close()
    except Exception as e:
        print("Neighborhood enrichment failed:", e)

    try:
        genes = adata[:, adata.var["highly_variable"]].var_names[:sp_cfg["moran_n_genes"]]
        sq.gr.spatial_autocorr(adata, mode="moran", genes=genes, n_perms=sp_cfg["n_perms"], n_jobs=1)
        adata.uns["moranI"].to_csv(REPORT / "spatial_moranI_top_genes.csv")
    except Exception as e:
        print("Moran's I failed:", e)
else:
    print("No spatial coordinates found; spatial plots/statistics skipped.")

adata.write_h5ad(OUT, compression="gzip")
print(f"Wrote {OUT}")
