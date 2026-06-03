from pathlib import Path
import scanpy as sc
import matplotlib.pyplot as plt
from utils import load_yaml, ensure_dirs

config = load_yaml("config/project_config.yaml")
markers = load_yaml("config/marker_genes.yaml")
ensure_dirs(config)

IN = Path("data/processed/lca_human_qc_scanpy.h5ad")
OUT = Path("data/processed/lca_human_marker_scored.h5ad")
FIG = Path(config["paths"]["figures_dir"])
TAB = Path(config["paths"]["reports_dir"])

adata = sc.read_h5ad(IN)
score_cols = []

for program, genes in markers.items():
    present = [g for g in genes if g in adata.var_names]
    missing = [g for g in genes if g not in adata.var_names]
    print(f"{program}: {len(present)} present, {len(missing)} missing")
    if len(present) >= 2:
        col = f"score_{program}"
        sc.tl.score_genes(adata, gene_list=present, score_name=col, use_raw=True)
        score_cols.append(col)

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
    plt.savefig(FIG / "03_umap_marker_scores.png", bbox_inches="tight")
    plt.close()

marker_panel = []
for genes in markers.values():
    marker_panel.extend(genes)
marker_panel = list(dict.fromkeys([g for g in marker_panel if g in adata.var_names]))[:80]

if marker_panel:
    sc.pl.dotplot(adata, marker_panel, groupby="leiden_scanpy", standard_scale="var", show=False)
    plt.savefig(FIG / "04_marker_dotplot_by_leiden.png", bbox_inches="tight")
    plt.close()

if score_cols:
    summary = adata.obs.groupby("leiden_scanpy")[score_cols].mean()
    summary.to_csv(TAB / "marker_score_means_by_leiden.csv")
    print(summary)

adata.write_h5ad(OUT, compression="gzip")
print(f"Wrote {OUT}")
