from pathlib import Path
import squidpy as sq
import matplotlib.pyplot as plt

FIG = Path("figures")
REPORT = Path("reports")
FIG.mkdir(exist_ok=True)
REPORT.mkdir(exist_ok=True)

adata = sq.datasets.visium_hne_adata()
adata.obs["dataset_note"] = "Squidpy built-in Visium H&E smoke test; not liver biology"

import matplotlib.pyplot as plt

sq.pl.spatial_scatter(adata, color="cluster")
plt.savefig("figures/06_spatial_squidpy_smoketest.png", dpi=150, bbox_inches="tight")
plt.close()
plt.savefig(FIG / "spatial_smoketest_clusters.png", bbox_inches="tight")
plt.close()

sq.gr.spatial_neighbors(adata)
sq.gr.nhood_enrichment(adata, cluster_key="cluster", n_perms=100)
sq.pl.nhood_enrichment(adata, cluster_key="cluster", show=False)
plt.savefig(FIG / "spatial_smoketest_neighborhood_enrichment.png", bbox_inches="tight")
plt.close()

genes = adata.var_names[:500]
sq.gr.spatial_autocorr(adata, mode="moran", genes=genes, n_perms=100, n_jobs=1)
adata.uns["moranI"].to_csv(REPORT / "spatial_smoketest_moranI.csv")
adata.write_h5ad("data/processed/spatial_smoketest_squidpy.h5ad", compression="gzip")
print("Finished Squidpy smoke test")
