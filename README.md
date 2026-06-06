# Liver atlas scVI + Squidpy demo

A reproducible Python/Podman workflow for public liver single-cell and spatial transcriptomics data using Scanpy, scVI-tools, and Squidpy.

## Why this project exists

This project was built as an interview-preparation portfolio item for computational liver biology. It demonstrates that I can build a reproducible workflow, learn a new tissue system, reason about liver cell types, connect single-cell and spatial transcriptomics analysis, and communicate limitations honestly.

## Biological question

Can public liver atlas data recover broad liver cell compartments and spatially organized programs in a way that provides a useful starting point for liver disease and liver cancer research?

## Workflow

1. Download public Liver Cell Atlas resources.
2. Import single-cell data into AnnData.
3. Perform single-cell QC and Scanpy baseline analysis.
4. Score broad liver marker programs.
5. Train an scVI model on the main raw-count single-cell dataset.
6. Use the scVI latent representation for UMAP and clustering.
7. Import Visium-like spatial expression data.
8. Attach GEO-derived Visium tissue coordinates.
9. Use Squidpy for generic spatial marker plotting, spatial neighbors, neighborhood enrichment, and Moran's I.
10. Run an independent processed HCC holdout overview using the same marker programs without forcing scVI on non-count data.


## Quick start

Run these from the host. The Makefile executes Python inside the container, so you do not need host pandas/Scanpy/scVI for these commands.

```bash
make build
make gpu-test
make download
make inspect
make import-sc
make qc
make markers
make scvi
```

For spatial analysis:

```bash
make import-visium
make add-visium-coords
make spatial
```

For the processed HCC holdout, first prepare:

```bash
mkdir -p data/raw/holdout
gunzip -c data/raw/GSE156625/GSE156625_HCCscanpyobj.h5ad.gz > data/raw/holdout/holdout_liver.h5ad
```

Then run:

```bash
make holdout
```

For Squidpy installation/code-path checks only:

```bash
make spatial-smoketesttest
```

The smoke test proves the Squidpy code path but is not used for liver biological interpretation.

## Current outputs

Representative committed figures include:

```text
figures/01_qc_violin_before_filtering.png
figures/02_scanpy_umap_overview.png
figures/03_umap_marker_scores.png
figures/04_marker_dotplot_by_leiden.png
figures/05_scvi_umap_overview.png
figures/06_scvi_umap_marker_scores.png
figures/07_spatial_expression_umap.png
figures/08_spatial_marker_scores_generic.png
figures/09_spatial_neighborhood_enrichment.png
figures/10_holdout_marker_overview.png
```

## Caveats

- Marker scores are broad sanity checks, not definitive expert annotation.
- Visium spots are not single cells.
- The spatial workflow currently uses generic coordinates; image overlay may fall back because `adata.uns["spatial"]` is not fully built.
- scVI latent space is useful but not directly interpretable.
- scVI should only be trained on an appropriate raw-count matrix.
- The GSE156625 HCC holdout object is preprocessed and lacks a count layer, so the holdout workflow reuses the existing embedding and marker programs rather than forcing scVI.
- Raw data and large processed objects are not committed to Git.
