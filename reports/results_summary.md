# Results summary

## Project aim

This project demonstrates a reproducible Python/Podman workflow for public liver single-cell and spatial transcriptomics data using Scanpy, scVI-tools, and Squidpy.

## Single-cell baseline

The baseline Scanpy workflow imports public liver single-cell data, computes QC metrics, filters low-quality cells, normalizes expression, identifies highly variable genes, runs PCA, builds a neighbor graph, computes UMAP, and clusters cells with Leiden.

## Marker-based biological sanity checks

Broad marker programs were used to identify likely hepatocyte, cholangiocyte, endothelial, macrophage/Kupffer-like, stellate/fibroblast, T-cell, NK-cell, B-cell, and plasma-cell regions. These scores are sanity checks, not definitive annotations.

## scVI analysis

The scVI workflow trains a count-based latent variable model on the raw count layer. The learned latent representation is used for neighbors, UMAP, and Leiden clustering. Results are interpreted with marker programs rather than treating latent space as directly biological.

## Spatial analysis

The Squidpy workflow imports Visium-like spatial data when available, computes spot-level marker scores, visualizes programs on tissue coordinates, builds a spatial neighbor graph, estimates neighborhood enrichment, and calculates spatial autocorrelation for highly variable genes.

## Holdout analysis

The independent GSE156625 HCC holdout can be placed at `data/raw/holdout/holdout_liver.h5ad`. The validated object is already processed: it contains existing PCA/UMAP embeddings, louvain clusters, and patient/tumor metadata, but no count layer. The holdout workflow therefore reuses the existing embedding and applies the same marker programs as simple mean-expression scores instead of forcing raw-count preprocessing or scVI training on transformed values. This demonstrates workflow reuse and honest limitation handling rather than one-dataset tuning.

## Limitations

- Public data only.
- No raw FASTQ alignment in this sprint.
- Marker scoring is a broad sanity check, not expert annotation.
- Visium spots are not single cells.
- Image-overlay spatial plotting is not complete until `adata.uns["spatial"]` is built properly.
- scVI latent dimensions are not directly interpretable and should only be trained on an appropriate count matrix.
- The first GSE156625 HCC holdout is preprocessed and lacks a raw count layer, so it is used for marker-program reuse rather than scVI.
- Final biological conclusions would require expert review and experimental validation.
