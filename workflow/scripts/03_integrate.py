"""Normalize, find variable genes, reduce dimensions and integrate.

Interferon beta changes the expression of hundreds of genes in every cell type.
Without correction, treated and control cells of the same type sit apart on
the UMAP. Harmony adjusts the PCA so that cells group by type and not by
treatment. We only use the corrected PCA for clustering and the UMAP. All
differential expression later on is done on the raw counts.
"""

import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import harmonypy
import scanpy as sc

sys.stdout = sys.stderr = open(snakemake.log[0], "w")
cfg = snakemake.config
seed = cfg["seed"]

adata = sc.read_h5ad(snakemake.input[0])

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Pick variable genes in each condition and keep genes that show up in both
sc.pp.highly_variable_genes(adata, n_top_genes=cfg["n_top_genes"], batch_key="condition", flavor="seurat")
print(f"Variable genes: {adata.var['highly_variable'].sum()}")

sc.pp.pca(adata, n_comps=50, mask_var="highly_variable", random_state=seed)

# UMAP before integration, kept only for the comparison figure
sc.pp.neighbors(adata, n_pcs=cfg["n_pcs"], random_state=seed)
sc.tl.umap(adata, random_state=seed)
adata.obsm["X_umap_uncorrected"] = adata.obsm["X_umap"].copy()

harmony = harmonypy.run_harmony(adata.obsm["X_pca"], adata.obs, "condition", random_state=seed)
corrected = np.asarray(harmony.Z_corr)
# Older harmonypy versions return dimensions by cells, newer ones cells by dimensions
if corrected.shape[0] != adata.n_obs:
    corrected = corrected.T
adata.obsm["X_pca_harmony"] = corrected
sc.pp.neighbors(adata, use_rep="X_pca_harmony", n_pcs=cfg["n_pcs"], random_state=seed)
sc.tl.umap(adata, random_state=seed)

fig, axes = plt.subplots(2, 2, figsize=(12, 10))
for row, (basis, name) in enumerate([("X_umap_uncorrected", "Before Harmony"), ("X_umap", "After Harmony")]):
    sc.pl.embedding(adata, basis=basis, color="condition", ax=axes[row, 0], show=False,
                    title=f"{name}: condition", frameon=False)
    sc.pl.embedding(adata, basis=basis, color="published_label", ax=axes[row, 1], show=False,
                    title=f"{name}: published cell type", frameon=False)
fig.tight_layout()
fig.savefig(snakemake.output.umap, dpi=150)
plt.close(fig)

adata.write_h5ad(snakemake.output.adata, compression="gzip")
