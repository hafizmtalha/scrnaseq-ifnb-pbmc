"""Cluster the cells and give each cluster a cell type name.

We name clusters from well known PBMC marker genes. For each cell we score
every marker set, then give each cluster the cell type with the highest
average score. At the end we compare our labels with the labels the
original authors published, as an independent check.
"""

import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import seaborn as sns

sys.stdout = sys.stderr = open(snakemake.log[0], "w")
cfg = snakemake.config

adata = sc.read_h5ad(snakemake.input[0])

sc.tl.leiden(adata, resolution=cfg["leiden_resolution"], flavor="igraph",
             n_iterations=2, random_state=cfg["seed"], key_added="leiden")
print(adata.obs["leiden"].value_counts().sort_index())

MARKERS = {
    "CD4 T cells": ["IL7R", "CCR7", "LTB", "CD3E"],
    "CD8 T cells": ["CD8A", "CD8B", "CD3E"],
    "NK cells": ["GNLY", "NKG7", "PRF1", "KLRD1"],
    "B cells": ["MS4A1", "CD79A", "CD79B", "CD19"],
    "CD14 monocytes": ["CD14", "LYZ", "S100A8", "S100A9"],
    "FCGR3A monocytes": ["FCGR3A", "MS4A7", "LST1"],
    "Dendritic cells": ["FCER1A", "CLEC10A", "CD1C"],
    "Plasmacytoid DCs": ["LILRA4", "IL3RA", "CLEC4C"],
    "Megakaryocytes": ["PPBP", "PF4", "GNG11"],
}
MARKERS = {k: [g for g in v if g in adata.var_names] for k, v in MARKERS.items()}

for cell_type, genes in MARKERS.items():
    sc.tl.score_genes(adata, genes, score_name=f"score_{cell_type}", random_state=cfg["seed"])

score_cols = [f"score_{c}" for c in MARKERS]
cluster_scores = adata.obs.groupby("leiden", observed=True)[score_cols].mean()
cluster_scores.columns = list(MARKERS)
# Highly expressed marker sets (like the monocyte one) give big raw scores
# everywhere. Scaling each cell type column first means a cluster gets the
# label whose markers stand out most in that cluster compared with the others.
z = (cluster_scores - cluster_scores.mean()) / cluster_scores.std()
cluster_to_type = z.idxmax(axis=1)
print(cluster_to_type)

adata.obs["cell_type"] = adata.obs["leiden"].map(cluster_to_type).astype("category")

# Scores per cluster, so a reader can see why each name was picked
fig, ax = plt.subplots(figsize=(9, 0.4 * len(z) + 2))
sns.heatmap(z, cmap="RdBu_r", center=0, ax=ax, cbar_kws={"label": "Scaled marker score"})
ax.set_yticklabels([f"{c}  ({cluster_to_type[c]})" for c in z.index], rotation=0)
ax.set_ylabel("Leiden cluster (assigned type)")
ax.set_title("Marker scores for each cluster")
fig.tight_layout()
fig.savefig(snakemake.output.scores, dpi=150)
plt.close(fig)

# UMAPs
for key, path, title in [("leiden", snakemake.output.clusters, "Leiden clusters"),
                         ("cell_type", snakemake.output.umap, "Cell types")]:
    fig, ax = plt.subplots(figsize=(8, 6))
    sc.pl.umap(adata, color=key, ax=ax, show=False, frameon=False, title=title,
               legend_loc="on data" if key == "leiden" else "right margin", legend_fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)

# Dot plot of the marker genes by cell type
dp = sc.pl.dotplot(adata, MARKERS, groupby="cell_type", standard_scale="var",
                   show=False, return_fig=True)
dp.savefig(snakemake.output.dotplot, dpi=150, bbox_inches="tight")
plt.close("all")

# Top marker genes per cluster found from the data itself
sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")
markers = sc.get.rank_genes_groups_df(adata, group=None)
markers = markers[markers["pvals_adj"] < 0.05].groupby("group", observed=True).head(25)
markers["cell_type"] = markers["group"].map(cluster_to_type)
markers.to_csv(snakemake.output.markers, sep="\t", index=False)

# How well do our labels agree with the published ones?
table = pd.crosstab(adata.obs["cell_type"], adata.obs["published_label"], normalize="index")
fig, ax = plt.subplots(figsize=(9, 6))
sns.heatmap(table, cmap="Blues", annot=True, fmt=".2f", ax=ax, cbar=False)
ax.set_xlabel("Published label (Kang et al. 2018)")
ax.set_ylabel("Our label")
ax.set_title("Share of each of our cell types given each published label")
fig.tight_layout()
fig.savefig(snakemake.output.compare, dpi=150)
plt.close(fig)

print(table.round(2))

counts = pd.crosstab([adata.obs["cell_type"]], adata.obs["condition"])
counts.to_csv(snakemake.output.counts, sep="\t")
print(counts)

adata.write_h5ad(snakemake.output.adata, compression="gzip")
