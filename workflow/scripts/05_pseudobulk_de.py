"""Differential expression between stimulated and control cells, per cell type.

Testing single cells as if they were independent samples gives far too many
false positives, because cells from the same donor are not independent.
Instead we add up the raw counts of all cells from one donor, one condition
and one cell type. This gives a "pseudobulk" sample that behaves like a
normal bulk RNA seq sample. We then run DESeq2 (through PyDESeq2) on these
samples with donor in the model, so each donor is compared with itself.
"""

import sys

import decoupler as dc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

sys.stdout = sys.stderr = open(snakemake.log[0], "w")
cfg = snakemake.config
PADJ, LFC = 0.05, 1.0

adata = sc.read_h5ad(snakemake.input[0])

pdata = dc.pp.pseudobulk(adata, sample_col="sample", groups_col="cell_type",
                         layer="counts", mode="sum")
dc.pp.filter_samples(pdata, min_cells=cfg["min_cells_pseudobulk"], min_counts=1000)
print(pdata)

# Overview of all pseudobulk samples
pview = pdata.copy()
sc.pp.normalize_total(pview, target_sum=1e6)
sc.pp.log1p(pview)
sc.pp.pca(pview, n_comps=10)
fig, axes = plt.subplots(1, 2, figsize=(13, 5))
sc.pl.pca(pview, color="cell_type", ax=axes[0], show=False, title="Pseudobulk samples by cell type", size=80)
sc.pl.pca(pview, color="condition", ax=axes[1], show=False, title="Pseudobulk samples by condition", size=80)
fig.tight_layout()
fig.savefig(snakemake.output.pca, dpi=150)
plt.close(fig)

results = []
for cell_type in pdata.obs["cell_type"].cat.categories:
    sub = pdata[pdata.obs["cell_type"] == cell_type].copy()
    n_ctrl = (sub.obs["condition"] == "ctrl").sum()
    n_stim = (sub.obs["condition"] == "stim").sum()
    if min(n_ctrl, n_stim) < 3:
        print(f"Skipping {cell_type}: only {n_ctrl} control and {n_stim} stimulated samples")
        continue

    # Keep genes with enough reads in enough samples
    dc.pp.filter_by_expr(sub, group="condition", min_count=10, min_total_count=15)
    sub.X = np.round(sub.X).astype(int)
    print(f"{cell_type}: {sub.n_obs} samples, {sub.n_vars} genes")

    dds = DeseqDataSet(adata=sub, design="~donor + condition", refit_cooks=True,
                       n_cpus=snakemake.threads, quiet=True)
    dds.deseq2()
    stats = DeseqStats(dds, contrast=["condition", "stim", "ctrl"], n_cpus=snakemake.threads, quiet=True)
    stats.summary()
    res = stats.results_df.reset_index(names="gene")
    res.insert(0, "cell_type", cell_type)
    results.append(res)

de = pd.concat(results, ignore_index=True)
de = de.sort_values(["cell_type", "padj"])
de.to_csv(snakemake.output.all_results, sep="\t", index=False)

sig = de[(de["padj"] < PADJ) & (de["log2FoldChange"].abs() > LFC)]
summary = (sig.assign(direction=np.where(sig["log2FoldChange"] > 0, "up", "down"))
              .groupby(["cell_type", "direction"]).size().unstack(fill_value=0)
              .reindex(columns=["up", "down"], fill_value=0))
summary["total"] = summary.sum(axis=1)
summary = summary.sort_values("total", ascending=False)
summary.to_csv(snakemake.output.summary, sep="\t")
print(summary)

# Bar chart of DE gene counts
fig, ax = plt.subplots(figsize=(8, 5))
summary[["up", "down"]].plot.barh(stacked=True, ax=ax, color=["#C44E52", "#4C72B0"])
ax.invert_yaxis()
ax.set_xlabel(f"Number of genes (adjusted p < {PADJ}, |log2 fold change| > {LFC:g})")
ax.set_ylabel("")
ax.set_title("Genes changed by interferon beta in each cell type")
ax.legend(["Up in stim", "Down in stim"])
fig.tight_layout()
fig.savefig(snakemake.output.bar, dpi=150)
plt.close(fig)

# Volcano plot for CD14 monocytes, the cell type with the strongest response
mono = de[de["cell_type"] == "CD14 monocytes"].dropna(subset=["padj"]).copy()
mono["neglog10"] = -np.log10(mono["padj"].clip(lower=1e-300))
colour = np.where((mono["padj"] < PADJ) & (mono["log2FoldChange"] > LFC), "#C44E52",
                  np.where((mono["padj"] < PADJ) & (mono["log2FoldChange"] < -LFC), "#4C72B0", "lightgrey"))
fig, ax = plt.subplots(figsize=(8, 6))
ax.scatter(mono["log2FoldChange"], mono["neglog10"], c=colour, s=6, alpha=0.7, linewidths=0)
for _, r in mono.nsmallest(15, "padj").iterrows():
    ax.annotate(r["gene"], (r["log2FoldChange"], r["neglog10"]), fontsize=8,
                xytext=(3, 3), textcoords="offset points")
ax.axvline(-LFC, ls="--", c="grey", lw=0.8)
ax.axvline(LFC, ls="--", c="grey", lw=0.8)
ax.axhline(-np.log10(PADJ), ls="--", c="grey", lw=0.8)
ax.set_xlabel("log2 fold change (stim vs ctrl)")
ax.set_ylabel("-log10 adjusted p value")
ax.set_title("CD14 monocytes: pseudobulk DESeq2")
fig.tight_layout()
fig.savefig(snakemake.output.volcano, dpi=150)
plt.close(fig)

# Heatmap: fold changes of the top genes from each cell type, across all cell types
top_genes = (de.dropna(subset=["padj"]).sort_values("padj").groupby("cell_type").head(6)["gene"].unique())
lfc = de[de["gene"].isin(top_genes)].pivot(index="gene", columns="cell_type", values="log2FoldChange")
lfc = lfc.fillna(0)
g = sns.clustermap(lfc, cmap="RdBu_r", center=0, vmin=-6, vmax=6, figsize=(9, 0.25 * len(lfc) + 3),
                   cbar_kws={"label": "log2 fold change"}, col_cluster=True, row_cluster=True)
g.ax_heatmap.set_xlabel("")
g.ax_heatmap.set_ylabel("")
g.figure.suptitle("Top genes from each cell type (0 means the gene was not tested there)", y=1.01)
g.savefig(snakemake.output.heatmap, dpi=150, bbox_inches="tight")
plt.close("all")
