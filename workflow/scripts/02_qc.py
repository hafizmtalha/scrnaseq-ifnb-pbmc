"""Remove low quality cells and droplets that hold two cells.

We use three kinds of evidence:
  1. Simple counts per cell (number of genes found, share of mitochondrial reads)
  2. Demuxlet calls from the original study. Cells from two donors in one droplet
     show a mix of genotypes and are flagged as doublets.
  3. Scrublet, which finds doublets from expression alone. This also catches
     doublets made of two cells from the same donor, which demuxlet cannot see.
"""

import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc

sys.stdout = sys.stderr = open(snakemake.log[0], "w")
cfg = snakemake.config

adata = sc.read_h5ad(snakemake.input[0])
adata.var["mt"] = adata.var_names.str.startswith("MT-")
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], percent_top=None, log1p=False, inplace=True)


def qc_plot(a, path, title):
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    metrics = [("n_genes_by_counts", "Genes per cell"),
               ("total_counts", "UMIs per cell"),
               ("pct_counts_mt", "Mitochondrial reads (%)")]
    for ax, (col, label) in zip(axes, metrics):
        sc.pl.violin(a, col, groupby="condition", ax=ax, show=False, stripplot=False)
        ax.set_title(label)
        ax.set_xlabel("")
        ax.set_ylabel("")
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


qc_plot(adata, snakemake.output.before, f"Before filtering ({adata.n_obs:,} cells)")

steps = [("Starting cells", adata.n_obs)]

keep = (adata.obs["n_genes_by_counts"] >= cfg["min_genes"]) & \
       (adata.obs["n_genes_by_counts"] <= cfg["max_genes"]) & \
       (adata.obs["pct_counts_mt"] <= cfg["max_pct_mito"])
adata = adata[keep].copy()
steps.append(("After gene count and mitochondrial filters", adata.n_obs))

adata = adata[adata.obs["demuxlet"] == "singlet"].copy()
steps.append(("After removing demuxlet doublets and unclear cells", adata.n_obs))

# Run Scrublet on each sample on its own, since doublets only form within a run
sc.pp.scrublet(adata, batch_key="condition", random_state=cfg["seed"])
print(adata.obs["predicted_doublet"].value_counts())
adata = adata[~adata.obs["predicted_doublet"]].copy()
steps.append(("After removing Scrublet doublets", adata.n_obs))

sc.pp.filter_genes(adata, min_cells=cfg["min_cells_per_gene"])

qc_plot(adata, snakemake.output.after, f"After filtering ({adata.n_obs:,} cells)")

summary = pd.DataFrame(steps, columns=["step", "cells"])
summary.to_csv(snakemake.output.table, sep="\t", index=False)
print(summary)
print(f"Genes kept: {adata.n_vars}")

# Keep the raw counts safe. Later steps change adata.X but DE needs raw counts.
adata.layers["counts"] = adata.X.copy()
adata.write_h5ad(snakemake.output.adata, compression="gzip")
