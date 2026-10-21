"""Turn gene level results into pathway level results.

We use two public gene set collections through decoupler:
  * PROGENy, which gives the genes that respond when 14 signalling pathways
    are switched on. This tells us which signals the cells are receiving.
  * MSigDB hallmark gene sets, 50 well defined biological processes.

For each cell type we take the DESeq2 test statistic of every gene and ask,
with a simple linear model (ulm), whether the genes in a set lean up or down
more than chance.
"""

import sys

import decoupler as dc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import seaborn as sns

sys.stdout = sys.stderr = open(snakemake.log[0], "w")

de = pd.read_csv(snakemake.input.de, sep="\t")
stat = de.pivot(index="cell_type", columns="gene", values="stat").fillna(0)

progeny = dc.op.progeny(organism="human", top=500)
hallmark = dc.op.hallmark(organism="human")

tables = []
scores = {}
for name, net in [("PROGENy", progeny), ("Hallmark", hallmark)]:
    es, pv = dc.mt.ulm(data=stat, net=net)
    scores[name] = es
    long = es.reset_index(names="cell_type").melt(id_vars="cell_type", var_name="pathway", value_name="score")
    pvl = pv.reset_index(names="cell_type").melt(id_vars="cell_type", var_name="pathway", value_name="pvalue")
    long = long.merge(pvl, on=["cell_type", "pathway"])
    long.insert(0, "collection", name)
    tables.append(long)

pd.concat(tables).to_csv(snakemake.output.table, sep="\t", index=False)

# Keep the hallmark sets that change the most across cell types
hm = scores["Hallmark"]
hm = hm[hm.abs().max().sort_values(ascending=False).index[:15]]
hm.columns = [c.replace("HALLMARK_", "").replace("_", " ").title() for c in hm.columns]

fig, axes = plt.subplots(1, 2, figsize=(17, 6), gridspec_kw={"width_ratios": [14, 15]})
for ax, (title, data) in zip(axes, [("PROGENy pathway activity", scores["PROGENy"]),
                                     ("Top 15 hallmark gene sets", hm)]):
    sns.heatmap(data, cmap="RdBu_r", center=0, ax=ax, cbar_kws={"label": "Activity score (stim vs ctrl)"})
    ax.set_title(title)
    ax.set_ylabel("")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=60, ha="right")
fig.tight_layout()
fig.savefig(snakemake.output.heatmap, dpi=150)
plt.close(fig)

# Score each single cell for the interferon alpha response genes
adata = sc.read_h5ad(snakemake.input.adata)
ifn_genes = hallmark.loc[hallmark["source"] == "INTERFERON_ALPHA_RESPONSE", "target"]
ifn_genes = [g for g in ifn_genes if g in adata.var_names]
sc.tl.score_genes(adata, ifn_genes, score_name="interferon_score", random_state=0)

fig, axes = plt.subplots(1, 2, figsize=(12, 5))
vmax = adata.obs["interferon_score"].quantile(0.99)
for ax, cond in zip(axes, ["ctrl", "stim"]):
    sub = adata[adata.obs["condition"] == cond]
    sc.pl.umap(sub, color="interferon_score", ax=ax, show=False, frameon=False,
               vmin=0, vmax=vmax, cmap="viridis", title=f"Interferon alpha response score: {cond}")
fig.tight_layout()
fig.savefig(snakemake.output.umap, dpi=150)
plt.close(fig)
