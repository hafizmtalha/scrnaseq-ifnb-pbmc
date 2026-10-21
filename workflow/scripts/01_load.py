"""Read the two count matrices from GEO and join them into one AnnData object.

The GEO files are in the classic 10x format: a sparse matrix, a list of
barcodes and a list of genes. The authors also give a metadata table with
the donor of each cell (found with demuxlet), their cell type labels and a
flag for droplets that hold more than one cell.
"""

import sys

import anndata as ad
import pandas as pd
import scanpy as sc

sys.stdout = sys.stderr = open(snakemake.log[0], "w")

raw = "data/raw"
genes = pd.read_csv(f"{raw}/genes.gz", sep="\t", header=None, names=["gene_id", "gene_name"])

parts = []
for condition in ["ctrl", "stim"]:
    # The matrix on GEO has genes as rows, so we flip it to cells as rows
    a = sc.read_mtx(f"{raw}/{condition}_matrix.gz").T
    barcodes = pd.read_csv(f"{raw}/{condition}_barcodes.gz", header=None)[0].values
    a.obs["barcode"] = barcodes
    a.obs_names = [f"{condition}_{b}" for b in barcodes]
    # Use the Ensembl ids for now because a few gene names appear twice
    a.var_names = genes["gene_id"].values
    a.obs["condition"] = condition
    parts.append(a)
    print(f"{condition}: {a.n_obs} cells")

adata = ad.concat(parts)
adata.var["gene_id"] = adata.var_names
adata.var_names = genes.set_index("gene_id").loc[adata.var_names, "gene_name"].values
adata.var_names_make_unique()

# The metadata table lists cells in the same order as the two barcode files.
# Some barcodes appear in both conditions, so the authors renamed a few of them
# (for example AAAC...-1 became AAAC...-11). We match rows by position and
# check that the barcodes agree before the dash.
meta = pd.read_csv(f"{raw}/metadata.gz", sep="\t", index_col=0)
assert len(meta) == adata.n_obs
same = meta.index.str.split("-").str[0] == adata.obs["barcode"].str.split("-").str[0]
assert same.all(), "metadata and barcodes are not in the same order"
assert (meta["stim"].values == adata.obs["condition"].values).all()

adata.obs["donor"] = meta["ind"].astype(str).values
adata.obs["published_label"] = meta["cell"].values
adata.obs["demuxlet"] = meta["multiplets"].values

# Sample label (donor plus condition)
adata.obs["sample"] = adata.obs["donor"] + "_" + adata.obs["condition"]

for col in ["condition", "donor", "sample", "published_label", "demuxlet"]:
    adata.obs[col] = adata.obs[col].astype("category")

print(adata)
print(adata.obs["demuxlet"].value_counts())
adata.write_h5ad(snakemake.output[0], compression="gzip")
