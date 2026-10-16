#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os as os
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
import anndata as ad
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import seaborn as sns
import statsmodels.stats.multitest as smm
import statsmodels.api as sm
import statsmodels.formula.api as smf
import adjustText
from scipy.stats import pearsonr, spearmanr
from scipy.stats import chi2_contingency
import bipca
from bipca import plotting
from bipca.utils import stabilize_matrix
from bipca.experiments import new_svd as SVD
from openTSNE import TSNE
from multiprocessing import Pool
import matplotlib.pyplot as plt
from tqdm.contrib.concurrent import thread_map,process_map
from functools import partial
import torch
from threadpoolctl import threadpool_limits


# In[2]:


celltypes = ['04_Crnde_Lmx1a', '27_Lmx1a_Ttc6', '32_Gucy2c_Ttc6', "36_Ebf2_Crhbp"]
adata = sc.read("LMX_bipca.h5ad")
adata


# In[2]:


adata = sc.read("allgenes_bipca_hypomap.h5ad")
adata


# In[4]:


ax = sc.pl.embedding(
    adata,
    basis="tsne",
    color="detailed_celltype_annotation",
    legend_loc="on data",
    legend_fontsize="xx-small",
    return_fig=True)
ax=plt.gca()
texts = ax.texts
adjustText.adjust_text(texts, 
                       ax=ax,
                       arrowprops=dict(arrowstyle="-", color='black', lw=0.5))
plt.show()


# In[5]:


marker_genes = ["Lmx1a", "Lmx1b", "Irx3", "Irx5"]
sc.pl.embedding(
    adata,
    basis="tsne",
    color=marker_genes,
    layer="scran_normalization", 
    color_map="YlOrBr", 
    save="_lmx_irx.pdf", 
    vmax=1.5)


# In[4]:


adata.obs['plot'] = adata.obs['detailed_celltype_annotation'].where(
    adata.obs['detailed_celltype_annotation'].isin(celltypes)
)

sc.pl.embedding(
    adata,
    basis="tsne",
    color="plot",
    legend_loc="on data",
    legend_fontsize="xx-small",
    save="_lmx_subset_large.pdf"
)


# In[11]:


gpcr_genes = [
    "Slc17a6", #VGlut2
    "Slc17a7", #VGlut1 LHA https://www.cell.com/neuron/fulltext/S0896-6273%2821%2900693-0
    "Gad1",
    "Gad2",
    "Slc32a1", #Vgat
    "Hcrt", #Ox, should be Vglut2 positive
    "Pmch", #MCH should be Gad1 positive
    "Nts", #LH
    "Gal",
    "Crh",
    "Pdyn", 
    "Tac1", 
    "Tac2", 
    "Pnoc", #LH
    "Hdc", #Histidine decarboxylase, Tuberomammillary
    "Slc18a2", #Vmat2, Tuberomammillary
    "Lepr", #Premammillary reproductive
    "Gucy2c",
    "Slc6a3", #Dopamine
    "Th",
    "Ddc"
    ]
sc.pl.embedding(adata,basis="tsne",
                color=gpcr_genes, 
                color_map="YlOrBr", 
                vmax=1.5,
                save="_lmx_classification_markers.pdf")


# In[9]:


gpcr_genes = [
    "Lepr",
    "Insr",
    "Ghsr", #Ghrelin
    "Glp1r",
    "Gipr",
    "Mc4r",
    "Mc3r",
    "Npy1r",
    "Npy2r",
    "Npy5r",
    "Mchr1",
    "Hcrtr1", #Orexin
    "Hcrtr2", 
    "Htr2c", #Serotonin
    "Cnr1", #Cannabinoid
    "Oprm1", #Opioid µ
    "Hrh1", #Histamine
    "Hrh2",
    "Hrh3",
    "Calcr",
    "Ramp1", 
    "Ramp3"
    ]
sc.pl.embedding(adata,basis="tsne",
                color=gpcr_genes, 
                color_map="YlOrBr", 
                vmax=1.5,
                save="_lmx_GPCRs.pdf")


# In[5]:


gpcr_genes = [
    "Adcyap1", "Nos1", "Bdnf", "Esr1", "Lepr", "Nr5a1", "Tac1", "Rprm", "Pdyn", "Sst", "Hpcal1", "Gal"]
sc.pl.embedding(adata,basis="tsne",
                color=gpcr_genes, 
                color_map="YlOrBr", 
                vmax=1.5,
                save="_vmh_classification_markers.pdf")


# In[6]:


celltypes = ['05_Qrfpr_Pdzrn3', '00_Hmcn2_Gpr149', '03_Rai14_Acvr1c', "02_Fezf1_Lama1", "09_Nfia_Ccbe1"]
adata.obs['plot'] = adata.obs['detailed_celltype_annotation'].where(
    adata.obs['detailed_celltype_annotation'].isin(celltypes)
)

sc.pl.embedding(
    adata,
    basis="tsne",
    color="plot",
    legend_loc="on data",
    legend_fontsize="xx-small",
    save="_vmh_subset_large.pdf"
)


# In[5]:


adata = adata[adata.obs['plot'].notna()].copy()
adata


# In[6]:


print(adata.X[:5, :5].toarray())


# In[7]:


gene= "Glp1r"
if gene in adata.var_names:
    print(f"✔️ Found {gene!r} in adata.var_names (n_vars = {adata.n_vars})")
else:
    print(f"❌ {gene!r} NOT in adata.var_names (n_vars = {adata.n_vars})")


# In[8]:


X_sub,index = stabilize_matrix(adata.X,threshold=100)
X_sub


# In[9]:


adata_sub = sc.AnnData(X=X_sub,obs=adata.obs.iloc[index[0],:],var = adata.var.iloc[index[1],:])
adata_sub.X


# In[10]:


# gene-name indices aligned to the two matrices
genes1 = pd.Index(adata_sub.var_names)   # length should be 15916
genes2 = pd.Index(adata.var_names)   # length should be 22180

# set differences
only_in_1 = genes1.difference(genes2)   # present in matrix 1, absent in matrix 2
only_in_2 = genes2.difference(genes1)   # present in matrix 2, absent in matrix 1

# save to CSV
only_in_1.to_series().to_csv("genes_only_in_matrix1.csv", index=False)
only_in_2.to_series().to_csv("genes_only_in_matrix2.csv", index=False)

# (optional) common genes, if you need them
common = genes1.intersection(genes2)
common.to_series().to_csv("genes_common.csv", index=False)

print(f"Only in matrix1: {len(only_in_1)}")
print(f"Only in matrix2: {len(only_in_2)}")
print(f"Common: {len(common)}")


# In[11]:


X0 = adata.X
n_genes = X0.shape[1]
kept_mask = np.zeros(n_genes, dtype=bool)
kept_mask[index[1]] = True
excluded_mask = ~kept_mask

def col_vector_to_1d(arr):
    # turn (1,n) sparse or ndarray into 1D
    if sparse.issparse(arr):
        return np.asarray(arr).ravel()
    return np.asarray(arr).ravel()

is_sparse = sparse.issparse(X0)

# total counts per gene
gene_sum = col_vector_to_1d(X0.sum(axis=0))
# number of cells with >0 counts per gene
if is_sparse:
    n_cells_expr = col_vector_to_1d((X0 > 0).sum(axis=0))
else:
    n_cells_expr = (X0 > 0).sum(axis=0)

# per-gene max count (use efficient sparse max)
def max_per_gene(X):
    if sparse.issparse(X):
        # work in CSC so columns are contiguous in memory
        Xc = X.tocsc()
        n = Xc.shape[1]
        out = np.zeros(n, dtype=Xc.data.dtype if Xc.data.size else float)
        indptr = Xc.indptr
        data = Xc.data
        for j in range(n):
            if indptr[j] < indptr[j+1]:
                out[j] = data[indptr[j]:indptr[j+1]].max()
            else:
                out[j] = 0.0  # empty column -> max is 0
        return out
    else:
        # dense
        return np.asarray(X).max(axis=0).ravel()

gene_max = max_per_gene(X0)
print("gene_max shape:", gene_max.shape)  # should be (n_genes,)

# per-gene mean and variance
if is_sparse:
    # means efficiently
    gene_mean = gene_sum / X0.shape[0]
    # variance: compute E[x^2] - (E[x])^2
    x2_sum = col_vector_to_1d(X0.power(2).sum(axis=0))
    gene_var = (x2_sum / X0.shape[0]) - gene_mean**2
    gene_var = np.clip(gene_var, 0, None)
else:
    gene_mean = X0.mean(axis=0)
    gene_var  = X0.var(axis=0)


# 1) Coerce everything to flat 1D arrays of length n_genes
def to_1d(a):
    a = np.asarray(a)
    return a.reshape(-1)  # handles (n,), (1,n), (n,1)

gene_sum      = to_1d(gene_sum)
gene_mean     = to_1d(gene_mean)
gene_var      = to_1d(gene_var)
gene_max      = to_1d(gene_max)
n_cells_expr  = to_1d(n_cells_expr)
kept_mask     = to_1d(kept_mask).astype(bool)
excluded_mask = to_1d(excluded_mask).astype(bool)

# 2) Use the gene names that align with X0 (NOT the filtered adata)
# Replace `adata_before` with whatever holds your original var_names for X0
gene_names = np.array(adata.var_names)          # length must equal n_genes
assert len(gene_names) == n_genes, "gene_names length doesn't match X0.shape[1]"

# 3) (Optional) sanity check all lengths before building the DataFrame
for name, arr in [
    ("gene_sum", gene_sum), ("gene_mean", gene_mean), ("gene_var", gene_var),
    ("gene_max", gene_max), ("n_cells_expr", n_cells_expr),
    ("kept_mask", kept_mask), ("excluded_mask", excluded_mask),
]:
    assert arr.shape[0] == n_genes, f"{name} has length {arr.shape[0]} != {n_genes}"

# 4) Build the metrics table with the gene names as the index
metrics = pd.DataFrame(
    {
        "kept": kept_mask,
        "excluded": excluded_mask,
        "sum": gene_sum,
        "mean": gene_mean,
        "var": gene_var,
        "max": gene_max,
        "n_cells_expr": n_cells_expr,
    },
    index=gene_names
)

# Save if you like
metrics.to_csv("stabilize_matrix_gene_metrics.csv")
print(metrics.shape)


# In[12]:


# (A) Histogram of per-gene MAX with a vertical line at your threshold=100
plt.figure(figsize=(8,5))
plt.hist(np.log1p(metrics.loc[kept_mask, "max"]), bins=20, alpha=0.6, label="Kept")
plt.hist(np.log1p(metrics.loc[excluded_mask, "max"]), bins=20, alpha=0.6, label="Excluded")
plt.xlabel("log1p(Max count per gene)")
plt.ylabel("Number of genes")
plt.title("Distribution of per-gene MAX (Kept vs Excluded)")
plt.legend()
plt.tight_layout()
plt.show()

# (B) Mean–Variance plot (log–log) highlighting excluded genes
plt.figure(figsize=(6,6))
plt.scatter(np.log1p(metrics.loc[kept_mask, "mean"]),
            np.log1p(metrics.loc[kept_mask, "var"]),
            s=8, alpha=0.3, label="Kept")
plt.scatter(np.log1p(metrics.loc[excluded_mask, "mean"]),
            np.log1p(metrics.loc[excluded_mask, "var"]),
            s=10, alpha=0.7, label="Excluded")
plt.xlabel("log1p(mean)")
plt.ylabel("log1p(variance)")
plt.title("Mean–Variance (Kept vs Excluded)")
plt.legend()
plt.tight_layout()
plt.show()

# (C) Expression prevalence vs mean (how many cells express each gene)
plt.figure(figsize=(8,5))
plt.scatter(metrics.loc[kept_mask, "n_cells_expr"], metrics.loc[kept_mask, "mean"],
            s=8, alpha=0.3, label="Kept")
plt.scatter(metrics.loc[excluded_mask, "n_cells_expr"], metrics.loc[excluded_mask, "mean"],
            s=10, alpha=0.7, label="Excluded")
plt.xlabel("# cells with >0 expression")
plt.ylabel("Mean expression")
plt.title("Prevalence vs Mean (Kept vs Excluded)")
plt.legend()
plt.tight_layout()
plt.show()

# --- 5) If you want a CSV of excluded genes only (with metrics) ---
metrics.loc[metrics["excluded"]].to_csv("excluded_genes_with_metrics.csv")
print(f"Excluded genes: {metrics['excluded'].sum()} / {len(metrics)}")


# In[13]:


adata = sc.AnnData(X=X_sub,obs=adata.obs.iloc[index[0],:],var = adata.var.iloc[index[1],:])


# In[14]:


n_threads = 8
torch.set_num_threads(n_threads)
with threadpool_limits(limits=n_threads):
    op = bipca.BiPCA(n_components=-1,seed=42) # get the BiPCA operator, here n_components=-1 is to do full SVD
    Z = op.fit_transform(X_sub.toarray()) # take the count data as input and run bipca
    op.get_plotting_spectrum() # get the fitting parameters


# In[15]:


op.write_to_adata(adata) # save the results back to anndata


# In[16]:


plotting.MP_histograms_from_bipca(adata,median=False,bins=100) 
plt.savefig('bipca_LMX_filtered.pdf')


# In[17]:


adata.obsm['pcs'] = SVD(adata.layers['Z_biwhite'],adata.uns["bipca"]['rank'])
sc.pp.neighbors(adata, use_rep = "pcs",method="gauss", n_pcs=adata.uns["bipca"]['rank'])
sc.tl.leiden(adata, resolution=1, random_state=42, flavor="leidenalg", n_iterations=-1)
adata.obsm['tsne'] = np.array(TSNE(n_jobs=32).fit(adata.obsm['pcs']))


# In[18]:


sc.pl.embedding(adata,basis="tsne",color="leiden", legend_loc="on data", legend_fontsize="xx-small", save="_leiden_LMX.pdf")


# In[19]:


adata.write("raw_LMX_bipca.h5ad")
adata


# In[21]:


adata.obs['n_counts'] = adata.X.sum(axis=1)
adata.obs['n_genes'] = (adata.X > 0).sum(axis=1)
scales_counts = sc.pp.normalize_total(adata, target_sum=None, inplace=False)
# log1p transform
adata.layers["log1p_norm"] = sc.pp.log1p(scales_counts["X"], copy=True)
fig, axes = plt.subplots(1, 2, figsize=(10, 5))
p1 = sns.histplot(adata.obs["n_counts"], bins=100, kde=False, ax=axes[0])
axes[0].set_title("Total counts")
p2 = sns.histplot(adata.layers["log1p_norm"].sum(1), bins=100, kde=False, ax=axes[1])
axes[1].set_title("Shifted logarithm")
plt.show()


# In[22]:


os.environ['R_HOME'] = '/.conda/envs/pca_env/lib/R'
import logging
import anndata2ri
import rpy2.rinterface_lib.callbacks as rcb
import rpy2.robjects as ro
rcb.logger.setLevel(logging.ERROR)
ro.pandas2ri.activate()
anndata2ri.activate()

get_ipython().run_line_magic('load_ext', 'rpy2.ipython')


# In[23]:


get_ipython().run_cell_magic('R', '', 'library(scran)\nlibrary(BiocParallel)\n')


# In[24]:


adata.layers["counts"] = adata.X
adata


# In[25]:


from scipy.sparse import issparse
data_mat = adata.X.T
if issparse(data_mat):
    if data_mat.nnz > 2**31 - 1:
        data_mat = data_mat.tocoo()
    else:
        data_mat = data_mat.tocsc()
ro.globalenv["data_mat"] = data_mat
ro.globalenv["input_groups"] = adata.obs["leiden"]


# In[26]:


get_ipython().run_cell_magic('R', '-o size_factors', '\nsize_factors = sizeFactors(\n    computeSumFactors(\n        SingleCellExperiment(\n            list(counts=data_mat)), \n            clusters = input_groups,\n            min.mean = 0.1,\n            BPPARAM = MulticoreParam()\n    )\n)\n')


# In[27]:


from scipy.sparse import csr_matrix
adata.obs["size_factors"] = size_factors
scran = adata.X / adata.obs["size_factors"].values[:, None]
scran = scran.tocsc()
adata.layers["scran_normalization"] = csr_matrix(sc.pp.log1p(scran))


# In[28]:


fig, axes = plt.subplots(1, 2, figsize=(10, 5))
p1 = sns.histplot(adata.obs["n_counts"], bins=100, kde=False, ax=axes[0])
axes[0].set_title("Total counts")
p2 = sns.histplot(
    adata.layers["scran_normalization"].sum(1), bins=100, kde=False, ax=axes[1]
)
axes[1].set_title("log1p with Scran estimated size factors")
plt.show()


# In[29]:


adata.write("raw_LMX_bipca.h5ad")
adata


# In[2]:


bdata = sc.read("raw_LMX_bipca.h5ad")


# In[3]:


###adata here is the original subset from the whole atlas
# filtered object (has embeddings in obsm)
adata_filt = bdata                         #ran BiPCA
# PRE-filtered object (same cells, more genes)
adata_full = sc.read("allgenes_bipca_hypomap.h5ad")

# 1) align cells (order & intersection) so n_obs matches
common_cells = adata_filt.obs_names.intersection(adata_full.obs_names)
adata_filt   = adata_filt[common_cells].copy()
adata_full   = adata_full[common_cells].copy()

# 2) choose what goes into raw.X (use the full scran-normalized matrix)
if "scran_normalization" in adata_full.layers:
    adata_full.X = adata_full.layers["scran_normalization"]
# (optional) make sure it's sparse to save RAM
if not sparse.issparse(adata_full.X):
    adata_full.X = sparse.csr_matrix(adata_full.X)

# 3) attach the full matrix (with all genes) as .raw
adata_filt.raw = adata_full


# In[4]:


# --- Plot previously excluded genes ---
# These genes can be absent from adata_filt.var_names as long as they exist in adata_filt.raw.var_names
sc.pl.embedding(
    adata_filt,
    basis="tsne",                                # or "pcs", "U_biwhite" if registered as embeddings
    color=["Glp1r"],
    use_raw=True                                 # <- key bit
)


# In[58]:


marker_genes = ["Lmx1a", "Lmx1b", "Irx3", "Irx5"]
sc.pl.embedding(adata_filt,
                basis="tsne",
                color=marker_genes,
                layer="scran_normalization",
                color_map="YlOrBr", 
                vmax=1.5)


# In[5]:


sc.pl.embedding(
    adata_filt,
    basis="tsne",
    color="detailed_celltype_annotation",
    legend_loc="on data",
    legend_fontsize="xx-small"
)


# In[52]:


###This part can't be saved so it's only ran once###
adata.X = adata.layers["scran_normalization"]
sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon", key_added="dea_leiden")


# In[53]:


sc.tl.dendrogram(adata, groupby="leiden", use_rep = "pcs", n_pcs=adata.uns["bipca"]['rank'])
sc.pl.rank_genes_groups_dotplot(adata, groupby="leiden", standard_scale="var", n_genes=5, key="dea_leiden", save="markers_LMX.pdf")


# In[54]:


sc.tl.filter_rank_genes_groups(
    adata,
    min_in_group_fraction=0.2,
    max_out_group_fraction=0.2,
    key="dea_leiden",
    key_added="dea_leiden_filtered"
)
#To preserve the original structure of adata.uns[‘rank_genes_groups’], filtered genes are set to NaN.
#This needs to be removed because it raises an error while saving to h5ad.


# In[55]:


sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby="leiden",
    standard_scale="var",
    n_genes=5,
    key="dea_leiden_filtered",
    save="markers_LMX_filtered.pdf"
)


# In[60]:


df = sc.get.rank_genes_groups_df(adata, group=None, key = "dea_leiden_filtered", log2fc_min=2)
df.to_csv("LMX_leiden_filtered.csv")
###This part can't be saved so it's only ran once###


# In[6]:


#Markers >20% in cluster <20% outside cluster >2LogFchange
cl_annotation = adata_filt.obs["leiden"]
cl_annotation["0"] = "00_Sox6_male"
cl_annotation["1"] = "01_Slc6a3_female"
cl_annotation["2"] = "02_Calcr"
cl_annotation["3"] = "03_Fgf10"
cl_annotation["4"] = "04_Nts"
cl_annotation["5"] = "05_Cubn"
cl_annotation["6"] = "06_Acvr1c_Npy2r"
cl_annotation["7"] = "07_Ebf1"
cl_annotation["8"] = "08_Gfra1"
cl_annotation["9"] = "09_Cntnap3"
cl_annotation["10"] = "10_Trpc3"
cl_annotation["11"] = "11_Rxfp1"
cl_annotation["12"] = "12_Cdh23"
cl_annotation["13"] = "13_Eya1"
cl_annotation["14"] = "14_Ebf3"
cl_annotation["15"] = "15_Nfia_Nfib"
cl_annotation["16"] = "drop"
cl_annotation["17"] = "drop"
adata_filt.obs["subset_celltype_annotation"] = adata_filt.obs.leiden.map(cl_annotation)


# In[7]:


sc.pl.embedding(adata_filt,
                basis="tsne",
                color="subset_celltype_annotation", 
                legend_loc="on data", 
                legend_fontsize="xx-small", 
                save="_LMX_celltypes.pdf")


# In[8]:


adata = adata_filt[adata_filt.obs["subset_celltype_annotation"]!="drop"]
sc.pl.embedding(adata,basis="tsne",color="subset_celltype_annotation", legend_loc="on data", legend_fontsize="xx-small")


# In[9]:


adata.write("LMX_bipca.h5ad")
adata


# In[2]:


adata = sc.read("LMX_bipca.h5ad")
adata


# In[3]:


gpcr_genes = [
    "Glp1r",
    "Mc4r",
    "Npy1r",
    "Npy2r",
    "Calcr",
    "Oxtr",
    "Ghsr",
    "Gipr",
    "Lepr"]
sc.pl.embedding(adata,basis="tsne",
                color=gpcr_genes, 
                color_map="YlOrBr", 
                vmax=1.5,
                use_raw=True)


# In[10]:


table = pd.crosstab(adata.obs["sex"], adata.obs["subset_celltype_annotation"])
# Proportion of each cell type *per animal*
row_props = table.div(table.sum(axis=1), axis=0)
# Proportion of each animal *per cell type*
col_props = table.div(table.sum(axis=0), axis=1)

chi2, p, dof, expected = chi2_contingency(table)
residuals = table - expected
n_rows, n_cols = table.shape
fig_width = max(10, n_cols * 0.5)
fig_height = max(6, n_rows * 0.5)
plt.figure(figsize=(fig_width, fig_height))
sns.heatmap(residuals, cmap="coolwarm", center=0, annot=True, fmt=".0f")
plt.title("Cell Type Overrepresentation (Observed - Expected)")
plt.ylabel("Animal")
plt.xlabel("Cell Type")
plt.tight_layout()
plt.show()


# In[11]:


std_residuals = (table - expected) / (expected**0.5)
n_rows, n_cols = table.shape
fig_width = max(10, n_cols * 0.5)
fig_height = max(6, n_rows * 0.5)
plt.figure(figsize=(fig_width, fig_height))
sns.heatmap(std_residuals, cmap="coolwarm", center=0, annot=True, fmt=".0f")
plt.title("Cell Type Standardized Overrepresentation (Observed - Expected)")
plt.ylabel("Animal")
plt.xlabel("Cell Type")
plt.tight_layout()
plt.show()


# In[12]:


ct_col  = "subset_celltype_annotation"
sex_col = "sex"
tab = pd.crosstab(adata.obs[sex_col], adata.obs[ct_col])

# expected counts under independence
N = tab.values.sum()
expected = np.outer(tab.sum(1).values, tab.sum(0).values) / N
expected = pd.DataFrame(expected, index=tab.index, columns=tab.columns)

# Pearson standardized residuals: (obs - exp) / sqrt(exp)
# (guard against tiny expected counts)
eps = 1e-9
resid = (tab - expected) / np.sqrt(expected + eps)   # DataFrame same shape as tab
# turn (celltype, sex) residual grid into a lookup table
look = (
    resid.stack()
         .rename("sex_resid")
         .reset_index()
         .rename(columns={"level_0": ct_col, "level_1": sex_col})
)

# join back to obs
adata.obs = adata.obs.merge(
    look, on=[ct_col, sex_col], how="left"
)


# In[13]:


sc.pl.embedding(
    adata,
    basis="tsne",
    color="sex_resid",
    cmap="coolwarm",
    vcenter=0,                  # 0 = no enrichment; neg under-, pos over-representation
    size=5
)


# In[14]:


from pandas.api.types import is_categorical_dtype
# make a single signed score: positive = male overrep, negative = female overrep
sex_bias = resid.loc['male'].subtract(resid.loc['female'], fill_value=0)/2
# sex_bias is a Series indexed by cell type

# map that score back to each cell via its cell type
adata.obs['sex_bias'] = adata.obs[ct_col].map(sex_bias)
col = "sex_bias"
# Coerce to numeric floats (and drop any categorical dtype)
adata.obs[col] = pd.to_numeric(adata.obs[col], errors="coerce")
# (if it was a categorical dtype, the line above removes that;
#  otherwise, you can explicitly do:)
if is_categorical_dtype(adata.obs[col]):
    adata.obs[col] = adata.obs[col].astype(float)


# In[15]:


from matplotlib.colors import TwoSlopeNorm
sc.pl.embedding(
    adata,
    basis='tsne',
    color='sex_bias',
    cmap='PiYG',                    
    vcenter=0,
    size=6,
    frameon=False,
    title='Pearson residuals (male→green, female→violet)',
    save="LMX_male_female.pdf"
)


# In[16]:


col = "detailed_celltype_annotation"
selected = ['04_Crnde_Lmx1a','27_Lmx1a_Ttc6','32_Gucy2c_Ttc6','36_Ebf2_Crhbp']
selected_colors = {'04_Crnde_Lmx1a': '#008941',
 '27_Lmx1a_Ttc6': '#ffaa92',
 '32_Gucy2c_Ttc6': '#000035',
 '36_Ebf2_Crhbp': '#0aa6d8'}
dtype = pd.CategoricalDtype(categories=selected, ordered=True)
adata.obs[col] = adata.obs[col].astype(dtype)
cats_now = list(adata.obs[col].cat.categories)
adata.uns[f"{col}_colors"] = [selected_colors.get(c, "#bdbdbd") for c in cats_now]


# In[17]:


sc.pl.embedding(adata,
                basis="tsne",
                color="detailed_celltype_annotation", 
                legend_loc="on data", 
                legend_fontsize="xx-small",
                save="LMX_subset_detailed.pdf")


# In[12]:


sc.pl.embedding(adata,basis="tsne",
                color="Lepr", 
                color_map="YlOrBr", 
                vmax=1.5)


# In[18]:


adata.write("LMX_bipca.h5ad")
#adata = sc.read("LMX_bipca.h5ad")
adata


# In[3]:


gpcr_genes = [
    "Slc17a6", #VGlut2
    "Slc17a7", #VGlut1 LHA https://www.cell.com/neuron/fulltext/S0896-6273%2821%2900693-0
    "Gad1",
    "Gad2",
    "Slc32a1", #Vgat
    "Hcrt", #Ox, should be Vglut2 positive
    "Pmch", #MCH should be Gad1 positive
    "Nts", #LH
    "Gal",
    "Crh",
    "Pdyn", 
    "Tac1", 
    "Tac2", 
    "Pnoc", #LH
    "Hdc", #Histidine decarboxylase, Tuberomammillary
    "Slc18a2", #Vmat2, Tuberomammillary
    "Lepr", #Premammillary reproductive
    "Gucy2c",
    "Slc6a3", #Dopamine
    "Th",
    "Ddc"
    ]
sc.pl.embedding(adata,basis="tsne",
                color=gpcr_genes, 
                color_map="YlOrBr", 
                vmax=1.5,
                use_raw=True,
                save="_LMX_subset_classifiers")


# In[4]:


gpcr_genes = [
    "Lepr",
    "Insr",
    "Ghsr", #Ghrelin
    "Glp1r",
    "Gipr",
    "Mc4r",
    "Mc3r",
    "Npy1r",
    "Npy2r",
    "Npy5r",
    "Mchr1",
    "Hcrtr1", #Orexin
    "Hcrtr2", 
    "Htr2c", #Serotonin
    "Cnr1", #Cannabinoid
    "Oprm1", #Opioid µ
    "Hrh1", #Histamine
    "Hrh2",
    "Hrh3",
    "Calcr",
    "Ramp1", 
    "Ramp3"
    ]
sc.pl.embedding(adata,basis="tsne",
                color=gpcr_genes, 
                color_map="YlOrBr", 
                vmax=1.5,
                use_raw=True,
                save="_LMX_subset_receptors")


# In[5]:


bmi_genes = [
    "Ptprn",
    "Mllt11",
    "Arhgdig",
    "Lingo2",
    "Kcnb2",
    "Mast1",
    "Aldoa",
    "Dnm1",
    "Efna5",
    "Calcr",
    "Syt4",
    "Dlk1",
    "Atcay", 
    "Pacsin1",
    "Elavl4",
    "Cd151",
    "Slc35g2",
    "Rac3",
    "Moap1",
    "Sez6l2"
    ]
sc.pl.embedding(adata,basis="tsne",
                color=gpcr_genes, 
                color_map="YlOrBr", 
                vmax=1.5,
                use_raw=True,
                save="_LMX_subset_BMI")


# In[4]:


sc.pl.embedding(
    adata,
    basis="tsne",                       
    color="Yengo_BMI",
    cmap="RdBu_r",
    save="_LMX_subset_scDRS"
)


# In[3]:


import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Define the columns of interest
c_levels = ['C1_HM', 'C2_HM', 'C3_HM', 'C4_HM']

fig, axes = plt.subplots(2, 2, figsize=(18, 14))
axes = axes.flatten()

for i, col in enumerate(c_levels):
    counts = adata.obs[col].value_counts()
    pct = (counts / counts.sum() * 100).round(2)

    df = pd.DataFrame({'count': counts, 'percent': pct})

    # Plot top 20 labels
    top = df.head(20)
    top['percent'].plot(kind='barh', ax=axes[i], color='steelblue', edgecolor='white')
    axes[i].set_title(f'{col}', fontsize=13, fontweight='bold')
    axes[i].set_xlabel('% of cells')
    axes[i].set_ylabel('')
    axes[i].invert_yaxis()

    # Annotate bars with counts
    for j, (pct_val, cnt) in enumerate(zip(top['percent'], top['count'])):
        axes[i].text(pct_val + 0.3, j, f'n={cnt}', va='center', fontsize=8)

plt.suptitle('Enriched Human Hypomap Labels in Lmx1a Subset', fontsize=15, fontweight='bold', y=1.01)
plt.tight_layout()
plt.savefig('LMX_hypomap_enrichment.pdf', bbox_inches='tight')
plt.show()


# In[4]:


# Also print a clean summary table for each level
for col in c_levels:
    counts = adata.obs[col].value_counts()
    pct = (counts / counts.sum() * 100).round(2)
    df = pd.DataFrame({'count': counts, 'percent (%)': pct}).head(15)
    print(f"\n{'='*50}")
    print(f"  {col}")
    print(f"{'='*50}")
    print(df.to_string())


# In[5]:


# Check what the detailed_celltype_annotation values are
celltypes = adata.obs['detailed_celltype_annotation'].unique()
# For each C level, build a crosstab: rows = HM label, cols = detailed_celltype_annotation
# Normalize by column (within each celltype) to get % composition

fig, axes = plt.subplots(2, 2, figsize=(22, 18))
axes = axes.flatten()

for i, col in enumerate(c_levels):
    # Crosstab: counts
    ct = pd.crosstab(adata.obs[col], adata.obs['detailed_celltype_annotation'])

    # Keep only top N labels overall to keep plot readable
    top_labels = adata.obs[col].value_counts().head(20).index
    ct = ct.loc[ct.index.isin(top_labels)]

    # Normalize within each cell type (column-wise %)
    ct_pct = ct.div(ct.sum(axis=0), axis=1) * 100

    # Sort rows by mean abundance across cell types
    ct_pct = ct_pct.loc[ct_pct.mean(axis=1).sort_values(ascending=False).index]

    sns.heatmap(
        ct_pct,
        ax=axes[i],
        cmap='YlOrRd',
        linewidths=0.3,
        linecolor='lightgrey',
        annot=True,
        fmt='.1f',
        annot_kws={'size': 7},
        cbar_kws={'label': '% within cell type'}
    )
    axes[i].set_title(f'{col}', fontsize=13, fontweight='bold')
    axes[i].set_xlabel('Detailed cell type annotation', fontsize=10)
    axes[i].set_ylabel('Human Hypomap Label', fontsize=10)
    axes[i].tick_params(axis='x', rotation=45, labelsize=9)
    axes[i].tick_params(axis='y', rotation=0, labelsize=8)

plt.suptitle('Human Hypomap Label Abundance per Detailed Cell Type\n(% within each cell type annotation)',
             fontsize=14, fontweight='bold', y=1.01)
plt.tight_layout()
plt.savefig('LMX_hypomap_by_celltype_heatmap.pdf', bbox_inches='tight')
plt.show()


# In[6]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Check the distribution of sex_bias scores first
print("sex_bias summary:")
print(adata.obs['sex_bias'].describe())
print(f"\nNaN count: {adata.obs['sex_bias'].isna().sum()}")


# In[15]:


# Visualise the distribution to pick sensible thresholds
fig, ax = plt.subplots(1, 1, figsize=(8, 4))
adata.obs['sex_bias'].hist(bins=50, ax=ax, color='steelblue', edgecolor='white')
ax.axvline(0, color='black', linestyle='--', linewidth=1)
# Common Pearson residual thresholds: ±2 (similar to z-score)
ax.axvline(2, color='firebrick', linestyle='--', linewidth=1, label='threshold +2')
ax.axvline(-2, color='royalblue', linestyle='--', linewidth=1, label='threshold -2')
ax.set_xlabel('Sex bias score (Pearson residual)')
ax.set_ylabel('Cell count')
ax.set_title('Distribution of sex bias scores')
ax.legend()
plt.tight_layout()
plt.show()


# In[16]:


# --- Define thresholds (adjust if needed after inspecting the histogram) ---
threshold = 10.0  # Pearson residuals: |r| > 2 is conventionally meaningful

def classify_sex_bias(score, thr=threshold):
    if pd.isna(score):
        return 'neutral'
    elif score > thr:
        return 'male'
    elif score < -thr:
        return 'female'
    else:
        return 'neutral'

adata.obs['sex_bias_category'] = adata.obs['sex_bias'].apply(classify_sex_bias)

# Summary
print("Sex bias category counts:")
print(adata.obs['sex_bias_category'].value_counts())
print("\nBy detailed_celltype_annotation:")
print(pd.crosstab(adata.obs['detailed_celltype_annotation'], 
                  adata.obs['sex_bias_category']))


# In[17]:


# Visualise: stacked bar per cell type + scatter of underlying scores
fig, axes = plt.subplots(1, 2, figsize=(16, 5))

palette = {'male': '#4878CF', 'female': '#D65F5F', 'neutral': '#A8A8A8'}

# --- Left: stacked bar showing proportion per detailed cell type ---
ct = pd.crosstab(
    adata.obs['detailed_celltype_annotation'],
    adata.obs['sex_bias_category'],
    normalize='index'
) * 100

# Ensure consistent column order
for cat in ['male', 'neutral', 'female']:
    if cat not in ct.columns:
        ct[cat] = 0
ct = ct[['male', 'neutral', 'female']]

ct.plot(kind='bar', stacked=True, ax=axes[0],
        color=[palette[c] for c in ct.columns],
        edgecolor='white', linewidth=0.5)
axes[0].set_title('Sex bias composition\nper cell type', fontsize=12, fontweight='bold')
axes[0].set_ylabel('% of cells')
axes[0].set_xlabel('')
axes[0].tick_params(axis='x', rotation=45)
axes[0].legend(title='Sex bias', bbox_to_anchor=(1.01, 1), loc='upper left')

# --- Right: boxplot of raw sex_bias scores per cell type, coloured by category median ---
order = adata.obs.groupby('detailed_celltype_annotation')['sex_bias'] \
               .median().sort_values().index.tolist()

sns.boxplot(
    data=adata.obs,
    x='detailed_celltype_annotation', y='sex_bias',
    order=order,
    ax=axes[1],
    palette='coolwarm',
    linewidth=0.8,
    flierprops=dict(marker='o', markersize=2, alpha=0.3)
)
axes[1].axhline(threshold, color='firebrick', linestyle='--', linewidth=1, label=f'+{threshold}')
axes[1].axhline(-threshold, color='royalblue', linestyle='--', linewidth=1, label=f'-{threshold}')
axes[1].axhline(0, color='black', linestyle='-', linewidth=0.8)
axes[1].set_title('Raw sex bias scores\nper cell type', fontsize=12, fontweight='bold')
axes[1].set_ylabel('Sex bias (Pearson residual)')
axes[1].set_xlabel('')
axes[1].tick_params(axis='x', rotation=45)
axes[1].legend(title='Threshold')

plt.suptitle('Sex Bias in Lmx1a Neuron Subtypes', fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig('LMX_sex_bias_categories.pdf', bbox_inches='tight')
plt.show()


# In[18]:


# Check what the detailed_celltype_annotation values are
celltypes = adata.obs['sex_bias_category'].unique()
# For each C level, build a crosstab: rows = HM label, cols = detailed_celltype_annotation
# Normalize by column (within each celltype) to get % composition

fig, axes = plt.subplots(2, 2, figsize=(22, 18))
axes = axes.flatten()

for i, col in enumerate(c_levels):
    # Crosstab: counts
    ct = pd.crosstab(adata.obs[col], adata.obs['sex_bias_category'])

    # Keep only top N labels overall to keep plot readable
    top_labels = adata.obs[col].value_counts().head(20).index
    ct = ct.loc[ct.index.isin(top_labels)]

    # Normalize within each cell type (column-wise %)
    ct_pct = ct.div(ct.sum(axis=0), axis=1) * 100

    # Sort rows by mean abundance across cell types
    ct_pct = ct_pct.loc[ct_pct.mean(axis=1).sort_values(ascending=False).index]

    sns.heatmap(
        ct_pct,
        ax=axes[i],
        cmap='YlOrRd',
        linewidths=0.3,
        linecolor='lightgrey',
        annot=True,
        fmt='.1f',
        annot_kws={'size': 7},
        cbar_kws={'label': '% within category'}
    )
    axes[i].set_title(f'{col}', fontsize=13, fontweight='bold')
    axes[i].set_xlabel('Sex bias category', fontsize=10)
    axes[i].set_ylabel('Human Hypomap Label', fontsize=10)
    axes[i].tick_params(axis='x', rotation=45, labelsize=9)
    axes[i].tick_params(axis='y', rotation=0, labelsize=8)

plt.suptitle('Human Hypomap Label Abundance by sex\n(% within each category)',
             fontsize=14, fontweight='bold', y=1.01)
plt.tight_layout()
plt.savefig('LMX_hypomap_by_sex_heatmap.pdf', bbox_inches='tight')
plt.show()


# In[22]:


import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Define the columns of interest
c_levels = ['region_HM']

# Create a single subplot - unpack properly
fig, ax = plt.subplots(figsize=(18, 14))

# Process the data
col = c_levels[0]
counts = adata.obs[col].value_counts()
pct = (counts / counts.sum() * 100).round(2)

df = pd.DataFrame({'count': counts, 'percent': pct})

# Plot top 20 labels
top = df.head(20)
top['percent'].plot(kind='barh', ax=ax, color='steelblue', edgecolor='white')
ax.set_title(f'{col}', fontsize=13, fontweight='bold')
ax.set_xlabel('% of cells')
ax.set_ylabel('')
ax.invert_yaxis()

# Annotate bars with counts
for j, (pct_val, cnt) in enumerate(zip(top['percent'], top['count'])):
    ax.text(pct_val + 0.3, j, f'n={cnt}', va='center', fontsize=8)

plt.suptitle('Enriched Human Hypomap Regions in Lmx1a Subset', fontsize=15, fontweight='bold', y=1.01)
plt.tight_layout()
plt.savefig('LMX_hypomap_enrichment_region.pdf', bbox_inches='tight')
plt.show()


# In[23]:


# --- Left: stacked bar showing proportion per detailed cell type ---
ct = pd.crosstab(
    adata.obs['region_HM'],
    adata.obs['sex_bias_category'],
    normalize='index'
) * 100

# Ensure consistent column order
for cat in ['male', 'neutral', 'female']:
    if cat not in ct.columns:
        ct[cat] = 0
ct = ct[['male', 'neutral', 'female']]
# Check what the detailed_celltype_annotation values are
celltypes = adata.obs['sex_bias_category'].unique()

# For each C level, build a crosstab: rows = HM label, cols = detailed_celltype_annotation
# Normalize by column (within each celltype) to get % composition
fig, axes = plt.subplots(nrows=len(c_levels), figsize=(18, 14))

# Handle single subplot case - axes won't be an array
if len(c_levels) == 1:
    axes = [axes]

for i, col in enumerate(c_levels):
    # Crosstab: counts
    ct = pd.crosstab(adata.obs[col], adata.obs['sex_bias_category'])

    # Keep only top N labels overall to keep plot readable
    top_labels = adata.obs[col].value_counts().head(20).index
    ct = ct.loc[ct.index.isin(top_labels)]

    # Normalize within each cell type (column-wise %)
    ct_pct = ct.div(ct.sum(axis=0), axis=1) * 100

    # Sort rows by mean abundance across cell types
    ct_pct = ct_pct.loc[ct_pct.mean(axis=1).sort_values(ascending=False).index]

    sns.heatmap(
        ct_pct,
        ax=axes[i],
        cmap='YlOrRd',
        linewidths=0.3,
        linecolor='lightgrey',
        annot=True,
        fmt='.1f',
        annot_kws={'size': 7},
        cbar_kws={'label': '% within category'}
    )
    axes[i].set_title(f'{col}', fontsize=13, fontweight='bold')
    axes[i].set_xlabel('Sex bias category', fontsize=10)
    axes[i].set_ylabel('Human Hypomap Label', fontsize=10)
    axes[i].tick_params(axis='x', rotation=45, labelsize=9)
    axes[i].tick_params(axis='y', rotation=0, labelsize=8)

plt.suptitle('Human Hypomap Region Abundance by sex\n(% within each category)',
             fontsize=14, fontweight='bold', y=1.01)
plt.tight_layout()
plt.savefig('LMX_hypomap_by_sex_heatmap_region.pdf', bbox_inches='tight')
plt.show()


# In[7]:


"""
Cross-species validation of Lmx1a BMI enrichment — CORRECTED
Uses LDSC results throughout. Group A split by LDSC nominal (p<0.05).
Sankey order used for all C4 plots.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns
import scanpy as sc

UNIFIED_PATH  = "hypomap_cellect_unified.csv"
EFFECTOR_PATH = "genetics_effector.txt"
EXOME_PATH    = "genetics_exome.txt"
GWAS          = 'BMI_Yengo2018'
BONF_C4       = 0.05 / 415
BONF_C3       = 0.05 / 156
BONF_EXOME    = 0.05 / 426
NOM_C4        = 0.05            # LDSC nominal threshold for Group A split

OUTDIR = "figures_hypomap_corrected"
os.makedirs(OUTDIR, exist_ok=True)

sns.set_style("whitegrid")
plt.rcParams.update({'font.size':9,'axes.titlesize':11,'axes.labelsize':10,
                     'figure.dpi':150,'savefig.dpi':300,'savefig.bbox':'tight',
                     'pdf.fonttype':42,'ps.fonttype':42})

MOUSE_CELLTYPES = ['04_Crnde_Lmx1a','27_Lmx1a_Ttc6','32_Gucy2c_Ttc6','36_Ebf2_Crhbp']
MOUSE_COLORS    = {'04_Crnde_Lmx1a':'#008941','27_Lmx1a_Ttc6':'#ffaa92',
                   '32_Gucy2c_Ttc6':'#000035','36_Ebf2_Crhbp':'#0aa6d8'}

GROUP_A_C4 = ['C4-226','C4-231','C4-247','C4-248','C4-229',
               'C4-223','C4-224','C4-252','C4-228','C4-230']
GROUP_B_C4 = ['C4-33','C4-407','C4-409','C4-411']
LMX_C4_IDS = GROUP_A_C4 + GROUP_B_C4

LMX_C3_IDS = ['C3-92','C3-84','C3-139','C3-29','C3-94','C3-83','C3-82','C3-26',
               'C3-137','C3-85','C3-28','C3-24','C3-93','C3-88','C3-115','C3-95']

GROUP_A_C3       = ['C3-92','C3-93','C3-84','C3-94','C3-82','C3-85','C3-83']
CRNDE_C3         = ['C3-95']
GABA_NEIGHBOR_C3 = ['C3-28']
GROUP_B_C3       = ['C3-29','C3-137','C3-88','C3-139','C3-26','C3-115','C3-24']

# ── Simple colour scheme: Group A nom / Group A not-nom / Group B ─────────────
COL_A_NOM  = '#8B0000'   # Group A LDSC nominal (p<0.05)
COL_A_NS   = '#E8A8A8'   # Group A not nominal
COL_B      = '#95a5a6'   # Group B (LMX1A-)
COL_HL     = '#8B0000'   # highlight border colour

# ── 1. Load unified CSV ───────────────────────────────────────────────────────
print("Loading unified CELLECT results …")
unified = pd.read_csv(UNIFIED_PATH)

def get_bmi(method, level):
    return (unified[(unified['gwas']==GWAS) & (unified['method']==method) &
                    (unified['level']==level) & (unified['is_neuronal'])]
            .set_index('cluster_id').copy())

# Primary: LDSC throughout
ldsc_c4  = get_bmi('new_LDSC',  'C4')
ldsc_c3  = get_bmi('new_LDSC',  'C3')
magma_c3 = get_bmi('new_MAGMA', 'C3')   # kept for C3 Bonf-sig markers only

print(f"  new_LDSC C4: {len(ldsc_c4)} clusters  "
      f"({ldsc_c4['bonf_sig'].sum()} Bonf-sig  "
      f"{(ldsc_c4['pvalue']<NOM_C4).sum()} nominal)")

# ── 2. Build C4 summary table ─────────────────────────────────────────────────
records = []
for cid in LMX_C4_IDS:
    if cid not in ldsc_c4.index: continue
    r = ldsc_c4.loc[cid]
    records.append({
        'cluster_id':  cid,
        'annotation':  r['annotation_full'],
        'group':       'A' if cid in GROUP_A_C4 else 'B',
        'beta':        r['beta'],
        'pvalue':      r['pvalue'],
        'log10p':      r['log10p'],
        'bonf_sig':    r['bonf_sig'],
        'nom_sig':     r['pvalue'] < NOM_C4,
        'n_effector':  r['n_effector_genes'],
        'has_exome':   r['has_exome_sig'],
        'exome_genes': r['exome_sig_genes'],
    })

df_c4 = pd.DataFrame(records)
df_c4['is_grp_a'] = df_c4['cluster_id'].isin(GROUP_A_C4)

# Group A LDSC nominal clusters → effector gene analysis (expect 7/10)
grp_a_nom = df_c4[df_c4['is_grp_a'] & df_c4['nom_sig']]['cluster_id'].tolist()
print(f"  Group A LDSC nominal (p<0.05): {len(grp_a_nom)} / {len(GROUP_A_C4)}")
print(f"  → {grp_a_nom}")
# ── After computing grp_a_nom, add FDR computation ────────────────────────────
from statsmodels.stats.multitest import multipletests

all_pvals = ldsc_c4['pvalue'].values
_, fdr_pvals, _, _ = multipletests(all_pvals, method='fdr_bh')
ldsc_c4['fdr'] = fdr_pvals

fdr_sig_mask = fdr_pvals < 0.05
if fdr_sig_mask.any():
    fdr_cutoff_p   = all_pvals[fdr_sig_mask].max()
    fdr_cutoff_log = -np.log10(fdr_cutoff_p)
else:
    fdr_cutoff_p = fdr_cutoff_log = None
    print("No clusters pass FDR<0.05")

# Update grp_a_nom to use FDR instead of nominal p=0.05
grp_a_nom = [c for c in GROUP_A_C4
             if c in ldsc_c4.index and ldsc_c4.loc[c, 'fdr'] < 0.05]
print(f"  Group A FDR<0.05: {len(grp_a_nom)} / {len(GROUP_A_C4)}")
print(f"  → {grp_a_nom}")
# ── 3. Cell counts + Sankey order ─────────────────────────────────────────────
print("\nExtracting label transfer counts …")
if 'C4_short' not in adata.obs.columns:
    adata.obs['C4_short'] = adata.obs['C4_HM'].str.extract(r'(C4-\d+)')
if 'C3_short' not in adata.obs.columns:
    adata.obs['C3_short'] = adata.obs['C3_HM'].str.extract(r'(C3-\d+)')

lmx_c4_counts = adata.obs['C4_short'].value_counts()
lmx_c3_counts = adata.obs['C3_short'].value_counts()

# Sankey order: Group A desc by mouse count, then Group B desc by mouse count
c4_order_sankey = (
    sorted(GROUP_A_C4, key=lambda c: -lmx_c4_counts.get(c, 0)) +
    sorted(GROUP_B_C4, key=lambda c: -lmx_c4_counts.get(c, 0))
)

df_c4['n_mouse_cells'] = df_c4['cluster_id'].map(lmx_c4_counts).fillna(0).astype(int)

spec_map = {
    'C4-226':100.0,'C4-231':100.0,'C4-247': 99.3,'C4-248': 98.3,
    'C4-229': 97.2,'C4-223': 97.0,'C4-224': 90.6,'C4-252': 88.4,
    'C4-228': 61.1,'C4-230': 59.0,'C4-33':  27.7,'C4-407': 22.4,
    'C4-409': 22.2,'C4-411': 16.0,
}
df_c4['pct_lmx1a_specific'] = df_c4['cluster_id'].map(spec_map).fillna(0)

# ── 4. Effector genes — Group A LDSC nominal only ─────────────────────────────
print("Loading effector genes …")
eff = pd.read_csv(EFFECTOR_PATH, sep='\t')
eff['cluster_id'] = (eff['annotation']
                     .str.replace(r'^C4_','C4-',regex=True)
                     .str.replace(r'^C3_','C3-',regex=True))

eff_nom = eff[eff['cluster_id'].isin(grp_a_nom)].copy()
gene_summary = (eff_nom.groupby('gene_symbol')
                .agg(n_clusters=('cluster_id','nunique'),
                     mean_esmu=('esmu','mean'),
                     min_magma_p=('magma_gene_pval','min'))
                .sort_values(['n_clusters','min_magma_p'], ascending=[False,True]))
top_effectors = gene_summary[gene_summary['n_clusters'] >= 3].head(20)
print(f"  Effector genes in ≥3 Group A nominal clusters: {len(top_effectors)}")

# ── 5. Exome convergence ──────────────────────────────────────────────────────
print("Loading exome results …")
exome = pd.read_csv(EXOME_PATH, sep='\t')
exome['cluster_id_list'] = exome['cluster'].apply(
    lambda s: [x.strip().replace('C4_','C4-').replace('C3_','C3-')
               for x in str(s).split(',')])
exome_records = []
for _, row in exome.iterrows():
    cl = set(row['cluster_id_list'])
    overlap_nom = cl & set(grp_a_nom)
    overlap_all = cl & set(LMX_C4_IDS)
    if overlap_all:
        exome_records.append({
            'gene':          row['gene'],
            'p_exome':       row['p_exome'],
            'mask':          row['mask_exome'],
            'gwas_signal':   row['GWAS_signal'],
            'n_nom':         len(overlap_nom),
            'n_lmx_total':   len(overlap_all),
            'nom_clusters':  ','.join(sorted(overlap_nom)),
        })
df_exome_lmx = pd.DataFrame(exome_records).sort_values('p_exome')

# ══════════════════════════════════════════════════════════════════════════════
#  FIGURES
# ══════════════════════════════════════════════════════════════════════════════

def bar_color(cid):
    if cid in set(GROUP_B_C4):  return COL_B
    if cid in set(grp_a_nom):   return COL_A_NOM   # now FDR-based
    return COL_A_NS

# ── FIGURE A: C4 LDSC bar chart — Sankey order, flipped ──────────────────────
fig, ax = plt.subplots(figsize=(6, 6.5))

# Sankey order bottom-to-top (reversed for horizontal bar)
plot_order = list(reversed(c4_order_sankey))
df_c4_idx  = df_c4.set_index('cluster_id')

colors = [bar_color(c) for c in plot_order]
edges  = ['white' for _ in plot_order]   # no highlights
lws    = [0.4     for _ in plot_order]
vals   = [ldsc_c4.loc[c,'log10p'] if c in ldsc_c4.index else 0.0 for c in plot_order]

ax.barh(range(len(plot_order)), vals,
        color=colors, edgecolor=edges, linewidth=lws, height=0.72)

if fdr_cutoff_log is not None:
    ax.axvline(fdr_cutoff_log, color='steelblue', ls='--', lw=1)

ax.set_yticks(range(len(plot_order)))
ax.set_yticklabels(plot_order, fontsize=8)
for lbl, cid in zip(ax.get_yticklabels(), plot_order):
    if cid in set(GROUP_B_C4):
        lbl.set_color('#7f8c8d')
    # no bold/colour for C4-229 or C4-248

for i, cid in enumerate(plot_order):
    spec = spec_map.get(cid, 0)
    v    = vals[i]
    ax.text(v + 0.03, i, f'{spec:.0f}%', va='center', fontsize=6, color='grey')

ax.set_xlabel('−log₁₀(new_LDSC p-value)  [Yengo2018 BMI]')
ax.set_title('BMI enrichment of human HYPOMAP C4 clusters\n'
             'corresponding to mouse Lmx1a neurons  (LDSC)',
             fontweight='bold')
ax.legend(handles=[
    Patch(fc=COL_A_NOM, label=f'Group A  LMX1A+  FDR<0.05 ({len(grp_a_nom)}/10)'),
    Patch(fc=COL_A_NS,  label='Group A  LMX1A+  not significant'),
    Patch(fc=COL_B,     label='Group B  LMX1A−'),
    plt.Line2D([0],[0], color='steelblue', ls='--', lw=1,
               label=f'FDR<0.05 (p={fdr_cutoff_p:.2f})' if fdr_cutoff_log else 'FDR<0.05 (none)'),
], fontsize=7.5, loc='lower right', framealpha=0.9)
plt.tight_layout()
plt.savefig(f'{OUTDIR}/Fig_A_LDSC_BMI_C4.pdf')
plt.savefig(f'{OUTDIR}/Fig_A_LDSC_BMI_C4.png')
print("  Saved Fig A")

# ── FIGURE B: Context histogram ───────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(6, 3.5))
ax.hist(ldsc_c4['log10p'], bins=60, color='#D0D0D0', edgecolor='white',
        linewidth=0.3, label=f'All C4 (n={len(ldsc_c4)})')
ax.hist([ldsc_c4.loc[c,'log10p'] for c in grp_a_nom if c in ldsc_c4.index],
        bins=10, color=COL_A_NOM, alpha=0.85, edgecolor='white', linewidth=0.3,
        label=f'Group A nominal LDSC ({len(grp_a_nom)})')
ax.hist([ldsc_c4.loc[c,'log10p'] for c in GROUP_A_C4
         if c in ldsc_c4.index and c not in grp_a_nom],
        bins=5, color=COL_A_NS, alpha=0.85, edgecolor='white', linewidth=0.3,
        label='Group A not nominal')
ax.hist([ldsc_c4.loc[c,'log10p'] for c in GROUP_B_C4 if c in ldsc_c4.index],
        bins=5, color=COL_B, alpha=0.85, edgecolor='white', linewidth=0.3,
        label=f'Group B LMX1A− ({len(GROUP_B_C4)})')
ax.axvline(-np.log10(NOM_C4), color='grey', ls=':', lw=0.8)
ax.set_xlabel('−log₁₀(new_LDSC p)  |  Yengo2018 BMI')
ax.set_ylabel('Number of C4 clusters')
ax.set_title('Group A LMX1A+ clusters enriched for BMI signal\n'
             'Group B LMX1A− scatter with the bulk', fontweight='bold')
ax.legend(fontsize=8)
plt.tight_layout()
plt.savefig(f'{OUTDIR}/Fig_B_context_histogram.pdf')
plt.savefig(f'{OUTDIR}/Fig_B_context_histogram.png')
print("  Saved Fig B")

# ── FIGURE C: Effector gene heatmap — Group A LDSC nominal only ───────────────
fig, ax = plt.subplots(figsize=(8, 6))

top_genes = top_effectors.index.tolist()[:20]
mat = pd.DataFrame(0.0, index=top_genes, columns=grp_a_nom)
for _, row in eff_nom[eff_nom['gene_symbol'].isin(top_genes)].iterrows():
    if row['cluster_id'] in mat.columns:
        mat.loc[row['gene_symbol'], row['cluster_id']] = row['esmu']

gene_order = mat.apply(lambda r: (r>0).sum(), axis=1).sort_values(ascending=False).index
mat = mat.loc[gene_order]

col_labels = [
    f"{c}\n{ldsc_c4.loc[c,'annotation_full'][4:].replace('_',' ')[:25]}\n"
    f"p={ldsc_c4.loc[c,'pvalue']:.1e}"
    for c in mat.columns
]
sns.heatmap(mat, ax=ax, cmap='YlOrRd', linewidths=0.4, linecolor='lightgrey',
            cbar_kws={'label':'ESµ specificity','shrink':0.7},
            xticklabels=col_labels, yticklabels=True,
            mask=(mat==0), vmin=0, vmax=0.8)
ax.set_title('Top effector genes in Group A LMX1A+ clusters\n'
             '(LDSC nominal p<0.05, Yengo2018)', fontweight='bold', pad=12)
ax.set_xlabel('Human C4 cluster  (p = new_LDSC)', fontsize=9)
ax.set_ylabel('BMI effector gene', fontsize=9)
ax.tick_params(axis='x', rotation=45, labelsize=6.5)
plt.tight_layout()
plt.savefig(f'{OUTDIR}/Fig_C_effector_heatmap.pdf')
plt.savefig(f'{OUTDIR}/Fig_C_effector_heatmap.png')
print("  Saved Fig C")

# ── FIGURE D: Exome convergence ───────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(5, 5))
df_ex = df_exome_lmx[df_exome_lmx['p_exome'].notna()].copy()
df_ex['log10p'] = -np.log10(df_ex['p_exome'].clip(lower=1e-50))
df_ex = df_ex.sort_values('log10p', ascending=False).head(20)

KNOWN_OBESITY = {'MC4R','PCSK1','POMC','CALCR','BSN','CORO1A'}
colors_ex = ['#8B0000' if g in KNOWN_OBESITY else '#4878CF' for g in df_ex['gene']]
ax.barh(range(len(df_ex)), df_ex['log10p'], color=colors_ex,
        edgecolor='white', linewidth=0.5, height=0.7)
ax.set_yticks(range(len(df_ex)))
ax.set_yticklabels(df_ex['gene'], fontsize=8)
ax.invert_yaxis()
ax.axvline(-np.log10(0.05),    color='grey',     ls=':', lw=0.8)
ax.axvline(-np.log10(BONF_EXOME), color='black', ls='--', lw=1,
           label=f'Exome-wide sig (p={BONF_EXOME:.1e})')
for i, (_, r) in enumerate(df_ex.iterrows()):
    ax.text(r['log10p']+0.1, i,
            f"nom={r['n_nom']}  tot={r['n_lmx_total']}",
            va='center', fontsize=6.5, color='grey')
ax.set_xlabel('−log₁₀(rare variant exome p-value)')
ax.set_title('Rare variant evidence for effector genes\n'
             'in Group A LMX1A+ nominal clusters', fontweight='bold')
ax.legend(handles=[Patch(fc='#8B0000', label='Known obesity gene'),
                   Patch(fc='#4878CF', label='Other effector gene')],
          fontsize=7.5, loc='lower right')
plt.tight_layout()
plt.savefig(f'{OUTDIR}/Fig_D_exome_convergence.pdf')
plt.savefig(f'{OUTDIR}/Fig_D_exome_convergence.png')
print("  Saved Fig D")

# ── FIGURE E: C3 cluster assignment bar chart ─────────────────────────────────
fig, ax = plt.subplots(figsize=(8, 4.5))
c3_total = lmx_c3_counts.sum()
c3_top   = [c for c in LMX_C3_IDS if c in lmx_c3_counts.index]
c3_pcts  = [lmx_c3_counts[c]/c3_total*100 for c in c3_top]

# Colour by corrected C3 classification
def c3_bar_col(cid):
    if cid in set(GROUP_A_C3) | set(CRNDE_C3): return COL_A_NOM
    if cid in set(GABA_NEIGHBOR_C3):            return '#8e44ad'
    return COL_B

bar_cols_c3 = [c3_bar_col(c) for c in c3_top]
ax.bar(range(len(c3_top)), c3_pcts, color=bar_cols_c3,
       edgecolor='white', linewidth=0.5)
ax.set_xticks(range(len(c3_top)))
ax.set_xticklabels(c3_top, rotation=45, ha='right', fontsize=8)
for lbl, cid in zip(ax.get_xticklabels(), c3_top):
    if cid in {'C3-92','C3-84','C3-95'}:
        lbl.set_fontweight('bold')
for i, cid in enumerate(c3_top):
    if cid in magma_c3.index and magma_c3.loc[cid,'bonf_sig']:
        ax.text(i, c3_pcts[i]+0.2, '***', ha='center', fontsize=8, color='black')
ax.set_ylabel('% of mouse Lmx1a cells mapping to C3 cluster')
ax.set_title('Human HYPOMAP C3 cluster assignment\nof mouse Lmx1a neurons',
             fontweight='bold')
ax.legend(handles=[
    Patch(fc=COL_A_NOM,   label='Group A LMX1A+  /  CRNDE'),
    Patch(fc='#8e44ad',   label='GABA-neighbour'),
    Patch(fc=COL_B,       label='Group B LMX1A−'),
], fontsize=8, loc='upper right')
plt.tight_layout()
plt.savefig(f'{OUTDIR}/Fig_E_mouse_human_mapping.pdf')
plt.savefig(f'{OUTDIR}/Fig_E_mouse_human_mapping.png')
print("  Saved Fig E")

# ── FIGURE F: Drug target genes ───────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(5.5, 3.5))
drug_genes  = ['CALCR','GIPR','PCSK1','BSN','BDNF','GLP1R','LEPR','MC4R','LMX1B','FOXP2']
drug_labels = ['CALCR (amylin)','GIPR (tirzepatide)','PCSK1 (prohormone)',
               'BSN (synaptic)','BDNF (neurotrophin)','GLP1R (semaglutide)',
               'LEPR (leptin)','MC4R (setmelanotide)','LMX1B (TF)','FOXP2 (TF)']
DRUG_APPROVED = {'CALCR','GIPR','PCSK1','GLP1R','MC4R','LEPR'}
drug_data = []
for gene, label in zip(drug_genes, drug_labels):
    hits = eff_nom[eff_nom['gene_symbol']==gene]
    drug_data.append({'gene':gene,'label':label,
                      'n_nom': hits['cluster_id'].nunique() if not hits.empty else 0,
                      'mean_esmu': hits['esmu'].mean() if not hits.empty else 0,
                      'approved': gene in DRUG_APPROVED})
df_drug = pd.DataFrame(drug_data).sort_values('n_nom', ascending=True)
colors_d = ['#8B0000' if r['approved'] else '#4878CF' for _, r in df_drug.iterrows()]
ax.barh(range(len(df_drug)), df_drug['n_nom'],
        color=colors_d, edgecolor='white', linewidth=0.5, height=0.65)
ax.set_yticks(range(len(df_drug)))
ax.set_yticklabels(df_drug['label'], fontsize=8)
ax.set_xlabel('# Group A LMX1A+ nominal clusters with gene as effector')
ax.set_title('Obesity-relevant genes as effectors\n'
             'in Group A LMX1A+ nominal clusters', fontweight='bold')
for i, (_, r) in enumerate(df_drug.iterrows()):
    if r['mean_esmu'] > 0:
        ax.text(r['n_nom']+0.05, i, f'ESµ={r["mean_esmu"]:.2f}',
                va='center', fontsize=7, color='grey')
ax.legend(handles=[Patch(fc='#8B0000', label='Approved/known drug target'),
                   Patch(fc='#4878CF', label='Relevant TF or pathway gene')],
          fontsize=7.5, loc='lower right')
plt.tight_layout()
plt.savefig(f'{OUTDIR}/Fig_F_drug_targets.pdf')
plt.savefig(f'{OUTDIR}/Fig_F_drug_targets.png')
print("  Saved Fig F")

# ── Summary table ─────────────────────────────────────────────────────────────
print("\n" + "="*85)
print("SUMMARY: C4 Lmx1a clusters — LDSC BMI evidence")
print("="*85)
print(f"{'Cluster':<10} {'Group':>6} {'LMX1A%':>7} {'p_LDSC':>10} "
      f"{'nom':>5} {'bonf':>5}  exome_genes")
print("-"*70)
for cid in c4_order_sankey:
    r = df_c4[df_c4['cluster_id']==cid]
    if r.empty: continue
    r = r.iloc[0]
    grp  = 'A' if r['is_grp_a'] else 'B'
    nom  = '✓' if r['nom_sig']  else ''
    bonf = '✓' if r['bonf_sig'] else ''
    spec = f"{spec_map.get(cid,0):.0f}%"
    pv   = f"{r['pvalue']:.3e}"
    eg   = str(r['exome_genes'])[:20] if r['has_exome'] else ''
    print(f"  {cid:<10} {grp:>6} {spec:>7} {pv:>10} {nom:>5} {bonf:>5}  {eg}")

print(f"\nGroup A nominal LDSC ({len(grp_a_nom)} clusters): {grp_a_nom}")
print(f"\nAll figures saved to {OUTDIR}/")


# In[11]:


# ── Panel C standalone: ESµ heatmap with explicit filtering ───────────────────
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

EFFECTOR_PATH = "genetics_effector.txt"
OUTDIR        = "figures_hypomap_corrected"

# ── Filtering toggles ─────────────────────────────────────────────────────────
ESM_MIN      = 0.0
MAGMA_P_MAX  = 0.05
MIN_CLUSTERS = 3
MAX_GENES    = 30      # None = show all

CLUSTERS = grp_a_nom  # ['C4-224', 'C4-226', 'C4-231', 'C4-248']

# ── 1. Load and filter ────────────────────────────────────────────────────────
eff = pd.read_csv(EFFECTOR_PATH, sep='\t')
eff['cluster_id'] = eff['annotation'].str.replace(r'^C4_', 'C4-', regex=True)

eff4 = eff[
    eff['cluster_id'].isin(CLUSTERS) &
    (eff['esmu']            >= ESM_MIN) &
    (eff['magma_gene_pval'] <  MAGMA_P_MAX)
].copy()

gene_stats = (eff4.groupby('gene_symbol')
              .agg(n_clusters  =('cluster_id', 'nunique'),
                   max_esmu    =('esmu', 'max'),
                   min_magma_p =('magma_gene_pval', 'min'))
              .reset_index()
              .query('n_clusters >= @MIN_CLUSTERS')
              .sort_values(['n_clusters', 'min_magma_p'],
                           ascending=[False, True])
              .reset_index(drop=True))

if MAX_GENES:
    gene_stats = gene_stats.head(MAX_GENES)

genes = gene_stats['gene_symbol'].tolist()
print(f"Filter: ESµ≥{ESM_MIN}  MAGMA p<{MAGMA_P_MAX}  "
      f"n_clusters≥{MIN_CLUSTERS}  → {len(genes)} genes")

# ── 2. Build ESµ matrix ───────────────────────────────────────────────────────
mat = (eff4[eff4['gene_symbol'].isin(genes)]
       .pivot_table(index='gene_symbol', columns='cluster_id',
                    values='esmu', aggfunc='max')
       .reindex(index=genes, columns=CLUSTERS))

# ── 3. Plot ───────────────────────────────────────────────────────────────────
fig_h = max(4, len(genes) * 0.28 + 1.5)
fig, ax = plt.subplots(figsize=(5, fig_h))

col_labels = [
    f"{c}\n{ldsc_c4.loc[c,'annotation_full'][4:].replace('_',' ')[:22]}\n"
    f"p={ldsc_c4.loc[c,'pvalue']:.1e}"
    for c in CLUSTERS
]

sns.heatmap(mat, ax=ax,
            cmap='viridis',
            vmin=0, vmax=1.0,
            linewidths=0.4, linecolor='lightgrey',
            mask=mat.isna(),
            xticklabels=col_labels,
            yticklabels=True,
            cbar_kws={'label': 'ESµ specificity', 'shrink': 0.5})

ax.set_xticklabels(ax.get_xticklabels(), rotation=35, ha='right', fontsize=8)
ax.set_yticklabels(ax.get_yticklabels(), fontsize=8)
ax.set_xlabel('C4 cluster  (FDR<0.05, CELLECT-LDSC  Yengo2018)', fontsize=9)
ax.set_ylabel('BMI effector gene', fontsize=9)
ax.set_title(
    f'BMI effector genes — Group A LMX1A+ clusters\n'
    f'ESµ≥{ESM_MIN}  |  MAGMA gene p<{MAGMA_P_MAX}  |  '
    f'≥{MIN_CLUSTERS} clusters  |  n={len(genes)} genes',
    fontweight='bold', fontsize=9, loc='left', pad=8)

plt.tight_layout()
plt.savefig(f'{OUTDIR}/Fig_panelC_heatmap.pdf', bbox_inches='tight')
plt.savefig(f'{OUTDIR}/Fig_panelC_heatmap.png', bbox_inches='tight')
plt.show()
print(f"  Saved Panel C heatmap  ({len(genes)} genes)")


# In[13]:


# ── STRING network plot: effector genes for C4-226/231/248/224 ────────────────
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
import networkx as nx
import requests
import matplotlib.patheffects as mpe

EFFECTOR_PATH  = "genetics_effector.txt"
CLUSTERS       = ['C4-226','C4-231','C4-248','C4-224']
CLUSTER_COLORS = {
    'C4-248': '#e41a1c',   # red
    'C4-231': '#377eb8',   # blue
    'C4-226': '#4daf4a',   # green
    'C4-224': '#ff7f00',   # orange
}

ESM_THRESHOLD  = 0.3     # min ESmu to include
NOM_P          = 0.05    # nominal MAGMA gene p-value threshold
STRING_SCORE   = 400     # STRING interaction confidence (0-1000)
SPECIES        = 9606    # human

# ── 1. Filter effector genes ──────────────────────────────────────────────────
eff = pd.read_csv(EFFECTOR_PATH, sep='\t')
eff['cluster_id'] = eff['annotation'].str.replace(r'^C4_','C4-', regex=True)

sub = eff[eff['cluster_id'].isin(CLUSTERS)].copy()
filt = sub[(sub['esmu'] >= ESM_THRESHOLD) &
           (sub['magma_gene_pval'] < NOM_P)].copy()

# Per gene: max ESmu across clusters, min MAGMA p, cluster membership
gene_stats = (filt.groupby('gene_symbol')
              .agg(max_esmu=('esmu','max'),
                   min_magma_p=('magma_gene_pval','min'),
                   n_clusters=('cluster_id','nunique'),
                   clusters=('cluster_id', lambda x: sorted(set(x))))
              .reset_index())

genes = gene_stats['gene_symbol'].tolist()
print(f"Querying STRING with {len(genes)} genes: {genes}")

# ── 2. STRING API calls ───────────────────────────────────────────────────────
base = "https://string-db.org/api/json"

# 2a. Map identifiers → STRING internal IDs
r_map = requests.post(f"{base}/get_string_ids", data={
    'identifiers': '\r'.join(genes),
    'species':     SPECIES,
    'limit':       1,
    'caller_identity': 'lmx1a_analysis'
})
r_map.raise_for_status()
id_map  = {d['queryItem']: d['stringId'] for d in r_map.json()}
string_ids = list(id_map.values())
print(f"Mapped {len(string_ids)}/{len(genes)} genes to STRING IDs")

# 2b. Fetch network interactions
r_net = requests.post(f"{base}/network", data={
    'identifiers':    '\r'.join(string_ids),
    'species':        SPECIES,
    'required_score': STRING_SCORE,
    'caller_identity':'lmx1a_analysis'
})
r_net.raise_for_status()
interactions = r_net.json()
print(f"Interactions returned: {len(interactions)}")

# ── 3. Build NetworkX graph ───────────────────────────────────────────────────
# Reverse ID map: stringId → gene_symbol
rev_map = {v: k for k, v in id_map.items()}

G = nx.Graph()
for gene in genes:
    G.add_node(gene)

for edge in interactions:
    g1 = rev_map.get(edge['stringId_A'])
    g2 = rev_map.get(edge['stringId_B'])
    score = edge['score']
    if g1 and g2 and g1 != g2:
        G.add_edge(g1, g2, weight=score)

# Remove isolated nodes if desired (keep all for now)
print(f"Graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

# ── 4. Node attributes ────────────────────────────────────────────────────────
gs_idx = gene_stats.set_index('gene_symbol')

node_esmu   = {g: gs_idx.loc[g,'max_esmu']   if g in gs_idx.index else 0.5
               for g in G.nodes()}
node_pval   = {g: gs_idx.loc[g,'min_magma_p'] if g in gs_idx.index else 0.05
               for g in G.nodes()}
node_clust  = {g: gs_idx.loc[g,'clusters']   if g in gs_idx.index else []
               for g in G.nodes()}

# Node size: −log10(MAGMA gene p), scaled
min_sz, max_sz = 80, 600
pvals_arr = np.array([max(node_pval[g], 1e-12) for g in G.nodes()])
log_p     = -np.log10(pvals_arr)
sz_norm   = (log_p - log_p.min()) / (log_p.max() - log_p.min() + 1e-9)
node_sizes = {g: min_sz + sz_norm[i]*(max_sz-min_sz)
              for i, g in enumerate(G.nodes())}

# Node colour: max ESmu → YlOrRd
cmap_node = plt.colormaps['YlOrRd']
esmu_vals = np.array([node_esmu[g] for g in G.nodes()])
esmu_norm = (esmu_vals - ESM_THRESHOLD) / (1.0 - ESM_THRESHOLD + 1e-9)
node_colors = {g: cmap_node(esmu_norm[i]) for i, g in enumerate(G.nodes())}

# Edge width: STRING score
edge_weights = [G[u][v]['weight']/1000 * 3 for u, v in G.edges()]

# ── 5. Layout ─────────────────────────────────────────────────────────────────
# Spring layout seeded for reproducibility; weight by STRING score
pos = nx.spring_layout(G, weight='weight', seed=42, k=1.8)

# ── 6. Figure ─────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(10, 8))

# Edges
for u, v in G.edges():
    x0, y0 = pos[u]
    x1, y1 = pos[v]
    score   = G[u][v]['weight']
    lw      = max(0.5, score / 1000 * 5)   # scale by STRING score
    ax.plot([x0, x1], [y0, y1],
            color='#555555', alpha=0.55,
            linewidth=lw, zorder=1,
            solid_capstyle='round')

strong_edges = [(u, v) for u, v in G.edges() if G[u][v]['weight'] >= 700]
edge_labels  = {(u, v): f"{int(G[u][v]['weight'])}"
                for u, v in strong_edges}
for u, v in G.edges():
    if G[u][v]['weight'] >= 700:
        mx = (pos[u][0] + pos[v][0]) / 2
        my = (pos[u][1] + pos[v][1]) / 2
        ax.text(mx, my, str(int(G[u][v]['weight'])),
                ha='center', va='center', fontsize=5,
                color='#333333', zorder=2)
# Nodes — one scatter per node for individual colours/sizes
for gene in G.nodes():
    ax.scatter(*pos[gene],
               s=node_sizes[gene],
               c=[node_colors[gene]],
               edgecolors='white',
               linewidths=0.8,
               zorder=4)

# Labels
for gene, (x, y) in pos.items():
    ax.text(x, y, gene,
            ha='center', va='center',
            fontsize=7.5, fontweight='bold',
            color='white', zorder=6,
            path_effects=[
                mpe.withStroke(linewidth=2.5, foreground='black')
            ])

# ── Cluster membership: outer ring coloured by which cluster(s) ──────────────
for gene in G.nodes():
    clust_list = node_clust.get(gene, [])
    if len(clust_list) == 1:
        ring_col = CLUSTER_COLORS[clust_list[0]]
        ax.scatter(*pos[gene],
                   s=node_sizes[gene]*1.5,
                   c='none',
                   edgecolors=ring_col,
                   linewidths=2.5,
                   zorder=3)
    elif len(clust_list) > 1:
        # Multiple clusters: black thick ring
        ax.scatter(*pos[gene],
                   s=node_sizes[gene]*1.5,
                   c='none',
                   edgecolors='black',
                   linewidths=2.5,
                   zorder=3)

# ── Colourbar: ESmu ──────────────────────────────────────────────────────────
sm = plt.cm.ScalarMappable(cmap=cmap_node,
                            norm=mcolors.Normalize(vmin=ESM_THRESHOLD, vmax=1.0))
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax, shrink=0.55, pad=0.02)
cbar.set_label('Max ESµ specificity', fontsize=9)

# ── Legend ────────────────────────────────────────────────────────────────────
legend_handles = []
for cid, col in CLUSTER_COLORS.items():
    legend_handles.append(
        mpatches.Patch(fc='none', ec=col, lw=2.5, label=f'{cid} specific'))
legend_handles.append(
    mpatches.Patch(fc='none', ec='black', lw=2.5, label='Multiple clusters'))
legend_handles.append(
    mpatches.Patch(fc='none', ec='none', lw=0, label='─────────────'))

# Size legend entries
for label, p_eg in [('p=1e-3 (small)', 1e-3), ('p=1e-9 (large)', 1e-9)]:
    lp   = -np.log10(p_eg)
    lp_range = log_p.max() - log_p.min()
    sz_n = (lp - log_p.min()) / (lp_range + 1e-9)
    sz   = min_sz + sz_n*(max_sz-min_sz)
    sz   = max(sz, 1.0)   # guard against zero/negative
    legend_handles.append(
        mpatches.Circle((0,0), radius=np.sqrt(sz/np.pi),
                        fc='#cccccc', ec='white',
                        label=f'Node size: MAGMA {label}'))

ax.legend(handles=legend_handles, fontsize=7.5,
          loc='lower left', framealpha=0.9,
          title='Ring = cluster  |  Fill = ESµ  |  Size = MAGMA −log₁₀(p)',
          title_fontsize=7)

ax.set_title(
    f'STRING functional network — BMI effector genes\n'
    f'Clusters: {", ".join(CLUSTERS)}  |  '
    f'ESµ≥{ESM_THRESHOLD}  |  MAGMA gene p<{NOM_P}  |  '
    f'STRING score≥{STRING_SCORE} (low-medium confidence)',
    fontweight='bold', fontsize=11
)
ax.axis('off')
plt.tight_layout()
plt.savefig('figures_hypomap_corrected/STRING_effector_network.pdf',
            bbox_inches='tight')
plt.savefig('figures_hypomap_corrected/STRING_effector_network.png',
            bbox_inches='tight', dpi=200)
plt.show()

# ── Summary ───────────────────────────────────────────────────────────────────
print("\nGene summary:")
print(gene_stats.sort_values('max_esmu', ascending=False)
      [['gene_symbol','max_esmu','min_magma_p','n_clusters','clusters']]
      .to_string(index=False))
print(f"\nIsolated nodes (no STRING interactions): "
      f"{[g for g in G.nodes() if G.degree(g)==0]}")


# In[6]:


# ─────────────────────────────────────────────────────────────────────────────
# Panels B + C — consolidated effector-gene figure
# Lenient cutoffs: every gene appearing as effector in any of the 4 FDR clusters
# Run after the existing Fig_A / Fig_C cells in step6
# ─────────────────────────────────────────────────────────────────────────────
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
import matplotlib.patheffects as mpe
import networkx as nx
import requests
from statsmodels.stats.multitest import multipletests

# ── Config ────────────────────────────────────────────────────────────────────
UNIFIED_PATH   = "hypomap_cellect_unified.csv"
EFFECTOR_PATH  = "genetics_effector.txt"
OUTDIR         = "figures_hypomap_corrected"
GWAS           = "BMI_Yengo2018"
GROUP_A_C4     = ['C4-226','C4-231','C4-247','C4-248','C4-229',
                  'C4-223','C4-224','C4-252','C4-228','C4-230']
CLUSTER_COLORS = {                       # same palette as your STRING figure
    'C4-248': '#e41a1c',                 # red
    'C4-231': '#377eb8',                 # blue
    'C4-226': '#4daf4a',                 # green
    'C4-224': '#ff7f00',                 # orange
}

# Lenient toggles — change here if you want to tighten later
ESM_FLOOR      = 0.3       # Smu cutoff
MAGMA_P_MAX    = 1.0       # MAGMA cutoff
STRING_SCORE   = 400       # confidence
SPECIES        = 9606
MAX_GENES_HEAT = None      # cap rows in dot-heatmap; None = show all

os.makedirs(OUTDIR, exist_ok=True)
plt.rcParams.update({'font.size':9, 'axes.titlesize':11, 'axes.labelsize':10,
                     'figure.dpi':150, 'savefig.dpi':300, 'savefig.bbox':'tight',
                     'pdf.fonttype':42, 'ps.fonttype':42})

# ── 1. Load CELLECT unified CSV → recompute FDR-sig Group A clusters ──────────
unified = pd.read_csv(UNIFIED_PATH)
ldsc_c4 = (unified[(unified['gwas']==GWAS) & (unified['method']=='new_LDSC') &
                   (unified['level']=='C4') & (unified['is_neuronal'])]
           .set_index('cluster_id').copy())

_, fdr_pvals, _, _ = multipletests(ldsc_c4['pvalue'].values, method='fdr_bh')
ldsc_c4['fdr'] = fdr_pvals

CLUSTERS = [c for c in GROUP_A_C4
            if c in ldsc_c4.index and ldsc_c4.loc[c, 'fdr'] < 0.05]
# Order columns by ascending p-value (most significant left-most)
CLUSTERS = sorted(CLUSTERS, key=lambda c: ldsc_c4.loc[c, 'pvalue'])
print(f"FDR<0.05 Group A clusters: {CLUSTERS}")

# ── 2. Load effector genes → filter to the 4 clusters ────────────────────────
eff = pd.read_csv(EFFECTOR_PATH, sep='\t')
eff['cluster_id'] = eff['annotation'].str.replace(r'^C4_', 'C4-', regex=True)

eff4 = eff[eff['cluster_id'].isin(CLUSTERS) &
           (eff['esmu']            >= ESM_FLOOR) &
           (eff['magma_gene_pval'] <= MAGMA_P_MAX)].copy()

# Per-gene summary across the 4 clusters
gene_stats = (eff4.groupby('gene_symbol')
              .agg(n_clusters    =('cluster_id', 'nunique'),
                   max_esmu      =('esmu', 'max'),
                   min_magma_p   =('magma_gene_pval', 'min'),
                   clusters      =('cluster_id', lambda x: sorted(set(x))))
              .reset_index())

# ESmu matrix [gene × cluster] — NaN where gene is not an effector in that cluster
esmu_mat = (eff4.pivot_table(index='gene_symbol', columns='cluster_id',
                             values='esmu', aggfunc='max')
                .reindex(columns=CLUSTERS))

# Order rows: sharing block desc, then max ESmu desc within block
order = (gene_stats.sort_values(['n_clusters','max_esmu'],
                                ascending=[False, False])
                   ['gene_symbol'].tolist())
if MAX_GENES_HEAT:
    order = order[:MAX_GENES_HEAT]
esmu_mat   = esmu_mat.loc[order]
gene_stats = gene_stats.set_index('gene_symbol').loc[order]

print(f"  Total effector genes across 4 clusters: {len(order)}")
print(f"  Sharing distribution: "
      f"{gene_stats['n_clusters'].value_counts().sort_index(ascending=False).to_dict()}")

# ──────────────────────────────────────────────────────────────────────────────
# PANEL B — dot-heatmap (gene × cluster)
#   colour = ESmu in cluster   |   size = MAGMA −log10(p) gene-level
# ──────────────────────────────────────────────────────────────────────────────
n_genes = len(order)
fig_h   = max(6, n_genes * 0.18 + 1.5)
fig, (ax, ax_n) = plt.subplots(1, 2, figsize=(7.5, fig_h),
                                gridspec_kw={'width_ratios':[1, 0.06],
                                             'wspace':0.05})

# Size scaling for dots from MAGMA −log10(p)
log_p_gene  = -np.log10(gene_stats['min_magma_p'].clip(lower=1e-300).values)
size_min, size_max = 25, 320
sz_norm = (log_p_gene - log_p_gene.min()) / (log_p_gene.max() - log_p_gene.min() + 1e-9)
node_sizes  = size_min + sz_norm * (size_max - size_min)

cmap        = plt.colormaps['YlOrRd']
cnorm       = mcolors.Normalize(vmin=0.2, vmax=1.0)

# Place dots
for yi, gene in enumerate(order):
    g_size = node_sizes[yi]
    for xi, cl in enumerate(CLUSTERS):
        v = esmu_mat.loc[gene, cl]
        if pd.isna(v):
            ax.scatter(xi, yi, s=g_size*0.4, marker='x',
                       color='#dddddd', linewidths=0.6, zorder=2)
        else:
            ax.scatter(xi, yi, s=g_size, c=[cmap(cnorm(v))],
                       edgecolors='black', linewidths=0.4, zorder=3)

# Sharing-block separators
boundaries = []
prev = gene_stats.iloc[0]['n_clusters']
for yi in range(1, n_genes):
    cur = gene_stats.iloc[yi]['n_clusters']
    if cur != prev:
        ax.axhline(yi - 0.5, color='#bbbbbb', lw=0.5, ls='--', zorder=1)
        boundaries.append(yi)
        prev = cur

ax.set_xlim(-0.5, len(CLUSTERS)-0.5)
ax.set_ylim(-0.5, n_genes-0.5)
ax.invert_yaxis()
ax.set_xticks(range(len(CLUSTERS)))
ax.set_xticklabels([f"{c}\np={ldsc_c4.loc[c,'pvalue']:.1e}" for c in CLUSTERS],
                   fontsize=8)
# Colour x tick labels by cluster
for i, c in enumerate(CLUSTERS):
    ax.get_xticklabels()[i].set_color(CLUSTER_COLORS.get(c, 'black'))
    ax.get_xticklabels()[i].set_fontweight('bold')
ax.set_yticks(range(n_genes))
ax.set_yticklabels(order, fontsize=7)
ax.tick_params(axis='y', length=0)
ax.set_xlabel('Human C4 cluster (FDR<0.05, BMI-LDSC Yengo2018)', fontsize=9)
ax.set_title('Panel B — BMI effector genes × LMX1A+ posterior-glutamatergic clusters\n'
             '(colour = ESµ in cluster, size = gene-level MAGMA −log₁₀p, '
             '× = not an effector)', fontweight='bold', fontsize=10, loc='left', pad=10)

# Right-hand annotation: # clusters
n_cl_arr = gene_stats['n_clusters'].values
ax_n.imshow(n_cl_arr.reshape(-1,1), aspect='auto', cmap='Greys',
            vmin=0, vmax=4, interpolation='none')
ax_n.set_xticks([0]); ax_n.set_xticklabels(['n cl'], fontsize=7)
ax_n.set_yticks([])
for yi, n in enumerate(n_cl_arr):
    ax_n.text(0, yi, str(int(n)), ha='center', va='center',
              fontsize=6.5, color='white' if n>=3 else 'black')

# Colour bar for ESmu
cax = fig.add_axes([0.93, 0.55, 0.018, 0.30])
cb  = plt.colorbar(plt.cm.ScalarMappable(norm=cnorm, cmap=cmap),
                   cax=cax)
cb.set_label('ESµ specificity', fontsize=8)
cb.ax.tick_params(labelsize=7)

# Size legend (bottom-right)
size_legend_ax = fig.add_axes([0.92, 0.15, 0.06, 0.18])
size_legend_ax.axis('off')
for i, p_ex in enumerate([1e-3, 1e-7, 1e-12]):
    lp = -np.log10(p_ex)
    sn = (lp - log_p_gene.min()) / (log_p_gene.max() - log_p_gene.min() + 1e-9)
    sn = np.clip(sn, 0, 1)
    sz = size_min + sn*(size_max-size_min)
    size_legend_ax.scatter(0.2, 0.85 - i*0.30, s=sz,
                           c='#cccccc', edgecolors='black', linewidths=0.4)
    size_legend_ax.text(0.55, 0.85 - i*0.30, f'p={p_ex:.0e}',
                        va='center', fontsize=7)
size_legend_ax.set_title('MAGMA gene p', fontsize=7.5, loc='left', pad=2)
size_legend_ax.set_xlim(0,1); size_legend_ax.set_ylim(0,1)

plt.savefig(f'{OUTDIR}/Fig_panelB_dot_heatmap.pdf', bbox_inches='tight')
plt.savefig(f'{OUTDIR}/Fig_panelB_dot_heatmap.png', bbox_inches='tight')
plt.show()
print("  Saved Panel B")

# ──────────────────────────────────────────────────────────────────────────────
# PANEL C — STRING functional network on the same gene set (lenient)
# ──────────────────────────────────────────────────────────────────────────────
genes = order   # same set as Panel B for consistency

# 2a. Map symbols → STRING IDs
r_map = requests.post("https://string-db.org/api/json/get_string_ids", data={
    'identifiers': '\r'.join(genes),
    'species': SPECIES, 'limit': 1,
    'caller_identity': 'lmx1a_panelC'
})
r_map.raise_for_status()
id_map = {d['queryItem']: d['stringId'] for d in r_map.json()}
print(f"Mapped {len(id_map)}/{len(genes)} genes to STRING")

# 2b. Fetch interactions
r_net = requests.post("https://string-db.org/api/json/network", data={
    'identifiers': '\r'.join(id_map.values()),
    'species': SPECIES, 'required_score': STRING_SCORE,
    'caller_identity': 'lmx1a_panelC'
})
r_net.raise_for_status()
interactions = r_net.json()
rev_map = {v: k for k, v in id_map.items()}

G = nx.Graph()
G.add_nodes_from(genes)
for e in interactions:
    g1, g2 = rev_map.get(e['stringId_A']), rev_map.get(e['stringId_B'])
    if g1 and g2 and g1 != g2:
        G.add_edge(g1, g2, weight=e['score'])

print(f"STRING graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

# Split connected vs isolated
isolated = [g for g in G.nodes() if G.degree(g) == 0]
connected = [g for g in G.nodes() if G.degree(g) > 0]
G_conn = G.subgraph(connected).copy()
print(f"  Connected: {len(connected)}   Isolated: {len(isolated)}")

# Node attributes
node_esmu  = {g: gene_stats.loc[g,'max_esmu']    for g in G.nodes()}
node_pval  = {g: gene_stats.loc[g,'min_magma_p'] for g in G.nodes()}
node_clst  = {g: gene_stats.loc[g,'clusters']    for g in G.nodes()}

log_p_all  = -np.log10(np.array([max(node_pval[g],1e-300) for g in G.nodes()]))
sz_n       = (log_p_all - log_p_all.min()) / (log_p_all.max() - log_p_all.min() + 1e-9)
node_size  = {g: 80 + sz_n[i]*620 for i,g in enumerate(G.nodes())}


# Layout: connected component on a spring layout, isolated nodes in a tray below
pos = nx.spring_layout(G_conn, weight='weight', seed=42, k=2.0/np.sqrt(len(connected)))
# rescale connected layout to roughly [-1,1]×[-0.2,1]
xs = np.array([p[0] for p in pos.values()]); ys = np.array([p[1] for p in pos.values()])
def rescale(p):
    x = 2*(p[0]-xs.min())/(xs.max() - xs.min()+1e-9) - 1
    y = 1.2*(p[1]-ys.min())/(ys.max() - ys.min()+1e-9) - 0.1
    return (x, y)
pos = {g: rescale(p) for g,p in pos.items()}

# Place isolated nodes in a row at y = -0.7
n_iso = len(isolated)
if n_iso:
    iso_xs = np.linspace(-1, 1, n_iso) if n_iso > 1 else np.array([0.0])
    for g, x in zip(isolated, iso_xs):
        pos[g] = (x, -0.7)

# Draw
fig, ax = plt.subplots(figsize=(11, 9))

# Edges
for u, v, d in G_conn.edges(data=True):
    ax.plot([pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]],
            color='#666666', alpha=0.5, lw=max(0.4, d['weight']/1000*4),
            zorder=1, solid_capstyle='round')

# Hub detection for label emphasis
deg = dict(G.degree())
hub_thresh = np.percentile([deg[g] for g in connected], 75) if connected else 1

# Nodes
for g in G.nodes():
    cl_list = node_clst[g]
    ring_col = CLUSTER_COLORS[cl_list[0]] if len(cl_list)==1 else 'black'
    ax.scatter(*pos[g], s=node_size[g]*1.55,
               facecolors='none', edgecolors=ring_col,
               linewidths=2.0, zorder=3)
    ax.scatter(*pos[g], s=node_size[g],
               c=[cmap(cnorm(node_esmu[g]))],
               edgecolors='white', linewidths=0.6, zorder=4)

# Labels — slightly bigger for hubs
for g, (x, y) in pos.items():
    fs = 8.5 if deg.get(g,0) >= hub_thresh else 6.5
    fw = 'bold' if deg.get(g,0) >= hub_thresh else 'normal'
    ax.text(x, y, g, ha='center', va='center',
            fontsize=fs, fontweight=fw, color='white', zorder=6,
            path_effects=[mpe.withStroke(linewidth=2.0, foreground='black')])

# Tray separator + label
if n_iso:
    ax.axhline(-0.45, color='grey', lw=0.5, ls=':', zorder=0)
    ax.text(-1.05, -0.7, f'No STRING\ninteraction\nscore≥{STRING_SCORE}',
            ha='right', va='center', fontsize=7, color='grey', style='italic')

# Cluster ring legend
ring_handles = [mpatches.Patch(fc='none', ec=col, lw=2,
                               label=f'{cl} specific')
                for cl,col in CLUSTER_COLORS.items()]
ring_handles.append(mpatches.Patch(fc='none', ec='black', lw=2,
                                   label='≥2 clusters'))
leg1 = ax.legend(handles=ring_handles, loc='upper left', fontsize=7.5,
                 title='Ring colour', title_fontsize=8, framealpha=0.9)
ax.add_artist(leg1)

# Size legend
size_handles = []
for p_ex, label in [(1e-3,'p=1e−3'), (1e-7,'p=1e−7'), (1e-12,'p=1e−12')]:
    lp = -np.log10(p_ex)
    sn = np.clip((lp - log_p_all.min())/(log_p_all.max() - log_p_all.min()+1e-9), 0, 1)
    s  = 80 + sn*620
    size_handles.append(plt.scatter([],[], s=s, c='#cccccc',
                                    edgecolors='black', linewidths=0.4,
                                    label=label))
ax.legend(handles=size_handles, loc='upper right', fontsize=7.5,
          title='Node size: MAGMA gene p', title_fontsize=8,
          labelspacing=1.6, framealpha=0.9)

# Colourbar for ESmu
cax = fig.add_axes([0.92, 0.30, 0.015, 0.28])
cb  = plt.colorbar(plt.cm.ScalarMappable(norm=cnorm, cmap=cmap), cax=cax)
cb.set_label('Max ESµ across the 4 clusters', fontsize=8)
cb.ax.tick_params(labelsize=7)

ax.set_xlim(-1.15, 1.15); ax.set_ylim(-0.95, 1.15)
ax.set_title(
    f'Panel C — STRING functional network of BMI effector genes\n'
    f'{len(genes)} genes (effector in ≥1 of {", ".join(CLUSTERS)})  |  '
    f'STRING score≥{STRING_SCORE}',
    fontweight='bold', fontsize=10, loc='left', pad=10)
ax.axis('off')
plt.tight_layout()
plt.savefig(f'{OUTDIR}/Fig_panelC_string_network.pdf', bbox_inches='tight')
plt.savefig(f'{OUTDIR}/Fig_panelC_string_network.png', bbox_inches='tight')
plt.show()
print("  Saved Panel C")

# ── Quick summary printout ────────────────────────────────────────────────────
print("\n=== Effector gene summary across the 4 FDR clusters ===")
print(f"{'gene':<12} {'n_cl':>4} {'max_ESµ':>8} {'minMAGMA p':>11}  clusters")
for g in order[:25]:
    r = gene_stats.loc[g]
    print(f"  {g:<12} {int(r['n_clusters']):>4} {r['max_esmu']:>8.2f} "
          f"{r['min_magma_p']:>11.2e}  {','.join(r['clusters'])}")
if len(order) > 25:
    print(f"  ... and {len(order)-25} more")


# In[4]:


"""
Yengo_BMI scDRS score × BMI effector gene analysis — CORRECTED
===============================================================
Uses Group A effector genes (from unified CSV / genetics_effector.txt)
and corrected cluster highlights (C4-229 = strongest Bonf-sig Lmx1a hit).

Key corrections vs previous version:
  - BMI_GENES derived from Group A effectors (≥3 clusters) in correct
    mouse gene name format, not an arbitrary hardcoded list
  - HIGHLIGHT_C4 updated to C4-229/C4-248 (not C4-248/C4-411)
  - Figure 4: Group A Bonf-sig vs Group B comparison
  - Figure colour scheme consistent with corrected tier labels
"""

import os
import numpy as np
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import scanpy as sc
import scipy.sparse as sp
import warnings
warnings.filterwarnings("ignore")

# ── CONFIG ────────────────────────────────────────────────────────────────────
EFFECTOR_PATH = "genetics_effector.txt"
OUTDIR        = "figures_yengo_effectors_corrected"
os.makedirs(OUTDIR, exist_ok=True)

# Corrected cluster sets
GROUP_A_C4 = ['C4-226','C4-231','C4-247','C4-248','C4-229',
               'C4-223','C4-224','C4-252','C4-228','C4-230']
GROUP_B_C4 = ['C4-33','C4-407','C4-409','C4-411']

# Updated highlights: C4-229 = strongest Bonf-sig (p=3.1e-08, 97% specific)
#                     C4-248 = canonical LMX1A+ (p=1.3e-03, 98% specific)
HIGHLIGHT_C4  = {'C4-229', 'C4-248'}
HIGHLIGHT_C4_COLORS = {'C4-229': '#8B0000', 'C4-248': '#C44E52'}

MOUSE_COLORS = {
    '04_Crnde_Lmx1a': '#008941',
    '27_Lmx1a_Ttc6':  '#ffaa92',
    '32_Gucy2c_Ttc6': '#000035',
    '36_Ebf2_Crhbp':  '#0aa6d8',
}

sns.set_style("whitegrid")
plt.rcParams.update({'font.size':9,'axes.titlesize':11,'axes.labelsize':10,
                     'figure.dpi':150,'savefig.dpi':300,'savefig.bbox':'tight',
                     'pdf.fonttype':42,'ps.fonttype':42})

# ── 1. Derive BMI_GENES from Group A effectors ────────────────────────────────
print("Deriving BMI gene list from Group A effectors …")
eff = pd.read_csv(EFFECTOR_PATH, sep='\t')
eff['cluster_id'] = eff['annotation'].str.replace(r'^C4_','C4-',regex=True)

eff_a = eff[eff['cluster_id'].isin(GROUP_A_C4)].copy()
gene_summary_full = (eff_a.groupby('gene_symbol')
                     .agg(n_clusters=('cluster_id','nunique'),
                          mean_esmu=('esmu','mean'),
                          min_magma_p=('magma_gene_pval','min'))
                     .sort_values(['n_clusters','min_magma_p'],
                                  ascending=[False,True]))

# Keep genes in ≥3 Group A clusters — convert to mouse gene name format
# (human ELAVL4 → mouse Elavl4: first letter upper, rest lower)
top_human_genes = gene_summary_full[gene_summary_full['n_clusters']>=3].head(25)
BMI_GENES_HUMAN = top_human_genes.index.tolist()
BMI_GENES_MOUSE = [g[0].upper()+g[1:].lower() for g in BMI_GENES_HUMAN]

# Build human → mouse lookup for annotation
gene_h2m = dict(zip(BMI_GENES_HUMAN, BMI_GENES_MOUSE))
gene_m2h = {v:k for k,v in gene_h2m.items()}

print(f"  Group A effectors (≥3 clusters): {len(BMI_GENES_MOUSE)}")
print(f"  Genes: {BMI_GENES_MOUSE}")

# ── 2. Data already loaded as adata (LMX_bipca.h5ad) ─────────────────────────
print(f"\nadata: {adata.n_obs} cells × {adata.n_vars} genes")

# Expression layer
for layer in ['log1p_norm','scran_normalization',None]:
    if layer is None or layer in adata.layers:
        expr_layer = layer
        break
print(f"  Using expression layer: {expr_layer or 'X'}")

# C4/C3 short IDs
if 'C4_short' not in adata.obs.columns:
    adata.obs['C4_short'] = adata.obs['C4_HM'].str.extract(r'(C4-\d+)')
if 'C3_short' not in adata.obs.columns:
    adata.obs['C3_short'] = adata.obs['C3_HM'].str.extract(r'(C3-\d+)')

# Cell-type column
ct_col = ('subset_celltype_annotation'
          if 'subset_celltype_annotation' in adata.obs.columns
          else 'detailed_celltype_annotation')
print(f"  Cell-type column: {ct_col}")

# ── 3. Extract expression matrix ──────────────────────────────────────────────
genes_present = [g for g in BMI_GENES_MOUSE if g in adata.var_names]
genes_missing  = [g for g in BMI_GENES_MOUSE if g not in adata.var_names]
print(f"\n  Mouse effector genes found: {len(genes_present)}/{len(BMI_GENES_MOUSE)}")
if genes_missing:
    print(f"  Missing: {genes_missing}")

adata_sub = adata[:, genes_present]
X_expr    = adata_sub.layers[expr_layer] if expr_layer else adata_sub.X
if sp.issparse(X_expr):
    X_expr = X_expr.toarray()
df_expr = pd.DataFrame(X_expr, index=adata.obs_names, columns=genes_present)

# ── 4. Metadata ───────────────────────────────────────────────────────────────
obs    = adata.obs.copy()
yengo  = obs['Yengo_BMI'].astype(float)

# Group A/B labels
obs['c4_group'] = 'Other'
obs.loc[obs['C4_short'].isin(GROUP_A_C4), 'c4_group'] = 'Group A'
obs.loc[obs['C4_short'].isin(GROUP_B_C4), 'c4_group'] = 'Group B'
obs.loc[obs['C4_short'].isin(HIGHLIGHT_C4),'c4_group'] = 'Highlight'

# ── 5. Module score ───────────────────────────────────────────────────────────
print("\n[1] Computing module score …")
z_expr       = (df_expr - df_expr.mean()) / (df_expr.std() + 1e-8)
module_score = z_expr.mean(axis=1)
obs['module_score'] = module_score.values
obs['Yengo_BMI']    = yengo.values

r_mod, p_mod = stats.spearmanr(yengo, module_score)
print(f"  Module score vs Yengo_BMI: ρ={r_mod:.3f}, p={p_mod:.2e}")

# ── 6. Gene-level correlations ────────────────────────────────────────────────
print("\n[2] Gene-level Spearman correlations …")
corr_records = []
for gene in genes_present:
    expr   = df_expr[gene].values
    pct_ex = (expr > 0).mean()
    if pct_ex < 0.02: continue
    r, p = stats.spearmanr(yengo, expr)
    corr_records.append({'gene_mouse':gene,
                         'gene_human': gene_m2h.get(gene, gene),
                         'rho':r, 'pvalue':p,
                         'pct_expressed':pct_ex,
                         'n_clusters_grpA': gene_summary_full.loc[
                             gene_m2h.get(gene,''), 'n_clusters']
                         if gene_m2h.get(gene,'') in gene_summary_full.index else 0})
df_corr = pd.DataFrame(corr_records).sort_values('rho', ascending=False)
bonf_thr = 0.05 / max(len(df_corr), 1)
df_corr['bonf_sig'] = df_corr['pvalue'] < bonf_thr

# ── 7. High vs low Yengo cells ────────────────────────────────────────────────
print("\n[3] High vs Low Yengo_BMI differential expression …")
q75, q25 = yengo.quantile(0.75), yengo.quantile(0.25)
high_mask, low_mask = yengo >= q75, yengo <= q25

de_records = []
for gene in genes_present:
    hi = df_expr.loc[high_mask, gene].values
    lo = df_expr.loc[low_mask,  gene].values
    st, p = stats.mannwhitneyu(hi, lo, alternative='two-sided')
    fc = (hi.mean()+0.01)/(lo.mean()+0.01)
    de_records.append({'gene_mouse':gene,
                       'gene_human':gene_m2h.get(gene,gene),
                       'log2fc':np.log2(fc), 'pvalue':p,
                       'mean_high':hi.mean(), 'mean_low':lo.mean()})
df_de = pd.DataFrame(de_records).sort_values('log2fc', ascending=False)
df_de['pval_bonf']   = df_de['pvalue'] * len(df_de)
df_de['significant'] = df_de['pval_bonf'] < 0.05

# ── 8. Cluster-level summaries ────────────────────────────────────────────────
print("\n[4] Cluster-level analysis …")
top_c4 = obs['C4_short'].value_counts().head(15).index.tolist()
c4_summary = (obs[obs['C4_short'].isin(top_c4)]
              .groupby('C4_short')
              .agg(n_cells=('C4_short','count'),
                   mean_yengo=('Yengo_BMI','mean'),
                   median_yengo=('Yengo_BMI','median'),
                   mean_module=('module_score','mean'))
              .reset_index()
              .sort_values('mean_yengo', ascending=False))
c4_summary['is_highlight'] = c4_summary['C4_short'].isin(HIGHLIGHT_C4)
c4_summary['group'] = c4_summary['C4_short'].map(
    lambda c: 'A' if c in GROUP_A_C4 else 'B' if c in GROUP_B_C4 else '?')

ct_summary = (obs.groupby(ct_col)
              .agg(n_cells=(ct_col,'count'),
                   mean_yengo=('Yengo_BMI','mean'),
                   mean_module=('module_score','mean'))
              .reset_index()
              .sort_values('mean_yengo', ascending=False))

# ══════════════════════════════════════════════════════════════════════════════
#                        FIGURE GENERATION
# ══════════════════════════════════════════════════════════════════════════════

# ── FIGURE 1: Gene correlation bar chart ─────────────────────────────────────
fig, ax = plt.subplots(figsize=(6, 7))

df_p = df_corr.sort_values('rho').copy()
colors = ['#8B0000' if r>0 else '#4878CF' for r in df_p['rho']]

ax.barh(range(len(df_p)), df_p['rho'],
        color=colors, edgecolor='white', linewidth=0.5, height=0.72)

for i, (_, row) in enumerate(df_p.iterrows()):
    if row['bonf_sig']:
        ax.text(row['rho']+(0.003 if row['rho']>0 else -0.003),
                i, '★', va='center',
                ha='left' if row['rho']>0 else 'right',
                fontsize=8, color='black')

ax.axvline(0, color='black', linewidth=0.8)
ax.set_yticks(range(len(df_p)))
# Labels: mouse gene + human gene + n_clusters_grpA
yticklabels = [f"{r['gene_mouse']}  [{r['gene_human']}  n={int(r['n_clusters_grpA'])}]"
               for _, r in df_p.iterrows()]
ax.set_yticklabels(yticklabels, fontsize=7.5)
ax.set_xlabel("Spearman ρ with Yengo_BMI scDRS score")
ax.set_title("Group A effector gene correlations\nwith per-cell Yengo_BMI scDRS score",
             fontweight='bold')
ax.legend(handles=[
    mpatches.Patch(fc='#8B0000', label='Positive ρ'),
    mpatches.Patch(fc='#4878CF', label='Negative ρ'),
], fontsize=8, loc='lower right')
ax.text(0.98, 0.02, f'n={len(obs)} cells  |  ★=Bonf.sig',
        transform=ax.transAxes, ha='right', fontsize=7, color='grey')
plt.tight_layout()
plt.savefig(f'{OUTDIR}/Fig1_gene_correlations.pdf')
plt.savefig(f'{OUTDIR}/Fig1_gene_correlations.png')
print("\nSaved Fig1")

# ── FIGURE 2: Module score vs Yengo_BMI + cell-type violins ──────────────────
fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

# Left: hexbin scatter
ax = axes[0]
hb = ax.hexbin(yengo, module_score, gridsize=45, cmap='YlOrRd',
               mincnt=1, linewidths=0.2)
fig.colorbar(hb, ax=ax, label='Cell count')
m, b = np.polyfit(yengo, module_score, 1)
xr = np.linspace(yengo.min(), yengo.max(), 200)
ax.plot(xr, m*xr+b, color='#8B0000', lw=2)
ax.set_xlabel('Yengo_BMI scDRS score')
ax.set_ylabel('Group A effector module score\n(mean z-scored expression)')
ax.set_title(f'Module score vs scDRS\nρ={r_mod:.3f}  p={p_mod:.2e}',
             fontweight='bold')
ax.text(0.02, 0.98,
        f'Genes: {len(genes_present)} mouse orthologs\nof Group A effectors (≥3 clusters)',
        transform=ax.transAxes, va='top', fontsize=7.5, color='grey',
        bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='grey', alpha=0.8))

# Right: violins per mouse cell type
ax = axes[1]
ct_order = ct_summary[ct_col].tolist()
groups_ct = [obs.loc[obs[ct_col]==ct, 'Yengo_BMI'].values for ct in ct_order]
vp = ax.violinplot(groups_ct, positions=range(len(ct_order)),
                   showmedians=True, showextrema=False)
for body, ct in zip(vp['bodies'], ct_order):
    body.set_facecolor(MOUSE_COLORS.get(ct,'#999999')); body.set_alpha(0.75)
vp['cmedians'].set_color('black'); vp['cmedians'].set_linewidth(1.5)
ax.set_xticks(range(len(ct_order)))
ax.set_xticklabels(ct_order, rotation=35, ha='right', fontsize=7.5)
ax.set_ylabel('Yengo_BMI scDRS score')
ax.set_title('Yengo_BMI score by\nmouse cell type', fontweight='bold')
for i, (ct, grp) in enumerate(zip(ct_order, groups_ct)):
    ax.text(i, ax.get_ylim()[0], f'n={len(grp)}',
            ha='center', va='bottom', fontsize=6, color='grey')
plt.tight_layout()
plt.savefig(f'{OUTDIR}/Fig2_module_scatter_celltype.pdf')
plt.savefig(f'{OUTDIR}/Fig2_module_scatter_celltype.png')
print("Saved Fig2")

# ── FIGURE 3: DE heatmap (high vs low Yengo_BMI) ─────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))

# Left: log2FC bar chart
ax = axes[0]
df_de_p = df_de.sort_values('log2fc', ascending=True).copy()
bar_c = ['#8B0000' if (r['significant'] and r['log2fc']>0) else
         '#4878CF' if (r['significant'] and r['log2fc']<0) else
         '#BBBBBB' for _, r in df_de_p.iterrows()]
ax.barh(range(len(df_de_p)), df_de_p['log2fc'],
        color=bar_c, edgecolor='white', linewidth=0.5, height=0.72)
ax.axvline(0, color='black', lw=0.8)
ax.set_yticks(range(len(df_de_p)))
ax.set_yticklabels([f"{r['gene_mouse']}  [{r['gene_human']}]"
                    for _, r in df_de_p.iterrows()], fontsize=7.5)
ax.set_xlabel('log₂(FC)  high vs low Yengo_BMI cells (Q3 vs Q1)')
ax.set_title('Effector gene upregulation\nin high-scDRS cells', fontweight='bold')
ax.legend(handles=[
    mpatches.Patch(fc='#8B0000', label='Higher in high-BMI (Bonf.sig)'),
    mpatches.Patch(fc='#4878CF', label='Lower in high-BMI (Bonf.sig)'),
    mpatches.Patch(fc='#BBBBBB', label='Not significant'),
], fontsize=7, loc='lower right')

# Right: heatmap — genes × mouse cell types
ax = axes[1]
ct_list = ct_summary[ct_col].tolist()
hmap = pd.DataFrame(index=genes_present, columns=ct_list, dtype=float)
for ct in ct_list:
    mask_ct = obs[ct_col] == ct
    hmap[ct] = df_expr.loc[mask_ct].mean(axis=0)
hmap_z = hmap.sub(hmap.mean(axis=1),axis=0).div(hmap.std(axis=1)+1e-8,axis=0)
gene_order = df_corr.sort_values('rho',ascending=False)['gene_mouse'].tolist()
hmap_z = hmap_z.loc[[g for g in gene_order if g in hmap_z.index]]
# Rename index to human names
hmap_z.index = [f"{g}  [{gene_m2h.get(g,g)}]" for g in hmap_z.index]

sns.heatmap(hmap_z, ax=ax, cmap='RdBu_r', center=0, linewidths=0.4,
            linecolor='lightgrey', vmin=-2, vmax=2,
            cbar_kws={'label':'Mean expr (z-score)','shrink':0.7})
ax.set_title('Group A effector genes\nby mouse cell type', fontweight='bold')
ax.tick_params(axis='x', rotation=40, labelsize=7.5)
ax.tick_params(axis='y', labelsize=7)
for lbl in ax.get_xticklabels():
    lbl.set_color(MOUSE_COLORS.get(lbl.get_text(),'black'))
    lbl.set_fontweight('bold')
plt.tight_layout()
plt.savefig(f'{OUTDIR}/Fig3_DE_heatmap.pdf')
plt.savefig(f'{OUTDIR}/Fig3_DE_heatmap.png')
print("Saved Fig3")

# ── FIGURE 4: C4 cluster Yengo_BMI — Group A/B comparison ───────────────────
fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

# Left: mean Yengo_BMI per C4 cluster, coloured by group
ax = axes[0]
df_c4p = c4_summary.sort_values('mean_yengo', ascending=True).copy()

def c4_col_fn(r):
    if r['C4_short'] == 'C4-229':   return '#8B0000'   # Tier 1 top
    if r['C4_short'] == 'C4-248':   return '#C44E52'   # Tier 3 canonical
    if r['group']    == 'A':        return '#E07B54'   # Other Group A
    if r['group']    == 'B':        return '#95a5a6'   # Group B
    return '#cccccc'

colors_c4 = [c4_col_fn(r) for _, r in df_c4p.iterrows()]
edges_c4  = ['black' if r['is_highlight'] else 'white' for _, r in df_c4p.iterrows()]
lws_c4    = [1.2     if r['is_highlight'] else 0.4     for _, r in df_c4p.iterrows()]

ax.barh(range(len(df_c4p)), df_c4p['mean_yengo'],
        color=colors_c4, edgecolor=edges_c4, linewidth=lws_c4, height=0.7)
ax.set_yticks(range(len(df_c4p)))
ax.set_yticklabels(df_c4p['C4_short'], fontsize=8)
for lbl, r in zip(ax.get_yticklabels(), df_c4p.itertuples()):
    if r.is_highlight:
        lbl.set_fontweight('bold'); lbl.set_color('#8B0000')
    elif r.group == 'B':
        lbl.set_color('#7f8c8d')
for i, (_, r) in enumerate(df_c4p.iterrows()):
    ax.text(r['mean_yengo']+0.002, i,
            f"n={r['n_cells']}  Grp{r['group']}",
            va='center', fontsize=6, color='grey')
ax.set_xlabel('Mean Yengo_BMI scDRS score')
ax.set_title('Mean BMI genetic enrichment\nper human C4 cluster (top 15)',
             fontweight='bold')
ax.legend(handles=[
    mpatches.Patch(fc='#8B0000', ec='black', lw=1, label='C4-229 (Tier1 Bonf-sig)'),
    mpatches.Patch(fc='#C44E52', ec='black', lw=1, label='C4-248 (Tier3 LMX1A+)'),
    mpatches.Patch(fc='#E07B54', label='Other Group A'),
    mpatches.Patch(fc='#95a5a6', label='Group B (contaminant)'),
], fontsize=7.5, loc='lower right')

# Right: violins — Group A Bonf-sig vs Group B vs no-label
ax = axes[1]
TIER1_C4_BONFSIG = ['C4-229','C4-247','C4-252','C4-223','C4-224','C4-228','C4-230']
grp_tier1 = obs.loc[obs['C4_short'].isin(TIER1_C4_BONFSIG), 'Yengo_BMI'].values
grp_tier3 = obs.loc[obs['C4_short'].isin(['C4-226','C4-231','C4-248']),
                    'Yengo_BMI'].values
grp_b     = obs.loc[obs['C4_short'].isin(GROUP_B_C4),  'Yengo_BMI'].values

groups_vp = [g for g in [grp_tier1, grp_tier3, grp_b] if len(g) > 0]
labels_vp = ['Group A\nBonf-sig\n(Tier1)',
              'Group A\nnominal\n(Tier3)',
              'Group B\ncontaminant'][:len(groups_vp)]
colors_vp = ['#8B0000','#E8A8A8','#95a5a6'][:len(groups_vp)]

vp = ax.violinplot(groups_vp, positions=range(len(groups_vp)),
                   showmedians=True, showextrema=False)
for body, col in zip(vp['bodies'], colors_vp):
    body.set_facecolor(col); body.set_alpha(0.8)
vp['cmedians'].set_color('black'); vp['cmedians'].set_linewidth(1.5)

# Pairwise MWU tests
if len(groups_vp) >= 2:
    mwu_t1_b, p_t1_b = stats.mannwhitneyu(grp_tier1, grp_b,
                                           alternative='two-sided')
    ax.text(0.5, ax.get_ylim()[1]*0.98,
            f'Tier1 vs B: p={p_t1_b:.2e}',
            ha='center', fontsize=7.5, color='grey',
            transform=ax.get_xaxis_transform())

ax.set_xticks(range(len(labels_vp)))
ax.set_xticklabels(labels_vp, fontsize=9)
ax.set_ylabel('Yengo_BMI scDRS score')
ax.set_title('Yengo_BMI: Group A Bonf-sig\nvs Group A nominal vs Group B',
             fontweight='bold')
for i, grp in enumerate(groups_vp):
    ax.text(i, ax.get_ylim()[0], f'n={len(grp)}',
            ha='center', va='bottom', fontsize=7.5, color='grey')
plt.tight_layout()
plt.savefig(f'{OUTDIR}/Fig4_C4_yengo.pdf')
plt.savefig(f'{OUTDIR}/Fig4_C4_yengo.png')
print("Saved Fig4")

# ── FIGURE 5: Top 6 gene distributions ───────────────────────────────────────
fig, axes = plt.subplots(2, 3, figsize=(10, 7))
axes = axes.flatten()
top6 = df_corr.sort_values('rho', ascending=False).head(6)['gene_mouse'].tolist()
bins_all = np.linspace(0, df_expr[top6].values.max(), 40)

for i, gene in enumerate(top6):
    ax = axes[i]
    hi_v = df_expr.loc[high_mask, gene].values
    lo_v = df_expr.loc[low_mask,  gene].values
    ax.hist(lo_v, bins=bins_all, color='#4878CF', alpha=0.6, density=True,
            label=f'Low BMI Q1 (n={low_mask.sum()})', edgecolor='white',
            linewidth=0.3)
    ax.hist(hi_v, bins=bins_all, color='#8B0000', alpha=0.6, density=True,
            label=f'High BMI Q3 (n={high_mask.sum()})', edgecolor='white',
            linewidth=0.3)
    r_g = df_corr.loc[df_corr['gene_mouse']==gene,'rho'].values[0]
    p_g = df_corr.loc[df_corr['gene_mouse']==gene,'pvalue'].values[0]
    n_cl = int(df_corr.loc[df_corr['gene_mouse']==gene,'n_clusters_grpA'].values[0])
    ax.set_title(f"{gene}  [{gene_m2h.get(gene,gene)}]\n"
                 f"ρ={r_g:.3f}  p={p_g:.2e}  n_clusters={n_cl}",
                 fontweight='bold', fontsize=8.5)
    ax.set_xlabel('log-normalised expression')
    ax.set_ylabel('Density')
    if i == 0:
        ax.legend(fontsize=7)

plt.suptitle('Top Group A effector genes:\nhigh vs low Yengo_BMI cells',
             fontweight='bold', y=1.01)
plt.tight_layout()
plt.savefig(f'{OUTDIR}/Fig5_top_gene_distributions.pdf')
plt.savefig(f'{OUTDIR}/Fig5_top_gene_distributions.png')
print("Saved Fig5")

# ── Summary tables ────────────────────────────────────────────────────────────
print("\n" + "="*80)
print("GENE-LEVEL CORRELATIONS (Yengo_BMI scDRS vs Group A effector expression)")
print("="*80)
print(f"Bonferroni threshold: {bonf_thr:.2e}  ({len(df_corr)} genes tested)")
print("-"*80)
print(f"{'Mouse':12} {'Human':12} {'ρ':>8} {'p-value':>12} "
      f"{'Bonf':>6} {'%exp':>6} {'n_clust':>7}")
print("-"*80)
for _, r in df_corr.iterrows():
    sig = '***' if r['bonf_sig'] else ''
    print(f"{r['gene_mouse']:<12} {r['gene_human']:<12} {r['rho']:>+.4f} "
          f"{r['pvalue']:>12.3e} {sig:>6} "
          f"{r['pct_expressed']*100:>5.1f}% {int(r['n_clusters_grpA']):>7}")

print(f"\nModule score overall: ρ={r_mod:.4f}  p={p_mod:.3e}")
print(f"\nAll figures saved to {OUTDIR}/")


# In[10]:


# ── Distribution of C4_HM label prediction scores (fixed) ───────────────────
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# ── 1. Find the prediction-score column ─────────────────────────────────────
score_candidates = [c for c in adata.obs.columns
                    if ('C4' in c or 'c4' in c) and 'score' in c.lower()]
if not score_candidates:
    score_candidates = [c for c in adata.obs.columns if 'score' in c.lower()]

print("Candidate prediction-score columns found:", score_candidates)
SCORE_COL = score_candidates[0]   # adjust index if needed
print(f"Using: '{SCORE_COL}'")

obs = adata.obs.copy()
obs[SCORE_COL] = pd.to_numeric(obs[SCORE_COL], errors='coerce')

# ── 2. Extract short cluster ID  e.g. "C4-248 Post-1 ..." → "C4-248" ────────
obs['C4_short'] = obs['C4_HM'].str.extract(r'(C4-\d+)')

# ── 3. Top N by cell count, using the short ID ───────────────────────────────
TOP_N = 20
top_ids = obs['C4_short'].value_counts().head(TOP_N).index.tolist()
obs_top = obs[obs['C4_short'].isin(top_ids)].copy()

# Order by median prediction score (descending)
order = (obs_top.groupby('C4_short')[SCORE_COL]
                .median()
                .sort_values(ascending=False)
                .index.tolist())

fig, axes = plt.subplots(1, 2, figsize=(18, 6))

# ── Left: violin per C4 cluster ──────────────────────────────────────────────
ax = axes[0]
sns.violinplot(
    data=obs_top,
    x='C4_short', y=SCORE_COL,
    order=order,
    ax=ax,
    palette='viridis',
    inner='box',
    linewidth=0.6,
    cut=0,
)
# Annotate n cells per cluster
for i, cid in enumerate(order):
    n = (obs_top['C4_short'] == cid).sum()
    ax.text(i, ax.get_ylim()[0] - 0.02, f'n={n}',
            ha='center', va='top', fontsize=6.5, color='grey')

ax.set_xticklabels(order, rotation=45, ha='right', fontsize=8)
ax.set_xlabel('Human Hypomap C4 cluster', fontsize=10)
ax.set_ylabel('Label prediction score', fontsize=10)
ax.set_title(f'C4_HM label transfer prediction score\nLmx1a neuron subset',
             fontweight='bold')

# ── Right: overall histogram ─────────────────────────────────────────────────
ax = axes[1]
ax.hist(obs[SCORE_COL].dropna(), bins=60, color='steelblue',
        edgecolor='white', linewidth=0.4)
med = obs[SCORE_COL].median()
q25 = obs[SCORE_COL].quantile(0.25)
q75 = obs[SCORE_COL].quantile(0.75)
ax.axvline(med, color='firebrick', linestyle='--', lw=1.5,
           label=f'Median = {med:.3f}')
ax.axvline(q25, color='#aaa', linestyle=':', lw=1, label=f'Q25 = {q25:.3f}')
ax.axvline(q75, color='#aaa', linestyle=':', lw=1, label=f'Q75 = {q75:.3f}')
ax.set_xlabel('Label prediction score', fontsize=10)
ax.set_ylabel('Number of cells', fontsize=10)
ax.set_title('Overall distribution of C4_HM\nlabel prediction scores', fontweight='bold')
ax.legend(fontsize=8)

plt.tight_layout()
plt.savefig('LMX_C4_prediction_score_distribution.pdf', bbox_inches='tight')
plt.savefig('LMX_C4_prediction_score_distribution.png', bbox_inches='tight')
plt.show()

print(f"\nPrediction score summary (all cells):")
print(obs[SCORE_COL].describe().round(4))
print(f"\nMedian prediction score per C4 cluster (top {TOP_N}):")
print(obs_top.groupby('C4_short')[SCORE_COL]
             .agg(['median','mean','count'])
             .sort_values('median', ascending=False).round(4))


# In[12]:


# ── C4_HM prediction score × mouse Lmx1a population ────────────────────────
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

MOUSE_CELLTYPES = ['04_Crnde_Lmx1a', '27_Lmx1a_Ttc6', '32_Gucy2c_Ttc6', '36_Ebf2_Crhbp']
MOUSE_COLORS = {
    '04_Crnde_Lmx1a': '#008941',
    '27_Lmx1a_Ttc6':  '#ffaa92',
    '32_Gucy2c_Ttc6': '#000035',
    '36_Ebf2_Crhbp':  '#0aa6d8',
}
HIGHLIGHT_C4 = {'C4-248', 'C4-411'}

# ── Build working obs ────────────────────────────────────────────────────────
obs2 = adata.obs.copy()
obs2[SCORE_COL]   = pd.to_numeric(obs2[SCORE_COL], errors='coerce')
obs2['C4_short']  = obs2['C4_HM'].str.extract(r'(C4-\d+)')
obs2['mouse_pop'] = obs2['detailed_celltype_annotation'].astype(str)

# Keep only the 4 populations of interest
obs2 = obs2[obs2['mouse_pop'].isin(MOUSE_CELLTYPES)].copy()

# Top C4 clusters overall (across all 4 populations)
TOP_N = 15
top_c4 = obs2['C4_short'].value_counts().head(TOP_N).index.tolist()
obs2_top = obs2[obs2['C4_short'].isin(top_c4)].copy()

# Order C4 clusters by overall median prediction score
c4_order = (obs2_top.groupby('C4_short')[SCORE_COL]
                    .median()
                    .sort_values(ascending=False)
                    .index.tolist())

# ════════════════════════════════════════════════════════════════════════════
# FIGURE 1 — Violin grid: one row per mouse population, columns = C4 clusters
# ════════════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(4, 1, figsize=(16, 14), sharex=True, sharey=True)

for ax, pop in zip(axes, MOUSE_CELLTYPES):
    sub = obs2_top[obs2_top['mouse_pop'] == pop]

    # Only plot C4 clusters that actually have cells from this population
    pop_c4_order = [c for c in c4_order if c in sub['C4_short'].values]

    palette = {c: '#8B0000' if c in HIGHLIGHT_C4 else MOUSE_COLORS[pop]
               for c in pop_c4_order}

    sns.violinplot(
        data=sub,
        x='C4_short', y=SCORE_COL,
        order=pop_c4_order,
        ax=ax,
        palette=palette,
        inner='box',
        linewidth=0.5,
        cut=0,
        scale='width',
    )
    # Annotate n per cluster
    for i, cid in enumerate(pop_c4_order):
        n = (sub['C4_short'] == cid).sum()
        ax.text(i, ax.get_ylim()[0], f'n={n}',
                ha='center', va='bottom', fontsize=6, color='grey')

    ax.set_ylabel(f'{pop}\nprediction score', fontsize=8, fontweight='bold',
                  color=MOUSE_COLORS[pop])
    ax.set_xlabel('')
    ax.set_title('')
    ax.tick_params(axis='x', labelsize=8, rotation=30)
    # Highlight the x-tick labels for C4-248 and C4-411
    for lbl in ax.get_xticklabels():
        if lbl.get_text() in HIGHLIGHT_C4:
            lbl.set_color('#8B0000')
            lbl.set_fontweight('bold')

axes[-1].set_xlabel('Human Hypomap C4 cluster', fontsize=10)
fig.suptitle('C4_HM label transfer prediction score\nstratified by mouse Lmx1a population',
             fontsize=13, fontweight='bold', y=1.01)
plt.tight_layout()
plt.savefig('LMX_C4_score_by_mouse_pop_violin.pdf', bbox_inches='tight')
plt.savefig('LMX_C4_score_by_mouse_pop_violin.png', bbox_inches='tight')
plt.show()


# ════════════════════════════════════════════════════════════════════════════
# FIGURE 2 — Heatmap: mean prediction score  (C4 cluster × mouse population)
#            + stacked bar showing % cell composition per C4 cluster
# ════════════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 2, figsize=(14, 6),
                         gridspec_kw={'width_ratios': [2, 1]})

# Left: heatmap of mean prediction score
pivot_score = obs2_top.pivot_table(
    index='C4_short', columns='mouse_pop',
    values=SCORE_COL, aggfunc='mean'
)[MOUSE_CELLTYPES].loc[c4_order]   # enforce consistent ordering

# Highlight rows
row_colors = ['#8B0000' if c in HIGHLIGHT_C4 else '#444444' for c in pivot_score.index]

ax = axes[0]
sns.heatmap(
    pivot_score,
    ax=ax,
    cmap='YlOrRd',
    linewidths=0.4,
    linecolor='lightgrey',
    annot=True,
    fmt='.2f',
    annot_kws={'size': 8},
    vmin=0, vmax=1,
    cbar_kws={'label': 'Mean prediction score', 'shrink': 0.7},
    xticklabels=True,
    yticklabels=True,
)
ax.set_xlabel('Mouse Lmx1a population', fontsize=10)
ax.set_ylabel('Human Hypomap C4 cluster', fontsize=10)
ax.set_title('Mean C4_HM prediction score\nper mouse population × C4 cluster',
             fontweight='bold')
ax.tick_params(axis='x', rotation=30, labelsize=8)

# Colour y-tick labels: highlight C4-248 and C4-411
for lbl, color in zip(ax.get_yticklabels(), row_colors):
    lbl.set_color(color)
    if lbl.get_text() in HIGHLIGHT_C4:
        lbl.set_fontweight('bold')

# Colour x-tick labels by mouse colour
for lbl in ax.get_xticklabels():
    lbl.set_color(MOUSE_COLORS.get(lbl.get_text(), 'black'))
    lbl.set_fontweight('bold')

# Right: stacked bar — % of each C4 cluster's cells from each mouse population
pivot_pct = (obs2_top.groupby(['C4_short', 'mouse_pop'])
                     .size()
                     .unstack(fill_value=0)
                     .reindex(c4_order)[MOUSE_CELLTYPES])
pivot_pct_norm = pivot_pct.div(pivot_pct.sum(axis=1), axis=0) * 100

ax = axes[1]
bottom = np.zeros(len(pivot_pct_norm))
for pop in MOUSE_CELLTYPES:
    vals = pivot_pct_norm[pop].values
    bars = ax.barh(range(len(c4_order)), vals, left=bottom,
                   color=MOUSE_COLORS[pop], label=pop, edgecolor='white', linewidth=0.4)
    bottom += vals

ax.set_yticks(range(len(c4_order)))
ax.set_yticklabels(c4_order, fontsize=8)
for lbl in ax.get_yticklabels():
    if lbl.get_text() in HIGHLIGHT_C4:
        lbl.set_color('#8B0000')
        lbl.set_fontweight('bold')

ax.set_xlim(0, 100)
ax.set_xlabel('% of cells', fontsize=10)
ax.set_title('Cell composition\nper C4 cluster', fontweight='bold')
ax.legend(title='Mouse population', fontsize=7, title_fontsize=7,
          bbox_to_anchor=(1.01, 1), loc='upper left')

plt.tight_layout()
plt.savefig('LMX_C4_score_heatmap_composition.pdf', bbox_inches='tight')
plt.savefig('LMX_C4_score_heatmap_composition.png', bbox_inches='tight')
plt.show()


# ════════════════════════════════════════════════════════════════════════════
# FIGURE 3 — Strip/box plot: prediction score per mouse population,
#            faceted by whether the assigned C4 cluster is a highlight or not
# ════════════════════════════════════════════════════════════════════════════
obs2['c4_group'] = obs2['C4_short'].apply(
    lambda x: 'C4-248 / C4-411\n(top Lmx1a)' if x in HIGHLIGHT_C4 else 'Other C4 clusters'
)

fig, ax = plt.subplots(figsize=(9, 5))
sns.boxplot(
    data=obs2,
    x='mouse_pop', y=SCORE_COL,
    hue='c4_group',
    order=MOUSE_CELLTYPES,
    hue_order=['C4-248 / C4-411\n(top Lmx1a)', 'Other C4 clusters'],
    palette={'C4-248 / C4-411\n(top Lmx1a)': '#8B0000',
             'Other C4 clusters': '#aaaaaa'},
    linewidth=0.8,
    flierprops=dict(marker='o', markersize=2, alpha=0.3),
    ax=ax,
)

ax.set_xlabel('')
ax.set_ylabel('Label prediction score', fontsize=10)
ax.set_title('Prediction score: top C4 clusters vs others\nstratified by mouse Lmx1a population',
             fontweight='bold')
ax.set_xticklabels(MOUSE_CELLTYPES, rotation=20, ha='right', fontsize=9)
for lbl, pop in zip(ax.get_xticklabels(), MOUSE_CELLTYPES):
    lbl.set_color(MOUSE_COLORS[pop])
    lbl.set_fontweight('bold')

ax.legend(title='C4 cluster group', fontsize=8, title_fontsize=8)
plt.tight_layout()
plt.savefig('LMX_C4_score_highlight_vs_other_by_pop.pdf', bbox_inches='tight')
plt.savefig('LMX_C4_score_highlight_vs_other_by_pop.png', bbox_inches='tight')
plt.show()

# ── Summary table ────────────────────────────────────────────────────────────
print("\nMean prediction score per mouse population × C4 cluster group:")
print(obs2.groupby(['mouse_pop', 'c4_group'])[SCORE_COL]
          .agg(['mean', 'median', 'count']).round(3).to_string())


# In[18]:


# ── Hypomap label hierarchy consistency check ────────────────────────────────
# For each cell: do C4 → C3 → C2 → C1 assignments follow a valid path
# in the original hypomap cluster tree?

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

EXCEL = "41586_2024_8504_MOESM3_ESM.xlsx"

# ── 1. Load the cluster tree edgelist and annotation sheet ───────────────────
df_edges = pd.read_excel(EXCEL, sheet_name="5_cluster_tree_edgelist")
df_ann   = pd.read_excel(EXCEL, sheet_name="6_cluster_tree_annotations")

print("Edgelist columns:", df_edges.columns.tolist())
print("Edgelist sample:")
print(df_edges.head(10).to_string(index=False))
print(f"\nTotal edges: {len(df_edges)}")

# ── 2. Build parent → children and child → parent maps ──────────────────────
# Edgelist likely has 'from'/'to' or 'parent'/'child' columns
# Detect column names
edge_cols = df_edges.columns.tolist()
print(f"\nEdge columns: {edge_cols}")

# Normalise to 'parent' / 'child'
if 'from' in edge_cols and 'to' in edge_cols:
    df_edges = df_edges.rename(columns={'from': 'parent', 'to': 'child'})
elif 'parent' not in edge_cols:
    # Use first two columns
    df_edges.columns = ['parent', 'child'] + list(df_edges.columns[2:])

# Build lookup dicts using cluster_name (since obs columns store names)
child_to_parent = dict(zip(df_edges['child'], df_edges['parent']))
parent_to_children = df_edges.groupby('parent')['child'].apply(list).to_dict()

print(f"\nChild→parent entries: {len(child_to_parent)}")
print(f"\nSample child→parent:")
for k, v in list(child_to_parent.items())[:5]:
    print(f"  {k}  →  {v}")

# ── 3. Check what format the obs label columns use ───────────────────────────
adata_obs = adata.obs.copy()
c_cols = ['C1_HM', 'C2_HM', 'C3_HM', 'C4_HM']

print("\nSample obs labels (first 3 cells):")
print(adata_obs[c_cols].head(3).to_string())

print("\nSample tree node names:")
print(list(child_to_parent.keys())[:5])
print(list(child_to_parent.values())[:5])


# In[19]:


# ── Full hierarchy consistency analysis (adjusted) ────────────────────────────

import re

# ── 3. Build the parent map directly from the short IDs in from/to ───────────
# No name→ID conversion needed — the tree already uses short IDs
df_edges_work = df_edges.copy()

# Rename if not already done
if 'from' in df_edges_work.columns:
    df_edges_work = df_edges_work.rename(columns={'from': 'parent', 'to': 'child'})

# parent/child columns already contain short IDs like "C3-84", "C2-30"
child_id_to_parent_id = dict(zip(df_edges_work['child'], df_edges_work['parent']))

print("Sample child→parent (short IDs):")
for k, v in list(child_id_to_parent_id.items())[4:12]:
    print(f"  {k}  →  {v}")

# ── 4. Extract short ID from obs full-name labels ─────────────────────────────
# obs: "C3-84 Post-1 GLU-1 ZNF676" → "C3-84"
def extract_id(label_str, level):
    if pd.isna(label_str) or str(label_str).strip() == '':
        return np.nan
    m = re.match(rf'({level}-\d+)', str(label_str).strip())
    return m.group(1) if m else np.nan

obs = adata.obs.copy()
for level in ['C1', 'C2', 'C3', 'C4']:
    col = f'{level}_HM'
    if col in obs.columns:
        obs[f'{level}_id'] = obs[col].apply(lambda x: extract_id(x, level))

# Sanity check extraction
print("\nExtraction check (first 5 cells):")
print(obs[['C1_id','C2_id','C3_id','C4_id']].head(5).to_string())

# Verify a known pair is in the map
print(f"\nC3-84 parent in tree: {child_id_to_parent_id.get('C3-84', 'NOT FOUND')}")
print(f"C4-223 parent in tree: {child_id_to_parent_id.get('C4-223', 'NOT FOUND')}")

# ── 5. Consistency check per cell ────────────────────────────────────────────
def is_consistent(child_id, parent_id):
    """Check if parent_id is the actual tree parent of child_id."""
    if pd.isna(child_id) or pd.isna(parent_id):
        return np.nan
    expected = child_id_to_parent_id.get(str(child_id))
    if expected is None:
        return np.nan   # child not in tree
    return str(parent_id) == str(expected)

results = []
for _, row in obs.iterrows():
    c1 = row.get('C1_id', np.nan)
    c2 = row.get('C2_id', np.nan)
    c3 = row.get('C3_id', np.nan)
    c4 = row.get('C4_id', np.nan)

    c4_c3 = is_consistent(c4, c3)
    c3_c2 = is_consistent(c3, c2)
    c2_c1 = is_consistent(c2, c1)

    assessable = [v for v in [c4_c3, c3_c2, c2_c1] if v is not np.nan and not (isinstance(v, float) and np.isnan(v))]
    n_violations    = sum(1 for v in assessable if v is False)
    fully_consistent = all(v is True for v in assessable) if assessable else np.nan

    results.append({
        'C1_id': c1, 'C2_id': c2, 'C3_id': c3, 'C4_id': c4,
        'C4_C3_ok':       c4_c3,
        'C3_C2_ok':       c3_c2,
        'C2_C1_ok':       c2_c1,
        'exp_C3_for_C4':  child_id_to_parent_id.get(str(c4), np.nan) if pd.notna(c4) else np.nan,
        'exp_C2_for_C3':  child_id_to_parent_id.get(str(c3), np.nan) if pd.notna(c3) else np.nan,
        'exp_C1_for_C2':  child_id_to_parent_id.get(str(c2), np.nan) if pd.notna(c2) else np.nan,
        'n_violations':   n_violations,
        'fully_consistent': fully_consistent,
        'mouse_pop':  row.get('detailed_celltype_annotation', 'unknown'),
        'subset_ann': row.get('subset_celltype_annotation',  'unknown'),
    })

df_res = pd.DataFrame(results, index=obs.index)

# ── 6. Summary statistics ─────────────────────────────────────────────────────
total = len(df_res)

print("\n" + "=" * 65)
print("HIERARCHY CONSISTENCY SUMMARY")
print("=" * 65)

for pair, col in [('C4→C3', 'C4_C3_ok'),
                  ('C3→C2', 'C3_C2_ok'),
                  ('C2→C1', 'C2_C1_ok')]:
    n_ok  = (df_res[col] == True).sum()
    n_bad = (df_res[col] == False).sum()
    n_na  = df_res[col].isna().sum()
    pct   = n_ok / (n_ok + n_bad) * 100 if (n_ok + n_bad) > 0 else 0
    print(f"\n  {pair}:  consistent={n_ok:>5}  violated={n_bad:>5}  "
          f"unassessable={n_na:>5}  ({pct:.1f}% consistent)")

n_full_ok = df_res['fully_consistent'].sum()
n_any_bad = (df_res['n_violations'] > 0).sum()
print(f"\n  Fully consistent (all 3 pairs): "
      f"{n_full_ok:.0f} / {total} cells ({n_full_ok/total*100:.1f}%)")
print(f"  At least 1 violation:           "
      f"{n_any_bad} / {total} cells ({n_any_bad/total*100:.1f}%)")

# Quick check: are violations real or extraction artefacts?
viol_sample = df_res[df_res['C4_C3_ok'] == False][
    ['C4_id','C3_id','exp_C3_for_C4']].head(10)
print(f"\nSample C4→C3 violations (C4 assigned, C3 assigned, C3 expected):")
print(viol_sample.to_string(index=False))

# ── 7. Figures ────────────────────────────────────────────────────────────────
MOUSE_COLORS = {
    '04_Crnde_Lmx1a': '#008941',
    '27_Lmx1a_Ttc6':  '#ffaa92',
    '32_Gucy2c_Ttc6': '#000035',
    '36_Ebf2_Crhbp':  '#0aa6d8',
}

fig = plt.figure(figsize=(18, 14))
gs  = fig.add_gridspec(2, 3, hspace=0.45, wspace=0.4)

# Panel A: overall consistency per pair
ax = fig.add_subplot(gs[0, 0])
pairs = ['C4→C3', 'C3→C2', 'C2→C1']
ok_pct, vl_pct, na_pct = [], [], []
for col in ['C4_C3_ok', 'C3_C2_ok', 'C2_C1_ok']:
    n_ok  = (df_res[col] == True).sum()
    n_bad = (df_res[col] == False).sum()
    n_na  = df_res[col].isna().sum()
    ok_pct.append(n_ok  / total * 100)
    vl_pct.append(n_bad / total * 100)
    na_pct.append(n_na  / total * 100)

x = np.arange(len(pairs))
ax.bar(x, ok_pct, color='#2ecc71', label='Consistent',    edgecolor='white')
ax.bar(x, vl_pct, bottom=ok_pct,
       color='#e74c3c', label='Violated',    edgecolor='white')
ax.bar(x, na_pct, bottom=np.array(ok_pct)+np.array(vl_pct),
       color='#bdc3c7', label='Unassessable', edgecolor='white')
ax.set_xticks(x)
ax.set_xticklabels(pairs, fontsize=10)
ax.set_ylabel('% of cells')
ax.set_ylim(0, 105)
ax.set_title('Label hierarchy consistency\nper parent-child pair', fontweight='bold')
ax.legend(fontsize=7.5)
for xi, (ok, vl) in enumerate(zip(ok_pct, vl_pct)):
    if vl > 0.5:
        ax.text(xi, ok + vl/2, f'{vl:.1f}%', ha='center',
                va='center', fontsize=9, color='white', fontweight='bold')

# Panel B: violations by mouse population
ax = fig.add_subplot(gs[0, 1])
pop_data = (df_res.groupby('mouse_pop')
            .agg(n_cells          =('n_violations','count'),
                 pct_any_violation=('n_violations',
                                    lambda x: (x>0).sum()/len(x)*100))
            .reset_index()
            .sort_values('pct_any_violation', ascending=True))
colors_pop = [MOUSE_COLORS.get(p, '#999999') for p in pop_data['mouse_pop']]
ax.barh(range(len(pop_data)), pop_data['pct_any_violation'],
        color=colors_pop, edgecolor='white', linewidth=0.5, height=0.6)
ax.set_yticks(range(len(pop_data)))
ax.set_yticklabels(pop_data['mouse_pop'], fontsize=8)
for i, (_, r) in enumerate(pop_data.iterrows()):
    ax.text(r['pct_any_violation']+0.3, i, f"n={r['n_cells']:.0f}",
            va='center', fontsize=7, color='grey')
ax.set_xlabel('% cells with ≥1 hierarchy violation')
ax.set_title('Violations by mouse\nLmx1a population', fontweight='bold')

# Panel C: violation rate for top 20 C4 clusters
ax = fig.add_subplot(gs[0, 2])
c4_sum = (df_res.groupby('C4_id')
          .agg(n_cells    =('C4_C3_ok','count'),
               n_viol_c4c3=('C4_C3_ok', lambda x: (x==False).sum()),
               n_viol_c3c2=('C3_C2_ok', lambda x: (x==False).sum()),
               n_viol_c2c1=('C2_C1_ok', lambda x: (x==False).sum()))
          .reset_index())
c4_sum['pct_viol'] = (
    (c4_sum['n_viol_c4c3'] + c4_sum['n_viol_c3c2'] + c4_sum['n_viol_c2c1']) /
    (c4_sum['n_cells'] * 3) * 100)
top20_c4 = c4_sum.nlargest(20, 'n_cells').sort_values('pct_viol', ascending=True)
hl_cols  = ['#8B0000' if c in {'C4-248','C4-411'} else '#4878CF'
            for c in top20_c4['C4_id']]
ax.barh(range(len(top20_c4)), top20_c4['pct_viol'],
        color=hl_cols, edgecolor='white', linewidth=0.4, height=0.72)
ax.set_yticks(range(len(top20_c4)))
ax.set_yticklabels(top20_c4['C4_id'], fontsize=7.5)
ax.set_xlabel('% parent-child checks violated')
ax.set_title('Violation rate\ntop 20 C4 clusters by cell count', fontweight='bold')
ax.legend(handles=[mpatches.Patch(fc='#8B0000', label='C4-248 / C4-411'),
                   mpatches.Patch(fc='#4878CF', label='Other C4 clusters')],
          fontsize=7)

# Panel D: C4×C3 co-occurrence heatmap
ax = fig.add_subplot(gs[1, :2])
top15_c4 = c4_sum.nlargest(15, 'n_cells')['C4_id'].tolist()
sub      = df_res[df_res['C4_id'].isin(top15_c4)].copy()

pivot_ok  = (sub[sub['C4_C3_ok']==True]
             .groupby(['C4_id','C3_id']).size().unstack(fill_value=0))
pivot_bad = (sub[sub['C4_C3_ok']==False]
             .groupby(['C4_id','C3_id']).size().unstack(fill_value=0))
all_c3    = sorted(set(pivot_ok.columns) | set(pivot_bad.columns))
pivot_ok  = pivot_ok.reindex(columns=all_c3, fill_value=0)
pivot_bad = pivot_bad.reindex(columns=all_c3, fill_value=0)
total_pair= pivot_ok + pivot_bad
frac_ok   = pivot_ok.div(total_pair.replace(0, np.nan))
frac_ok   = frac_ok.loc[:, total_pair.sum(axis=0) > 0]
frac_ok   = frac_ok.reindex(top15_c4).fillna(np.nan)

sns.heatmap(frac_ok, ax=ax, cmap='RdYlGn', vmin=0, vmax=1,
            linewidths=0.3, linecolor='lightgrey',
            cbar_kws={'label': 'Fraction consistent', 'shrink': 0.6},
            mask=frac_ok.isna(), annot=True, fmt='.2f',
            annot_kws={'size': 6.5})

# Black border on expected diagonal (true parent)
for i, c4 in enumerate(top15_c4):
    exp_c3 = child_id_to_parent_id.get(c4)
    if exp_c3 and exp_c3 in frac_ok.columns:
        j = list(frac_ok.columns).index(exp_c3)
        ax.add_patch(plt.Rectangle((j, i), 1, 1, fill=False,
                                   edgecolor='black', lw=2, zorder=5))
ax.set_title('C4→C3 label consistency\n'
             '(fraction consistent per C4×C3 pair  |  □ = expected parent)',
             fontweight='bold')
ax.set_xlabel('Assigned C3 label', fontsize=9)
ax.set_ylabel('Assigned C4 label', fontsize=9)
ax.tick_params(axis='x', rotation=45, labelsize=7.5)
ax.tick_params(axis='y', rotation=0,  labelsize=7.5)

# Panel E: per-cell violation count distribution
ax = fig.add_subplot(gs[1, 2])
counts = df_res['n_violations'].value_counts().sort_index()
ax.bar(counts.index, counts.values / total * 100,
       color=['#2ecc71','#f39c12','#e74c3c','#8e44ad'][:len(counts)],
       edgecolor='white', linewidth=0.5, width=0.6)
ax.set_xticks(counts.index)
ax.set_xticklabels([f'{int(i)} violation{"s" if i!=1 else ""}'
                    for i in counts.index], fontsize=9)
ax.set_ylabel('% of cells')
ax.set_title('Per-cell hierarchy\nviolation count', fontweight='bold')
for n_viol, cnt in counts.items():
    ax.text(n_viol, cnt/total*100 + 0.5, f'{cnt/total*100:.1f}%',
            ha='center', fontsize=9, fontweight='bold')

plt.suptitle('Hypomap C1–C4 label hierarchy consistency — Lmx1a neuron subset',
             fontweight='bold', fontsize=13, y=1.01)
plt.savefig('LMX_hierarchy_consistency.pdf', bbox_inches='tight')
plt.savefig('LMX_hierarchy_consistency.png', bbox_inches='tight')
plt.show()

# Export violations
violations = df_res[df_res['n_violations'] > 0][
    ['C1_id','C2_id','C3_id','C4_id',
     'exp_C3_for_C4','exp_C2_for_C3','exp_C1_for_C2',
     'C4_C3_ok','C3_C2_ok','C2_C1_ok',
     'n_violations','mouse_pop']].copy()
print(f"\n{len(violations)} cells with violations — sample:")
print(violations.head(15).to_string(index=False))
violations.to_csv('LMX_hierarchy_violations.csv')
print("Saved LMX_hierarchy_violations.csv")


# In[5]:


# ── Plot top 2 C4 and C3 labels on the Lmx1a subset tSNE ─────────────────────

TOP_C4 = ['C4-248', 'C4-411']
TOP_C3 = ['C3-92',  'C3-84']

# Extract short IDs
adata.obs['C4_short'] = adata.obs['C4_HM'].str.extract(r'(C4-\d+)')
adata.obs['C3_short'] = adata.obs['C3_HM'].str.extract(r'(C3-\d+)')

# Create filtered columns — NaN for everything outside the top labels
adata.obs['C4_top2'] = adata.obs['C4_short'].where(
    adata.obs['C4_short'].isin(TOP_C4))

adata.obs['C3_top2'] = adata.obs['C3_short'].where(
    adata.obs['C3_short'].isin(TOP_C3))

# Assign fixed colors
adata.uns['C4_top2_colors'] = ['#FF4A46', '#1CE6FF']   # C4-248, C4-411
adata.uns['C3_top2_colors'] = ['#FFFF00', '#006FA6']   # C3-92,  C3-84

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

sc.pl.embedding(adata, basis='tsne', color='C4_top2',
                legend_loc='on data', legend_fontsize='small',
                title='Top 2 C4 Hypomap labels\n(C4-248, C4-411)',
                na_color='#eeeeee', ax=axes[0], show=False)

sc.pl.embedding(adata, basis='tsne', color='C3_top2',
                legend_loc='on data', legend_fontsize='small',
                title='Top 2 C3 Hypomap labels\n(C3-92, C3-84)',
                na_color='#eeeeee', ax=axes[1], show=False)

plt.tight_layout()
plt.savefig('LMX_top2_C4_C3_labels.pdf', bbox_inches='tight')
plt.savefig('LMX_top2_C4_C3_labels.png', bbox_inches='tight')
plt.show()


# In[ ]:




