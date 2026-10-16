#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os as os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.sparse import csr_matrix

from multiprocessing import Pool
import matplotlib.pyplot as plt
from tqdm.contrib.concurrent import thread_map,process_map
from functools import partial
import torch
from threadpoolctl import threadpool_limits

import bipca
from bipca import plotting
from bipca.utils import stabilize_matrix
from bipca.experiments import new_svd as SVD
from openTSNE import TSNE


# In[3]:


adata = sc.read_h5ad("7plex_raw_RNA_filtered.h5ad")
adata


# In[5]:


adata.layers["raw_filtered"] = adata.X
adata.layers


# In[6]:


adata.var = pd.DataFrame(index=adata.var_names)
adata


# In[7]:


# Remove specific columns from adata.var
columns_to_remove = ["Index.10x", "nCount_HTO", "nFeature_HTO", "Index.HTO", "HTO_mcl_maxID", "HTO_mcl_secondID", "HTO_mcl_margin", "HTO_mcl_classification", "nCount_spliced", "nCount_unspliced", "nFeature_spliced", "nFeature_unspliced", "nCount_SCT", "nFeature_SCT", "SCT_snn_res.0.8", "integrated_snn_res.0.8"]  
adata.obs = adata.obs.drop(columns=columns_to_remove)
adata.obs


# In[8]:


sc.pl.embedding(adata,basis="X_umap",color="seurat_clusters", legend_loc="on data", legend_fontsize="xx-small")


# In[9]:


metadata_df = pd.read_csv('7plex_original_seurat_clusters.csv', index_col=1)
metadata_df = metadata_df[["original.clusters", "louvian.clusters"]]
metadata_df


# In[10]:


adata.obs = adata.obs.join(metadata_df, how='left')
adata.obs


# In[42]:


adata.write("raw_filtered_cleaned.h5ad")


# In[11]:


#adata = sc.read_h5ad("raw_filtered_cleaned.h5ad")
adata


# In[12]:


X_sub,index = stabilize_matrix(adata.X,threshold=100)


# In[13]:


X_sub


# In[14]:


index


# In[15]:


adata = sc.AnnData(X=X_sub,obs=adata.obs.iloc[index[0],:],var = adata.var.iloc[index[1],:])


# In[16]:


n_threads = 32
torch.set_num_threads(n_threads)
with threadpool_limits(limits=n_threads):
    op = bipca.BiPCA(n_components=-1,seed=42) # get the BiPCA operator, here n_components=-1 is to do full SVD
    Z = op.fit_transform(X_sub.toarray()) # take the count data as input and run bipca
    op.get_plotting_spectrum() # get the fitting parameters


# In[17]:


op.write_to_adata(adata) # save the results back to anndata


# In[18]:


# get MP
plotting.MP_histograms_from_bipca(adata,median=False,bins=100) 
plt.savefig('bipca_MP_raw_filtered.pdf')


# In[19]:


adata.obsm['pcs'] = SVD(adata.layers['Z_biwhite'],adata.uns["bipca"]['rank'])


# In[20]:


sc.pp.neighbors(adata, use_rep = "pcs",method="gauss", n_pcs=adata.uns["bipca"]['rank'])


# In[21]:


sc.tl.leiden(adata, resolution=1, random_state=42, flavor="leidenalg", n_iterations=-1)


# In[22]:


adata.obsm['tsne'] = np.array(TSNE(n_jobs=32).fit(adata.obsm['pcs']))


# In[23]:


sc.pl.embedding(adata,basis="tsne",color="leiden", legend_loc="on data", legend_fontsize="xx-small", save="_leiden_raw_filtered.pdf")


# In[24]:


adata.write("raw_filtered_bipca.h5ad")
adata


# In[4]:


adata = sc.read_h5ad("raw_filtered_bipca.h5ad")
adata


# In[5]:


# Calculate the number of genes with non-zero expression per cell
adata.obs['n_genes'] = (adata.X > 0).sum(axis=1)
sc.pl.scatter(adata, "nFeature_RNA", "n_genes")


# In[6]:


# Calculate the total counts per cell
adata.obs['n_counts'] = adata.X.sum(axis=1)
sc.pl.scatter(adata, "nCount_RNA", "n_counts")


# In[7]:


sns.histplot(adata.obs["n_counts"], bins=100, kde=False)


# In[9]:


scales_counts = sc.pp.normalize_total(adata, target_sum=None, inplace=False)
# log1p transform
adata.layers["log1p_norm"] = sc.pp.log1p(scales_counts["X"], copy=True)
fig, axes = plt.subplots(1, 2, figsize=(10, 5))
p1 = sns.histplot(adata.obs["n_counts"], bins=100, kde=False, ax=axes[0])
axes[0].set_title("Total counts")
p2 = sns.histplot(adata.layers["log1p_norm"].sum(1), bins=100, kde=False, ax=axes[1])
axes[1].set_title("Shifted logarithm")
plt.show()


# In[23]:


adata
print(adata.X[:5, :5].toarray())


# In[14]:


metadata_df = pd.read_csv('7plex_metadata.csv', index_col=0)
metadata_df = metadata_df[["minor.cell.type", "original.cell.type", "cellex.minor", "cellex.orig"]]
metadata_df


# In[15]:


adata.obs = adata.obs.join(metadata_df, how='left')
adata.obs


# In[17]:


sc.pl.embedding(adata,basis="tsne",color="nCount_RNA",layer="Z_biwhite", color_map="turbo", save="_nCount_RNA_filtered.pdf")


# In[18]:


marker_genes = ["Rbfox3","Syp","Th","Chat","Gad1","Slc17a6","Gfap","Cx3cr1","Sox10","Olig1","Vim","Col23a1"]
sc.pl.embedding(adata,basis="tsne",color=marker_genes,layer="Z_biwhite", save="_markers_raw_filtered.pdf")


# In[20]:


sc.pl.embedding(adata,basis="tsne",color="minor.cell.type", legend_loc="on data", legend_fontsize="xx-small", save="_minorCT_filtered.pdf")


# In[21]:


sc.pl.embedding(adata,basis="tsne",color="original.cell.type", legend_loc="on data", legend_fontsize="xx-small", save="_origCT_filtered.pdf")


# In[24]:


adata.write("raw_filtered_bipca.h5ad")
adata


# In[3]:


adata = sc.read_h5ad("raw_filtered_bipca.h5ad")
adata
print(adata.X[:5, :5].toarray())


# In[4]:


adata.layers["counts"] = adata.X
adata


# In[5]:


os.environ['R_HOME'] = '/.conda/envs/pca_env/lib/R'
import logging
import anndata2ri
import rpy2.rinterface_lib.callbacks as rcb
import rpy2.robjects as ro
rcb.logger.setLevel(logging.ERROR)
ro.pandas2ri.activate()
anndata2ri.activate()

get_ipython().run_line_magic('load_ext', 'rpy2.ipython')


# In[6]:


get_ipython().run_cell_magic('R', '', 'library(scran)\nlibrary(BiocParallel)\n')


# In[7]:


# Preliminary clustering for differentiated normalisation
adata_pp = adata.copy()
sc.pp.normalize_total(adata_pp)
sc.pp.log1p(adata_pp)
sc.pp.pca(adata_pp, n_comps=15)
sc.pp.neighbors(adata_pp)
sc.tl.leiden(adata_pp, key_added="groups", flavor="igraph", n_iterations=2)


# In[8]:


from scipy.sparse import issparse
data_mat = adata_pp.X.T
# convert to CSC if possible. See https://github.com/MarioniLab/scran/issues/70
if issparse(data_mat):
    if data_mat.nnz > 2**31 - 1:
        data_mat = data_mat.tocoo()
    else:
        data_mat = data_mat.tocsc()
ro.globalenv["data_mat"] = data_mat
ro.globalenv["input_groups"] = adata_pp.obs["groups"]


# In[9]:


del adata_pp


# In[10]:


get_ipython().run_cell_magic('R', '-o size_factors', '\nsize_factors = sizeFactors(\n    computeSumFactors(\n        SingleCellExperiment(\n            list(counts=data_mat)), \n            clusters = input_groups,\n            min.mean = 0.1,\n            BPPARAM = MulticoreParam()\n    )\n)\n')


# In[18]:


adata.obs["size_factors"] = size_factors
scran = adata.X / adata.obs["size_factors"].values[:, None]
scran = scran.tocsc()
adata.layers["scran_normalization"] = csr_matrix(sc.pp.log1p(scran))


# In[21]:


fig, axes = plt.subplots(1, 2, figsize=(10, 5))
p1 = sns.histplot(adata.obs["n_counts"], bins=100, kde=False, ax=axes[0])
axes[0].set_title("Total counts")
p2 = sns.histplot(
    adata.layers["scran_normalization"].sum(1), bins=100, kde=False, ax=axes[1]
)
axes[1].set_title("log1p with Scran estimated size factors")
plt.show()


# In[22]:


adata.write("raw_filtered_bipca.h5ad")
adata


# In[2]:


adata = sc.read("raw_filtered_bipca.h5ad", cache=True)
adata


# In[5]:


analytic_pearson = sc.experimental.pp.normalize_pearson_residuals(adata, inplace=False)
adata.layers["analytic_pearson_residuals"] = csr_matrix(analytic_pearson["X"])
adata


# In[7]:


fig, axes = plt.subplots(1, 2, figsize=(10, 5))
p1 = sns.histplot(adata.obs["n_counts"], bins=100, kde=False, ax=axes[0])
axes[0].set_title("Total counts")
p2 = sns.histplot(
    adata.layers["analytic_pearson_residuals"].sum(1), bins=100, kde=False, ax=axes[1]
)
axes[1].set_title("Analytic Pearson residuals")
plt.show()


# In[9]:


adata.write("raw_filtered_bipca.h5ad")
adata


# In[3]:


print(adata.X[:5, :5].toarray())


# In[4]:


os.environ['R_HOME'] = '/.conda/envs/pca_env/lib/R'
import logging
import anndata2ri
import rpy2.rinterface_lib.callbacks as rcb
import rpy2.robjects as ro
rcb.logger.setLevel(logging.ERROR)
ro.pandas2ri.activate()
anndata2ri.activate()

get_ipython().run_line_magic('load_ext', 'rpy2.ipython')


# In[17]:


get_ipython().run_cell_magic('R', '', 'library(scry)\nlibrary(SummarizedExperiment)\n')


# In[30]:


bdata = adata.copy()
bdata.layers.clear()
ro.globalenv["adata"] = bdata


# In[31]:


get_ipython().run_cell_magic('R', '', 'sce = devianceFeatureSelection(adata, assay="X")\n')


# In[32]:


binomial_deviance = ro.r("rowData(sce)$binomial_deviance").T


# In[33]:


idx = binomial_deviance.argsort()[-4000:]
mask = np.zeros(adata.var_names.shape, dtype=bool)
mask[idx] = True

adata.var["highly_deviant"] = mask
adata.var["binomial_deviance"] = binomial_deviance


# In[34]:


sc.pp.highly_variable_genes(adata, layer="scran_normalization")


# In[35]:


ax = sns.scatterplot(
    data=adata.var, x="means", y="dispersions", hue="highly_deviant", s=5
)
ax.set_xlim(None, 1.5)
ax.set_ylim(None, 3)
plt.show()


# In[37]:


adata.write("raw_filtered_bipca.h5ad")
adata


# In[38]:


###BACKUP created here MAR6


# In[ ]:




