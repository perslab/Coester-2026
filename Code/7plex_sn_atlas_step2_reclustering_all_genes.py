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

import statsmodels.stats.multitest as smm
from scipy.stats import pearsonr, spearmanr


# In[3]:


adata = sc.read_h5ad("7plex_raw_RNA.h5ad")
adata


# In[4]:


gene= "Dbh"
if gene in adata.var_names:
    print(f"✔️ Found {gene!r} in adata.var_names (n_vars = {adata.n_vars})")
else:
    print(f"❌ {gene!r} NOT in adata.var_names (n_vars = {adata.n_vars})")


# In[8]:


adata.layers["counts"] = adata.X
adata.layers


# In[9]:


adata.var = pd.DataFrame(index=adata.var_names)
adata


# In[10]:


# Remove specific columns from adata.var
columns_to_remove = ["Index.10x", "nCount_HTO", "nFeature_HTO", "Index.HTO", "HTO_mcl_maxID", "HTO_mcl_secondID", "HTO_mcl_margin", "HTO_mcl_classification", "nCount_spliced", "nCount_unspliced", "nFeature_spliced", "nFeature_unspliced", "nCount_SCT", "nFeature_SCT", "SCT_snn_res.0.8", "integrated_snn_res.0.8"]  
adata.obs = adata.obs.drop(columns=columns_to_remove)
adata.obs


# In[11]:


sc.pl.embedding(adata,basis="X_umap",color="seurat_clusters", legend_loc="on data", legend_fontsize="xx-small")


# In[12]:


metadata_df = pd.read_csv('filtered_bipca_annotated.csv', index_col=0)
metadata_df


# In[17]:


keep_cells = list(set(adata.obs.index).intersection(set(metadata_df.index)))
print(len(keep_cells))


# In[19]:


subset_adata = adata[adata.obs_names.isin(keep_cells)].copy()
subset_adata
del adata


# In[24]:


metadata_df = metadata_df[["manual_celltype_annotation", "predicted_brain_region", "predicted_cell_type", "detailed_celltype_annotation"]]
subset_adata.obs = subset_adata.obs.join(metadata_df, how='left')
subset_adata.obs


# In[25]:


subset_adata.write("allgenes_bipca_annotated.h5ad")


# In[8]:


adata = sc.read_h5ad("allgenes_bipca_annotated.h5ad")
adata


# In[5]:


X_sub,index = stabilize_matrix(adata.X,threshold=100)


# In[6]:


X_sub


# In[7]:


index


# In[8]:


adata = sc.AnnData(X=X_sub,obs=adata.obs.iloc[index[0],:],var = adata.var.iloc[index[1],:])


# In[9]:


if gene in adata.var_names:
    print(f"✔️ Found {gene!r} in adata.var_names (n_vars = {adata.n_vars})")
else:
    print(f"❌ {gene!r} NOT in adata.var_names (n_vars = {adata.n_vars})")


# In[13]:


n_threads = 16
torch.set_num_threads(n_threads)
with threadpool_limits(limits=n_threads):
    op = bipca.BiPCA(n_components=-1,seed=42) # get the BiPCA operator, here n_components=-1 is to do full SVD
    Z = op.fit_transform(X_sub.toarray()) # take the count data as input and run bipca
    op.get_plotting_spectrum() # get the fitting parameters


# In[14]:


op.write_to_adata(adata) # save the results back to anndata


# In[15]:


# get MP
plotting.MP_histograms_from_bipca(adata,median=False,bins=100) 
plt.savefig('bipca_MP.pdf')


# In[16]:


adata.obsm['pcs'] = SVD(adata.layers['Z_biwhite'],adata.uns["bipca"]['rank'])


# In[17]:


sc.pp.neighbors(adata, use_rep = "pcs",method="gauss", n_pcs=adata.uns["bipca"]['rank'])


# In[18]:


sc.tl.leiden(adata, resolution=1, random_state=42, flavor="leidenalg", n_iterations=-1)


# In[19]:


adata.obsm['tsne'] = np.array(TSNE(n_jobs=16).fit(adata.obsm['pcs']))


# In[20]:


sc.pl.embedding(adata,basis="tsne",color="leiden", legend_loc="on data", legend_fontsize="xx-small", save="_leiden_7plex.pdf")


# In[ ]:


adata.write("allgenes_bipca_annotated.h5ad")
adata


# In[22]:


sc.pl.embedding(adata,basis="tsne",color="manual_celltype_annotation", legend_loc="on data", legend_fontsize="xx-small", save="_manual_7plex.pdf")


# In[23]:


sc.pl.embedding(adata,basis="tsne",color="detailed_celltype_annotation", legend_loc="on data", legend_fontsize="xx-small", save="_detailed_7plex.pdf")


# In[3]:


adata = sc.read_h5ad("allgenes_bipca_annotated.h5ad")
adata


# In[4]:


sns.histplot(adata.obs["nCount_RNA"], bins=100, kde=False)


# In[5]:


adata
print(adata.X[:5, :5].toarray())


# In[28]:


sc.pl.embedding(adata,basis="tsne",color="n_genes",layer="Z_biwhite", color_map="turbo", save="_7plex_nCount.pdf")


# In[6]:


adata.layers["counts"] = adata.X
adata


# In[7]:


os.environ['R_HOME'] = '/.conda/envs/pca_env/lib/R'
import logging
import anndata2ri
import rpy2.rinterface_lib.callbacks as rcb
import rpy2.robjects as ro
rcb.logger.setLevel(logging.ERROR)
ro.pandas2ri.activate()
anndata2ri.activate()

get_ipython().run_line_magic('load_ext', 'rpy2.ipython')


# In[8]:


get_ipython().run_cell_magic('R', '', 'library(scran)\nlibrary(BiocParallel)\n')


# In[10]:


from scipy.sparse import issparse
bdata = adata.copy()
data_mat = bdata.X.T
# convert to CSC if possible. See https://github.com/MarioniLab/scran/issues/70
if issparse(data_mat):
    if data_mat.nnz > 2**31 - 1:
        data_mat = data_mat.tocoo()
    else:
        data_mat = data_mat.tocsc()
ro.globalenv["data_mat"] = data_mat
ro.globalenv["input_groups"] = adata.obs["leiden"]
del bdata


# In[11]:


get_ipython().run_cell_magic('R', '-o size_factors', '\nsize_factors = sizeFactors(\n    computeSumFactors(\n        SingleCellExperiment(\n            list(counts=data_mat)), \n            clusters = input_groups,\n            min.mean = 0.1,\n            BPPARAM = MulticoreParam(16)\n    )\n)\n')


# In[17]:


adata.obs["size_factors"] = size_factors
scran = adata.X / adata.obs["size_factors"].values[:, None]
scran = scran.tocsc()
adata.layers["scran_normalization"] = csr_matrix(sc.pp.log1p(scran))


# In[18]:


fig, axes = plt.subplots(1, 2, figsize=(10, 5))
p1 = sns.histplot(adata.obs["nCount_RNA"], bins=100, kde=False, ax=axes[0])
axes[0].set_title("Total counts")
p2 = sns.histplot(
    adata.layers["scran_normalization"].sum(1), bins=100, kde=False, ax=axes[1]
)
axes[1].set_title("log1p with Scran estimated size factors")
plt.show()


# In[15]:


adata.write("allgenes_bipca_annotated.h5ad")
#adata = sc.read("allgenes_bipca_annotated.h5ad", cache=True)
adata


# In[19]:


os.environ['R_HOME'] = '/.conda/envs/pca_env/lib/R'
import logging
import anndata2ri
import rpy2.rinterface_lib.callbacks as rcb
import rpy2.robjects as ro
rcb.logger.setLevel(logging.ERROR)
ro.pandas2ri.activate()
anndata2ri.activate()

get_ipython().run_line_magic('load_ext', 'rpy2.ipython')


# In[20]:


get_ipython().run_cell_magic('R', '', 'library(scry)\nlibrary(SummarizedExperiment)\n')


# In[21]:


print(adata.X[:5, :5].toarray())


# In[22]:


bdata = adata.copy()
bdata.layers.clear()
ro.globalenv["adata"] = bdata


# In[23]:


del bdata


# In[24]:


get_ipython().run_cell_magic('R', '', 'sce = devianceFeatureSelection(adata, assay="X")\n')


# In[25]:


binomial_deviance = ro.r("rowData(sce)$binomial_deviance").T


# In[26]:


idx = binomial_deviance.argsort()[-4000:]
mask = np.zeros(adata.var_names.shape, dtype=bool)
mask[idx] = True

adata.var["highly_deviant"] = mask
adata.var["binomial_deviance"] = binomial_deviance


# In[27]:


sc.pp.highly_variable_genes(adata, layer="scran_normalization")


# In[28]:


ax = sns.scatterplot(
    data=adata.var, x="means", y="dispersions", hue="highly_deviant", s=5
)
ax.set_xlim(None, 1.5)
ax.set_ylim(None, 3)
plt.show()


# In[3]:


#adata.write("allgenes_bipca_annotated.h5ad")
adata = sc.read("allgenes_bipca_annotated.h5ad", cache=True)
adata


# In[30]:


counts = pd.DataFrame(adata.layers["scran_normalization"].toarray(),
                  index=adata.obs_names,
                  columns=adata.var_names)
print(counts)


# In[31]:


meta1 = adata.obs.manual_celltype_annotation.to_frame()
meta2 = adata.obs.detailed_celltype_annotation.to_frame()


# In[32]:


counts.to_pickle("allgenes_counts.pkl")
meta1.to_pickle("allgenes_manual.pkl")
meta2.to_pickle("allgenes_detailed.pkl")


# In[29]:


metadata_df = pd.read_csv('prioritization.csv')
metadata_dict = dict(zip(metadata_df["annotation"], metadata_df["pvalue"]))
metadata_dict


# In[20]:


adata.obs["cellect_bmi_manual"] = adata.obs["manual_celltype_annotation"].map(metadata_dict)
adata.obs["cellect_bmi_manual"] = pd.to_numeric(adata.obs["cellect_bmi_manual"], errors="coerce")

adata.obs["cellect_bmi_detailed"] = adata.obs["detailed_celltype_annotation"].map(metadata_dict)
adata.obs["cellect_bmi_detailed"] = pd.to_numeric(adata.obs["cellect_bmi_detailed"], errors="coerce")


# In[36]:


sc.pl.embedding(adata,basis="tsne",color="cellect_bmi_manual",layer="Z_biwhite", color_map="cividis_r")
sc.pl.embedding(adata,basis="tsne",color="cellect_bmi_detailed",layer="Z_biwhite", color_map="cividis_r")


# In[34]:


detailed = metadata_df["specificity_id"] == "detailed"
manual = metadata_df["specificity_id"] == "manual"
# Split the DataFrame
detailed_df = metadata_df[detailed].copy()
manual_df = metadata_df[manual].copy()
# Apply BH-correction
detailed_df["pvalue_BH"] = smm.multipletests(detailed_df["pvalue"], method="fdr_bh")[1]
manual_df["pvalue_BH"] = smm.multipletests(manual_df["pvalue"], method="fdr_bh")[1]


# In[37]:


metadata_dict = dict(zip(manual_df["annotation"], manual_df["pvalue_BH"]))
adata.obs["cellect_bmi_manual_BH"] = adata.obs["manual_celltype_annotation"].map(metadata_dict)
adata.obs["cellect_bmi_manual_BH"] = pd.to_numeric(adata.obs["cellect_bmi_manual_BH"], errors="coerce")
metadata_dict = dict(zip(detailed_df["annotation"], detailed_df["pvalue_BH"]))
adata.obs["cellect_bmi_detailed_BH"] = adata.obs["detailed_celltype_annotation"].map(metadata_dict)
adata.obs["cellect_bmi_detailed_BH"] = pd.to_numeric(adata.obs["cellect_bmi_detailed_BH"], errors="coerce")


# In[38]:


sc.pl.embedding(adata,basis="tsne",color="cellect_bmi_manual_BH",layer="Z_biwhite", color_map="cividis_r")
sc.pl.embedding(adata,basis="tsne",color="cellect_bmi_detailed_BH",layer="Z_biwhite", color_map="cividis_r")


# In[39]:


adata.obs["manual_significant_BH"] = adata.obs["cellect_bmi_manual_BH"] < 0.05
adata.obs["detailed_significant_BH"] = adata.obs["cellect_bmi_detailed_BH"] < 0.05


# In[40]:


sc.pl.embedding(adata,basis="tsne",color="manual_significant_BH",layer="Z_biwhite")
sc.pl.embedding(adata,basis="tsne",color="detailed_significant_BH",layer="Z_biwhite")


# In[41]:


#adata.write("allgenes_bipca_annotated.h5ad")
adata = sc.read("allgenes_bipca_annotated.h5ad", cache=True)
adata


# In[ ]:


adata.obs["sort"] = pd.to_numeric(adata.obs["sort"], errors="coerce")


# In[54]:


sc.pl.embedding(adata,basis="tsne",color="sort", color_map="pink_r")


# In[57]:


marker_genes = ["v147"]
sc.pl.embedding(adata,basis="tsne",color=marker_genes,layer="scran_normalization", color_map="pink_r", vmax=1)


# In[79]:


gene_name = "v147"
gene_expression = adata[:, gene_name].layers["scran_normalization"].toarray().flatten()

sort_values = adata.obs["sort"].values

pearson_corr, pearson_p = pearsonr(gene_expression, sort_values)
spearman_corr, spearman_p = spearmanr(gene_expression, sort_values)

print(f"Pearson Correlation: {pearson_corr}, p-value: {pearson_p}")
print(f"Spearman Correlation: {spearman_corr}, p-value: {spearman_p}")
plot_data = pd.DataFrame({
    'Gene Expression': gene_expression,
    'Sort Values': sort_values
})

# Plot using seaborn
sns.scatterplot(data=plot_data, x='Gene Expression', y='Sort Values')
plt.title(f'{gene_name} Expression vs Sort Values\nPearson: {pearson_corr:.2f}, Spearman: {spearman_corr:.2f}')
plt.xlabel(f'{gene_name} Expression')
plt.ylabel('Sort Values')
plt.show()


# In[ ]:




