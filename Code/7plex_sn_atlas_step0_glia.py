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
os.chdir("/projects/perslab/people/qbx224/7plex_paper/atlas/")


# In[2]:


adata = sc.read_h5ad("glia_bipca_subset.h5ad")
adata


# In[3]:


print(adata.X[:5, :5].toarray()) #Check if counts are the active layer


# In[4]:


###BIPCA workflow
X_sub, index = stabilize_matrix(adata.X, threshold=100)
# Subset the AnnData object
adata_subset = adata[index[0], index[1]].copy()
# Replace the data matrix with the filtered one
adata_subset.X = X_sub
adata_subset


# In[5]:


adata = adata_subset.copy()
del adata_subset
adata


# In[6]:


n_threads = 16
torch.set_num_threads(n_threads)
with threadpool_limits(limits=n_threads):
    op = bipca.BiPCA(n_components=-1,seed=42) # get the BiPCA operator, here n_components=-1 is to do full SVD
    Z = op.fit_transform(X_sub.toarray()) # take the count data as input and run bipca
    op.get_plotting_spectrum() # get the fitting parameters


# In[7]:


op.write_to_adata(adata) # save the results back to anndata


# In[8]:


# get MP
plotting.MP_histograms_from_bipca(adata,median=False,bins=100) 
plt.savefig('bipca_MP_glia.pdf')


# In[9]:


adata.obsm['pcs'] = SVD(adata.layers['Z_biwhite'],adata.uns["bipca"]['rank'])


# In[10]:


sc.pp.neighbors(adata, use_rep = "pcs",method="gauss", n_pcs=adata.uns["bipca"]['rank'])


# In[11]:


sc.tl.leiden(adata, resolution=1, random_state=42, flavor="leidenalg", n_iterations=-1)


# In[12]:


adata.obsm['tsne'] = np.array(TSNE(n_jobs=32).fit(adata.obsm['pcs']))


# In[13]:


sc.pl.embedding(adata,basis="tsne",color="leiden", legend_loc="on data", legend_fontsize="xx-small", save="_leiden_glia.pdf")


# In[14]:


adata.write("glia_bipca_annotated.h5ad")
adata


# In[3]:


adata = sc.read_h5ad("glia_bipca_annotated.h5ad")
adata


# In[15]:


sc.pl.embedding(adata,basis="tsne",color="n_genes",layer="Z_biwhite", color_map="turbo", save="_nCount_glia.pdf")


# In[16]:


sc.pl.embedding(adata,basis="tsne",color="minor.cell.type", legend_loc="on data", legend_fontsize="xx-small", save="_minorCT_glia.pdf")


# In[17]:


sc.pl.embedding(adata,basis="tsne",color="original.cell.type", legend_loc="on data", legend_fontsize="xx-small", save="_origCT_glia.pdf")


# In[18]:


sc.pl.embedding(adata,basis="tsne",color="manual_celltype_annotation", legend_loc="on data", legend_fontsize="xx-small", save="_celltype_glia.pdf")


# In[19]:


adata.layers["counts"] = adata.X


# In[20]:


fig, axes = plt.subplots(1, 2, figsize=(10, 5))
p1 = sns.histplot(adata.obs["n_counts"], bins=100, kde=False, ax=axes[0])
axes[0].set_title("Total counts")
p2 = sns.histplot(
    adata.layers["scran_normalization"].sum(1), bins=100, kde=False, ax=axes[1]
)
axes[1].set_title("log1p with Scran estimated size factors")
plt.show()


# In[21]:


ax = sns.scatterplot(
    data=adata.var, x="means", y="dispersions", hue="highly_deviant", s=5
)
ax.set_xlim(None, 1.5)
ax.set_ylim(None, 3)
plt.show()


# In[22]:


adata.X = adata.layers["scran_normalization"]
sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon", key_added="dea_leiden")


# In[23]:


sc.tl.dendrogram(adata, groupby="leiden", use_rep = "pcs", n_pcs=adata.uns["bipca"]['rank'])
sc.pl.rank_genes_groups_dotplot(adata, groupby="leiden", standard_scale="var", n_genes=5, key="dea_leiden", save="markers_glia.pdf")


# In[24]:


sc.tl.filter_rank_genes_groups(
    adata,
    min_in_group_fraction=0.2,
    max_out_group_fraction=0.2,
    key="dea_leiden",
    key_added="dea_leiden_filtered"
)#To preserve the original structure of adata.uns[‘rank_genes_groups’], filtered genes are set to NaN.


# In[25]:


sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby="leiden",
    standard_scale="var",
    n_genes=5,
    key="dea_leiden_filtered",
    save="markers_glia_filtered.pdf"
)


# In[26]:


df = sc.get.rank_genes_groups_df(adata, group=None, key = "dea_leiden_filtered", log2fc_min=3)
df.to_csv("dea_glia_filtered.csv")


# In[34]:


cl_annotation = adata.obs["leiden"]
cl_annotation["0"] = "Oligodendrocyte"
cl_annotation["1"] = "Astrocyte"
cl_annotation["2"] = "Oligodendrocyte"
cl_annotation["3"] = "Tanycyte"
cl_annotation["4"] = "Astrocyte"
cl_annotation["5"] = "OPC"
cl_annotation["6"] = "Astrocyte"
cl_annotation["7"] = "Microglia"
cl_annotation["8"] = "Pituitary_Rxfp1"
cl_annotation["9"] = "VLMC"
cl_annotation["10"] = "Pituicyte"
cl_annotation["11"] = "Pituitary_Pappa2"
cl_annotation["12"] = "Pituitary_Mtnr1a"
cl_annotation["13"] = "Oligodendrocyte"
cl_annotation["14"] = "Pituitary_Pax7"
cl_annotation["15"] = "Pituitary_Gnrhr"
cl_annotation["16"] = "Endothelial_cell"
cl_annotation["17"] = "Astrocyte"
cl_annotation["18"] = "Epithelial-like_cell"
cl_annotation["19"] = "VLMC"
cl_annotation["20"] = "Neuron-like"
cl_annotation["21"] = "OPC"
cl_annotation["22"] = "Ependymal_cell"
cl_annotation["23"] = "Pituitary_Crhr1"
cl_annotation["24"] = "Oligodendrocyte"
cl_annotation["25"] = "Undefined"
cl_annotation["26"] = "Schwann_cell"
cl_annotation["27"] = "Undefined"
cl_annotation["28"] = "Pericyte"
cl_annotation["29"] = "OPC"
cl_annotation["30"] = "Undefined"
cl_annotation["31"] = "Undefined"
cl_annotation["32"] = "Undefined"
cl_annotation["33"] = "OPC"
cl_annotation["34"] = "Neuron-like"

adata.obs["glia_celltype"] = adata.obs.leiden.map(cl_annotation)
sc.pl.embedding(adata,basis="tsne",color="glia_celltype", legend_loc="on data", legend_fontsize="xx-small", save="_glia_clusters.pdf")


# In[43]:


sc.tl.dendrogram(adata, groupby="leiden", use_rep = "pcs", n_pcs=adata.uns["bipca"]['rank'])
sc.pl.correlation_matrix(adata, "glia_celltype")


# In[32]:


marker_genes = ["Gnrhr"]
sc.pl.embedding(adata,basis="tsne",color=marker_genes,layer="Z_biwhite")


# In[36]:


glia = adata[adata.obs["glia_celltype"]!="Neuron-like"].copy()
sc.pl.embedding(glia,basis="tsne",color="glia_celltype", legend_loc="on data", legend_fontsize="xx-small")


# In[37]:


del glia.uns['dea_leiden_filtered']
glia.write("glia_bipca_annotated.h5ad")
glia


# In[39]:


metadata = glia.obs
metadata


# In[40]:


metadata.to_csv("glia_bipca_annotated.csv")


# In[ ]:




