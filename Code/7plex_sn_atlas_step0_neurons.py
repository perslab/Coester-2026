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


# In[2]:


adata = sc.read_h5ad("neurons_bipca_subset.h5ad")
adata


# In[3]:


print(adata.X[:5, :5].toarray())#Check if counts are the active layer


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
plt.savefig('bipca_MP_neurons.pdf')


# In[9]:


adata.obsm['pcs'] = SVD(adata.layers['Z_biwhite'],adata.uns["bipca"]['rank'])


# In[10]:


sc.pp.neighbors(adata, use_rep = "pcs",method="gauss", n_pcs=adata.uns["bipca"]['rank'])


# In[11]:


sc.tl.leiden(adata, resolution=1, random_state=42, flavor="leidenalg", n_iterations=-1)


# In[12]:


adata.obsm['tsne'] = np.array(TSNE(n_jobs=32).fit(adata.obsm['pcs']))


# In[13]:


sc.pl.embedding(adata,basis="tsne",color="leiden", legend_loc="on data", legend_fontsize="xx-small", save="_leiden_neurons.pdf")


# In[14]:


adata.write("neurons_bipca_annotated.h5ad")
adata


# In[2]:


adata = sc.read_h5ad("neurons_bipca_annotated.h5ad")
adata


# In[15]:


sc.pl.embedding(adata,basis="tsne",color="n_genes",layer="Z_biwhite", color_map="turbo", save="_nCount_neurons.pdf")


# In[16]:


marker_genes = ["Pomc", "Agrp", "Kiss1", "Sst", "Lepr", "Calcr", "Nr5a1", "Rgs16", "Oxt", "Hdc", "Th", "Chat", "Gad1", "Slc17a6"]
sc.pl.embedding(adata,basis="tsne",color=marker_genes,layer="Z_biwhite", save="_markers_neurons.pdf")


# In[17]:


sc.pl.embedding(adata,basis="tsne",color="minor.cell.type", legend_loc="on data", legend_fontsize="xx-small", save="_minorCT_neurons.pdf")


# In[18]:


sc.pl.embedding(adata,basis="tsne",color="original.cell.type", legend_loc="on data", legend_fontsize="xx-small", save="_origCT_neurons.pdf")


# In[19]:


sc.pl.embedding(adata,basis="tsne",color="manual_celltype_annotation", legend_loc="on data", legend_fontsize="xx-small", save="_celltype_neurons.pdf")


# In[20]:


adata.layers["counts"] = adata.X


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


ax = sns.scatterplot(
    data=adata.var, x="means", y="dispersions", hue="highly_deviant", s=5
)
ax.set_xlim(None, 1.5)
ax.set_ylim(None, 3)
plt.show()


# In[23]:


adata.layers["counts"] = adata.X
adata.X = adata.layers["scran_normalization"]
sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon", key_added="dea_leiden")


# In[24]:


sc.tl.dendrogram(adata, groupby="leiden", use_rep = "pcs", n_pcs=adata.uns["bipca"]['rank'])
sc.pl.rank_genes_groups_dotplot(adata, groupby="leiden", standard_scale="var", n_genes=5, key="dea_leiden", save="markers_neurons.pdf")


# In[25]:


sc.tl.filter_rank_genes_groups(
    adata,
    min_in_group_fraction=0.2,
    max_out_group_fraction=0.2,
    key="dea_leiden",
    key_added="dea_leiden_filtered"
)#To preserve the original structure of adata.uns[‘rank_genes_groups’], filtered genes are set to NaN.


# In[26]:


sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby="leiden",
    standard_scale="var",
    n_genes=5,
    key="dea_leiden_filtered",
    save="markers_neurons_filtered.pdf"
)


# In[27]:


df = sc.get.rank_genes_groups_df(adata, group=None, key = "dea_leiden_filtered", log2fc_min=1.5)
df.to_csv("dea_neurons_filtered.csv")


# In[28]:


cl_annotation = adata.obs["leiden"]
cl_annotation["0"] = "00_Hmcn2_Gpr149"
cl_annotation["1"] = "01_Undefined"
cl_annotation["2"] = "02_Fezf1_Lama1"
cl_annotation["3"] = "03_Rai14_Acvr1c"
cl_annotation["4"] = "04_Crnde_Lmx1a"
cl_annotation["5"] = "05_Qrfpr_Pdzrn3"
cl_annotation["6"] = "06_Tbx3_Npy1r"
cl_annotation["7"] = "07_Nfib_Six3"
cl_annotation["8"] = "08_Lef1_Rspo2"
cl_annotation["9"] = "09_Nfia_Ccbe1"
cl_annotation["10"] = "10_Npy_Agrp"
cl_annotation["11"] = "11_Tafa4_Kank4"
cl_annotation["12"] = "12_Ntng2_Onecut1"
cl_annotation["13"] = "13_Pth2r_Otp"
cl_annotation["14"] = "14_Lhx1_Th"
cl_annotation["15"] = "15_Onecut1_Slc32a1"
cl_annotation["16"] = "16_Gal_Sox6"
cl_annotation["17"] = "17_Dkk2_Ebf3"
cl_annotation["18"] = "18_Ust_Bcl11b"
cl_annotation["19"] = "19_Bace2_Tspan18"
cl_annotation["20"] = "20_Tmem114_Npsr1"
cl_annotation["21"] = "21_Bnc2_Glp1r"
cl_annotation["22"] = "22_Crabp1_Vgll3"
cl_annotation["23"] = "23_Ust_Nkd2"
cl_annotation["24"] = "24_Nts_Ecel1"
cl_annotation["25"] = "25_Nxph4_Ebf3"
cl_annotation["26"] = "26_Dach2_Npsr1"
cl_annotation["27"] = "27_Lmx1a_Ttc6"
cl_annotation["28"] = "28_Kiss1_Slc38a3"
cl_annotation["29"] = "29_Dnah11_Pth2r"
cl_annotation["30"] = "30_Trh_Sim1"
cl_annotation["31"] = "31_Lhx9_Npvf"
cl_annotation["32"] = "32_Gucy2c_Ttc6"
cl_annotation["33"] = "33_Chat_Ghrh"
cl_annotation["34"] = "34_Pla2r1_Sim1"
cl_annotation["35"] = "35_Lhx9_Zic1"
cl_annotation["36"] = "36_Ebf2_Crhbp"
cl_annotation["37"] = "37_Qrfp_Qrfprl"
cl_annotation["38"] = "38_Hdc_Prph"
cl_annotation["39"] = "39_Nfib_Nfia"
cl_annotation["40"] = "40_Gpr50_Tac2"
cl_annotation["41"] = "41_Otx1_Zic1"
cl_annotation["42"] = "42_Foxb1_Myo5b"
cl_annotation["43"] = "Glia-like"
cl_annotation["44"] = "Glia-like"

adata.obs["neuronal_celltype"] = adata.obs.leiden.map(cl_annotation)
sc.pl.embedding(adata,basis="tsne",color="neuronal_celltype", legend_loc="on data", legend_fontsize="xx-small", save="_neuronal_clusters.pdf")


# In[29]:


marker_genes = ["Foxb1"]
sc.pl.embedding(adata,basis="tsne",color=marker_genes,layer="Z_biwhite")


# In[30]:


sc.pl.correlation_matrix(adata, "neuronal_celltype")


# In[31]:


neurons = adata[adata.obs["neuronal_celltype"]!="Glia-like"].copy()
sc.pl.embedding(neurons,basis="tsne",color="neuronal_celltype", legend_loc="on data", legend_fontsize="xx-small")


# In[32]:


del neurons.uns['dea_leiden_filtered']
neurons.write("neurons_bipca_annotated.h5ad")
neurons


# In[33]:


metadata = neurons.obs
metadata


# In[34]:


metadata.to_csv("neurons_bipca_annotated.csv")


# In[ ]:




