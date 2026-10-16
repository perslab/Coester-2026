#!/usr/bin/env python
# coding: utf-8

# In[2]:


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


# In[2]:


adata = sc.read("raw_filtered_bipca.h5ad", cache=True)
adata
#adata.filename


# In[3]:


adata.X = adata.layers["scran_normalization"]
sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon", key_added="dea_leiden")


# In[3]:


sc.tl.dendrogram(adata, groupby="leiden", use_rep = "pcs", n_pcs=adata.uns["bipca"]['rank'])
sc.pl.rank_genes_groups_dotplot(adata, groupby="leiden", standard_scale="var", n_genes=5, key="dea_leiden", save="markers_raw.pdf")


# In[4]:


sc.tl.filter_rank_genes_groups(
    adata,
    min_in_group_fraction=0.2,
    max_out_group_fraction=0.2,
    key="dea_leiden",
    key_added="dea_leiden_filtered"
)
#To preserve the original structure of adata.uns[‘rank_genes_groups’], filtered genes are set to NaN.
#This needs to be removed because it raises an error while saving to h5ad.


# In[5]:


sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby="leiden",
    standard_scale="var",
    n_genes=5,
    key="dea_leiden_filtered",
    save="markers_raw_filtered.pdf"
)


# In[6]:


sc.get.rank_genes_groups_df(adata, group="17", key = "dea_leiden_filtered", log2fc_min=2)
##Using the filtered list will turn filtered gene names into NAs to preserve the array dimensions


# In[7]:


df = sc.get.rank_genes_groups_df(adata, group=None, key = "dea_leiden_filtered", log2fc_min=2)
df.to_csv("dea_leiden_filtered.csv")


# In[14]:


#Markers >20% in cluster <20% outside cluster >2LogFchange
cl_annotation = adata.obs["leiden"]
cl_annotation["0"] = "00_Hmcn2_Rai14"
cl_annotation["1"] = "01_Astrocyte"
cl_annotation["2"] = "02_Sim1_Chrm2"
cl_annotation["3"] = "03_Oligodendrocyte"
cl_annotation["4"] = "04_Ttc6_Lmx1a"
cl_annotation["5"] = "05_Fezf1_Lama1"
cl_annotation["6"] = "06_Trhde_Ecel1"
cl_annotation["7"] = "07_Crabp1_Tbx3"
cl_annotation["8"] = "08_Qrfpr_Alk"
cl_annotation["9"] = "09_Tanycyte"
cl_annotation["10"] = "10_Lef1_Rspo2"
cl_annotation["11"] = "11_Slc17a8_Lgr5"
cl_annotation["12"] = "12_Tafa4_Kank4"
cl_annotation["13"] = "drop_mito"
cl_annotation["14"] = "14_OPC"
cl_annotation["15"] = "15_Bcl11b_Cdh7"
cl_annotation["16"] = "16_Ghrh_Ebf1"
cl_annotation["17"] = "17_Npy_Agrp"
cl_annotation["18"] = "18_Pou1f1_Rxfp1"
cl_annotation["19"] = "19_Ntng2_Onecut1"
cl_annotation["20"] = "20_Tmem114_Npsr1"
cl_annotation["21"] = "21_Ebf3_Nxph4"
cl_annotation["22"] = "22_Dkk2_Lhx1"
cl_annotation["23"] = "23_Microglia"
cl_annotation["24"] = "24_Pth2r_Tacr1"
cl_annotation["25"] = "25_Pomc_Npy1r"##POMC is not a real neuronal marker in this combined dataset
cl_annotation["26"] = "26_Lhx1_Th"
cl_annotation["27"] = "27_Pmfbp1_Onecut1"
cl_annotation["28"] = "28_Nts"
cl_annotation["29"] = "29_Lhx9_Qrfp"
cl_annotation["30"] = "30_Bace2_Ros1"
cl_annotation["31"] = "31_Hdc_Wif1"
cl_annotation["32"] = "32_VLMC"
cl_annotation["33"] = "33_Kiss1_Tac2"
cl_annotation["34"] = "34_Fibroblast"
cl_annotation["35"] = "35_Dnah11_Pth2r"
cl_annotation["36"] = "drop_immune"
cl_annotation["37"] = "drop_histone"
cl_annotation["38"] = "38_Sim1_Pla2r1"
cl_annotation["39"] = "39_Lhx9_Zic1"
cl_annotation["40"] = "40_Corticotropic_cell"
cl_annotation["41"] = "41_Endothelial_cell"
cl_annotation["42"] = "42_Gnrhr_Fshb"
cl_annotation["43"] = "43_Gpr50"
cl_annotation["44"] = "44_VLMC"
cl_annotation["45"] = "45_OPC"
cl_annotation["46"] = "drop_unknown"
cl_annotation["47"] = "47_Ependymal_cell"
cl_annotation["48"] = "drop_unknown"
cl_annotation["49"] = "49_Otx1_Igf1"
cl_annotation["50"] = "50_Schwann_cell"
cl_annotation["51"] = "51_Oligodendrocyte"
cl_annotation["52"] = "52_Endothelial_cell"
cl_annotation["53"] = "53_Astrocyte"
adata.obs["manual_celltype_annotation"] = adata.obs.leiden.map(cl_annotation)


# In[15]:


sc.pl.embedding(adata,basis="tsne",color="manual_celltype_annotation", legend_loc="on data", legend_fontsize="xx-small", save="_manual_celltypes.pdf")


# In[10]:


cl_annotation = adata.obs["leiden"]
cl_annotation["0"] = "VMH"
cl_annotation["1"] = "Glia"
cl_annotation["2"] = "LHA"
cl_annotation["3"] = "Glia"
cl_annotation["4"] = "PrM"
cl_annotation["5"] = "VMH"
cl_annotation["6"] = "HYP"
cl_annotation["7"] = "ARC"
cl_annotation["8"] = "VMH"
cl_annotation["9"] = "Glia"
cl_annotation["10"] = "PoH"
cl_annotation["11"] = "PrM"
cl_annotation["12"] = "DMH"
cl_annotation["13"] = "drop"
cl_annotation["14"] = "Glia"
cl_annotation["15"] = "HYP"
cl_annotation["16"] = "ARC"
cl_annotation["17"] = "ARC"
cl_annotation["18"] = "PIT"
cl_annotation["19"] = "DMH"
cl_annotation["20"] = "PrM"
cl_annotation["21"] = "MaM"
cl_annotation["22"] = "MaM"
cl_annotation["23"] = "Glia"
cl_annotation["24"] = "ARC"
cl_annotation["25"] = "ARC"
cl_annotation["26"] = "SCH"
cl_annotation["27"] = "LHA"
cl_annotation["28"] = "PoH"
cl_annotation["29"] = "LHA"
cl_annotation["30"] = "HYP"
cl_annotation["31"] = "LHA"
cl_annotation["32"] = "Glia"
cl_annotation["33"] = "ARC"
cl_annotation["34"] = "Glia"
cl_annotation["35"] = "HYP"
cl_annotation["36"] = "drop"
cl_annotation["37"] = "drop"
cl_annotation["38"] = "PVH"
cl_annotation["39"] = "THA"
cl_annotation["40"] = "Glia"
cl_annotation["41"] = "Glia"
cl_annotation["42"] = "PIT"
cl_annotation["43"] = "DMH"
cl_annotation["44"] = "Glia"
cl_annotation["45"] = "Glia"
cl_annotation["46"] = "drop"
cl_annotation["47"] = "Glia"
cl_annotation["48"] = "drop"
cl_annotation["49"] = "DMH"
cl_annotation["50"] = "Glia"
cl_annotation["51"] = "Glia"
cl_annotation["52"] = "Glia"
cl_annotation["53"] = "Glia"
adata.obs["predicted_brain_region"] = adata.obs.leiden.map(cl_annotation)


# In[11]:


sc.pl.embedding(adata,basis="tsne",color="predicted_brain_region", legend_fontsize="xx-small", save="_predicted_region.pdf")


# In[16]:


cl_annotation = adata.obs["leiden"]
cl_annotation["0"] = "Neuron"
cl_annotation["1"] = "Glia"
cl_annotation["2"] = "Neuron"
cl_annotation["3"] = "Glia"
cl_annotation["4"] = "Neuron"
cl_annotation["5"] = "Neuron"
cl_annotation["6"] = "Neuron"
cl_annotation["7"] = "Neuron"
cl_annotation["8"] = "Neuron"
cl_annotation["9"] = "Glia"
cl_annotation["10"] = "Neuron"
cl_annotation["11"] = "Neuron"
cl_annotation["12"] = "Neuron"
cl_annotation["13"] = "drop"
cl_annotation["14"] = "Glia"
cl_annotation["15"] = "Neuron"
cl_annotation["16"] = "Neuron"
cl_annotation["17"] = "Neuron"
cl_annotation["18"] = "Glia"
cl_annotation["19"] = "Neuron"
cl_annotation["20"] = "Neuron"
cl_annotation["21"] = "Neuron"
cl_annotation["22"] = "Neuron"
cl_annotation["23"] = "Glia"
cl_annotation["24"] = "Neuron"
cl_annotation["25"] = "Neuron"
cl_annotation["26"] = "Neuron"
cl_annotation["27"] = "Neuron"
cl_annotation["28"] = "Neuron"
cl_annotation["29"] = "Neuron"
cl_annotation["30"] = "Neuron"
cl_annotation["31"] = "Neuron"
cl_annotation["32"] = "Glia"
cl_annotation["33"] = "Neuron"
cl_annotation["34"] = "Glia"
cl_annotation["35"] = "Neuron"
cl_annotation["36"] = "drop"
cl_annotation["37"] = "drop"
cl_annotation["38"] = "Neuron"
cl_annotation["39"] = "Neuron"
cl_annotation["40"] = "Glia"
cl_annotation["41"] = "Glia"
cl_annotation["42"] = "Glia"
cl_annotation["43"] = "Neuron"
cl_annotation["44"] = "Glia"
cl_annotation["45"] = "Glia"
cl_annotation["46"] = "drop"
cl_annotation["47"] = "Glia"
cl_annotation["48"] = "drop"
cl_annotation["49"] = "Neuron"
cl_annotation["50"] = "Glia"
cl_annotation["51"] = "Glia"
cl_annotation["52"] = "Glia"
cl_annotation["53"] = "Glia"
adata.obs["predicted_cell_type"] = adata.obs.leiden.map(cl_annotation)


# In[17]:


sc.pl.embedding(adata,basis="tsne",color="predicted_cell_type", legend_fontsize="xx-small", save="_predicted_celltype.pdf")


# In[14]:


del adata.uns['dea_leiden_filtered']
adata.write("raw_filtered_bipca_annotated.h5ad")
adata


# In[24]:


layers_to_remove = ["Y_biwhite", "Z_biwhite", "analytic_pearson_residuals", "counts", "log1p_norm"]
for layer in layers_to_remove:
    if layer in adata.layers:
        del adata.layers[layer]


# In[25]:


print(adata.X[:5, :5].toarray())


# In[26]:


neurons = adata[adata.obs["predicted_cell_type"]=="Neuron"]
sc.pl.embedding(neurons,basis="tsne",color="manual_celltype_annotation", legend_loc="on data", legend_fontsize="xx-small")


# In[27]:


glia = adata[adata.obs["predicted_cell_type"]=="Glia"]
sc.pl.embedding(glia,basis="tsne",color="manual_celltype_annotation", legend_loc="on data", legend_fontsize="xx-small")


# In[28]:


adata


# In[29]:


neurons.write("neurons_bipca_subset.h5ad")
glia.write("glia_bipca_subset.h5ad")
del neurons
del glia


# In[18]:


sc.tl.dendrogram(adata, groupby="leiden", use_rep = "pcs", n_pcs=adata.uns["bipca"]['rank'])
sc.pl.dendrogram(adata, "manual_celltype_annotation")


# In[30]:


sc.pl.correlation_matrix(adata, "manual_celltype_annotation", save="_manual_celltype.pdf")


# In[ ]:


###Now adding in the more detailed clustering for glia & neuron
adata = sc.read("raw_filtered_bipca_annotated.h5ad", cache=True)
adata


# In[3]:


adata.X = adata.layers["counts"]
print(adata.X[:5, :5].toarray())
adata


# In[4]:


layers_to_remove = ["analytic_pearson_residuals", "log1p_norm", "Y_biwhite", "Z_biwhite"]
for layer in layers_to_remove:
    if layer in adata.layers:
        del adata.layers[layer]


# In[5]:


# Pulling the detailed annotations into adata.obs
neuron_df = pd.read_csv("neurons_bipca_annotated.csv", index_col=0)
neuron_df = neuron_df[["neuronal_celltype"]]
glia_df = pd.read_csv("glia_bipca_annotated.csv", index_col=0)
glia_df = glia_df[["glia_celltype"]]
metadata_df = neuron_df.join(glia_df, how="outer")
metadata_df['detailed_celltype_annotation'] = metadata_df['neuronal_celltype'].combine_first(metadata_df['glia_celltype'])
adata.obs = adata.obs.join(metadata_df, how='left')
adata.obs


# In[6]:


sc.pl.embedding(adata,basis="tsne",color="manual_celltype_annotation", legend_loc="on data", legend_fontsize="xx-small")


# In[7]:


# Removing cells that were omitted in the annotation process
adata = adata[~adata.obs['detailed_celltype_annotation'].isna(), :]


# In[15]:


sc.tl.dendrogram(adata, groupby="leiden", use_rep = "pcs", n_pcs=adata.uns["bipca"]['rank'])
adata


# In[21]:


sc.pl.dendrogram(adata, "detailed_celltype_annotation")


# In[20]:


sc.pl.correlation_matrix(adata, "detailed_celltype_annotation", save="_detailed_celltype.pdf")


# In[13]:


sc.pl.embedding(adata,basis="tsne",color="manual_celltype_annotation", legend_loc="on data", legend_fontsize="xx-small", save="_manual_celltype_annotation.pdf")


# In[14]:


sc.pl.embedding(adata,basis="tsne",color="detailed_celltype_annotation", legend_loc="on data", legend_fontsize="xx-small", save="_detailed_celltype_annotation.pdf")


# In[3]:


#adata.write("filtered_bipca_full_annotation.h5ad")
adata = sc.read("filtered_bipca_full_annotation.h5ad", cache=True)


# In[23]:


metadata = adata.obs
metadata.to_csv("filtered_bipca_annotated.csv")


# In[4]:


pcs = adata.obsm["pcs"]
df = pd.DataFrame(pcs, index=adata.obs_names)
df.to_csv("filtered_bipca_annotated_PCS.csv")


# In[5]:


tsne = adata.obsm["tsne"]
df = pd.DataFrame(tsne, index=adata.obs_names)
df.to_csv("filtered_bipca_annotated_TSNE.csv")


# In[4]:


adata


# In[11]:


counts = pd.DataFrame(adata.layers["scran_normalization"].toarray(),
                  index=adata.obs_names,
                  columns=adata.var_names)
print(counts)


# In[7]:


meta1 = adata.obs.manual_celltype_annotation.to_frame()
meta2 = adata.obs.detailed_celltype_annotation.to_frame()


# In[12]:


counts.to_pickle("7plex_counts.pkl")


# In[10]:


meta1.to_pickle("7plex_meta1.pkl")
meta2.to_pickle("7plex_meta2.pkl")


# In[ ]:


cellex_df = pd.read_csv('7plex_original_seurat_clusters.csv')

