#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os as os
import sys
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
import scdrs
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import seaborn as sns
import statsmodels.stats.multitest as smm
import statsmodels.api as sm
import statsmodels.formula.api as smf
import adjustText
from scipy.stats import pearsonr, spearmanr


# In[2]:


celltypes = ['00_Hmcn2_Gpr149', '02_Fezf1_Lama1', '03_Rai14_Acvr1c',
       '04_Crnde_Lmx1a', '05_Qrfpr_Pdzrn3', '06_Tbx3_Npy1r', '07_Nfib_Six3',
       '08_Lef1_Rspo2', '09_Nfia_Ccbe1', '10_Npy_Agrp', '11_Tafa4_Kank4',
       '12_Ntng2_Onecut1', '13_Pth2r_Otp', '14_Lhx1_Th', '15_Onecut1_Slc32a1',
       '16_Gal_Sox6', '17_Dkk2_Ebf3', '18_Ust_Bcl11b', '19_Bace2_Tspan18',
       '20_Tmem114_Npsr1', '21_Bnc2_Glp1r', '22_Crabp1_Vgll3', '23_Ust_Nkd2',
       '24_Nts_Ecel1', '25_Nxph4_Ebf3', '26_Dach2_Npsr1', '27_Lmx1a_Ttc6',
       '28_Kiss1_Slc38a3', '29_Dnah11_Pth2r', '30_Trh_Sim1', '31_Lhx9_Npvf',
       '32_Gucy2c_Ttc6', '33_Chat_Ghrh', '34_Pla2r1_Sim1', '36_Ebf2_Crhbp',
       '37_Qrfp_Qrfprl', '38_Hdc_Prph', '39_Nfib_Nfia', '40_Gpr50_Tac2',
       '41_Otx1_Zic1', '42_Foxb1_Myo5b', 'Pituitary_Mtnr1a', 'Pituitary_Pax7',
       'Pituitary_Rxfp1']
# ── Load and filter pattern for all downstream plots ─────────────────────────
df = pd.read_csv("hypomap_cellect_unified.csv")

# Key filters
NEURONAL_C4     = (df['level']=='C4') & (df['is_neuronal'])
NEURONAL_C3     = (df['level']=='C3') & (df['is_neuronal'])
YENGO           = df['gwas'] == 'BMI_Yengo2018'
PAPER           = df['method'] == 'paper_MAGMA'
NEW_MAGMA       = df['method'] == 'new_MAGMA'
NEW_LDSC        = df['method'] == 'new_LDSC'

# ── Master cluster definitions — use these in all downstream plots ────────────

# C4
GROUP_A_C4       = ['C4-226','C4-231','C4-247','C4-248','C4-229',
                    'C4-223','C4-224','C4-252','C4-228','C4-230']
GROUP_B_C4       = ['C4-33','C4-407','C4-409','C4-411']
# C4-7 excluded entirely
LMX_C4_ALL       = GROUP_A_C4 + GROUP_B_C4

# C3
GROUP_A_C3       = ['C3-92','C3-93','C3-84','C3-94','C3-82','C3-85','C3-83']
CRNDE_C3         = ['C3-95']
GABA_NEIGHBOR_C3 = ['C3-28']
GROUP_B_C3       = ['C3-29','C3-137','C3-88','C3-139','C3-26','C3-115','C3-24']
LMX_C3_ALL       = GROUP_A_C3 + CRNDE_C3 + GABA_NEIGHBOR_C3 + GROUP_B_C3

# Highlights (for scatter emphasis)
HL_C4 = {'C4-229','C4-248'}   # Tier 1 Bonf-sig + canonical LMX1A+
HL_C3 = {'C3-84','C3-93'}     # strongest BMI + highest LMX1A

# Colour scheme
C3_GROUP_COLORS = {
    'Group A':        '#8B0000',   # genuine LMX1A+
    'CRNDE':          '#e67e22',   # cross-species marker
    'GABA_neighbor':  '#8e44ad',   # specific but LMX1A-negative
    'Group B':        '#95a5a6',   # contaminant
}
C4_GROUP_COLORS = {
    'Group A':  '#8B0000',
    'Group B':  '#95a5a6',
}


# In[ ]:


adata = ad.read_h5ad("human_HYPOMAP_snRNASeq.h5ad")
adata


# In[11]:


sc.pl.embedding(
    adata,
    basis="X_umap",
    color="celltype_annotation")


# In[5]:


sc.pl.embedding(
    adata,
    basis="X_umap",
    color="region")


# In[13]:


marker_genes = ['LMX1A', 'LMX1B', 'IRX3', 'IRX5']
sc.pl.embedding(
    adata,
    basis="X_umap",
    color=marker_genes,
    layer="norm", 
    color_map="YlOrBr", 
    save="_hypomap_lmx_irx.pdf", 
    vmax=1.5)


# In[ ]:


adata.layers["counts"] = adata.raw.X
adata.layers["norm"] = adata.X
adata


# In[ ]:


adata.X = adata.layers["counts"]
print(adata.X[:5, :5].toarray())


# In[9]:


adata.raw = None
adata.layers = None
adata


# In[ ]:


adata.write("hypomap_slim.h5ad")


# In[3]:


adata = sc.read("hypomap_slim.h5ad")
adata


# In[19]:


print(adata.X[:5, :5].toarray())
adata.layers["counts"] = adata.X
adata


# In[20]:


# Row sums of X should match nCount_RNA if X is truly raw counts
row_sums = np.array(adata.X.sum(axis=1)).flatten()
correlation = np.corrcoef(row_sums, adata.obs["nCount_RNA"])[0, 1]
print(f"Correlation between X row sums and nCount_RNA: {correlation:.6f}")
# Expect ~1.0 if X == raw counts


# In[21]:


# The target_sum is recoverable if you check median library size
# Standard pipelines use 1e4; verify by checking what makes sense
print(f"Median nCount_RNA: {adata.obs['nCount_RNA'].median():.0f}")
# After normalize_total(target_sum=1e4), each cell sums to 10,000 before log1p
target_sum = 1e4  # adjust if your pipeline used a different value (e.g. 1e6 for CPM)


# In[22]:


import scipy.sparse as sp

size_factors = adata.obs["nCount_RNA"].values

if sp.issparse(adata.X):
    normed = adata.X.multiply(target_sum / size_factors[:, None])
    adata.X = sp.csr_matrix(normed)   # ← convert COO → CSR
else:
    adata.X = adata.X / size_factors[:, None] * target_sum

sc.pp.log1p(adata)

adata.layers["norm"] = adata.X.copy()
print(adata)


# In[23]:


row0_norm = np.expm1(np.array(adata.layers["norm"][0, :].todense()).flatten())
print(f"norm row 0 expm1 sum (should ≈ 10000): {row0_norm.sum():.1f}")

row0_counts = np.array(adata.layers["counts"][0, :].todense()).flatten().sum()
print(f"counts row 0 sum vs nCount_RNA: {row0_counts:.0f} vs {adata.obs['nCount_RNA'].iloc[0]:.0f}")


# In[3]:


#adata.write("hypomap_slim_norm.h5ad")
adata = sc.read("hypomap_slim_norm.h5ad")
adata


# In[5]:


counts = pd.DataFrame(adata.layers["norm"].toarray(),
                  index=adata.obs_names,
                  columns=adata.var_names)
print(counts)
counts.to_parquet("hypomap_counts.parquet")


# In[24]:


meta1 = adata.obs.C1_named.to_frame()
meta2 = adata.obs.C2_named.to_frame()
meta3 = adata.obs.C3_named.to_frame()
meta4 = adata.obs.C4_named.to_frame()


# In[25]:


meta1.to_parquet("hypomap_c1.parquet")
meta2.to_parquet("hypomap_c2.parquet")
meta3.to_parquet("hypomap_c3.parquet")
meta4.to_parquet("hypomap_c4.parquet")


# In[4]:


adata_mouse = sc.read("allgenes_bipca_hypomap.h5ad")
adata_mouse


# In[11]:


adata.obs.to_csv("hypomap_metadata.csv", index=True)  


# In[ ]:


###Transferred to "7plex_Seurat_hypomap" in R Kernel for label transfer and gene orthologue matching


# In[ ]:


adata = sc.read("allgenes_bipca_downstream.h5ad")
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


# Create a frequency table of cell type annotations
celltype_counts = adata.obs["detailed_celltype_annotation"].value_counts()

# Convert to DataFrame for better formatting
df = pd.DataFrame({
    'cell_type': celltype_counts.index,
    'count': celltype_counts.values
})

# Export to CSV
df.to_csv('celltype_counts_7plex_detailed.csv', index=False)
df


# In[6]:


marker_genes = ["Lmx1a", "Lmx1b", "Irx3", "Irx5"]
sc.pl.embedding(adata,basis="tsne",color=marker_genes,layer="scran_normalization", color_map="YlOrBr", save="_lmx_irx.pdf", vmax=1.5)


# In[7]:


marker_genes = ["Calcr", "Agrp", "Pomc", "Kiss1"]
sc.pl.embedding(adata,basis="tsne",color=marker_genes,layer="scran_normalization", color_map="YlOrBr", save="_arcuate_markers.pdf", vmax=1.5)


# In[5]:


meta = pd.read_csv("transferred_human_hypolabels.csv", index_col="cell_id")
# Align the metadata with AnnData cells
meta = meta.loc[adata.obs_names.intersection(meta.index)]
meta


# In[6]:


adata.obs


# In[7]:


# Add to AnnData.obs
adata.obs = adata.obs.join(meta, how="left")


# In[8]:


ax = sc.pl.embedding(
    adata,
    basis="tsne",
    color="celltype_HM",
    legend_loc="on data",
    legend_fontsize="xx-small",
    return_fig=True)
ax=plt.gca()
texts = ax.texts
adjustText.adjust_text(texts, 
                       ax=ax,
                       arrowprops=dict(arrowstyle="-", color='black', lw=0.5))
plt.show()


# In[9]:


adata.obs["region_HM"] = adata.obs["region_HM"].replace("Thalamaus", "Thalamus")
ax = sc.pl.embedding(
    adata,
    basis="tsne",
    color="region_HM",
    legend_loc="on data",
    legend_fontsize="xx-small",
    return_fig=True)
ax=plt.gca()
texts = ax.texts
adjustText.adjust_text(texts, 
                       ax=ax,
                       arrowprops=dict(arrowstyle="-", color='black', lw=0.5))
plt.show()


# In[10]:


ax = sc.pl.embedding(
    adata,
    basis="tsne",
    color="C1_HM",
    legend_loc="on data",
    legend_fontsize="xx-small",
    return_fig=True)
ax=plt.gca()
texts = ax.texts
adjustText.adjust_text(texts, 
                       ax=ax,
                       arrowprops=dict(arrowstyle="-", color='black', lw=0.5))
plt.show()


# In[11]:


ax = sc.pl.embedding(
    adata,
    basis="tsne",
    color="C2_HM",
    legend_loc="on data",
    legend_fontsize="xx-small",
    return_fig=True)
ax=plt.gca()
texts = ax.texts
adjustText.adjust_text(texts, 
                       ax=ax,
                       arrowprops=dict(arrowstyle="-", color='black', lw=0.5))
plt.show()


# In[12]:


ax = sc.pl.embedding(
    adata,
    basis="tsne",
    color="predicted_spatial_region_filtered_norm",
    legend_fontsize="xx-small",
    return_fig=True)
ax=plt.gca()
texts = ax.texts
adjustText.adjust_text(texts, 
                       ax=ax,
                       arrowprops=dict(arrowstyle="-", color='black', lw=0.5))
plt.show()


# In[13]:


adata.obs["region_HM_filtered"] = adata.obs["region_HM"].copy()
adata.obs["region_HM_filtered"] = adata.obs["region_HM_filtered"].where(adata.obs["score_region_HM"] > 0.8, np.nan)


# In[14]:


#Compute the total number of cells in each region.
region_total = adata.obs.groupby("region_HM_filtered").size().rename("total_cells")
region_total_df = region_total.reset_index()  # convert to DataFrame
#Compute counts for each celltype-region combination.
counts = (
    adata.obs
    .groupby(["detailed_celltype_annotation", "region_HM_filtered"])
    .size()
    .rename("celltype_count")
    .reset_index()
)
#Merge the total count per region into the celltype-region counts.
counts = counts.merge(region_total_df, on="region_HM_filtered", how="left")
#Compute normalized frequency:
counts["normalized"] = counts["celltype_count"] / counts["total_cells"]
# For each cell type, select the region with the highest normalized value.
# The following finds the index of the maximum normalized value per celltype.
idx = counts.groupby("detailed_celltype_annotation")["normalized"].idxmax()
max_counts = counts.loc[idx]

# Create a mapping from each cell type to its predicted region.
predicted_map = dict(zip(max_counts["detailed_celltype_annotation"], max_counts["region_HM_filtered"]))

# Step 6. Apply this mapping to adata.obs to create the new column.
adata.obs["region_HM_filtered_norm"] = adata.obs["detailed_celltype_annotation"].map(predicted_map)


# In[15]:


sc.pl.embedding(
    adata,
    basis="tsne",
    color="region_HM_filtered_norm"
)


# In[16]:


# Select all columns containing "score" in their name
score_cols = ['score_celltype_HM',
 'score_region_HM',
 'score_C1_HM',
 'score_C2_HM',
 'score_C3_HM',
 'score_C4_HM']

# Convert to a long-form dataframe for plotting
scores_df = adata.obs[score_cols].melt(var_name="ScoreType", value_name="Score")

g = sns.displot(
    data=scores_df, 
    x="Score", 
    col="ScoreType", 
    bins=100, 
    facet_kws={"sharex": True, "sharey": False},
    color="steelblue"
)
g.set(xlim=(0, 1))
plt.suptitle("Per-column Distribution of TransferData Scores", y=1.02)
plt.show()


# In[17]:


adata.obs["celltype_HM_filtered"] = adata.obs["celltype_HM"].copy()
adata.obs["celltype_HM_filtered"] = adata.obs["celltype_HM_filtered"].where(adata.obs["score_celltype_HM"] > 0.95, np.nan)

adata.obs["C1_HM_filtered"] = adata.obs["C1_HM"].copy()
adata.obs["C1_HM_filtered"] = adata.obs["C1_HM_filtered"].where(adata.obs["score_C1_HM"] > 0.95, np.nan)

adata.obs["C2_HM_filtered"] = adata.obs["C2_HM"].copy()
adata.obs["C2_HM_filtered"] = adata.obs["C2_HM_filtered"].where(adata.obs["score_C2_HM"] > 0.95, np.nan)

adata.obs["C3_HM_filtered"] = adata.obs["C3_HM"].copy()
adata.obs["C3_HM_filtered"] = adata.obs["C3_HM_filtered"].where(adata.obs["score_C3_HM"] > 0.8, np.nan)

adata.obs["C4_HM_filtered"] = adata.obs["C4_HM"].copy()
adata.obs["C4_HM_filtered"] = adata.obs["C4_HM_filtered"].where(adata.obs["score_C4_HM"] > 0.8, np.nan)


# In[18]:


sum_col = "celltype_HM_filtered"
new_col = f"{sum_col}_norm"
region_total = adata.obs.groupby(sum_col).size().rename("total_cells")
region_total_df = region_total.reset_index() 
counts = (
    adata.obs
    .groupby(["detailed_celltype_annotation", sum_col])
    .size()
    .rename("celltype_count")
    .reset_index()
)
counts = counts.merge(region_total_df, on=sum_col, how="left")
counts["normalized"] = counts["celltype_count"] / counts["total_cells"]
idx = counts.groupby("detailed_celltype_annotation")["normalized"].idxmax()
max_counts = counts.loc[idx]
predicted_map = dict(zip(max_counts["detailed_celltype_annotation"], max_counts[sum_col]))
adata.obs[new_col] = adata.obs["detailed_celltype_annotation"].map(predicted_map)


# In[19]:


sc.pl.embedding(
    adata,
    basis="tsne",
    color=new_col
)


# In[20]:


sc.pl.embedding(
    adata,
    basis="tsne",
    color="neurotransmitter"
)


# In[21]:


sum_col = "C1_HM_filtered"
new_col = f"{sum_col}_norm"
region_total = adata.obs.groupby(sum_col).size().rename("total_cells")
region_total_df = region_total.reset_index() 
counts = (
    adata.obs
    .groupby(["detailed_celltype_annotation", sum_col])
    .size()
    .rename("celltype_count")
    .reset_index()
)
counts = counts.merge(region_total_df, on=sum_col, how="left")
counts["normalized"] = counts["celltype_count"] / counts["total_cells"]
idx = counts.groupby("detailed_celltype_annotation")["normalized"].idxmax()
max_counts = counts.loc[idx]
predicted_map = dict(zip(max_counts["detailed_celltype_annotation"], max_counts[sum_col]))
adata.obs[new_col] = adata.obs["detailed_celltype_annotation"].map(predicted_map)
sc.pl.embedding(
    adata,
    basis="tsne",
    color=new_col
)


# In[22]:


sum_col = "C2_HM_filtered"
new_col = f"{sum_col}_norm"
region_total = adata.obs.groupby(sum_col).size().rename("total_cells")
region_total_df = region_total.reset_index() 
counts = (
    adata.obs
    .groupby(["detailed_celltype_annotation", sum_col])
    .size()
    .rename("celltype_count")
    .reset_index()
)
counts = counts.merge(region_total_df, on=sum_col, how="left")
counts["normalized"] = counts["celltype_count"] / counts["total_cells"]
idx = counts.groupby("detailed_celltype_annotation")["normalized"].idxmax()
max_counts = counts.loc[idx]
predicted_map = dict(zip(max_counts["detailed_celltype_annotation"], max_counts[sum_col]))
adata.obs[new_col] = adata.obs["detailed_celltype_annotation"].map(predicted_map)
sc.pl.embedding(
    adata,
    basis="tsne",
    color=new_col
)


# In[23]:


adata.obs["C2_HM_NT"] = np.select(
    [
        adata.obs["C2_HM_filtered_norm"].str.contains("GABA-GLU", case=False, na=False),
        adata.obs["C2_HM_filtered_norm"].str.contains("GABA", case=False, na=False),
        adata.obs["C2_HM_filtered_norm"].str.contains("GLU", case=False, na=False)
    ],
    ["GABA_Glut", "GABA", "Glut"],
    default="Glia"
)

sc.pl.embedding(
    adata,
    basis="tsne",
    color="C2_HM_NT"
)


# In[24]:


sum_col = "C3_HM_filtered"
new_col = f"{sum_col}_norm"
region_total = adata.obs.groupby(sum_col).size().rename("total_cells")
region_total_df = region_total.reset_index() 
counts = (
    adata.obs
    .groupby(["detailed_celltype_annotation", sum_col])
    .size()
    .rename("celltype_count")
    .reset_index()
)
counts = counts.merge(region_total_df, on=sum_col, how="left")
counts["normalized"] = counts["celltype_count"] / counts["total_cells"]
idx = counts.groupby("detailed_celltype_annotation")["normalized"].idxmax()
max_counts = counts.loc[idx]
predicted_map = dict(zip(max_counts["detailed_celltype_annotation"], max_counts[sum_col]))
adata.obs[new_col] = adata.obs["detailed_celltype_annotation"].map(predicted_map)
sc.pl.embedding(
    adata,
    basis="tsne",
    color=new_col
)


# In[25]:


adata.obs["C3_HM_NT"] = np.select(
    [
        adata.obs["C3_HM_filtered_norm"].str.contains("GABA-GLU", case=False, na=False),
        adata.obs["C3_HM_filtered_norm"].str.contains("GABA", case=False, na=False),
        adata.obs["C3_HM_filtered_norm"].str.contains("GLU", case=False, na=False)
    ],
    ["GABA_Glut", "GABA", "Glut"],
    default="Glia"
)

sc.pl.embedding(
    adata,
    basis="tsne",
    color="C3_HM_NT"
)


# In[26]:


sum_col = "C4_HM_filtered"
new_col = f"{sum_col}_norm"
region_total = adata.obs.groupby(sum_col).size().rename("total_cells")
region_total_df = region_total.reset_index() 
counts = (
    adata.obs
    .groupby(["detailed_celltype_annotation", sum_col])
    .size()
    .rename("celltype_count")
    .reset_index()
)
counts = counts.merge(region_total_df, on=sum_col, how="left")
counts["normalized"] = counts["celltype_count"] / counts["total_cells"]
idx = counts.groupby("detailed_celltype_annotation")["normalized"].idxmax()
max_counts = counts.loc[idx]
predicted_map = dict(zip(max_counts["detailed_celltype_annotation"], max_counts[sum_col]))
adata.obs[new_col] = adata.obs["detailed_celltype_annotation"].map(predicted_map)
sc.pl.embedding(
    adata,
    basis="tsne",
    color=new_col
)


# In[40]:


adata.obs["C4_HM_NT"] = np.select(
    [
        adata.obs["C4_HM_filtered_norm"].str.contains("GABA-GLU", case=False, na=False),
        adata.obs["C4_HM_filtered_norm"].str.contains("GABA", case=False, na=False),
        adata.obs["C4_HM_filtered_norm"].str.contains("GLU", case=False, na=False)
    ],
    ["GABA_Glut", "GABA", "Glut"],
    default="Glia"
)

sc.pl.embedding(
    adata,
    basis="tsne",
    color="C4_HM_NT"
)


# In[3]:


#adata.write("allgenes_bipca_hypomap.h5ad")
adata = sc.read("allgenes_bipca_hypomap.h5ad", cache=True)
adata


# In[16]:


#scdrs compute-score     --h5ad-file /projects/perslab/people/qbx224/7plex_paper/atlas/allgenes_bipca_hypomap.h5ad     --h5ad-species mouse     --gs-file /projects/perslab/data/misc/Yengo_BMI.gs     --gs-species human     --flag-filter-data True     --flag-raw-count True     --flag-return-ctrl-raw-score False     --flag-return-ctrl-norm-score True     --out-folder /proj
#ects/perslab/people/qbx224/7plex_paper/atlas/scDRS/
#******************************************************************************
#* Single-cell disease relevance score (scDRS)
#* Version 1.0.4
#* Martin Jinye Zhang and Kangcheng Hou
#* HSPH / Broad Institute / UCLA
#* MIT License
#******************************************************************************
#Call: scdrs compute-score \
#--h5ad-file /projects/perslab/people/qbx224/7plex_paper/atlas/allgenes_bipca_hypomap.h5ad \
#--h5ad-species mmusculus \
#--cov-file None \
#--gs-file /projects/perslab/data/misc/Yengo_BMI.gs \
#--gs-species hsapiens \
#--ctrl-match-opt mean_var \
#--weight-opt vs \
#--adj-prop None \
#--flag-filter-data True \
#--flag-raw-count True \
#--n-ctrl 1000 \
#--min-genes 250 \
#--min-cells 50 \
#--flag-return-ctrl-raw-score False \
#--flag-return-ctrl-norm-score True \
#--out-folder /projects/perslab/people/qbx224/7plex_paper/atlas/scDRS/

#Loading data:
#/home/qbx224/scDRS/scdrs/util.py:90: FutureWarning: Use sc.pp.normalize_total instead
#  sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
#/home/qbx224/.conda/envs/pca_env/lib/python3.10/site-packages/scanpy/preprocessing/_simple.py:588: FutureWarning: Use sc.pp.normalize_total instead
#  normalize_per_cell(
#--h5ad-file loaded: n_cell=178431, n_gene=22176 (sys_time=350.4s)
#First 3 cells: ['AAACCCAAGCTCGTGC_1_1', 'AAACCCAAGGAAAGGT_1_1', 'AAACCCAAGTCCGCGT_1_1']
#First 5 genes: ['Xkr4', 'Gm1992', 'Gm19938', 'Gm37381', 'Rp1']
#--gs-file loaded: n_trait=1 (sys_time=354.3s)
#Print info for first 3 traits:
#First 3 elements for 'Yengo_BMI': ['Nupr1', 'Cln3', 'Tmem219'], [10.0, 10.0, 10.0]

#Preprocessing:

#Computing scDRS score:
#Computing control scores: 100%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 1000/1000 [01:12<00:00, 13.86it/s]
#Trait=Yengo_BMI, n_gene=793: 0/178431 FDR<0.1 cells, 3/178431 FDR<0.2 cells (sys_time=948.6s)


# In[ ]:


scdrs perform-downstream \
    --h5ad-file /projects/perslab/people/qbx224/7plex_paper/atlas/allgenes_bipca_hypomap.h5ad \
    --score-file /projects/perslab/people/qbx224/7plex_paper/atlas/scDRS/Yengo_BMI.full_score.gz \
    --out-folder /projects/perslab/people/qbx224/7plex_paper/atlas/scDRS/ \
    --group-analysis detailed_celltype_annotation,hash.mcl.ID \
    --corr-analysis cellect_bmi_log10,viral_zscore_animal \
    --gene-analysis\
    --flag-filter-data True\
    --flag-raw-count True


# In[18]:


drs = pd.read_csv("scDRS/Yengo_BMI.full_score.gz", sep="\t", index_col=0)
drs


# In[22]:


print("Cells in adata:", adata.n_obs)
print("Cells in drs:", drs.shape[0])
print("Overlap:", len(set(adata.obs_names) & set(drs.index)))


# In[23]:


adata.obs["Yengo_BMI"] = (
    drs["norm_score"]
    .reindex(adata.obs_names)  # reorder to match adata, NaN if missing
    .astype(float)
)


# In[25]:


sc.pl.embedding(
    adata,
    basis="tsne",                       
    color="Yengo_BMI",
    cmap="RdBu_r"
)


# In[26]:


adata.write("allgenes_bipca_hypomap.h5ad")


# In[3]:


import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import statsmodels.formula.api as smf
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, pearsonr
from statsmodels.stats.multitest import multipletests

# =========================================================
# Settings
# =========================================================
animal_col = "hash.mcl.ID"
area_col = "area"
celltype_col = "detailed_celltype_annotation"
scdrs_col = "Yengo_BMI"
viral_gene = "v147"
norm_layer = "scran_normalization"

adata = sc.read("allgenes_bipca_hypomap.h5ad", cache=True)

# If adata is not already in memory, try these files in order
adata_candidates = [
    "allgenes_bipca_hypomap.h5ad",
    "allgenes_bipca_downstream.h5ad",
]

# Optional scDRS score file if scdrs_col is not already in adata.obs
scdrs_score_candidates = [
    "Yengo_BMI.full_score.gz",
    "scDRS/Yengo_BMI.full_score.gz",
]

# Glia background labels used earlier in your notebook
glial_labels = ["Glia"]

# Pseudobulk filters
min_cells_per_pb = 10          # minimum cells in an animal x celltype pseudobulk
min_animals_per_celltype = 3   # keep cell types observed in at least this many animals
min_animals_per_area = 3       # for simple animal-level correlations

# Region order from earlier notebook
region_order = ["BST", "CEA", "LHA", "PAG", "PBN", "PVH", "PVT"]


# =========================================================
# Helpers
# =========================================================
def to_dense_1d(x):
    if sp.issparse(x):
        x = x.toarray()
    x = np.asarray(x)
    return x.ravel()

def safe_standardize(s):
    s = pd.Series(s, copy=True)
    sd = s.std(ddof=0)
    if np.isclose(sd, 0) or np.isnan(sd):
        return pd.Series(np.zeros(len(s)), index=s.index)
    return (s - s.mean()) / sd

def fit_cluster_wls(formula, data, weight_col, cluster_col):
    fit = smf.wls(formula=formula, data=data, weights=data[weight_col]).fit(
        cov_type="cluster",
        cov_kwds={"groups": data[cluster_col]}
    )
    return fit

def extract_area_slopes(fit, areas, suffix):
    rows = []
    for area in areas:
        candidates = [
            f"C({area_col})[{area}]:{suffix}",
            f"C({area_col})[T.{area}]:{suffix}",
        ]
        found = None
        for nm in candidates:
            if nm in fit.params.index:
                found = nm
                break

        if found is None:
            pattern = re.compile(rf"^C\({re.escape(area_col)}\)\[(T\.)?{re.escape(area)}\]:{re.escape(suffix)}$")
            matches = [nm for nm in fit.params.index if pattern.match(nm)]
            if len(matches) == 1:
                found = matches[0]

        if found is None:
            rows.append({
                "area": area,
                "beta": np.nan,
                "se": np.nan,
                "pvalue": np.nan,
                "ci_low": np.nan,
                "ci_high": np.nan,
            })
        else:
            ci = fit.conf_int().loc[found]
            rows.append({
                "area": area,
                "beta": float(fit.params[found]),
                "se": float(fit.bse[found]),
                "pvalue": float(fit.pvalues[found]),
                "ci_low": float(ci.iloc[0]),
                "ci_high": float(ci.iloc[1]),
            })
    return pd.DataFrame(rows)

def compute_viral_zscore_from_glia(adata, animal_col, neuro_col="neurotransmitter", glial_labels=("Glia",), gene="v147", layer="scran_normalization"):
    if gene not in adata.var_names:
        raise ValueError(f"{gene} not found in adata.var_names")

    gene_idx = adata.var_names.get_loc(gene)
    z_scores = np.full(adata.n_obs, np.nan, dtype=float)

    for animal_id in adata.obs[animal_col].astype(str).unique():
        animal_mask = (adata.obs[animal_col].astype(str) == str(animal_id)).values
        ad_animal = adata[animal_mask].copy()

        if neuro_col not in ad_animal.obs.columns:
            raise ValueError(f"{neuro_col} not found in adata.obs; cannot compute glia-based z-score.")

        glia_mask = ad_animal.obs[neuro_col].isin(glial_labels).values
        if glia_mask.sum() == 0:
            continue

        expr = to_dense_1d(ad_animal.layers[layer][:, gene_idx])
        glia_expr = expr[glia_mask]

        mean_bg = np.mean(glia_expr)
        std_bg = np.std(glia_expr)

        if np.isclose(std_bg, 0):
            continue

        z = (expr - mean_bg) / std_bg
        z_scores[animal_mask] = z

    return z_scores


# =========================================================
# 1) Load or reuse AnnData
# =========================================================
try:
    adata
    print("Using existing `adata` object already in memory.")
except NameError:
    adata = None
    for p in adata_candidates:
        if Path(p).exists():
            print(f"Loading {p}")
            adata = sc.read(p)
            break
    if adata is None:
        raise FileNotFoundError(
            "No `adata` object in memory and none of these files were found:\n"
            + "\n".join(adata_candidates)
        )

print("\nADATA")
print(adata)
print("\nobs columns (first 30):")
print(list(adata.obs.columns[:30]))


# =========================================================
# 2) Ensure required columns exist
# =========================================================
required_obs = [animal_col, area_col, celltype_col]
missing_obs = [c for c in required_obs if c not in adata.obs.columns]
if missing_obs:
    raise ValueError(f"Missing required columns in adata.obs: {missing_obs}")

if viral_gene not in adata.var_names:
    raise ValueError(f"{viral_gene} not found in adata.var_names")

if norm_layer not in adata.layers:
    raise ValueError(f"{norm_layer} not found in adata.layers")


# =========================================================
# 3) Add scDRS score if missing
# =========================================================
if scdrs_col not in adata.obs.columns:
    score_path = next((p for p in scdrs_score_candidates if Path(p).exists()), None)
    if score_path is None:
        raise FileNotFoundError(
            f"`{scdrs_col}` not present in adata.obs and no score file found in:\n"
            + "\n".join(scdrs_score_candidates)
        )

    print(f"\nLoading scDRS scores from {score_path}")
    drs = pd.read_csv(score_path, sep="\t", index_col=0)
    if "norm_score" not in drs.columns:
        raise ValueError("Expected a `norm_score` column in scDRS score file.")

    adata.obs[scdrs_col] = (
        drs["norm_score"]
        .reindex(adata.obs_names)
        .astype(float)
    )

print(f"\nscDRS column `{scdrs_col}` present. Non-null cells: {adata.obs[scdrs_col].notna().sum()}")


# =========================================================
# 4) Add v147 expression and viral z-score if needed
# =========================================================
gene_idx = adata.var_names.get_loc(viral_gene)
adata.obs["v147_expr"] = to_dense_1d(adata.layers[norm_layer][:, gene_idx])

if "viral_zscore_animal" not in adata.obs.columns:
    print("\nComputing `viral_zscore_animal` from glial background per animal...")
    adata.obs["viral_zscore_animal"] = compute_viral_zscore_from_glia(
        adata=adata,
        animal_col=animal_col,
        neuro_col="neurotransmitter",
        glial_labels=glial_labels,
        gene=viral_gene,
        layer=norm_layer,
    )

adata.obs["v147_positive_z2"] = adata.obs["viral_zscore_animal"] > 2

print("\nv147 summary:")
print(adata.obs[["v147_expr", "viral_zscore_animal", "v147_positive_z2"]].describe(include="all"))


# =========================================================
# 5) Sanity checks on animal-region design
# =========================================================
design_check = (
    adata.obs[[animal_col, area_col]]
    .drop_duplicates()
    .groupby(animal_col)[area_col]
    .nunique()
)

if (design_check > 1).any():
    bad_animals = design_check[design_check > 1]
    raise ValueError(
        "Some animals map to more than one traced region, which breaks the assumed design:\n"
        f"{bad_animals}"
    )

animal_area = (
    adata.obs[[animal_col, area_col]]
    .drop_duplicates()
    .copy()
)

print("\nAnimals per area:")
print(animal_area[area_col].value_counts().reindex(region_order).dropna())


# =========================================================
# 6) Build pseudobulk tables
#    A. Primary: animal x celltype
#    B. Sensitivity: animal only
# =========================================================
obs_df = adata.obs[[animal_col, area_col, celltype_col, scdrs_col, "v147_expr", "viral_zscore_animal", "v147_positive_z2"]].copy()

# Drop rows without scDRS
obs_df = obs_df.loc[obs_df[scdrs_col].notna()].copy()

# Primary pseudobulk: animal x region x cell type
pb = (
    obs_df
    .groupby([animal_col, area_col, celltype_col], observed=True)
    .agg(
        n_cells=(scdrs_col, "size"),
        mean_scdrs=(scdrs_col, "mean"),
        median_scdrs=(scdrs_col, "median"),
        mean_v147=("v147_expr", "mean"),
        median_v147=("v147_expr", "median"),
        mean_viral_z=("viral_zscore_animal", "mean"),
        pct_v147_pos=("v147_positive_z2", "mean"),
    )
    .reset_index()
)

# Filter low-cell pseudobulks
pb = pb.loc[pb["n_cells"] >= min_cells_per_pb].copy()

# Keep cell types seen in enough animals
ct_counts = pb.groupby(celltype_col)[animal_col].nunique()
keep_ct = ct_counts[ct_counts >= min_animals_per_celltype].index
pb = pb.loc[pb[celltype_col].isin(keep_ct)].copy()

# Sensitivity pseudobulk: animal only
pb_animal = (
    obs_df
    .groupby([animal_col, area_col], observed=True)
    .agg(
        n_cells=(scdrs_col, "size"),
        mean_scdrs=(scdrs_col, "mean"),
        median_scdrs=(scdrs_col, "median"),
        mean_v147=("v147_expr", "mean"),
        median_v147=("v147_expr", "median"),
        mean_viral_z=("viral_zscore_animal", "mean"),
        pct_v147_pos=("v147_positive_z2", "mean"),
    )
    .reset_index()
)

print("\nPSEUDOBULK TABLES")
print("animal x celltype rows:", len(pb))
print("animal-only rows:", len(pb_animal))
print("cell types kept:", pb[celltype_col].nunique())


# =========================================================
# 7) Standardize predictors WITHIN REGION
#    so each slope is interpretable within a traced region
# =========================================================
for metric in ["mean_v147", "mean_viral_z", "pct_v147_pos"]:
    pb[f"{metric}_z_area"] = (
        pb.groupby(area_col, observed=True)[metric]
        .transform(safe_standardize)
    )
    pb_animal[f"{metric}_z_area"] = (
        pb_animal.groupby(area_col, observed=True)[metric]
        .transform(safe_standardize)
    )

# Make area categorical in desired order if possible
present_areas = [a for a in region_order if a in pb[area_col].unique()]
pb[area_col] = pd.Categorical(pb[area_col], categories=present_areas, ordered=True)
pb_animal[area_col] = pd.Categorical(pb_animal[area_col], categories=present_areas, ordered=True)


# =========================================================
# 8) PRIMARY MODEL
#    mean BMI scDRS ~ region-specific v147 slope + celltype FE
#    weighted by cells per pseudobulk
#    SE clustered by animal
# =========================================================
primary_results = []

# Model A: mean_v147
formula_mean_v147 = (
    f"mean_scdrs ~ 0 + C({area_col}) + C({celltype_col}) + "
    f"0 + C({area_col}):mean_v147_z_area"
)
fit_mean_v147 = fit_cluster_wls(
    formula=formula_mean_v147,
    data=pb,
    weight_col="n_cells",
    cluster_col=animal_col,
)

res_mean_v147 = extract_area_slopes(fit_mean_v147, present_areas, "mean_v147_z_area")
res_mean_v147["predictor"] = "mean_v147"
primary_results.append(res_mean_v147)

# Model B: pct_v147_pos
formula_pct_v147 = (
    f"mean_scdrs ~ 0 + C({area_col}) + C({celltype_col}) + "
    f"0 + C({area_col}):pct_v147_pos_z_area"
)
fit_pct_v147 = fit_cluster_wls(
    formula=formula_pct_v147,
    data=pb,
    weight_col="n_cells",
    cluster_col=animal_col,
)

res_pct_v147 = extract_area_slopes(fit_pct_v147, present_areas, "pct_v147_pos_z_area")
res_pct_v147["predictor"] = "pct_v147_pos"
primary_results.append(res_pct_v147)

# Model C: mean_viral_z (sensitivity, closer to earlier notebook)
formula_mean_vz = (
    f"mean_scdrs ~ 0 + C({area_col}) + C({celltype_col}) + "
    f"0 + C({area_col}):mean_viral_z_z_area"
)
fit_mean_vz = fit_cluster_wls(
    formula=formula_mean_vz,
    data=pb,
    weight_col="n_cells",
    cluster_col=animal_col,
)

res_mean_vz = extract_area_slopes(fit_mean_vz, present_areas, "mean_viral_z_z_area")
res_mean_vz["predictor"] = "mean_viral_z"
primary_results.append(res_mean_vz)

primary_results_df = pd.concat(primary_results, ignore_index=True)

# FDR within each predictor across the 7 areas
primary_results_df["fdr_by_predictor"] = np.nan
for pred in primary_results_df["predictor"].unique():
    idx = primary_results_df["predictor"] == pred
    primary_results_df.loc[idx, "fdr_by_predictor"] = multipletests(
        primary_results_df.loc[idx, "pvalue"].fillna(1.0),
        method="fdr_bh"
    )[1]

# Global FDR across all area x predictor tests
primary_results_df["fdr_global"] = multipletests(
    primary_results_df["pvalue"].fillna(1.0),
    method="fdr_bh"
)[1]


# =========================================================
# 9) SENSITIVITY: simple animal-level correlations within each region
# =========================================================
animal_corr_rows = []

for area in present_areas:
    df_area = pb_animal.loc[pb_animal[area_col] == area].copy()
    if len(df_area) < min_animals_per_area:
        continue

    for predictor in ["mean_v147", "mean_viral_z", "pct_v147_pos"]:
        x = df_area[predictor].astype(float).values
        y = df_area["mean_scdrs"].astype(float).values

        try:
            sp = spearmanr(x, y)
            sp_r = float(sp.statistic if hasattr(sp, "statistic") else sp[0])
            sp_p = float(sp.pvalue if hasattr(sp, "pvalue") else sp[1])
        except Exception:
            sp_r, sp_p = np.nan, np.nan

        try:
            pe = pearsonr(x, y)
            pe_r = float(pe.statistic if hasattr(pe, "statistic") else pe[0])
            pe_p = float(pe.pvalue if hasattr(pe, "pvalue") else pe[1])
        except Exception:
            pe_r, pe_p = np.nan, np.nan

        animal_corr_rows.append({
            "area": area,
            "predictor": predictor,
            "n_animals": len(df_area),
            "spearman_r": sp_r,
            "spearman_p": sp_p,
            "pearson_r": pe_r,
            "pearson_p": pe_p,
        })

animal_corr_df = pd.DataFrame(animal_corr_rows)
if not animal_corr_df.empty:
    animal_corr_df["fdr_spearman"] = np.nan
    animal_corr_df["fdr_pearson"] = np.nan
    for pred in animal_corr_df["predictor"].unique():
        idx = animal_corr_df["predictor"] == pred
        animal_corr_df.loc[idx, "fdr_spearman"] = multipletests(
            animal_corr_df.loc[idx, "spearman_p"].fillna(1.0),
            method="fdr_bh"
        )[1]
        animal_corr_df.loc[idx, "fdr_pearson"] = multipletests(
            animal_corr_df.loc[idx, "pearson_p"].fillna(1.0),
            method="fdr_bh"
        )[1]


# =========================================================
# 10) Print results
# =========================================================
print("\n" + "=" * 100)
print("PRIMARY PSEUDOBULK MODEL")
print("Outcome: mean BMI scDRS per animal x celltype pseudobulk")
print("Model: mean_scdrs ~ area-specific v147 slope + celltype FE")
print("Weights: n_cells per pseudobulk")
print("SE: clustered by animal")
print("=" * 100)

for pred in ["mean_v147", "mean_viral_z", "pct_v147_pos"]:
    tmp = (
        primary_results_df[primary_results_df["predictor"] == pred]
        .sort_values("pvalue")
        .copy()
    )
    print(f"\n--- Predictor: {pred} ---")
    print(
        tmp[[
            "area", "beta", "se", "ci_low", "ci_high", "pvalue",
            "fdr_by_predictor", "fdr_global"
        ]].to_string(index=False)
    )

print("\n" + "=" * 100)
print("ANIMAL-LEVEL SENSITIVITY CORRELATIONS")
print("Outcome: mean BMI scDRS per animal")
print("=" * 100)

if animal_corr_df.empty:
    print("No animal-level correlations computed.")
else:
    for pred in ["mean_v147", "mean_viral_z", "pct_v147_pos"]:
        tmp = animal_corr_df[animal_corr_df["predictor"] == pred].sort_values("spearman_p")
        print(f"\n--- Predictor: {pred} ---")
        print(
            tmp[[
                "area", "n_animals", "spearman_r", "spearman_p", "fdr_spearman",
                "pearson_r", "pearson_p", "fdr_pearson"
            ]].to_string(index=False)
        )


# =========================================================
# 11) Save outputs
# =========================================================
pb.to_csv("pseudobulk_animal_celltype_scdrs_v147.csv", index=False)
pb_animal.to_csv("pseudobulk_animal_scdrs_v147.csv", index=False)
primary_results_df.to_csv("region_specific_scdrs_v147_pseudobulk_model_results.csv", index=False)
if not animal_corr_df.empty:
    animal_corr_df.to_csv("region_specific_scdrs_v147_animal_level_correlations.csv", index=False)

print("\nSaved:")
print(" - pseudobulk_animal_celltype_scdrs_v147.csv")
print(" - pseudobulk_animal_scdrs_v147.csv")
print(" - region_specific_scdrs_v147_pseudobulk_model_results.csv")
if not animal_corr_df.empty:
    print(" - region_specific_scdrs_v147_animal_level_correlations.csv")


# =========================================================
# 12) Plot primary results
# =========================================================
fig, axes = plt.subplots(1, 3, figsize=(16, 4), sharey=False)

for ax, pred in zip(axes, ["mean_v147", "mean_viral_z", "pct_v147_pos"]):
    tmp = primary_results_df[primary_results_df["predictor"] == pred].copy()
    tmp[area_col] = pd.Categorical(tmp["area"], categories=present_areas, ordered=True)
    tmp = tmp.sort_values("area")

    ax.errorbar(
        x=np.arange(len(tmp)),
        y=tmp["beta"],
        yerr=1.96 * tmp["se"],
        fmt="o",
        capsize=4
    )
    ax.axhline(0, linestyle="--", linewidth=1)
    ax.set_xticks(np.arange(len(tmp)))
    ax.set_xticklabels(tmp["area"], rotation=45, ha="right")
    ax.set_title(pred)
    ax.set_ylabel("Slope: mean BMI scDRS per 1 SD predictor")
    for i, row in enumerate(tmp.itertuples(index=False)):
        ax.text(i, row.beta, f"p={row.pvalue:.3g}", fontsize=8, ha="center", va="bottom")

plt.tight_layout()
plt.show()


# =========================================================
# 13) Optional: show the animals per region and top rows
# =========================================================
print("\nAnimals per region in pseudobulk:")
print(pb[[animal_col, area_col]].drop_duplicates()[area_col].value_counts().reindex(region_order).dropna())

print("\nTop pseudobulk rows:")
print(
    pb.sort_values(["mean_scdrs", "mean_v147"], ascending=[False, False]).head(20).to_string(index=False)
)


# In[5]:


# ── Analysis: label specificity + LMX1A expression overlap ───────────────────
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from scipy.sparse import issparse

MOUSE_CELLTYPES = ['04_Crnde_Lmx1a', '27_Lmx1a_Ttc6',
                   '32_Gucy2c_Ttc6', '36_Ebf2_Crhbp']
LMX_C4_IDS = ['C4-248','C4-411','C4-33','C4-252','C4-407','C4-226','C4-229',
               'C4-223','C4-247','C4-230','C4-7','C4-224','C4-409','C4-228','C4-231']
LMX_C3_IDS = ['C3-92','C3-84','C3-139','C3-29','C3-94','C3-83','C3-82','C3-26',
               'C3-137','C3-85','C3-28','C3-24','C3-93','C3-88','C3-115','C3-95']

# ── 0. Prepare obs frames ─────────────────────────────────────────────────────

# MOUSE — obs already in memory; add C4/C3 short IDs if not present
obs_mouse = adata_mouse.obs.copy()

# C3 labels
if 'C3_short' not in obs_mouse.columns:
    if 'C3_HM' in obs_mouse.columns:
        obs_mouse['C3_short'] = obs_mouse['C3_HM'].str.extract(r'(C3-\d+)')
    else:
        obs_mouse['C3_short'] = pd.NA
        print("WARNING: C3_HM not found in mouse obs")

# C4 labels — sanity check showed C4_HM is absent; try several fallbacks
if 'C4_short' not in obs_mouse.columns:
    if 'C4_HM' in obs_mouse.columns:
        obs_mouse['C4_short'] = obs_mouse['C4_HM'].str.extract(r'(C4-\d+)')
        print("C4_short extracted from C4_HM")
    elif 'celltype_HM' in obs_mouse.columns:
        test = obs_mouse['celltype_HM'].str.extract(r'(C4-\d+)')
        if test.notna().any().any():
            obs_mouse['C4_short'] = test
            print("C4_short extracted from celltype_HM")
        else:
            obs_mouse['C4_short'] = pd.NA
            print("WARNING: no C4-style labels found in mouse obs — "
                  "re-run the label transfer cell that adds C4_HM first")
    else:
        obs_mouse['C4_short'] = pd.NA
        print("WARNING: C4 labels not available — "
              "Analysis 1 / Panel A / B will be empty")

print(f"Mouse C4_short non-null: {obs_mouse['C4_short'].notna().sum()}")
print(f"Mouse C3_short non-null: {obs_mouse['C3_short'].notna().sum()}")

obs_mouse['is_lmx1a_pop'] = obs_mouse['detailed_celltype_annotation'].isin(MOUSE_CELLTYPES)

# HUMAN — obs already has C4_short and C3_short from earlier cells ✓
obs_human = adata.obs.copy()
# Verify they exist; add from C3/C4 named columns if somehow absent
for level in ['C3','C4']:
    short_col = f'{level}_short'
    if short_col not in obs_human.columns:
        named_col = f'{level}_named' if f'{level}_named' in obs_human.columns else level
        if named_col in obs_human.columns:
            obs_human[short_col] = obs_human[named_col].str.extract(
                rf'({level}-\d+)')
            print(f"Added {short_col} from {named_col}")
        else:
            obs_human[short_col] = pd.NA
            print(f"WARNING: {short_col} unavailable for human obs")

print(f"Human C4_short non-null: {obs_human['C4_short'].notna().sum()}")
print(f"Human C3_short non-null: {obs_human['C3_short'].notna().sum()}")

# ── 1. Load unified CELLECT CSV ───────────────────────────────────────────────
cellect = pd.read_csv("hypomap_cellect_unified.csv")

cellect_wide = (cellect[
    (cellect['gwas']      == 'BMI_Yengo2018') &
    (cellect['level']     == 'C4') &
    (cellect['is_neuronal'])
].pivot_table(
    index=['cluster_id','has_exome_sig','n_effector_genes'],
    columns='method',
    values=['pvalue','log10p','bonf_sig','rank']
).reset_index())

cellect_wide.columns = [
    '_'.join(c).strip('_') if c[1] else c[0]
    for c in cellect_wide.columns
]
cellect_wide = cellect_wide.rename(columns={'cluster_id': 'C4_id'})

# ── 2. ANALYSIS 1: Mouse label specificity ────────────────────────────────────
print("\n" + "=" * 70)
print("ANALYSIS 1: C4 label specificity to mouse Lmx1a populations")
print("=" * 70)

spec_rows = []
for c4id in LMX_C4_IDS:
    cells_with_label = obs_mouse[obs_mouse['C4_short'] == c4id]
    n_total  = len(cells_with_label)
    n_lmx1a  = cells_with_label['is_lmx1a_pop'].sum()
    pct      = n_lmx1a / n_total * 100 if n_total > 0 else 0.0
    spec_rows.append({
        'C4_id':     c4id,
        'n_total':   n_total,
        'n_lmx1a':   n_lmx1a,
        'n_other':   n_total - n_lmx1a,
        'pct_lmx1a': pct,
        'specific':  pct >= 50,
    })

df_spec = pd.DataFrame(spec_rows).sort_values('pct_lmx1a', ascending=False)
print(df_spec.to_string(index=False))

non_lmx = obs_mouse[
    obs_mouse['C4_short'].isin(LMX_C4_IDS) & ~obs_mouse['is_lmx1a_pop']
]
print("\nNon-Lmx1a cell types receiving Lmx1a C4 labels:")
print(non_lmx['detailed_celltype_annotation'].value_counts().head(20).to_string())

# ── 3. ANALYSIS 2: LMX1A expression in human hypomap ─────────────────────────
print("\n" + "=" * 70)
print("ANALYSIS 2: LMX1A expression across human C4 clusters")
print("=" * 70)

human_gene = None
for candidate in ['LMX1A', 'Lmx1a', 'lmx1a']:
    if candidate in adata.var_names:
        human_gene = candidate
        break
if human_gene is None:
    matches = [g for g in adata.var_names if g.upper() == 'LMX1A']
    if matches:
        human_gene = matches[0]

if human_gene:
    print(f"Using gene: {human_gene}")
    gene_idx = list(adata.var_names).index(human_gene)
    X_norm   = adata.layers['norm']          # confirmed layer name ✓
    lmx1a_expr = (np.asarray(X_norm[:, gene_idx].todense()).ravel()
                  if issparse(X_norm)
                  else np.asarray(X_norm[:, gene_idx]).ravel())
    obs_human['LMX1A_expr'] = lmx1a_expr

    # Per-C4-cluster expression stats
    c4_lmx1a = (obs_human.groupby('C4_short')['LMX1A_expr']
                .agg(mean_expr='mean',
                     pct_expr=lambda x: (x > 0).mean() * 100,
                     n_cells='count')
                .reset_index()
                .sort_values('mean_expr', ascending=False)
                .reset_index(drop=True))
    c4_lmx1a['rank_pos'] = range(len(c4_lmx1a))   # for Panel C

    #expr_thresh = c4_lmx1a['mean_expr'].quantile(0.90)
    expr_thresh = 0.5
    pct_thresh  = 10.0
    c4_lmx1a['high_LMX1A'] = (
        (c4_lmx1a['mean_expr'] > expr_thresh) &
        (c4_lmx1a['pct_expr']  > pct_thresh)
    )

    print(f"LMX1A expression threshold (90th pct): {expr_thresh:.4f}")
    print(f"C4 clusters with high LMX1A: {c4_lmx1a['high_LMX1A'].sum()}")

    c4_lmx1a_focus = c4_lmx1a[c4_lmx1a['C4_short'].isin(LMX_C4_IDS)]
    print("\nLMX1A in Lmx1a C4 clusters (human):")
    print(c4_lmx1a_focus.sort_values('mean_expr', ascending=False).to_string(index=False))
else:
    print("LMX1A not found in human adata.var_names")
    obs_human['LMX1A_expr'] = 0.0
    c4_lmx1a = pd.DataFrame(columns=['C4_short','mean_expr','pct_expr',
                                      'n_cells','rank_pos','high_LMX1A'])
    expr_thresh = 0.0

# ── 4. ANALYSIS 3: Combined overlap table ────────────────────────────────────
print("\n" + "=" * 70)
print("ANALYSIS 3: Label transfer ∩ LMX1A expression ∩ BMI enrichment")
print("=" * 70)

overlap_rows = []
for c4id in LMX_C4_IDS:
    human_cells = obs_human[obs_human['C4_short'] == c4id]
    if len(human_cells) == 0:
        overlap_rows.append({'C4_id': c4id, 'n_human': 0,
                              'mean_expr': np.nan, 'pct_expr': np.nan})
        continue
    expr = human_cells['LMX1A_expr']
    overlap_rows.append({
        'C4_id':     c4id,
        'n_human':   len(human_cells),
        'mean_expr': expr.mean(),
        'pct_expr':  (expr > 0).mean() * 100,
    })

df_overlap   = pd.DataFrame(overlap_rows)
df_combined  = (df_spec
                .merge(df_overlap[['C4_id','n_human','mean_expr','pct_expr']],
                       on='C4_id', how='left')
                .rename(columns={'mean_expr':'human_LMX1A_mean',
                                 'pct_expr': 'human_LMX1A_pct'}))

# Merge BMI enrichment from unified CSV
df_combined_full = df_combined.merge(cellect_wide, on='C4_id', how='left')

# ── 5. Figures ────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(16, 13))
COL_HL  = '#8B0000'
COL_LMX = '#4878CF'

# ── Panel A: Mouse label specificity bar chart ────────────────────────────────
ax = axes[0, 0]
df_sp = df_spec.sort_values('pct_lmx1a', ascending=True)
colors_sp = [COL_HL if c in {'C4-248','C4-411'} else COL_LMX
             for c in df_sp['C4_id']]
ax.barh(range(len(df_sp)), df_sp['pct_lmx1a'],
        color=colors_sp, edgecolor='white', linewidth=0.4, height=0.72)
ax.axvline(50, color='black', ls='--', lw=1, label='50% threshold')
ax.axvline(50, color='grey',  ls=':',  lw=0.8)
for i, (_, r) in enumerate(df_sp.iterrows()):
    ax.text(r['pct_lmx1a'] + 0.5, i,
            f"{int(r['n_lmx1a'])}/{int(r['n_total'])}",
            va='center', fontsize=6.5, color='grey')
ax.set_yticks(range(len(df_sp)))
ax.set_yticklabels(df_sp['C4_id'], fontsize=8)
for lbl, c4 in zip(ax.get_yticklabels(), df_sp['C4_id']):
    if c4 in {'C4-248','C4-411'}:
        lbl.set_color(COL_HL); lbl.set_fontweight('bold')
ax.set_xlabel('% cells with C4 label that are Lmx1a mouse population')
ax.set_title('C4 label specificity to\nmouse Lmx1a neurons', fontweight='bold')
ax.set_xlim(0, 110); ax.legend(fontsize=8)

# ── Panel B: Non-Lmx1a contamination ─────────────────────────────────────────
ax = axes[0, 1]
top_contam = non_lmx['detailed_celltype_annotation'].value_counts().head(15)
ax.barh(range(len(top_contam)), top_contam.values,
        color='#95a5a6', edgecolor='white', linewidth=0.4, height=0.72)
ax.set_yticks(range(len(top_contam)))
ax.set_yticklabels(top_contam.index, fontsize=7.5)
ax.set_xlabel('Number of cells')
ax.set_title('Non-Lmx1a mouse cell types\nalso receiving Lmx1a C4 labels',
             fontweight='bold')
for i, (_, n) in enumerate(zip(top_contam.index, top_contam.values)):
    ax.text(n + 1, i, str(n), va='center', fontsize=7, color='grey')

# ── Panel C: LMX1A expression rank across all C4 clusters ────────────────────
ax = axes[1, 0]
if len(c4_lmx1a) > 0:
    ax.scatter(c4_lmx1a['rank_pos'], c4_lmx1a['mean_expr'],
               c='#cccccc', s=15, alpha=0.5, zorder=2,
               label=f'All C4 (n={len(c4_lmx1a)})')

    # Overlay Lmx1a clusters — all present since we plot all clusters
    lmx_in_ranked = c4_lmx1a[c4_lmx1a['C4_short'].isin(LMX_C4_IDS)]
    for _, row in lmx_in_ranked.iterrows():
        col = COL_HL if row['C4_short'] in {'C4-248','C4-411'} else '#E07B54'
        ax.scatter(row['rank_pos'], row['mean_expr'],
                   c=col, s=70, zorder=4, edgecolors='white', linewidths=0.5)
        ax.annotate(row['C4_short'],
                    xy=(row['rank_pos'], row['mean_expr']),
                    xytext=(4, 3), textcoords='offset points',
                    fontsize=6, color=col)

    ax.axhline(expr_thresh, color='black', ls='--', lw=1,
               label=f'90th pct ({expr_thresh:.3f})')

ax.set_xlabel('C4 cluster rank (by LMX1A expression, 1=highest)')
ax.set_ylabel('Mean LMX1A expression (norm layer)')
ax.set_title('LMX1A expression across human\nhypomap C4 clusters', fontweight='bold')
ax.legend(handles=[
    mpatches.Patch(fc=COL_HL,    label='C4-248 / C4-411'),
    mpatches.Patch(fc='#E07B54', label='Other Lmx1a C4s'),
    mpatches.Patch(fc='#cccccc', label='All other C4 clusters'),
], fontsize=7.5)

# ── Panel D: 2D — specificity vs LMX1A, coloured by BMI enrichment ───────────
ax = axes[1, 1]

# BMI colour = paper_MAGMA log10p (most informative for colour axis)
df_combined_full['log10p_pm'] = df_combined_full['log10p_paper_MAGMA'].fillna(0)

if df_combined_full['human_LMX1A_mean'].notna().any():
    sc = ax.scatter(
        df_combined_full['pct_lmx1a'],
        df_combined_full['human_LMX1A_mean'],
        c=df_combined_full['log10p_pm'],
        cmap='YlOrRd',
        vmin=0,
        vmax=df_combined_full['log10p_pm'].quantile(0.95) + 0.01,
        s=[max(40, n / 50) for n in df_combined_full['n_human'].fillna(0)],
        edgecolors=['black' if c in {'C4-248','C4-411'} else 'white'
                    for c in df_combined_full['C4_id']],
        linewidths=[1.2 if c in {'C4-248','C4-411'} else 0.4
                    for c in df_combined_full['C4_id']],
        zorder=3, alpha=0.9
    )
    plt.colorbar(sc, ax=ax, label='−log₁₀(p)  paper MAGMA (Yengo2018)')

    ax.axvline(50,          color='black', ls='--', lw=1, alpha=0.5)
    ax.axhline(expr_thresh, color='grey',  ls='--', lw=1, alpha=0.5)

    for _, r in df_combined_full.iterrows():
        if pd.notna(r['human_LMX1A_mean']):
            ax.annotate(r['C4_id'],
                        xy=(r['pct_lmx1a'], r['human_LMX1A_mean']),
                        xytext=(4, 2), textcoords='offset points',
                        fontsize=6.5,
                        color=COL_HL if r['C4_id'] in {'C4-248','C4-411'}
                              else '#333333')

    # Quadrant annotations
    xlim = ax.get_xlim(); ylim = ax.get_ylim()
    yspan = ylim[1] - ylim[0]
    ax.text(82,  ylim[1] - yspan*0.05, 'Specific + LMX1A high',
            fontsize=7.5, color='#2ecc71', fontweight='bold', va='top')
    ax.text(2,   ylim[1] - yspan*0.05, 'Non-specific + LMX1A high',
            fontsize=7.5, color='#e67e22', fontweight='bold', va='top')
    ax.text(82,  ylim[0] + yspan*0.02, 'Specific + LMX1A low',
            fontsize=7.5, color='#3498db', fontweight='bold')

ax.set_xlabel('% cells with label from Lmx1a mouse population', fontsize=9)
ax.set_ylabel('Mean LMX1A expression in human cluster', fontsize=9)
ax.set_title('Label specificity vs human LMX1A expression\n'
             '(colour = BMI enrichment, bubble = n human cells)',
             fontweight='bold')

plt.suptitle('Mouse Lmx1a label specificity × human LMX1A expression × BMI enrichment',
             fontweight='bold', fontsize=13, y=1.01)
plt.tight_layout()
plt.savefig('LMX_label_specificity_human_expression.pdf', bbox_inches='tight')
plt.savefig('LMX_label_specificity_human_expression.png', bbox_inches='tight')
plt.show()

# ── 6. Human UMAP: LMX1A expression + cluster highlights ─────────────────────
if 'X_umap' in adata.obsm and human_gene:
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    umap = adata.obsm['X_umap']

    ax = axes[0]
    sc_expr = ax.scatter(umap[:, 0], umap[:, 1],
                         c=obs_human['LMX1A_expr'],
                         cmap='YlOrRd', s=0.3, alpha=0.6,
                         vmin=0,
                         vmax=np.percentile(obs_human['LMX1A_expr'], 99))
    plt.colorbar(sc_expr, ax=ax, label='LMX1A expression (norm)')
    ax.set_title('Human HYPOMAP UMAP\nLMX1A expression', fontweight='bold')
    ax.axis('off')

    ax = axes[1]
    ax.scatter(umap[:, 0], umap[:, 1], c='#eeeeee', s=0.2, alpha=0.3)
    for c4id, col in [('C4-248', COL_HL), ('C4-411', '#C44E52')]:
        mask = obs_human['C4_short'] == c4id
        ax.scatter(umap[mask, 0], umap[mask, 1],
                   c=col, s=2, alpha=0.8, label=c4id)
    ax.set_title('Human HYPOMAP UMAP\nC4-248/411',
                 fontweight='bold')
    ax.legend(markerscale=5, fontsize=8)
    ax.axis('off')

    plt.tight_layout()
    plt.savefig('LMX_human_umap_expression.pdf', bbox_inches='tight')
    plt.savefig('LMX_human_umap_expression.png', bbox_inches='tight')
    plt.show()

# ── 7. Combined summary table ─────────────────────────────────────────────────
print("\n" + "=" * 95)
print("COMBINED SUMMARY: mouse label specificity + human LMX1A + CELLECT BMI")
print("=" * 95)
print(f"{'C4':^10} {'%Lmx1a':>8} {'Spec':>5} {'LMX1A_mean':>11} "
      f"{'p_paperM':>11} {'p_newM':>9} {'p_LDSC':>8} "
      f"{'ExomeSig':>9} {'nEff':>5}")
print("-" * 80)
for _, r in df_combined_full.sort_values('pct_lmx1a', ascending=False).iterrows():
    spec = '✓' if r['specific'] else ' '
    hmean = f"{r['human_LMX1A_mean']:.4f}" if pd.notna(r.get('human_LMX1A_mean')) else 'N/A'
    pp  = f"{r.get('pvalue_paper_MAGMA',np.nan):.3e}" if pd.notna(r.get('pvalue_paper_MAGMA')) else 'N/A'
    pnm = f"{r.get('pvalue_new_MAGMA',  np.nan):.3e}" if pd.notna(r.get('pvalue_new_MAGMA'))   else 'N/A'
    pl  = f"{r.get('pvalue_new_LDSC',   np.nan):.3e}" if pd.notna(r.get('pvalue_new_LDSC'))    else 'N/A'
    exm = '✓' if r.get('has_exome_sig', False) else ' '
    neff = int(r.get('n_effector_genes', 0))
    print(f"{r['C4_id']:^10} {r['pct_lmx1a']:>7.1f}% {spec:>5} "
          f"{hmean:>11} {pp:>11} {pnm:>9} {pl:>8} {exm:>9} {neff:>5}")


# In[7]:


import matplotlib.patheffects as mpe
from adjustText import adjust_text

fig, ax = plt.subplots(figsize=(9, 7))

# ── Colour by biological group (no BMI) ──────────────────────────────────────
def point_color(c4_id):
    if c4_id in set(GROUP_A_C4):   return C4_GROUP_COLORS['Group A']   # '#8B0000'
    if c4_id in set(GROUP_B_C4):   return C4_GROUP_COLORS['Group B']   # '#95a5a6'
    return '#cccccc'

colors = [point_color(c) for c in df_combined_full['C4_id']]

if df_combined_full['human_LMX1A_mean'].notna().any():
    ax.scatter(
        df_combined_full['pct_lmx1a'],
        df_combined_full['human_LMX1A_mean'],
        c=colors,
        s=80,
        edgecolors='white',
        linewidths=0.4,
        zorder=3, alpha=0.9
    )

    ax.axvline(50,          color='black', ls='--', lw=1, alpha=0.5)
    ax.axhline(expr_thresh, color='grey',  ls='--', lw=1, alpha=0.5)

    # ── Labels ────────────────────────────────────────────────────────────────
    texts = []
    for _, r in df_combined_full.iterrows():
        if pd.notna(r['human_LMX1A_mean']):
            t = ax.text(
                r['pct_lmx1a'], r['human_LMX1A_mean'],
                r['C4_id'],
                fontsize=6.5,
                color=(C4_GROUP_COLORS['Group A']
                       if r['C4_id'] in set(GROUP_A_C4) else '#555555'),
                zorder=5
            )
            texts.append(t)
    adjust_text(
        texts, ax=ax,
        expand=(1.15, 1.3),
        arrowprops=dict(arrowstyle='-', color='#aaaaaa', lw=0.5),
    )

    # ── Quadrant labels ───────────────────────────────────────────────────────
    ax.autoscale(tight=False)
    xlim  = ax.get_xlim();  ylim  = ax.get_ylim()
    xspan = xlim[1]-xlim[0]; yspan = ylim[1]-ylim[0]
    pad_x = xspan*0.02;      pad_y = yspan*0.03

    for qx, qy, label, ha, col in [
        (xlim[1]-pad_x, ylim[1]-0.7, 'Specific + LMX1A high',     'right', '#000000'),
        (xlim[0]+pad_x, ylim[1]-0.7, 'Non-specific + LMX1A high', 'left',  '#000000'),
        (xlim[1]-pad_x, ylim[0]+pad_y,'Specific + LMX1A low',     'right', '#000000'),
    ]:
        ax.text(qx, qy, label, fontsize=7.5, color=col, fontweight='bold',
                ha=ha,
                va='top' if qy > (ylim[0]+yspan/2) else 'bottom',
                path_effects=[mpe.withStroke(linewidth=2, foreground='white')],
                clip_on=True)

    # ── Legend ────────────────────────────────────────────────────────────────
    ax.legend(handles=[
        mpatches.Patch(fc=C4_GROUP_COLORS['Group A'],
                       label=f'Group A — genuine LMX1A+  (n={len(GROUP_A_C4)})'),
        mpatches.Patch(fc=C4_GROUP_COLORS['Group B'],
                       label=f'Group B — contaminant labels (n={len(GROUP_B_C4)})'),
    ], fontsize=8, loc='upper left', framealpha=0.9)

ax.set_xlabel('% cells with C4 label from mouse Lmx1a population', fontsize=9)
ax.set_ylabel('Mean normalized LMX1A expression in human cluster', fontsize=9)
ax.set_title('C4 label specificity vs human LMX1A expression\n'
             'Group A (genuine LMX1A+) vs Group B (contaminant)',
             fontweight='bold')
plt.tight_layout()
plt.savefig('LMX_label_specificity_C4.pdf', bbox_inches='tight')
plt.savefig('LMX_label_specificity_C4.png', bbox_inches='tight')
plt.show()


# In[28]:


# ── Sankey + LMX1A bar — C4 label assignment from mouse Lmx1a neurons ─────────
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import PathPatch
from matplotlib.path import Path

# ── 1. Compute cross-tabulation from mouse obs ────────────────────────────────
MOUSE_CELLTYPES = ['04_Crnde_Lmx1a','27_Lmx1a_Ttc6',
                   '32_Gucy2c_Ttc6','36_Ebf2_Crhbp']
MOUSE_COLORS    = {
    '04_Crnde_Lmx1a': '#008941',
    '27_Lmx1a_Ttc6':  '#ffaa92',
    '32_Gucy2c_Ttc6': '#000035',
    '36_Ebf2_Crhbp':  '#0aa6d8',
}

# Extract C4_short from C4_HM if not already present
if 'C4_short' not in adata_mouse.obs.columns:
    adata_mouse.obs['C4_short'] = (adata_mouse.obs['C4_HM']
                                   .str.extract(r'(C4-\d+)'))
    print(f"C4_short extracted: {adata_mouse.obs['C4_short'].notna().sum()} non-null")

# Cross-tab: mouse celltype rows × C4 cluster columns
lmx_obs = adata_mouse.obs[
    adata_mouse.obs['detailed_celltype_annotation'].isin(MOUSE_CELLTYPES)
].copy()
lmx_obs['C4_short'] = lmx_obs['C4_short'].astype(str)

# Only clusters present in our master sets
ALL_C4 = LMX_C4_ALL   # GROUP_A_C4 + GROUP_B_C4

xtab = (lmx_obs[lmx_obs['C4_short'].isin(ALL_C4)]
        .groupby(['detailed_celltype_annotation','C4_short'])
        .size()
        .unstack(fill_value=0))

# Ensure all C4 clusters are present
for c in ALL_C4:
    if c not in xtab.columns:
        xtab[c] = 0

# Order C4 clusters: Group A descending total, then Group B descending
c4_totals  = xtab.sum(axis=0)
c4_order   = (sorted(GROUP_A_C4, key=lambda c: -c4_totals.get(c,0)) +
               sorted(GROUP_B_C4, key=lambda c: -c4_totals.get(c,0)))
c4_order   = [c for c in c4_order if c in xtab.columns]
xtab       = xtab[c4_order]

# LMX1A expression for bar chart (from df_combined_full or recompute)
lmx1a_mean = {}
for cid in c4_order:
    mask = obs_human['C4_short'].astype(str) == cid
    if mask.sum() > 0:
        lmx1a_mean[cid] = obs_human.loc[mask, 'LMX1A_expr'].mean()
    else:
        lmx1a_mean[cid] = 0.0

# ── 2. Layout constants ───────────────────────────────────────────────────────

fig, (ax_s, ax_b) = plt.subplots(
    1, 2, figsize=(16, 8),
    gridspec_kw={'width_ratios': [3.5, 0.6], 'wspace': 0.01}  # was [2.5, 1]
)

N_CT  = len(MOUSE_CELLTYPES)     # 4 left nodes
N_C4  = len(c4_order)            # right nodes
X_L   = 0.0                      # x position of left (mouse) column
X_R   = 1.0                      # x position of right (C4) column
CURVE = 0.5                      # bezier control point offset

# ── 3. Compute node y-positions ───────────────────────────────────────────────

# Left: mouse cell types — equal height, evenly spaced
ct_totals   = xtab.sum(axis=1)   # total cells per mouse celltype
grand_total = ct_totals.sum()

# Use proportional height
ct_heights  = ct_totals / grand_total
ct_y_bottom = {}
y = 1.0
for ct in MOUSE_CELLTYPES:
    h = ct_heights.get(ct, 0)
    ct_y_bottom[ct] = y - h
    y -= h
ct_y_top = {ct: ct_y_bottom[ct] + ct_heights.get(ct, 0)
             for ct in MOUSE_CELLTYPES}

# Right: C4 clusters — proportional height + small gap between Group A and B
GAP   = 0.015
c4_heights = c4_totals[c4_order] / grand_total

y = 1.0
c4_y_bottom = {}
for i, cid in enumerate(c4_order):
    # Add gap before first Group B cluster
    if i > 0 and cid == GROUP_B_C4[0] and c4_order[i-1] in GROUP_A_C4:
        y -= GAP
    h = c4_heights[cid]
    c4_y_bottom[cid] = y - h
    y -= h
c4_y_top = {cid: c4_y_bottom[cid] + c4_heights[cid] for cid in c4_order}

# ── 4. Draw Sankey flows ──────────────────────────────────────────────────────

# Track current fill position within each node
ct_fill  = {ct: ct_y_top[ct]   for ct in MOUSE_CELLTYPES}
c4_fill  = {cid: c4_y_top[cid] for cid in c4_order}

def bezier_band(x0, y0_bot, y0_top, x1, y1_bot, y1_top, color, alpha=0.45):
    """Draw a filled bezier band (Sankey flow)."""
    cx = (x0 + x1) / 2
    verts = [
        (x0, y0_bot),
        (cx, y0_bot),
        (cx, y1_bot),
        (x1, y1_bot),
        (x1, y1_top),
        (cx, y1_top),
        (cx, y0_top),
        (x0, y0_top),
        (x0, y0_bot),
    ]
    codes = [Path.MOVETO,
             Path.CURVE4, Path.CURVE4, Path.CURVE4,
             Path.LINETO,
             Path.CURVE4, Path.CURVE4, Path.CURVE4,
             Path.CLOSEPOLY]
    # Bezier needs 3 control points; approximate with 4-point (CURVE4)
    # Rewrite with proper control points
    verts = [
        (x0, y0_bot),                          # start bottom-left
        (cx, y0_bot),  (cx, y1_bot), (x1, y1_bot),   # bezier to bottom-right
        (x1, y1_top),                          # top-right
        (cx, y1_top),  (cx, y0_top), (x0, y0_top),   # bezier to top-left
        (x0, y0_bot),                          # close
    ]
    path  = Path(verts, codes)
    patch = PathPatch(path, facecolor=color, edgecolor='none', alpha=alpha)
    ax_s.add_patch(patch)

# Iterate mouse celltypes × C4 clusters
for ct in MOUSE_CELLTYPES:
    col = MOUSE_COLORS[ct]
    for cid in c4_order:
        n = xtab.loc[ct, cid] if ct in xtab.index else 0
        if n == 0:
            continue
        frac = n / grand_total
        # Left node: fill downward from current position
        y0_top = ct_fill[ct]
        y0_bot = y0_top - frac
        ct_fill[ct] = y0_bot
        # Right node: fill downward
        y1_top = c4_fill[cid]
        y1_bot = y1_top - frac
        c4_fill[cid] = y1_bot
        bezier_band(X_L, y0_bot, y0_top, X_R, y1_bot, y1_top,
                    color=col, alpha=0.45)

# ── 5. Draw node bars ─────────────────────────────────────────────────────────
BAR_W = 0.04

# Left nodes — mouse cell types
for ct in MOUSE_CELLTYPES:
    yb = ct_y_bottom[ct]
    yt = ct_y_top[ct]
    ax_s.add_patch(mpatches.FancyArrow(
        0, 0, 0, 0, width=0))   # dummy; use Rectangle instead
    ax_s.add_patch(plt.Rectangle(
        (X_L - BAR_W/2, yb), BAR_W, yt - yb,
        facecolor=MOUSE_COLORS[ct], edgecolor='white', lw=0.5, zorder=4
    ))
    pct = ct_totals.get(ct, 0) / grand_total * 100
    ax_s.text(X_L - BAR_W/2 - 0.02, (yb+yt)/2,
              f'{ct.replace("_"," ")}\n({pct:.1f}%)',
              ha='right', va='center', fontsize=10, fontweight='bold',
              color=MOUSE_COLORS[ct])

# Right nodes — C4 clusters
for cid in c4_order:
    yb  = c4_y_bottom[cid]
    yt  = c4_y_top[cid]
    col = (C4_GROUP_COLORS['Group A'] if cid in set(GROUP_A_C4)
           else C4_GROUP_COLORS['Group B'])
    ax_s.add_patch(plt.Rectangle(
        (X_R - BAR_W/2, yb), BAR_W, yt - yb,
        facecolor=col, edgecolor='white', lw=0.5, zorder=4
    ))
    pct = c4_totals[cid] / grand_total * 100
    ax_s.text(X_R + BAR_W/2 + 0.02, (yb+yt)/2,
              f'{cid}  {pct:.1f}%',
              ha='left', va='center', fontsize=10,
              color=col, fontweight='bold' if cid in set(GROUP_A_C4) else 'normal')

ax_s.set_xlim(-0.55, 1.35)
ax_s.set_ylim(-0.02, 1.02)
ax_s.axis('off')
ax_s.set_title('Mouse Lmx1a neuron → human C4 label assignment',
               fontweight='bold', fontsize=11)

# Legend: mouse cell types
leg_handles = [mpatches.Patch(fc=MOUSE_COLORS[ct], label=ct.replace('_',' '))
               for ct in MOUSE_CELLTYPES]
leg_handles += [
    mpatches.Patch(fc=C4_GROUP_COLORS['Group A'], label='Group A  (LMX1A+)'),
    mpatches.Patch(fc=C4_GROUP_COLORS['Group B'], label='Group B  (LMX1A-)'),
]
#ax_s.legend(handles=leg_handles, fontsize=7.5,
#            loc='lower left', bbox_to_anchor=(0, 0), framealpha=0.9,
#            title='Mouse cell type / C4 group', title_fontsize=7.5)

# ── 6. Bar chart: human LMX1A expression ─────────────────────────────────────
THIN_BAR    = 0.012
x_max_data  = max(bar_vals) * 1.15   # xlim ceiling

# Group B values are ~0.001–0.005 — invisible on this scale.
# Show them as a hollow diamond marker at x=0 with a "≈0" label instead.
NEAR_ZERO   = 0.05 * max(bar_vals)   # threshold for "invisible" bar

grp_a_nodes = [c for c in c4_order if c in set(GROUP_A_C4)]
grp_b_nodes = [c for c in c4_order if c in set(GROUP_B_C4)]

# Group A bars (genuine values)
for cid, yc, val, col in zip(c4_order, node_centres, bar_vals, bar_cols):
    if val >= NEAR_ZERO:
        ax_b.barh(yc, val, height=THIN_BAR, color=col,
                  edgecolor='white', linewidth=0.3)
    else:
        # Near-zero: diamond marker
        ax_b.scatter(0, yc, marker='D', s=18,
                     color=C4_GROUP_COLORS['Group B'],
                     edgecolors='white', linewidths=0.3, zorder=4)

ax_b.axvline(0, color='black', lw=0.5)

# ── Y limits flush with Sankey ────────────────────────────────────────────────
y_top_all = c4_y_top[c4_order[0]]
y_bot_all = c4_y_bottom[c4_order[-1]]
ax_b.set_ylim(y_bot_all - 0.005, y_top_all + 0.005)

# ── Group brackets — pure data coordinates, tight to bars ────────────────────
grp_a_y_top = c4_y_top[grp_a_nodes[0]]
grp_a_y_bot = c4_y_bottom[grp_a_nodes[-1]]
grp_b_y_top = c4_y_top[grp_b_nodes[0]]
grp_b_y_bot = c4_y_bottom[grp_b_nodes[-1]]

bracket_x    = x_max_data * 1.01   # just past xlim
bracket_arm  = x_max_data * 0.03   # horizontal arm length
text_x       = bracket_x + bracket_arm + x_max_data * 0.01

for y_top, y_bot, label, col in [
    (grp_a_y_top, grp_a_y_bot, 'Group A\n(LMX1A+)', C4_GROUP_COLORS['Group A']),
    (grp_b_y_top, grp_b_y_bot, 'Group B\n(LMX1A−)', C4_GROUP_COLORS['Group B']),
]:
    y_mid = (y_top + y_bot) / 2
    # Bracket: [ shape in data coords
    ax_b.plot([bracket_x, bracket_x + bracket_arm,
               bracket_x + bracket_arm, bracket_x],
              [y_top,     y_top,
               y_bot,     y_bot],
              color=col, lw=1.2, clip_on=False)
    ax_b.text(text_x, y_mid, label,
              ha='left', va='center', fontsize=11, fontweight='bold',
              color=col, clip_on=False)

ax_b.set_xlim(0, x_max_data)
ax_b.set_xlabel('Mean normalized LMX1A expression', fontsize=11)
ax_b.set_yticks([])
ax_b.tick_params(axis='x', labelsize=11)
ax_b.set_title('Human LMX1A expression', fontweight='bold', fontsize=11)
ax_b.spines['top'].set_visible(False)
ax_b.spines['right'].set_visible(False)
ax_b.spines['left'].set_visible(False)

plt.savefig('LMX_C4_sankey_LMX1A.pdf', bbox_inches='tight')
plt.savefig('LMX_C4_sankey_LMX1A.png', bbox_inches='tight')
plt.show()


# In[10]:


# Add significance markers to the summary table
cellect = pd.read_csv("hypomap_cellect_unified.csv")
bonf = 0.05/415
print("Lmx1a C4 clusters — paper_MAGMA significance (corrected mapping):")
print(f"Bonferroni threshold: {bonf:.2e}")
paper_lmx = cellect[
    (cellect['method']=='paper_MAGMA') &
    (cellect['cluster_id'].isin(LMX_C4_IDS)) &
    (cellect['is_neuronal'])
].sort_values('pvalue')
for _, r in paper_lmx.iterrows():
    sig = '***' if r['pvalue'] < bonf else '(nominal)' if r['pvalue'] < 0.05 else ''
    lmx1a_pct = df_combined_full.loc[
        df_combined_full['C4_id']==r['cluster_id'], 'pct_lmx1a'].values
    lmx1a_pct = f"{lmx1a_pct[0]:.0f}%" if len(lmx1a_pct) > 0 else 'N/A'
    print(f"  {r['cluster_id']:<10}  p={r['pvalue']:.3e}  {sig:<12}  "
          f"mouse_spec={lmx1a_pct}")


# Group A — TRUE Lmx1a (high mouse specificity + high human LMX1A + paper Bonf-sig):
#   C4-229: 97% specific, 90% LMX1A+, p=3.1e-08 ***
#   C4-252: 88% specific, 88% LMX1A+, p=1.0e-06 ***
#   C4-247: 99% specific, 86% LMX1A+, p=1.3e-06 ***
#   C4-224: 91% specific, 89% LMX1A+, p=4.5e-06 ***
#   C4-223: 97% specific, 82% LMX1A+, p=4.5e-06 ***
#   C4-228: 61% specific, 64% LMX1A+, p=4.6e-06 ***  ← intermediate
#   C4-230: 59% specific, 78% LMX1A+, p=4.9e-06 ***  ← intermediate
#   C4-226: 100% specific, 80% LMX1A+, p=1.4e-04 (marginal)
#   C4-231: 100% specific, 88% LMX1A+, p=7.2e-04 (nominal)
#   C4-248: 98% specific, 90% LMX1A+, p=1.3e-03 (nominal)
# 
# Group B — CONTAMINANT labels (low specificity + near-zero LMX1A):
#   C4-411: 16% specific, 0.8% LMX1A+, p=3.8e-04 (nominal)
#   C4-407: 22% specific, 0.5% LMX1A+, p=6.1e-05 ***  ← sig but biologically wrong
#   C4-409: 22% specific, 0.7% LMX1A+, p=1.6e-02 (nominal)
#   C4-33:  28% specific, 0.1% LMX1A+, p=9.4e-03 (nominal)
#   C4-7:    2% specific, 0.3% LMX1A+, p=3.2e-02 (nominal) ← REMOVE

# In[16]:


# ── Analysis: C3 label specificity + LMX1A expression + BMI overlap ──────────
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from scipy.sparse import issparse

MOUSE_CELLTYPES = ['04_Crnde_Lmx1a', '27_Lmx1a_Ttc6',
                   '32_Gucy2c_Ttc6', '36_Ebf2_Crhbp']

LMX_C3_IDS = ['C3-92','C3-84','C3-139','C3-29','C3-94','C3-83','C3-82','C3-26',
               'C3-137','C3-85','C3-28','C3-24','C3-93','C3-88','C3-115','C3-95']

# ── 0. Prepare obs frames ─────────────────────────────────────────────────────
obs_mouse = adata_mouse.obs.copy()
if 'C3_short' not in obs_mouse.columns:
    if 'C3_HM' in obs_mouse.columns:
        obs_mouse['C3_short'] = obs_mouse['C3_HM'].str.extract(r'(C3-\d+)')
    else:
        obs_mouse['C3_short'] = pd.NA
        print("WARNING: C3_HM not found in mouse obs")

print(f"Mouse C3_short non-null: {obs_mouse['C3_short'].notna().sum()}")
obs_mouse['is_lmx1a_pop'] = obs_mouse['detailed_celltype_annotation'].isin(MOUSE_CELLTYPES)

obs_human = adata.obs.copy()
if 'C3_short' not in obs_human.columns:
    obs_human['C3_short'] = obs_human['C3'].str.extract(r'(C3-\d+)')
print(f"Human C3_short non-null: {obs_human['C3_short'].notna().sum()}")

# ── 1. Load C3 BMI results from unified CSV ───────────────────────────────────
# Paper has no CELLEX-based C3 — use new_MAGMA as primary
cellect = pd.read_csv("hypomap_cellect_unified.csv")
GWAS = 'BMI_Yengo2018'
BONF_C3 = 0.05 / 156

c3_magma = (cellect[(cellect['gwas']==GWAS) & (cellect['method']=='new_MAGMA') &
                    (cellect['level']=='C3') & cellect['is_neuronal']]
            .copy())
c3_ldsc  = (cellect[(cellect['gwas']==GWAS) & (cellect['method']=='new_LDSC') &
                    (cellect['level']=='C3') & cellect['is_neuronal']]
            .copy())

# Build flat wide table for C3
cellect_wide_c3 = (c3_magma[['cluster_id','has_exome_sig','n_effector_genes',
                               'pvalue','log10p','bonf_sig']]
                   .rename(columns={'pvalue':    'pvalue_new_MAGMA',
                                     'log10p':    'log10p_new_MAGMA',
                                     'bonf_sig':  'bonf_sig_new_MAGMA'})
                   .merge(
                       c3_ldsc[['cluster_id','pvalue','log10p','bonf_sig']]
                       .rename(columns={'pvalue':   'pvalue_new_LDSC',
                                         'log10p':   'log10p_new_LDSC',
                                         'bonf_sig': 'bonf_sig_new_LDSC'}),
                       on='cluster_id', how='left')
                   .rename(columns={'cluster_id': 'C3_id'}))

print(f"\nnew_MAGMA C3 Bonf-sig: {cellect_wide_c3['bonf_sig_new_MAGMA'].sum()} / 156")
print(f"new_LDSC  C3 Bonf-sig: {cellect_wide_c3['bonf_sig_new_LDSC'].sum()}  / 156")

# ── 2. ANALYSIS 1: Mouse label specificity at C3 ─────────────────────────────
print("\n" + "=" * 70)
print("ANALYSIS 1: C3 label specificity to mouse Lmx1a populations")
print("=" * 70)

spec_rows = []
for c3id in LMX_C3_IDS:
    cells = obs_mouse[obs_mouse['C3_short'] == c3id]
    n_total  = len(cells)
    n_lmx1a  = cells['is_lmx1a_pop'].sum()
    pct      = n_lmx1a / n_total * 100 if n_total > 0 else 0.0
    # Annotation: what mouse cell types receive this label?
    top_ct = (cells['detailed_celltype_annotation']
              .value_counts().head(3).index.tolist())
    spec_rows.append({
        'C3_id':     c3id,
        'n_total':   n_total,
        'n_lmx1a':   n_lmx1a,
        'n_other':   n_total - n_lmx1a,
        'pct_lmx1a': pct,
        'specific':  pct >= 50,
        'top_celltypes': ' | '.join(top_ct),
    })

df_spec_c3 = pd.DataFrame(spec_rows).sort_values('pct_lmx1a', ascending=False)
print(df_spec_c3[['C3_id','n_total','n_lmx1a','pct_lmx1a',
                   'specific','top_celltypes']].to_string(index=False))

non_lmx_c3 = obs_mouse[
    obs_mouse['C3_short'].isin(LMX_C3_IDS) & ~obs_mouse['is_lmx1a_pop']
]
print("\nNon-Lmx1a cell types receiving Lmx1a C3 labels (top 20):")
print(non_lmx_c3['detailed_celltype_annotation'].value_counts().head(20).to_string())

# ── 3. ANALYSIS 2: LMX1A expression per C3 cluster in human data ─────────────
print("\n" + "=" * 70)
print("ANALYSIS 2: LMX1A expression across human C3 clusters")
print("=" * 70)

human_gene = None
for candidate in ['LMX1A', 'Lmx1a', 'lmx1a']:
    if candidate in adata.var_names:
        human_gene = candidate; break
if human_gene is None:
    matches = [g for g in adata.var_names if g.upper() == 'LMX1A']
    if matches: human_gene = matches[0]

if human_gene:
    print(f"Using gene: {human_gene}")
    gene_idx   = list(adata.var_names).index(human_gene)
    X_norm     = adata.layers['norm']
    lmx1a_expr = (np.asarray(X_norm[:, gene_idx].todense()).ravel()
                  if issparse(X_norm)
                  else np.asarray(X_norm[:, gene_idx]).ravel())
    obs_human['LMX1A_expr'] = lmx1a_expr

    c3_lmx1a = (obs_human.groupby('C3_short', observed=True)['LMX1A_expr']
                .agg(mean_expr='mean',
                     pct_expr=lambda x: (x > 0).mean() * 100,
                     n_cells='count')
                .reset_index()
                .sort_values('mean_expr', ascending=False)
                .reset_index(drop=True))
    c3_lmx1a['rank_pos'] = range(len(c3_lmx1a))

    expr_thresh = 0.5   # same threshold as C4 analysis
    pct_thresh  = 10.0
    c3_lmx1a['high_LMX1A'] = ((c3_lmx1a['mean_expr'] > expr_thresh) &
                                (c3_lmx1a['pct_expr']  > pct_thresh))

    c3_lmx1a_focus = c3_lmx1a[c3_lmx1a['C3_short'].isin(LMX_C3_IDS)]
    print(f"\nLMX1A threshold: mean>{expr_thresh}  pct>{pct_thresh}%")
    print(f"C3 clusters with high LMX1A: {c3_lmx1a['high_LMX1A'].sum()}")
    print("\nLMX1A in Lmx1a C3 clusters (human):")
    print(c3_lmx1a_focus.sort_values('mean_expr', ascending=False)
          [['C3_short','mean_expr','pct_expr','n_cells','high_LMX1A']]
          .to_string(index=False))
else:
    print("LMX1A not found")
    obs_human['LMX1A_expr'] = 0.0
    c3_lmx1a = pd.DataFrame(columns=['C3_short','mean_expr','pct_expr',
                                      'n_cells','rank_pos','high_LMX1A'])
    expr_thresh = 0.0

# ── 4. ANALYSIS 3: Combined overlap table ────────────────────────────────────
print("\n" + "=" * 70)
print("ANALYSIS 3: Label specificity ∩ LMX1A expression ∩ BMI enrichment")
print("=" * 70)

overlap_rows = []
for c3id in LMX_C3_IDS:
    hcells = obs_human[obs_human['C3_short'] == c3id]
    if len(hcells) == 0:
        overlap_rows.append({'C3_id':c3id,'n_human':0,
                              'mean_expr':np.nan,'pct_expr':np.nan})
        continue
    expr = hcells['LMX1A_expr']
    overlap_rows.append({'C3_id':c3id,'n_human':len(hcells),
                          'mean_expr':expr.mean(),
                          'pct_expr':(expr>0).mean()*100})

df_overlap_c3 = pd.DataFrame(overlap_rows)
df_combined_c3 = (df_spec_c3
                  .merge(df_overlap_c3[['C3_id','n_human',
                                         'mean_expr','pct_expr']],
                         on='C3_id', how='left')
                  .rename(columns={'mean_expr':'human_LMX1A_mean',
                                   'pct_expr': 'human_LMX1A_pct'}))

df_combined_full_c3 = (df_combined_c3
                       .merge(cellect_wide_c3, on='C3_id', how='left'))

# Print clean combined table
print(f"\n{'C3_id':<10} {'%mouse':>7} {'spec':>5} {'LMX1A_mean':>11} "
      f"{'LMX1A_pct':>10} {'p_MAGMA':>10} {'Bonf_M':>7} "
      f"{'p_LDSC':>10} {'Bonf_L':>7}")
print("-" * 90)
for _, r in df_combined_full_c3.sort_values('pct_lmx1a', ascending=False).iterrows():
    spec  = '✓' if r['specific']  else ' '
    hmean = f"{r['human_LMX1A_mean']:.3f}" if pd.notna(r['human_LMX1A_mean']) else 'N/A'
    hpct  = f"{r['human_LMX1A_pct']:.1f}%"  if pd.notna(r['human_LMX1A_pct'])  else 'N/A'
    pm    = f"{r['pvalue_new_MAGMA']:.2e}"   if pd.notna(r.get('pvalue_new_MAGMA')) else 'N/A'
    bm    = '***' if r.get('bonf_sig_new_MAGMA', False) else ''
    pl    = f"{r['pvalue_new_LDSC']:.2e}"    if pd.notna(r.get('pvalue_new_LDSC'))  else 'N/A'
    bl    = '***' if r.get('bonf_sig_new_LDSC', False) else ''
    print(f"{r['C3_id']:<10} {r['pct_lmx1a']:>6.1f}% {spec:>5} "
          f"{hmean:>11} {hpct:>10} {pm:>10} {bm:>7} {pl:>10} {bl:>7}")

# ── 5. Figures ────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(16, 13))
COL_HL    = '#8B0000'
COL_LMX   = '#4878CF'

# ── Panel A: Mouse label specificity bar chart ────────────────────────────────
ax   = axes[0, 0]
df_sp = df_spec_c3.sort_values('pct_lmx1a', ascending=True)
colors_sp = [COL_HL if c in {'C3-84','C3-92'} else COL_LMX
             for c in df_sp['C3_id']]
ax.barh(range(len(df_sp)), df_sp['pct_lmx1a'],
        color=colors_sp, edgecolor='white', linewidth=0.4, height=0.72)
ax.axvline(50, color='black', ls='--', lw=1, label='50% threshold')
for i, (_, r) in enumerate(df_sp.iterrows()):
    ax.text(r['pct_lmx1a']+0.5, i,
            f"{int(r['n_lmx1a'])}/{int(r['n_total'])}",
            va='center', fontsize=6.5, color='grey')
ax.set_yticks(range(len(df_sp)))
ax.set_yticklabels(df_sp['C3_id'], fontsize=8)
for lbl, c3 in zip(ax.get_yticklabels(), df_sp['C3_id']):
    if c3 in {'C3-84','C3-92'}:
        lbl.set_color(COL_HL); lbl.set_fontweight('bold')
ax.set_xlabel('% cells with C3 label that are Lmx1a mouse population')
ax.set_title('C3 label specificity to\nmouse Lmx1a neurons', fontweight='bold')
ax.set_xlim(0, 110); ax.legend(fontsize=8)

# ── Panel B: Non-Lmx1a contamination ─────────────────────────────────────────
ax = axes[0, 1]
top_contam = non_lmx_c3['detailed_celltype_annotation'].value_counts().head(15)
ax.barh(range(len(top_contam)), top_contam.values,
        color='#95a5a6', edgecolor='white', linewidth=0.4, height=0.72)
ax.set_yticks(range(len(top_contam)))
ax.set_yticklabels(top_contam.index, fontsize=7.5)
ax.set_xlabel('Number of cells')
ax.set_title('Non-Lmx1a mouse cell types\nalso receiving Lmx1a C3 labels',
             fontweight='bold')
for i, (_, n) in enumerate(zip(top_contam.index, top_contam.values)):
    ax.text(n+1, i, str(n), va='center', fontsize=7, color='grey')

# ── Panel C: LMX1A expression rank across all C3 clusters ────────────────────
ax = axes[1, 0]
if len(c3_lmx1a) > 0:
    ax.scatter(c3_lmx1a['rank_pos'], c3_lmx1a['mean_expr'],
               c='#cccccc', s=15, alpha=0.5, zorder=2,
               label=f'All C3 (n={len(c3_lmx1a)})')
    for _, row in c3_lmx1a[c3_lmx1a['C3_short'].isin(LMX_C3_IDS)].iterrows():
        col = COL_HL if row['C3_short'] in {'C3-84','C3-92','C3-95'} else '#E07B54'
        ax.scatter(row['rank_pos'], row['mean_expr'],
                   c=col, s=70, zorder=4,
                   edgecolors='white', linewidths=0.5)
        ax.annotate(row['C3_short'],
                    xy=(row['rank_pos'], row['mean_expr']),
                    xytext=(4,3), textcoords='offset points',
                    fontsize=6, color=col)
    ax.axhline(expr_thresh, color='black', ls='--', lw=1,
               label=f'Threshold ({expr_thresh})')

ax.set_xlabel('C3 cluster rank (by LMX1A expression, 1=highest)')
ax.set_ylabel('Mean LMX1A expression (norm layer)')
ax.set_title('LMX1A expression across human\nhypomap C3 clusters', fontweight='bold')
ax.legend(handles=[
    mpatches.Patch(fc=COL_HL,    label='C3-84/92/95 (top Lmx1a)'),
    mpatches.Patch(fc='#E07B54', label='Other Lmx1a C3s'),
    mpatches.Patch(fc='#cccccc', label='All other C3 clusters'),
], fontsize=7.5)

# ── Panel D: Specificity vs LMX1A coloured by BMI ────────────────────────────
ax = axes[1, 1]
df_combined_full_c3['log10p_pm'] = df_combined_full_c3['log10p_new_MAGMA'].fillna(0)

if df_combined_full_c3['human_LMX1A_mean'].notna().any():
    sc = ax.scatter(
        df_combined_full_c3['pct_lmx1a'],
        df_combined_full_c3['human_LMX1A_mean'],
        c=df_combined_full_c3['log10p_pm'],
        cmap='YlOrRd',
        vmin=0,
        vmax=df_combined_full_c3['log10p_pm'].quantile(0.95)+0.01,
        s=[max(40, n/20) for n in df_combined_full_c3['n_human'].fillna(0)],
        edgecolors=['black' if c in {'C3-84','C3-92','C3-95'} else 'white'
                    for c in df_combined_full_c3['C3_id']],
        linewidths=[1.2 if c in {'C3-84','C3-92','C3-95'} else 0.4
                    for c in df_combined_full_c3['C3_id']],
        zorder=3, alpha=0.9
    )
    plt.colorbar(sc, ax=ax, label='−log₁₀(p)  new_MAGMA (Yengo2018)')

    ax.axvline(50,          color='black', ls='--', lw=1, alpha=0.5,
               label='50% mouse specificity')
    ax.axhline(expr_thresh, color='grey',  ls='--', lw=1, alpha=0.5,
               label=f'LMX1A threshold ({expr_thresh})')

    for _, r in df_combined_full_c3.iterrows():
        if pd.notna(r['human_LMX1A_mean']):
            ax.annotate(r['C3_id'],
                        xy=(r['pct_lmx1a'], r['human_LMX1A_mean']),
                        xytext=(4,2), textcoords='offset points',
                        fontsize=6.5,
                        color=COL_HL if r['C3_id'] in {'C3-84','C3-92','C3-95'}
                              else '#333333')

    xlim = ax.get_xlim(); ylim = ax.get_ylim(); yspan = ylim[1]-ylim[0]
    ax.text(55,  ylim[1]-yspan*0.05, 'Specific + LMX1A high',
            fontsize=7.5, color='#2ecc71', fontweight='bold', va='top')
    ax.text(2,   ylim[1]-yspan*0.05, 'Non-specific + LMX1A high',
            fontsize=7.5, color='#e67e22', fontweight='bold', va='top')
    ax.text(55,  ylim[0]+yspan*0.02, 'Specific + LMX1A low',
            fontsize=7.5, color='#3498db', fontweight='bold')

    ax.legend(fontsize=7.5, loc='upper left')

ax.set_xlabel('% cells with C3 label from Lmx1a mouse population', fontsize=9)
ax.set_ylabel('Mean LMX1A expression in human cluster', fontsize=9)
ax.set_title('C3 label specificity vs human LMX1A expression\n'
             '(colour = BMI enrichment, bubble = n human cells)',
             fontweight='bold')

plt.suptitle('Mouse Lmx1a C3 label specificity × human LMX1A expression × BMI enrichment\n'
             'BMI: new_MAGMA (paper has no CELLEX-based C3)',
             fontweight='bold', fontsize=13, y=1.01)
plt.tight_layout()
plt.savefig('LMX_C3_label_specificity_human_expression.pdf', bbox_inches='tight')
plt.savefig('LMX_C3_label_specificity_human_expression.png', bbox_inches='tight')
plt.show()

# ── 6. Clean Group A/B assignment based on BOTH criteria ─────────────────────
print("\n" + "=" * 80)
print("C3 GROUP ASSIGNMENT: specificity ≥50% AND human LMX1A mean >0.5")
print("=" * 80)
for _, r in df_combined_full_c3.sort_values('pct_lmx1a', ascending=False).iterrows():
    spec_ok = r['pct_lmx1a'] >= 50
    lmx_ok  = pd.notna(r['human_LMX1A_mean']) and r['human_LMX1A_mean'] > expr_thresh
    if spec_ok and lmx_ok:
        grp = 'Group A  (genuine LMX1A+)'
    elif spec_ok:
        grp = 'Group A-low  (specific, LMX1A low in human)'
    else:
        grp = 'Group B  (contaminant / non-LMX1A)'
    bm    = ('***M' if r.get('bonf_sig_new_MAGMA', False)
             else 'nomM' if r.get('pvalue_new_MAGMA', 1) < 0.05
             else '')
    bl    = ('***L' if r.get('bonf_sig_new_LDSC', False)
             else 'nomL' if r.get('pvalue_new_LDSC', 1) < 0.05
             else '')
    hmean = f"{r['human_LMX1A_mean']:.3f}" if pd.notna(r['human_LMX1A_mean']) else 'N/A'
    print(f"  {r['C3_id']:<10}  spec={r['pct_lmx1a']:>5.1f}%  "
          f"LMX1A={hmean:>6}  "
          f"BMI: {bm:<6} {bl:<6}  → {grp}")


# In[17]:


# ── Sankey + LMX1A bar — C3 label assignment from mouse Lmx1a neurons ─────────
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import PathPatch
from matplotlib.path import Path

# ── 1. Compute cross-tabulation from mouse obs ────────────────────────────────
MOUSE_CELLTYPES = ['04_Crnde_Lmx1a','27_Lmx1a_Ttc6',
                   '32_Gucy2c_Ttc6','36_Ebf2_Crhbp']
MOUSE_COLORS    = {
    '04_Crnde_Lmx1a': '#008941',
    '27_Lmx1a_Ttc6':  '#ffaa92',
    '32_Gucy2c_Ttc6': '#000035',
    '36_Ebf2_Crhbp':  '#0aa6d8',
}

# C3-level group definitions
# Group A (LMX1A+): genuine Group A + CRNDE (both express LMX1A in human)
GRP_A_C3_PLOT = list(GROUP_A_C3) + list(CRNDE_C3)
# Group B (LMX1A-): GABA-neighbour + contaminants
GRP_B_C3_PLOT = list(GABA_NEIGHBOR_C3) + list(GROUP_B_C3)
ALL_C3_PLOT   = GRP_A_C3_PLOT + GRP_B_C3_PLOT

# Reuse same colour scheme as C4 plot
C3_GROUP_PLOT_COLORS = {
    'Group A': C3_GROUP_COLORS['Group A'],   # '#8B0000'
    'Group B': '#95a5a6',
}

if 'C3_short' not in adata_mouse.obs.columns:
    adata_mouse.obs['C3_short'] = (adata_mouse.obs['C3_HM']
                                   .str.extract(r'(C3-\d+)'))
    print(f"C3_short extracted: {adata_mouse.obs['C3_short'].notna().sum()} non-null")

lmx_obs = adata_mouse.obs[
    adata_mouse.obs['detailed_celltype_annotation'].isin(MOUSE_CELLTYPES)
].copy()
lmx_obs['C3_short'] = lmx_obs['C3_short'].astype(str)

xtab = (lmx_obs[lmx_obs['C3_short'].isin(ALL_C3_PLOT)]
        .groupby(['detailed_celltype_annotation','C3_short'])
        .size()
        .unstack(fill_value=0))

for c in ALL_C3_PLOT:
    if c not in xtab.columns:
        xtab[c] = 0

c3_totals = xtab.sum(axis=0)
c3_order  = (sorted(GRP_A_C3_PLOT, key=lambda c: -c3_totals.get(c, 0)) +
              sorted(GRP_B_C3_PLOT, key=lambda c: -c3_totals.get(c, 0)))
c3_order  = [c for c in c3_order if c in xtab.columns]
xtab      = xtab[c3_order]

lmx1a_mean = {}
for cid in c3_order:
    mask = obs_human['C3_short'].astype(str) == cid
    if mask.sum() > 0:
        lmx1a_mean[cid] = obs_human.loc[mask, 'LMX1A_expr'].mean()
    else:
        lmx1a_mean[cid] = 0.0

# ── 2. Layout constants ───────────────────────────────────────────────────────
fig, (ax_s, ax_b) = plt.subplots(
    1, 2, figsize=(16, 8),
    gridspec_kw={'width_ratios': [3.5, 0.6], 'wspace': 0.01}
)

X_L, X_R = 0.0, 1.0

# ── 3. Node y-positions ───────────────────────────────────────────────────────
ct_totals   = xtab.sum(axis=1)
grand_total = ct_totals.sum()
ct_heights  = ct_totals / grand_total

ct_y_bottom = {}
y = 1.0
for ct in MOUSE_CELLTYPES:
    h = ct_heights.get(ct, 0)
    ct_y_bottom[ct] = y - h
    y -= h
ct_y_top = {ct: ct_y_bottom[ct] + ct_heights.get(ct, 0) for ct in MOUSE_CELLTYPES}

GAP          = 0.015
c3_heights   = c3_totals[c3_order] / grand_total
y            = 1.0
c3_y_bottom  = {}
for i, cid in enumerate(c3_order):
    if i > 0 and cid == GRP_B_C3_PLOT[0] and c3_order[i-1] in GRP_A_C3_PLOT:
        y -= GAP
    h = c3_heights[cid]
    c3_y_bottom[cid] = y - h
    y -= h
c3_y_top = {cid: c3_y_bottom[cid] + c3_heights[cid] for cid in c3_order}

# ── 4. Draw Sankey flows ──────────────────────────────────────────────────────
ct_fill = {ct:  ct_y_top[ct]  for ct in MOUSE_CELLTYPES}
c3_fill = {cid: c3_y_top[cid] for cid in c3_order}

def bezier_band(x0, y0_bot, y0_top, x1, y1_bot, y1_top, color, alpha=0.45):
    cx    = (x0 + x1) / 2
    verts = [(x0, y0_bot), (cx, y0_bot), (cx, y1_bot), (x1, y1_bot),
             (x1, y1_top), (cx, y1_top), (cx, y0_top), (x0, y0_top),
             (x0, y0_bot)]
    codes = [Path.MOVETO,
             Path.CURVE4, Path.CURVE4, Path.CURVE4,
             Path.LINETO,
             Path.CURVE4, Path.CURVE4, Path.CURVE4,
             Path.CLOSEPOLY]
    ax_s.add_patch(PathPatch(Path(verts, codes),
                             facecolor=color, edgecolor='none', alpha=alpha))

for ct in MOUSE_CELLTYPES:
    col = MOUSE_COLORS[ct]
    for cid in c3_order:
        n = xtab.loc[ct, cid] if ct in xtab.index else 0
        if n == 0: continue
        frac   = n / grand_total
        y0_top = ct_fill[ct];   y0_bot = y0_top - frac; ct_fill[ct] = y0_bot
        y1_top = c3_fill[cid];  y1_bot = y1_top - frac; c3_fill[cid] = y1_bot
        bezier_band(X_L, y0_bot, y0_top, X_R, y1_bot, y1_top,
                    color=col, alpha=0.45)

# ── 5. Draw node bars ─────────────────────────────────────────────────────────
BAR_W = 0.04

for ct in MOUSE_CELLTYPES:
    yb = ct_y_bottom[ct]; yt = ct_y_top[ct]
    ax_s.add_patch(plt.Rectangle(
        (X_L - BAR_W/2, yb), BAR_W, yt - yb,
        facecolor=MOUSE_COLORS[ct], edgecolor='white', lw=0.5, zorder=4))
    pct = ct_totals.get(ct, 0) / grand_total * 100
    ax_s.text(X_L - BAR_W/2 - 0.02, (yb+yt)/2,
              f'{ct.replace("_"," ")}\n({pct:.1f}%)',
              ha='right', va='center', fontsize=10, fontweight='bold',
              color=MOUSE_COLORS[ct])

for cid in c3_order:
    yb  = c3_y_bottom[cid]; yt = c3_y_top[cid]
    col = (C3_GROUP_PLOT_COLORS['Group A'] if cid in set(GRP_A_C3_PLOT)
           else C3_GROUP_PLOT_COLORS['Group B'])
    ax_s.add_patch(plt.Rectangle(
        (X_R - BAR_W/2, yb), BAR_W, yt - yb,
        facecolor=col, edgecolor='white', lw=0.5, zorder=4))
    pct = c3_totals[cid] / grand_total * 100
    ax_s.text(X_R + BAR_W/2 + 0.02, (yb+yt)/2,
              f'{cid}  {pct:.1f}%',
              ha='left', va='center', fontsize=10,
              color=col,
              fontweight='bold' if cid in set(GRP_A_C3_PLOT) else 'normal')

ax_s.set_xlim(-0.55, 1.35)
ax_s.set_ylim(-0.02, 1.02)
ax_s.axis('off')
ax_s.set_title('Mouse Lmx1a neuron → human C3 label assignment',
               fontweight='bold', fontsize=11)

# ── 6. Bar chart: human LMX1A expression ─────────────────────────────────────
node_centres = [(c3_y_bottom[c] + c3_y_top[c]) / 2 for c in c3_order]
bar_vals     = [lmx1a_mean.get(c, 0) for c in c3_order]
bar_cols     = [C3_GROUP_PLOT_COLORS['Group A'] if c in set(GRP_A_C3_PLOT)
                else C3_GROUP_PLOT_COLORS['Group B']
                for c in c3_order]

THIN_BAR   = 0.012
x_max_data = max(bar_vals) * 1.15
NEAR_ZERO  = 0.05 * max(bar_vals)

grp_a_nodes = [c for c in c3_order if c in set(GRP_A_C3_PLOT)]
grp_b_nodes = [c for c in c3_order if c in set(GRP_B_C3_PLOT)]

for cid, yc, val, col in zip(c3_order, node_centres, bar_vals, bar_cols):
    if val >= NEAR_ZERO:
        ax_b.barh(yc, val, height=THIN_BAR, color=col,
                  edgecolor='white', linewidth=0.3)
    else:
        ax_b.scatter(0, yc, marker='D', s=18,
                     color=C3_GROUP_PLOT_COLORS['Group B'],
                     edgecolors='white', linewidths=0.3, zorder=4)

ax_b.axvline(0, color='black', lw=0.5)

y_top_all = c3_y_top[c3_order[0]]
y_bot_all = c3_y_bottom[c3_order[-1]]
ax_b.set_ylim(y_bot_all - 0.005, y_top_all + 0.005)

grp_a_y_top  = c3_y_top[grp_a_nodes[0]]
grp_a_y_bot  = c3_y_bottom[grp_a_nodes[-1]]
grp_b_y_top  = c3_y_top[grp_b_nodes[0]]
grp_b_y_bot  = c3_y_bottom[grp_b_nodes[-1]]
bracket_x    = x_max_data * 1.01
bracket_arm  = x_max_data * 0.03
text_x       = bracket_x + bracket_arm + x_max_data * 0.01

for y_top, y_bot, label, col in [
    (grp_a_y_top, grp_a_y_bot, 'Group A\n(LMX1A+)', C3_GROUP_PLOT_COLORS['Group A']),
    (grp_b_y_top, grp_b_y_bot, 'Group B\n(LMX1A−)', C3_GROUP_PLOT_COLORS['Group B']),
]:
    y_mid = (y_top + y_bot) / 2
    ax_b.plot([bracket_x, bracket_x + bracket_arm,
               bracket_x + bracket_arm, bracket_x],
              [y_top, y_top, y_bot, y_bot],
              color=col, lw=1.2, clip_on=False)
    ax_b.text(text_x, y_mid, label,
              ha='left', va='center', fontsize=11, fontweight='bold',
              color=col, clip_on=False)

ax_b.set_xlim(0, x_max_data)
ax_b.set_xlabel('Mean normalized LMX1A expression', fontsize=11)
ax_b.set_yticks([])
ax_b.tick_params(axis='x', labelsize=11)
ax_b.set_title('Human LMX1A expression', fontweight='bold', fontsize=11)
ax_b.spines['top'].set_visible(False)
ax_b.spines['right'].set_visible(False)
ax_b.spines['left'].set_visible(False)

plt.savefig('LMX_C3_sankey_LMX1A.pdf', bbox_inches='tight')
plt.savefig('LMX_C3_sankey_LMX1A.png', bbox_inches='tight')
plt.show()


# In[9]:


import scanpy as sc
TOP_C4 = ['C4-248', 'C4-411']
TOP_C3 = ['C3-92',  'C3-84']

# Extract short IDs
adata_mouse.obs['C4_short'] = adata_mouse.obs['C4_HM'].str.extract(r'(C4-\d+)')
adata_mouse.obs['C3_short'] = adata_mouse.obs['C3_HM'].str.extract(r'(C3-\d+)')

# Create filtered columns — NaN for everything outside the top labels
adata_mouse.obs['C4_top2'] = adata_mouse.obs['C4_short'].where(
    adata_mouse.obs['C4_short'].isin(TOP_C4))

adata_mouse.obs['C3_top2'] = adata_mouse.obs['C3_short'].where(
    adata_mouse.obs['C3_short'].isin(TOP_C3))

# Assign fixed colors
adata_mouse.uns['C4_top2_colors'] = ['#FF4A46', '#1CE6FF']   # C4-248, C4-411
adata_mouse.uns['C3_top2_colors'] = ['#FFFF00', '#006FA6']   # C3-92,  C3-84

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

sc.pl.embedding(adata_mouse, basis='tsne', color='C4_top2',
                legend_loc='on data', legend_fontsize='small',
                title='Top 2 C4 Hypomap labels\n(C4-248, C4-411)',
                na_color='#eeeeee', ax=axes[0], show=False)

sc.pl.embedding(adata_mouse, basis='tsne', color='C3_top2',
                legend_loc='on data', legend_fontsize='small',
                title='Top 2 C3 Hypomap labels\n(C3-92, C3-84)',
                na_color='#eeeeee', ax=axes[1], show=False)

plt.tight_layout()
plt.savefig('allgenes_top2_C4_C3_labels.pdf', bbox_inches='tight')
plt.savefig('allgenes_top2_C4_C3_labels.png', bbox_inches='tight')
plt.show()


# In[32]:


# ── CELLECT BMI prioritization — corrected with unified CSV + master constants ─
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

UNIFIED_PATH = "hypomap_cellect_unified.csv"
GWAS         = 'BMI_Yengo2018'

# ── Master constants (already loaded — referenced here for clarity) ────────────
# C4
# GROUP_A_C4, GROUP_B_C4, LMX_C4_ALL, HL_C4
# C3
# GROUP_A_C3, CRNDE_C3, GABA_NEIGHBOR_C3, GROUP_B_C3, LMX_C3_ALL, HL_C3

LMX_C4_IDS = LMX_C4_ALL
LMX_C3_IDS = LMX_C3_ALL

BONF_C4    = 0.05 / 415
BONF_C3    = 0.05 / 156

# Colours — from master constants
COL_GRP_A_C4   = '#1a5276'   # C4 genuine Lmx1a
COL_GRP_A_C3   = C3_GROUP_COLORS['Group A']        # '#8B0000'
COL_CRNDE      = C3_GROUP_COLORS['CRNDE']           # '#e67e22'
COL_GABA_NB    = C3_GROUP_COLORS['GABA_neighbor']   # '#8e44ad'
COL_GRP_B      = '#95a5a6'
COL_HL         = '#8B0000'
COL_SIG        = '#E8A87C'
COL_NEG        = '#bdc3c7'

# ── 1. Load unified CSV ───────────────────────────────────────────────────────
df = pd.read_csv(UNIFIED_PATH)

def get(method, level):
    return (df[(df['gwas']   == GWAS) &
               (df['method'] == method) &
               (df['level']  == level) &
               (df['is_neuronal'])]
            .copy())

paper_c4 = get('paper_MAGMA', 'C4')
new_m_c4 = get('new_MAGMA',   'C4')
new_l_c4 = get('new_LDSC',    'C4')
new_m_c3 = get('new_MAGMA',   'C3')
new_l_c3 = get('new_LDSC',    'C3')

# C4 flags
for df_sub in [paper_c4, new_m_c4, new_l_c4]:
    df_sub['is_lmx']   = df_sub['cluster_id'].isin(LMX_C4_IDS)
    df_sub['is_grp_a'] = df_sub['cluster_id'].isin(GROUP_A_C4)
    df_sub['is_grp_b'] = df_sub['cluster_id'].isin(GROUP_B_C4)
    df_sub['is_hl']    = df_sub['cluster_id'].isin(HL_C4)

# C3 flags — full 4-way classification
for df_sub in [new_m_c3, new_l_c3]:
    df_sub['is_lmx']     = df_sub['cluster_id'].isin(LMX_C3_IDS)
    df_sub['is_grp_a']   = df_sub['cluster_id'].isin(GROUP_A_C3)
    df_sub['is_crnde']   = df_sub['cluster_id'].isin(CRNDE_C3)
    df_sub['is_gaba_nb'] = df_sub['cluster_id'].isin(GABA_NEIGHBOR_C3)
    df_sub['is_grp_b']   = df_sub['cluster_id'].isin(GROUP_B_C3)
    df_sub['is_hl']      = df_sub['cluster_id'].isin(HL_C3)

# ── 2. Sanity check printout ──────────────────────────────────────────────────
print("=" * 70)
print("SANITY CHECK: three-method CELLECT C4 comparison")
print("=" * 70)

for label, sub, bonf in [
    ('paper_MAGMA (bootstrapped ESmu)', paper_c4, BONF_C4),
    ('new_MAGMA   (single-shot ESmu)',  new_m_c4, BONF_C4),
    ('new_LDSC    (single-shot ESmu)',  new_l_c4, BONF_C4),
]:
    sig = (sub['pvalue'] < bonf).sum()
    top = sub.nsmallest(1,'pvalue').iloc[0]
    print(f"\n  {label}")
    print(f"  Bonf-sig: {sig:>3} / {len(sub)}  |  "
          f"top: {top['cluster_id']}  p={top['pvalue']:.2e}")
    print(f"  Top annotation: {top['annotation_full']}")

print(f"""
KEY CORRECTIONS:
  ✓ C4 Group A ({len(GROUP_A_C4)} clusters): ≥59% mouse-specific + human LMX1A >0.5
    Strongest: C4-229 p=3.10e-08 *** | 97% spec | 90% LMX1A+
  ✓ C4 Group B ({len(GROUP_B_C4)} clusters): <28% mouse-specific + near-zero LMX1A
    C4-7 excluded entirely (1.8% specific — noise label)
""")

print("=" * 70)
print("C3-LEVEL CELLECT BMI RESULTS (new_MAGMA — not in paper)")
print("=" * 70)
print(f"\n  C3 Group A ({len(GROUP_A_C3)}): genuine LMX1A+ Post-1 GLU  "
      f"[{', '.join(GROUP_A_C3)}]")
print(f"  C3 CRNDE  ({len(CRNDE_C3)}):    {CRNDE_C3}  — cross-species marker")
print(f"  C3 GABA-nb({len(GABA_NEIGHBOR_C3)}):    {GABA_NEIGHBOR_C3}"
      f"  — specific but LMX1A-negative")
print(f"  C3 Group B ({len(GROUP_B_C3)}): contaminants  "
      f"[incl. C3-115 Mid-2 GLU, C3-28 corrected]")
print(f"\n  Bonf-sig new_MAGMA: {(new_m_c3['pvalue']<BONF_C3).sum()}")
print(f"  Bonf-sig new_LDSC:  {(new_l_c3['pvalue']<BONF_C3).sum()}")

print(f"\nTop 10 C3 (new_MAGMA):")
print(new_m_c3.nsmallest(10,'pvalue')[
    ['cluster_id','annotation_full','pvalue','log10p','bonf_sig']
    ].to_string(index=False))

lmx_c3 = new_m_c3[new_m_c3['is_lmx']].sort_values('pvalue')
print(f"\nLmx1a C3 clusters (new_MAGMA):")
print(lmx_c3[['cluster_id','annotation_full','pvalue','log10p',
               'bonf_sig']].to_string(index=False))

# ── 3. Colour helpers ─────────────────────────────────────────────────────────
def c4_bar_color(row):
    if row['cluster_id'] in HL_C4:   return COL_HL
    if row['is_grp_b']:               return COL_GRP_B
    if row['bonf_sig']:               return COL_SIG
    if row['beta'] < 0:               return COL_NEG
    return COL_GRP_A_C4

def c3_bar_color(row):
    """4-way: Group A / CRNDE / GABA-neighbour / Group B."""
    cid = row['cluster_id']
    if cid in HL_C3:              return COL_HL
    if cid in GROUP_A_C3:
        return COL_SIG if row['bonf_sig'] else COL_GRP_A_C3
    if cid in CRNDE_C3:           return COL_CRNDE
    if cid in GABA_NEIGHBOR_C3:   return COL_GABA_NB
    if row['beta'] < 0:           return COL_NEG
    return COL_GRP_B

# ── 4. Figures ────────────────────────────────────────────────────────────────
sns.set_style("whitegrid")
plt.rcParams.update({'font.size':9,'axes.titlesize':11,'figure.dpi':150,
                     'savefig.dpi':300,'pdf.fonttype':42})

fig = plt.figure(figsize=(20, 15))
gs  = fig.add_gridspec(2, 3, hspace=0.48, wspace=0.42)

# ── Panel A: C4 histogram ─────────────────────────────────────────────────────
ax = fig.add_subplot(gs[0, 0])
ax.hist(paper_c4['log10p'], bins=50, color='#D0D0D0',
        edgecolor='white', linewidth=0.3, label=f'All C4 (n={len(paper_c4)})')
ax.hist(paper_c4[paper_c4['is_grp_a']]['log10p'], bins=15,
        color=COL_GRP_A_C4, alpha=0.85, edgecolor='white', linewidth=0.3,
        label=f'Group A genuine ({len(GROUP_A_C4)})')
ax.hist(paper_c4[paper_c4['is_grp_b']]['log10p'], bins=8,
        color=COL_GRP_B, alpha=0.85, edgecolor='white', linewidth=0.3,
        label=f'Group B contaminant ({len(GROUP_B_C4)})')
ax.axvline(-np.log10(BONF_C4), color='black', ls='--', lw=1,
           label=f'Bonf (p={BONF_C4:.1e})')
for _, r in paper_c4[paper_c4['cluster_id'].isin(HL_C4)].iterrows():
    ax.axvline(r['log10p'], color=COL_HL, lw=2, zorder=5)
    ax.annotate(r['cluster_id'],
                xy=(r['log10p'], ax.get_ylim()[1]*0.78),
                fontsize=7, fontweight='bold', color=COL_HL, ha='center',
                bbox=dict(boxstyle='round,pad=0.2', fc='white',
                          ec=COL_HL, alpha=0.9))
ax.set_xlabel('−log₁₀(p)', fontsize=9)
ax.set_ylabel('Number of C4 clusters')
ax.set_title('BMI enrichment — all C4\n(paper_MAGMA, bootstrapped ESmu)',
             fontweight='bold')
ax.legend(fontsize=7.5)
ax.text(0.98, 0.98,
        'Bonf-sig:\n  paper_MAGMA: 291\n  new_MAGMA:   240\n  new_LDSC:     10',
        transform=ax.transAxes, ha='right', va='top', fontsize=7,
        color='grey', style='italic')

# ── Panel B: C3 histogram — 4-way coloured ────────────────────────────────────
ax = fig.add_subplot(gs[0, 1])
ax.hist(new_m_c3['log10p'], bins=40, color='#D0D0D0',
        edgecolor='white', linewidth=0.3,
        label=f'All C3 (n={len(new_m_c3)})')
ax.hist(new_m_c3[new_m_c3['is_grp_a']]['log10p'], bins=10,
        color=COL_GRP_A_C3, alpha=0.85, edgecolor='white', linewidth=0.3,
        label=f'Group A genuine ({len(GROUP_A_C3)})')
ax.hist(new_m_c3[new_m_c3['is_crnde']]['log10p'], bins=5,
        color=COL_CRNDE, alpha=0.85, edgecolor='white', linewidth=0.3,
        label='CRNDE cross-species (1)')
ax.hist(new_m_c3[new_m_c3['is_gaba_nb']]['log10p'], bins=3,
        color=COL_GABA_NB, alpha=0.85, edgecolor='white', linewidth=0.3,
        label='GABA neighbour (1)')
for _, r in new_m_c3[new_m_c3['is_hl']].iterrows():
    ax.axvline(r['log10p'], color=COL_HL, lw=2.5, zorder=5)
    ax.annotate(r['cluster_id'],
                xy=(r['log10p'], ax.get_ylim()[1]*0.85),
                fontsize=8, fontweight='bold', color=COL_HL, ha='center',
                bbox=dict(boxstyle='round,pad=0.2', fc='white',
                          ec=COL_HL, alpha=0.9))
ax.axvline(-np.log10(BONF_C3), color='black', ls='--', lw=1,
           label=f'Bonf (p={BONF_C3:.1e})')
ax.axvline(-np.log10(0.05), color='grey', ls=':', lw=0.8)
ax.set_xlabel('−log₁₀(p)', fontsize=9)
ax.set_ylabel('Number of C3 clusters')
ax.set_title('BMI enrichment — all C3\n(new_MAGMA; paper has no CELLEX C3)',
             fontweight='bold')
ax.legend(fontsize=7.5)
ax.text(0.98, 0.98,
        'Bonf-sig:\n  new_MAGMA: 80\n  new_LDSC:   5',
        transform=ax.transAxes, ha='right', va='top', fontsize=7,
        color='grey', style='italic')

# ── Panel C: Lmx1a C3 bar chart — MAGMA + LDSC with 4-way colour ─────────────
ax = fig.add_subplot(gs[0, 2])
lmx_c3_sorted = lmx_c3.sort_values('log10p', ascending=True).copy()
lmx_c3_ldsc   = new_l_c3[new_l_c3['is_lmx']].set_index('cluster_id')

bar_y     = np.arange(len(lmx_c3_sorted))
bar_cols  = [c3_bar_color(r) for _, r in lmx_c3_sorted.iterrows()]
bar_edges = ['black' if r['cluster_id'] in HL_C3 else 'white'
              for _, r in lmx_c3_sorted.iterrows()]

ax.barh(bar_y - 0.2, lmx_c3_sorted['log10p'],
        color=bar_cols, edgecolor=bar_edges, linewidth=0.9,
        height=0.38, label='new_MAGMA')

ldsc_vals = [lmx_c3_ldsc.loc[r['cluster_id'], 'log10p']
             if r['cluster_id'] in lmx_c3_ldsc.index else 0
             for _, r in lmx_c3_sorted.iterrows()]
ax.barh(bar_y + 0.2, ldsc_vals,
        color=[c+'88' for c in bar_cols],
        edgecolor='white', linewidth=0.4,
        height=0.38, hatch='//', label='new_LDSC')

ax.axvline(-np.log10(BONF_C3), color='black', ls='--', lw=1,
           label=f'Bonf (p={BONF_C3:.1e})')
ax.axvline(-np.log10(0.05), color='grey', ls=':', lw=0.8)
ax.set_yticks(bar_y)
ax.set_yticklabels(lmx_c3_sorted['cluster_id'], fontsize=8)
for lbl, r in zip(ax.get_yticklabels(), lmx_c3_sorted.itertuples()):
    if r.cluster_id in HL_C3:
        lbl.set_fontweight('bold'); lbl.set_color(COL_HL)
    elif r.cluster_id in CRNDE_C3:
        lbl.set_color(COL_CRNDE)
    elif r.cluster_id in GABA_NEIGHBOR_C3:
        lbl.set_color(COL_GABA_NB)
    elif r.cluster_id in GROUP_B_C3:
        lbl.set_color('#7f8c8d')

ax.set_xlabel('−log₁₀(p)')
ax.set_title('Lmx1a C3 clusters\nBMI enrichment (MAGMA + LDSC)', fontweight='bold')
ax.legend(fontsize=7.5, loc='lower right')

legend_c3 = [
    mpatches.Patch(fc=COL_HL,       ec='black', lw=1,
                   label=f'Highlighted ({", ".join(HL_C3)})'),
    mpatches.Patch(fc=COL_SIG,
                   label='Group A Bonf-sig'),
    mpatches.Patch(fc=COL_GRP_A_C3,
                   label=f'Group A nominal ({len(GROUP_A_C3)})'),
    mpatches.Patch(fc=COL_CRNDE,
                   label='CRNDE cross-species (C3-95)'),
    mpatches.Patch(fc=COL_GABA_NB,
                   label='GABA neighbour (C3-28)'),
    mpatches.Patch(fc=COL_GRP_B,
                   label=f'Group B contaminant ({len(GROUP_B_C3)})'),
]
ax.legend(handles=legend_c3, fontsize=6.5, loc='lower right')

# ── Panel D: Lmx1a C4 bar chart ───────────────────────────────────────────────
ax = fig.add_subplot(gs[1, 0])
lmx_c4_paper = paper_c4[paper_c4['is_lmx']].sort_values('log10p',
                                                          ascending=True).copy()
bar_cols_c4  = [c4_bar_color(r) for _, r in lmx_c4_paper.iterrows()]
bar_edges_c4 = ['black' if r['cluster_id'] in HL_C4 else 'white'
                for _, r in lmx_c4_paper.iterrows()]

ax.barh(range(len(lmx_c4_paper)), lmx_c4_paper['log10p'],
        color=bar_cols_c4, edgecolor=bar_edges_c4, linewidth=0.9, height=0.72)
ax.axvline(-np.log10(BONF_C4), color='black', ls='--', lw=1,
           label=f'Bonf (p={BONF_C4:.1e})')
ax.axvline(-np.log10(0.05), color='grey', ls=':', lw=0.8)
ax.set_yticks(range(len(lmx_c4_paper)))
ax.set_yticklabels(lmx_c4_paper['cluster_id'], fontsize=8)
for lbl, r in zip(ax.get_yticklabels(), lmx_c4_paper.itertuples()):
    if r.cluster_id in HL_C4:
        lbl.set_fontweight('bold'); lbl.set_color(COL_HL)
    elif r.cluster_id in GROUP_B_C4:
        lbl.set_color('#7f8c8d')

ax.set_xlabel('−log₁₀(p)')
ax.set_title('Lmx1a C4 clusters\nBMI enrichment (paper_MAGMA)', fontweight='bold')
ax.legend(handles=[
    mpatches.Patch(fc=COL_HL,      ec='black', lw=1,
                   label=f'Highlighted ({", ".join(HL_C4)})'),
    mpatches.Patch(fc=COL_SIG,                  label='Bonf-sig Group A'),
    mpatches.Patch(fc=COL_GRP_A_C4,             label='Nominal Group A'),
    mpatches.Patch(fc=COL_GRP_B,                label='Group B contaminant'),
    mpatches.Patch(fc=COL_NEG,                   label='Negative beta'),
], fontsize=6.5, loc='lower right')

# ── Panel E: Three-method C4 comparison ──────────────────────────────────────
ax = fig.add_subplot(gs[1, 1:])
lmx_m = new_m_c4[new_m_c4['is_lmx']].set_index('cluster_id')
lmx_l = new_l_c4[new_l_c4['is_lmx']].set_index('cluster_id')

sorted_ids = (lmx_c4_paper.sort_values('log10p', ascending=False)
              ['cluster_id'].tolist())
n = len(sorted_ids)
x_p, x_nm, x_nl = np.arange(n)*3, np.arange(n)*3+0.85, np.arange(n)*3+1.70
w = 0.78

def get_log10p(idx, cid):
    return idx.loc[cid,'log10p'] if cid in idx.index else 0.0

for i, cid in enumerate(sorted_ids):
    base = COL_GRP_B if cid in GROUP_B_C4 else COL_GRP_A_C4
    edge = 'black' if cid in HL_C4 else 'white'
    lw   = 1.2 if cid in HL_C4 else 0.3
    lp_p  = lmx_c4_paper[lmx_c4_paper['cluster_id']==cid]['log10p'].values[0]
    ax.bar(x_p[i],  lp_p,               width=w, color=base, alpha=1.0,
           edgecolor=edge, linewidth=lw)
    ax.bar(x_nm[i], get_log10p(lmx_m,cid), width=w, color=base, alpha=0.65,
           edgecolor=edge, linewidth=lw, hatch='//')
    ax.bar(x_nl[i], get_log10p(lmx_l,cid), width=w, color=base, alpha=0.35,
           edgecolor=edge, linewidth=lw, hatch='xx')

ax.axhline(-np.log10(BONF_C4), color='black', ls='--', lw=1)
ax.axhline(-np.log10(0.05),    color='grey',  ls=':', lw=0.8)
ax.set_xticks(x_p + 1.27)
ax.set_xticklabels(sorted_ids, rotation=45, ha='right', fontsize=7.5)
for lbl, cid in zip(ax.get_xticklabels(), sorted_ids):
    if cid in HL_C4:       lbl.set_fontweight('bold'); lbl.set_color(COL_HL)
    elif cid in GROUP_B_C4: lbl.set_color('#7f8c8d')
ax.set_ylabel('−log₁₀(p)', fontsize=9)
ax.set_title('Lmx1a C4 — three-method BMI comparison (Yengo2018)\n'
             'Group A (blue) = genuine LMX1A+   |   Group B (grey) = contaminant',
             fontweight='bold')
ax.legend(handles=[
    mpatches.Patch(fc=COL_GRP_A_C4, alpha=1.0,
                   label='paper_MAGMA bootstrapped (291 sig)'),
    mpatches.Patch(fc=COL_GRP_A_C4, alpha=0.65, hatch='//',
                   label='new_MAGMA single-shot (240 sig)'),
    mpatches.Patch(fc=COL_GRP_A_C4, alpha=0.35, hatch='xx',
                   label='new_LDSC single-shot (10 sig)'),
    mpatches.Patch(fc=COL_GRP_B,    label='Group B contaminant'),
], fontsize=7.5, loc='upper right')

plt.suptitle(
    'CELLECT BMI prioritization — Lmx1a Hypomap labels (CORRECTED)\n'
    f'C4 (n=415 neuronal) and C3 (n=156)  |  GWAS: {GWAS}\n'
    'Unified CSV source: hypomap_cellect_unified.csv',
    fontweight='bold', fontsize=12, y=1.01
)
plt.savefig('LMX_CELLECT_C3_C4_BMI_corrected.pdf', bbox_inches='tight')
plt.savefig('LMX_CELLECT_C3_C4_BMI_corrected.png', bbox_inches='tight')
plt.show()

# ── 5. Final summary ──────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("KEY FINDINGS — CORRECTED")
print("=" * 70)
top_paper = paper_c4.nsmallest(1,'pvalue').iloc[0]
top_newm  = new_m_c4.nsmallest(1,'pvalue').iloc[0]
top_newl  = new_l_c4.nsmallest(1,'pvalue').iloc[0]

print(f"""
Top C4 hits:
  paper_MAGMA: {top_paper['cluster_id']}  p={top_paper['pvalue']:.2e}
  new_MAGMA:   {top_newm['cluster_id']}   p={top_newm['pvalue']:.2e}
  new_LDSC:    {top_newl['cluster_id']}   p={top_newl['pvalue']:.2e}

C4 Group A ({len(GROUP_A_C4)} clusters — genuine LMX1A+):
  Bonf-sig paper_MAGMA: {paper_c4[paper_c4['cluster_id'].isin(GROUP_A_C4)]['bonf_sig'].sum()} / {len(GROUP_A_C4)}
  Strongest: C4-229 p=3.10e-08 *** | 97% spec | 90% LMX1A+

C4 Group B ({len(GROUP_B_C4)} clusters — contaminant):
  All nominal in paper_MAGMA  |  C4-7 excluded (1.8% specific)

C3 Group A ({len(GROUP_A_C3)} clusters — genuine LMX1A+ Post-1 GLU):
  {', '.join(GROUP_A_C3)}
  Bonf-sig new_MAGMA: {new_m_c3[new_m_c3['cluster_id'].isin(GROUP_A_C3)]['bonf_sig'].sum()} / {len(GROUP_A_C3)}

C3 CRNDE (C3-95):  MAGMA p={new_m_c3[new_m_c3['cluster_id']=='C3-95']['pvalue'].values[0]:.3e} ***
  34.6% mouse-specific BUT human LMX1A=1.23 — cross-species marker

C3 GABA neighbour (C3-28):  MAGMA p={new_m_c3[new_m_c3['cluster_id']=='C3-28']['pvalue'].values[0]:.3e} ***
  58% mouse-specific BUT human LMX1A=0.001 — Lmx1a-adjacent Pre-1 GABA
  NOT LMX1A+ in humans — separate category in all plots

C3 Group B ({len(GROUP_B_C3)} clusters — contaminants):
  C3-115: 5.5% specific — Mid-2 GLU contaminator (BMI signal unrelated to Lmx1a)
  C3-24:  1.6% specific — noise label
""")


# In[43]:


from statsmodels.stats.multitest import multipletests

# ── Compute FDR threshold from all C4 LDSC p-values ──────────────────────────
all_pvals   = ldsc_idx['pvalue'].values
_, fdr_pvals, _, _ = multipletests(all_pvals, method='fdr_bh')
ldsc_idx    = ldsc_idx.copy()
ldsc_idx['fdr'] = fdr_pvals

# FDR threshold as -log10(p) of the largest raw p that still passes FDR<0.05
fdr_sig_mask = fdr_pvals < 0.05
if fdr_sig_mask.any():
    fdr_cutoff_p    = all_pvals[fdr_sig_mask].max()
    fdr_cutoff_log  = -np.log10(fdr_cutoff_p)
else:
    fdr_cutoff_log  = None
    print("No clusters pass FDR<0.05")
# ── Panel D replot: Lmx1a C4 — new_LDSC, Sankey order, simplified colours ────
fig, ax = plt.subplots(figsize=(5, 6))

ldsc_idx = new_l_c4.set_index('cluster_id')

# Reverse c4_order so most abundant is at top
c4_order_flipped = list(reversed(c4_order))

bar_y_vals  = []
bar_colors  = []
bar_edges   = []
bar_lws     = []
tick_labels = []
HL_C4 = []
for cid in c4_order_flipped:
    lp   = ldsc_idx.loc[cid, 'log10p'] if cid in ldsc_idx.index else 0.0
    col  = (C4_GROUP_COLORS['Group A'] if cid in set(GROUP_A_C4)
            else C4_GROUP_COLORS['Group B'])
    edge = 'black' if cid in HL_C4 else 'white'
    lw   = 1.0     if cid in HL_C4 else 0.3

    bar_y_vals.append(lp)
    bar_colors.append(col)
    bar_edges.append(edge)
    bar_lws.append(lw)
    tick_labels.append(cid)

y_pos = np.arange(len(c4_order_flipped))

ax.barh(y_pos, bar_y_vals,
        color=bar_colors, edgecolor=bar_edges,
        linewidth=bar_lws, height=0.72)
ax.axvline(-np.log10(0.05), color='grey', ls=':', lw=0.8)
if fdr_cutoff_log is not None:
    ax.axvline(fdr_cutoff_log, color='steelblue', ls='--', lw=1)
ax.legend(handles=[
    mpatches.Patch(fc=C4_GROUP_COLORS['Group A'], label='Group A (LMX1A+)'),
    mpatches.Patch(fc=C4_GROUP_COLORS['Group B'], label='Group B (LMX1A−)'),
    plt.Line2D([0],[0], color='steelblue', ls='--', lw=1,
               label=f'FDR<0.05 (p={fdr_cutoff_p:.2f})' if fdr_cutoff_log else 'FDR<0.05 (none)'),
    plt.Line2D([0],[0], color='grey', ls=':', lw=0.8, label='Nominal p=0.05'),
], fontsize=8, loc='lower right')
ax.set_yticks(y_pos)
ax.set_yticklabels(tick_labels, fontsize=8)
for lbl, cid in zip(ax.get_yticklabels(), c4_order_flipped):
    if cid in HL_C4:
        lbl.set_fontweight('bold'); lbl.set_color(C4_GROUP_COLORS['Group A'])
    elif cid in set(GROUP_B_C4):
        lbl.set_color('#7f8c8d')
ax.set_xlabel('−log₁₀(p)  LDSC  (Yengo2018)', fontsize=9)
ax.set_title('Lmx1a C4 clusters\nBMI enrichment  (CELLECT-LDSC)',
             fontweight='bold')
plt.tight_layout()
plt.savefig('LMX_C4_LDSC_sankey_order.pdf', bbox_inches='tight')
plt.savefig('LMX_C4_LDSC_sankey_order.png', bbox_inches='tight')
plt.show()


# In[19]:


# ── C3 MAGMA vs LDSC comparison — Yengo2018, master constants ────────────────
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.stats import spearmanr, pearsonr

UNIFIED_PATH = "hypomap_cellect_unified.csv"
GWAS         = 'BMI_Yengo2018'
BONF_C4      = 0.05 / 415
BONF_C3      = 0.05 / 156

# Master constants already loaded:
# GROUP_A_C3, CRNDE_C3, GABA_NEIGHBOR_C3, GROUP_B_C3, LMX_C3_ALL
# HL_C3, C3_GROUP_COLORS
# GROUP_A_C4, GROUP_B_C4, LMX_C4_ALL, HL_C4

# Convenience colour refs
COL_HL_LMX = C3_GROUP_COLORS['Group A']        # '#8B0000'
COL_CRNDE  = C3_GROUP_COLORS['CRNDE']           # '#e67e22'
COL_GABA   = C3_GROUP_COLORS['GABA_neighbor']   # '#8e44ad'
COL_LMX    = '#1a5276'                           # other Group A
COL_GRP_B  = '#95a5a6'

# Ensure set types for membership tests with | operator
_HL_C3      = set(HL_C3)
_CRNDE_C3   = set(CRNDE_C3)
_GABA_C3    = set(GABA_NEIGHBOR_C3)
_GROUP_A_C3 = set(GROUP_A_C3)
_GROUP_B_C3 = set(GROUP_B_C3)

# ── 1. Load ───────────────────────────────────────────────────────────────────
unified = pd.read_csv(UNIFIED_PATH)

def get(method, level):
    return (unified[
        (unified['method']    == method) &
        (unified['gwas']      == GWAS) &
        (unified['level']     == level) &
        (unified['is_neuronal'])
    ].copy().set_index('cluster_id'))

magma_c3 = get('new_MAGMA', 'C3')
ldsc_c3  = get('new_LDSC',  'C3')
magma_c4 = get('new_MAGMA', 'C4')
ldsc_c4  = get('new_LDSC',  'C4')

# ── 2. Build comparison tables ────────────────────────────────────────────────
def build_comp_c3(m_df, l_df):
    comp = pd.DataFrame({
        'p_magma':       m_df['pvalue'],
        'beta_magma':    m_df['beta'],
        'log10p_magma':  m_df['log10p'],
        'p_ldsc':        l_df['pvalue'],
        'beta_ldsc':     l_df['beta'],
        'log10p_ldsc':   l_df['log10p'],
        'annotation':    m_df['annotation_full'],
    }).dropna(subset=['p_magma','p_ldsc'])

    comp['is_lmx']     = comp.index.isin(LMX_C3_ALL)
    comp['is_grp_a']   = comp.index.isin(GROUP_A_C3)
    comp['is_crnde']   = comp.index.isin(CRNDE_C3)
    comp['is_gaba_nb'] = comp.index.isin(GABA_NEIGHBOR_C3)
    comp['is_grp_b']   = comp.index.isin(GROUP_B_C3)
    comp['is_hl']      = comp.index.isin(HL_C3)

    comp['bonf_magma']   = comp['p_magma'] < BONF_C3
    comp['bonf_ldsc']    = comp['p_ldsc']  < BONF_C3
    comp['nom_magma']    = comp['p_magma'] < 0.05
    comp['nom_ldsc']     = comp['p_ldsc']  < 0.05
    comp['nom_both']     = comp['nom_magma'] & comp['nom_ldsc']
    comp['signed_magma'] = np.sign(comp['beta_magma']) * comp['log10p_magma']
    comp['signed_ldsc']  = np.sign(comp['beta_ldsc'])  * comp['log10p_ldsc']
    return comp

def build_comp_c4(m_df, l_df):
    comp = pd.DataFrame({
        'p_magma':       m_df['pvalue'],
        'beta_magma':    m_df['beta'],
        'log10p_magma':  m_df['log10p'],
        'p_ldsc':        l_df['pvalue'],
        'beta_ldsc':     l_df['beta'],
        'log10p_ldsc':   l_df['log10p'],
    }).dropna(subset=['p_magma','p_ldsc'])
    comp['is_lmx']   = comp.index.isin(LMX_C4_ALL)
    comp['is_grp_a'] = comp.index.isin(GROUP_A_C4)
    comp['is_grp_b'] = comp.index.isin(GROUP_B_C4)
    comp['is_hl']    = comp.index.isin(HL_C4)
    comp['signed_magma'] = np.sign(comp['beta_magma']) * comp['log10p_magma']
    comp['signed_ldsc']  = np.sign(comp['beta_ldsc'])  * comp['log10p_ldsc']
    return comp

c3_comp = build_comp_c3(magma_c3, ldsc_c3)
c4_comp = build_comp_c4(magma_c4, ldsc_c4)

# ── 3. Correlations ───────────────────────────────────────────────────────────
r_sp,  p_sp  = spearmanr(c3_comp['log10p_magma'], c3_comp['log10p_ldsc'])
r_sg,  p_sg  = spearmanr(c3_comp['signed_magma'],  c3_comp['signed_ldsc'])

print("=" * 72)
print(f"C3 new_MAGMA vs new_LDSC  |  GWAS: {GWAS}")
print("=" * 72)
print(f"\n  new_MAGMA C3 Bonf-sig: {magma_c3['bonf_sig'].sum()} / {len(magma_c3)}")
print(f"  new_LDSC  C3 Bonf-sig: {ldsc_c3['bonf_sig'].sum()}  / {len(ldsc_c3)}")
dir_ag = (np.sign(c3_comp['signed_magma'])==np.sign(c3_comp['signed_ldsc'])).mean()*100
print(f"  Direction agree:       {dir_ag:.0f}%")
print(f"  Spearman (unsigned):   ρ={r_sp:.3f}  p={p_sp:.2e}")
print(f"  Spearman (signed):     ρ={r_sg:.3f}  p={p_sg:.2e}")

print(f"\nLmx1a C3 clusters — new_MAGMA vs new_LDSC:")
print(f"{'C3_id':<10} {'group':<12} {'p_MAGMA':>10} {'p_LDSC':>10} "
      f"{'Bonf_M':>7} {'Bonf_L':>7} {'nom_both':>9}")
print("-" * 68)
lmx_c3_comp = c3_comp[c3_comp['is_lmx']].sort_values('p_magma')
for idx, r in lmx_c3_comp.iterrows():
    grp = ('Group A'   if r['is_grp_a'] else
           'CRNDE'     if r['is_crnde'] else
           'GABA-nb'   if r['is_gaba_nb'] else 'Group B')
    bm = '***' if r['bonf_magma'] else ''
    bl = '***' if r['bonf_ldsc']  else ''
    nb = '✓'   if r['nom_both']   else ''
    print(f"{idx:<10} {grp:<12} {r['p_magma']:>10.3e} {r['p_ldsc']:>10.3e} "
          f"{bm:>7} {bl:>7} {nb:>9}")

# ── 4. Dot colour helper ──────────────────────────────────────────────────────
def dot_color(idx, comp_df):
    if idx in _HL_C3:    return COL_HL_LMX
    if idx in _CRNDE_C3: return COL_CRNDE
    if idx in _GABA_C3:  return COL_GABA
    if idx in _GROUP_A_C3: return COL_LMX
    return COL_GRP_B


# ── 5. Figures ────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(18, 7))

# ── Panel A: C3 scatter — MAGMA vs LDSC ──────────────────────────────────────
ax = axes[0]
ax.scatter(c3_comp['log10p_magma'], c3_comp['log10p_ldsc'],
           c='#dddddd', s=18, alpha=0.55, zorder=2)

for idx, r in c3_comp[c3_comp['is_lmx']].iterrows():
    col = dot_color(idx, c3_comp)
    sz = 110 if idx in _HL_C3 else 90 if idx in (_CRNDE_C3 | _GABA_C3) else 55
    ax.scatter(r['log10p_magma'], r['log10p_ldsc'],
               c=col, s=sz, zorder=5, edgecolors='white', linewidths=0.5)
    if idx in (_HL_C3 | _CRNDE_C3 | _GABA_C3):
        ax.annotate(idx, xy=(r['log10p_magma'], r['log10p_ldsc']),
                    xytext=(5,3), textcoords='offset points',
                    fontsize=7.5, fontweight='bold', color=col)

m_fit, b_fit = np.polyfit(c3_comp['log10p_magma'], c3_comp['log10p_ldsc'], 1)
xr = np.linspace(0, c3_comp['log10p_magma'].max(), 200)
ax.plot(xr, m_fit*xr+b_fit, color='#e74c3c', lw=1.5, ls='--',
        label=f'OLS slope={m_fit:.2f}')

bx = -np.log10(BONF_C3)
if bx <= c3_comp['log10p_magma'].max()*1.1:
    ax.axvline(bx, color='steelblue', ls='--', lw=0.9, alpha=0.8,
               label=f'Bonf MAGMA')
if bx <= c3_comp['log10p_ldsc'].max()*1.1:
    ax.axhline(bx, color='#c0392b', ls='--', lw=0.9, alpha=0.8,
               label=f'Bonf LDSC')
ax.axhline(-np.log10(0.05), color='grey', ls=':', lw=0.7, alpha=0.6)
ax.axvline(-np.log10(0.05), color='grey', ls=':', lw=0.7, alpha=0.6)
ax.set_xlim(0, c3_comp['log10p_magma'].max()*1.06)
ax.set_ylim(0, c3_comp['log10p_ldsc'].max()*1.06)
ax.set_xlabel(f'−log₁₀(p)  new_MAGMA  (Bonf-sig: {magma_c3["bonf_sig"].sum()})',
              fontsize=10)
ax.set_ylabel(f'−log₁₀(p)  new_LDSC   (Bonf-sig: {ldsc_c3["bonf_sig"].sum()})',
              fontsize=10)
ax.set_title(f'C3  new_MAGMA vs new_LDSC\n'
             f'Spearman ρ={r_sp:.3f}  p={p_sp:.2e}  (n={len(c3_comp)})',
             fontweight='bold')
ax.legend(handles=[
    mpatches.Patch(fc=COL_HL_LMX, label=f'Group A highlight ({", ".join(HL_C3)})'),
    mpatches.Patch(fc=COL_LMX,    label=f'Group A other ({len(GROUP_A_C3)})'),
    mpatches.Patch(fc=COL_CRNDE,  label='CRNDE C3-95'),
    mpatches.Patch(fc=COL_GABA,   label='GABA-nb C3-28'),
    mpatches.Patch(fc=COL_GRP_B,  label=f'Group B ({len(GROUP_B_C3)})'),
], fontsize=7.5, loc='upper left')

# ── Panel B: C3 ranked landscape ─────────────────────────────────────────────
ax = axes[1]
c3_ranked = (magma_c3.reset_index()
             .sort_values('log10p', ascending=False)
             .reset_index(drop=True))
c3_ranked['rank_pos'] = range(len(c3_ranked))

ax.scatter(c3_ranked['rank_pos'], c3_ranked['log10p'],
           c='#dddddd', s=15, alpha=0.55, zorder=2)

for _, row_m in c3_ranked.iterrows():
    cid = row_m['cluster_id']
    if cid in ldsc_c3.index:
        lp_scaled = ldsc_c3.loc[cid,'log10p'] * (
            c3_ranked['log10p'].max() / max(ldsc_c3['log10p'].max(), 1e-9))
        ax.scatter(row_m['rank_pos'], lp_scaled,
                   marker='^', c='#aaaaaa', s=15, alpha=0.4, zorder=2)

for _, r in c3_ranked[c3_ranked['cluster_id'].isin(LMX_C3_ALL)].iterrows():
    cid = r['cluster_id']
    col = dot_color(cid, c3_comp)
    sz = 100 if cid in (_HL_C3 | _CRNDE_C3 | _GABA_C3) else 50
    ax.scatter(r['rank_pos'], r['log10p'],
               c=col, s=sz, zorder=5, edgecolors='white', linewidths=0.5)
    if cid in (_HL_C3 | _CRNDE_C3 | _GABA_C3):
        ax.annotate(cid, xy=(r['rank_pos'], r['log10p']),
                    xytext=(3,3), textcoords='offset points',
                    fontsize=7.5, fontweight='bold', color=col)

ax.axhline(-np.log10(BONF_C3), color='black', ls='--', lw=1,
           label=f'MAGMA Bonf (p={BONF_C3:.1e})')
ax.axhline(-np.log10(0.05), color='grey', ls=':', lw=0.8)
ax.set_xlabel('C3 cluster rank (1=most significant, MAGMA)', fontsize=10)
ax.set_ylabel('−log₁₀(p)  new_MAGMA\n▲ = new_LDSC scaled', fontsize=9)
ax.set_title(f'C3 BMI landscape\n'
             f'{magma_c3["bonf_sig"].sum()} MAGMA  |  '
             f'{ldsc_c3["bonf_sig"].sum()} LDSC Bonf-sig',
             fontweight='bold')
ax.legend(fontsize=7.5, handles=[
    mpatches.Patch(fc=COL_HL_LMX, label='Group A highlight'),
    mpatches.Patch(fc=COL_LMX,    label='Group A other'),
    mpatches.Patch(fc=COL_CRNDE,  label='CRNDE C3-95'),
    mpatches.Patch(fc=COL_GABA,   label='GABA-nb C3-28'),
    mpatches.Patch(fc='#aaaaaa',  label='LDSC scaled (▲)'),
])

# ── Panel C: Lmx1a C3 bar chart — dual y-axis MAGMA + LDSC ──────────────────
ax1   = axes[2]
# Use LMX_C3_ALL in defined group order for interpretable layout
focus_ids = (list(GROUP_A_C3) + list(CRNDE_C3) +
             list(GABA_NEIGHBOR_C3) + list(GROUP_B_C3))
focus_labels = ['C3-95\nCRNDE ★' if c=='C3-95' else c for c in focus_ids]

def safe(idx_df, cid):
    return idx_df.loc[cid,'log10p'] if cid in idx_df.index else 0.0

lp_m = [safe(magma_c3, c) for c in focus_ids]
lp_l = [safe(ldsc_c3,  c) for c in focus_ids]

def bar_col_m(cid):
    if cid in _HL_C3:    return COL_HL_LMX
    if cid in _GROUP_A_C3: return COL_LMX
    if cid in _CRNDE_C3: return COL_CRNDE
    if cid in _GABA_C3:  return COL_GABA
    return COL_GRP_B

def bar_col_l(cid):
    """Darker version of MAGMA colour for LDSC bars."""
    mapping = {COL_HL_LMX:'#5c0000', COL_LMX:'#0d2e4d',
                COL_CRNDE:'#c0742a', COL_GABA:'#5b2c6f',
                COL_GRP_B:'#717d7e'}
    return mapping.get(bar_col_m(cid), '#333333')

x = np.arange(len(focus_ids)); w = 0.38
cols_m = [bar_col_m(c) for c in focus_ids]
cols_l = [bar_col_l(c) for c in focus_ids]

# Group separator positions
sep_positions = [len(GROUP_A_C3)-0.5,
                 len(GROUP_A_C3)+len(CRNDE_C3)-0.5,
                 len(GROUP_A_C3)+len(CRNDE_C3)+len(GABA_NEIGHBOR_C3)-0.5]

ax1.bar(x - w/2, lp_m, width=w, color=cols_m,
        edgecolor='white', linewidth=0.4, label='new_MAGMA (left axis)')
ax1.set_ylabel('−log₁₀(p)  new_MAGMA', fontsize=9, color='#1a5276')
ax1.tick_params(axis='y', labelcolor='#1a5276')
ax1.axhline(-np.log10(BONF_C3), color='steelblue', ls='--', lw=1,
            label=f'MAGMA Bonf')
ax1.axhline(-np.log10(0.05), color='steelblue', ls=':', lw=0.7, alpha=0.6)

ax2 = ax1.twinx()
ax2.bar(x + w/2, lp_l, width=w, color=cols_l,
        edgecolor='white', linewidth=0.4, alpha=0.75,
        hatch='//', label='new_LDSC (right axis)')
ax2.set_ylabel('−log₁₀(p)  new_LDSC', fontsize=9, color='#c0392b')
ax2.tick_params(axis='y', labelcolor='#c0392b')
ax2.axhline(-np.log10(BONF_C3), color='#c0392b', ls='--', lw=1)

# Group separator lines + labels
for sp in sep_positions:
    ax1.axvline(sp, color='black', lw=0.8, ls='-', alpha=0.25)
for pos, lbl in [(len(GROUP_A_C3)/2-0.5,       'Group A'),
                  (len(GROUP_A_C3)+0.5,          'CRNDE'),
                  (len(GROUP_A_C3)+1.5,          'GABA-nb'),
                  (len(GROUP_A_C3)+len(CRNDE_C3)+len(GABA_NEIGHBOR_C3)
                   + len(GROUP_B_C3)/2 - 0.5,    'Group B')]:
    ax1.text(pos, ax1.get_ylim()[1]*0.92, lbl, ha='center',
             fontsize=7, color='grey', style='italic')

ax1.set_xticks(x)
ax1.set_xticklabels(focus_labels, rotation=45, ha='right', fontsize=7.5)
for lbl, cid in zip(ax1.get_xticklabels(), focus_ids):
    if cid in (_HL_C3 | _CRNDE_C3 | _GABA_C3):
        lbl.set_fontweight('bold')
    if cid in GROUP_B_C3:
        lbl.set_color('#7f8c8d')

ax1.set_title(f'Lmx1a C3 BMI enrichment by group\n'
              f'new_MAGMA (left) vs new_LDSC (right)  |  {GWAS}',
              fontweight='bold')

lines1, labs1 = ax1.get_legend_handles_labels()
lines2, labs2 = ax2.get_legend_handles_labels()
ax1.legend(lines1+lines2, labs1+labs2, fontsize=7, loc='upper right')

plt.suptitle(
    f'C3-level CELLECT BMI  —  new_MAGMA vs new_LDSC  ({GWAS})\n'
    f'4-way C3 classification: Group A | CRNDE | GABA-neighbour | Group B',
    fontweight='bold', fontsize=12, y=1.01
)
plt.tight_layout()
plt.savefig('LMX_C3_MAGMA_vs_LDSC.pdf', bbox_inches='tight')
plt.savefig('LMX_C3_MAGMA_vs_LDSC.png', bbox_inches='tight')
plt.show()

# ── 6. Summary ────────────────────────────────────────────────────────────────
print("\n" + "=" * 72)
print(f"SUMMARY: new_MAGMA vs new_LDSC  |  C3 Lmx1a  |  {GWAS}")
print("=" * 72)
print(f"""
Method power (Bonf threshold {BONF_C3:.2e}):
  new_MAGMA:  {magma_c3['bonf_sig'].sum()} / {len(magma_c3)} Bonf-sig
  new_LDSC:   {ldsc_c3['bonf_sig'].sum()}  / {len(ldsc_c3)} Bonf-sig
  Spearman ρ (unsigned): {r_sp:.3f}  p={p_sp:.2e}
  Direction agree: {dir_ag:.0f}%

Group A genuine LMX1A+ ({len(GROUP_A_C3)} clusters):""")
for cid in sorted(GROUP_A_C3,
                  key=lambda c: c3_comp.loc[c,'p_magma']
                  if c in c3_comp.index else 1):
    if cid not in c3_comp.index: continue
    r = c3_comp.loc[cid]
    bm = '***' if r['bonf_magma'] else ('nom' if r['nom_magma'] else ' ')
    bl = '***' if r['bonf_ldsc']  else ('nom' if r['nom_ldsc']  else ' ')
    print(f"  {cid:<10}  MAGMA={r['p_magma']:.3e} {bm}  "
          f"LDSC={r['p_ldsc']:.3e} {bl}")

print(f"\nCRNDE (C3-95) — cross-species marker:")
r95 = c3_comp.loc['C3-95']
print(f"  MAGMA={r95['p_magma']:.3e} {'***' if r95['bonf_magma'] else 'nom'}  "
      f"LDSC={r95['p_ldsc']:.3e} {'***' if r95['bonf_ldsc'] else 'nom'}")

print(f"\nGABA-neighbour (C3-28) — specific but LMX1A-negative in human:")
r28 = c3_comp.loc['C3-28']
print(f"  MAGMA={r28['p_magma']:.3e} {'***' if r28['bonf_magma'] else 'nom'}  "
      f"LDSC={r28['p_ldsc']:.3e} {'nom' if r28['nom_ldsc'] else ' '}")
print(f"  → BMI signal is real but NOT via LMX1A mechanism")

print(f"\nGroup B contaminants ({len(GROUP_B_C3)}):")
print(f"  {', '.join(GROUP_B_C3)}")
print(f"  C3-115: 5.5% mouse-specific, Mid-2 GLU contaminator")
print(f"  C3-28 corrected to GABA-neighbour category")


# In[21]:


from statsmodels.stats.multitest import multipletests
import pandas as pd
import numpy as np

# ── Load LDSC C3 results from unified CSV ─────────────────────────────────────
unified  = pd.read_csv("hypomap_cellect_unified.csv")
GWAS     = 'BMI_Yengo2018'
BONF_C3  = 0.05 / 156

new_l_c3 = (unified[(unified['gwas']      == GWAS) &
                     (unified['method']    == 'new_LDSC') &
                     (unified['level']     == 'C3') &
                     (unified['is_neuronal'])]
            .copy())

print(f"new_l_c3: {len(new_l_c3)} rows  "
      f"({(new_l_c3['pvalue'] < BONF_C3).sum()} Bonf-sig)")

# ── Compute FDR threshold from all C3 LDSC p-values ──────────────────────────
ldsc_idx_c3  = new_l_c3.set_index('cluster_id')
all_pvals_c3 = ldsc_idx_c3['pvalue'].values
_, fdr_pvals_c3, _, _ = multipletests(all_pvals_c3, method='fdr_bh')
ldsc_idx_c3  = ldsc_idx_c3.copy()
ldsc_idx_c3['fdr'] = fdr_pvals_c3

fdr_sig_mask_c3 = fdr_pvals_c3 < 0.05
if fdr_sig_mask_c3.any():
    fdr_cutoff_p_c3   = all_pvals_c3[fdr_sig_mask_c3].max()
    fdr_cutoff_log_c3 = -np.log10(fdr_cutoff_p_c3)
else:
    fdr_cutoff_log_c3 = None
    print("No C3 clusters pass FDR<0.05")

# ── C3 bar chart: Sankey order, simplified colours ────────────────────────────
fig, ax = plt.subplots(figsize=(5, 6))

# Sankey order (Group A desc by mouse count, then CRNDE, GABA-nb, Group B)
# reversed so most abundant is at top
c3_order_flipped = list(reversed(c3_order))

bar_y_vals  = []
bar_colors  = []
bar_edges   = []
bar_lws     = []
tick_labels = []
HL_C3_PLOT  = []   # no highlights for this simplified version

def c3_bar_col_simple(cid):
    if cid in set(GRP_A_C3_PLOT):  return C3_GROUP_PLOT_COLORS['Group A']
    return C3_GROUP_PLOT_COLORS['Group B']

for cid in c3_order_flipped:
    lp   = ldsc_idx_c3.loc[cid, 'log10p'] if cid in ldsc_idx_c3.index else 0.0
    col  = c3_bar_col_simple(cid)
    edge = 'black' if cid in HL_C3_PLOT else 'white'
    lw   = 1.0     if cid in HL_C3_PLOT else 0.3
    bar_y_vals.append(lp)
    bar_colors.append(col)
    bar_edges.append(edge)
    bar_lws.append(lw)
    tick_labels.append(cid)

y_pos = np.arange(len(c3_order_flipped))
ax.barh(y_pos, bar_y_vals,
        color=bar_colors, edgecolor=bar_edges,
        linewidth=bar_lws, height=0.72)

ax.axvline(-np.log10(0.05), color='grey', ls=':', lw=0.8)
if fdr_cutoff_log_c3 is not None:
    ax.axvline(fdr_cutoff_log_c3, color='steelblue', ls='--', lw=1)

ax.legend(handles=[
    mpatches.Patch(fc=C3_GROUP_PLOT_COLORS['Group A'], label='Group A (LMX1A+)'),
    mpatches.Patch(fc=C3_GROUP_PLOT_COLORS['Group B'], label='Group B (LMX1A−)'),
    plt.Line2D([0],[0], color='steelblue', ls='--', lw=1,
               label=f'FDR<0.05 (p={fdr_cutoff_p_c3:.2f})'
                     if fdr_cutoff_log_c3 else 'FDR<0.05 (none)'),
    plt.Line2D([0],[0], color='grey', ls=':', lw=0.8, label='Nominal p=0.05'),
], fontsize=8, loc='lower right')

ax.set_yticks(y_pos)
ax.set_yticklabels(tick_labels, fontsize=8)
for lbl, cid in zip(ax.get_yticklabels(), c3_order_flipped):
    if cid in set(GRP_B_C3_PLOT):
        lbl.set_color('#7f8c8d')

ax.set_xlabel('−log₁₀(p)  LDSC  (Yengo2018)', fontsize=9)
ax.set_title('Lmx1a C3 clusters\nBMI enrichment  (CELLECT-LDSC)',
             fontweight='bold')
plt.tight_layout()
plt.savefig('LMX_C3_LDSC_sankey_order.pdf', bbox_inches='tight')
plt.savefig('LMX_C3_LDSC_sankey_order.png', bbox_inches='tight')
plt.show()


# In[5]:


# ── C3 tSNE: top-2 | expanded = corrected master classification ───────────────
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

UNIFIED_PATH = "hypomap_cellect_unified.csv"
GWAS         = 'BMI_Yengo2018'
BONF_C3      = 0.05 / 156

# Master constants already loaded — derive evidence tiers from them
# ── Evidence tiers (from master constants + BMI results) ─────────────────────
# Group A Bonf-sig MAGMA: C3-93, C3-84, C3-83, C3-82, C3-94
GRP_A_BONFSIG = ['C3-93','C3-84','C3-83','C3-82','C3-94']
# Group A nominal only: C3-85, C3-92
GRP_A_NOMINAL = ['C3-85','C3-92']
# CRNDE and GABA-neighbour: own categories
# CRNDE_C3 = ['C3-95'], GABA_NEIGHBOR_C3 = ['C3-28']
# Globally Bonf-sig, non-Lmx1a (LDSC confirmed)
GLOBAL_SIG_C3 = ['C3-118','C3-62','C3-123','C3-63','C3-100']

ALL_HIGHLIGHT_C3 = (GRP_A_BONFSIG + list(CRNDE_C3) +
                    list(GABA_NEIGHBOR_C3) + GRP_A_NOMINAL + GLOBAL_SIG_C3)

# ── Colour map aligned with master C3_GROUP_COLORS ───────────────────────────
C3_TSNE_COLORS = {
    # CRNDE — cross-species marker
    'C3-95':  C3_GROUP_COLORS['CRNDE'],           # '#e67e22'
    # Group A Bonf-sig — red gradient (darkest = HL_C3)
    'C3-93':  '#8B0000',   # HL_C3 — highest LMX1A + strongest BMI
    'C3-84':  '#B22222',   # HL_C3 — ZNF676
    'C3-83':  '#CD5C5C',
    'C3-82':  '#E07B7B',
    'C3-94':  '#EBA8A8',
    # GABA-neighbour — purple (specific but LMX1A-negative)
    'C3-28':  C3_GROUP_COLORS['GABA_neighbor'],   # '#8e44ad'
    # Group A nominal — muted orange
    'C3-85':  '#D4690A',
    'C3-92':  '#F0AD88',
    # Global Bonf-sig non-Lmx1a — blue gradient
    'C3-118': '#1a5276',
    'C3-62':  '#2e86c1',
    'C3-123': '#5dade2',
    'C3-63':  '#a9cce3',
    'C3-100': '#d6eaf8',
}

# ── Load BMI results for legend p-values ─────────────────────────────────────
unified  = pd.read_csv(UNIFIED_PATH)
magma_c3 = (unified[(unified['gwas']==GWAS) & (unified['method']=='new_MAGMA') &
                    (unified['level']=='C3') & unified['is_neuronal']]
            .set_index('cluster_id'))
ldsc_c3  = (unified[(unified['gwas']==GWAS) & (unified['method']=='new_LDSC') &
                    (unified['level']=='C3') & unified['is_neuronal']]
            .set_index('cluster_id'))

def plab(cid, df):
    return f'p={df.loc[cid,"pvalue"]:.1e}' if cid in df.index else ''

def siglab(cid):
    pm = magma_c3.loc[cid,'pvalue'] if cid in magma_c3.index else 1.0
    pl = ldsc_c3.loc[cid,'pvalue']  if cid in ldsc_c3.index  else 1.0
    bm = pm < BONF_C3; bl = pl < BONF_C3
    nm = pm < 0.05;    nl = pl < 0.05
    if bm and bl:  return '*** both'
    if bm and nl:  return '*** M  nom L'
    if bm:         return '*** M'
    if nm and nl:  return 'nom both'
    if nm:         return 'nom M'
    return ''

# ── Prepare mouse obs ─────────────────────────────────────────────────────────
if 'C3_short' not in adata_mouse.obs.columns:
    adata_mouse.obs['C3_short'] = (adata_mouse.obs['C3_HM']
                                   .str.extract(r'(C3-\d+)'))

# Top-2: HL_C3 = {'C3-84','C3-93'}
top2 = sorted(HL_C3)   # ['C3-84','C3-93']
if 'C3_top2' not in adata_mouse.obs.columns:
    adata_mouse.obs['C3_top2'] = adata_mouse.obs['C3_short'].where(
        adata_mouse.obs['C3_short'].isin(top2))
    adata_mouse.obs['C3_top2'] = pd.Categorical(
        adata_mouse.obs['C3_top2'], categories=top2, ordered=True)
    adata_mouse.uns['C3_top2_colors'] = ['#8B0000','#B22222']

# Expanded
adata_mouse.obs['C3_expanded'] = adata_mouse.obs['C3_short'].where(
    adata_mouse.obs['C3_short'].isin(ALL_HIGHLIGHT_C3))
adata_mouse.obs['C3_expanded'] = pd.Categorical(
    adata_mouse.obs['C3_expanded'],
    categories=ALL_HIGHLIGHT_C3, ordered=True)
adata_mouse.uns['C3_expanded_colors'] = [C3_TSNE_COLORS[c]
                                          for c in ALL_HIGHLIGHT_C3]

# ── Plot ──────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(20, 8))

sc.pl.embedding(adata_mouse, basis='tsne', color='C3_top2',
                legend_loc='on data', legend_fontsize='small',
                title=f'Top 2 C3 labels\n'
                      f'(C3-93 CBLN2, C3-84 ZNF676)',
                na_color='#eeeeee', ax=axes[0], show=False)

sc.pl.embedding(adata_mouse, basis='tsne', color='C3_expanded',
                legend_loc='on data', legend_fontsize=5,
                title='C3 BMI-enriched Lmx1a labels — corrected classification\n'
                      '(new_MAGMA + new_LDSC | Yengo2018)',
                na_color='#eeeeee', ax=axes[1], show=False)

ax = axes[1]
legend_el = [
    # CRNDE
    mpatches.Patch(fc='none', ec='none',
                   label='── CRNDE: cross-species marker ──'),
    mpatches.Patch(fc=C3_GROUP_COLORS['CRNDE'],
                   label=f'C3-95  CRNDE    {plab("C3-95",magma_c3)}  '
                         f'{siglab("C3-95")}  ★'),
    # Group A Bonf-sig
    mpatches.Patch(fc='none', ec='none',
                   label='── Group A: genuine LMX1A+ (MAGMA Bonf-sig) ──'),
    mpatches.Patch(fc='#8B0000',
                   label=f'C3-93  CBLN2    {plab("C3-93",magma_c3)}  '
                         f'{siglab("C3-93")}  ★'),
    mpatches.Patch(fc='#B22222',
                   label=f'C3-84  ZNF676   {plab("C3-84",magma_c3)}  '
                         f'{siglab("C3-84")}  ★'),
    mpatches.Patch(fc='#CD5C5C',
                   label=f'C3-83  ZEB2     {plab("C3-83",magma_c3)}  '
                         f'{siglab("C3-83")}'),
    mpatches.Patch(fc='#E07B7B',
                   label=f'C3-82  EBF2     {plab("C3-82",magma_c3)}  '
                         f'{siglab("C3-82")}'),
    mpatches.Patch(fc='#EBA8A8',
                   label=f'C3-94  THSD7B   {plab("C3-94",magma_c3)}  '
                         f'{siglab("C3-94")}'),
    # GABA-neighbour
    mpatches.Patch(fc='none', ec='none',
                   label='── GABA-neighbour: specific, LMX1A-negative ──'),
    mpatches.Patch(fc=C3_GROUP_COLORS['GABA_neighbor'],
                   label=f'C3-28  KCNH8    {plab("C3-28",magma_c3)}  '
                         f'{siglab("C3-28")}  (Pre-1 GABA)'),
    # Group A nominal
    mpatches.Patch(fc='none', ec='none',
                   label='── Group A: genuine LMX1A+ (nominal) ──'),
    mpatches.Patch(fc='#D4690A',
                   label=f'C3-85  EBF2     {plab("C3-85",magma_c3)}  '
                         f'{siglab("C3-85")}'),
    mpatches.Patch(fc='#F0AD88',
                   label=f'C3-92  AGTR1    {plab("C3-92",magma_c3)}  '
                         f'{siglab("C3-92")}'),
    # Global contrast
    mpatches.Patch(fc='none', ec='none',
                   label='── Global Bonf-sig non-Lmx1a (LDSC) ──'),
    mpatches.Patch(fc='#1a5276',
                   label=f'C3-118 CALB1    {plab("C3-118",ldsc_c3)}  LDSC***'),
    mpatches.Patch(fc='#2e86c1',
                   label=f'C3-62  ADAMTSL1 {plab("C3-62",ldsc_c3)}  LDSC***'),
    mpatches.Patch(fc='#5dade2',
                   label=f'C3-123 RGS22    {plab("C3-123",ldsc_c3)}  LDSC***'),
    mpatches.Patch(fc='#a9cce3',
                   label=f'C3-63  TRHDE    {plab("C3-63",ldsc_c3)}  LDSC***'),
    mpatches.Patch(fc='#d6eaf8', ec='#aaaaaa',
                   label=f'C3-100 OTP      {plab("C3-100",ldsc_c3)}  LDSC***'),
]
ax.legend(handles=legend_el, fontsize=6.5, loc='lower left',
          bbox_to_anchor=(1.01, 0), framealpha=0.95,
          title=('★ = HL_C3 highlight or CRNDE cross-species\n'
                 '*** = Bonferroni sig  nom = p<0.05\n'
                 'M = new_MAGMA  L = new_LDSC\n'
                 'C3-115 removed (Group B: 5.5% specific)\n'
                 'C3-29/137 removed (Group B: <32% specific)'),
          title_fontsize=6.5)

plt.suptitle('Human C3 Hypomap labels on mouse Lmx1a tSNE — corrected classification\n'
             'Group A (LMX1A+) | CRNDE | GABA-neighbour | Global contrast\n'
             'GWAS: BMI_Yengo2018',
             fontweight='bold', fontsize=12, y=1.01)
plt.tight_layout()
plt.savefig('LMX_C3_tsne_expanded_corrected.pdf', bbox_inches='tight')
plt.savefig('LMX_C3_tsne_expanded_corrected.png', bbox_inches='tight')
plt.show()

# ── Cell count summary ────────────────────────────────────────────────────────
MOUSE_CELLTYPES = ['04_Crnde_Lmx1a','27_Lmx1a_Ttc6',
                   '32_Gucy2c_Ttc6','36_Ebf2_Crhbp']
print("\nCell counts per highlighted C3 cluster (mouse atlas):")
print(f"{'C3':^10} {'total':>8} {'Lmx1a':>8} {'%':>6}  group         evidence")
print("-" * 65)
for cid in ALL_HIGHLIGHT_C3:
    total = (adata_mouse.obs['C3_short'] == cid).sum()
    lmx   = (adata_mouse.obs[
                adata_mouse.obs['detailed_celltype_annotation']
                .isin(MOUSE_CELLTYPES)
             ]['C3_short'] == cid).sum()
    pct   = lmx/total*100 if total > 0 else 0
    grp   = ('GrpA-Bonfsig'  if cid in GRP_A_BONFSIG  else
             'CRNDE'          if cid in set(CRNDE_C3)   else
             'GABA-nb'        if cid in set(GABA_NEIGHBOR_C3) else
             'GrpA-nominal'   if cid in GRP_A_NOMINAL  else
             'Global')
    ev    = siglab(cid)
    print(f"  {cid:<10} {total:>6}  {lmx:>6}  {pct:>5.1f}%  "
          f"[{grp:<14}]  {ev}")


# In[6]:


# ── C4 tSNE: top-2 | expanded = corrected master classification ───────────────
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

UNIFIED_PATH = "hypomap_cellect_unified.csv"
GWAS         = 'BMI_Yengo2018'
BONF_C4      = 0.05 / 415

# ── Load BMI results ──────────────────────────────────────────────────────────
unified  = pd.read_csv(UNIFIED_PATH)

def get_c4(method):
    return (unified[(unified['gwas']==GWAS) & (unified['method']==method) &
                    (unified['level']=='C4') & (unified['is_neuronal'])]
            .set_index('cluster_id'))

paper_c4 = get_c4('paper_MAGMA')
magma_c4 = get_c4('new_MAGMA')
ldsc_c4  = get_c4('new_LDSC')

# ── Evidence tiers — derived from master GROUP_A_C4 ──────────────────────────
# Tier 1: paper Bonf-sig + new_MAGMA Bonf-sig
TIER1_LMX_C4 = ['C4-229','C4-247','C4-252']
# Tier 2: paper Bonf-sig only
TIER2_LMX_C4 = ['C4-223','C4-224','C4-228','C4-230']
# Tier 3: Group A nominal (paper + LDSC)
TIER3_LMX_C4 = ['C4-226','C4-231','C4-248']
# Master GROUP_B_C4 already loaded: ['C4-33','C4-407','C4-409','C4-411']
# Globally Bonf-sig paper + LDSC, non-Lmx1a contrast
GLOBAL_SIG_C4 = ['C4-333','C4-334','C4-150','C4-337','C4-152',
                  'C4-338','C4-356','C4-279','C4-285','C4-363']

ALL_HIGHLIGHT_C4 = (TIER1_LMX_C4 + TIER2_LMX_C4 + TIER3_LMX_C4 +
                    GROUP_B_C4 + GLOBAL_SIG_C4)

# ── Colour map (named C4_TSNE_COLORS to avoid shadowing C4_GROUP_COLORS) ─────
C4_TSNE_COLORS = {
    # Tier 1 — darkest red (HL_C4 = {'C4-229','C4-248'})
    'C4-229': C4_GROUP_COLORS['Group A'],   # '#8B0000'
    'C4-247': '#B22222',
    'C4-252': '#CD5C5C',
    # Tier 2 — orange gradient
    'C4-223': '#D4690A',
    'C4-224': '#E07E30',
    'C4-228': '#E89560',
    'C4-230': '#F0AD88',
    # Tier 3 — light red/pink; C4-248 is HL but weaker signal
    'C4-226': '#C44E52',
    'C4-231': '#DA7F7F',
    'C4-248': '#E8A8A8',
    # Group B — grey gradient (master C4_GROUP_COLORS['Group B'] = '#95a5a6')
    'C4-407': C4_GROUP_COLORS['Group B'],   # '#95a5a6'
    'C4-33':  '#aab7b8',
    'C4-409': '#bdc3c7',
    'C4-411': '#d5d8dc',
    # Global non-Lmx1a contrast — blue gradient
    'C4-333': '#1a5276', 'C4-334': '#1f618d', 'C4-150': '#2874a6',
    'C4-337': '#2e86c1', 'C4-152': '#3498db', 'C4-338': '#5dade2',
    'C4-356': '#7fb3d3', 'C4-279': '#85c1e9', 'C4-285': '#aed6f1',
    'C4-363': '#d6eaf8',
}

# ── Prepare mouse obs ─────────────────────────────────────────────────────────
if 'C4_short' not in adata_mouse.obs.columns:
    for col in ['C4_HM','celltype_HM']:
        if col in adata_mouse.obs.columns:
            extracted = adata_mouse.obs[col].str.extract(r'(C4-\d+)')
            if extracted.notna().any().any():
                adata_mouse.obs['C4_short'] = extracted
                print(f"C4_short extracted from {col}")
                break
    else:
        adata_mouse.obs['C4_short'] = pd.NA
        print("WARNING: C4_short not available")

print(f"Mouse C4_short non-null: {adata_mouse.obs['C4_short'].notna().sum()}")

# Top-2: HL_C4 = {'C4-229','C4-248'}
TOP2_C4 = ['C4-229','C4-248']
adata_mouse.obs['C4_top2'] = adata_mouse.obs['C4_short'].where(
    adata_mouse.obs['C4_short'].isin(TOP2_C4))
adata_mouse.obs['C4_top2'] = pd.Categorical(
    adata_mouse.obs['C4_top2'], categories=TOP2_C4, ordered=True)
adata_mouse.uns['C4_top2_colors'] = [C4_TSNE_COLORS['C4-229'],
                                      C4_TSNE_COLORS['C4-248']]

# Expanded
adata_mouse.obs['C4_expanded'] = adata_mouse.obs['C4_short'].where(
    adata_mouse.obs['C4_short'].isin(ALL_HIGHLIGHT_C4))
adata_mouse.obs['C4_expanded'] = pd.Categorical(
    adata_mouse.obs['C4_expanded'],
    categories=ALL_HIGHLIGHT_C4, ordered=True)
adata_mouse.uns['C4_expanded_colors'] = [C4_TSNE_COLORS[c]
                                          for c in ALL_HIGHLIGHT_C4]

# ── Legend helpers ────────────────────────────────────────────────────────────
def plab(cid, df):
    return f'p={df.loc[cid,"pvalue"]:.1e}' if cid in df.index else ''

def tier_label(cid):
    pp = paper_c4.loc[cid,'pvalue'] if cid in paper_c4.index else 1.0
    pm = magma_c4.loc[cid,'pvalue'] if cid in magma_c4.index else 1.0
    pl = ldsc_c4.loc[cid,'pvalue']  if cid in ldsc_c4.index  else 1.0
    if pp < BONF_C4 and pm < BONF_C4: return '***P+M'
    if pp < BONF_C4 and pl < 0.05:    return '***P nomL'
    if pp < BONF_C4:                   return '***P'
    if pp < 0.05 and pl < 0.05:       return 'nom P+L'
    return 'ns'

# ── Plot ──────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(20, 8))

sc.pl.embedding(adata_mouse, basis='tsne', color='C4_top2',
                legend_loc='on data', legend_fontsize='small',
                title='Top 2 C4 labels\n'
                      'C4-229 ZNF676/CDH23 (strongest Bonf-sig, 97% spec)\n'
                      'C4-248 AGTR1/EPB41L4A (canonical LMX1A+, 98% spec)',
                na_color='#eeeeee', ax=axes[0], show=False)

sc.pl.embedding(adata_mouse, basis='tsne', color='C4_expanded',
                legend_loc='on data', legend_fontsize=4,
                title='C4 BMI-enriched Lmx1a labels — corrected classification\n'
                      'Evidence tiers + global contrast  |  GWAS: BMI_Yengo2018',
                na_color='#eeeeee', ax=axes[1], show=False)

ax = axes[1]
legend_el = [
    mpatches.Patch(fc='none', ec='none',
                   label='── Tier 1: Bonf-sig paper + new_MAGMA ──'),
    mpatches.Patch(fc=C4_TSNE_COLORS['C4-229'],
                   label=f'C4-229  ZNF676/CDH23    {plab("C4-229",paper_c4)}'
                         f'  {tier_label("C4-229")}  ★'),
    mpatches.Patch(fc=C4_TSNE_COLORS['C4-247'],
                   label=f'C4-247  AGTR1/EBF2      {plab("C4-247",paper_c4)}'
                         f'  {tier_label("C4-247")}'),
    mpatches.Patch(fc=C4_TSNE_COLORS['C4-252'],
                   label=f'C4-252  THSD7B/PTH2R    {plab("C4-252",paper_c4)}'
                         f'  {tier_label("C4-252")}'),
    mpatches.Patch(fc='none', ec='none',
                   label='── Tier 2: Bonf-sig paper only ──'),
    mpatches.Patch(fc=C4_TSNE_COLORS['C4-223'],
                   label=f'C4-223  ZNF676/EBF2     {plab("C4-223",paper_c4)}'
                         f'  {tier_label("C4-223")}'),
    mpatches.Patch(fc=C4_TSNE_COLORS['C4-224'],
                   label=f'C4-224  ZNF676/NXPH2    {plab("C4-224",paper_c4)}'
                         f'  {tier_label("C4-224")}'),
    mpatches.Patch(fc=C4_TSNE_COLORS['C4-228'],
                   label=f'C4-228  ZNF676/POU6F2   {plab("C4-228",paper_c4)}'
                         f'  {tier_label("C4-228")}'),
    mpatches.Patch(fc=C4_TSNE_COLORS['C4-230'],
                   label=f'C4-230  EBF2/IL1RAPL2   {plab("C4-230",paper_c4)}'
                         f'  {tier_label("C4-230")}'),
    mpatches.Patch(fc='none', ec='none',
                   label='── Tier 3: Group A nominal (paper + LDSC) ──'),
    mpatches.Patch(fc=C4_TSNE_COLORS['C4-226'],
                   label=f'C4-226  ZNF676/IRX6     {plab("C4-226",paper_c4)}'
                         f'  {tier_label("C4-226")}'),
    mpatches.Patch(fc=C4_TSNE_COLORS['C4-231'],
                   label=f'C4-231  EBF2/PLEKHH2    {plab("C4-231",paper_c4)}'
                         f'  {tier_label("C4-231")}'),
    mpatches.Patch(fc=C4_TSNE_COLORS['C4-248'],
                   label=f'C4-248  AGTR1/EPB41L4A  {plab("C4-248",paper_c4)}'
                         f'  {tier_label("C4-248")}  ★'),
    mpatches.Patch(fc='none', ec='none',
                   label='── Group B: contaminant labels ──'),
    mpatches.Patch(fc=C4_TSNE_COLORS['C4-407'],
                   label=f'C4-407  ADGRL4/IGF2BP2  {plab("C4-407",paper_c4)}'
                         f'  22% spec'),
    mpatches.Patch(fc=C4_TSNE_COLORS['C4-33'],
                   label=f'C4-33   EYA2/CXCL14     {plab("C4-33",paper_c4)}'
                         f'  28% spec'),
    mpatches.Patch(fc=C4_TSNE_COLORS['C4-409'],
                   label=f'C4-409  KIT/MCTP2        {plab("C4-409",paper_c4)}'
                         f'  22% spec'),
    mpatches.Patch(fc=C4_TSNE_COLORS['C4-411'],
                   label=f'C4-411  SALL3/PVALB      {plab("C4-411",paper_c4)}'
                         f'  16% spec'),
    mpatches.Patch(fc='none', ec='none',
                   label='── Global: Bonf-sig paper + LDSC (non-Lmx1a) ──'),
    mpatches.Patch(fc=C4_TSNE_COLORS['C4-333'],
                   label=f'C4-333  CALB1/ONECUT1   {plab("C4-333",paper_c4)}  ***P+L'),
    mpatches.Patch(fc=C4_TSNE_COLORS['C4-334'],
                   label=f'C4-334  CALB1/TINAG     {plab("C4-334",paper_c4)}  ***P+L'),
    mpatches.Patch(fc=C4_TSNE_COLORS['C4-150'],
                   label=f'C4-150  ADAMTSL1/HMCN1  {plab("C4-150",paper_c4)}  ***P+L'),
    mpatches.Patch(fc=C4_TSNE_COLORS['C4-337'],
                   label=f'C4-337  CALB1/SATB2     {plab("C4-337",paper_c4)}  ***P+L'),
    mpatches.Patch(fc=C4_TSNE_COLORS['C4-152'],
                   label=f'C4-152  TRHDE/CALCR     {plab("C4-152",paper_c4)}  ***P+L'),
    mpatches.Patch(fc=C4_TSNE_COLORS['C4-338'],
                   label=f'C4-338  CALB1/RASSF6    {plab("C4-338",paper_c4)}  ***P+L'),
    mpatches.Patch(fc=C4_TSNE_COLORS['C4-356'],
                   label=f'C4-356  RGS22/GPC3      {plab("C4-356",paper_c4)}  ***P+L'),
    mpatches.Patch(fc=C4_TSNE_COLORS['C4-279'],
                   label=f'C4-279  CBLN2/ONECUT3   {plab("C4-279",paper_c4)}  ***P+L'),
    mpatches.Patch(fc=C4_TSNE_COLORS['C4-285'],
                   label=f'C4-285  OTP/NXPH2       {plab("C4-285",paper_c4)}  ***P+L'),
    mpatches.Patch(fc=C4_TSNE_COLORS['C4-363'], ec='#aaaaaa',
                   label=f'C4-363  PPP1R17/ADGRG6  {plab("C4-363",paper_c4)}  ***P+L'),
]
ax.legend(handles=legend_el, fontsize=6.5, loc='lower left',
          bbox_to_anchor=(1.01, 0), framealpha=0.95,
          title=('★ = HL_C4 highlight\n'
                 '***P = paper Bonf-sig  ***M = new_MAGMA Bonf-sig\n'
                 '***L = new_LDSC Bonf-sig  nom = p<0.05\n'
                 'spec = % mouse Lmx1a specificity\n'
                 'C4-7 excluded (1.8% specific — noise label)'),
          title_fontsize=6.5)

plt.suptitle('Human C4 Hypomap labels on mouse Lmx1a tSNE — corrected classification\n'
             'Group A tiers | Group B contaminants | Global contrast\n'
             'GWAS: BMI_Yengo2018',
             fontweight='bold', fontsize=12, y=1.01)
plt.tight_layout()
plt.savefig('LMX_C4_tsne_expanded_corrected.pdf', bbox_inches='tight')
plt.savefig('LMX_C4_tsne_expanded_corrected.png', bbox_inches='tight')
plt.show()

# ── Cell count summary ────────────────────────────────────────────────────────
MOUSE_CELLTYPES = ['04_Crnde_Lmx1a','27_Lmx1a_Ttc6',
                   '32_Gucy2c_Ttc6','36_Ebf2_Crhbp']
print("\nCell counts per highlighted C4 cluster (mouse atlas):")
print(f"{'C4':^10} {'total':>7} {'Lmx1a':>7} {'%':>6}  tier     p_paper    annotation")
print("-" * 80)
for tier, ids in [('Tier1',  TIER1_LMX_C4),
                   ('Tier2',  TIER2_LMX_C4),
                   ('Tier3',  TIER3_LMX_C4),
                   ('GrpB',   GROUP_B_C4),
                   ('Global', GLOBAL_SIG_C4[:5])]:
    for cid in ids:
        total = (adata_mouse.obs['C4_short'] == cid).sum()
        lmx   = (adata_mouse.obs[
                    adata_mouse.obs['detailed_celltype_annotation']
                    .isin(MOUSE_CELLTYPES)
                 ]['C4_short'] == cid).sum()
        pct  = lmx/total*100 if total > 0 else 0
        ann  = paper_c4.loc[cid,'annotation_full'] if cid in paper_c4.index else ''
        pp   = paper_c4.loc[cid,'pvalue']          if cid in paper_c4.index else np.nan
        print(f"  {cid:<10} {total:>6}  {lmx:>6}  {pct:>5.1f}%  "
              f"[{tier:<7}] p={pp:.2e}  {ann[4:40]}")


# In[7]:


# ── C3 + C4 on human HYPOMAP UMAP — corrected master classification ───────────
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Master constants already loaded. Local colour maps use prefixed names
# to avoid shadowing C3_GROUP_COLORS / C4_GROUP_COLORS.

# ── C3 tier definitions — corrected ──────────────────────────────────────────
# Group A Bonf-sig MAGMA (genuine LMX1A+)
GRP_A_BONFSIG_C3 = ['C3-93','C3-84','C3-83','C3-82','C3-94']
# Group A nominal
GRP_A_NOMINAL_C3 = ['C3-85','C3-92']
# CRNDE_C3 = ['C3-95']       — cross-species marker
# GABA_NEIGHBOR_C3 = ['C3-28'] — specific but LMX1A-negative
GLOBAL_SIG_C3    = ['C3-118','C3-62','C3-123','C3-63','C3-100']

ALL_C3  = (GRP_A_BONFSIG_C3 + list(CRNDE_C3) +
           list(GABA_NEIGHBOR_C3) + GRP_A_NOMINAL_C3 + GLOBAL_SIG_C3)
TOP2_C3 = sorted(HL_C3)   # ['C3-84','C3-93']

C3_UMAP_COLORS = {
    # Group A Bonf-sig — red gradient; HL_C3 darkest
    'C3-93':  C3_GROUP_COLORS['Group A'],   # '#8B0000'  HL
    'C3-84':  '#B22222',                    # HL
    'C3-83':  '#CD5C5C',
    'C3-82':  '#E07B7B',
    'C3-94':  '#EBA8A8',
    # CRNDE — orange
    'C3-95':  C3_GROUP_COLORS['CRNDE'],     # '#e67e22'
    # GABA-neighbour — purple
    'C3-28':  C3_GROUP_COLORS['GABA_neighbor'],  # '#8e44ad'
    # Group A nominal — muted orange
    'C3-85':  '#D4690A',
    'C3-92':  '#F0AD88',
    # Global non-Lmx1a — blue gradient
    'C3-118': '#1a5276', 'C3-62': '#2e86c1', 'C3-123': '#5dade2',
    'C3-63':  '#a9cce3', 'C3-100':'#d6eaf8',
}

# ── C4 tier definitions — same as corrected tSNE cell ────────────────────────
TIER1_LMX_C4  = ['C4-229','C4-247','C4-252']
TIER2_LMX_C4  = ['C4-223','C4-224','C4-228','C4-230']
TIER3_LMX_C4  = ['C4-226','C4-231','C4-248']
# GROUP_B_C4 from master constants
GLOBAL_SIG_C4 = ['C4-333','C4-334','C4-150','C4-337','C4-152',
                  'C4-338','C4-356','C4-279','C4-285','C4-363']
ALL_C4  = (TIER1_LMX_C4 + TIER2_LMX_C4 + TIER3_LMX_C4 +
           GROUP_B_C4 + GLOBAL_SIG_C4)
TOP2_C4 = sorted(HL_C4)   # ['C4-229','C4-248']

C4_UMAP_COLORS = {
    'C4-229': C4_GROUP_COLORS['Group A'],   # '#8B0000'  HL
    'C4-247': '#B22222',
    'C4-252': '#CD5C5C',
    'C4-223': '#D4690A', 'C4-224': '#E07E30',
    'C4-228': '#E89560', 'C4-230': '#F0AD88',
    'C4-226': '#C44E52', 'C4-231': '#DA7F7F',
    'C4-248': '#E8A8A8',
    'C4-407': C4_GROUP_COLORS['Group B'],   # '#95a5a6'
    'C4-33':  '#aab7b8', 'C4-409': '#bdc3c7', 'C4-411': '#d5d8dc',
    'C4-333': '#1a5276', 'C4-334': '#1f618d', 'C4-150': '#2874a6',
    'C4-337': '#2e86c1', 'C4-152': '#3498db', 'C4-338': '#5dade2',
    'C4-356': '#7fb3d3', 'C4-279': '#85c1e9', 'C4-285': '#aed6f1',
    'C4-363': '#d6eaf8',
}

# ── Verify short IDs on human adata ──────────────────────────────────────────
for level, col_short, col_named in [('C3','C3_short','C3'),
                                      ('C4','C4_short','C4')]:
    if col_short not in adata.obs.columns:
        adata.obs[col_short] = adata.obs[col_named].str.extract(
            rf'({level}-\d+)')
        print(f"Built {col_short} from {col_named}")

print(f"Human C3_short non-null: {adata.obs['C3_short'].notna().sum()}")
print(f"Human C4_short non-null: {adata.obs['C4_short'].notna().sum()}")

# ── Build categorical columns ─────────────────────────────────────────────────
def make_cat(obs_col, ids, color_map):
    adata.obs[obs_col] = adata.obs[obs_col.replace('_human','_short')
                                   if '_human' in obs_col
                                   else obs_col].where(
        adata.obs[obs_col.replace('_human','_short')
                  if '_human' in obs_col else obs_col].isin(ids))
    adata.obs[obs_col] = pd.Categorical(
        adata.obs[obs_col], categories=ids, ordered=True)
    adata.uns[f'{obs_col}_colors'] = [color_map[c] for c in ids]

# C3 top-2
adata.obs['C3_top2_human'] = adata.obs['C3_short'].where(
    adata.obs['C3_short'].isin(TOP2_C3))
adata.obs['C3_top2_human'] = pd.Categorical(
    adata.obs['C3_top2_human'], categories=TOP2_C3, ordered=True)
adata.uns['C3_top2_human_colors'] = [C3_UMAP_COLORS[c] for c in TOP2_C3]

# C3 expanded
adata.obs['C3_expanded_human'] = adata.obs['C3_short'].where(
    adata.obs['C3_short'].isin(ALL_C3))
adata.obs['C3_expanded_human'] = pd.Categorical(
    adata.obs['C3_expanded_human'], categories=ALL_C3, ordered=True)
adata.uns['C3_expanded_human_colors'] = [C3_UMAP_COLORS[c] for c in ALL_C3]

# C4 top-2
adata.obs['C4_top2_human'] = adata.obs['C4_short'].where(
    adata.obs['C4_short'].isin(TOP2_C4))
adata.obs['C4_top2_human'] = pd.Categorical(
    adata.obs['C4_top2_human'], categories=TOP2_C4, ordered=True)
adata.uns['C4_top2_human_colors'] = [C4_UMAP_COLORS[c] for c in TOP2_C4]

# C4 expanded
adata.obs['C4_expanded_human'] = adata.obs['C4_short'].where(
    adata.obs['C4_short'].isin(ALL_C4))
adata.obs['C4_expanded_human'] = pd.Categorical(
    adata.obs['C4_expanded_human'], categories=ALL_C4, ordered=True)
adata.uns['C4_expanded_human_colors'] = [C4_UMAP_COLORS[c] for c in ALL_C4]

# ── Load BMI results for legend p-values ─────────────────────────────────────
unified  = pd.read_csv("hypomap_cellect_unified.csv")
GWAS     = 'BMI_Yengo2018'
BONF_C3  = 0.05 / 156
BONF_C4  = 0.05 / 415

def get_c(method, level):
    return (unified[(unified['gwas']==GWAS) & (unified['method']==method) &
                    (unified['level']==level) & unified['is_neuronal']]
            .set_index('cluster_id'))

paper_c4 = get_c('paper_MAGMA','C4')
magma_c3 = get_c('new_MAGMA',  'C3')
ldsc_c3  = get_c('new_LDSC',   'C3')
ldsc_c4  = get_c('new_LDSC',   'C4')

def plab(cid, df):
    return f'p={df.loc[cid,"pvalue"]:.1e}' if cid in df.index else ''

# ── Legend builders — corrected 4-way C3 classification ──────────────────────
def c3_legend():
    return [
        # Group A Bonf-sig
        mpatches.Patch(fc='none', ec='none',
                       label='── Group A genuine LMX1A+ (MAGMA Bonf-sig) ──'),
        mpatches.Patch(fc=C3_UMAP_COLORS['C3-93'],
                       label=f'C3-93  CBLN2    {plab("C3-93",magma_c3)}'
                             f'  ***M nomL  ★'),
        mpatches.Patch(fc=C3_UMAP_COLORS['C3-84'],
                       label=f'C3-84  ZNF676   {plab("C3-84",magma_c3)}'
                             f'  ***M nomL  ★'),
        mpatches.Patch(fc=C3_UMAP_COLORS['C3-83'],
                       label=f'C3-83  ZEB2     {plab("C3-83",magma_c3)}  ***M'),
        mpatches.Patch(fc=C3_UMAP_COLORS['C3-82'],
                       label=f'C3-82  EBF2     {plab("C3-82",magma_c3)}  ***M'),
        mpatches.Patch(fc=C3_UMAP_COLORS['C3-94'],
                       label=f'C3-94  THSD7B   {plab("C3-94",magma_c3)}  ***M'),
        # CRNDE
        mpatches.Patch(fc='none', ec='none',
                       label='── CRNDE: cross-species marker ──'),
        mpatches.Patch(fc=C3_UMAP_COLORS['C3-95'],
                       label=f'C3-95  CRNDE    {plab("C3-95",magma_c3)}'
                             f'  ***M nomL  ★'),
        # GABA-neighbour
        mpatches.Patch(fc='none', ec='none',
                       label='── GABA-neighbour: specific, LMX1A-negative ──'),
        mpatches.Patch(fc=C3_UMAP_COLORS['C3-28'],
                       label=f'C3-28  KCNH8    {plab("C3-28",magma_c3)}'
                             f'  ***M nomL  (Pre-1 GABA)'),
        # Group A nominal
        mpatches.Patch(fc='none', ec='none',
                       label='── Group A genuine LMX1A+ (nominal) ──'),
        mpatches.Patch(fc=C3_UMAP_COLORS['C3-85'],
                       label=f'C3-85  EBF2     {plab("C3-85",magma_c3)}  nom both'),
        mpatches.Patch(fc=C3_UMAP_COLORS['C3-92'],
                       label=f'C3-92  AGTR1    {plab("C3-92",magma_c3)}  nom both'),
        # Global contrast
        mpatches.Patch(fc='none', ec='none',
                       label='── Global Bonf-sig non-Lmx1a (LDSC) ──'),
        mpatches.Patch(fc='#1a5276',
                       label=f'C3-118 CALB1    {plab("C3-118",ldsc_c3)}  ***L'),
        mpatches.Patch(fc='#2e86c1',
                       label=f'C3-62  ADAMTSL1 {plab("C3-62",ldsc_c3)}  ***L'),
        mpatches.Patch(fc='#5dade2',
                       label=f'C3-123 RGS22    {plab("C3-123",ldsc_c3)}  ***L'),
        mpatches.Patch(fc='#a9cce3',
                       label=f'C3-63  TRHDE    {plab("C3-63",ldsc_c3)}  ***L'),
        mpatches.Patch(fc='#d6eaf8', ec='#aaa',
                       label=f'C3-100 OTP      {plab("C3-100",ldsc_c3)}  ***L'),
    ]

def c4_legend():
    return [
        mpatches.Patch(fc='none', ec='none',
                       label='── Tier 1: Bonf-sig paper + new_MAGMA ──'),
        mpatches.Patch(fc=C4_UMAP_COLORS['C4-229'],
                       label=f'C4-229  ZNF676/CDH23    {plab("C4-229",paper_c4)}'
                             f'  ***P+M  ★'),
        mpatches.Patch(fc=C4_UMAP_COLORS['C4-247'],
                       label=f'C4-247  AGTR1/EBF2      {plab("C4-247",paper_c4)}'
                             f'  ***P+M'),
        mpatches.Patch(fc=C4_UMAP_COLORS['C4-252'],
                       label=f'C4-252  THSD7B/PTH2R    {plab("C4-252",paper_c4)}'
                             f'  ***P+M'),
        mpatches.Patch(fc='none', ec='none',
                       label='── Tier 2: Bonf-sig paper only ──'),
        mpatches.Patch(fc=C4_UMAP_COLORS['C4-223'],
                       label=f'C4-223  ZNF676/EBF2     {plab("C4-223",paper_c4)}'
                             f'  ***P'),
        mpatches.Patch(fc=C4_UMAP_COLORS['C4-224'],
                       label=f'C4-224  ZNF676/NXPH2    {plab("C4-224",paper_c4)}'
                             f'  ***P'),
        mpatches.Patch(fc=C4_UMAP_COLORS['C4-228'],
                       label=f'C4-228  ZNF676/POU6F2   {plab("C4-228",paper_c4)}'
                             f'  ***P'),
        mpatches.Patch(fc=C4_UMAP_COLORS['C4-230'],
                       label=f'C4-230  EBF2/IL1RAPL2   {plab("C4-230",paper_c4)}'
                             f'  ***P'),
        mpatches.Patch(fc='none', ec='none',
                       label='── Tier 3: Group A nominal ──'),
        mpatches.Patch(fc=C4_UMAP_COLORS['C4-226'],
                       label=f'C4-226  ZNF676/IRX6     {plab("C4-226",paper_c4)}'
                             f'  nom P+L'),
        mpatches.Patch(fc=C4_UMAP_COLORS['C4-231'],
                       label=f'C4-231  EBF2/PLEKHH2    {plab("C4-231",paper_c4)}'
                             f'  nom P+L'),
        mpatches.Patch(fc=C4_UMAP_COLORS['C4-248'],
                       label=f'C4-248  AGTR1/EPB41L4A  {plab("C4-248",paper_c4)}'
                             f'  nom P+L  ★'),
        mpatches.Patch(fc='none', ec='none',
                       label='── Group B: contaminant labels ──'),
        mpatches.Patch(fc=C4_UMAP_COLORS['C4-407'],
                       label=f'C4-407  ADGRL4/IGF2BP2  {plab("C4-407",paper_c4)}'
                             f'  22% spec'),
        mpatches.Patch(fc=C4_UMAP_COLORS['C4-33'],
                       label=f'C4-33   EYA2/CXCL14     {plab("C4-33",paper_c4)}'
                             f'  28% spec'),
        mpatches.Patch(fc=C4_UMAP_COLORS['C4-409'],
                       label=f'C4-409  KIT/MCTP2        {plab("C4-409",paper_c4)}'
                             f'  22% spec'),
        mpatches.Patch(fc=C4_UMAP_COLORS['C4-411'],
                       label=f'C4-411  SALL3/PVALB      {plab("C4-411",paper_c4)}'
                             f'  16% spec'),
        mpatches.Patch(fc='none', ec='none',
                       label='── Global: Bonf-sig paper + LDSC (non-Lmx1a) ──'),
        mpatches.Patch(fc='#1a5276',
                       label=f'C4-333  CALB1/ONECUT1   {plab("C4-333",paper_c4)}  ***P+L'),
        mpatches.Patch(fc='#1f618d',
                       label=f'C4-334  CALB1/TINAG     {plab("C4-334",paper_c4)}  ***P+L'),
        mpatches.Patch(fc='#2874a6',
                       label=f'C4-150  ADAMTSL1/HMCN1  {plab("C4-150",paper_c4)}  ***P+L'),
        mpatches.Patch(fc='#2e86c1',
                       label=f'C4-337  CALB1/SATB2     {plab("C4-337",paper_c4)}  ***P+L'),
        mpatches.Patch(fc='#3498db',
                       label=f'C4-152  TRHDE/CALCR     {plab("C4-152",paper_c4)}  ***P+L'),
        mpatches.Patch(fc='#5dade2',
                       label=f'C4-338  CALB1/RASSF6    {plab("C4-338",paper_c4)}  ***P+L'),
        mpatches.Patch(fc='#7fb3d3',
                       label=f'C4-356  RGS22/GPC3      {plab("C4-356",paper_c4)}  ***P+L'),
        mpatches.Patch(fc='#85c1e9',
                       label=f'C4-279  CBLN2/ONECUT3   {plab("C4-279",paper_c4)}  ***P+L'),
        mpatches.Patch(fc='#aed6f1',
                       label=f'C4-285  OTP/NXPH2       {plab("C4-285",paper_c4)}  ***P+L'),
        mpatches.Patch(fc='#d6eaf8', ec='#aaa',
                       label=f'C4-363  PPP1R17/ADGRG6  {plab("C4-363",paper_c4)}  ***P+L'),
    ]

# ── Figure: 2 rows × 2 cols ───────────────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(22, 18))
LEG_KW = dict(fontsize=6.5, loc='lower left',
              bbox_to_anchor=(1.01, 0), framealpha=0.95)

# Row 0: C3
sc.pl.embedding(adata, basis='X_umap', color='C3_top2_human',
                legend_loc='on data', legend_fontsize='small',
                title='C3 top-2 labels — human HYPOMAP\n'
                      f'C3-93 CBLN2  |  C3-84 ZNF676',
                na_color='#eeeeee', ax=axes[0, 0], show=False)
axes[0, 0].set_xlabel('UMAP 1'); axes[0, 0].set_ylabel('UMAP 2')

sc.pl.embedding(adata, basis='X_umap', color='C3_expanded_human',
                legend_loc='on data', legend_fontsize=4,
                title='C3 BMI-enriched Lmx1a labels — human HYPOMAP\n'
                      'Group A | CRNDE | GABA-nb | Global  |  '
                      'GWAS: BMI_Yengo2018',
                na_color='#eeeeee', ax=axes[0, 1], show=False)
axes[0, 1].set_xlabel('UMAP 1'); axes[0, 1].set_ylabel('UMAP 2')
axes[0, 1].legend(
    handles=c3_legend(),
    title=('★ = HL_C3 highlight or CRNDE cross-species\n'
           '***M = new_MAGMA Bonf-sig  ***L = new_LDSC\n'
           'nom = p<0.05  |  C3-115/29/137 → Group B\n'
           'C3-28 → GABA-neighbour (specific, LMX1A-neg)'),
    title_fontsize=6.5, **LEG_KW)

# Row 1: C4
sc.pl.embedding(adata, basis='X_umap', color='C4_top2_human',
                legend_loc='on data', legend_fontsize='small',
                title='C4 top-2 labels — human HYPOMAP\n'
                      'C4-229 ZNF676/CDH23 (Bonf-sig)  |  '
                      'C4-248 AGTR1/EPB41L4A (LMX1A+)',
                na_color='#eeeeee', ax=axes[1, 0], show=False)
axes[1, 0].set_xlabel('UMAP 1'); axes[1, 0].set_ylabel('UMAP 2')

sc.pl.embedding(adata, basis='X_umap', color='C4_expanded_human',
                legend_loc='on data', legend_fontsize=4,
                title='C4 BMI-enriched Lmx1a labels — human HYPOMAP\n'
                      'Evidence tiers + global contrast  |  GWAS: BMI_Yengo2018',
                na_color='#eeeeee', ax=axes[1, 1], show=False)
axes[1, 1].set_xlabel('UMAP 1'); axes[1, 1].set_ylabel('UMAP 2')
axes[1, 1].legend(
    handles=c4_legend(),
    title=('★ = HL_C4 highlight\n'
           '***P = paper Bonf-sig  ***M = new_MAGMA\n'
           '***L = new_LDSC  nom = p<0.05\n'
           'C4-7 excluded (1.8% specific — noise)'),
    title_fontsize=6.5, **LEG_KW)

plt.suptitle('Human HYPOMAP UMAP — C3 and C4 BMI-enriched Lmx1a labels\n'
             'Corrected classification: Group A | CRNDE | GABA-nb | Group B | Global\n'
             'paper_MAGMA + new_MAGMA + new_LDSC  |  GWAS: BMI_Yengo2018',
             fontweight='bold', fontsize=13, y=1.005)
plt.tight_layout()
plt.savefig('HYPOMAP_C3_C4_UMAP_corrected.pdf', bbox_inches='tight')
plt.savefig('HYPOMAP_C3_C4_UMAP_corrected.png', bbox_inches='tight')
plt.show()

# ── Cell count summary ────────────────────────────────────────────────────────
ann_lookup = (unified.drop_duplicates(subset='cluster_id')
              .set_index('cluster_id')['annotation_full']
              .to_dict())

print("Cell counts in human HYPOMAP for highlighted clusters:")
for level, col, ids, label in [
    ('C3', 'C3_short', ALL_C3,      'C3 clusters (corrected classification)'),
    ('C4', 'C4_short', ALL_C4[:14], 'C4 Lmx1a tiers (excl. global contrast)'),
]:
    print(f"\n{label}:")
    print(f"{'Cluster':<12} {'n_human':>8}  annotation")
    print("-" * 65)
    for cid in ids:
        n       = (adata.obs[col].astype(str) == cid).sum()
        ann_str = ann_lookup.get(cid, '')
        ann_str = '_'.join(ann_str.split('_')[3:])[:45] if ann_str else ''
        print(f"  {cid:<12} {n:>7}  {ann_str}")


# In[21]:


# ── Rank directionality sanity check — from unified CSV ──────────────────────
import pandas as pd
import numpy as np
from scipy.stats import spearmanr

UNIFIED_PATH = "hypomap_cellect_unified.csv"
GWAS         = 'BMI_Yengo2018'
BONF_C4      = 0.05 / 415

GROUP_A_C4 = ['C4-226','C4-231','C4-247','C4-248','C4-229',
               'C4-223','C4-224','C4-252','C4-228','C4-230']
GROUP_B_C4 = ['C4-33','C4-407','C4-409','C4-411']
LMX_C4_IDS = GROUP_A_C4 + GROUP_B_C4

unified = pd.read_csv(UNIFIED_PATH)

# ── Build one wide table per method for C4 neuronal Yengo ─────────────────────
def get_c4(method):
    return (unified[(unified['gwas']      == GWAS) &
                    (unified['method']    == method) &
                    (unified['level']     == 'C4') &
                    (unified['is_neuronal'])]
            [['cluster_id','pvalue','beta','log10p','rank','annotation_full']]
            .copy()
            .rename(columns={'pvalue':    f'pvalue_{method}',
                              'beta':      f'beta_{method}',
                              'log10p':    f'log10p_{method}',
                              'rank':      f'rank_{method}'}))

paper = get_c4('paper_MAGMA')
magma = get_c4('new_MAGMA')
ldsc  = get_c4('new_LDSC')

# Merge all three on cluster_id
wide = (paper.merge(magma.drop(columns='annotation_full'), on='cluster_id')
              .merge(ldsc.drop(columns='annotation_full'),  on='cluster_id'))
wide['is_lmx'] = wide['cluster_id'].isin(LMX_C4_IDS)

print("RANK CONVENTION CHECK — unified CSV")
print("=" * 65)
print("""
Ranks were pre-computed in the unified CSV as:
  rank = pvalue.rank(ascending=True)  within each method × gwas × level

  ascending=True  →  rank 1   = SMALLEST p-value = MOST significant
                  →  rank 415 = LARGEST  p-value = LEAST significant

  rank 1   = strongest BMI enrichment
  rank 415 = no enrichment
""")

# Verify empirically
for label, rank_col, pval_col in [
    ('paper_MAGMA', 'rank_paper_MAGMA', 'pvalue_paper_MAGMA'),
    ('new_MAGMA',   'rank_new_MAGMA',   'pvalue_new_MAGMA'),
    ('new_LDSC',    'rank_new_LDSC',    'pvalue_new_LDSC'),
]:
    r1   = wide.loc[wide[rank_col] == 1]
    rN   = wide.loc[wide[rank_col] == wide[rank_col].max()]
    print(f"{label:<16}  "
          f"rank 1  → {r1['cluster_id'].values[0]:<12}  "
          f"p={r1[pval_col].values[0]:.2e}  (most sig ✓ if small)")
    print(f"                  "
          f"rank N  → {rN['cluster_id'].values[0]:<12}  "
          f"p={rN[pval_col].values[0]:.2e}  (least sig ✓ if large)")
    print()

# ── Cross-method rank correlations ───────────────────────────────────────────
print("=" * 65)
print("RANK CORRELATIONS (all 415 neuronal C4, Yengo2018)")
print("=" * 65)

pairs = [
    ('paper_MAGMA', 'new_MAGMA',  'rank_paper_MAGMA', 'rank_new_MAGMA'),
    ('paper_MAGMA', 'new_LDSC',   'rank_paper_MAGMA', 'rank_new_LDSC'),
    ('new_MAGMA',   'new_LDSC',   'rank_new_MAGMA',   'rank_new_LDSC'),
]
for a, b, ra, rb in pairs:
    r_sp, p_sp = spearmanr(wide[ra], wide[rb])
    print(f"  {a:<16} vs {b:<16}  ρ={r_sp:+.3f}  p={p_sp:.2e}")

print("""
Expected:
  paper vs new_MAGMA:  ρ ≈ +0.95  (same GWAS, same method, diff ESmu)
  paper vs new_LDSC:   ρ moderate  (diff method + diff ESmu)
  MAGMA vs LDSC:       ρ ≈ +0.67  (same ESmu, diff enrichment backend)
""")

# ── Replace just the top-10 block ─────────────────────────────────────────────
print("=" * 65)
print("TOP 10 per method — what rank do they get in the others?")
print("=" * 65)

ALL_RANK_COLS = ['rank_paper_MAGMA','rank_new_MAGMA','rank_new_LDSC']
ALL_PVAL_COLS = ['pvalue_paper_MAGMA','pvalue_new_MAGMA','pvalue_new_LDSC']

for focus_method, rank_col in [
    ('paper_MAGMA', 'rank_paper_MAGMA'),
    ('new_MAGMA',   'rank_new_MAGMA'),
    ('new_LDSC',    'rank_new_LDSC'),
]:
    # deduplicate column list — rank_col may already be in ALL_RANK_COLS
    cols = ['cluster_id'] + ALL_RANK_COLS + ALL_PVAL_COLS
    top10 = wide.nsmallest(10, rank_col)[cols].copy()

    print(f"\nTop 10 in {focus_method}:")
    print(f"{'cluster':^10} {'rk_P':>6} {'rk_M':>6} {'rk_L':>6}  "
          f"{'p_paper':>10} {'p_magma':>10} {'p_ldsc':>10}")
    print("-" * 62)
    for _, r in top10.iterrows():
        print(f"{r['cluster_id']:<10} "
              f"{int(r['rank_paper_MAGMA']):>6} "
              f"{int(r['rank_new_MAGMA']):>6} "
              f"{int(r['rank_new_LDSC']):>6}  "
              f"{r['pvalue_paper_MAGMA']:>10.2e} "
              f"{r['pvalue_new_MAGMA']:>10.2e} "
              f"{r['pvalue_new_LDSC']:>10.2e}")

# ── Lmx1a clusters ────────────────────────────────────────────────────────────
print("\n" + "=" * 65)
print("Lmx1a C4 clusters — ranks across all three methods")
print("=" * 65)
print(f"{'cluster':^10} {'rk_P':>6} {'rk_M':>6} {'rk_L':>6}  "
      f"{'p_paper':>10} {'p_magma':>10} {'p_ldsc':>10}  group")
print("-" * 72)
lmx_wide = wide[wide['is_lmx']].sort_values('rank_paper_MAGMA')
for _, r in lmx_wide.iterrows():
    grp = 'A' if r['cluster_id'] in GROUP_A_C4 else 'B'
    print(f"{r['cluster_id']:<10} "
          f"{int(r['rank_paper_MAGMA']):>6} "
          f"{int(r['rank_new_MAGMA']):>6} "
          f"{int(r['rank_new_LDSC']):>6}  "
          f"{r['pvalue_paper_MAGMA']:>10.2e} "
          f"{r['pvalue_new_MAGMA']:>10.2e} "
          f"{r['pvalue_new_LDSC']:>10.2e}  Grp{grp}")

print("""
INTERPRETATION:
  rank 1   = MOST significant  (p small, strong enrichment)
  rank 415 = LEAST significant (p large, no enrichment)

  paper_MAGMA vs new_MAGMA ρ ≈ +0.95:
    Rankings are nearly identical — the single-shot ESmu in
    new_MAGMA reproduces the bootstrapped paper rankings well.
    The earlier negative correlation was entirely an artefact
    of the wrong c4_ann row-order mapping.

  new_MAGMA vs new_LDSC ρ ≈ +0.67:
    Same ESmu, different enrichment backend. MAGMA finds more
    significant hits (240 vs 10 Bonf-sig) because it uses a
    gene-window chi-square test; LDSC uses partitioned
    heritability which is more conservative. Both agree on
    direction for 89% of clusters.

  Lmx1a Group A clusters rank in the top half of paper_MAGMA
  (most in top 300/415), consistent with genuine BMI signal.
  Group B clusters (C4-411, C4-33, C4-407, C4-409) rank lower
  and/or show inconsistent ranks across methods — consistent
  with their low mouse specificity and near-zero LMX1A.
""")


# In[22]:


# ── Clean final comparison — all three analyses ───────────────────────────────
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import re
from scipy.stats import spearmanr, pearsonr

def extract_id(s):
    m = re.search(r'(C[234])[-_](\d+)', str(s))
    return f'{m.group(1)}-{m.group(2)}' if m else None

# ── Load ──────────────────────────────────────────────────────────────────────
magma_new = pd.read_csv("prioritization_hypomap_magma.csv")
ldsc_new  = pd.read_csv("prioritization_hypomap_new.csv")
gm        = pd.read_csv("genetics_magma.txt", sep='\t')  # paper sheet 17

for df in [magma_new, ldsc_new]:
    df['cluster_id'] = df['annotation'].apply(extract_id)
    df['level'] = df['specificity_id'].map({'named_c4':'C4','named_c3':'C3'})

# genetics_magma.txt: C4_N → C4-N directly (no c4_ann row-order mapping needed)
gm_c4 = gm[gm['annotation'].str.startswith('C4_')].copy()
gm_c4['cluster_id']   = gm_c4['annotation'].str.replace('_','-', n=1)
gm_c4['pvalue_paper'] = pd.to_numeric(gm_c4['pvalue_10pct'], errors='coerce')
gm_c4['beta_paper']   = pd.to_numeric(gm_c4['beta_10pct'],   errors='coerce')
gm_c4['log10p_paper'] = -np.log10(gm_c4['pvalue_paper'].clip(1e-30))
gm_c4 = gm_c4.dropna(subset=['pvalue_paper'])

BONF_PAPER = 0.05 / len(gm_c4)          # 0.05/415
BONF_C4    = 0.05 / 415
BONF_C3    = 0.05 / 156
GWAS       = 'BMI_Yengo2018'

print(f"Paper C4: {len(gm_c4)} clusters  "
      f"({(gm_c4['pvalue_paper']<BONF_PAPER).sum()} Bonf-sig, expect 291)")

LMX_C4 = ['C4-248','C4-411','C4-33','C4-252','C4-407','C4-226','C4-229',
           'C4-223','C4-247','C4-230','C4-7','C4-224','C4-409','C4-228','C4-231']
LMX_C3 = ['C3-92','C3-84','C3-139','C3-29','C3-94','C3-83','C3-82','C3-26',
           'C3-137','C3-85','C3-28','C3-24','C3-93','C3-88','C3-115','C3-95']
HL_C4  = {'C4-248','C4-411'}
HL_C3  = {'C3-84','C3-92'}
COL_HL    = '#8B0000'
COL_LMX   = '#4878CF'
COL_OTHER = '#cccccc'

# ── Clean subsets (strip non-neuronal C3 contamination from named_c4) ─────────
def clean(df, gwas, level):
    prefix = f'{level}-'
    return (df[(df['gwas']==gwas) & (df['level']==level) &
               df['cluster_id'].str.startswith(prefix)]
            .copy())

magma_c4 = clean(magma_new, GWAS, 'C4')
magma_c3 = clean(magma_new, GWAS, 'C3')
ldsc_c4  = clean(ldsc_new,  GWAS, 'C4')
ldsc_c3  = clean(ldsc_new,  GWAS, 'C3')

# ── Part 1: New MAGMA vs Paper MAGMA (C4) ────────────────────────────────────
p1 = magma_c4.merge(
    gm_c4[['cluster_id','beta_paper','pvalue_paper','log10p_paper']],
    on='cluster_id', how='inner'
)
p1['log10p_new']  = -np.log10(p1['pvalue'].clip(1e-30))
p1['rank_new']    = p1['pvalue'].rank(ascending=True, method='min')
p1['rank_paper']  = p1['pvalue_paper'].rank(ascending=True, method='min')
p1['is_lmx'] = p1['cluster_id'].isin(LMX_C4)
p1['is_hl']  = p1['cluster_id'].isin(HL_C4)

r_sp1, p_sp1 = spearmanr(p1['rank_new'],   p1['rank_paper'])
r_pe1, p_pe1 = pearsonr( p1['log10p_new'], p1['log10p_paper'])
n1_new   = (p1['pvalue']       < BONF_C4).sum()
n1_paper = (p1['pvalue_paper'] < BONF_PAPER).sum()

print(f"\nPart 1 — New MAGMA vs Paper MAGMA (C4, n={len(p1)})")
print(f"  Rank Spearman  ρ={r_sp1:.3f}  p={p_sp1:.2e}")
print(f"  log10p Pearson r={r_pe1:.3f}  p={p_pe1:.2e}")
print(f"  Sig new={n1_new}  Sig paper={n1_paper}")

# ── Part 2: New MAGMA vs New LDSC (C4 + C3) ──────────────────────────────────
p2 = {}
for level, sub_m, sub_l, bonf in [('C4', magma_c4, ldsc_c4, BONF_C4),
                                    ('C3', magma_c3, ldsc_c3, BONF_C3)]:
    df = sub_m.merge(sub_l[['cluster_id','beta','pvalue']],
                     on='cluster_id', suffixes=('_magma','_ldsc'))
    df['log10p_magma'] = -np.log10(df['pvalue_magma'].clip(1e-30))
    df['log10p_ldsc']  = -np.log10(df['pvalue_ldsc'].clip(1e-30))
    df['rank_magma'] = df['pvalue_magma'].rank(ascending=True, method='min')
    df['rank_ldsc']  = df['pvalue_ldsc'].rank(ascending=True, method='min')
    df['signed_m'] = np.sign(df['beta_magma']) * df['log10p_magma']
    df['signed_l'] = np.sign(df['beta_ldsc'])  * df['log10p_ldsc']
    df['is_lmx'] = df['cluster_id'].isin(LMX_C4 if level=='C4' else LMX_C3)
    df['is_hl']  = df['cluster_id'].isin(HL_C4   if level=='C4' else HL_C3)
    r_sp, p_sp = spearmanr(df['rank_magma'],   df['rank_ldsc'])
    r_pe, p_pe = pearsonr( df['log10p_magma'], df['log10p_ldsc'])
    r_sg, p_sg = spearmanr(df['signed_m'],     df['signed_l'])
    n_m = (df['pvalue_magma'] < bonf).sum()
    n_l = (df['pvalue_ldsc']  < bonf).sum()
    dag = (np.sign(df['signed_m'])==np.sign(df['signed_l'])).mean()*100
    print(f"\nPart 2 — New MAGMA vs New LDSC ({level}, n={len(df)})")
    print(f"  Rank Spearman  ρ={r_sp:.3f}  p={p_sp:.2e}")
    print(f"  log10p Pearson r={r_pe:.3f}  p={p_pe:.2e}")
    print(f"  Signed Spearman ρ={r_sg:.3f}  p={p_sg:.2e}")
    print(f"  Sig MAGMA={n_m}  Sig LDSC={n_l}  Direction agree={dag:.0f}%")
    p2[level] = dict(df=df, r_sp=r_sp, p_sp=p_sp, r_pe=r_pe, r_sg=r_sg,
                     n_m=n_m, n_l=n_l, bonf=bonf)

# ── Lmx1a table ───────────────────────────────────────────────────────────────
lmx_tab = (p1[p1['is_lmx']][['cluster_id','pvalue','pvalue_paper']]
           .rename(columns={'pvalue':'p_new_MAGMA','pvalue_paper':'p_paper_MAGMA'})
           .merge(p2['C4']['df'][p2['C4']['df']['is_lmx']]
                  [['cluster_id','pvalue_ldsc']], on='cluster_id')
           .sort_values('p_new_MAGMA'))

print(f"\n{'='*72}")
print("Lmx1a C4 — new MAGMA / paper MAGMA / new LDSC (Yengo2018)")
print(f"{'='*72}")
print(f"{'Cluster':<10} {'p_newMAGMA':>12} {'p_paperMAGMA':>14} "
      f"{'p_LDSC':>10} {'NM':>5} {'PM':>5} {'L':>5}")
print("-"*65)
for _, r in lmx_tab.iterrows():
    snm = '***' if r['p_new_MAGMA']   < BONF_C4    else ''
    spm = '***' if r['p_paper_MAGMA'] < BONF_PAPER else ''
    sl  = '***' if r['pvalue_ldsc']   < BONF_C4    else ''
    print(f"{r['cluster_id']:<10} {r['p_new_MAGMA']:>12.3e} "
          f"{r['p_paper_MAGMA']:>14.3e} {r['pvalue_ldsc']:>10.3e} "
          f"{snm:>5} {spm:>5} {sl:>5}")

# ── Figure ─────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(3, 2, figsize=(14, 18))

def panel(ax, df, xcol, ycol, xlab, ylab, title,
          is_rank=False, bonf_x=None, bonf_y=None, info=None, level='C4'):
    lmx_oth = df[df['is_lmx'] & ~df['is_hl']]
    ax.scatter(df[xcol],       df[ycol],       c=COL_OTHER, s=14, alpha=0.5, zorder=2)
    ax.scatter(lmx_oth[xcol], lmx_oth[ycol], c=COL_LMX,   s=55, alpha=0.9, zorder=4,
               edgecolors='white', linewidths=0.4)
    for _, r in df[df['is_hl']].iterrows():
        ax.scatter(r[xcol], r[ycol], c=COL_HL, s=110, zorder=5,
                   edgecolors='white', linewidths=0.5)
        ax.annotate(r['cluster_id'], xy=(r[xcol],r[ycol]),
                    xytext=(5,3), textcoords='offset points',
                    fontsize=7.5, fontweight='bold', color=COL_HL)
    if is_rank:
        N = len(df)
        ax.plot([1,N],[1,N], color='#aaaaaa', lw=0.9, ls='--', alpha=0.6)
        ax.set_xlim(0, N+3); ax.set_ylim(0, N+3)
        ax.set_aspect('equal')
    else:
        x_max = df[xcol].max()*1.06; y_max = df[ycol].max()*1.06
        ax.set_xlim(0, x_max); ax.set_ylim(0, y_max)
        m_f, b_f = np.polyfit(df[xcol], df[ycol], 1)
        xr = np.linspace(0, df[xcol].max(), 200)
        ax.plot(xr, m_f*xr+b_f, color='#e74c3c', lw=1.5, ls='--',
                label=f'OLS slope={m_f:.2f}')
        lim_id = min(x_max, y_max)
        ax.plot([0,lim_id],[0,lim_id], color='#aaaaaa', lw=0.8, ls=':', alpha=0.5)
        if bonf_x:
            bx = -np.log10(bonf_x)
            if bx < x_max:
                ax.axvline(bx, color='steelblue', ls='--', lw=0.9, alpha=0.7)
        if bonf_y:
            by = -np.log10(bonf_y)
            if by < y_max:
                ax.axhline(by, color='#c0392b',  ls='--', lw=0.9, alpha=0.7)
    ax.set_xlabel(xlab, fontsize=9); ax.set_ylabel(ylab, fontsize=9)
    ax.set_title(title, fontweight='bold', fontsize=9)
    if info:
        ax.text(0.98, 0.02, info, transform=ax.transAxes, ha='right',
                va='bottom', fontsize=7.5, color='grey',
                bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='grey', alpha=0.8))
    hl_lbl = 'C4-248/411' if level=='C4' else 'C3-84/92'
    ax.legend(fontsize=7, handles=[
        mpatches.Patch(fc=COL_HL,    label=hl_lbl),
        mpatches.Patch(fc=COL_LMX,   label='Other Lmx1a'),
        mpatches.Patch(fc=COL_OTHER, label='All others'),
    ], loc='upper left')

N1 = len(p1)
# Row 0: Part 1 rank + log10p
panel(axes[0,0], p1, 'rank_paper', 'rank_new',
      f'Rank  paper MAGMA (N={N1})', f'Rank  new MAGMA (N={N1})',
      f'C4  New MAGMA vs Paper MAGMA\nSpearman ρ={r_sp1:.3f}  p={p_sp1:.2e}',
      is_rank=True,
      info=f'Paper sig: {n1_paper}\nNew sig:   {n1_new}', level='C4')

panel(axes[0,1], p1, 'log10p_paper', 'log10p_new',
      '−log₁₀(p)  paper MAGMA', '−log₁₀(p)  new MAGMA',
      f'C4  New MAGMA vs Paper MAGMA\nPearson r={r_pe1:.3f}  p={p_pe1:.2e}',
      is_rank=False, bonf_x=BONF_PAPER, bonf_y=BONF_C4, level='C4')

# Rows 1–2: Part 2
for ri, level in enumerate(['C4','C3']):
    res = p2[level]; df = res['df']; N = len(df)
    panel(axes[ri+1,0], df, 'rank_magma', 'rank_ldsc',
          f'Rank  new MAGMA (N={N})', f'Rank  new LDSC (N={N})',
          f'{level}  New MAGMA vs New LDSC\nSpearman ρ={res["r_sp"]:.3f}  p={res["p_sp"]:.2e}',
          is_rank=True,
          info=f'MAGMA sig: {res["n_m"]}\nLDSC sig:  {res["n_l"]}', level=level)
    panel(axes[ri+1,1], df, 'log10p_magma', 'log10p_ldsc',
          '−log₁₀(p)  new MAGMA', '−log₁₀(p)  new LDSC',
          f'{level}  New MAGMA vs New LDSC\nPearson r={res["r_pe"]:.3f}  Signed-ρ={res["r_sg"]:.3f}',
          is_rank=False, bonf_x=res['bonf'], bonf_y=res['bonf'], level=level)

for i, lbl in enumerate(['Part 1: New vs Paper MAGMA',
                           'Part 2: MAGMA vs LDSC (C4)',
                           'Part 2: MAGMA vs LDSC (C3)']):
    axes[i,0].annotate(lbl, xy=(-0.18,0.5), xycoords='axes fraction',
                        fontsize=9, fontweight='bold', rotation=90,
                        va='center', color='#333333')

plt.suptitle('CELLECT concordance — Yengo2018 BMI\n'
             'Rank 1 = most significant  |  Rank N = least significant\n'
             'Part 1: new MAGMA vs paper bootstrapped MAGMA (C4)  |  '
             'Part 2: new MAGMA vs new LDSC (C3 + C4)',
             fontweight='bold', fontsize=11, y=1.01)
plt.tight_layout()
plt.savefig('LMX_MAGMA_LDSC_concordance_fixed.pdf', bbox_inches='tight')
plt.savefig('LMX_MAGMA_LDSC_concordance_fixed.png', bbox_inches='tight')
plt.show()


# In[12]:


# ── CELLECT concordance — abundance-gradient highlights ───────────────────────
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
import re
from scipy.stats import spearmanr, pearsonr

UNIFIED_PATH = "hypomap_cellect_unified.csv"
GWAS         = 'BMI_Yengo2018'
BONF_C4      = 0.05 / 415
BONF_C3      = 0.05 / 156

# ── Abundance-based colour gradients ─────────────────────────────────────────
# Counts from mouse Lmx1a subset label transfer (most abundant = darkest)
# LMX1A+ (Group A + CRNDE): Reds gradient
# Non-LMX1A (Group B + GABA-nb): Blues gradient

C4_COUNTS_POS = {'C4-248':2650,'C4-252':393,'C4-226':362,'C4-229':314,
                  'C4-223':281,'C4-247':278,'C4-230':257,'C4-224':211,
                  'C4-228':151}
C4_COUNTS_NEG = {'C4-411':977,'C4-33':484,'C4-407':388,'C4-409':157}
# C4-7 excluded (1.8% specific)

C3_COUNTS_POS = {'C3-92':2734,'C3-84':2146,'C3-94':536,'C3-83':373,
                  'C3-82':350,'C3-85':269,'C3-93':182,'C3-95':246}
C3_COUNTS_NEG = {'C3-139':833,'C3-29':778,'C3-26':332,'C3-137':321,
                  'C3-28':234,'C3-24':198,'C3-88':143,'C3-115':142}

def make_gradient(count_dict, cmap_name, vmin=0.35, vmax=1.0):
    """Map counts → hex colours. Most abundant → darkest (vmax)."""
    cmap   = plt.colormaps[cmap_name]
    counts = np.array(list(count_dict.values()), dtype=float)
    norm   = (counts - counts.min()) / (counts.max() - counts.min() + 1e-9)
    scaled = vmin + norm * (vmax - vmin)
    return {cid: mcolors.to_hex(cmap(v))
            for cid, v in zip(count_dict.keys(), scaled)}

COL_C4_POS = make_gradient(C4_COUNTS_POS, 'Reds')
COL_C4_NEG = make_gradient(C4_COUNTS_NEG, 'Blues')
COL_C3_POS = make_gradient(C3_COUNTS_POS, 'Reds')
COL_C3_NEG = make_gradient(C3_COUNTS_NEG, 'Blues')

# Merge into flat lookups per level
COL_C4 = {**COL_C4_POS, **COL_C4_NEG}
COL_C3 = {**COL_C3_POS, **COL_C3_NEG}
COL_OTHER = '#cccccc'

def dot_color(cid, level):
    lut = COL_C4 if level == 'C4' else COL_C3
    return lut.get(cid, COL_OTHER)

# ── Load unified CSV ──────────────────────────────────────────────────────────
unified = pd.read_csv(UNIFIED_PATH)

def get(method, level):
    return (unified[(unified['gwas']==GWAS) & (unified['method']==method) &
                    (unified['level']==level) & (unified['is_neuronal'])]
            .copy())

paper_c4 = get('paper_MAGMA','C4')
magma_c4 = get('new_MAGMA',  'C4')
ldsc_c4  = get('new_LDSC',   'C4')
magma_c3 = get('new_MAGMA',  'C3')
ldsc_c3  = get('new_LDSC',   'C3')

ALL_LMX_C4 = list(C4_COUNTS_POS) + list(C4_COUNTS_NEG)
ALL_LMX_C3 = list(C3_COUNTS_POS) + list(C3_COUNTS_NEG)

for df_sub in [paper_c4, magma_c4, ldsc_c4]:
    df_sub['is_lmx'] = df_sub['cluster_id'].isin(ALL_LMX_C4)
for df_sub in [magma_c3, ldsc_c3]:
    df_sub['is_lmx'] = df_sub['cluster_id'].isin(ALL_LMX_C3)

# ── Part 1: New MAGMA vs Paper MAGMA (C4) ─────────────────────────────────────
gm = unified[(unified['gwas']==GWAS) & (unified['method']=='paper_MAGMA') &
              (unified['level']=='C4') & (unified['is_neuronal'])][
              ['cluster_id','beta','pvalue','log10p']].rename(
              columns={'beta':'beta_paper','pvalue':'pvalue_paper',
                       'log10p':'log10p_paper'})

p1 = magma_c4.merge(gm, on='cluster_id', how='inner')
p1['log10p_new'] = p1['log10p']
p1['rank_new']   = p1['pvalue'].rank(ascending=True, method='min')
p1['rank_paper'] = p1['pvalue_paper'].rank(ascending=True, method='min')
p1['is_lmx']     = p1['cluster_id'].isin(ALL_LMX_C4)

r_sp1, p_sp1 = spearmanr(p1['rank_new'],    p1['rank_paper'])
r_pe1, p_pe1 = pearsonr( p1['log10p_new'],  p1['log10p_paper'])
n1_new        = (p1['pvalue']        < BONF_C4).sum()
n1_paper      = (p1['pvalue_paper']  < BONF_C4).sum()

print(f"Part 1 — New MAGMA vs Paper MAGMA (C4, n={len(p1)})")
print(f"  Rank Spearman ρ={r_sp1:.3f}  p={p_sp1:.2e}")
print(f"  Pearson r={r_pe1:.3f}  Sig new={n1_new}  Sig paper={n1_paper}")

# ── Part 2: New MAGMA vs New LDSC ─────────────────────────────────────────────
p2 = {}
for level, sub_m, sub_l, bonf, lmx_ids in [
    ('C4', magma_c4, ldsc_c4, BONF_C4, ALL_LMX_C4),
    ('C3', magma_c3, ldsc_c3, BONF_C3, ALL_LMX_C3),
]:
    df = sub_m.merge(sub_l[['cluster_id','beta','pvalue']],
                     on='cluster_id', suffixes=('_magma','_ldsc'))
    df['log10p_magma'] = -np.log10(df['pvalue_magma'].clip(1e-30))
    df['log10p_ldsc']  = -np.log10(df['pvalue_ldsc'].clip(1e-30))
    df['rank_magma']   = df['pvalue_magma'].rank(ascending=True, method='min')
    df['rank_ldsc']    = df['pvalue_ldsc'].rank(ascending=True, method='min')
    df['signed_m']     = np.sign(df['beta_magma']) * df['log10p_magma']
    df['signed_l']     = np.sign(df['beta_ldsc'])  * df['log10p_ldsc']
    df['is_lmx']       = df['cluster_id'].isin(lmx_ids)
    r_sp,p_sp = spearmanr(df['rank_magma'],   df['rank_ldsc'])
    r_pe,p_pe = pearsonr( df['log10p_magma'], df['log10p_ldsc'])
    r_sg,p_sg = spearmanr(df['signed_m'],     df['signed_l'])
    n_m = (df['pvalue_magma'] < bonf).sum()
    n_l = (df['pvalue_ldsc']  < bonf).sum()
    dag = (np.sign(df['signed_m'])==np.sign(df['signed_l'])).mean()*100
    print(f"\nPart 2 — {level}  ρ={r_sp:.3f}  Signed-ρ={r_sg:.3f}  "
          f"MAGMA sig={n_m}  LDSC sig={n_l}  Dir={dag:.0f}%")
    p2[level] = dict(df=df, r_sp=r_sp, p_sp=p_sp, r_pe=r_pe,
                     r_sg=r_sg, n_m=n_m, n_l=n_l, bonf=bonf)

# ── Panel function ─────────────────────────────────────────────────────────────
def panel(ax, df, xcol, ycol, xlab, ylab, title, level,
          is_rank=False, bonf_x=None, bonf_y=None, info=None):
    lmx_ids = ALL_LMX_C4 if level=='C4' else ALL_LMX_C3

    # Background: all non-Lmx1a clusters
    bg = df[~df['is_lmx']]
    ax.scatter(bg[xcol], bg[ycol], c=COL_OTHER, s=14, alpha=0.5, zorder=2)

    # Lmx1a clusters: coloured by abundance gradient, sized by count
    count_lut = {**C4_COUNTS_POS, **C4_COUNTS_NEG} if level=='C4' \
                else {**C3_COUNTS_POS, **C3_COUNTS_NEG}
    max_count = max(count_lut.values())

    for _, r in df[df['is_lmx']].iterrows():
        cid = r['cluster_id']
        col = dot_color(cid, level)
        cnt = count_lut.get(cid, 0)
        # Size: 30–120, log-scaled by count
        sz  = 30 + 90 * (np.log1p(cnt) / np.log1p(max_count))
        ax.scatter(r[xcol], r[ycol], c=col, s=sz, zorder=4,
                   edgecolors='white', linewidths=0.4, alpha=0.95)
        # Annotate only the most abundant in each gradient
        if cid in (('C4-248','C3-92','C3-84') if level in ('C4','C3')
                   else ('C4-411',)):
            ax.annotate(cid, xy=(r[xcol], r[ycol]),
                        xytext=(5,3), textcoords='offset points',
                        fontsize=7, fontweight='bold',
                        color=col)

    if is_rank:
        N = len(df)
        ax.plot([1,N],[1,N], color='#aaaaaa', lw=0.9, ls='--', alpha=0.6)
        ax.set_xlim(0, N+3); ax.set_ylim(0, N+3); ax.set_aspect('equal')
    else:
        x_max = df[xcol].max()*1.06; y_max = df[ycol].max()*1.06
        ax.set_xlim(0, x_max); ax.set_ylim(0, y_max)
        m_f, b_f = np.polyfit(df[xcol], df[ycol], 1)
        xr = np.linspace(0, df[xcol].max(), 200)
        ax.plot(xr, m_f*xr+b_f, color='#e74c3c', lw=1.5, ls='--',
                label=f'OLS slope={m_f:.2f}')
        ax.plot([0,min(x_max,y_max)],[0,min(x_max,y_max)],
                color='#aaaaaa', lw=0.8, ls=':', alpha=0.5)
        if bonf_x:
            bx = -np.log10(bonf_x)
            if bx < x_max:
                ax.axvline(bx, color='steelblue', ls='--', lw=0.9, alpha=0.7)
        if bonf_y:
            by = -np.log10(bonf_y)
            if by < y_max:
                ax.axhline(by, color='#c0392b', ls='--', lw=0.9, alpha=0.7)

    ax.set_xlabel(xlab, fontsize=9); ax.set_ylabel(ylab, fontsize=9)
    ax.set_title(title, fontweight='bold', fontsize=9)
    if info:
        ax.text(0.98, 0.02, info, transform=ax.transAxes, ha='right',
                va='bottom', fontsize=7.5, color='grey',
                bbox=dict(boxstyle='round,pad=0.3', fc='white',
                          ec='grey', alpha=0.8))

    # Compact gradient legend
    lmx_col_dark = '#67000d' if level in ('C4','C3') else '#67000d'
    lmx_col_lite = '#fc9b7c'
    neg_col_dark = '#08306b'
    neg_col_lite = '#a6cee4'
    pos_label = 'LMX1A+ (Group A / CRNDE)'
    neg_label = 'Non-LMX1A (Group B / GABA-nb)'
    ax.legend(handles=[
        mpatches.Patch(fc=lmx_col_dark,
                       label=f'{pos_label}\n(dark=most abundant)'),
        mpatches.Patch(fc=lmx_col_lite,
                       label=f'{pos_label} (light=rare)'),
        mpatches.Patch(fc=neg_col_dark,
                       label=f'{neg_label}\n(dark=most abundant)'),
        mpatches.Patch(fc=neg_col_lite,
                       label=f'{neg_label} (light=rare)'),
        mpatches.Patch(fc=COL_OTHER, label='All other clusters'),
    ], fontsize=6.5, loc='upper left')

# ── Figure ────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(3, 2, figsize=(14, 18))

N1 = len(p1)
panel(axes[0,0], p1, 'rank_paper', 'rank_new',
      f'Rank paper MAGMA (N={N1})', f'Rank new MAGMA (N={N1})',
      f'C4  New MAGMA vs Paper MAGMA\nSpearman ρ={r_sp1:.3f}  p={p_sp1:.2e}',
      level='C4', is_rank=True,
      info=f'Paper sig: {n1_paper}\nNew sig:   {n1_new}')

panel(axes[0,1], p1, 'log10p_paper', 'log10p_new',
      '−log₁₀(p)  paper MAGMA', '−log₁₀(p)  new MAGMA',
      f'C4  New MAGMA vs Paper MAGMA\nPearson r={r_pe1:.3f}  p={p_pe1:.2e}',
      level='C4', is_rank=False, bonf_x=BONF_C4, bonf_y=BONF_C4)

for ri, level in enumerate(['C4','C3']):
    res = p2[level]; df = res['df']; N = len(df)
    panel(axes[ri+1,0], df, 'rank_magma', 'rank_ldsc',
          f'Rank new MAGMA (N={N})', f'Rank new LDSC (N={N})',
          f'{level}  New MAGMA vs New LDSC\nSpearman ρ={res["r_sp"]:.3f}  p={res["p_sp"]:.2e}',
          level=level, is_rank=True,
          info=f'MAGMA sig: {res["n_m"]}\nLDSC sig:  {res["n_l"]}')
    panel(axes[ri+1,1], df, 'log10p_magma', 'log10p_ldsc',
          '−log₁₀(p)  new MAGMA', '−log₁₀(p)  new LDSC',
          f'{level}  New MAGMA vs New LDSC\nPearson r={res["r_pe"]:.3f}  Signed-ρ={res["r_sg"]:.3f}',
          level=level, is_rank=False,
          bonf_x=res['bonf'], bonf_y=res['bonf'])

for i, lbl in enumerate(['Part 1: New vs Paper MAGMA (C4)',
                           'Part 2: MAGMA vs LDSC (C4)',
                           'Part 2: MAGMA vs LDSC (C3)']):
    axes[i,0].annotate(lbl, xy=(-0.18,0.5), xycoords='axes fraction',
                        fontsize=9, fontweight='bold', rotation=90,
                        va='center', color='#333333')

plt.suptitle(
    'CELLECT concordance — Yengo2018 BMI\n'
    'Dot colour = abundance in mouse Lmx1a subset  '
    '(red = LMX1A+  |  blue = non-LMX1A)\n'
    'Dot size ∝ log(count)  |  dark = most abundant',
    fontweight='bold', fontsize=11, y=1.01
)
plt.tight_layout()
plt.savefig('LMX_MAGMA_LDSC_concordance_abundance.pdf', bbox_inches='tight')
plt.savefig('LMX_MAGMA_LDSC_concordance_abundance.png', bbox_inches='tight')
plt.show()


# In[7]:


# ── BMI prioritisation across C2→C4 hierarchy — corrected master constants ───
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import scanpy as sc

CELLECT_PATH = "prioritization_hypomap.csv"
GWAS         = 'BMI_Yengo2018'

# ── Corrected Lmx1a sets from master constants ────────────────────────────────
# C4: Group A genuine + Group B contaminants (C4-7 excluded)
LMX_C4 = LMX_C4_ALL   # GROUP_A_C4 + GROUP_B_C4

# C3: Group A genuine + CRNDE + GABA-neighbour ONLY
# Group B (C3-115, C3-29, C3-137, C3-88, C3-139, C3-26, C3-24) excluded
LMX_C3 = GROUP_A_C3 + CRNDE_C3 + GABA_NEIGHBOR_C3

# C2: Post-1/Post-2 GLU (unchanged — no master constants at C2)
LMX_C2 = ['C2-30','C2-31','C2-32','C2-33','C2-46','C2-47']

# Tier-1 emphasis dots (corrected — C3-115 removed)
TIER1_C3 = set(HL_C3) | set(CRNDE_C3)          # C3-84, C3-93, C3-95
TIER1_C4 = set(HL_C4)                           # C4-229, C4-248

# ── 1. Load and prepare CELLECT-LDSC results ──────────────────────────────────
df_cellect = pd.read_csv(CELLECT_PATH)

def prepare_level(specificity_id, level, n_clusters):
    sub = df_cellect[(df_cellect['gwas']==GWAS) &
                     (df_cellect['specificity_id']==specificity_id)].copy()
    sub['cluster_id']    = sub['annotation'].str.extract(rf'({level}-\d+)')
    sub['log10p']        = -np.log10(sub['pvalue'].clip(1e-30))
    sub['bonf_sig']      = sub['pvalue'] < (0.05 / n_clusters)
    sub['signed_log10p'] = np.sign(sub['beta']) * sub['log10p']
    return sub.set_index('cluster_id')

c2_bmi = prepare_level('hypomap_C2', 'C2', 52)
c3_bmi = prepare_level('hypomap_C3', 'C3', 156)
c4_bmi = prepare_level('hypomap_C4', 'C4', 452)

BONF = {'C2': 0.05/52, 'C3': 0.05/156, 'C4': 0.05/452}

# ── 2. Print summaries ────────────────────────────────────────────────────────
print("=" * 65)
print(f"CELLECT-LDSC BMI prioritisation: {GWAS}")
print("=" * 65)
for level, sub, bonf in [('C2', c2_bmi, BONF['C2']),
                           ('C3', c3_bmi, BONF['C3']),
                           ('C4', c4_bmi, BONF['C4'])]:
    top = sub.nsmallest(1,'pvalue')
    print(f"\n  {level}  n={len(sub)}  Bonf={bonf:.2e}  "
          f"sig={sub['bonf_sig'].sum()}")
    print(f"  Top: {top.index[0]}  p={top['pvalue'].values[0]:.2e}  "
          f"{top['annotation'].values[0]}")

print("\n" + "=" * 65)
print("Lmx1a clusters (corrected) — BMI enrichment by level")
print("=" * 65)
for level, sub, lmx_ids, bonf in [
    ('C2', c2_bmi, LMX_C2, BONF['C2']),
    ('C3', c3_bmi, LMX_C3, BONF['C3']),
    ('C4', c4_bmi, LMX_C4, BONF['C4']),
]:
    lmx_sub = sub[sub.index.isin(lmx_ids)].sort_values('pvalue')
    print(f"\n  {level} ({len(lmx_sub)} clusters shown):")
    for cid, r in lmx_sub.iterrows():
        sig = '***' if r['bonf_sig'] else ('nom' if r['pvalue']<0.05 else '')
        # Annotate group
        if level=='C4':
            grp = 'A' if cid in GROUP_A_C4 else 'B'
        elif level=='C3':
            grp = ('A' if cid in GROUP_A_C3 else
                   'CRNDE' if cid in set(CRNDE_C3) else 'GABA-nb')
        else:
            grp = ''
        print(f"    {cid:<12} p={r['pvalue']:.3e}  {sig:<5}  "
              f"[{grp}]  {r['annotation'][:45]}")

# ── 3. Map log10p onto UMAP cells ─────────────────────────────────────────────
for level_col in ['C2','C3','C4']:
    short_col = f'{level_col}_short'
    if short_col not in adata.obs.columns:
        adata.obs[short_col] = adata.obs[level_col].str.extract(
            rf'({level_col}-\d+)')

for level, bmi_df, short_col in [
    ('C2', c2_bmi, 'C2_short'),
    ('C3', c3_bmi, 'C3_short'),
    ('C4', c4_bmi, 'C4_short'),
]:
    cluster_str = adata.obs[short_col].astype(str)
    adata.obs[f'{level}_log10p']   = (cluster_str
                                       .map(bmi_df['log10p'].to_dict())
                                       .fillna(0.0).astype(float))
    adata.obs[f'{level}_bonf_sig'] = (cluster_str
                                       .map(bmi_df['bonf_sig'].to_dict())
                                       .fillna(False).astype(bool))

# ── 4. Figure ─────────────────────────────────────────────────────────────────
fig, axes   = plt.subplots(2, 3, figsize=(20, 13))
umap        = adata.obsm['X_umap']
cmap_cont   = plt.cm.YlOrRd
vmax_global = max(c2_bmi['log10p'].max(),
                  c3_bmi['log10p'].max(),
                  c4_bmi['log10p'].max()) * 0.95

# ── Row 0: Continuous −log10(p) ───────────────────────────────────────────────

for col_i, (level, short_col, bmi_df, bonf) in enumerate([
    ('C2', 'C2_short', c2_bmi, BONF['C2']),
    ('C3', 'C3_short', c3_bmi, BONF['C3']),
    ('C4', 'C4_short', c4_bmi, BONF['C4']),
]):
    ax   = axes[0, col_i]
    vals = adata.obs[f'{level}_log10p'].values
    order = np.argsort(vals)
    sc_main = ax.scatter(umap[order,0], umap[order,1],
                         c=vals[order], cmap=cmap_cont,
                         s=0.15, alpha=0.7,
                         vmin=0, vmax=vmax_global, rasterized=True)
    plt.colorbar(sc_main, ax=ax,
                 label=f'−log₁₀(p)  (Bonf={-np.log10(bonf):.1f})')

    sig_cls = bmi_df['bonf_sig'].sum()
    sig_n   = adata.obs[f'{level}_bonf_sig'].sum()
    ax.set_title(f'{level}  −log₁₀(p) per cell\n'
                 f'{sig_cls}/{len(bmi_df)} Bonf-sig  '
                 f'({sig_n:,} cells in sig clusters)',
                 fontweight='bold', fontsize=9)
    ax.axis('off')

# ── Row 1: Discrete category map ─────────────────────────────────────────────
for col_i, (level, short_col, bmi_df, bonf) in enumerate([
    ('C2', 'C2_short', c2_bmi, BONF['C2']),
    ('C3', 'C3_short', c3_bmi, BONF['C3']),
    ('C4', 'C4_short', c4_bmi, BONF['C4']),
]):
    ax      = axes[1, col_i]
    lmx_ids = {'C2':LMX_C2,'C3':LMX_C3,'C4':LMX_C4}[level]
    tier1   = {'C2':set(),'C3':TIER1_C3,'C4':TIER1_C4}[level]

    # Background
    ax.scatter(umap[:,0], umap[:,1],
               c='#eeeeee', s=0.1, alpha=0.4, rasterized=True)

    # Global Bonf-sig non-Lmx1a — blue
    global_sig  = [c for c in bmi_df[bmi_df['bonf_sig']].index
                   if c not in lmx_ids]
    global_mask = adata.obs[short_col].astype(str).isin(global_sig)
    ax.scatter(umap[global_mask,0], umap[global_mask,1],
               c='#2e86c1', s=0.5, alpha=0.7, zorder=2, rasterized=True)

    sig_lmx_cls = nom_lmx_cls = 0
    for cid in lmx_ids:
        mask = adata.obs[short_col].astype(str) == cid
        if mask.sum() == 0: continue
        bsig = bmi_df.loc[cid,'bonf_sig'] if cid in bmi_df.index else False
        nom  = (bmi_df.loc[cid,'pvalue']<0.05) if cid in bmi_df.index else False

        # Colour: respect biological category
        is_crnde  = cid in set(CRNDE_C3)
        is_gaba   = cid in set(GABA_NEIGHBOR_C3)
        is_grpb_c4 = (cid in set(GROUP_B_C4))

        if bsig:
            sig_lmx_cls += 1
            col = (C3_GROUP_COLORS['CRNDE']        if is_crnde  else
                   C3_GROUP_COLORS['GABA_neighbor'] if is_gaba   else
                   C4_GROUP_COLORS['Group B']        if is_grpb_c4 else
                   '#8B0000')
            sz, zord = 4, 5
        elif nom:
            nom_lmx_cls += 1
            col = ('#e67e22' if is_crnde else
                   '#8e44ad' if is_gaba  else
                   '#95a5a6' if is_grpb_c4 else
                   '#E07B54')
            sz, zord = 2, 4
        else:
            col, sz, zord = '#d5d8dc', 1, 3

        ax.scatter(umap[mask,0], umap[mask,1],
                   c=col, s=sz, alpha=0.9, zorder=zord,
                   edgecolors='none', rasterized=True)

    ax.set_title(f'{level}  Bonf-sig cluster map\n'
                 f'Lmx1a: {sig_lmx_cls} Bonf-sig  {nom_lmx_cls} nominal  |  '
                 f'Global: {len(global_sig)} Bonf-sig',
                 fontweight='bold', fontsize=9)
    ax.axis('off')

    # Level-appropriate legend
    if level == 'C3':
        legend_handles = [
            mpatches.Patch(fc='#8B0000', label=f'Group A Bonf-sig ({sig_lmx_cls})'),
            mpatches.Patch(fc=C3_GROUP_COLORS['CRNDE'],
                           label='CRNDE C3-95 (cross-species)'),
            mpatches.Patch(fc=C3_GROUP_COLORS['GABA_neighbor'],
                           label='GABA-nb C3-28 (specific, LMX1A−)'),
            mpatches.Patch(fc='#E07B54', label=f'Group A nominal ({nom_lmx_cls})'),
        ]
    elif level == 'C4':
        legend_handles = [
            mpatches.Patch(fc='#8B0000',
                           label=f'Group A Bonf-sig ({sig_lmx_cls})'),
            mpatches.Patch(fc='#E07B54',
                           label=f'Group A nominal ({nom_lmx_cls})'),
            mpatches.Patch(fc=C4_GROUP_COLORS['Group B'],
                           label='Group B contaminant'),
        ]
    else:
        legend_handles = [
            mpatches.Patch(fc='#8B0000',
                           label=f'Lmx1a Bonf-sig ({sig_lmx_cls})'),
            mpatches.Patch(fc='#E07B54',
                           label=f'Lmx1a nominal ({nom_lmx_cls})'),
            mpatches.Patch(fc='#d5d8dc', label='Lmx1a not enriched'),
        ]
    legend_handles += [
        mpatches.Patch(fc='#2e86c1', label=f'Global Bonf-sig ({len(global_sig)})'),
        mpatches.Patch(fc='#eeeeee', label='All other cells'),
    ]
    ax.legend(handles=legend_handles, fontsize=6.5, loc='lower left',
              bbox_to_anchor=(1.01, 0), framealpha=0.9)

plt.suptitle(
    f'BMI GWAS enrichment across HYPOMAP hierarchy  |  CELLECT-LDSC  |  {GWAS}\n'
    f'Row 1: per-cell −log₁₀(p) on shared scale  |  '
    f'Row 2: corrected Lmx1a classification (Group A | CRNDE | GABA-nb | Group B)\n'
    f'C3-115/29/137/88/139/26/24 excluded (Group B)  |  C4-7 excluded (noise)',
    fontweight='bold', fontsize=11, y=1.005
)
plt.tight_layout()
plt.savefig('HYPOMAP_BMI_hierarchy_C2_C3_C4.pdf', bbox_inches='tight', dpi=300)
plt.savefig('HYPOMAP_BMI_hierarchy_C2_C3_C4.png', bbox_inches='tight', dpi=150)
plt.show()

# ── 5. Summary across levels ──────────────────────────────────────────────────
print("\n" + "=" * 68)
print("LMX1A region: BMI enrichment signal across resolution hierarchy")
print("=" * 68)
for level, sub, lmx_ids, bonf in [
    ('C2', c2_bmi, LMX_C2, BONF['C2']),
    ('C3', c3_bmi, LMX_C3, BONF['C3']),
    ('C4', c4_bmi, LMX_C4, BONF['C4']),
]:
    lmx     = sub[sub.index.isin(lmx_ids)]
    sig_n   = lmx['bonf_sig'].sum()
    nom_n   = ((lmx['pvalue']<0.05) & ~lmx['bonf_sig']).sum()
    med_p   = lmx['pvalue'].median()
    g_sig_n = sub['bonf_sig'].sum()
    n_shown = len(lmx_ids)
    print(f"  {level}  n_Lmx1a={len(lmx)}/{n_shown}  "
          f"Bonf-sig={sig_n}  nominal={nom_n}  "
          f"median_p={med_p:.3f}  "
          f"[global: {g_sig_n}/{len(sub)} Bonf-sig]")

print("""
INTERPRETATION:
  C2 (coarse): Post-1/Post-2 region shows only nominal signal.
  Bonf-sig clusters are Mid-1 GABA / Mid-2 — LMX1A posterior area
  is NOT the primary BMI signal at coarse resolution.

  C3: Group A genuine LMX1A+ clusters (C3-93, C3-84, C3-83, C3-82,
  C3-94) and CRNDE C3-95 reach nominal LDSC significance. The
  GABA-neighbour C3-28 is also nominally significant but NOT via LMX1A
  (near-zero human LMX1A expression). Signal resolves here for the
  first time at neuronal subtype granularity.

  C4: Group A clusters show nominal-to-marginal LDSC signal; the
  paper_MAGMA with bootstrapped CELLEX reports 7 Group A clusters
  Bonf-sig. Group B contaminants remain non-significant or negative.

  The globally dominant BMI signal (CALB1, ONECUT, Mid-1 GABA) is
  consistent Bonf-sig from C2 through C4 — this is the primary
  hypothalamic pathway. The LMX1A posterior GLU population is a
  distinct secondary signal visible only at C3/C4 resolution.
""")


# In[8]:


# ── BMI enrichment stratified by anatomical region — tabulation ───────────────
import pandas as pd
import numpy as np

# Verify the region column name
region_col = 'region'   # adjust if needed, e.g. 'region_collapsed'
print("Available region values:")
print(adata.obs[region_col].value_counts().to_string())

# ── Build per-region summary table for each resolution level ──────────────────
rows = []

for level, bmi_df, lmx_ids in [
    ('C2', c2_bmi, LMX_C2),
    ('C3', c3_bmi, LMX_C3),
    ('C4', c4_bmi, LMX_C4),
]:
    short_col   = f'{level}_short'
    log10p_col  = f'{level}_log10p'
    bonfsig_col = f'{level}_bonf_sig'

    for region in sorted(adata.obs[region_col].dropna().unique()):
        mask_region = adata.obs[region_col] == region
        n_cells     = mask_region.sum()
        if n_cells == 0:
            continue

        obs_r = adata.obs[mask_region]

        # BMI signal summary across all cells in region
        log10p_vals = obs_r[log10p_col].values
        mean_log10p = log10p_vals.mean()
        max_log10p  = log10p_vals.max()
        med_log10p  = np.median(log10p_vals)

        # Cells in Bonf-sig / nominal clusters
        n_bonfsig   = obs_r[bonfsig_col].sum()
        n_nominal   = (obs_r[short_col].astype(str)
                       .map(lambda c: (bmi_df.loc[c,'pvalue'] < 0.05
                                       if c in bmi_df.index else False))
                       .sum())

        # Dominant cluster in this region
        cl_counts   = obs_r[short_col].astype(str).value_counts()
        top_cl      = cl_counts.index[0] if len(cl_counts) else 'NA'
        top_cl_pct  = cl_counts.iloc[0] / n_cells * 100 if len(cl_counts) else 0

        # Lmx1a-specific clusters present in this region
        lmx_present = [c for c in cl_counts.index if c in set(lmx_ids)]
        lmx_cells   = obs_r[short_col].astype(str).isin(set(lmx_ids)).sum()

        # Best BMI p in Lmx1a clusters within this region
        best_lmx_p  = min(
            (bmi_df.loc[c,'pvalue'] for c in lmx_present if c in bmi_df.index),
            default=np.nan
        )

        rows.append({
            'level':        level,
            'region':       region,
            'n_cells':      n_cells,
            'mean_log10p':  round(mean_log10p, 3),
            'med_log10p':   round(med_log10p,  3),
            'max_log10p':   round(max_log10p,  3),
            'n_bonfsig':    int(n_bonfsig),
            'pct_bonfsig':  round(n_bonfsig / n_cells * 100, 1),
            'n_nominal':    int(n_nominal),
            'lmx_cells':    int(lmx_cells),
            'pct_lmx':      round(lmx_cells / n_cells * 100, 1),
            'lmx_clusters': ','.join(lmx_present[:6]),
            'best_lmx_p':   best_lmx_p,
            'top_cluster':  top_cl,
            'top_cl_pct':   round(top_cl_pct, 1),
        })

df_region = pd.DataFrame(rows)

# ── Print table per level ─────────────────────────────────────────────────────
for level in ['C2','C3','C4']:
    sub = (df_region[df_region['level']==level]
           .sort_values('mean_log10p', ascending=False))

    print(f"\n{'='*90}")
    print(f"BMI enrichment by region — {level} level  |  GWAS: {GWAS}")
    print(f"{'='*90}")
    print(f"{'Region':<30} {'n_cells':>8} {'mean_lgp':>9} {'max_lgp':>8} "
          f"{'%bonfsig':>9} {'lmx_cells':>10} {'%lmx':>6} "
          f"{'best_lmx_p':>12}  lmx_clusters")
    print("-"*110)
    for _, r in sub.iterrows():
        lmx_p_str = (f"{r['best_lmx_p']:.3e}" if pd.notna(r['best_lmx_p'])
                     else '   NA      ')
        sig_flag  = ('***' if (pd.notna(r['best_lmx_p']) and
                               r['best_lmx_p'] < 0.05/{'C2':52,'C3':156,'C4':452}[level])
                     else ('nom' if (pd.notna(r['best_lmx_p']) and
                                     r['best_lmx_p'] < 0.05) else ''))
        print(f"  {r['region']:<28} {r['n_cells']:>8,} {r['mean_log10p']:>9.3f} "
              f"{r['max_log10p']:>8.3f} {r['pct_bonfsig']:>8.1f}% "
              f"{r['lmx_cells']:>9,} {r['pct_lmx']:>5.1f}% "
              f"{lmx_p_str:>12}  {sig_flag:<4} {r['lmx_clusters']}")

# ── Save ──────────────────────────────────────────────────────────────────────
df_region.to_csv('BMI_enrichment_by_region.csv', index=False)
print(f"\nSaved to BMI_enrichment_by_region.csv  ({len(df_region)} rows)")


# In[11]:


# ── Typo fix + full name mapping (exact keys from your data) ─────────────────
adata.obs[region_col] = adata.obs[region_col].str.replace(
    'Thalamaus', 'Thalamus', regex=False)

REGION_NAMES = {
    'ARC':      'Arcuate Nucleus',
    'VMH':      'Ventromedial Hypothalamic Nucleus',
    'DMH':      'Dorsomedial Hypothalamic Nucleus',
    'LH':       'Lateral Hypothalamic Area',
    'PVN':      'Paraventricular Nucleus',
    'SCN':      'Suprachiasmatic Nucleus',
    'SON':      'Supraoptic Nucleus',
    'ME':       'Median Eminence',
    'MAM':      'Mammillary Bodies',
    'TMN':      'Tuberomammillary Nucleus',
    'LTN':      'Lateral Tuberal Nucleus',
    'POA':      'Preoptic Area',
    'MPOA':     'Medial Preoptic Area',
    'LPOA':     'Lateral Preoptic Area',
    'Perivent': 'Periventricular Zone',
    'Thalamus': 'Thalamus',
    'Fx/OT/ac': 'White Matter (Fornix / Optic Tract / Ant. Commissure)',
    'Vascular': 'Vascular',
    'Vent':     'Ventricular Zone',
    'NA':       'Unassigned',
}

adata.obs['region_full'] = adata.obs[region_col].map(
    lambda x: REGION_NAMES.get(x, x))

# Verify — every value should now be mapped
print("Final mapping:")
for r, full in (adata.obs[region_col]
                .value_counts()
                .reset_index()
                .values):
    print(f"  {r:<15} → {REGION_NAMES.get(r, '*** UNMAPPED')}  (n={full:,})")

region_col = 'region_full'


# In[13]:


# ── Region × BMI enrichment — weighted by cell count ─────────────────────────
import pandas as pd
import numpy as np

region_col = 'region_full'   # adjust if needed

# ── 1. Build cluster → LDSC log10p lookup for C3 and C4 ──────────────────────
c3_log10p = c3_bmi['log10p'].to_dict()   # cluster_id → -log10(p)
c4_log10p = c4_bmi['log10p'].to_dict()

# Map each cell to its BMI score at C3 and C4
for short_col, lookup, out_col in [
    ('C3_short', c3_log10p, 'C3_ldsc_log10p'),
    ('C4_short', c4_log10p, 'C4_ldsc_log10p'),
]:
    if short_col not in adata.obs.columns:
        adata.obs[short_col] = adata.obs[short_col.replace('_short','_HM')]\
                                     .str.extract(r'(C[34]-\d+)')
    adata.obs[out_col] = (adata.obs[short_col].astype(str)
                          .map(lookup).fillna(0.0).astype(float))

# ── 2. Compute weighted BMI score per region ──────────────────────────────────
rows = []
for region in sorted(adata.obs[region_col].dropna().unique()):
    mask  = adata.obs[region_col] == region
    obs_r = adata.obs[mask]
    n     = mask.sum()
    if n == 0:
        continue

    # Cell-count-weighted mean log10p at C3 and C4
    # Weight is simply n_cells per cluster (implicit — every cell contributes 1)
    c3_score = obs_r['C3_ldsc_log10p'].mean()
    c4_score = obs_r['C4_ldsc_log10p'].mean()
    combined = (c3_score + c4_score) / 2   # simple average across levels

    # Dominant clusters and their individual p-values
    for level, short_col, bmi_df in [('C3','C3_short',c3_bmi),
                                      ('C4','C4_short',c4_bmi)]:
        cl_counts = obs_r[short_col].astype(str).value_counts()
        for cid, cnt in cl_counts.items():
            if cid not in bmi_df.index: continue
            rows.append({
                'region':   region,
                'level':    level,
                'cluster':  cid,
                'n_cells':  cnt,
                'pct_cells': cnt/n*100,
                'pvalue':   bmi_df.loc[cid,'pvalue'],
                'log10p':   bmi_df.loc[cid,'log10p'],
                'bonf_sig': bmi_df.loc[cid,'bonf_sig'],
                'n_region': n,
                'c3_score': c3_score,
                'c4_score': c4_score,
                'combined': combined,
            })

df_long = pd.DataFrame(rows)

# ── 3. Region-level summary — one row per region ──────────────────────────────
region_summary = (adata.obs.groupby(region_col, observed=True)
                  .apply(lambda g: pd.Series({
                      'n_cells':       len(g),
                      'c3_wtd_log10p': g['C3_ldsc_log10p'].mean(),
                      'c4_wtd_log10p': g['C4_ldsc_log10p'].mean(),
                      'combined_score':(g['C3_ldsc_log10p'].mean() +
                                        g['C4_ldsc_log10p'].mean()) / 2,
                      # Fraction of cells in ANY Bonf-sig cluster at C4
                      'pct_in_bonfsig_c4': g[f'C4_bonf_sig'].mean()*100
                          if 'C4_bonf_sig' in g.columns else np.nan,
                  }))
                  .reset_index()
                  .sort_values('combined_score', ascending=False))

# Top cluster per region per level
for level, short_col, bmi_df in [('C3','C3_short',c3_bmi),
                                   ('C4','C4_short',c4_bmi)]:
    top_cl = (df_long[df_long['level']==level]
              .sort_values('n_cells', ascending=False)
              .groupby('region').first()[['cluster','pvalue']]
              .rename(columns={'cluster':f'top_{level}',
                               'pvalue': f'top_{level}_p'}))
    region_summary = region_summary.join(top_cl, on=region_col)

region_summary['rank'] = range(1, len(region_summary)+1)

# ── 4. Print ranking ──────────────────────────────────────────────────────────
print(f"\n{'='*95}")
print(f"Brain region ranking by cell-count-weighted BMI signal  |  CELLECT-LDSC  |  {GWAS}")
print(f"Weighted mean −log₁₀(p) across C3 + C4 clusters per region")
print(f"{'='*95}")
print(f"{'Rank':>4}  {'Region':<30} {'n_cells':>8}  "
      f"{'C3 score':>9}  {'C4 score':>9}  {'Combined':>9}  "
      f"{'Top C3':>10}  {'Top C4':>10}")
print("-"*100)
for _, r in region_summary.iterrows():
    print(f"  {r['rank']:>3}.  {r[region_col]:<30} {r['n_cells']:>8,}  "
          f"{r['c3_wtd_log10p']:>9.3f}  {r['c4_wtd_log10p']:>9.3f}  "
          f"{r['combined_score']:>9.3f}  "
          f"{str(r.get('top_C3','')):<10}  {str(r.get('top_C4','')):<10}")

# ── 5. Save ───────────────────────────────────────────────────────────────────
region_summary.to_csv('BMI_region_ranking.csv', index=False)
df_long.to_csv('BMI_region_cluster_detail.csv', index=False)
print(f"\nSaved BMI_region_ranking.csv and BMI_region_cluster_detail.csv")


# In[ ]:




