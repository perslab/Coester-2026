#!/usr/bin/env python
# coding: utf-8

# In[19]:


import os as os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.stats.multitest as smm
import statsmodels.api as sm
import statsmodels.formula.api as smf
import adjustText
from scipy.stats import pearsonr, spearmanr


# In[3]:


adata = sc.read("allgenes_bipca_annotated.h5ad")
adata


# In[5]:


adata.file.close()


# In[3]:


sc.pl.embedding(adata,basis="tsne",color="leiden", legend_loc="on data", legend_fontsize="xx-small")
sc.pl.embedding(adata,basis="tsne",color="manual_celltype_annotation", legend_loc="on data", legend_fontsize="xx-small")
sc.pl.embedding(adata,basis="tsne",color="detailed_celltype_annotation", legend_loc="on data", legend_fontsize="xx-small")


# In[6]:


marker_genes = ["Lepr"]
sc.pl.embedding(adata,basis="tsne",color=marker_genes,layer="scran_normalization", color_map="pink_r")


# In[ ]:


##This saves the EMBEDDINGS for import to Seurat
pcs = adata.obsm["pcs"]
df = pd.DataFrame(pcs, index=adata.obs_names)
df.to_csv("allgenes_bipca_annotated_PCS.csv")
tsne = adata.obsm["tsne"]
df = pd.DataFrame(tsne, index=adata.obs_names)
df.to_csv("allgenes_bipca_annotated_TSNE.csv")


# In[6]:


adata.varm["V_biwhite"]


# In[9]:


##This saves the LOADINGS for import to Seurat
pca_loadings = adata.varm["V_biwhite"]  # Typical Scanpy storage
# Convert to DataFrame (genes x PCs)
loadings_df = pd.DataFrame(
    pca_loadings,
    index=adata.var_names,  # Gene names
    columns=[f"PC_{i+1}" for i in range(pca_loadings.shape[1])]
)
loadings_df


# In[10]:


# Save to CSV
loadings_df.to_csv("allgenes_bipca_annotated_PCS_loadings.csv")


# In[4]:


print(adata.uns["manual_celltype_annotation_colors"])


# In[5]:


print(adata.obs["manual_celltype_annotation"].cat.categories)


# In[6]:


print(adata.uns["detailed_celltype_annotation_colors"])


# In[7]:


print(adata.obs["detailed_celltype_annotation"].cat.categories)


# In[3]:


adata.obs["sort"] = pd.to_numeric(adata.obs["sort"], errors="coerce")


# In[4]:


sc.pl.embedding(adata,basis="tsne",color="sort", color_map="pink_r")


# In[5]:


marker_genes = ["v147"]
sc.pl.embedding(adata,basis="tsne",color=marker_genes,layer="scran_normalization", color_map="pink_r", vmax=1)


# In[11]:


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

sns.regplot(data=plot_data, x='Gene Expression', y='Sort Values', scatter_kws={'alpha':0.5})
plt.title(f'{gene_name} Expression vs Sort Values\nPearson: {pearson_corr:.2f}, Spearman: {spearman_corr:.2f}')
plt.xlabel(f'{gene_name} Expression')
plt.ylabel('Sort Values')
plt.show()


# In[14]:


data = pd.DataFrame({
    'Gene_Expression': gene_expression,
    'Sort_Values': sort_values,
    'Animal': adata.obs['hash.mcl.ID']
})

# Fit a linear mixed-effects model
model = smf.mixedlm("Gene_Expression ~ Sort_Values", data, groups=data["Animal"])
result = model.fit()

# Print the results
print(result.summary())


# In[28]:


adata_sort1 = adata[adata.obs['sort'] == 1]
adata_sort1


# In[32]:


gene_name = "v147"
gene_expression = adata_sort1[:, gene_name].layers["scran_normalization"].toarray().flatten()
plt.hist(gene_expression, bins=50, color='skyblue', edgecolor='black')
plt.ylim(0,35000)
plt.title('Expression of v147')
plt.xlabel('Expression Level')
plt.ylabel('Number of Cells')
plt.show()


# In[36]:


adata_sort2 = adata[adata.obs['sort'] < 1]
gene_name = "v147"
gene_expression = adata_sort2[:, gene_name].layers["scran_normalization"].toarray().flatten()
plt.hist(gene_expression, bins=50, color='skyblue', edgecolor='black')
plt.ylim(0,35000)
plt.title('Expression of v147')
plt.xlabel('Expression Level')
plt.ylabel('Number of Cells')
plt.show()


# In[37]:


gene_name = "v147"
gene_series = sc.get.obs_df(adata, keys=gene_name)
adata_subset = adata[gene_series[gene_name] > 0, :]
gene_expression = adata_subset[:, gene_name].layers["scran_normalization"].toarray().flatten()
plt.hist(gene_expression, bins=50, color='skyblue', edgecolor='black')
plt.title('Expression of v147')
plt.xlabel('Expression Level')
plt.ylabel('Number of Cells')
plt.show()


# In[46]:


# 1. Create a DataFrame with expression and condition
df = pd.DataFrame({
    'Expression': adata_subset[:, gene_name].layers["scran_normalization"].toarray().flatten(), 
    'Condition': adata_subset.obs['sort']  # The column to color by
})

sns.jointplot(
    data=df,
    x='Condition',
    y='Expression',
    kind='scatter',
    alpha=0.3,
    marginal_kws={'bins': 30}
)

plt.suptitle('Expression vs. Continuous Condition', y=1.02)
plt.show()


# In[54]:


# Create a jointplot
g = sns.jointplot(
    data=df,
    x='Condition',
    y='Expression',
    kind='scatter',
    alpha=0.5,
    marginal_kws={'bins': 30}
)

from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
from sklearn.pipeline import make_pipeline
# Fit polynomial regression (degree=2 for quadratic)
X = df['Condition'].values.reshape(-1, 1)
y = df['Expression'].values
model = make_pipeline(PolynomialFeatures(degree=2), LinearRegression())
model.fit(X, y)

# Predict and plot
x_range = np.linspace(df['Condition'].min(), df['Condition'].max(), 100)
y_pred = model.predict(x_range.reshape(-1, 1))
g.ax_joint.plot(x_range, y_pred, color='red', label='Quadratic Fit', lw=2)


# In[5]:


gene_name = "v147"
gene_expression = adata[:, gene_name].X.toarray().flatten() 
adata.obs["AAV"] = np.where(gene_expression > 0, True, False)
sc.pl.embedding(adata,basis="tsne",color="AAV")


# In[60]:


from scipy.stats import pointbiserialr

# Extract data
binary = adata.obs["AAV"].astype(int)  # Convert to 0/1
continuous = adata.obs["sort"]

# Calculate correlation
corr, p_value = pointbiserialr(binary, continuous)
print(f"Point-Biserial Correlation: {corr:.3f}, p-value: {p_value:.3g}")


# In[65]:


sns.histplot(
    data=adata.obs,
    x="sort",
    hue="AAV",
    kde=True,
    alpha=0.5,
    palette={True: "red", False: "gray"}
)
plt.title("Distribution by Binary Group")
plt.show()


# In[6]:


# Calculate the ratio of 'True' AAV per celltype
aav_ratio_per_celltype = (
    adata.obs.groupby("detailed_celltype_annotation")["AAV"]  # Group by celltype
    .mean()                              # Mean of True (1) and False (0) gives ratio
    .rename("AAV_ratio")                 # Rename the series
)

# Map the ratio back to the original obs DataFrame
adata.obs["AAV_ratio"] = adata.obs["detailed_celltype_annotation"].map(aav_ratio_per_celltype)
adata.obs["AAV_ratio"] = pd.to_numeric(adata.obs["AAV_ratio"], errors="coerce")


# In[7]:


sc.pl.embedding(adata,basis="tsne",color="AAV_ratio", color_map="inferno")


# In[5]:


# Calculate AAV ratio per celltype AND per area
aav_ratio_per_celltype_area = (
    adata.obs.groupby(["detailed_celltype_annotation", "area"])["AAV"]  # Group by both
    .mean()                              # Mean of True (1) and False (0) gives ratio
    .reset_index()                       # Convert to DataFrame
)
aav_ratio_per_celltype_area


# In[8]:


areas = ["BST", "CEA", "LHA", "PAG", "PBN", "PVH", "PVT"]
for area in areas:
    # Filter ratios for the current area
    area_ratios = (
        aav_ratio_per_celltype_area
        .loc[aav_ratio_per_celltype_area["area"] == area]
        .set_index("detailed_celltype_annotation")["AAV"]
        .rename(f"AAV_ratio_{area}")
    )
    # Map to adata.obs
    adata.obs[f"AAV_ratio_{area}"] = adata.obs["detailed_celltype_annotation"].map(area_ratios)
adata.obs


# In[9]:


sc.pl.embedding(adata,basis="tsne",color="AAV_ratio_BST", color_map="inferno")
sc.pl.embedding(adata,basis="tsne",color="AAV_ratio_CEA", color_map="inferno")
sc.pl.embedding(adata,basis="tsne",color="AAV_ratio_LHA", color_map="inferno")
sc.pl.embedding(adata,basis="tsne",color="AAV_ratio_PAG", color_map="inferno")
sc.pl.embedding(adata,basis="tsne",color="AAV_ratio_PBN", color_map="inferno")
sc.pl.embedding(adata,basis="tsne",color="AAV_ratio_PVH", color_map="inferno")
sc.pl.embedding(adata,basis="tsne",color="AAV_ratio_PVT", color_map="inferno")


# In[10]:


adata.obs.columns


# In[11]:


sc.pl.embedding(adata,basis="tsne",color="predicted_cell_type")


# In[12]:


# Group by predicted_cell_type + area, then compute mean AAV ratio
aav_ratio_per_group = (
    adata.obs.groupby(["predicted_cell_type", "area"])["AAV"]
    .mean()
    .unstack(level="predicted_cell_type")  # Pivot to get Glia/Neuron as columns
    .rename(columns={"Glia": "glia_ratio", "Neuron": "neuron_ratio"})
)
aav_ratio_per_group


# In[13]:


adata.obs['AAV_ratio_BST_norm'] = (adata.obs['AAV_ratio_BST'])*0.021718
adata.obs['AAV_ratio_CEA_norm'] = (adata.obs['AAV_ratio_CEA'])*0.103428
adata.obs['AAV_ratio_LHA_norm'] = (adata.obs['AAV_ratio_LHA'])*0.116275
adata.obs['AAV_ratio_PAG_norm'] = (adata.obs['AAV_ratio_PAG'])*0.018910
adata.obs['AAV_ratio_PBN_norm'] = (adata.obs['AAV_ratio_PBN'])*0.021776
adata.obs['AAV_ratio_PVH_norm'] = (adata.obs['AAV_ratio_PVH'])*0.148016
adata.obs['AAV_ratio_PVT_norm'] = (adata.obs['AAV_ratio_PVT'])*0.021016


# In[14]:


sc.pl.embedding(adata,basis="tsne",color="AAV_ratio_BST_norm", color_map="inferno")
sc.pl.embedding(adata,basis="tsne",color="AAV_ratio_CEA_norm", color_map="inferno")
sc.pl.embedding(adata,basis="tsne",color="AAV_ratio_LHA_norm", color_map="inferno")
sc.pl.embedding(adata,basis="tsne",color="AAV_ratio_PAG_norm", color_map="inferno")
sc.pl.embedding(adata,basis="tsne",color="AAV_ratio_PBN_norm", color_map="inferno")
sc.pl.embedding(adata,basis="tsne",color="AAV_ratio_PVH_norm", color_map="inferno")
sc.pl.embedding(adata,basis="tsne",color="AAV_ratio_PVT_norm", color_map="inferno")


# In[9]:


metadata = adata.obs
metadata.to_csv("allgenes_bipca_annotated_obs.csv")
adata


# In[3]:


#adata.write("allgenes_bipca_annotated.h5ad")
#adata.write("allgenes_bipca_annotated_APR02.h5ad")
adata = sc.read("allgenes_bipca_annotated.h5ad", cache=True)
adata


# In[4]:


sc.pl.embedding(adata,basis="tsne",color="sex")


# In[5]:


sc.pl.embedding(adata,basis="tsne",color="area")


# In[15]:


sc.pl.embedding(adata,basis="tsne",color="detailed_celltype_annotation", legend_loc="on data", legend_fontsize="xx-small")


# In[46]:


ax = sc.pl.embedding(
    adata,
    basis="tsne",
    color="detailed_celltype_annotation",
    legend_loc="on data",
    legend_fontsize="xx-small",
    return_fig=True
)
# Obtain all text objects from the Axes (these are the on-data labels)
ax=plt.gca()
texts = ax.texts
adjustText.adjust_text(texts, 
                       ax=ax,
                       arrowprops=dict(arrowstyle="-", color='black', lw=0.5)
                      )

#plt.show()
plt.savefig("detailed_repel.pdf")


# In[26]:


ax = sc.pl.embedding(
    adata,
    basis="tsne",
    color="detailed_celltype_annotation",
    legend_loc="on data",
    legend_fontsize="xx-small",
    return_fig=True
)


# In[51]:


adata.obs.columns


# In[52]:


freq_table = adata.obs["detailed_celltype_annotation"].value_counts(dropna=False).reset_index()
freq_table.columns = ["detailed_celltype_annotation", "frequency"]
freq_table.to_csv("frequency_table.csv", index=False)


# In[10]:


sc.pl.embedding(adata,basis="tsne",color="cellect_bmi_detailed_BH",layer="Z_biwhite", color_map="Spectral")


# In[8]:


adata.obs["cellect_bmi_log10"] = -np.log10(adata.obs["cellect_bmi_detailed"])


# In[11]:


sc.pl.embedding(adata,basis="tsne",color="cellect_bmi_log10",layer="Z_biwhite", color_map="Spectral_r")


# In[14]:


metadata_df = pd.read_csv("7plex_region_labeltransfer_meta.csv", index_col=0)
metadata_df


# In[15]:


metadata_df = metadata_df[["region_collapsed", "score_region_collapsed"]]
adata.obs = adata.obs.join(metadata_df, how='left')
adata.obs


# In[17]:


colors_regions = [
"#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
"#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#4FC601", "#8FB0FF",
"#997D87", "#5A0007", "#809693", "#FF1A59", "#1B4400", "#BF5650", "#3B5DFF",
"#4A3B53"
]
adata.uns["region_collapsed_colors"]=colors_regions
sc.pl.embedding(adata,basis="tsne",
                color="region_collapsed",
               save="_7plex_region_unfiltered.pdf")


# In[17]:


freq_table = adata.obs["region_collapsed"].value_counts(dropna=False).reset_index()
freq_table


# In[4]:


# Group by 'celltype' and compute the most common 'region' for each cell type.
# If there's a tie, mode() returns all modes; here we take the first one.
predicted_map = adata.obs.groupby("detailed_celltype_annotation")["region_collapsed"].agg(lambda x: x.mode().iloc[0])
# Map the predicted region back to the full dataframe based on the cell type.
adata.obs["predicted_spatial_region"] = adata.obs["detailed_celltype_annotation"].map(predicted_map)


# In[16]:


adata.uns["predicted_spatial_region_colors"]=colors_regions
ax = sc.pl.embedding(
    adata,
    basis="tsne",
    color="predicted_spatial_region",
    save="_7plex_region.pdf"
)


# In[9]:


#Compute the total number of cells in each region.
region_total = adata.obs.groupby("region_collapsed").size().rename("total_cells")
region_total_df = region_total.reset_index()  # convert to DataFrame
region_total_df


# In[11]:


#Compute counts for each celltype-region combination.
counts = (
    adata.obs
    .groupby(["detailed_celltype_annotation", "region_collapsed"])
    .size()
    .rename("celltype_count")
    .reset_index()
)
#Merge the total count per region into the celltype-region counts.
counts = counts.merge(region_total_df, on="region_collapsed", how="left")
#Compute normalized frequency:
counts["normalized"] = counts["celltype_count"] / counts["total_cells"]
counts


# In[12]:


# For each cell type, select the region with the highest normalized value.
# The following finds the index of the maximum normalized value per celltype.
idx = counts.groupby("detailed_celltype_annotation")["normalized"].idxmax()
max_counts = counts.loc[idx]

# Create a mapping from each cell type to its predicted region.
predicted_map = dict(zip(max_counts["detailed_celltype_annotation"], max_counts["region_collapsed"]))

# Step 6. Apply this mapping to adata.obs to create the new column.
adata.obs["predicted_spatial_region_norm"] = adata.obs["detailed_celltype_annotation"].map(predicted_map)

# Optionally, view the counts and mapping:
print(max_counts)
print(predicted_map)


# In[26]:


colors_regions_n = [
"#FFFF00", "#1CE6FF", "#FF4A46", "#008941", "#A30059",
"#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#4FC601", "#8FB0FF",
"#997D87", "#809693", "#FF1A59", "#1B4400", "#BF5650", "#3B5DFF",
"#4A3B53"
]
adata.uns["predicted_spatial_region_norm_colors"]=colors_regions_n
sc.pl.embedding(
    adata,
    basis="tsne",
    color="predicted_spatial_region_norm",
    save="_7plex_region_norm.pdf"
)


# In[18]:


adata.obs["region_collapsed_filtered"] = adata.obs["region_collapsed"].copy()
adata.obs["region_collapsed_filtered"] = adata.obs["region_collapsed_filtered"].where(adata.obs["score_region_collapsed"] > 0.4, np.nan)
adata.uns["region_collapsed_filtered_colors"] = colors_regions
sc.pl.embedding(adata,basis="tsne",
                color="region_collapsed_filtered")


# In[19]:


predicted_map = adata.obs.groupby("detailed_celltype_annotation")["region_collapsed_filtered"].agg(lambda x: x.mode().iloc[0])
# Map the predicted region back to the full dataframe based on the cell type.
adata.obs["predicted_spatial_region_filtered"] = adata.obs["detailed_celltype_annotation"].map(predicted_map)
adata.uns["predicted_spatial_region_filtered_colors"]=colors_regions
ax = sc.pl.embedding(
    adata,
    basis="tsne",
    color="predicted_spatial_region_filtered",
    save="_7plex_region.pdf"
)


# In[20]:


#Compute the total number of cells in each region.
region_total = adata.obs.groupby("region_collapsed_filtered").size().rename("total_cells")
region_total_df = region_total.reset_index()  # convert to DataFrame
#Compute counts for each celltype-region combination.
counts = (
    adata.obs
    .groupby(["detailed_celltype_annotation", "region_collapsed_filtered"])
    .size()
    .rename("celltype_count")
    .reset_index()
)
#Merge the total count per region into the celltype-region counts.
counts = counts.merge(region_total_df, on="region_collapsed_filtered", how="left")
#Compute normalized frequency:
counts["normalized"] = counts["celltype_count"] / counts["total_cells"]
# For each cell type, select the region with the highest normalized value.
# The following finds the index of the maximum normalized value per celltype.
idx = counts.groupby("detailed_celltype_annotation")["normalized"].idxmax()
max_counts = counts.loc[idx]

# Create a mapping from each cell type to its predicted region.
predicted_map = dict(zip(max_counts["detailed_celltype_annotation"], max_counts["region_collapsed_filtered"]))

# Step 6. Apply this mapping to adata.obs to create the new column.
adata.obs["predicted_spatial_region_filtered_norm"] = adata.obs["detailed_celltype_annotation"].map(predicted_map)


# In[22]:


print(max_counts)
print(predicted_map)


# In[23]:


colors_regions_f = [
"#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", 
"#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#4FC601", "#8FB0FF",
"#997D87", "#5A0007", "#FF1A59", "#1B4400", "#BF5650", "#3B5DFF",
"#4A3B53"
]
adata.uns["predicted_spatial_region_filtered_norm_colors"]=colors_regions_f
sc.pl.embedding(
    adata,
    basis="tsne",
    color="predicted_spatial_region_filtered_norm",
    save="_7plex_region_norm_04.pdf"
)


# In[2]:


#adata.write("allgenes_bipca_annotated.h5ad")
adata = sc.read("allgenes_bipca_annotated.h5ad", cache=True)
adata


# In[3]:


sc.pl.embedding(
    adata,
    basis="tsne",
    color="leiden",
    legend_loc="on data",
    legend_fontsize="xx-small"
)
sc.pl.embedding(
    adata,
    basis="tsne",
    color="detailed_celltype_annotation",
    legend_loc="on data",
    legend_fontsize="xx-small"
)


# In[4]:


cell_plot = ["41_Otx1_Zic1"]

# Create a temporary column that only keeps the desired cell types; others become NaN.
adata.obs['plot'] = adata.obs['detailed_celltype_annotation'].where(
    adata.obs['detailed_celltype_annotation'].isin(cell_plot)
)

# Plot using the new column.
sc.pl.embedding(
    adata,
    basis="tsne",
    color="plot",
    legend_loc="on data",
    legend_fontsize="xx-small"
)


# In[5]:


neurotransmitter_markers = {
    "GABAergic": ["Gad1", "Gad2", "Slc32a1"],
    "Glutamatergic": ["Slc17a7", "Slc17a6", "Slc17a8"],
    "Cholinergic": ["Chat", "Slc18a3"],
    "Dopaminergic": ["Th", "Ddc", "Slc6a3"],
    "Histaminergic": ["Hdc", "Slc18a2"]
}
all_markers = list({gene for genes in neurotransmitter_markers.values() for gene in genes})

sc.pl.dotplot(
    adata,
    var_names=neurotransmitter_markers,
    groupby='detailed_celltype_annotation',
    layer="scran_normalization",
    dendrogram=True,
    color_map='Reds',
    save="neurotransmitters_7plex.pdf"
)


# In[6]:


sc.pl.tsne(adata, 
           color=["Gad2", "Slc17a6", "Chat", "Th", "Ddc", "Hdc", "Cartpt"], 
           cmap='pink_r', 
           layer="scran_normalization", 
           vmax=1, 
          save="_neurotransmitters_7plex.pdf")


# In[7]:


cl_annotation = adata.obs["detailed_celltype_annotation"]
cl_annotation["00_Hmcn2_Gpr149"]    = "Glut"
cl_annotation["01_Undefined"]       = "GABA_Glut"
cl_annotation["02_Fezf1_Lama1"]     = "Glut"
cl_annotation["03_Rai14_Acvr1c"]    = "Glut"
cl_annotation["04_Crnde_Lmx1a"]     = "Glut"
cl_annotation["05_Qrfpr_Pdzrn3"]    = "Glut"
cl_annotation["06_Tbx3_Npy1r"]      = "GABA"
cl_annotation["07_Nfib_Six3"]       = "GABA"
cl_annotation["08_Lef1_Rspo2"]      = "GABA"
cl_annotation["09_Nfia_Ccbe1"]      = "Glut"
cl_annotation["10_Npy_Agrp"]        = "GABA"
cl_annotation["11_Tafa4_Kank4"]     = "Glut"
cl_annotation["12_Ntng2_Onecut1"]   = "Glut"
cl_annotation["13_Pth2r_Otp"]       = "GABA"
cl_annotation["14_Lhx1_Th"]         = "GABA"
cl_annotation["15_Onecut1_Slc32a1"] = "GABA"
cl_annotation["16_Gal_Sox6"]        = "GABA"
cl_annotation["17_Dkk2_Ebf3"]       = "Glut"
cl_annotation["18_Ust_Bcl11b"]      = "Glut"
cl_annotation["19_Bace2_Tspan18"]   = "Glut"
cl_annotation["20_Tmem114_Npsr1"]   = "GABA"
cl_annotation["21_Bnc2_Glp1r"]      = "GABA"
cl_annotation["22_Crabp1_Vgll3"]    = "GABA"
cl_annotation["23_Ust_Nkd2"]        = "Glut"
cl_annotation["24_Nts_Ecel1"]       = "GABA"
cl_annotation["25_Nxph4_Ebf3"]      = "Glut"
cl_annotation["26_Dach2_Npsr1"]     = "GABA"
cl_annotation["27_Lmx1a_Ttc6"]      = "Glut"
cl_annotation["28_Kiss1_Slc38a3"]   = "Glut"
cl_annotation["29_Dnah11_Pth2r"]    = "GABA"
cl_annotation["30_Trh_Sim1"]        = "Glut"
cl_annotation["31_Lhx9_Npvf"]       = "Glut"
cl_annotation["32_Gucy2c_Ttc6"]     = "Glut"
cl_annotation["33_Chat_Ghrh"]       = "Chat"
cl_annotation["34_Pla2r1_Sim1"]     = "Glut"
cl_annotation["35_Lhx9_Zic1"]       = "Glut"
cl_annotation["36_Ebf2_Crhbp"]      = "Glut"
cl_annotation["37_Qrfp_Qrfprl"]     = "Glut"
cl_annotation["38_Hdc_Prph"]        = "Hist"
cl_annotation["39_Nfib_Nfia"]       = "GABA"
cl_annotation["40_Gpr50_Tac2"]      = "GABA"
cl_annotation["41_Otx1_Zic1"]       = "Glut"
cl_annotation["42_Foxb1_Myo5b"]     = "Glut"
#GLIA
cl_annotation["Astrocyte"]           = "Glia"
cl_annotation["Ependymal_cell"]      = "Glia"
cl_annotation["Endothelial_cell"]    = "Glia"
cl_annotation["Epithelial-like_cell"]= "Glia"
cl_annotation["Microglia"]           = "Glia"
cl_annotation["Oligodendrocyte"]     = "Glia"
cl_annotation["OPC"]                 = "Glia"
cl_annotation["Pericyte"]            = "Glia"
cl_annotation["Pituicyte"]           = "Glia"
cl_annotation["Pituitary_Crhr1"]     = "Pituitary"
cl_annotation["Pituitary_Gnrhr"]     = "Pituitary"
cl_annotation["Pituitary_Mtnr1a"]    = "Pituitary"
cl_annotation["Pituitary_Pappa2"]    = "Pituitary"
cl_annotation["Pituitary_Pax7"]      = "Pituitary"
cl_annotation["Pituitary_Rxfp1"]     = "Pituitary"
cl_annotation["Schwann_cell"]        = "Glia"
cl_annotation["Tanycyte"]            = "Glia"
cl_annotation["Undefined"]           = "Glia"
cl_annotation["VLMC"]                = "Glia"
adata.obs["neurotransmitter"] = adata.obs.detailed_celltype_annotation.map(cl_annotation)


# In[9]:


sc.pl.embedding(
    adata,
    basis="tsne",
    color="neurotransmitter",
    save="_neurotransmitters_7plex_umap.pdf"
)


# In[10]:


from scipy.stats import chi2_contingency
table = pd.crosstab(adata.obs["hash.mcl.ID"], adata.obs["detailed_celltype_annotation"])
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


from scipy.stats import norm
from statsmodels.stats.multitest import multipletests
# Two-tailed p-values for each residual
pvals = 2 * (1 - norm.cdf(np.abs(std_residuals)))
pvals_df = pd.DataFrame(pvals, index=table.index, columns=table.columns)
# Flatten and correct
pvals_flat = pvals_df.values.flatten()
_, pvals_adj, _, _ = multipletests(pvals_flat, method='fdr_bh')

# Reshape back to original shape
pvals_adj_df = pd.DataFrame(pvals_adj.reshape(pvals_df.shape), index=table.index, columns=table.columns)
# Significance mask
significant_mask = pvals_adj_df < 0.001

# Masked residuals (non-significant values set to NaN)
masked_residuals = std_residuals.where(significant_mask)
plt.figure(figsize=(max(10, std_residuals.shape[1] * 0.5),
                    max(6, std_residuals.shape[0] * 0.5)))

sns.heatmap(masked_residuals, annot=True, fmt=".0f", cmap="coolwarm", center=0,
            linewidths=0.5, linecolor='grey', cbar_kws={'label': 'Standardized Residual'})

plt.title("Significant Over-/Underrepresentation (FDR < 0.001)")
plt.xlabel("Cell Type")
plt.ylabel("Animal")
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.show()


# In[13]:


table = pd.crosstab(adata.obs["area"], adata.obs["detailed_celltype_annotation"])
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
plt.ylabel("Area")
plt.xlabel("Cell Type")
plt.tight_layout()
plt.show()


# In[14]:


std_residuals = (table - expected) / (expected**0.5)
n_rows, n_cols = table.shape
fig_width = max(10, n_cols * 0.5)
fig_height = max(6, n_rows * 0.5)
plt.figure(figsize=(fig_width, fig_height))
sns.heatmap(std_residuals, cmap="coolwarm", center=0, annot=True, fmt=".0f")
plt.title("Cell Type Standardized Overrepresentation (Observed - Expected)")
plt.ylabel("Area")
plt.xlabel("Cell Type")
plt.tight_layout()
plt.show()


# In[15]:


# Two-tailed p-values for each residual
pvals = 2 * (1 - norm.cdf(np.abs(std_residuals)))
pvals_df = pd.DataFrame(pvals, index=table.index, columns=table.columns)
# Flatten and correct
pvals_flat = pvals_df.values.flatten()
_, pvals_adj, _, _ = multipletests(pvals_flat, method='fdr_bh')

# Reshape back to original shape
pvals_adj_df = pd.DataFrame(pvals_adj.reshape(pvals_df.shape), index=table.index, columns=table.columns)
# Significance mask
significant_mask = pvals_adj_df < 0.001

# Masked residuals (non-significant values set to NaN)
masked_residuals = std_residuals.where(significant_mask)
plt.figure(figsize=(max(10, std_residuals.shape[1] * 0.5),
                    max(6, std_residuals.shape[0] * 0.5)))

sns.heatmap(masked_residuals, annot=True, fmt=".0f", cmap="coolwarm", center=0,
            linewidths=0.5, linecolor='grey', cbar_kws={'label': 'Standardized Residual'})

plt.title("Significant Over-/Underrepresentation (FDR < 0.001)")
plt.xlabel("Cell Type")
plt.ylabel("Animal")
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.show()


# In[35]:


adata.write("allgenes_bipca_annotated.h5ad")
#adata = sc.read("allgenes_bipca_annotated.h5ad", cache=True)
adata


# In[16]:


###For analysis purposes a few cell types are removed
#01 and 35 for massive misrepresentation across batches
#Undefined for massive misrepresentation from a single animal
#Schwann, Epithelial-like and Pituitary Crhr1/Gnrhr/Pappa2 for less than 10 mappings across all spatial sections
#The remaining pituitary can be reasonably explained to be misrepresented because of dissection variation
downstream = ["01_Undefined", "35_Lhx9_Zic1", "Undefined", "Schwann_cell", "Epithelial-like_cell", "Pituitary_Crhr1", "Pituitary_Gnrhr", "Pituitary_Pappa2"]
mask = ~adata.obs["detailed_celltype_annotation"].isin(downstream)
subset = adata[mask].copy()
subset


# In[17]:


sc.pl.embedding(
    subset,
    basis="tsne",
    color="detailed_celltype_annotation",
    legend_loc="on data",
    legend_fontsize="xx-small",
    save="_downstream_7plex_umap.pdf"
)


# In[20]:


ax = sc.pl.embedding(
    subset,
    basis="tsne",
    color="detailed_celltype_annotation",
    legend_loc="on data",
    legend_fontsize="xx-small",
    return_fig=True
)
# Obtain all text objects from the Axes (these are the on-data labels)
ax=plt.gca()
texts = ax.texts
adjustText.adjust_text(texts, 
                       ax=ax,
                       arrowprops=dict(arrowstyle="-", color='black', lw=0.5)
                      )
#plt.show()
plt.savefig("downstream_7plex_repel.pdf")


# In[21]:


table = pd.crosstab(subset.obs["hash.mcl.ID"], subset.obs["detailed_celltype_annotation"])
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


# In[22]:


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


# In[23]:


table = pd.crosstab(subset.obs["area"], subset.obs["detailed_celltype_annotation"])
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
plt.ylabel("Area")
plt.xlabel("Cell Type")
plt.tight_layout()
plt.show()


# In[24]:


std_residuals = (table - expected) / (expected**0.5)
n_rows, n_cols = table.shape
fig_width = max(10, n_cols * 0.5)
fig_height = max(6, n_rows * 0.5)
plt.figure(figsize=(fig_width, fig_height))
sns.heatmap(std_residuals, cmap="coolwarm", center=0, annot=True, fmt=".0f")
plt.title("Cell Type Standardized Overrepresentation (Observed - Expected)")
plt.ylabel("Area")
plt.xlabel("Cell Type")
plt.tight_layout()
plt.show()


# In[25]:


table = pd.crosstab(subset.obs["sex"], subset.obs["detailed_celltype_annotation"])
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
plt.ylabel("Area")
plt.xlabel("Cell Type")
plt.tight_layout()
plt.show()


# In[33]:


abs_counts = pd.crosstab(subset.obs["detailed_celltype_annotation"], subset.obs["sex"])

# Plot absolute counts
fig1, ax1 = plt.subplots(figsize=(14, 6))
abs_counts.plot.bar(stacked=True, ax=ax1)
ax1.set_xlabel('Cell Type')
ax1.set_ylabel('Number of Cells')
ax1.set_title('Absolute Number of Cells per Cell Type by Sex')
plt.tight_layout()
plt.savefig("downstream_7plex_sex_abs.pdf")

rel_counts = abs_counts.div(abs_counts.sum(axis=1), axis=0)

# Plot relative proportions
fig2, ax2 = plt.subplots(figsize=(14, 6))
rel_counts.plot.bar(stacked=True, ax=ax2, legend=False)
ax2.set_xlabel('Cell Type')
ax2.set_ylabel('Proportion of Cells')
ax2.set_title('Proportion of Cells per Cell Type by Sex')
plt.tight_layout()
plt.savefig("downstream_7plex_sex_rel.pdf")


# In[34]:


abs_counts = pd.crosstab(subset.obs["detailed_celltype_annotation"], subset.obs["area"])

# Plot absolute counts
fig1, ax1 = plt.subplots(figsize=(14, 6))
abs_counts.plot.bar(stacked=True, ax=ax1)
ax1.set_xlabel('Cell Type')
ax1.set_ylabel('Number of Cells')
ax1.set_title('Absolute Number of Cells per Cell Type by Area')
plt.tight_layout()
plt.savefig("downstream_7plex_area_abs.pdf")

rel_counts = abs_counts.div(abs_counts.sum(axis=1), axis=0)

# Plot relative proportions
fig2, ax2 = plt.subplots(figsize=(14, 6))
rel_counts.plot.bar(stacked=True, ax=ax2, legend=False)
ax2.set_xlabel('Cell Type')
ax2.set_ylabel('Proportion of Cells')
ax2.set_title('Proportion of Cells per Cell Type by Area')
plt.tight_layout()
plt.savefig("downstream_7plex_area_rel.pdf")


# In[36]:


subset.write("allgenes_bipca_downstream.h5ad")
#subset = sc.read("allgenes_bipca_annotated.h5ad", cache=True)
subset


# In[ ]:




