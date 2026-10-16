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
import statsmodels.stats.multitest as smm
import statsmodels.api as sm
import statsmodels.formula.api as smf
import adjustText
from scipy.stats import pearsonr, spearmanr


# In[2]:


adata = sc.read("allgenes_bipca_downstream.h5ad")
adata


# In[3]:


adata.obs


# In[6]:


df = adata.obs.copy()
meta_df = pd.crosstab(
    df['hash.mcl.ID'],
    df['area'])
# Boolean mask of non-zero entries, then sum down each column
nonzero_per_area = (meta_df > 0).sum(axis=0)
print(nonzero_per_area)


# In[3]:


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
sc.pl.embedding(
    adata,
    basis="tsne",
    color="neurotransmitter")


# In[10]:


marker_genes = ["Lmx1a", "Lmx1b", "Irx3", "Irx5"]
sc.pl.embedding(adata,basis="tsne",color=marker_genes,layer="scran_normalization", color_map="YlOrBr", save="_lmx_irx.pdf", vmax=1.5)


# In[3]:


marker_genes = ["Lepr", "Calcr", "Glp1r", "Mc4r", "Gipr"]
sc.pl.embedding(adata,basis="tsne",color=marker_genes,layer="scran_normalization", color_map="YlOrBr", save="_lmx_receptors.pdf", vmax=1.5)


# In[4]:


marker_genes = ["Otp", "Meis1", "Lhx9", "Lhx2", "Npsr1", "Samd3"]
sc.pl.embedding(adata,basis="tsne",color=marker_genes,layer="scran_normalization", color_map="YlOrBr", save="_lmx_human_markers.pdf", vmax=1.5)


# In[4]:


# 1) grab your obs and compute the ordering by median
df = adata.obs
order = (
    df.groupby('detailed_celltype_annotation')['cellect_bmi_log10']
      .median()
      .sort_values()
      .index
      .tolist()
)
plt.figure(figsize=(12, 6))
sns.set_style("whitegrid")
sns.pointplot(
    data=df,
    x='detailed_celltype_annotation',
    y='cellect_bmi_log10',
    order=order,
    join=False,         # don’t connect the dots
    capsize=0.1,        # error bar cap size
    errwidth=1,         # error bar line width
    ci='sd'             # use ±1 SD for the error bars
)

# 3) add the intersect line at y=3.1
plt.axhline(3.1, color='red', linestyle='--', label='y = 3.1')
# polish up
plt.xticks(rotation=45, ha='right')
plt.xlabel('Cell Type')
plt.ylabel('log10(Cellect BMI)')
plt.title('Cellect BMI (log10) by Cell Type')
plt.tight_layout()
plt.savefig("downstream_7plex_CELLECT.pdf")


# In[5]:


abs_counts = pd.crosstab(adata.obs["neurotransmitter"], adata.obs["area"])

# Plot absolute counts
fig1, ax1 = plt.subplots(figsize=(14, 6))
abs_counts.plot.bar(stacked=True, ax=ax1)
ax1.set_xlabel('Cell Type')
ax1.set_ylabel('Number of Cells')
ax1.set_title('Absolute Number of Cells per Cell Type by Area')
plt.tight_layout()
plt.savefig("downstream_7plex_neurotransmitter_abs.pdf")

rel_counts = abs_counts.div(abs_counts.sum(axis=1), axis=0)

# Plot relative proportions
fig2, ax2 = plt.subplots(figsize=(14, 6))
rel_counts.plot.bar(stacked=True, ax=ax2, legend=False)
ax2.set_xlabel('Cell Type')
ax2.set_ylabel('Proportion of Cells')
ax2.set_title('Proportion of Cells per Cell Type by Area')
plt.tight_layout()
plt.savefig("downstream_7plex_neurotransmitter_rel.pdf")


# In[6]:


from scipy.stats import chi2_contingency
table = pd.crosstab(adata.obs["hash.mcl.ID"], adata.obs["neurotransmitter"])
# Proportion of each cell type *per animal*
row_props = table.div(table.sum(axis=1), axis=0)

# Proportion of each animal *per cell type*
col_props = table.div(table.sum(axis=0), axis=1)

chi2, p, dof, expected = chi2_contingency(table)
residuals = table - expected
std_residuals = (table - expected) / (expected**0.5)
n_rows, n_cols = table.shape
fig_width = max(10, n_cols * 0.5)
fig_height = max(6, n_rows * 0.5)
plt.figure(figsize=(fig_width, fig_height))
sns.heatmap(std_residuals, cmap="coolwarm", center=0, annot=True, fmt=".0f")
plt.title("Cell Type Standardized Overrepresentation (Observed - Expected)")
plt.ylabel("Animal")
plt.xlabel("Neurotransmitter")
plt.tight_layout()
plt.show()


# In[7]:


viral_gene = 'v147'
gene_idx = adata.var_names.get_loc(viral_gene)
glial_labels = ["Glia"]
# prepare container
stats = []
# Prepare array to hold z-scores
z_scores = np.zeros(adata.n_obs)
z_scores[:] = np.nan  # initialize with NaNs
# Loop over animals
for animal_id in adata.obs['hash.mcl.ID'].unique():
    # Subset cells from this animal
    animal_mask = adata.obs['hash.mcl.ID'] == animal_id
    adata_animal = adata[animal_mask]
    # Identify glial cells for this animal
    glia_mask = adata_animal.obs['neurotransmitter'].isin(glial_labels)
    n_glia    = glia_mask.sum()
    if glia_mask.sum() == 0:
        print(f"Warning: No glial cells for animal {animal_id}")
        continue

    # Extract expression for viral gene
    expr = adata_animal.layers["scran_normalization"][:, gene_idx]
    expr = expr.toarray().flatten() if hasattr(expr, 'toarray') else expr.flatten()
    glia_expr = expr[glia_mask]

    # Compute mean and std from glia
    mean_bg = glia_expr.mean()
    std_bg = glia_expr.std()
    # Add stats
    stats.append({
        'animal': animal_id,
        'n_glia': n_glia,
        'mean_bg': mean_bg,
        'std_bg':  std_bg
    })
    # Avoid division by zero
    if std_bg == 0:
        print(f"Warning: std is 0 for animal {animal_id}")
        continue

    # Z-score for all cells in this animal
    z = (expr - mean_bg) / std_bg

    # Assign to full dataset
    z_scores[animal_mask] = z

# Save into obs
adata.obs['viral_zscore_animal'] = z_scores
df_stats = pd.DataFrame(stats).set_index('animal')


# In[8]:


fig, ax = plt.subplots(figsize=(8,4))

# simple barplot with error bars
x = np.arange(len(df_stats))
means = df_stats['mean_bg'].values
errs  = df_stats['std_bg'].values

ax.bar(x, means, yerr=errs, capsize=5)
ax.set_xticks(x)
ax.set_xticklabels(df_stats.index, rotation=45, ha='right')
ax.set_xlabel('Animal ID')
ax.set_ylabel('Glial background expression (mean ± SD)')
ax.set_title(f'{viral_gene} z-background per Animal')
plt.tight_layout()
plt.savefig("7plex_virus_background_hash.pdf")


# In[9]:


viral_gene = 'v147'
gene_idx = adata.var_names.get_loc(viral_gene)
glial_labels = ["Glia"]
# prepare container
stats = []
# Prepare array to hold z-scores
z_scores = np.zeros(adata.n_obs)
z_scores[:] = np.nan  # initialize with NaNs
# Loop over animals
for animal_id in adata.obs['area'].unique():
    # Subset cells from this animal
    animal_mask = adata.obs['area'] == animal_id
    adata_animal = adata[animal_mask]
    # Identify glial cells for this animal
    glia_mask = adata_animal.obs['neurotransmitter'].isin(glial_labels)
    n_glia    = glia_mask.sum()
    if glia_mask.sum() == 0:
        print(f"Warning: No glial cells for animal {animal_id}")
        continue

    # Extract expression for viral gene
    expr = adata_animal.layers["scran_normalization"][:, gene_idx]
    expr = expr.toarray().flatten() if hasattr(expr, 'toarray') else expr.flatten()
    glia_expr = expr[glia_mask]

    # Compute mean and std from glia
    mean_bg = glia_expr.mean()
    std_bg = glia_expr.std()
    # Add stats
    stats.append({
        'area': animal_id,
        'n_glia': n_glia,
        'mean_bg': mean_bg,
        'std_bg':  std_bg
    })
    # Avoid division by zero
    if std_bg == 0:
        print(f"Warning: std is 0 for animal {animal_id}")
        continue

    # Z-score for all cells in this animal
    z = (expr - mean_bg) / std_bg

    # Assign to full dataset
    z_scores[animal_mask] = z

df_stats = pd.DataFrame(stats).set_index('area')

fig, ax = plt.subplots(figsize=(8,4))

# simple barplot with error bars
x = np.arange(len(df_stats))
means = df_stats['mean_bg'].values
errs  = df_stats['std_bg'].values

ax.bar(x, means, yerr=errs, capsize=5)
ax.set_xticks(x)
ax.set_xticklabels(df_stats.index, rotation=45, ha='right')
ax.set_xlabel('Area')
ax.set_ylabel('Glial background expression (mean ± SD)')
ax.set_title(f'{viral_gene} z-background per injection area')
plt.tight_layout()
plt.savefig("7plex_virus_background_area.pdf")


# In[10]:


plt.figure(figsize=(8, 5))
sns.histplot(adata.obs['viral_zscore_animal'], bins=100, kde=True, color='purple')
plt.xlabel('Viral gene Z-score')
plt.ylabel('Number of cells')
plt.title('Distribution of Viral Gene Z-scores')
plt.axvline(x=2, color='red', linestyle='--', label='Z = 2 threshold')
plt.legend()
plt.tight_layout()
plt.show()


# In[11]:


z_col = 'viral_zscore_animal'
# Compute the 95th percentile cutoff
z_cap = np.nanpercentile(adata.obs[z_col], 99)
adata.obs[f'{z_col}_capped'] = np.clip(adata.obs[z_col], a_min=None, a_max=z_cap)
plt.figure(figsize=(8, 5))
sns.histplot(adata.obs[f'{z_col}_capped'], bins=100, kde=True, color='purple')
plt.xlabel('Viral gene Z-score 99th percentile')
plt.ylabel('Number of cells')
plt.title('Distribution of Viral Gene Z-scores')
plt.axvline(x=2, color='red', linestyle='--', label='Z = 2 threshold')
plt.legend()
plt.tight_layout()
plt.show()


# In[12]:


adata.obs['virus_positive'] = adata.obs['viral_zscore_animal'] > 2
# Create a pivot table of counts
heatmap_data = (
    adata.obs[adata.obs['virus_positive']]
    .groupby(['detailed_celltype_annotation', 'area'], observed=False)
    .size()
    .unstack(fill_value=0)
)
# Total cells per group
total_counts = (
    adata.obs
    .groupby(['detailed_celltype_annotation', 'area'], observed=False)
    .size()
    .unstack(fill_value=0)
)

# Proportion = virus-positive / total
proportion_data = heatmap_data / total_counts

# Calculate dynamic figure height based on number of rows
n_rows = proportion_data.shape[0]
n_cols = proportion_data.shape[1]
row_height = 0.4  # adjust this to control spacing
fig_height = max(6, n_rows * row_height)
fig_width = max(8, n_cols * 0.5)

# Plot with adjusted size
plt.figure(figsize=(fig_width, fig_height))
sns.heatmap(proportion_data, cmap='rocket_r', linewidths=0.5, annot=True, fmt=".2f")
plt.title('Proportion of Virus-Positive Cells by Cell Type and Area')
plt.xlabel('Area')
plt.ylabel('Cell Type')
plt.tight_layout()

plt.savefig("Virus_pos_zscore.pdf")


# In[13]:


mask = proportion_data <= 0.09
cmap = sns.color_palette("rocket_r", as_cmap=True)
cmap.set_bad("white")                    # masked values = white

# 3) plot with grid‐lines
plt.figure(figsize=(fig_width, fig_height))
ax = sns.heatmap(
    proportion_data,
    mask=mask,
    cmap=cmap,
    linewidths=0.8,                       # thickness of your grid lines
    linecolor='gray',                     # color of grid lines
    annot=True,
    fmt=".2f",
    square=True                           # make each cell square
)
ax.set_title('Proportion of Virus-Positive Cells\nby Cell Type and Area')
ax.set_xlabel('Area')
ax.set_ylabel('Cell Type')
plt.tight_layout()
plt.savefig("Virus_pos_zscore_T.pdf")


# In[14]:


z_col = 'viral_zscore_animal'
# Compute the 95th percentile cutoff
z_cap = np.nanpercentile(adata.obs[z_col], 99.9)
adata.obs[f'{z_col}_capped'] = np.clip(adata.obs[z_col], a_min=None, a_max=z_cap)
plt.figure(figsize=(8, 5))
sns.histplot(data=adata.obs, 
             x='viral_zscore_animal_capped',
             kde=True,
             stat='density', 
             common_norm=False)
plt.axvline(2, color='red', linestyle='--', label='Z = 2 cutoff')
plt.title('Histogram of Viral Z-scores')
plt.xlabel('Z-score')
plt.ylabel('Density')
plt.tight_layout()
plt.show()


# In[15]:


z_col = 'viral_zscore_animal'
# Compute the 99th-percentile cap
z_cap = np.nanpercentile(adata.obs[z_col], 99)
adata.obs[f'{z_col}_capped'] = np.clip(adata.obs[z_col], None, z_cap)

plt.figure(figsize=(8, 5))
sns.histplot(
    data=adata.obs,
    x=f'{z_col}_capped',
    hue='neurotransmitter',
    bins=100,
    stat='density',        # or 'count' if you want absolute numbers
    common_norm=False,     # don’t normalize across hues
    multiple='stack',      # <— this stacks the bars
    element='bars',        # bars (default), not step
    kde=False              # drop the KDE for clean stacking
)
plt.axvline(2, color='red', linestyle='--', label='Z = 2 cutoff')
plt.title('Stacked Histogram of Viral Z-scores by Neurotransmitter')
plt.xlabel('Z-score (capped at 99th pct)')
plt.ylabel('Density')
plt.legend(title='Neurotransmitter', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.savefig("Histogram_neurotransmitter.pdf")


# In[16]:


z_col = 'viral_zscore_animal'
z_cap = np.nanpercentile(adata.obs[z_col], 99.9)
ax = sns.histplot(
    data=adata.obs,
    x=f'{z_col}',
    hue='neurotransmitter',
    bins=100,
    stat='density',
    common_norm=False,
    multiple='stack',
    element='bars',
    kde=False,
    binrange=(0, z_cap)    # <-- restrict histogram to [0, z_cap]
)
plt.axvline(2, color='red', linestyle='--', label='Z = 2 cutoff')
# no need for set_xlim if binrange covers it
plt.tight_layout()
plt.savefig("Histogram_neurotransmitter_0cap.pdf")


# In[17]:


df = adata.obs.copy()
rel_contrib = pd.crosstab(
    df['virus_positive'],
    df['neurotransmitter'],
    normalize='index')

rel_contrib.index = rel_contrib.index.map({False: 'Virus Negative', True: 'Virus Positive'})

fig, ax = plt.subplots(figsize=(8, 6))
rel_contrib.plot(
    kind='bar',
    stacked=True,
    ax=ax,
    width=0.7
)
ax.set_xlabel('')
ax.set_ylabel('Proportion of Cells')
ax.set_title('Neurotransmitter Composition\nby Virus Status')
ax.legend(title='Neurotransmitter', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.xticks(rotation=0)
plt.tight_layout()
plt.savefig("Barplot_neurotransmitter.pdf")


# In[18]:


heatmap_data = (
    adata.obs[adata.obs['virus_positive']]
    .groupby(['predicted_spatial_region_filtered_norm', 'area'], observed=False)
    .size()
    .unstack(fill_value=0))
total_counts = (
    adata.obs
    .groupby(['predicted_spatial_region_filtered_norm', 'area'], observed=False)
    .size()
    .unstack(fill_value=0))
proportion_data = heatmap_data / total_counts

n_rows = proportion_data.shape[0]
n_cols = proportion_data.shape[1]
row_height = 0.4
fig_height = max(6, n_rows * row_height)
fig_width = max(8, n_cols * 0.5)
plt.figure(figsize=(fig_width, fig_height))
sns.heatmap(proportion_data, cmap='rocket_r', linewidths=0.5, annot=True, fmt=".2f")
plt.title('Proportion of Virus-Positive Cells by Region and Injection Area')
plt.xlabel('Area')
plt.ylabel('Region')
plt.tight_layout()

plt.savefig("Virus_pos_zscore_regional.pdf")


# In[19]:


spatial_calls = pd.read_csv("projection_calls_spatial_o.csv", index_col=0)
spatial_calls = spatial_calls.T #intersect between called tracings spatial x sn labeltransfer 0.4 cutoff
spatial_calls_t = pd.read_csv("cells_calls_spatial_o.csv", index_col=1)
spatial_calls["total"] = spatial_calls_t["Freq"] #total numbers from seurat_merged$comp1 (filtered 0.4)
# spatial_calls contains the number of positive cells for each projection and each cell type
df_ratio = spatial_calls.div(spatial_calls["total"], axis=0) # relative values


# In[20]:


# projection calls from the single-cell data
heatmap_data = (
    adata.obs[adata.obs['virus_positive']]
    .groupby(['detailed_celltype_annotation', 'area'], observed=False)
    .size()
    .unstack(fill_value=0)
)
total_counts = (
    adata.obs
    .groupby(['detailed_celltype_annotation', 'area'], observed=False)
    .size()
    .unstack(fill_value=0)
)
# Proportion = virus-positive / total
proportion_data = heatmap_data / total_counts


# In[21]:


# Align on both axes, keeping only the intersection
df1_clean, df2_clean = proportion_data.align(df_ratio, join='inner')
df1_clean = df1_clean.filter(regex=r'\d', axis=0)#filtering out glia
df2_clean = df2_clean.filter(regex=r'\d', axis=0)#filtering out glia
# Now df1_clean.index == df2_clean.index and df1_clean.columns == df2_clean.columns
print("Common index   :", df1_clean.index)
print("Common columns :", df1_clean.columns)

# You can then compute e.g. the global correlation
s1 = df1_clean.stack()
s2 = df2_clean.stack()
s1, s2 = s1.align(s2, join="inner")
print("Pearson R:", s1.corr(s2))
global_spearman = s1.corr(s2, method="spearman")
print(f"Global Spearman correlation: {global_spearman:.3f}")


# In[22]:


per_area_corr = df1_clean.corrwith(df2_clean, axis=0)  # correlation column‐by‐column
print(per_area_corr)


# In[23]:


from scipy.stats import linregress

# Build your long‐form DataFrame as before
s1 = df1_clean.stack().rename("val1")
s2 = df2_clean.stack().rename("val2")
df_long = (
    pd.concat([s1, s2], axis=1)
      .reset_index()
      .rename(columns={"level_0": "celltype", "level_1": "region"})
)

regions = df_long["region"].unique()
cmap = plt.get_cmap("tab20", len(regions))
region_to_color = {r: cmap(i) for i, r in enumerate(regions)}

plt.figure(figsize=(8, 6))

for region, grp in df_long.groupby("region"):
    # Drop any rows where val1 or val2 is NA
    grp = grp.dropna(subset=["val1", "val2"])
    x = grp["val1"].values
    y = grp["val2"].values
    col = region_to_color[region]

    # scatter
    plt.scatter(x, y, color=col, alpha=0.6, s=20, label=region)

    # Only fit if there's variance in both x and y
    if len(x) > 1 and np.ptp(x) > 0 and np.ptp(y) > 0:
        try:
            lr = linregress(x, y)
        except Exception:
            continue

        # build line
        x_line = np.array([x.min(), x.max()])
        y_line = lr.intercept + lr.slope * x_line
        plt.plot(x_line, y_line, color=col, lw=1.5)

plt.legend(bbox_to_anchor=(1, 1), title="Region")
plt.xlabel("snRNAseq values")
plt.ylabel("spatial values")
plt.title("Scatter + Linear Fit by Region")
plt.tight_layout()
plt.savefig("scatterplot_correlation_spatial_projection.pdf")
#PAG, PBN, CEA = v36r
#PVT, BST, LHA = v147r


# In[24]:


spatial_calls = pd.read_csv("projection_calls_spatial_z.csv", index_col=0)
spatial_calls = spatial_calls.T #intersect between called tracings spatial x sn labeltransfer 0.4 cutoff
spatial_calls_t = pd.read_csv("cells_calls_spatial_z.csv", index_col=1)
spatial_calls["total"] = spatial_calls_t["Freq"] #total numbers from seurat_merged$comp1 (filtered 0.4)
# spatial_calls contains the number of positive cells for each projection and each cell type
df_ratio = spatial_calls.div(spatial_calls["total"], axis=0) # relative values


# In[25]:


# projection calls from the single-cell data
heatmap_data = (
    adata.obs[adata.obs['virus_positive']]
    .groupby(['detailed_celltype_annotation', 'area'], observed=False)
    .size()
    .unstack(fill_value=0)
)
total_counts = (
    adata.obs
    .groupby(['detailed_celltype_annotation', 'area'], observed=False)
    .size()
    .unstack(fill_value=0)
)
# Proportion = virus-positive / total
proportion_data = heatmap_data / total_counts


# In[26]:


# Align on both axes, keeping only the intersection
df1_clean, df2_clean = proportion_data.align(df_ratio, join='inner')
df1_clean = df1_clean.filter(regex=r'\d', axis=0)#filtering out glia
df2_clean = df2_clean.filter(regex=r'\d', axis=0)#filtering out glia
# Now df1_clean.index == df2_clean.index and df1_clean.columns == df2_clean.columns
print("Common index   :", df1_clean.index)
print("Common columns :", df1_clean.columns)

# You can then compute e.g. the global correlation
s1 = df1_clean.stack()
s2 = df2_clean.stack()
s1, s2 = s1.align(s2, join="inner")
print("Pearson R:", s1.corr(s2))
global_spearman = s1.corr(s2, method="spearman")
print(f"Global Spearman correlation: {global_spearman:.3f}")


# In[27]:


per_area_corr = df1_clean.corrwith(df2_clean, axis=0)  # correlation column‐by‐column
print(per_area_corr)


# In[28]:


# Build your long‐form DataFrame as before
s1 = df1_clean.stack().rename("val1")
s2 = df2_clean.stack().rename("val2")
df_long = (
    pd.concat([s1, s2], axis=1)
      .reset_index()
      .rename(columns={"level_0": "celltype", "level_1": "region"})
)

regions = df_long["region"].unique()
cmap = plt.get_cmap("tab20", len(regions))
region_to_color = {r: cmap(i) for i, r in enumerate(regions)}

plt.figure(figsize=(8, 6))

for region, grp in df_long.groupby("region"):
    # Drop any rows where val1 or val2 is NA
    grp = grp.dropna(subset=["val1", "val2"])
    x = grp["val1"].values
    y = grp["val2"].values
    col = region_to_color[region]

    # scatter
    plt.scatter(x, y, color=col, alpha=0.6, s=20, label=region)

    # Only fit if there's variance in both x and y
    if len(x) > 1 and np.ptp(x) > 0 and np.ptp(y) > 0:
        try:
            lr = linregress(x, y)
        except Exception:
            continue

        # build line
        x_line = np.array([x.min(), x.max()])
        y_line = lr.intercept + lr.slope * x_line
        plt.plot(x_line, y_line, color=col, lw=1.5)

plt.legend(bbox_to_anchor=(1, 1), title="Region")
plt.xlabel("snRNAseq values")
plt.ylabel("spatial values")
plt.title("Scatter + Linear Fit by Region")
plt.tight_layout()
plt.savefig("scatterplot_correlation_spatial_zscore.pdf")
#PAG, PBN, CEA = v36r
#PVT, BST, LHA = v147r


# In[29]:


spatial_calls = pd.read_csv("projection_calls_spatial_c.csv", index_col=0)
spatial_calls = spatial_calls.T #intersect between called tracings spatial x sn labeltransfer 0.4 cutoff
spatial_calls_t = pd.read_csv("cells_calls_spatial_c.csv", index_col=1)
spatial_calls["total"] = spatial_calls_t["Freq"] #total numbers from seurat_merged$comp1 (filtered 0.4)
# spatial_calls contains the number of positive cells for each projection and each cell type
df_ratio = spatial_calls.div(spatial_calls["total"], axis=0) # relative values
# projection calls from the single-cell data
heatmap_data = (
    adata.obs[adata.obs['virus_positive']]
    .groupby(['detailed_celltype_annotation', 'area'], observed=False)
    .size()
    .unstack(fill_value=0)
)
total_counts = (
    adata.obs
    .groupby(['detailed_celltype_annotation', 'area'], observed=False)
    .size()
    .unstack(fill_value=0)
)
# Proportion = virus-positive / total
proportion_data = heatmap_data / total_counts
# Align on both axes, keeping only the intersection
df1_clean, df2_clean = proportion_data.align(df_ratio, join='inner')
df1_clean = df1_clean.filter(regex=r'\d', axis=0)#filtering out glia
df2_clean = df2_clean.filter(regex=r'\d', axis=0)#filtering out glia
# Now df1_clean.index == df2_clean.index and df1_clean.columns == df2_clean.columns
print("Common index   :", df1_clean.index)
print("Common columns :", df1_clean.columns)

# You can then compute e.g. the global correlation
s1 = df1_clean.stack()
s2 = df2_clean.stack()
s1, s2 = s1.align(s2, join="inner")
print("Pearson R:", s1.corr(s2))
global_spearman = s1.corr(s2, method="spearman")
print(f"Global Spearman correlation: {global_spearman:.3f}")
per_area_corr = df1_clean.corrwith(df2_clean, axis=0)  # correlation column‐by‐column
print(per_area_corr)


# In[30]:


# Build your long‐form DataFrame as before
s1 = df1_clean.stack().rename("val1")
s2 = df2_clean.stack().rename("val2")
df_long = (
    pd.concat([s1, s2], axis=1)
      .reset_index()
      .rename(columns={"level_0": "celltype", "level_1": "region"})
)

regions = df_long["region"].unique()
cmap = plt.get_cmap("tab20", len(regions))
region_to_color = {r: cmap(i) for i, r in enumerate(regions)}

plt.figure(figsize=(8, 6))

for region, grp in df_long.groupby("region"):
    # Drop any rows where val1 or val2 is NA
    grp = grp.dropna(subset=["val1", "val2"])
    x = grp["val1"].values
    y = grp["val2"].values
    col = region_to_color[region]

    # scatter
    plt.scatter(x, y, color=col, alpha=0.6, s=20, label=region)

    # Only fit if there's variance in both x and y
    if len(x) > 1 and np.ptp(x) > 0 and np.ptp(y) > 0:
        try:
            lr = linregress(x, y)
        except Exception:
            continue

        # build line
        x_line = np.array([x.min(), x.max()])
        y_line = lr.intercept + lr.slope * x_line
        plt.plot(x_line, y_line, color=col, lw=1.5)

plt.legend(bbox_to_anchor=(1, 1), title="Region")
plt.xlabel("snRNAseq values")
plt.ylabel("spatial values")
plt.title("Scatter + Linear Fit by Region")
plt.tight_layout()
plt.savefig("scatterplot_correlation_spatial_coronal.pdf")
#PAG, PBN, CEA = v36r
#PVT, BST, LHA = v147r


# In[26]:


outdir = "area_pickles"
os.makedirs(outdir, exist_ok=True)

for area in adata.obs["area"].unique():
    #subset AnnData
    mask = adata.obs["area"] == area
    sub = adata[mask]

    #build DataFrames
    counts = pd.DataFrame(
        sub.layers["scran_normalization"].toarray(),
        index=sub.obs_names,
        columns=sub.var_names
    )
    meta1 = sub.obs["manual_celltype_annotation"].to_frame()
    meta2 = sub.obs["detailed_celltype_annotation"].to_frame()

    #dump to pickle
    counts.to_pickle(os.path.join(outdir, f"{area}_counts.pkl"))
    meta1.to_pickle(os.path.join(outdir, f"{area}_manual.pkl"))
    meta2.to_pickle(os.path.join(outdir, f"{area}_detailed.pkl"))

    print(f"Saved pickles for area={area!r}")


# In[3]:


indir = "area_cellex_results"
#Loop through each CSV, pull out the v147 row, and build a per‐area DataFrame
records = []
for fn in os.listdir(indir):
    if not fn.endswith("_detailed_cellex.csv"):
        continue
    area = fn.replace("_detailed_cellex.csv","")
    df   = pd.read_csv(os.path.join(indir, fn), index_col=0)
    if "v147" not in df.index:
        continue
    # extract the row as a Series of cell‐type → ESMU
    s = df.loc["v147"].rename(area)
    records.append(s)

#Combine into one DataFrame: rows=areas, cols=cell‐types
combined = pd.DataFrame(records)
#Make sure columns match your actual adata.obs categories
celltypes = sorted(adata.obs["detailed_celltype_annotation"].unique())
combined = combined.reindex(columns=celltypes, fill_value=0)
combined = combined.sort_index()
combined_swap = combined.T
print(combined)


# In[28]:


#Calculate dynamic figure height based on number of rows
n_rows = combined_swap.shape[0]
n_cols = combined_swap.shape[1]
row_height = 0.4  # adjust this to control spacing
fig_height = max(6, n_rows * row_height)
fig_width = max(8, n_cols * 0.5)
plt.figure(figsize=(fig_width, fig_height))
sns.heatmap(combined_swap, cmap='rocket_r', linewidths=0.5, annot=True, fmt=".2f")
plt.xlabel("Area")
plt.ylabel("Cell Type")
plt.title("v147 Enrichment (ESMU) by Area × Cell Type")
plt.tight_layout()
plt.savefig("ESmu_detailed_v147.pdf")


# In[33]:


# projection calls from the single-cell data
heatmap_data = (
    adata.obs[adata.obs['virus_positive']]
    .groupby(['detailed_celltype_annotation', 'area'], observed=False)
    .size()
    .unstack(fill_value=0)
)
total_counts = (
    adata.obs
    .groupby(['detailed_celltype_annotation', 'area'], observed=False)
    .size()
    .unstack(fill_value=0)
)
# Proportion = virus-positive / total
proportion_data = heatmap_data / total_counts
# Align on both axes, keeping only the intersection
df1_clean, df2_clean = proportion_data.align(combined_swap, join='inner')
df1_clean = df1_clean.filter(regex=r'\d', axis=0)#filtering out glia
df2_clean = df2_clean.filter(regex=r'\d', axis=0)#filtering out glia
# Now df1_clean.index == df2_clean.index and df1_clean.columns == df2_clean.columns
print("Common index   :", df1_clean.index)
print("Common columns :", df1_clean.columns)

# You can then compute e.g. the global correlation
s1 = df1_clean.stack()
s2 = df2_clean.stack()

s1.index.names = [None] * s1.index.nlevels
s2.index.names = [None] * s2.index.nlevels
s1, s2 = s1.align(s2, join="inner")
print("Pearson R:", s1.corr(s2))
global_spearman = s1.corr(s2, method="spearman")
print(f"Global Spearman correlation: {global_spearman:.3f}")
per_area_corr = df1_clean.corrwith(df2_clean, axis=0)  # correlation column‐by‐column
print(per_area_corr)


# In[34]:


# Build your long‐form DataFrame as before
s1 = df1_clean.stack().rename("val1")
s2 = df2_clean.stack().rename("val2")
s1 = s1.filter(regex=r'\d', axis=0)
s2 = s2.filter(regex=r'\d', axis=0)
df_long = (
    pd.concat([s1, s2], axis=1)
      .reset_index()
      .rename(columns={"level_0": "celltype", "level_1": "region"})
)

regions = df_long["region"].unique()
cmap = plt.get_cmap("tab20", len(regions))
region_to_color = {r: cmap(i) for i, r in enumerate(regions)}

plt.figure(figsize=(8, 6))

for region, grp in df_long.groupby("region"):
    # Drop any rows where val1 or val2 is NA
    grp = grp.dropna(subset=["val1", "val2"])
    x = grp["val1"].values
    y = grp["val2"].values
    col = region_to_color[region]

    # scatter
    plt.scatter(x, y, color=col, alpha=0.6, s=20, label=region)

    # Only fit if there's variance in both x and y
    if len(x) > 1 and np.ptp(x) > 0 and np.ptp(y) > 0:
        try:
            lr = linregress(x, y)
        except Exception:
            continue

        # build line
        x_line = np.array([x.min(), x.max()])
        y_line = lr.intercept + lr.slope * x_line
        plt.plot(x_line, y_line, color=col, lw=1.5)

plt.legend(bbox_to_anchor=(1, 1), title="Region")
plt.xlabel("z-score proportion values")
plt.ylabel("ESmu v147 values")
plt.title("Scatter + Linear Fit by Region")
plt.tight_layout()
plt.savefig("scatterplot_correlation_ESMu.pdf")
#PAG, PBN, CEA = v36r
#PVT, BST, LHA = v147r


# In[26]:


indir = "area_cellex_results"
runs  = ["spatial1", "spatial2", "spatial3"]
genes = ["p36-NLS-Cre", "p147-mCherry-2A-iCre"]

# map each (run → gene) to your desired area label
gene_area_map = {
    "spatial1": {
        "p36-NLS-Cre":          "PAG",
        "p147-mCherry-2A-iCre": "BST",
    },
    "spatial2": {
        "p36-NLS-Cre":          "PBN",
        "p147-mCherry-2A-iCre": "LHA",
    },
    "spatial3": {
        "p36-NLS-Cre":          "CEA",
        "p147-mCherry-2A-iCre": "PVT",
    },
}

records = []
for run in runs:
    fn = f"{run}_detailed_cellex.csv"
    path = os.path.join(indir, fn)
    df   = pd.read_csv(path, index_col=0)  # rows=genes, cols=cell‐types

    for gene in genes:
        if gene not in df.index:
            raise KeyError(f"Gene {gene!r} not found in {fn}")
        area = gene_area_map[run][gene]
        # extract that row (a Series: cell‐type → enrichment), name it by area
        s = df.loc[gene].rename(area)
        records.append(s)

# build combined DataFrame: rows=areas, cols=cell‐types
combined_spatial = pd.DataFrame(records)

# ensure you have all your cell‐types in the right order (pull from adata if needed)
celltypes = sorted(adata.obs["detailed_celltype_annotation"].unique())
combined_spatial = combined_spatial.reindex(columns=celltypes, fill_value=0)
combined_spatial_swap = combined_spatial.T
combined_spatial.to_csv("region_celltype_scores_cellex_spatial.csv")
print(combined_spatial)


# In[25]:


#Calculate dynamic figure height based on number of rows
n_rows = combined_spatial_swap.shape[0]
n_cols = combined_spatial_swap.shape[1]
row_height = 0.4  # adjust this to control spacing
fig_height = max(6, n_rows * row_height)
fig_width = max(8, n_cols * 0.5)
plt.figure(figsize=(fig_width, fig_height))
sns.heatmap(combined_spatial_swap, cmap='rocket_r', linewidths=0.5, annot=True, fmt=".2f")
plt.xlabel("Area")
plt.ylabel("Cell Type")
plt.title("v147 or v63 Enrichment (ESMU) by Injection Area × Cell Type")
plt.tight_layout()
plt.savefig("ESmu_detailed_spatial.pdf")


# In[37]:


# projection calls from the single-cell data
heatmap_data = (
    adata.obs[adata.obs['virus_positive']]
    .groupby(['detailed_celltype_annotation', 'area'], observed=False)
    .size()
    .unstack(fill_value=0)
)
total_counts = (
    adata.obs
    .groupby(['detailed_celltype_annotation', 'area'], observed=False)
    .size()
    .unstack(fill_value=0)
)
# Proportion = virus-positive / total
proportion_data = heatmap_data / total_counts
# Align on both axes, keeping only the intersection
df1_clean, df2_clean = proportion_data.align(combined_spatial_swap, join='inner')
df1_clean = df1_clean.filter(regex=r'\d', axis=0)#filtering out glia
df2_clean = df2_clean.filter(regex=r'\d', axis=0)#filtering out glia
# Now df1_clean.index == df2_clean.index and df1_clean.columns == df2_clean.columns
print("Common index   :", df1_clean.index)
print("Common columns :", df1_clean.columns)

# You can then compute e.g. the global correlation
s1 = df1_clean.stack()
s2 = df2_clean.stack()

s1.index.names = [None] * s1.index.nlevels
s2.index.names = [None] * s2.index.nlevels
s1, s2 = s1.align(s2, join="inner")
print("Pearson R:", s1.corr(s2))
global_spearman = s1.corr(s2, method="spearman")
print(f"Global Spearman correlation: {global_spearman:.3f}")
per_area_corr = df1_clean.corrwith(df2_clean, axis=0)  # correlation column‐by‐column
print(per_area_corr)


# In[38]:


# Build your long‐form DataFrame as before
s1 = df1_clean.stack().rename("val1")
s2 = df2_clean.stack().rename("val2")
s1 = s1.filter(regex=r'\d', axis=0)
s2 = s2.filter(regex=r'\d', axis=0)
df_long = (
    pd.concat([s1, s2], axis=1)
      .reset_index()
      .rename(columns={"level_0": "celltype", "level_1": "region"})
)

regions = df_long["region"].unique()
cmap = plt.get_cmap("tab20", len(regions))
region_to_color = {r: cmap(i) for i, r in enumerate(regions)}

plt.figure(figsize=(8, 6))

for region, grp in df_long.groupby("region"):
    # Drop any rows where val1 or val2 is NA
    grp = grp.dropna(subset=["val1", "val2"])
    x = grp["val1"].values
    y = grp["val2"].values
    col = region_to_color[region]

    # scatter
    plt.scatter(x, y, color=col, alpha=0.6, s=20, label=region)

    # Only fit if there's variance in both x and y
    if len(x) > 1 and np.ptp(x) > 0 and np.ptp(y) > 0:
        try:
            lr = linregress(x, y)
        except Exception:
            continue

        # build line
        x_line = np.array([x.min(), x.max()])
        y_line = lr.intercept + lr.slope * x_line
        plt.plot(x_line, y_line, color=col, lw=1.5)

plt.legend(bbox_to_anchor=(1, 1), title="Region")
plt.xlabel("z-score proportion values")
plt.ylabel("ESmu spatial values")
plt.title("Scatter + Linear Fit by Region")
plt.tight_layout()
plt.savefig("scatterplot_correlation_ESMu_spatial.pdf")


# In[39]:


# Align on both axes, keeping only the intersection
df1_clean, df2_clean = combined_swap.align(combined_spatial_swap, join='inner')
df1_clean = df1_clean.filter(regex=r'\d', axis=0)#filtering out glia
df2_clean = df2_clean.filter(regex=r'\d', axis=0)#filtering out glia
# Now df1_clean.index == df2_clean.index and df1_clean.columns == df2_clean.columns
print("Common index   :", df1_clean.index)
print("Common columns :", df1_clean.columns)

# You can then compute e.g. the global correlation
s1 = df1_clean.stack()
s2 = df2_clean.stack()

s1.index.names = [None] * s1.index.nlevels
s2.index.names = [None] * s2.index.nlevels
s1, s2 = s1.align(s2, join="inner")
print("Pearson R:", s1.corr(s2))
global_spearman = s1.corr(s2, method="spearman")
print(f"Global Spearman correlation: {global_spearman:.3f}")
per_area_corr = df1_clean.corrwith(df2_clean, axis=0)  # correlation column‐by‐column
print(per_area_corr)


# In[40]:


# Build your long‐form DataFrame as before
s1 = df1_clean.stack().rename("val1")
s2 = df2_clean.stack().rename("val2")
s1 = s1.filter(regex=r'\d', axis=0)
s2 = s2.filter(regex=r'\d', axis=0)
df_long = (
    pd.concat([s1, s2], axis=1)
      .reset_index()
      .rename(columns={"level_0": "celltype", "level_1": "region"})
)
regions = df_long["region"].unique()
cmap = plt.get_cmap("tab20", len(regions))
region_to_color = {r: cmap(i) for i, r in enumerate(regions)}

plt.figure(figsize=(8, 6))

for region, grp in df_long.groupby("region"):
    # Drop any rows where val1 or val2 is NA
    grp = grp.dropna(subset=["val1", "val2"])
    x = grp["val1"].values
    y = grp["val2"].values
    col = region_to_color[region]

    # scatter
    plt.scatter(x, y, color=col, alpha=0.6, s=20, label=region)

    # Only fit if there's variance in both x and y
    if len(x) > 1 and np.ptp(x) > 0 and np.ptp(y) > 0:
        try:
            lr = linregress(x, y)
        except Exception:
            continue

        # build line
        x_line = np.array([x.min(), x.max()])
        y_line = lr.intercept + lr.slope * x_line
        plt.plot(x_line, y_line, color=col, lw=1.5)

plt.legend(bbox_to_anchor=(1, 1), title="Region")
plt.xlabel("ESmu spatial values")
plt.ylabel("ESmu v147 values")
plt.title("Scatter + Linear Fit by Region")
plt.tight_layout()
plt.savefig("scatterplot_correlation_ESMu_ESMu.pdf")


# In[41]:


# Ensure both dataframes have the same projection-area columns
combined_swap_clean = combined_swap.filter(regex=r'\d', axis=0)#filtering out glia
combined_spatial_swap_clean = combined_spatial_swap.filter(regex=r'\d', axis=0)#filtering out glia
areas = combined_swap_clean.columns.intersection(combined_spatial_swap_clean.columns)
celltypes = combined_swap_clean.index.tolist()

# x locations for the cell types
x = np.arange(len(celltypes))
width = 0.35  # width of each bar

for area in areas:
    swap_vals    = combined_swap[area].reindex(celltypes).values
    spatial_vals = combined_spatial_swap[area].reindex(celltypes).values

    fig, ax = plt.subplots(figsize=(max(8, len(celltypes)*0.3), 6))
    ax.bar(x - width/2, swap_vals,    width, label='Single-cell',    alpha=0.8)
    ax.bar(x + width/2, spatial_vals, width, label='Spatial', alpha=0.8)

    ax.set_xticks(x)
    ax.set_xticklabels(celltypes, rotation=90, ha='right')
    ax.set_ylabel('Enrichment Score')
    ax.set_title(f'Enrichment by Cell Type for {area}')
    ax.legend(title='Dataset')

    plt.tight_layout()
    plt.show()


# In[43]:


# 1) Build a mapping from each cell-type to its neurotransmitter category
cluster2neu = (
    adata.obs[['detailed_celltype_annotation','neurotransmitter']]
        .drop_duplicates()
        .set_index('detailed_celltype_annotation')['neurotransmitter']
)

# 2) All the neurotransmitter categories we care about
neu_categories = adata.obs['neurotransmitter'].unique()

# 3) Loop over areas and count positives/negatives
records = []
threshold = 0.09

for area in combined_swap.columns:
    # which clusters are “positive” vs “negative”
    positive_clusters = combined_swap.index[combined_swap[area] > threshold]
    negative_clusters = combined_swap.index[combined_swap[area] <= threshold]

    # count how many of those clusters fall into each neurotransmitter
    pos_counts = cluster2neu.loc[positive_clusters].value_counts()
    neg_counts = cluster2neu.loc[negative_clusters].value_counts()

    for neu in neu_categories:
        records.append({
            'area': area,
            'neurotransmitter': neu,
            'positive': int(pos_counts.get(neu, 0)),
            'negative': int(neg_counts.get(neu, 0))
        })

# 4) Build your final crosstable DataFrame
df_counts = pd.DataFrame(records)

# Multi-indexed by area & neurotransmitter, with columns positive/negative
df_counts = df_counts.set_index(['area','neurotransmitter'])[['positive','negative']]

print(df_counts)


# In[45]:


mask = df_counts.index.get_level_values('neurotransmitter').isin(['GABA','Glut'])
df_filtered = df_counts[mask]

print(df_filtered)


# In[47]:


# turn df_counts into a “long” DataFrame
df_long = (
    df_filtered
      .reset_index()  # columns = ['area','neurotransmitter','positive','negative']
      .melt(
         id_vars=['area','neurotransmitter'],
         value_vars=['positive','negative'],
         var_name='status',
         value_name='count'
      )
)

# catplot: one subplot per area
g = sns.catplot(
    data=df_long,
    x='neurotransmitter',
    y='count',
    hue='status',
    col='area',
    kind='bar',
    col_wrap=4,       # e.g. wrap into 2 rows of 4
    height=4,
    aspect=1
)

# polish
g.set_titles("{col_name}")
g.set_axis_labels("Neurotransmitter", "Number of Cell‐types")
for ax in g.axes.flatten():
    ax.tick_params(axis='x', rotation=45)
g._legend.set_title("Status")
plt.tight_layout()
plt.show()


# In[51]:


df_frac = df_filtered.copy()
df_frac['total'] = df_frac['positive'] + df_frac['negative']
df_frac['pos_frac'] = df_frac['positive'] / df_frac['total']
df_frac['neg_frac'] = df_frac['negative'] / df_frac['total']

# keep only the two fraction columns
df_frac = df_frac[['pos_frac','neg_frac']].rename(columns={
    'pos_frac':'Positive',
    'neg_frac':'Negative'
})

# 2B) seaborn facet grid: all areas in one figure
# turn into long form
df_long = (
    df_frac
      .reset_index()
      .melt(
        id_vars=['area','neurotransmitter'],
        value_vars=['Positive','Negative'],
        var_name='Status',
        value_name='Fraction'
      )
)

g = sns.catplot(
    data=df_long,
    x='neurotransmitter',
    y='Fraction',
    hue='Status',
    col='area',
    kind='bar',
    col_wrap=4,
    height=4,
    aspect=1
)
g.set_titles('{col_name}')
g.set_axis_labels('Neurotransmitter','Fraction')
for ax in g.axes.flatten():
    ax.set_ylim(0,1)
    ax.tick_params(axis='x', rotation=45)
g._legend.set_title('Status')
plt.tight_layout()
plt.show()


# In[52]:


# 1) Build a mapping from each cell-type to its neurotransmitter category
cluster2neu = (
    adata.obs[['detailed_celltype_annotation','neurotransmitter']]
        .drop_duplicates()
        .set_index('detailed_celltype_annotation')['neurotransmitter']
)

# 2) All the neurotransmitter categories we care about
neu_categories = adata.obs['neurotransmitter'].unique()

# 3) Loop over areas and count positives/negatives
records = []
threshold = 0.09

for area in combined_spatial_swap.columns:
    # which clusters are “positive” vs “negative”
    positive_clusters = combined_spatial_swap.index[combined_spatial_swap[area] > threshold]
    negative_clusters = combined_spatial_swap.index[combined_spatial_swap[area] <= threshold]

    # count how many of those clusters fall into each neurotransmitter
    pos_counts = cluster2neu.loc[positive_clusters].value_counts()
    neg_counts = cluster2neu.loc[negative_clusters].value_counts()

    for neu in neu_categories:
        records.append({
            'area': area,
            'neurotransmitter': neu,
            'positive': int(pos_counts.get(neu, 0)),
            'negative': int(neg_counts.get(neu, 0))
        })

# 4) Build your final crosstable DataFrame
df_counts = pd.DataFrame(records)

# Multi-indexed by area & neurotransmitter, with columns positive/negative
df_counts = df_counts.set_index(['area','neurotransmitter'])[['positive','negative']]

print(df_counts)


# In[53]:


mask = df_counts.index.get_level_values('neurotransmitter').isin(['GABA','Glut'])
df_filtered = df_counts[mask]
# turn df_counts into a “long” DataFrame
df_long = (
    df_filtered
      .reset_index()  # columns = ['area','neurotransmitter','positive','negative']
      .melt(
         id_vars=['area','neurotransmitter'],
         value_vars=['positive','negative'],
         var_name='status',
         value_name='count'
      )
)

# catplot: one subplot per area
g = sns.catplot(
    data=df_long,
    x='neurotransmitter',
    y='count',
    hue='status',
    col='area',
    kind='bar',
    col_wrap=4,       # e.g. wrap into 2 rows of 4
    height=4,
    aspect=1
)

# polish
g.set_titles("{col_name}")
g.set_axis_labels("Neurotransmitter", "Number of Cell‐types")
for ax in g.axes.flatten():
    ax.tick_params(axis='x', rotation=45)
g._legend.set_title("Status")
plt.tight_layout()
plt.show()


# In[54]:


spatial_calls = pd.read_csv("projection_calls_spatial_z.csv", index_col=0)
spatial_calls = spatial_calls.T #intersect between called tracings spatial x sn labeltransfer 0.4 cutoff
spatial_calls_t = pd.read_csv("cells_calls_spatial_z.csv", index_col=1)
spatial_calls["total"] = spatial_calls_t["Freq"] #total numbers from seurat_merged$comp1 (filtered 0.4)
# spatial_calls contains the number of positive cells for each projection and each cell type
df_ratio = spatial_calls.div(spatial_calls["total"], axis=0) # relative values

# Align on both axes, keeping only the intersection
df1_clean, df2_clean = df_ratio.align(combined_spatial_swap, join='inner')
df1_clean = df1_clean.filter(regex=r'\d', axis=0)#filtering out glia
df2_clean = df2_clean.filter(regex=r'\d', axis=0)#filtering out glia
# Now df1_clean.index == df2_clean.index and df1_clean.columns == df2_clean.columns
print("Common index   :", df1_clean.index)
print("Common columns :", df1_clean.columns)

# You can then compute e.g. the global correlation
s1 = df1_clean.stack()
s2 = df2_clean.stack()

s1.index.names = [None] * s1.index.nlevels
s2.index.names = [None] * s2.index.nlevels
s1, s2 = s1.align(s2, join="inner")
print("Pearson R:", s1.corr(s2))
global_spearman = s1.corr(s2, method="spearman")
print(f"Global Spearman correlation: {global_spearman:.3f}")
per_area_corr = df1_clean.corrwith(df2_clean, axis=0)  # correlation column‐by‐column
print(per_area_corr)


# In[55]:


# Build your long‐form DataFrame as before
s1 = df1_clean.stack().rename("val1")
s2 = df2_clean.stack().rename("val2")
s1 = s1.filter(regex=r'\d', axis=0)
s2 = s2.filter(regex=r'\d', axis=0)
df_long = (
    pd.concat([s1, s2], axis=1)
      .reset_index()
      .rename(columns={"level_0": "celltype", "level_1": "region"})
)

regions = df_long["region"].unique()
cmap = plt.get_cmap("tab20", len(regions))
region_to_color = {r: cmap(i) for i, r in enumerate(regions)}

plt.figure(figsize=(8, 6))

for region, grp in df_long.groupby("region"):
    # Drop any rows where val1 or val2 is NA
    grp = grp.dropna(subset=["val1", "val2"])
    x = grp["val1"].values
    y = grp["val2"].values
    col = region_to_color[region]

    # scatter
    plt.scatter(x, y, color=col, alpha=0.6, s=20, label=region)

    # Only fit if there's variance in both x and y
    if len(x) > 1 and np.ptp(x) > 0 and np.ptp(y) > 0:
        try:
            lr = linregress(x, y)
        except Exception:
            continue

        # build line
        x_line = np.array([x.min(), x.max()])
        y_line = lr.intercept + lr.slope * x_line
        plt.plot(x_line, y_line, color=col, lw=1.5)

plt.legend(bbox_to_anchor=(1, 1), title="Region")
plt.xlabel("z-score proportion values")
plt.ylabel("ESmu v147 values")
plt.title("Scatter + Linear Fit by Region")
plt.tight_layout()
plt.savefig("scatterplot_correlation_ESMu_zscore_spatial.pdf")
#PAG, PBN, CEA = v36r
#PVT, BST, LHA = v147r


# In[3]:


neurons = [
    '00_Hmcn2_Gpr149',
    '02_Fezf1_Lama1',
    '03_Rai14_Acvr1c',
    '04_Crnde_Lmx1a',
    '05_Qrfpr_Pdzrn3',
    '06_Tbx3_Npy1r',
    '07_Nfib_Six3',
    '08_Lef1_Rspo2',
    '09_Nfia_Ccbe1',
    '10_Npy_Agrp',
    '11_Tafa4_Kank4',
    '12_Ntng2_Onecut1',
    '13_Pth2r_Otp',
    '14_Lhx1_Th',
    '15_Onecut1_Slc32a1',
    '16_Gal_Sox6',
    '17_Dkk2_Ebf3',
    '18_Ust_Bcl11b',
    '19_Bace2_Tspan18',
    '20_Tmem114_Npsr1',
    '21_Bnc2_Glp1r',
    '22_Crabp1_Vgll3',
    '23_Ust_Nkd2',
    '24_Nts_Ecel1',
    '25_Nxph4_Ebf3',
    '26_Dach2_Npsr1',
    '27_Lmx1a_Ttc6',
    '28_Kiss1_Slc38a3',
    '29_Dnah11_Pth2r',
    '30_Trh_Sim1',
    '31_Lhx9_Npvf',
    '32_Gucy2c_Ttc6',
    '34_Pla2r1_Sim1',
    '36_Ebf2_Crhbp',
    '37_Qrfp_Qrfprl',
    '38_Hdc_Prph',
    '39_Nfib_Nfia',
    '40_Gpr50_Tac2',
    '41_Otx1_Zic1',
    '42_Foxb1_Myo5b',
    'Pituitary_Mtnr1a',
    'Pituitary_Pax7',
    'Pituitary_Rxfp1']
glia = ['Microglia', 
        'OPC', 
        'Oligodendrocyte', 
        'Astrocyte', 
        'VLMC', 
        'Pituicyte', 
        'Ependymal_cell', 
        'Tanycyte', 
        'Pericyte', 
        'Endothelial_cell']


# In[22]:


indir = "detailed_cellex_results"
areas = ['PAG', 'PBN', 'PVT', 'BST', 'CEA', 'PVH', 'LHA']
#Loop through each CSV, pull out the v147 row, and build a per‐area DataFrame
records = []
for fn in os.listdir(indir):
    if not fn.endswith("_cellex.csv"):
        continue
    area = fn[:3]
    neur = fn[4:]
    neur = neur.replace("_cellex.csv","")
    df   = pd.read_csv(os.path.join(indir, fn), index_col=0)
    if "v147" not in df.index:
        value = 0.00
        records.append((area, neur, value))
        continue
    # extract the row as a Series of cell‐type → ESMU
    value = df.at["v147", neur]
    records.append((area, neur, value))

# Convert to tidy DataFrame
summary = pd.DataFrame(records, columns=["area", "neur", "value"])

# Pivot into wide format: rows=areas, cols=neurons, values=extracted numbers
result = summary.pivot(index="area", columns="neur", values="value")
result = result.sort_index()
combined_swap = result.T
print(result)


# In[23]:


#Calculate dynamic figure height based on number of rows
n_rows = combined_swap.shape[0]
n_cols = combined_swap.shape[1]
row_height = 0.4  # adjust this to control spacing
fig_height = max(6, n_rows * row_height)
fig_width = max(8, n_cols * 0.5)
plt.figure(figsize=(fig_width, fig_height))
sns.heatmap(combined_swap, cmap='rocket_r', linewidths=0.5, annot=True, fmt=".2f")
plt.xlabel("Area")
plt.ylabel("Cell Type")
plt.title("v147 Enrichment (ESMU) by Area × Individual Cell Type")
plt.tight_layout()
plt.savefig("ESmu_individual_v147.pdf")


# In[4]:


# Extract the v147 counts from the "counts" layer
v147_counts = adata.layers["counts"][:, adata.var_names == "v147"]

# Convert to a flat 1D array
v147_counts = v147_counts.A1 if hasattr(v147_counts, "A1") else v147_counts.toarray().ravel()

# Add to obs
adata.obs["v147_counts"] = v147_counts

# Group by celltype and compute the mean
avg_counts = (
    adata.obs
    .groupby("detailed_celltype_annotation")["v147_counts"]
    .mean()
    .sort_values(ascending=False)
)

# Plot as a bar chart
avg_counts.plot(kind="bar", figsize=(10,5))
plt.ylabel("Average v147 counts per cell")
plt.xlabel("Cell type")
plt.title("Average viral reads (v147) per cell type")
plt.tight_layout()
plt.show()


# In[4]:


indir = "area_cellex_results"
#Loop through each CSV, pull out the v147 row, and build a per‐area DataFrame
records = []
for fn in os.listdir(indir):
    if not fn.endswith("_detailed_cellex.csv"):
        continue
    area = fn.replace("_detailed_cellex.csv","")
    df   = pd.read_csv(os.path.join(indir, fn), index_col=0)
    if "v147" not in df.index:
        continue
    # extract the row as a Series of cell‐type → ESMU
    s = df.loc["v147"].rename(area)
    records.append(s)

#Combine into one DataFrame: rows=areas, cols=cell‐types
combined = pd.DataFrame(records)
#Make sure columns match your actual adata.obs categories
celltypes = sorted(adata.obs["detailed_celltype_annotation"].unique())
combined = combined.reindex(columns=celltypes, fill_value=0)
combined = combined.sort_index()
combined_swap = combined.T
print(combined)


# In[5]:


print(combined.index)   # should match adata.obs['area'] categories
print(combined.columns) # should match adata.obs['detailed_celltype_annotation'] categories


# In[8]:


# build one score per cell: value = table[area, celltype]
def lookup_score(row):
    area = row["area"]
    ct   = row["detailed_celltype_annotation"]
    try:
        return combined.loc[area, ct]
    except KeyError:
        # if some combos are missing, return NaN
        return np.nan

adata.obs["region_celltype_score"] = adata.obs.apply(lookup_score, axis=1)


# In[9]:


sc.pl.embedding(
    adata,
    basis="tsne",                       
    color="region_celltype_score",
    vmin=0, vmax=1,
    cmap="viridis",
    title="Region × celltype score"
)


# In[11]:


for region in combined.index:
    colname = f"score_{region}"
    # map celltype → value for this region
    m = combined.loc[region]
    adata.obs[colname] = adata.obs["detailed_celltype_annotation"].map(m)


# In[12]:


sc.pl.embedding(
    adata,
    basis="tsne",                       
    color="score_BST",
    vmin=0, vmax=1,
    cmap="viridis"
)


# In[15]:


score_col = "region_celltype_score"   # already in adata.obs
area_col  = "area"

# start with 'other' everywhere
adata.obs["area_highlight"] = "no_projection"

# cells that should be colored by their area
mask = adata.obs[score_col] > 0.1

# set those entries to the area name
adata.obs.loc[mask, "area_highlight"] = adata.obs.loc[mask, area_col]

# make it categorical (helps with consistent colors)
adata.obs["area_highlight"] = adata.obs["area_highlight"].astype("category")

# example: get unique areas from obs['area']
areas = [a for a in adata.obs[area_col].unique() if pd.notna(a)]

# pick N distinct colors
cmap = plt.get_cmap("tab10")
area_color_map = {area: cmap(i % 10) for i, area in enumerate(sorted(areas))}

# Add grey for 'other'
area_color_map["no_projection"] = "#bdbdbd"

# build the palette in the correct category order
cats = list(adata.obs["area_highlight"].cat.categories)
adata.uns["area_highlight_colors"] = [area_color_map.get(c, "#bdbdbd") for c in cats]


# In[17]:


sc.pl.embedding(
    adata,
    basis="tsne",                      
    color="area_highlight",            
    size=6,
    title="Area-specific highlight (score > 0.1)",
)


# In[22]:


# 1) Map areas → colors (in your specified order)
area_colors = {
    "BST": "#4D67C8",
    "CEA": "#A7EA51",
    "LHA": "#5DCEAF",
    "PAG": "#702FA0",
    "PBN": "#FFC001",
    "PVH": "#FF8021",
    "PVT": "#5DCBF3",
}
zero_color = "#BDBDBD"  # grey for score == 0

# 2) Build a new column with area name if score>0, else "zero"
score_col = "region_celltype_score"  # the column you already computed
new_col   = "area_score_color_label" # name for the new label column

adata.obs[new_col] = np.where(
    adata.obs[score_col] > 0.1,
    adata.obs["area"],      # label = area (BST, CEA, …)
    "no_projection"                  # label = "zero" if score == 0
)

# Make it categorical with a fixed ordering (areas + zero)
cats = list(area_colors.keys()) + ["no_projection"]
adata.obs[new_col] = adata.obs[new_col].astype(
    pd.CategoricalDtype(categories=cats, ordered=True)
)

# 3) Install the matching color palette in .uns for this column
palette = [area_colors[a] for a in area_colors.keys()] + [zero_color]
adata.uns[f"{new_col}_colors"] = palette


# In[23]:


sc.pl.embedding(
    adata,
    basis="tsne",               
    color=new_col,                  
    size=6,
    title="Area-specific coloring for score > 0.1"
)


# In[ ]:




