#!/usr/bin/env python
# coding: utf-8

# In[9]:


import os as os
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
import anndata as ad
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import seaborn as sns
import adjustText
from pathlib import Path


# In[7]:


adata = sc.read_h5ad("HYP_integrated_2plex_labelled.h5ad")
adata


# In[8]:


adata.layers["RNA"] = adata.X
adata.layers


# In[9]:


print(adata.X[:5, :5].toarray())#Check if counts are the active layer


# In[10]:


adata.var = pd.DataFrame(index=adata.var_names)
adata


# In[11]:


# Remove specific columns from adata.var
columns_to_remove = ["nUMI_hto", "nCount_SCT", "nFeature_SCT", "SCT_snn_res.0.8", "SCT_snn_res.0.8", "integrated_snn_res.0.8"]  
adata.obs = adata.obs.drop(columns=columns_to_remove)
adata.obs


# In[12]:


sc.pl.embedding(adata,basis="X_umap",color="seurat_clusters", legend_loc="on data", legend_fontsize="xx-small")


# In[14]:


ax = sc.pl.embedding(
    adata,
    basis="X_umap",
    color="detailed_celltype_annotation",
    legend_loc="on data",
    legend_fontsize="x-small",
    return_fig=True)
ax=plt.gca()
texts = ax.texts
adjustText.adjust_text(texts, 
                       ax=ax,
                       arrowprops=dict(arrowstyle="-", color='black', lw=0.5))
plt.savefig("testplot.pdf")


# In[15]:


adata.write("HYP_integrated_2plex_anndata.h5ad")
adata


# In[2]:


adata = sc.read_h5ad("HYP_integrated_2plex_anndata.h5ad")
adata


# In[3]:


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


# In[8]:


adata_pp = adata.copy()
from scipy.sparse import issparse
data_mat = adata_pp.X.T
# convert to CSC if possible. See https://github.com/MarioniLab/scran/issues/70
if issparse(data_mat):
    if data_mat.nnz > 2**31 - 1:
        data_mat = data_mat.tocoo()
    else:
        data_mat = data_mat.tocsc()
ro.globalenv["data_mat"] = data_mat
ro.globalenv["input_groups"] = adata_pp.obs["detailed_celltype_annotation"]


# In[9]:


get_ipython().run_cell_magic('R', '-o size_factors', '\nsize_factors = sizeFactors(\n    computeSumFactors(\n        SingleCellExperiment(\n            list(counts=data_mat)), \n            clusters = input_groups,\n            min.mean = 0.1,\n            BPPARAM = MulticoreParam()\n    )\n)\n')


# In[11]:


from scipy.sparse import csr_matrix
adata.obs["size_factors"] = size_factors
scran = adata.X / adata.obs["size_factors"].values[:, None]
scran = scran.tocsc()
adata.layers["scran_normalization"] = csr_matrix(sc.pp.log1p(scran))


# In[13]:


fig, axes = plt.subplots(1, 2, figsize=(10, 5))
p1 = sns.histplot(adata.obs["nCount_RNA"], bins=100, kde=False, ax=axes[0])
axes[0].set_title("Total counts")
p2 = sns.histplot(
    adata.layers["scran_normalization"].sum(1), bins=100, kde=False, ax=axes[1]
)
axes[1].set_title("log1p with Scran estimated size factors")
plt.show()


# In[14]:


adata.write("HYP_integrated_2plex_scran.h5ad")
adata


# In[10]:


# Retrograde (v36) & anterograde (v147) tracing from/to NTS
# Integration with the 7-plex reference atlas

# %% Cell 1 — Imports
import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.patches import Patch
import seaborn as sns
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

OUTDIR = "figures_dvc"
os.makedirs(OUTDIR, exist_ok=True)


# In[11]:


sns.set_style("whitegrid")
plt.rcParams.update({
    'font.size': 9, 'axes.titlesize': 11, 'axes.labelsize': 10,
    'figure.dpi': 150, 'savefig.dpi': 300, 'savefig.bbox': 'tight',
    'pdf.fonttype': 42, 'ps.fonttype': 42,
})

ct_col = 'detailed_celltype_annotation'
score_col = 'score_detailed_celltype_annotation'
lmx_types = ['04_Crnde_Lmx1a', '27_Lmx1a_Ttc6', '32_Gucy2c_Ttc6', '36_Ebf2_Crhbp']

# %% Cell 2 — Load data
adata_2p = sc.read("HYP_integrated_2plex_scran.h5ad")
adata_ref = sc.read("allgenes_bipca_annotated.h5ad")

print(f"2-plex: {adata_2p.n_obs} cells × {adata_2p.n_vars} genes")
print(f"Atlas:  {adata_ref.n_obs} cells × {adata_ref.n_vars} genes")
adata_2p


# In[4]:


# %% Cell 3 — Parse animal IDs and assign tracing groups
adata_2p.obs['animal_id'] = adata_2p.obs['sample_id'].str.extract(r'(\d+)').astype(int).values

# Animal assignments from surgery sheet (Biocenter_May23.xlsx)
# Exclude 188 (failed NTS injection, NTS−) and 147 (no SAC date)
nts_retro_ids  = [168, 176, 149, 159, 165]
nts_antero_ids = [184, 145, 186, 190, 157]
excluded_ids   = [188, 147]

adata_2p.obs['tracing_group'] = 'other'  # VMH/MM anterograde
adata_2p.obs.loc[adata_2p.obs['animal_id'].isin(nts_retro_ids),  'tracing_group'] = 'NTS_retro'
adata_2p.obs.loc[adata_2p.obs['animal_id'].isin(nts_antero_ids), 'tracing_group'] = 'NTS_antero'
adata_2p.obs.loc[adata_2p.obs['animal_id'].isin(excluded_ids),   'tracing_group'] = 'excluded'

print(adata_2p.obs['tracing_group'].value_counts())
print()
print(pd.crosstab(adata_2p.obs['tracing_group'], adata_2p.obs['area']))


# In[5]:


# %% Cell 4 — Subset to NTS retrograde & anterograde
retro  = adata_2p[adata_2p.obs['tracing_group'] == 'NTS_retro'].copy()
antero = adata_2p[adata_2p.obs['tracing_group'] == 'NTS_antero'].copy()

print(f"NTS-R: {retro.n_obs} cells from {retro.obs['animal_id'].nunique()} animals")
print(f"NTS-A: {antero.n_obs} cells from {antero.obs['animal_id'].nunique()} animals")


# In[6]:


# %% Cell 5 — Label transfer quality check
fig, axes = plt.subplots(1, 2, figsize=(12, 4))

for ax, (name, data) in zip(axes, [('NTS retrograde', retro), ('NTS anterograde', antero)]):
    scores = data.obs[score_col]
    ax.hist(scores, bins=50, color='steelblue', edgecolor='white', linewidth=0.3)
    ax.axvline(0.1, color='firebrick', ls='--', lw=1.5, label='threshold 0.1')
    above = (scores > 0.1).mean() * 100
    ax.set_title(f'{name}\nmedian={scores.median():.3f}, >{0.1}: {above:.0f}%')
    ax.set_xlabel('Prediction score')
    ax.legend()

plt.tight_layout()
plt.savefig(f'{OUTDIR}/QC_label_transfer_scores.pdf')
plt.show()


# In[7]:


# %% Cell 6 — Filter low-confidence labels (threshold 0.4)
retro_hq  = retro[retro.obs[score_col] > 0.1].copy()
antero_hq = antero[antero.obs[score_col] > 0.1].copy()

print(f"NTS-R after filtering: {retro_hq.n_obs}/{retro.n_obs} ({retro_hq.n_obs/retro.n_obs*100:.1f}%)")
print(f"NTS-A after filtering: {antero_hq.n_obs}/{antero.n_obs} ({antero_hq.n_obs/antero.n_obs*100:.1f}%)")


# In[8]:


# %% Cell 7 — Cell type enrichment: NTS-R vs atlas (Fisher's exact)
ref_cts   = adata_ref.obs[ct_col].value_counts()
retro_cts = retro_hq.obs[ct_col].value_counts()
all_types = sorted(set(ref_cts.index) | set(retro_cts.index))

n_retro = retro_hq.n_obs
n_ref   = adata_ref.n_obs

records = []
for ct in all_types:
    n_ct_r = retro_cts.get(ct, 0)
    n_ct_a = ref_cts.get(ct, 0)

    pct_r = n_ct_r / n_retro * 100
    pct_a = n_ct_a / n_ref * 100
    fold  = (pct_r / pct_a) if pct_a > 0 else np.inf

    table = np.array([[n_ct_r, n_retro - n_ct_r],
                      [n_ct_a, n_ref   - n_ct_a]])
    odds, pval = fisher_exact(table, alternative='greater')

    records.append({
        'cell_type': ct, 'n_retro': n_ct_r, 'pct_retro': pct_r,
        'n_atlas': n_ct_a, 'pct_atlas': pct_a,
        'fold_enrichment': fold, 'odds_ratio': odds, 'pvalue': pval,
        'is_lmx': ct in lmx_types,
    })

df_enrich = pd.DataFrame(records)
_, fdr, _, _ = multipletests(df_enrich['pvalue'], method='fdr_bh')
df_enrich['fdr'] = fdr
df_enrich['significant'] = df_enrich['fdr'] < 0.05
df_enrich = df_enrich.sort_values('fold_enrichment', ascending=False)

df_enrich.to_csv(f'{OUTDIR}/DVC_enrichment_table.csv', index=False)
print("Saved DVC_enrichment_table.csv")


# In[9]:


# %% Cell 8 — Print enrichment results
sig = df_enrich[(df_enrich['significant']) & (df_enrich['n_retro'] >= 5)]
print(f"Significantly enriched in DVC-projecting fraction (FDR<0.05, n≥5):\n")
print(f"{'Cell type':30s} {'n_retro':>7s} {'%_retro':>7s} {'%_atlas':>7s} {'fold':>6s} {'FDR':>10s}")
print("-" * 75)
for _, r in sig.iterrows():
    tag = " ★" if r['is_lmx'] else ""
    print(f"{r['cell_type']:30s} {r['n_retro']:>7.0f} {r['pct_retro']:>6.1f}% "
          f"{r['pct_atlas']:>6.1f}% {r['fold_enrichment']:>6.2f} {r['fdr']:>10.2e}{tag}")

# Total Lmx1a enrichment
lmx_r = retro_hq.obs[ct_col].isin(lmx_types).sum()
lmx_a = adata_ref.obs[ct_col].isin(lmx_types).sum()
print(f"\nTotal Lmx1a: {lmx_r}/{retro_hq.n_obs} ({lmx_r/retro_hq.n_obs*100:.1f}%) "
      f"vs atlas {lmx_a}/{adata_ref.n_obs} ({lmx_a/adata_ref.n_obs*100:.1f}%) "
      f"→ {(lmx_r/retro_hq.n_obs)/(lmx_a/adata_ref.n_obs):.2f}× enriched")


# In[10]:


# %% Cell 9 — Per-animal breakdown
ct_by_animal = pd.crosstab(retro_hq.obs['animal_id'], retro_hq.obs[ct_col])
ct_by_animal_pct = ct_by_animal.div(ct_by_animal.sum(axis=1), axis=0) * 100

lmx_cols = [c for c in lmx_types if c in ct_by_animal_pct.columns]

print("Per-animal Lmx1a representation:\n")
for animal in ct_by_animal_pct.index:
    lmx_pct = ct_by_animal_pct.loc[animal, lmx_cols].sum()
    n_cells = ct_by_animal.loc[animal].sum()
    sex = retro_hq.obs[retro_hq.obs['animal_id'] == animal]['sex'].iloc[0]
    print(f"  Animal {animal} ({sex}): {lmx_pct:.1f}% Lmx1a+ ({n_cells} cells total)")


# In[11]:


# %% Cell 10 — Figure: Enrichment bar chart
fig, ax = plt.subplots(figsize=(6, 8))

df_plot = df_enrich.nlargest(25, 'n_retro').sort_values('pct_retro', ascending=True)

colors = []
for _, r in df_plot.iterrows():
    if r['is_lmx'] and r['significant']:   colors.append('#8B0000')
    elif r['significant']:                  colors.append('#E8A87C')
    elif r['is_lmx']:                       colors.append('#C44E52')
    else:                                   colors.append('#C0C0C0')

y = range(len(df_plot))
ax.barh(y, df_plot['pct_retro'].values, color=colors,
        edgecolor='white', linewidth=0.5, height=0.7)

# Atlas proportion as black ticks
for i, (_, r) in enumerate(df_plot.iterrows()):
    ax.plot(r['pct_atlas'], i, 'k|', markersize=8, markeredgewidth=1.5, zorder=5)

# Fold annotation for significant enriched types
for i, (_, r) in enumerate(df_plot.iterrows()):
    if r['significant'] and r['fold_enrichment'] > 1:
        ax.text(r['pct_retro'] + 0.3, i, f'{r["fold_enrichment"]:.1f}×',
                va='center', fontsize=7, fontweight='bold',
                color='#8B0000' if r['is_lmx'] else '#666')

ax.set_yticks(y)
ax.set_yticklabels(df_plot['cell_type'].values, fontsize=7)
ax.set_xlabel('% of sorted GFP+ nuclei')
ax.set_title('Cell type composition of DVC-projecting\nhypothalamic neurons (NTS retrograde)',
             fontweight='bold')

legend_el = [
    Patch(fc='#8B0000', label='Lmx1a + FDR<0.05'),
    Patch(fc='#C44E52', label='Lmx1a (n.s.)'),
    Patch(fc='#E8A87C', label='FDR<0.05'),
    Patch(fc='#C0C0C0', label='Not enriched'),
    plt.Line2D([0], [0], color='k', marker='|', ls='None',
               markersize=8, markeredgewidth=1.5, label='Atlas proportion'),
]
ax.legend(handles=legend_el, fontsize=7, loc='lower right')
plt.tight_layout()
plt.savefig(f'{OUTDIR}/Fig_DVC_enrichment_bar.pdf')
plt.show()


# In[5]:


# %% Cell 10 (revised) — Enrichment bar chart, neurons only
OUTDIR = "figures_dvc"
df_enrich = pd.read_csv(f'{OUTDIR}/DVC_enrichment_table.csv')

# Cell types to exclude
glial_types = {
    'Astrocyte', 'Oligodendrocyte', 'OPC', 'Microglia',
    'VLMC', 'Ependymal_cell', 'Tanycyte', 'Pericyte',
    'Schwann_cell', 'Pituicyte', 'Epithelial-like_cell',
    'Endothelial_cell',
}
pituitary_types = {c for c in df_enrich['cell_type'] if 'Pituitary' in c}
exclude = glial_types | pituitary_types | {'01_Undefined', 'Undefined'}

df_enrich_neurons = df_enrich[~df_enrich['cell_type'].isin(exclude)].copy()

df_plot = df_enrich_neurons.nlargest(25, 'n_retro').sort_values('pct_retro', ascending=True)

# ── Figure ────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(6, 8))

colors = []
for _, r in df_plot.iterrows():
    if r['is_lmx'] and r['significant']:  colors.append('#8B0000')
    elif r['significant']:                 colors.append('#E8A87C')
    elif r['is_lmx']:                      colors.append('#C44E52')
    else:                                  colors.append('#C0C0C0')

y = range(len(df_plot))
ax.barh(y, df_plot['pct_retro'].values, color=colors,
        edgecolor='white', linewidth=0.5, height=0.7)

for i, (_, r) in enumerate(df_plot.iterrows()):
    ax.plot(r['pct_atlas'], i, 'k|', markersize=8, markeredgewidth=1.5, zorder=5)

for i, (_, r) in enumerate(df_plot.iterrows()):
    if r['significant'] and r['fold_enrichment'] > 1:
        ax.text(r['pct_retro'] + 0.3, i, f'{r["fold_enrichment"]:.1f}×',
                va='center', fontsize=7, fontweight='bold',
                color='#8B0000' if r['is_lmx'] else '#666')

ax.set_yticks(y)
ax.set_yticklabels(df_plot['cell_type'].values, fontsize=7)
ax.set_xlabel('% of sorted GFP+ nuclei')
ax.set_title('Cell type composition of DVC-projecting\nhypothalamic neurons (NTS retrograde)',
             fontweight='bold')

legend_el = [
    Patch(fc='#8B0000', label='Lmx1a + FDR<0.05'),
    Patch(fc='#C44E52', label='Lmx1a (n.s.)'),
    Patch(fc='#E8A87C', label='FDR<0.05'),
    Patch(fc='#C0C0C0', label='Not enriched'),
    plt.Line2D([0], [0], color='k', marker='|', ls='None',
               markersize=8, markeredgewidth=1.5, label='Atlas proportion'),
]
ax.legend(handles=legend_el, fontsize=7, loc='lower right')
plt.tight_layout()
plt.savefig(f'{OUTDIR}/Fig_DVC_enrichment_bar_neurons.pdf')
plt.show()


# In[12]:


# %% Cell 11 — Figure: Per-animal stacked bar
fig, ax = plt.subplots(figsize=(8, 4))

top_types = df_enrich.nlargest(10, 'n_retro')['cell_type'].tolist()
for lt in lmx_types:
    if lt not in top_types and lt in retro_hq.obs[ct_col].values:
        top_types.append(lt)

animal_comp = pd.crosstab(retro_hq.obs['animal_id'], retro_hq.obs[ct_col])
other_cols = [c for c in animal_comp.columns if c not in top_types]
if other_cols:
    animal_comp['Other'] = animal_comp[other_cols].sum(axis=1)
    animal_comp = animal_comp.drop(columns=other_cols)

animal_pct = animal_comp.div(animal_comp.sum(axis=1), axis=0) * 100

lmx_colors = {'04_Crnde_Lmx1a': '#8B0000', '27_Lmx1a_Ttc6': '#C44E52',
               '32_Gucy2c_Ttc6': '#E8A87C', '36_Ebf2_Crhbp': '#EDCBB1'}
palette = plt.cm.Set3(np.linspace(0, 1, 15))
ct_colors, ci = {}, 0
for ct in animal_pct.columns:
    if ct in lmx_colors:        ct_colors[ct] = lmx_colors[ct]
    elif ct == 'Other':          ct_colors[ct] = '#D0D0D0'
    else:                        ct_colors[ct] = palette[ci]; ci += 1

animal_pct.plot(kind='bar', stacked=True, ax=ax,
                color=[ct_colors.get(c, '#999') for c in animal_pct.columns],
                edgecolor='white', linewidth=0.3)

for i, animal in enumerate(animal_pct.index):
    sex = retro_hq.obs[retro_hq.obs['animal_id'] == animal]['sex'].iloc[0]
    ax.text(i, -3, f'({sex[0].upper()})', ha='center', fontsize=8, color='grey')

ax.set_ylabel('% of GFP+ nuclei')
ax.set_xlabel('Animal ID')
ax.set_title('Cell type composition per animal\n(NTS retrograde, GFP+ sorted)', fontweight='bold')
ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=7, ncol=1)
ax.set_xticklabels(animal_pct.index, rotation=0)
plt.tight_layout()
plt.savefig(f'{OUTDIR}/Fig_DVC_per_animal.pdf')
plt.show()


# In[13]:


# %% Cell 13 — Compute DVC row for projection heatmap
# Score = fold_enrichment / 5, capped at 1.0, zeroed if fold < 1
# This maps 5× enrichment → 1.0 on the ESµ scale

dvc_scores = {}
for _, r in df_enrich.iterrows():
    ct = r['cell_type']
    if r['n_retro'] >= 3 and r['pct_atlas'] > 0 and r['fold_enrichment'] > 1:
        dvc_scores[ct] = min(r['fold_enrichment'] / 5.0, 1.0)
    else:
        dvc_scores[ct] = 0.0

dvc_series = pd.Series(dvc_scores, name='DVC (NTS)')

print("Top DVC projection scores:\n")
for ct, score in dvc_series.nlargest(15).items():
    tag = " ★" if ct in lmx_types else ""
    print(f"  {ct:30s}  {score:.3f}{tag}")


# In[14]:


# %% Cell 14 — Build and save DVC heatmap row aligned to dendrogram
dendro_info = adata_ref.uns['dendrogram_detailed_celltype_annotation']
cell_type_order = dendro_info['categories_ordered']

dvc_heatmap_row = pd.DataFrame(
    [[dvc_series.get(ct, np.nan) for ct in cell_type_order]],
    index=['DVC (NTS)'],
    columns=cell_type_order
)

dvc_heatmap_row.to_csv(f'{OUTDIR}/DVC_heatmap_row.csv')
print(f"Saved DVC_heatmap_row.csv ({len(cell_type_order)} columns)")


# In[15]:


# %% Cell 15 — Plot standalone DVC heatmap row
fig, ax = plt.subplots(figsize=(18, 1.5))

cmap = plt.get_cmap('rocket_r')
cmap.set_bad(color='lightgray')
norm = Normalize(vmin=0, vmax=1)
arr = dvc_heatmap_row.values.astype(float)

ax.imshow(np.ma.masked_invalid(arr), aspect='auto', cmap=cmap,
          norm=norm, interpolation='nearest')

ax.set_xticks(np.arange(-.5, len(cell_type_order), 1), minor=True)
ax.grid(which='minor', color='gray', ls='-', lw=0.5, alpha=0.3)

for col_i, ct in enumerate(cell_type_order):
    val = arr[0, col_i]
    if not np.isnan(val) and val > 0.36:
        c = 'white' if norm(val) > 0.5 else 'black'
        ax.text(col_i, 0, '*', ha='center', va='center',
                fontsize=10, color=c, fontweight='bold')

ax.set_yticks([0])
ax.set_yticklabels(['DVC (NTS)'], fontsize=10, fontweight='bold')
ax.set_xticks(range(len(cell_type_order)))
ax.set_xticklabels(cell_type_order, rotation=90, ha='right', fontsize=7)

plt.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=ax,
             orientation='vertical', fraction=0.02, pad=0.02,
             label='Enrichment score')
plt.tight_layout()
plt.savefig(f'{OUTDIR}/Fig_DVC_heatmap_row.pdf')
plt.show()


# In[16]:


# Use only highly variable genes
adata_hvg = adata_ref[:, adata_ref.var['highly_deviant']].copy()

sc.tl.dendrogram(
    adata_hvg, 
    groupby="detailed_celltype_annotation"
)
# Get the dendrogram ordering
dendro_info = adata_hvg.uns['dendrogram_detailed_celltype_annotation']
cell_type_order = dendro_info['categories_ordered']
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
celltypes = sorted(adata_ref.obs["detailed_celltype_annotation"].unique())
combined = combined.reindex(columns=celltypes, fill_value=0)
combined = combined.sort_index()
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
celltypes = sorted(adata_ref.obs["detailed_celltype_annotation"].unique())
combined_spatial = combined_spatial.reindex(columns=celltypes, fill_value=0)
combined_spatial.to_csv('df_spatial.csv', index=False)
score_df = pd.DataFrame(index=combined.index, columns=cell_type_order, dtype=float)


# In[17]:


# %% Cell 16 — Integrate DVC row into existing combined heatmap

dvc_row = pd.read_csv('figures_dvc/DVC_heatmap_row.csv', index_col=0)
dvc_aligned = dvc_row.reindex(columns=score_df.columns, fill_value=np.nan)
score_df_ext = pd.concat([score_df, dvc_aligned])


# In[18]:

single_dataset_regions = {'PVH'}  # only snRNA available, do not divide by 2

for region in combined.index:
    for celltype in cell_type_order:
        snrna_val   = combined.loc[region, celltype] \
                      if celltype in combined.columns else np.nan
        spatial_val = combined_spatial.loc[region, celltype] \
                      if (region in combined_spatial.index
                          and celltype in combined_spatial.columns) else np.nan

        if np.isnan(snrna_val) and np.isnan(spatial_val):
            score_df.loc[region, celltype] = np.nan
        elif region in single_dataset_regions:
            # One dataset only — raw ESµ, scale already 0–1
            score_df.loc[region, celltype] = np.nansum([snrna_val, spatial_val])
        else:
            # Two datasets — sum then divide by 2, keeping scale 0–1
            score_df.loc[region, celltype] = np.nansum([snrna_val, spatial_val]) / 2

# Convert to float array for imshow
score_array = score_df.values.astype(float)

# Normalize to 0-1 linear scale
norm = Normalize(vmin=0, vmax=1)
cmap = plt.get_cmap('rocket_r')
cmap.set_bad(color='lightgray')

# Create figure
fig, ax = plt.subplots(figsize=(18, 4))

# Plot the combined score grid
im = ax.imshow(np.ma.masked_invalid(score_array),
               aspect='auto', cmap=cmap, norm=norm, interpolation='nearest')

# Add faint grid
ax.set_xticks(np.arange(-.5, len(cell_type_order), 1), minor=True)
ax.set_yticks(np.arange(-.5, len(score_df.index), 1), minor=True)
ax.grid(which='minor', color='gray', linestyle='-', linewidth=0.5, alpha=0.3)

# Asterisks: threshold > 0.36 across all regions
for row_i, region in enumerate(score_df.index):
    for col_i in range(len(cell_type_order)):
        val = score_array[row_i, col_i]
        if np.isnan(val):
            continue
        if val > 0.36:
            brightness = norm(val)
            text_color = 'white' if brightness > 0.5 else 'black'
            ax.text(col_i, row_i, '*', ha='center', va='center',
                    fontsize=10, color=text_color, fontweight='bold')

# Set labels
ax.set_yticks(range(len(score_df.index)))
ax.set_yticklabels(score_df.index, fontsize=12)
ax.set_xticks(range(len(cell_type_order)))
ax.set_xticklabels(cell_type_order, rotation=90, ha='right', fontsize=8)

# Add colorbar
cbar = plt.colorbar(
    ScalarMappable(norm=norm, cmap=cmap),
    ax=ax,
    orientation='vertical',
    fraction=0.02,
    pad=0.02
)
cbar.set_label('Combined projection score\n(mean ESµ)', fontsize=10)
cbar.set_ticks([0, 0.25, 0.5, 0.75, 1.0])
cbar.ax.tick_params(labelsize=9)

plt.tight_layout()
plt.savefig('region_celltype_annotation_combined_score_DVC.pdf', bbox_inches='tight', dpi=300)
plt.show()


# In[19]:

# %% Cell 17 — Anterograde composition (supplementary)
fig, ax = plt.subplots(figsize=(6, 6))

antero_cts = antero_hq.obs[ct_col].value_counts()
antero_pct = (antero_cts / antero_cts.sum() * 100).head(20)

antero_pct[::-1].plot(kind='barh', ax=ax, color='#4878CF',
                       edgecolor='white', linewidth=0.5)

for i, (ct, pct) in enumerate(antero_pct[::-1].items()):
    ax.text(pct + 0.2, i, f'n={antero_cts[ct]}', va='center', fontsize=7)

ax.set_xlabel('% of sorted GFP+ nuclei')
ax.set_title('Cell types receiving DVC→HYP input\n(NTS anterograde)', fontweight='bold')
plt.tight_layout()
plt.savefig(f'{OUTDIR}/Fig_DVC_anterograde.pdf')
plt.show()


# In[20]:


# %% Cell 18 — Summary
print("=" * 70)
print("DVC TRACING ANALYSIS SUMMARY")
print("=" * 70)

lmx_r = retro_hq.obs[ct_col].isin(lmx_types).sum()
lmx_a = adata_ref.obs[ct_col].isin(lmx_types).sum()

print(f"""
Animals used (NTS retrograde): {nts_retro_ids}
  Excluded: 188 (failed injection), 147 (no SAC)

Cells:
  NTS-R total:     {retro.n_obs}
  NTS-R (score>0.4): {retro_hq.n_obs}
  NTS-A total:     {antero.n_obs}
  NTS-A (score>0.4): {antero_hq.n_obs}

Lmx1a enrichment in DVC-projecting fraction:
  DVC retro: {lmx_r}/{retro_hq.n_obs} ({lmx_r/retro_hq.n_obs*100:.1f}%)
  Atlas:     {lmx_a}/{adata_ref.n_obs} ({lmx_a/adata_ref.n_obs*100:.1f}%)
  Fold:      {(lmx_r/retro_hq.n_obs)/(lmx_a/adata_ref.n_obs):.2f}×

Enriched cell types (FDR<0.05, n≥5): {len(sig)}
""")

for _, r in sig.iterrows():
    tag = " ★ Lmx1a" if r['is_lmx'] else ""
    print(f"  {r['cell_type']:30s} {r['fold_enrichment']:.1f}× (FDR={r['fdr']:.1e}){tag}")

print(f"\nAll outputs in: {OUTDIR}/")


# In[21]:


print(adata_2p.X[:5, :5].toarray())


# In[22]:


print(adata_ref.X[:5, :5].toarray())


# In[23]:


# ── Reference atlas: v147 expression at animal × lane level (scran-normalized) ─

def get_expr_scran(adata, gene):
    """Extract scran-normalized expression, falling back to .X if layer absent."""
    if gene not in adata.var_names:
        print(f"  WARNING: {gene} not found"); return None
    layer = adata[:, gene].layers['scran_normalization']
    x = layer.toarray().flatten() if hasattr(layer, 'toarray') else layer.flatten()
    return x

ref_v147  = get_expr_scran(adata_ref, 'v147')
ref_sort  = adata_ref.obs['sort'].values.astype(float)

ref_records = []
for (animal, lane), grp in adata_ref.obs.groupby(['hash.mcl.ID', 'orig.ident'],
                                                   observed=True):
    mask = adata_ref.obs_names.isin(grp.index)
    e    = ref_v147[mask]
    s    = ref_sort[mask].mean()

    ref_records.append({
        'animal':    animal,
        'lane':      lane,
        'sort_prob': s,
        'pct_pos':   (e > 0).mean() * 100,   # positivity unchanged by normalization
        'mean_expr': e.mean(),                 # now on scran scale
        'n_cells':   len(e),
    })

df_ref = pd.DataFrame(ref_records).sort_values('sort_prob')

bins   = [0, 0.05, 0.2, 0.5, 0.85, 0.99, 1.01]
labels = ['<5%', '5–20%', '20–50%', '50–85%', '85–99%', '100%']
df_ref['sort_bin'] = pd.cut(df_ref['sort_prob'], bins=bins, labels=labels, right=False)

print("Animal×lane combinations per sort bin (scran-normalized):")
print(df_ref.groupby('sort_bin', observed=True)['pct_pos']
           .agg(['median', 'mean', 'std', 'count']).round(2))
print(f"\nAnimal×lane combos with sort=1.0: "
      f"{(df_ref['sort_prob'] == 1.0).sum()} / {len(df_ref)}")

# ── Extract scran-normalized tracer expression from 2-plex ───────────────────

def get_expr_scran_2plex(adata, gene):
    """Try scran layer first, fall back to .X (raw counts) if absent."""
    if gene not in adata.var_names:
        print(f"  WARNING: {gene} not found"); return None
    if 'scran_normalization' in adata.layers:
        layer = adata[:, gene].layers['scran_normalization']
    else:
        print(f"  NOTE: scran_normalization layer not found in 2-plex, using .X")
        layer = adata[:, gene].X
    return layer.toarray().flatten() if hasattr(layer, 'toarray') else layer.flatten()

retro_expr  = get_expr_scran_2plex(retro_hq,  'v36')
antero_expr = get_expr_scran_2plex(antero_hq, 'v147')

print(f"NTS-R  v36  : {(retro_expr  > 0).mean()*100:.2f}% positive, "
      f"n={len(retro_expr)}")
print(f"NTS-A  v147 : {(antero_expr > 0).mean()*100:.2f}% positive, "
      f"n={len(antero_expr)}")
# ── Figure ────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(13, 5))

ax = axes[0]
sc = ax.scatter(df_ref['sort_prob'], df_ref['pct_pos'],
                c=df_ref['sort_prob'], cmap='Blues',
                s=30, alpha=0.7, edgecolors='k', linewidths=0.3, zorder=3,
                label='Atlas animal×lane')

# 2-plex positivity is still count>0, consistent with scran (zeros are preserved)
retro_pct  = (retro_expr  > 0).mean() * 100
antero_pct = (antero_expr > 0).mean() * 100

ax.scatter(1.0, retro_pct,  marker='*', s=300, color='firebrick', zorder=6,
           label=f'2-plex NTS-R / v36 ({retro_pct:.1f}%)')
ax.scatter(1.0, antero_pct, marker='*', s=300, color='darkorange', zorder=6,
           label=f'2-plex NTS-A / v147 ({antero_pct:.1f}%)')

z    = np.polyfit(df_ref['sort_prob'], df_ref['pct_pos'], 1)
xfit = np.linspace(0, 1, 200)
ax.plot(xfit, np.poly1d(z)(xfit), 'k--', lw=1.2, alpha=0.5,
        label=f'Atlas fit (slope={z[0]:.1f}%/unit sort)')

plt.colorbar(sc, ax=ax, label='Sort probability')
ax.set_xlabel('Sort probability (animal × lane)')
ax.set_ylabel('% tracer-positive nuclei')
ax.set_title('v147 atlas (scran-normalized) vs 2-plex\nanimal × lane resolution',
             fontweight='bold')
ax.legend(fontsize=8)
ax.set_xlim(-0.02, 1.08); ax.set_ylim(bottom=0)
ax.grid(alpha=0.3)

ax = axes[1]
order = [l for l in labels if l in df_ref['sort_bin'].cat.categories]
sns.boxplot(data=df_ref, x='sort_bin', y='pct_pos', order=order,
            palette='Blues', ax=ax, width=0.5, linewidth=0.8,
            flierprops=dict(marker='o', ms=4, alpha=0.5))
sns.stripplot(data=df_ref, x='sort_bin', y='pct_pos', order=order,
              color='k', size=3, alpha=0.4, ax=ax)

ax.axhline(retro_pct,  color='firebrick',  ls='--', lw=1.8,
           label=f'2-plex NTS-R / v36 ({retro_pct:.1f}%)')
ax.axhline(antero_pct, color='darkorange', ls=':',  lw=1.8,
           label=f'2-plex NTS-A / v147 ({antero_pct:.1f}%)')

if '100%' in df_ref['sort_bin'].cat.categories:
    med_100 = df_ref.loc[df_ref['sort_bin'] == '100%', 'pct_pos'].median()
    ax.text(order.index('100%'), med_100 + 1,
            f'median\n{med_100:.1f}%', ha='center', fontsize=8, color='steelblue')

ax.set_xlabel('Sort probability bin (atlas animal × lane)')
ax.set_ylabel('% tracer-positive nuclei')
ax.set_title('Atlas v147-retro by sort bin\n(scran-normalized, animal × lane)',
             fontweight='bold')
ax.legend(fontsize=8)
ax.grid(axis='y', alpha=0.3)

plt.suptitle('Tracer mRNA detection: 2-plex GFP-sorted cells vs atlas 100%-sorted lanes\n'
             '(scran-normalized counts)',
             fontweight='bold', fontsize=11)
plt.tight_layout()
plt.savefig(f'{OUTDIR}/Fig_tracer_expression_validation_scran.pdf', bbox_inches='tight', dpi=300)
plt.show()

# ── Direct comparison ─────────────────────────────────────────────────────────
print("\nDIRECT COMPARISON: 2-plex vs atlas sort=1.0 (scran-normalized)")
print("="*60)
atlas_100 = df_ref[df_ref['sort_prob'] == 1.0]['pct_pos']
print(f"Atlas sort=1.0 lanes : n={len(atlas_100)}, "
      f"median={atlas_100.median():.1f}%, "
      f"range={atlas_100.min():.1f}–{atlas_100.max():.1f}%")
print(f"2-plex NTS-R (v36)   : {retro_pct:.1f}%")
print(f"2-plex NTS-A (v147)  : {antero_pct:.2f}%")


# In[24]:


# %% Cell 20 — % tracer-positive nuclei by cell type
# Atlas: sort=1.0 animal×lane combos, v147-retro
# 2-plex NTS-R: v36  | 2-plex NTS-A: v147

# ── 1. Identify atlas sort=1.0 cells ─────────────────────────────────────────
atlas_100_mask = np.isin(
    adata_ref.obs_names,
    df_ref.loc[df_ref['sort_prob'] == 1.0]
           .apply(lambda r: adata_ref.obs[
               (adata_ref.obs['hash.mcl.ID'] == r['animal']) &
               (adata_ref.obs['orig.ident']  == r['lane'])
           ].index, axis=1)
           .explode().values
)
adata_atlas_100 = adata_ref[atlas_100_mask].copy()
print(f"Atlas sort=1.0 cells: {adata_atlas_100.n_obs} "
      f"from {adata_atlas_100.obs['hash.mcl.ID'].nunique()} animals, "
      f"{adata_atlas_100.obs['orig.ident'].nunique()} lanes")

# ── 2. Compute % positive per cell type for each dataset ─────────────────────
def pct_positive_by_celltype(adata, gene, ct_col='detailed_celltype_annotation',
                              layer=None, min_cells=10):
    """Return DataFrame: cell_type | n_cells | n_pos | pct_pos"""
    if gene not in adata.var_names:
        return pd.DataFrame()
    x = adata[:, gene].layers[layer] if layer and layer in adata.layers \
        else adata[:, gene].X
    expr = x.toarray().flatten() if hasattr(x, 'toarray') else np.asarray(x).flatten()

    records = []
    for ct, grp in adata.obs.groupby(ct_col):
        idx  = adata.obs_names.isin(grp.index)
        e    = expr[idx]
        n    = len(e)
        if n < min_cells:
            continue
        n_pos = (e > 0).sum()
        records.append({'cell_type': ct, 'n_cells': n,
                        'n_pos': n_pos, 'pct_pos': n_pos / n * 100})
    return pd.DataFrame(records).set_index('cell_type')

df_atlas = pct_positive_by_celltype(adata_atlas_100, 'v147',    layer='counts')
df_retro = pct_positive_by_celltype(retro_hq,        'v36',     ct_col=ct_col)
df_antero= pct_positive_by_celltype(antero_hq,       'v147',    ct_col=ct_col)

# Align to shared cell types
all_cts = sorted(set(df_atlas.index) | set(df_retro.index) | set(df_antero.index))

df_combined = pd.DataFrame({
    'Atlas sort=1.0 (v147-retro)': df_atlas['pct_pos'].reindex(all_cts),
    '2-plex NTS-R (v36-retro)':    df_retro['pct_pos'].reindex(all_cts),
    '2-plex NTS-A (v147-antero)':  df_antero['pct_pos'].reindex(all_cts),
})
df_n = pd.DataFrame({
    'n_atlas':  df_atlas['n_cells'].reindex(all_cts, fill_value=0),
    'n_retro':  df_retro['n_cells'].reindex(all_cts, fill_value=0),
    'n_antero': df_antero['n_cells'].reindex(all_cts, fill_value=0),
})

# Sort rows by atlas % positive (descending), put NaN at bottom
sort_key = df_combined['Atlas sort=1.0 (v147-retro)'].fillna(-1)
row_order = sort_key.sort_values(ascending=False).index
df_combined = df_combined.loc[row_order]
df_n        = df_n.loc[row_order]

print(f"\nCell types represented: {len(df_combined)}")
print(df_combined.round(2).to_string())

# ── 3. Figure: heatmap ────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(16, max(8, len(all_cts) * 0.28)),
                          gridspec_kw={'width_ratios': [3, 1, 1]})

cmap = plt.get_cmap('rocket_r')
cmap.set_bad('lightgray')

# Heatmap
arr = df_combined.values.astype(float)
im  = axes[0].imshow(np.ma.masked_invalid(arr), aspect='auto',
                      cmap=cmap, vmin=0, vmax=arr[~np.isnan(arr)].max(),
                      interpolation='nearest')

axes[0].set_xticks(range(3))
axes[0].set_xticklabels(df_combined.columns, rotation=25, ha='right', fontsize=9)
axes[0].set_yticks(range(len(row_order)))
axes[0].set_yticklabels(row_order, fontsize=7)
axes[0].set_title('% tracer-positive nuclei\nby cell type', fontweight='bold')

# Annotate cells
for ri in range(arr.shape[0]):
    for ci in range(arr.shape[1]):
        v = arr[ri, ci]
        if not np.isnan(v) and v > 0:
            tc = 'white' if v > arr[~np.isnan(arr)].max() * 0.6 else 'black'
            axes[0].text(ci, ri, f'{v:.1f}', ha='center', va='center',
                         fontsize=5.5, color=tc)

plt.colorbar(im, ax=axes[0], fraction=0.02, pad=0.01, label='% positive')

# N cells bar — atlas
axes[1].barh(range(len(row_order)), df_n['n_atlas'],
             color='steelblue', alpha=0.7, height=0.6)
axes[1].set_yticks([]); axes[1].set_xlabel('N cells')
axes[1].set_title('N cells\n(atlas)', fontsize=9, fontweight='bold')
axes[1].grid(axis='x', alpha=0.3)

# N cells bar — 2-plex retro
axes[2].barh(range(len(row_order)), df_n['n_retro'],
             color='firebrick', alpha=0.6, height=0.6, label='NTS-R')
axes[2].barh(range(len(row_order)), df_n['n_antero'],
             color='darkorange', alpha=0.6, height=0.6,
             left=df_n['n_retro'], label='NTS-A')
axes[2].set_yticks([]); axes[2].set_xlabel('N cells')
axes[2].set_title('N cells\n(2-plex)', fontsize=9, fontweight='bold')
axes[2].legend(fontsize=7); axes[2].grid(axis='x', alpha=0.3)

plt.suptitle('Tracer gene detection rate by cell type\n'
             'Atlas sort=1.0 lanes vs 2-plex GFP-sorted',
             fontweight='bold', fontsize=11)
plt.tight_layout()
plt.savefig(f'{OUTDIR}/Fig_tracer_pct_by_celltype.pdf', bbox_inches='tight', dpi=300)
plt.show()

# ── 4. Correlation: do atlas and 2-plex retro agree on which CTs detect well? ─
shared = df_combined[['Atlas sort=1.0 (v147-retro)',
                       '2-plex NTS-R (v36-retro)']].dropna()
if len(shared) >= 5:
    from scipy.stats import spearmanr
    r, p = spearmanr(shared.iloc[:, 0], shared.iloc[:, 1])
    print(f"\nSpearman correlation (atlas v147-retro vs 2-plex v36-retro "
          f"% positive per cell type): ρ={r:.3f}, p={p:.3e}, n={len(shared)} cell types")
    print("→ High ρ = cell-type-level detection is a property of the cell type "
          "(e.g. transcriptional activity, nuclear size), not the specific AAV")


# In[25]:


# %% Cell 21 — Stratify 2-plex by tracer expression level using glia as background

# ── 1. Define neuron vs glia labels ──────────────────────────────────────────
glial_types = {
    'Astrocyte', 'Oligodendrocyte', 'OPC', 'Microglia',
    'VLMC', 'Ependymal_cell', 'Tanycyte', 'Pericyte',
    'Schwann_cell', 'Pituicyte', 'Epithelial-like_cell',
    'Endothelial_cell', 'VLMC',
}
pituitary_types = {c for c in retro_hq.obs[ct_col].unique() if 'Pituitary' in c}
non_neuronal = glial_types | pituitary_types

def add_class(adata):
    adata.obs['cell_class'] = 'Neuron'
    adata.obs.loc[adata.obs[ct_col].isin(glial_types),     'cell_class'] = 'Glia'
    adata.obs.loc[adata.obs[ct_col].isin(pituitary_types), 'cell_class'] = 'Pituitary'

add_class(retro_hq)
add_class(antero_hq)

print("NTS-R cell class breakdown:")
print(retro_hq.obs['cell_class'].value_counts())
print("\nNTS-A cell class breakdown:")
print(antero_hq.obs['cell_class'].value_counts())

# ── 2. Extract expression ─────────────────────────────────────────────────────
retro_expr  = get_expr_scran(retro_hq,  'v36')
antero_expr = get_expr_scran(antero_hq, 'v147')

retro_hq.obs['tracer_expr']  = retro_expr
antero_hq.obs['tracer_expr'] = antero_expr

# ── 3. Derive background threshold from glia ──────────────────────────────────
# Use the 95th percentile of glial expression as a conservative background ceiling

for name, adata, expr, gene in [
    ('NTS-R', retro_hq,  retro_expr,  'v36'),
    ('NTS-A', antero_hq, antero_expr, 'v147'),
]:
    glia_mask   = adata.obs['cell_class'] == 'Glia'
    neuron_mask = adata.obs['cell_class'] == 'Neuron'
    glia_expr   = expr[glia_mask.values]
    neuron_expr = expr[neuron_mask.values]

    # Threshold: 95th percentile of glial expression among glia that are >0
    # (avoids zeros dragging threshold to 0)
    glia_pos = glia_expr[glia_expr > 0]
    thresh_95 = np.percentile(glia_pos, 95) if len(glia_pos) > 0 else 1.0
    thresh_99 = np.percentile(glia_pos, 99) if len(glia_pos) > 0 else 1.0
    # Also try: any expression > 1 count (simple integer threshold)
    thresh_count2 = 1.5  # i.e. counts >= 2

    print(f"\n{name} / {gene}")
    print(f"  Glia  : n={glia_mask.sum()}, "
          f"% >0: {(glia_expr>0).mean()*100:.1f}%, "
          f"95th pct of pos glia: {thresh_95:.2f}, "
          f"99th pct: {thresh_99:.2f}")
    print(f"  Neuron: n={neuron_mask.sum()}, "
          f"% >0: {(neuron_expr>0).mean()*100:.1f}%")
    for t in [thresh_95, thresh_99, thresh_count2]:
        n_pass = (neuron_expr > t).sum()
        print(f"    Threshold >{t:.1f}: "
              f"{n_pass}/{neuron_mask.sum()} neurons pass "
              f"({n_pass/neuron_mask.sum()*100:.1f}%)")

    adata.obs[f'tracer_thresh95'] = expr > thresh_95
    adata.obs[f'tracer_thresh99'] = expr > thresh_99

# ── 4. Figure A: Expression distributions, neurons vs glia ───────────────────
fig, axes = plt.subplots(2, 3, figsize=(15, 9))

for row, (name, adata, expr, gene) in enumerate([
    ('NTS-R', retro_hq,  retro_expr,  'v36'),
    ('NTS-A', antero_hq, antero_expr, 'v147'),
]):
    glia_expr   = expr[adata.obs['cell_class'].values == 'Glia']
    neuron_expr = expr[adata.obs['cell_class'].values == 'Neuron']
    glia_pos    = glia_expr[glia_expr > 0]
    thresh_95   = np.percentile(glia_pos, 95) if len(glia_pos) > 0 else 1.0

    # Panel 1: full distribution (all counts > 0)
    ax = axes[row, 0]
    bins = np.arange(0.5, max(expr.max(), 5) + 1.5, 1)
    ax.hist(neuron_expr[neuron_expr > 0], bins=bins, alpha=0.7,
            color='steelblue', label=f'Neurons (n={neuron_expr.shape[0]})',
            density=True, edgecolor='white')
    ax.hist(glia_expr[glia_expr > 0], bins=bins, alpha=0.7,
            color='coral', label=f'Glia (n={glia_expr.shape[0]})',
            density=True, edgecolor='white')
    ax.axvline(thresh_95, color='black', ls='--', lw=1.5,
               label=f'Glia 95th pct = {thresh_95:.1f}')
    ax.set_xlabel(f'{gene} count (positive cells only)')
    ax.set_ylabel('Density')
    ax.set_title(f'{name}: expression distribution\n(positive cells only)',
                 fontweight='bold')
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    # Panel 2: ECDF — lets you read off any threshold visually
    ax = axes[row, 1]
    for e, label, color in [
        (neuron_expr, 'Neurons', 'steelblue'),
        (glia_expr,   'Glia',    'coral'),
    ]:
        e_pos = np.sort(e[e > 0])
        if len(e_pos) == 0:
            continue
        ecdf = np.arange(1, len(e_pos) + 1) / len(e_pos)
        ax.plot(e_pos, ecdf, color=color, lw=1.5, label=label)
    ax.axvline(thresh_95, color='black', ls='--', lw=1.5,
               label=f'Glia 95th pct = {thresh_95:.1f}')
    ax.set_xlabel(f'{gene} count (positive cells only)')
    ax.set_ylabel('Cumulative proportion')
    ax.set_title(f'{name}: ECDF of expression\n(threshold selection)',
                 fontweight='bold')
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    # Panel 3: % tracer+ neurons per cell type, above glia 95th pct threshold
    ax = axes[row, 2]
    neuron_adata = adata[adata.obs['cell_class'] == 'Neuron'].copy()
    ct_records = []
    for ct, grp in neuron_adata.obs.groupby(ct_col, observed=True):
        idx  = neuron_adata.obs_names.isin(grp.index)
        e    = expr[adata.obs['cell_class'].values == 'Neuron'][idx]
        n    = len(e)
        if n < 10:
            continue
        ct_records.append({
            'cell_type': ct,
            'n_cells':   n,
            'pct_above_bg': (e > thresh_95).mean() * 100,
            'pct_pos_any':  (e > 0).mean() * 100,
        })
    df_ct = (pd.DataFrame(ct_records)
               .set_index('cell_type')
               .sort_values('pct_above_bg', ascending=True))

    colors = ['#8B0000' if ct in lmx_types else 'steelblue'
              for ct in df_ct.index]
    y = range(len(df_ct))
    ax.barh(y, df_ct['pct_above_bg'], color=colors,
            edgecolor='white', linewidth=0.4, height=0.7, label='Above bg threshold')
    ax.barh(y, df_ct['pct_pos_any'] - df_ct['pct_above_bg'],
            left=df_ct['pct_above_bg'], color='lightsteelblue',
            edgecolor='white', linewidth=0.4, height=0.7, label='Count>0 but ≤ threshold')

    ax.set_yticks(y)
    ax.set_yticklabels(df_ct.index, fontsize=7)
    ax.set_xlabel('% cells with tracer expression')
    ax.set_title(f'{name}: % tracer+ neurons by cell type\n'
                 f'(dark = above glia 95th pct bg, red = Lmx1a types)',
                 fontweight='bold')
    ax.legend(fontsize=7, loc='lower right')
    ax.grid(axis='x', alpha=0.3)

    # Add n label
    for i, (ct, row_data) in enumerate(df_ct.iterrows()):
        ax.text(df_ct['pct_pos_any'].max() + 0.3, i,
                f'n={int(row_data["n_cells"])}',
                va='center', fontsize=6, color='gray')

plt.suptitle('Tracer expression stratification: neurons vs glia background\n'
             'Background threshold = 95th percentile of glial positive cells',
             fontweight='bold', fontsize=12)
plt.tight_layout()
plt.savefig(f'{OUTDIR}/Fig_tracer_expr_stratified.pdf', bbox_inches='tight', dpi=300)
plt.show()


# In[3]:


adata = sc.read_h5ad("HYP_integrated_2plex_scran.h5ad")
adata


# In[5]:


# ── Filter 2-plex to NTS retrograde animals only ───────────────────────────
nts_retro_ids = [168, 176, 149, 159, 165, 188]

adata.obs['animal_id'] = (adata.obs['sample_id']
                               .str.extract(r'(\d+)')
                               .astype(int).values)

adata = adata[adata.obs['animal_id'].isin(nts_retro_ids)].copy()
print(f"NTS retrograde cells: {adata.n_obs} "
      f"from animals {sorted(adata.obs['animal_id'].unique())}")


# In[6]:


# ── Create parquets at each score threshold ─────────────
import os
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

ct_col    = 'detailed_celltype_annotation'
score_col = 'score_detailed_celltype_annotation'
outdir    = "dvc_threshold_pickles"
os.makedirs(outdir, exist_ok=True)

# Rename v36 → v147 before extracting counts
if 'v36' in adata.var_names and 'v147' not in adata.var_names:
    adata.var_names = adata.var_names.str.replace('^v36$', 'v147', regex=True)
elif 'v36' in adata.var_names and 'v147' in adata.var_names:
    idx_v36  = list(adata.var_names).index('v36')
    idx_v147 = list(adata.var_names).index('v147')
    for layer_name in list(adata.layers.keys()):
        layer = adata.layers[layer_name].tolil() \
                if sparse.issparse(adata.layers[layer_name]) \
                else adata.layers[layer_name]
        layer[:, idx_v147] = layer[:, idx_v147] + layer[:, idx_v36]
        adata.layers[layer_name] = sparse.csr_matrix(layer) \
                                   if sparse.issparse(adata.layers[layer_name]) \
                                   else layer
    adata = adata[:, adata.var_names != 'v36'].copy()

print(f"v147 in var_names: {'v147' in adata.var_names}")
print(f"Total cells: {adata.n_obs}")
print(f"Score range: {adata.obs[score_col].min():.3f} – "
      f"{adata.obs[score_col].max():.3f}\n")

# Cell types to exclude
exclude_cts = {
    'Epithelial-like_cell', 'Pituitary_Pappa2', '35_Lhx9_Zic1',
    '01_Undefined', 'Pituitary_Crhr1', 'Undefined', 'Schwann_cell',
    'Pituitary_Gnrhr'
}

thresholds = np.arange(0.1, 1.0, 0.1)

print(f"{'Threshold':>10s}  {'N cells':>8s}  {'N cell types':>13s}  "
      f"{'Min CT size':>12s}  {'v147 pos%':>10s}")
print("-" * 60)

for thresh in thresholds:
    label = f"thresh_{thresh:.1f}".replace('.', 'p')

    # Filter by score and excluded cell types
    mask = (adata.obs[score_col] >= thresh) & \
           (~adata.obs[ct_col].isin(exclude_cts))
    sub  = adata[mask].copy()

    if sub.n_obs == 0:
        print(f"{thresh:>10.1f}  {'EMPTY':>8s}")
        continue

    # Drop cell types with < 10 cells after filtering
    ct_counts  = sub.obs[ct_col].value_counts()
    valid_cts  = ct_counts[ct_counts >= 10].index
    sub        = sub[sub.obs[ct_col].isin(valid_cts)].copy()

    # v147 positivity check
    v147_expr  = sub[:, 'v147'].layers['scran_normalization']
    v147_expr  = v147_expr.toarray().flatten() if sparse.issparse(v147_expr) \
                 else v147_expr.flatten()
    pct_pos    = (v147_expr > 0).mean() * 100

    print(f"{thresh:>10.1f}  {sub.n_obs:>8d}  {sub.obs[ct_col].nunique():>13d}  "
          f"{ct_counts[valid_cts].min():>12d}  {pct_pos:>9.1f}%")

    # Build and save parquets
    mat = sub.layers['scran_normalization']
    counts_df = pd.DataFrame(
        mat.toarray() if sparse.issparse(mat) else mat,
        index   = sub.obs_names,
        columns = sub.var_names
    )
    meta_df = sub.obs[[ct_col]].copy()

    counts_df.to_parquet(os.path.join(outdir, f"{label}_counts.parquet"))
    meta_df.to_parquet(  os.path.join(outdir, f"{label}_detailed.parquet"))

print(f"\nSaved parquets for {len(thresholds)} thresholds to '{outdir}/'")


# In[8]:


# ── Notebook 1: Load and compare v147 ESmu across thresholds ─────────────────

indir = "dvc_threshold_cellex"

records = []
for label in thresholds:
    label = round(label, 1)
    label_fix = f"{label}".replace('.', 'p')
    path = os.path.join(indir, f"thresh_{label_fix}_cellex.csv")
    if not os.path.exists(path):
        print(f"Missing: {label}")
        continue
    df = pd.read_csv(path, index_col=0)
    if 'v147' not in df.index:
        print(f"v147 not in {label}")
        continue
    thresh_val = label
    records.append(df.loc['v147'].rename(thresh_val))

df_thresh = pd.DataFrame(records).sort_index()

# Only show cell types with signal at any threshold
has_signal = (df_thresh > 0).any(axis=0)
df_thresh  = df_thresh.loc[:, has_signal]

print(f"Cell types with signal at any threshold: {has_signal.sum()}")
print(df_thresh.round(3).to_string())

# ── Figure: heatmap of ESmu vs threshold ─────────────────────────────────────
import matplotlib.pyplot as plt
import numpy as np

# Sort cell types by score at lowest threshold
ct_order = df_thresh.loc[0.1].sort_values(ascending=False).index

fig, axes = plt.subplots(1, 2, figsize=(16, 6))

# Panel A: heatmap
ax  = axes[0]
arr = df_thresh[ct_order].values.astype(float)
im  = ax.imshow(arr, aspect='auto', cmap='rocket_r',
                vmin=0, vmax=1, interpolation='nearest')
ax.set_yticks(range(len(df_thresh)))
ax.set_yticklabels([f"≥{t:.1f}  (n={int(adata[adata.obs[score_col]>=t].n_obs)})"
                    for t in df_thresh.index], fontsize=8)
ax.set_xticks(range(len(ct_order)))
ax.set_xticklabels(ct_order, rotation=90, ha='right', fontsize=7)
ax.set_title('v147 ESmu by score threshold\n(rows = threshold, cols = cell types)',
             fontweight='bold')
ax.set_ylabel('Score threshold')
plt.colorbar(im, ax=ax, fraction=0.02, label='ESmu')

# Panel B: line plot for top cell types
ax = axes[1]
top_cts = df_thresh.max(axis=0).nlargest(10).index
for ct in top_cts:
    ax.plot(df_thresh.index, df_thresh[ct], marker='o', ms=4,
            lw=1.5, label=ct)
ax.set_xlabel('Score threshold')
ax.set_ylabel('ESmu (v147)')
ax.set_title('Top cell types: ESmu stability\nacross score thresholds',
             fontweight='bold')
ax.legend(fontsize=7, bbox_to_anchor=(1.02, 1), loc='upper left')
ax.grid(alpha=0.3)
ax.set_xticks(df_thresh.index)

plt.suptitle('DVC CELLEX sensitivity to label transfer score threshold',
             fontweight='bold', fontsize=11)
plt.tight_layout()
plt.savefig("Fig_DVC_threshold_sensitivity.pdf", bbox_inches='tight', dpi=300)
plt.show()


# In[ ]:




