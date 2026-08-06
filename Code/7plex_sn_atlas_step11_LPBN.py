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
import adjustText
os.chdir("/maps/projects/perslab/people/qbx224/7plex_paper/2plex")
from pathlib import Path

outdir = Path("/maps/projects/perslab/people/qbx224/7plex_paper/2plex/")
sc.settings.figdir = outdir  # <- controls where scanpy writes


# In[3]:


adata = sc.read_h5ad("HYP_integrated_LPBN_labelled.h5ad")
adata


# In[4]:


adata.layers["RNA"] = adata.X
adata.layers


# In[5]:


print(adata.X[:5, :5].toarray())#Check if counts are the active layer


# In[6]:


adata.var = pd.DataFrame(index=adata.var_names)
adata


# In[7]:


# Remove specific columns from adata.var
columns_to_remove = ["nUMI_hto", "nCount_SCT", "nFeature_SCT", "SCT_snn_res.0.8", "SCT_snn_res.0.8", "integrated_snn_res.0.8"]  
adata.obs = adata.obs.drop(columns=columns_to_remove)
adata.obs


# In[8]:


sc.pl.embedding(adata,basis="X_umap",color="seurat_clusters", legend_loc="on data", legend_fontsize="xx-small")


# In[9]:


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


# In[10]:


adata.write("HYP_integrated_LPBN_anndata.h5ad")
adata


# In[11]:


adata.layers["counts"] = adata.X
adata


# In[12]:


os.environ['R_HOME'] = '/home/qbx224/.conda/envs/pca_env/lib/R'
import logging
import anndata2ri
import rpy2.rinterface_lib.callbacks as rcb
import rpy2.robjects as ro
rcb.logger.setLevel(logging.ERROR)
ro.pandas2ri.activate()
anndata2ri.activate()

get_ipython().run_line_magic('load_ext', 'rpy2.ipython')


# In[13]:


get_ipython().run_cell_magic('R', '', 'library(scran)\nlibrary(BiocParallel)\n')


# In[14]:


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


# In[15]:


get_ipython().run_cell_magic('R', '-o size_factors', '\nsize_factors = sizeFactors(\n    computeSumFactors(\n        SingleCellExperiment(\n            list(counts=data_mat)), \n            clusters = input_groups,\n            min.mean = 0.1,\n            BPPARAM = MulticoreParam()\n    )\n)\n')


# In[16]:


from scipy.sparse import csr_matrix
adata.obs["size_factors"] = size_factors
scran = adata.X / adata.obs["size_factors"].values[:, None]
scran = scran.tocsc()
adata.layers["scran_normalization"] = csr_matrix(sc.pp.log1p(scran))


# In[17]:


fig, axes = plt.subplots(1, 2, figsize=(10, 5))
p1 = sns.histplot(adata.obs["nCount_RNA"], bins=100, kde=False, ax=axes[0])
axes[0].set_title("Total counts")
p2 = sns.histplot(
    adata.layers["scran_normalization"].sum(1), bins=100, kde=False, ax=axes[1]
)
axes[1].set_title("log1p with Scran estimated size factors")
plt.show()


# In[18]:


adata.write("HYP_integrated_LPBN_scran.h5ad")
adata


# In[19]:


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
# %% Cell 2 — Load data
OUTDIR = "figures_LPBN"
os.makedirs(OUTDIR, exist_ok=True)

ct_col    = 'detailed_celltype_annotation'
score_col = 'score_detailed_celltype_annotation'

adata_2p  = sc.read("HYP_integrated_LPBN_scran.h5ad")
adata_ref = sc.read("/projects/perslab/people/qbx224/7plex_paper/atlas/allgenes_bipca_annotated.h5ad")

print(f"LPBN 2-plex: {adata_2p.n_obs} cells × {adata_2p.n_vars} genes")
print(f"Atlas:       {adata_ref.n_obs} cells × {adata_ref.n_vars} genes")


# In[20]:


# %% Cell 3 — Annotate animals from surgery sheet
# All animals had LPBN retrograde with v147r
# Treatment groups from sort sheet
tx_833    = [52, 55, 63, 67, 50, 59, 62]
tx_saline = [54, 71]
tx_none   = [56, 68, 73, 58, 61]

# Animals re-sequenced on lane E7 (appear in both E3 and E7)
reseq_ids = [55, 56, 63, 67, 68, 73]

adata_2p.obs['animal_id'] = (adata_2p.obs['sample_id']
                               .str.extract(r'(\d+)')
                               .astype(int).values)

adata_2p.obs['treatment'] = 'unknown'
adata_2p.obs.loc[adata_2p.obs['animal_id'].isin(tx_833),    'treatment'] = '833'
adata_2p.obs.loc[adata_2p.obs['animal_id'].isin(tx_saline), 'treatment'] = 'saline'
adata_2p.obs.loc[adata_2p.obs['animal_id'].isin(tx_none),   'treatment'] = 'none'

adata_2p.obs['resequenced'] = adata_2p.obs['animal_id'].isin(reseq_ids)

print("\nAnimal × lane breakdown:")
print(pd.crosstab(adata_2p.obs['animal_id'], adata_2p.obs['old.ident']))
print("\nTreatment breakdown:")
print(adata_2p.obs['treatment'].value_counts())
print("\nRe-sequenced animals (2 lanes):")
print(adata_2p.obs[adata_2p.obs['resequenced']]
      .groupby('animal_id')['old.ident'].unique())


# In[21]:


# %% Cell 4 — Label transfer QC and filter
fig, ax = plt.subplots(figsize=(6, 4))
scores = adata_2p.obs[score_col]
ax.hist(scores, bins=50, color='steelblue', edgecolor='white', linewidth=0.3)
ax.axvline(0.1, color='firebrick', ls='--', lw=1.5, label='threshold 0.1')
above = (scores > 0.1).mean() * 100
ax.set_title(f'LPBN label transfer scores\nmedian={scores.median():.3f}, >0.1: {above:.0f}%')
ax.set_xlabel('Prediction score')
ax.legend()
plt.tight_layout()
plt.savefig(f'{OUTDIR}/QC_label_transfer_LPBN.pdf')
plt.show()

lpbn_hq = adata_2p[adata_2p.obs[score_col] > 0.1].copy()
print(f"\nAfter QC filter: {lpbn_hq.n_obs}/{adata_2p.n_obs} "
      f"({lpbn_hq.n_obs/adata_2p.n_obs*100:.1f}%)")


# In[22]:


# %% Cell 5 — Extract atlas PBN sort=1.0 cells for comparison
# Mirror the animal×lane approach established for the HYP dataset
ref_sort = adata_ref.obs['sort'].values.astype(float)

ref_records = []
for (animal, lane), grp in adata_ref.obs.groupby(
        ['hash.mcl.ID', 'orig.ident'], observed=True):
    mask = adata_ref.obs_names.isin(grp.index)
    s    = ref_sort[mask].mean()
    ref_records.append({'animal': animal, 'lane': lane,
                        'sort_prob': s, 'n_cells': mask.sum()})

df_ref_lanes = pd.DataFrame(ref_records)
sort100_combos = df_ref_lanes[df_ref_lanes['sort_prob'] == 1.0]

print(f"\nAtlas sort=1.0 animal×lane combos: {len(sort100_combos)}")
print(f"Unique animals: {sort100_combos['animal'].nunique()}")

# Subset atlas to sort=1.0 cells
atlas_100_idx = []
for _, row in sort100_combos.iterrows():
    idx = adata_ref.obs[
        (adata_ref.obs['hash.mcl.ID'] == row['animal']) &
        (adata_ref.obs['orig.ident']  == row['lane'])
    ].index
    atlas_100_idx.extend(idx.tolist())

adata_atlas_100 = adata_ref[adata_ref.obs_names.isin(atlas_100_idx)].copy()
print(f"Atlas sort=1.0 cells: {adata_atlas_100.n_obs}")


# In[23]:


# %% Cell 6 — v147 expression by cell type: LPBN 2-plex vs atlas sort=1.0
def get_expr_scran(adata, gene):
    if gene not in adata.var_names:
        print(f"  WARNING: {gene} not found"); return None
    if 'scran_normalization' in adata.layers:
        layer = adata[:, gene].layers['scran_normalization']
    else:
        print(f"  NOTE: scran layer absent, using .X")
        layer = adata[:, gene].X
    return layer.toarray().flatten() if hasattr(layer, 'toarray') else layer.flatten()

def celltype_v147_summary(adata, ct_col='detailed_celltype_annotation',
                           min_cells=10):
    expr = get_expr_scran(adata, 'v147')
    if expr is None:
        return None
    records = []
    for ct, grp in adata.obs.groupby(ct_col, observed=True):
        e = expr[adata.obs_names.isin(grp.index)]
        if len(e) < min_cells:
            continue
        records.append({
            'cell_type':   ct,
            'n_cells':     len(e),
            'pct_pos':     (e > 0).mean() * 100,
            'mean_expr':   e.mean(),
            'mean_pos':    e[e > 0].mean() if (e > 0).any() else 0,
        })
    return (pd.DataFrame(records)
              .set_index('cell_type')
              .sort_values('pct_pos', ascending=False))

df_lpbn  = celltype_v147_summary(lpbn_hq)
df_atlas = celltype_v147_summary(adata_atlas_100)

print("\nLPBN 2-plex — v147 by cell type:")
print(df_lpbn.round(3).to_string())
print("\nAtlas sort=1.0 — v147 by cell type:")
print(df_atlas.round(3).to_string())

# %% Cell 7 — Figure: side-by-side heatmap
all_cts = sorted(set(df_lpbn.index) | set(df_atlas.index))

df_compare = pd.DataFrame({
    'LPBN 2-plex (v147-retro, sort=100%)': 
        df_lpbn['pct_pos'].reindex(all_cts),
    'Atlas sort=1.0 (v147-retro)':         
        df_atlas['pct_pos'].reindex(all_cts),
})

# Sort by LPBN % positive
df_compare = df_compare.sort_values(
    'LPBN 2-plex (v147-retro, sort=100%)', ascending=False)

fig, axes = plt.subplots(1, 3, figsize=(16, max(8, len(all_cts) * 0.28)),
                          gridspec_kw={'width_ratios': [3, 1, 1]})

cmap = plt.get_cmap('rocket_r')
cmap.set_bad('lightgray')

arr = df_compare.values.astype(float)
vmax = np.nanmax(arr)
im = axes[0].imshow(np.ma.masked_invalid(arr), aspect='auto',
                     cmap=cmap, vmin=0, vmax=vmax, interpolation='nearest')

axes[0].set_xticks(range(2))
axes[0].set_xticklabels(df_compare.columns, rotation=20, ha='right', fontsize=9)
axes[0].set_yticks(range(len(df_compare)))
axes[0].set_yticklabels(df_compare.index, fontsize=7)
axes[0].set_title('% v147-positive nuclei\nby cell type', fontweight='bold')

for ri in range(arr.shape[0]):
    for ci in range(arr.shape[1]):
        v = arr[ri, ci]
        if not np.isnan(v) and v > 0:
            tc = 'white' if v > vmax * 0.6 else 'black'
            axes[0].text(ci, ri, f'{v:.1f}', ha='center', va='center',
                         fontsize=5.5, color=tc)

plt.colorbar(im, ax=axes[0], fraction=0.02, pad=0.01, label='% positive')

# N cells bars
for ax_i, (col, color, key) in enumerate([
    (axes[1], 'steelblue', df_atlas),
    (axes[2], 'firebrick', df_lpbn),
], start=1):
    n = key['n_cells'].reindex(df_compare.index, fill_value=0)
    axes[ax_i].barh(range(len(df_compare)), n,
                    color=color, alpha=0.7, height=0.6)
    axes[ax_i].set_yticks([])
    axes[ax_i].set_xlabel('N cells')
    axes[ax_i].set_title(
        'N cells\n(atlas)' if ax_i == 1 else 'N cells\n(LPBN)',
        fontsize=9, fontweight='bold')
    axes[ax_i].grid(axis='x', alpha=0.3)

plt.suptitle('v147 retrograde tracer detection by cell type\n'
             'LPBN 2-plex (sort=100%) vs atlas sort=1.0 lanes',
             fontweight='bold', fontsize=11)
plt.tight_layout()
plt.savefig(f'{OUTDIR}/Fig_v147_by_celltype_LPBN_vs_atlas.pdf',
            bbox_inches='tight', dpi=300)
plt.show()


# In[24]:


# %% Cell 8 — Spearman correlation between datasets
shared = df_compare.dropna()
if len(shared) >= 5:
    from scipy.stats import spearmanr
    r, p = spearmanr(shared.iloc[:, 0], shared.iloc[:, 1])
    print(f"\nSpearman ρ (LPBN vs atlas % positive per cell type): "
          f"ρ={r:.3f}, p={p:.3e}, n={len(shared)} cell types")


# In[25]:


# %% Cell — Rank-shift test: LPBN 2-plex vs atlas PBN animals
# Unit of analysis = per-animal cell type proportion
# Test = Mann-Whitney U on animal-level proportions (Wilcoxon rank-sum)

# ── 1. Identify PBN animals in the reference atlas ───────────────────────────
# The atlas 'area' column encodes the injection target
print("Atlas area values:")
print(adata_ref.obs['area'].value_counts())

# Subset to PBN animals regardless of sort probability
pbn_mask = adata_ref.obs['area'].str.contains('PBN', case=False, na=False)
adata_pbn = adata_ref[pbn_mask].copy()
print(f"\nAtlas PBN cells: {adata_pbn.n_obs} from "
      f"{adata_pbn.obs['hash.mcl.ID'].nunique()} animals")
print(adata_pbn.obs['hash.mcl.ID'].value_counts())

# ── 2. Compute per-animal cell type proportions ───────────────────────────────
def per_animal_proportions(adata, animal_col, ct_col='detailed_celltype_annotation'):
    """Returns DataFrame: rows=animals, cols=cell types, values=proportions."""
    records = {}
    for animal, grp in adata.obs.groupby(animal_col, observed=True):
        ct_counts = grp[ct_col].value_counts()
        records[animal] = ct_counts / ct_counts.sum()
    df = pd.DataFrame(records).T.fillna(0)
    return df

# LPBN 2-plex: animal_id already computed
props_lpbn = per_animal_proportions(lpbn_hq, 'animal_id')
props_pbn  = per_animal_proportions(adata_pbn, 'hash.mcl.ID')

print(f"\nLPBN 2-plex: {len(props_lpbn)} animals × {props_lpbn.shape[1]} cell types")
print(f"Atlas PBN:   {len(props_pbn)} animals × {props_pbn.shape[1]} cell types")

# ── 3. Mann-Whitney U on per-animal proportions ───────────────────────────────
from scipy.stats import mannwhitneyu

all_cts = sorted(set(props_lpbn.columns) | set(props_pbn.columns))

records = []
for ct in all_cts:
    lpbn_vals = props_lpbn[ct].values * 100 if ct in props_lpbn.columns \
                else np.zeros(len(props_lpbn))
    pbn_vals  = props_pbn[ct].values  * 100 if ct in props_pbn.columns  \
                else np.zeros(len(props_pbn))

    # Skip cell types absent in both
    if lpbn_vals.sum() == 0 and pbn_vals.sum() == 0:
        continue

    stat, pval = mannwhitneyu(lpbn_vals, pbn_vals, alternative='two-sided')

    records.append({
        'cell_type':       ct,
        'median_lpbn':     np.median(lpbn_vals),
        'median_pbn':      np.median(pbn_vals),
        'mean_lpbn':       np.mean(lpbn_vals),
        'mean_pbn':        np.mean(pbn_vals),
        'n_animals_lpbn':  len(lpbn_vals),
        'n_animals_pbn':   len(pbn_vals),
        'fold_change':     (np.mean(lpbn_vals) / np.mean(pbn_vals))
                           if np.mean(pbn_vals) > 0 else np.inf,
        'log2fc':          np.log2((np.mean(lpbn_vals) + 0.01) /
                                   (np.mean(pbn_vals)  + 0.01)),
        'mwu_stat':        stat,
        'pvalue':          pval,
    })

df_rank = pd.DataFrame(records)
_, fdr, _, _ = multipletests(df_rank['pvalue'], method='fdr_bh')
df_rank['fdr']         = fdr
df_rank['significant'] = df_rank['fdr'] < 0.05
df_rank['direction']   = np.where(df_rank['log2fc'] > 0,
                                   'Enriched in LPBN', 'Depleted in LPBN')

df_rank = df_rank.sort_values('log2fc', ascending=False)
df_rank.to_csv(f'{OUTDIR}/LPBN_vs_PBN_rankshift.csv', index=False)

# ── 4. Print results ──────────────────────────────────────────────────────────
print("\n" + "="*75)
print("ENRICHED in LPBN vs atlas PBN (FDR<0.05, animal-level MWU):")
print(f"{'Cell type':30s} {'med_LPBN':>9s} {'med_PBN':>9s} "
      f"{'log2FC':>7s} {'FDR':>10s}")
print("-"*75)
sig_up = df_rank[df_rank['significant'] & (df_rank['log2fc'] > 0)]
for _, r in sig_up.iterrows():
    print(f"{r['cell_type']:30s} {r['median_lpbn']:>8.2f}% {r['median_pbn']:>8.2f}% "
          f"{r['log2fc']:>7.2f} {r['fdr']:>10.2e}")

print("\n" + "="*75)
print("DEPLETED in LPBN vs atlas PBN (FDR<0.05):")
print(f"{'Cell type':30s} {'med_LPBN':>9s} {'med_PBN':>9s} "
      f"{'log2FC':>7s} {'FDR':>10s}")
print("-"*75)
sig_dn = df_rank[df_rank['significant'] & (df_rank['log2fc'] < 0)]
for _, r in sig_dn.iterrows():
    print(f"{r['cell_type']:30s} {r['median_lpbn']:>8.2f}% {r['median_pbn']:>8.2f}% "
          f"{r['log2fc']:>7.2f} {r['fdr']:>10.2e}")

# ── 5. Volcano plot ───────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(9, 7))

# Define glia for colouring
glial_types_hyp = {
    'Astrocyte', 'Oligodendrocyte', 'OPC', 'Microglia', 'VLMC',
    'Ependymal_cell', 'Tanycyte', 'Pericyte', 'Schwann_cell',
    'Pituicyte', 'Epithelial-like_cell', 'Endothelial_cell',
}

for _, r in df_rank.iterrows():
    is_glia = r['cell_type'] in glial_types_hyp
    if r['significant'] and not is_glia:
        color, zorder, size = ('#C0392B' if r['log2fc'] > 0 else '#2980B9'), 4, 60
    elif r['significant'] and is_glia:
        color, zorder, size = ('#E67E22' if r['log2fc'] > 0 else '#8E44AD'), 4, 60
    else:
        color, zorder, size = '#BBBBBB', 2, 25

    ax.scatter(r['log2fc'], -np.log10(r['fdr'] + 1e-10),
               color=color, s=size, alpha=0.8, zorder=zorder,
               edgecolors='white', linewidths=0.3)

# Label significant neurons
for _, r in df_rank[df_rank['significant'] &
                    ~df_rank['cell_type'].isin(glial_types_hyp)].iterrows():
    ax.annotate(r['cell_type'],
                xy=(r['log2fc'], -np.log10(r['fdr'] + 1e-10)),
                xytext=(4, 0), textcoords='offset points',
                fontsize=6.5, va='center',
                color='#C0392B' if r['log2fc'] > 0 else '#2980B9')

ax.axvline(0, color='k', lw=0.8, ls='--', alpha=0.5)
ax.axhline(-np.log10(0.05), color='grey', lw=0.8, ls=':', alpha=0.5,
           label='FDR=0.05')

ax.set_xlabel('log2 fold change (LPBN / atlas PBN)', fontsize=10)
ax.set_ylabel('-log10(FDR)', fontsize=10)
ax.set_title('Cell type rank-shift: LPBN 2-plex vs atlas PBN animals\n'
             '(Mann-Whitney U on per-animal proportions)',
             fontweight='bold')

legend_el = [
    Patch(fc='#C0392B', label='Neuron enriched in LPBN (FDR<0.05)'),
    Patch(fc='#2980B9', label='Neuron depleted in LPBN (FDR<0.05)'),
    Patch(fc='#E67E22', label='Glia enriched (FDR<0.05)'),
    Patch(fc='#8E44AD', label='Glia depleted (FDR<0.05)'),
    Patch(fc='#BBBBBB', label='Not significant'),
    plt.Line2D([0],[0], color='grey', ls=':', lw=1, label='FDR=0.05'),
]
ax.legend(handles=legend_el, fontsize=7, loc='upper left')
ax.grid(alpha=0.2)

plt.tight_layout()
plt.savefig(f'{OUTDIR}/Fig_LPBN_vs_PBN_volcano.pdf', bbox_inches='tight', dpi=300)
plt.show()

# ── 6. Stripchart for top hits: show per-animal distributions ─────────────────
top_hits = (df_rank[df_rank['significant']]
            .reindex(df_rank[df_rank['significant']]['log2fc']
                     .abs().sort_values(ascending=False).index)
            .head(16)['cell_type'].tolist())

n_cols = 4
n_rows = int(np.ceil(len(top_hits) / n_cols))
fig, axes = plt.subplots(n_rows, n_cols,
                          figsize=(n_cols * 3, n_rows * 3))
axes = axes.flatten()

for i, ct in enumerate(top_hits):
    ax = axes[i]
    lpbn_vals = props_lpbn[ct].values * 100 if ct in props_lpbn.columns \
                else np.zeros(len(props_lpbn))
    pbn_vals  = props_pbn[ct].values  * 100 if ct in props_pbn.columns  \
                else np.zeros(len(props_pbn))

    for xi, (vals, color, label) in enumerate([
        (pbn_vals,  '#2980B9', 'Atlas PBN'),
        (lpbn_vals, '#C0392B', 'LPBN 2-plex'),
    ]):
        jitter = np.random.uniform(-0.12, 0.12, len(vals))
        ax.scatter(np.full(len(vals), xi) + jitter, vals,
                   color=color, s=30, alpha=0.8, zorder=3,
                   edgecolors='k', linewidths=0.3)
        ax.plot([xi - 0.25, xi + 0.25],
                [np.median(vals), np.median(vals)],
                color='k', lw=2, zorder=4)

    r = df_rank[df_rank['cell_type'] == ct].iloc[0]
    ax.set_title(f"{ct}\nlog2FC={r['log2fc']:.2f}, FDR={r['fdr']:.2e}",
                 fontsize=7, fontweight='bold')
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['Atlas\nPBN', 'LPBN\n2-plex'], fontsize=7)
    ax.set_ylabel('% of nuclei', fontsize=7)
    ax.grid(axis='y', alpha=0.3)

# Hide unused panels
for j in range(len(top_hits), len(axes)):
    axes[j].set_visible(False)

plt.suptitle('Per-animal proportions for top rank-shifted cell types',
             fontweight='bold', fontsize=11)
plt.tight_layout()
plt.savefig(f'{OUTDIR}/Fig_LPBN_vs_PBN_stripchart.pdf', bbox_inches='tight', dpi=300)
plt.show()


# In[26]:


# %% Cell — v147 neuron vs glia expression in LPBN dataset

glial_types_hyp = {
    'Astrocyte', 'Oligodendrocyte', 'OPC', 'Microglia', 'VLMC',
    'Ependymal_cell', 'Tanycyte', 'Pericyte', 'Schwann_cell',
    'Pituicyte', 'Epithelial-like_cell', 'Endothelial_cell',
}
pituitary_types = {c for c in lpbn_hq.obs[ct_col].unique() if 'Pituitary' in c}
non_neuronal    = glial_types_hyp | pituitary_types

lpbn_hq.obs['cell_class'] = 'Neuron'
lpbn_hq.obs.loc[lpbn_hq.obs[ct_col].isin(glial_types_hyp),  'cell_class'] = 'Glia'
lpbn_hq.obs.loc[lpbn_hq.obs[ct_col].isin(pituitary_types),  'cell_class'] = 'Pituitary'

print("LPBN cell class breakdown:")
print(lpbn_hq.obs['cell_class'].value_counts())

# ── Extract v147 expression ───────────────────────────────────────────────────
v147_expr = get_expr_scran(lpbn_hq, 'v147')
if v147_expr is None:
    raise ValueError("v147 not found in lpbn_hq")

neuron_mask = lpbn_hq.obs['cell_class'].values == 'Neuron'
glia_mask   = lpbn_hq.obs['cell_class'].values == 'Glia'

neuron_expr = v147_expr[neuron_mask]
glia_expr   = v147_expr[glia_mask]

# Background threshold from glia
glia_pos  = glia_expr[glia_expr > 0]
thresh_95 = np.percentile(glia_pos, 95) if len(glia_pos) > 0 else 1.0
thresh_99 = np.percentile(glia_pos, 99) if len(glia_pos) > 0 else 1.0

print(f"\nGlia  : n={glia_mask.sum()}, "
      f"% >0 = {(glia_expr>0).mean()*100:.1f}%, "
      f"95th pct of pos = {thresh_95:.3f}, "
      f"99th pct = {thresh_99:.3f}")
print(f"Neuron: n={neuron_mask.sum()}, "
      f"% >0 = {(neuron_expr>0).mean()*100:.1f}%")
for t in [thresh_95, thresh_99]:
    n_pass = (neuron_expr > t).sum()
    print(f"  Threshold >{t:.3f}: {n_pass}/{neuron_mask.sum()} neurons "
          f"({n_pass/neuron_mask.sum()*100:.1f}%)")

# ── Figure: 3 panels ─────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# Panel 1: expression distribution (positive cells only)
ax = axes[0]
bins = np.linspace(0, max(v147_expr.max(), 1), 60)
ax.hist(neuron_expr[neuron_expr > 0], bins=bins, alpha=0.7,
        color='steelblue', label=f'Neurons (n={neuron_mask.sum()})',
        density=True, edgecolor='white')
ax.hist(glia_expr[glia_expr > 0], bins=bins, alpha=0.7,
        color='coral', label=f'Glia (n={glia_mask.sum()})',
        density=True, edgecolor='white')
ax.axvline(thresh_95, color='black', ls='--', lw=1.5,
           label=f'Glia 95th pct = {thresh_95:.3f}')
ax.set_xlabel('v147 scran-normalized expression\n(positive cells only)')
ax.set_ylabel('Density')
ax.set_title('LPBN: v147 expression distribution\n(positive cells only)',
             fontweight='bold')
ax.legend(fontsize=8)
ax.grid(alpha=0.3)

# Panel 2: ECDF
ax = axes[1]
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
           label=f'Glia 95th pct = {thresh_95:.3f}')
ax.set_xlabel('v147 scran-normalized expression\n(positive cells only)')
ax.set_ylabel('Cumulative proportion')
ax.set_title('LPBN: ECDF of v147 expression\n(threshold selection)',
             fontweight='bold')
ax.legend(fontsize=8)
ax.grid(alpha=0.3)

# Panel 3: % above background per cell type (neurons only)
ax = axes[2]
neuron_adata = lpbn_hq[lpbn_hq.obs['cell_class'] == 'Neuron'].copy()
ct_records = []
for ct, grp in neuron_adata.obs.groupby(ct_col, observed=True):
    idx = neuron_adata.obs_names.isin(grp.index)
    e   = v147_expr[neuron_mask][idx]
    if len(e) < 10:
        continue
    ct_records.append({
        'cell_type':    ct,
        'n_cells':      len(e),
        'pct_above_bg': (e > thresh_95).mean() * 100,
        'pct_pos_any':  (e > 0).mean() * 100,
    })

df_ct = (pd.DataFrame(ct_records)
           .set_index('cell_type')
           .sort_values('pct_above_bg', ascending=True))

y      = range(len(df_ct))
colors = ['#C0392B' if ct in sig_up['cell_type'].values else 'steelblue'
          for ct in df_ct.index]

ax.barh(y, df_ct['pct_above_bg'], color=colors,
        edgecolor='white', linewidth=0.4, height=0.7,
        label='Above bg (>glia 95th pct)')
ax.barh(y, df_ct['pct_pos_any'] - df_ct['pct_above_bg'],
        left=df_ct['pct_above_bg'],
        color='lightsteelblue', edgecolor='white', linewidth=0.4, height=0.7,
        label='Count>0 but ≤ threshold')

ax.set_yticks(y)
ax.set_yticklabels(df_ct.index, fontsize=7)
ax.set_xlabel('% cells with v147 expression')
ax.set_title('LPBN neurons: % v147+ by cell type\n'
             '(red = enriched vs atlas PBN)',
             fontweight='bold')
ax.legend(fontsize=7, loc='lower right')
ax.grid(axis='x', alpha=0.3)

# n label
for i, (ct, row) in enumerate(df_ct.iterrows()):
    ax.text(df_ct['pct_pos_any'].max() + 0.3, i,
            f'n={int(row["n_cells"])}',
            va='center', fontsize=6, color='gray')

plt.suptitle('v147 expression stratification: LPBN neurons vs glia background\n'
             'Background threshold = 95th percentile of glial positive cells',
             fontweight='bold', fontsize=11)
plt.tight_layout()
plt.savefig(f'{OUTDIR}/Fig_v147_neuron_glia_LPBN.pdf', bbox_inches='tight', dpi=300)
plt.show()

# ── Mann-Whitney summary ──────────────────────────────────────────────────────
from scipy.stats import mannwhitneyu
stat, pval = mannwhitneyu(neuron_expr, glia_expr, alternative='two-sided')
print(f"\nMann-Whitney U (all expression, neurons vs glia): "
      f"U={stat:.0f}, p={pval:.3e}")
stat2, pval2 = mannwhitneyu(neuron_expr[neuron_expr > 0],
                              glia_expr[glia_expr > 0],
                              alternative='two-sided')
print(f"Mann-Whitney U (positive cells only):             "
      f"U={stat2:.0f}, p={pval2:.3e}")


# In[2]:


adata = sc.read_h5ad("HYP_integrated_LPBN_scran.h5ad")
adata


# In[4]:


# ── Create parquets at each score threshold ─────────────
import os
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

ct_col    = 'detailed_celltype_annotation'
score_col = 'score_detailed_celltype_annotation'
outdir    = "lpbn_threshold_pickles"
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


# In[5]:


# ── Notebook 1: Load and compare v147 ESmu across thresholds ─────────────────

indir = "lpbn_threshold_cellex"

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

plt.suptitle('LPBN CELLEX sensitivity to label transfer score threshold',
             fontweight='bold', fontsize=11)
plt.tight_layout()
plt.savefig("Fig_LPBN_threshold_sensitivity.pdf", bbox_inches='tight', dpi=300)
plt.show()


# In[3]:


# Annotate animals from surgery sheet
# All animals had LPBN retrograde with v147r
# Treatment groups from sort sheet
tx_833    = [52, 55, 63, 67, 50, 59, 62]
tx_saline = [54, 71]
tx_none   = [56, 68, 73, 58, 61]

# Animals re-sequenced on lane E7 (appear in both E3 and E7)
reseq_ids = [55, 56, 63, 67, 68, 73]

adata.obs['animal_id'] = (adata.obs['sample_id']
                               .str.extract(r'(\d+)')
                               .astype(int).values)

adata.obs['treatment'] = 'unknown'
adata.obs.loc[adata.obs['animal_id'].isin(tx_833),    'treatment'] = '833'
adata.obs.loc[adata.obs['animal_id'].isin(tx_saline), 'treatment'] = 'saline'
adata.obs.loc[adata.obs['animal_id'].isin(tx_none),   'treatment'] = 'none'

adata.obs['resequenced'] = adata.obs['animal_id'].isin(reseq_ids)

print("\nAnimal × lane breakdown:")
print(pd.crosstab(adata.obs['animal_id'], adata.obs['old.ident']))
print("\nTreatment breakdown:")
print(adata.obs['treatment'].value_counts())
print("\nRe-sequenced animals (2 lanes):")
print(adata.obs[adata.obs['resequenced']]
      .groupby('animal_id')['old.ident'].unique())


# In[ ]:


import scanpy as sc
import pandas as pd
import numpy as np
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import warnings
warnings.filterwarnings('ignore')

# ── 1. Confirm we're using raw counts ────────────────────────────────────────
# Use the 'counts' layer (integers). Double-check it's actually integer counts.
assert np.allclose(adata.layers['counts'].data % 1, 0), \
    "counts layer doesn't look like integers — check your layer!"

# ── 2. Pseudobulk: aggregate per (animal_id × cell_type) ────────────────────
# Re-sequenced animals appear in two lanes but same animal_id → 
# summing across all cells from that animal is correct behaviour.

def make_pseudobulk(adata, groupby_cols, layer='counts'):
    """Sum raw counts over groupby_cols, return (count_df, meta_df)."""
    import scipy.sparse as sp

    groups = adata.obs[groupby_cols].astype(str).agg('__'.join, axis=1)
    unique_groups = groups.unique()

    rows = []
    meta_rows = []
    for grp in unique_groups:
        mask = groups == grp
        mat = adata.layers[layer][mask.values]  # .values converts Series → numpy bool array
        if sp.issparse(mat):
            mat = mat.toarray()
        rows.append(mat.sum(axis=0))
        # grab metadata from first cell (all same within group)
        meta_rows.append(adata.obs.loc[mask].iloc[0][groupby_cols].to_dict())

    counts_df = pd.DataFrame(
        np.vstack(rows),
        index=unique_groups,
        columns=adata.var_names
    )
    meta_df = pd.DataFrame(meta_rows, index=unique_groups)
    return counts_df, meta_df


pb_counts, pb_meta = make_pseudobulk(
    adata,
    groupby_cols=['animal_id', 'detailed_celltype_annotation'],
    layer='counts'
)

# Also pull n_cells per pseudobulk sample (useful for QC)
pb_ncells = (
    adata.obs
    .assign(_grp=lambda d: d['animal_id'].astype(str) + '__' + d['detailed_celltype_annotation'].astype(str))
    .groupby('_grp').size()
    .rename('n_cells')
)
pb_meta['n_cells'] = pb_ncells

print(f"Pseudobulk matrix: {pb_counts.shape[0]} samples × {pb_counts.shape[1]} genes")
print(pb_meta['detailed_celltype_annotation'].value_counts().head(10))


# In[7]:


# ── 3. DEG per cell type — 833 vs (saline + none) ───────────────────────────

MIN_CELLS_PER_SAMPLE  = 10
MIN_ANIMALS_PER_GROUP = 3

# Build treatment lookup once
animal_tx = (adata.obs[['animal_id','treatment']]
             .drop_duplicates('animal_id')
             .set_index('animal_id')['treatment'])

all_results = {}

for ct in pb_meta['detailed_celltype_annotation'].unique():
    mask     = pb_meta['detailed_celltype_annotation'] == ct
    ct_counts = pb_counts[mask].copy()
    ct_meta   = pb_meta[mask].copy()

    # ── QC: drop low-cell pseudobulk samples ─────────────────────────────
    keep      = ct_meta['n_cells'] >= MIN_CELLS_PER_SAMPLE
    ct_counts = ct_counts[keep]
    ct_meta   = ct_meta[keep].copy()

    # ── Map treatment, collapse saline+none → 'control' ──────────────────
    ct_meta['treatment'] = ct_meta['animal_id'].map(animal_tx)
    ct_meta['group'] = ct_meta['treatment'].map(
        {'833': '833', 'saline': 'control', 'none': 'control'}
    )

    # Drop any unknowns
    ct_meta   = ct_meta[ct_meta['group'].isin(['833', 'control'])]
    ct_counts = ct_counts.loc[ct_meta.index]

    # ── Require minimum animals per group ─────────────────────────────────
    group_counts = ct_meta['group'].value_counts()
    if not all(group_counts.get(g, 0) >= MIN_ANIMALS_PER_GROUP
               for g in ['833', 'control']):
        print(f"  [SKIP] {ct}: n={group_counts.to_dict()}")
        continue

    # ── Drop all-zero genes ───────────────────────────────────────────────
    expressed = ct_counts.sum(axis=0) > 0
    ct_counts = ct_counts.loc[:, expressed]

    print(f"\n[{ct}]  n_833={group_counts['833']}  "
          f"n_control={group_counts['control']}  genes={expressed.sum()}")

    try:
        dds = DeseqDataSet(
            counts         = ct_counts.astype(int),
            metadata       = ct_meta[['group']],
            design_factors = 'group',
            ref_level      = ['group', 'control'],
            refit_cooks    = True,
            quiet          = True
        )
        dds.deseq2()

        stat_res = DeseqStats(
            dds,
            contrast = ['group', '833', 'control'],
            quiet    = True
        )
        stat_res.summary()

        res = stat_res.results_df.copy()
        res['cell_type']  = ct
        res['n_833']      = group_counts['833']
        res['n_control']  = group_counts['control']
        all_results[ct]   = res

    except Exception as e:
        print(f"  [ERROR] {ct}: {e}")
        continue

# ── 4. Collate ───────────────────────────────────────────────────────────────
if not all_results:
    print("No cell types passed filters — check MIN_CELLS_PER_SAMPLE or group sizes")
else:
    results_df = pd.concat(all_results.values())
    sig = results_df[
        (results_df['padj'] < 0.05) &
        (results_df['log2FoldChange'].abs() > 0.5)
    ].sort_values(['cell_type', 'padj'])

    print("\n══════════ SIGNIFICANT DEGs (833 vs control) ══════════")
    print(sig.groupby('cell_type').size().rename('n_sig').sort_values(ascending=False))

    results_df.to_csv("pseudobulk_833_vs_control_all.csv")
    sig.to_csv("pseudobulk_833_vs_control_significant.csv")


# In[12]:


# 1. Pull the actual DEGs for your top hit to understand the biology
npsr1_degs = results_df[
    (results_df['cell_type'] == '26_Dach2_Npsr1') &
    (results_df['padj'] < 0.05) &
    (results_df['log2FoldChange'].abs() > 0.5)
].sort_values('log2FoldChange', ascending=False)

print(npsr1_degs[['log2FoldChange','padj']].to_string())

# 2. Run GO/pathway enrichment on those 29 genes
# (use gseapy or submit gene list to Enrichr)
import gseapy as gp
enr = gp.enrichr(
    gene_list   = npsr1_degs.index.tolist(),
    gene_sets   = ['GO_Biological_Process_2023', 'KEGG_2021_Human'],
    organism    = 'mouse',
    outdir      = 'enrichr_NPSR1_833'
)

# 3. Characterize the Undefined population
undef_degs = results_df[
    (results_df['cell_type'] == '01_Undefined') &
    (results_df['padj'] < 0.05)
].index.tolist()
print("Undefined marker genes:", undef_degs)


# In[15]:


# ── 0. Cell type order — from full-atlas cellex columns ──────────────────────
cellex_full    = pd.read_csv("all_cellex_results/adata_cellex.csv", index_col=0)
cell_type_order = list(cellex_full.columns)
celltypes       = cell_type_order  # same list, used interchangeably below

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.gridspec import GridSpec


# ── 1. v147 snRNA projection scores (one row per area from per-area csvs) ────
indir   = "/maps/projects/perslab/people/qbx224/7plex_paper/atlas/area_cellex_results"
records = []
for fn in os.listdir(indir):
    if not fn.endswith("_detailed_cellex.csv"):
        continue
    area = fn.replace("_detailed_cellex.csv", "")
    df   = pd.read_csv(os.path.join(indir, fn), index_col=0)
    if "v147" not in df.index:
        continue
    records.append(df.loc["v147"].rename(area))

combined = (pd.DataFrame(records)
              .reindex(columns=celltypes, fill_value=0)
              .sort_index())

print("snRNA areas:", combined.index.tolist())

# ── 2. Spatial projection scores (p36/p147 genes → area labels) ──────────────
indir        = "/maps/projects/perslab/people/qbx224/7plex_paper/atlas/area_cellex_results"
gene_area_map = {
    "spatial1": {"p36-NLS-Cre": "PAG", "p147-mCherry-2A-iCre": "BST"},
    "spatial2": {"p36-NLS-Cre": "PBN", "p147-mCherry-2A-iCre": "LHA"},
    "spatial3": {"p36-NLS-Cre": "CEA", "p147-mCherry-2A-iCre": "PVT"},
}
records_sp = []
for run, gene_map in gene_area_map.items():
    df = pd.read_csv(os.path.join(indir, f"{run}_detailed_cellex.csv"), index_col=0)
    for gene, area in gene_map.items():
        if gene not in df.index:
            raise KeyError(f"{gene!r} not found in {run}_detailed_cellex.csv")
        records_sp.append(df.loc[gene].rename(area))

combined_spatial = (pd.DataFrame(records_sp)
                      .reindex(columns=celltypes, fill_value=0))

print("Spatial areas:", combined_spatial.index.tolist())

# ── 3. Combined projection score (average snRNA + spatial where both exist) ───
single_dataset_regions = {'PVH', 'DVC'}   # snRNA only — do not halve

all_areas = sorted(set(combined.index) | set(combined_spatial.index))
score_df  = pd.DataFrame(index=all_areas, columns=cell_type_order, dtype=float)

for region in all_areas:
    for ct in cell_type_order:
        snrna   = combined.loc[region, ct]         if (region in combined.index
                                                        and ct in combined.columns)         else np.nan
        spatial = combined_spatial.loc[region, ct] if (region in combined_spatial.index
                                                        and ct in combined_spatial.columns) else np.nan
        if np.isnan(snrna) and np.isnan(spatial):
            score_df.loc[region, ct] = np.nan
        elif region in single_dataset_regions:
            score_df.loc[region, ct] = np.nansum([snrna, spatial])
        else:
            score_df.loc[region, ct] = np.nansum([snrna, spatial]) / 2

print("Score matrix:", score_df.shape, "| Areas:", score_df.index.tolist())

# ── 4. Receptor expression specificity from full-atlas cellex ─────────────────
receptor_genes = ['Calcr', 'Ramp1', 'Ramp2', 'Ramp3']
cellex_full    = pd.read_csv("all_cellex_results/adata_cellex.csv", index_col=0)
# rows = genes, cols = cell types
print("Cellex genes available (sample):", cellex_full.index[:5].tolist())
print("Cellex columns (sample):",         cellex_full.columns[:5].tolist())

receptors_present = [g for g in receptor_genes if g in cellex_full.index]
missing           = set(receptor_genes) - set(receptors_present)
if missing:
    print(f"WARNING: {missing} not found in cellex — will be skipped")

receptor_df = (cellex_full.loc[receptors_present]
                           .reindex(columns=cell_type_order, fill_value=0))

# ── 5. DEG-responsive cell types (from earlier analysis) ─────────────────────
responsive_celltypes = {
    '26_Dach2_Npsr1':   29,
    '01_Undefined':     17,
    '19_Bace2_Tspan18': 12,
    '05_Qrfpr_Pdzrn3':   2,
    '13_Pth2r_Otp':      1,
    '12_Ntng2_Onecut1':  1,
    '23_Ust_Nkd2':       1,
    'Astrocyte':          1,
}

# ── 6. Figure ─────────────────────────────────────────────────────────────────
norm      = Normalize(vmin=0, vmax=1)
cmap_rec  = plt.get_cmap('YlOrRd')
cmap_main = plt.get_cmap('rocket_r')
cmap_main.set_bad(color='lightgray')
sig_cmap  = plt.get_cmap('Oranges')
max_sig   = max(responsive_celltypes.values())

n_rec  = len(receptors_present)
n_area = len(score_df.index)

fig = plt.figure(figsize=(20, 4 + n_rec * 0.5))
gs  = GridSpec(2, 1,
               height_ratios=[n_rec, n_area],
               hspace=0.05,
               figure=fig)
ax_rec  = fig.add_subplot(gs[0])
ax_main = fig.add_subplot(gs[1])

# Receptor strip
ax_rec.imshow(receptor_df.values.astype(float),
              aspect='auto', cmap=cmap_rec, norm=norm, interpolation='nearest')
ax_rec.set_xticks(np.arange(-.5, len(cell_type_order), 1), minor=True)
ax_rec.set_yticks(np.arange(-.5, n_rec, 1), minor=True)
ax_rec.grid(which='minor', color='gray', linestyle='-', linewidth=0.5, alpha=0.3)
ax_rec.set_yticks(range(n_rec))
ax_rec.set_yticklabels(receptors_present, fontsize=9, fontstyle='italic')
ax_rec.set_xticks([])
ax_rec.set_title('Amylin receptor components (ESµ)  |  LPBN projection scores',
                 fontsize=11, pad=6)

cbar_rec = plt.colorbar(ScalarMappable(norm=norm, cmap=cmap_rec),
                         ax=ax_rec, fraction=0.015, pad=0.01)
cbar_rec.set_label('ESµ (receptor)', fontsize=8)
cbar_rec.set_ticks([0, 0.5, 1.0])
cbar_rec.ax.tick_params(labelsize=7)

# Main projection heatmap
score_array = score_df.values.astype(float)
ax_main.imshow(np.ma.masked_invalid(score_array),
               aspect='auto', cmap=cmap_main, norm=norm, interpolation='nearest')
ax_main.set_xticks(np.arange(-.5, len(cell_type_order), 1), minor=True)
ax_main.set_yticks(np.arange(-.5, n_area, 1), minor=True)
ax_main.grid(which='minor', color='gray', linestyle='-', linewidth=0.5, alpha=0.3)

for row_i, region in enumerate(score_df.index):
    for col_i, val in enumerate(score_array[row_i]):
        if np.isnan(val) or val <= 0.36:
            continue
        text_color = 'white' if norm(val) > 0.5 else 'black'
        ax_main.text(col_i, row_i, '*', ha='center', va='center',
                     fontsize=10, color=text_color, fontweight='bold')

ax_main.set_yticks(range(n_area))
ax_main.set_yticklabels(score_df.index, fontsize=11)
ax_main.set_xticks(range(len(cell_type_order)))

xticklabels = ax_main.set_xticklabels(cell_type_order, rotation=90,
                                        ha='right', fontsize=8)
for label in xticklabels:
    ct = label.get_text()
    if ct in responsive_celltypes:
        n_sig = responsive_celltypes[ct]
        label.set_color(sig_cmap(0.4 + 0.6 * (n_sig / max_sig)))
        label.set_fontweight('bold')

cbar_main = plt.colorbar(ScalarMappable(norm=norm, cmap=cmap_main),
                          ax=ax_main, fraction=0.015, pad=0.01)
cbar_main.set_label('Combined projection score\n(mean ESµ)', fontsize=9)
cbar_main.set_ticks([0, 0.25, 0.5, 0.75, 1.0])
cbar_main.ax.tick_params(labelsize=8)

legend_patches = [
    mpatches.Patch(color=sig_cmap(0.4 + 0.6 * (29 / max_sig)), label='833 DEGs: high (29)'),
    mpatches.Patch(color=sig_cmap(0.4 + 0.6 * (12 / max_sig)), label='833 DEGs: mid (12)'),
    mpatches.Patch(color=sig_cmap(0.4 + 0.6 * ( 1 / max_sig)), label='833 DEGs: low (1–2)'),
]
ax_main.legend(handles=legend_patches, loc='lower right', fontsize=8,
               title='DEG response to 833', title_fontsize=8, framealpha=0.85)

plt.savefig('region_celltype_833_receptor_projection.pdf', bbox_inches='tight', dpi=300)
plt.show()


# In[ ]:





# In[ ]:




