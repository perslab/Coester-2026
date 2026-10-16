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
from pathlib import Path


# In[2]:


adata = sc.read_h5ad("DVC_integrated_2plex_labelled.h5ad")
adata


# In[3]:


adata.layers["RNA"] = adata.X
adata.layers


# In[4]:


print(adata.X[:5, :5].toarray())#Check if counts are the active layer


# In[5]:


adata.var = pd.DataFrame(index=adata.var_names)
adata


# In[6]:


# Remove specific columns from adata.var
columns_to_remove = ["nUMI_hto", "nCount_SCT", "nFeature_SCT", "SCT_snn_res.0.8", "SCT_snn_res.0.8", "integrated_snn_res.0.8"]  
adata.obs = adata.obs.drop(columns=columns_to_remove)
adata.obs


# In[7]:


sc.pl.embedding(adata,basis="X_umap",color="seurat_clusters", legend_loc="on data", legend_fontsize="xx-small")


# In[8]:


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


# In[9]:


adata.write("DVC_integrated_2plex_anndata.h5ad")
adata


# In[10]:


print(adata.X[:5, :5].toarray())


# In[11]:


adata.layers["counts"] = adata.X
adata


# In[12]:


os.environ['R_HOME'] = '/.conda/envs/pca_env/lib/R'
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


# In[16]:


get_ipython().run_cell_magic('R', '-o size_factors', '\nsize_factors = sizeFactors(\n    computeSumFactors(\n        SingleCellExperiment(\n            list(counts = data_mat)),\n            clusters = NULL,   # let scran pool cells itself\n            min.mean = 0.1,\n            BPPARAM = MulticoreParam()\n    )\n)\n')


# In[17]:


from scipy.sparse import csr_matrix
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


# In[19]:


adata.write("DVC_integrated_2plex_anndata.h5ad")
adata


# In[22]:


# Imports
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
sns.set_style("whitegrid")
plt.rcParams.update({
    'font.size': 9, 'axes.titlesize': 11, 'axes.labelsize': 10,
    'figure.dpi': 150, 'savefig.dpi': 300, 'savefig.bbox': 'tight',
    'pdf.fonttype': 42, 'ps.fonttype': 42,
})

ct_col = 'detailed_celltype_annotation'
score_col = 'score_detailed_celltype_annotation'

#Load data
adata_2p = sc.read("DVC_integrated_2plex_anndata.h5ad")
# Load just the metadata instead of the full h5ad
ref_meta = pd.read_csv("DVC_mouse_ref_metadata.csv", index_col=0)
ref_meta = ref_meta.rename(columns={'cell.type': 'detailed_celltype_annotation'})

print(f"Reference: {len(ref_meta)} cells")
print(f"Columns: {ref_meta.columns.tolist()}")
print(ref_meta[ct_col].value_counts().head(10))

print(f"2-plex: {adata_2p.n_obs} cells × {adata_2p.n_vars} genes")
adata_2p


# In[23]:


# %% Cell 3 — Assign tracing groups from surgery sheet
# VMH ANTERO: AP -1.50, v36-1 anterograde from ventromedial hypothalamus
# MM  ANTERO: AP -2.80, v36-1 anterograde from mammillary bodies
# Exclusions from Section sheet:
#   VMH: 174 (no SAC), 180 (L-, R- failed), 163 (no SAC)
#   MM:  162 (L-, R+ marginal), 161 (no SAC), 167 (no SAC)

vmh_ids      = [166, 169, 177]          # VMH
mm_ids       = [170, 178, 153]          # MM
excluded_ids = [174, 180, 162]          # failures

adata_2p.obs['animal_id'] = (adata_2p.obs['sample_id']
                               .str.extract(r'(\d+)')
                               .astype(int).values)

adata_2p.obs['tracing_group'] = 'other'
adata_2p.obs.loc[adata_2p.obs['animal_id'].isin(vmh_ids),      'tracing_group'] = 'VMH_antero'
adata_2p.obs.loc[adata_2p.obs['animal_id'].isin(mm_ids),       'tracing_group'] = 'MM_antero'
adata_2p.obs.loc[adata_2p.obs['animal_id'].isin(excluded_ids), 'tracing_group'] = 'excluded'

print(adata_2p.obs['tracing_group'].value_counts())
print()
print(pd.crosstab(adata_2p.obs['tracing_group'], adata_2p.obs['sample_id']))


# In[24]:


# %% Cell 4 — Subset to each anterograde group
vmh  = adata_2p[adata_2p.obs['tracing_group'] == 'VMH_antero'].copy()
mm   = adata_2p[adata_2p.obs['tracing_group'] == 'MM_antero'].copy()

print(f"VMH antero: {vmh.n_obs} cells from {vmh.obs['animal_id'].nunique()} animals "
      f"({sorted(vmh.obs['animal_id'].unique())})")
print(f"MM  antero: {mm.n_obs} cells from {mm.obs['animal_id'].nunique()} animals "
      f"({sorted(mm.obs['animal_id'].unique())})")

# %% Cell 5 — Label transfer quality
fig, axes = plt.subplots(1, 2, figsize=(12, 4))
for ax, (name, data) in zip(axes, [('VMH anterograde', vmh), ('MM anterograde', mm)]):
    scores = data.obs[score_col]
    ax.hist(scores, bins=50, color='steelblue', edgecolor='white', linewidth=0.3)
    ax.axvline(0.1, color='firebrick', ls='--', lw=1.5, label='threshold 0.1')
    above = (scores > 0.1).mean() * 100
    ax.set_title(f'{name}\nmedian={scores.median():.3f}, >0.1: {above:.0f}%')
    ax.set_xlabel('Prediction score')
    ax.legend()
plt.tight_layout()
plt.savefig(f'{OUTDIR}/QC_label_transfer_VMH_MM.pdf')
plt.show()


# In[26]:


# %% Cell 6 — Filter low-confidence labels
vmh_hq = vmh[vmh.obs[score_col] > 0.1].copy()
mm_hq  = mm[mm.obs[score_col]  > 0.1].copy()

print(f"VMH after filtering: {vmh_hq.n_obs}/{vmh.n_obs} ({vmh_hq.n_obs/vmh.n_obs*100:.1f}%)")
print(f"MM  after filtering: {mm_hq.n_obs}/{mm.n_obs} ({mm_hq.n_obs/mm.n_obs*100:.1f}%)")

# %% Cell 7 — Enrichment function (reusable)
def compute_enrichment(query_adata, ref_meta, label,
                       min_n=5, ct_col='detailed_celltype_annotation'):
    """Fisher's exact enrichment — ref_meta is a plain DataFrame."""
    ref_cts   = ref_meta[ct_col].value_counts()
    query_cts = query_adata.obs[ct_col].value_counts()
    all_types = sorted(set(ref_cts.index) | set(query_cts.index))

    n_query = len(query_adata)
    n_ref   = len(ref_meta)

    records = []
    for ct in all_types:
        n_ct_q = query_cts.get(ct, 0)
        n_ct_r = ref_cts.get(ct, 0)
        pct_q  = n_ct_q / n_query * 100
        pct_r  = n_ct_r / n_ref   * 100
        fold   = (pct_q / pct_r) if pct_r > 0 else np.inf

        table      = np.array([[n_ct_q, n_query - n_ct_q],
                               [n_ct_r, n_ref   - n_ct_r]])
        odds, pval = fisher_exact(table, alternative='greater')

        records.append({
            'cell_type': ct, 'n_query': n_ct_q, 'pct_query': pct_q,
            'n_atlas': n_ct_r, 'pct_atlas': pct_r,
            'fold_enrichment': fold, 'odds_ratio': odds, 'pvalue': pval,
        })

    df = pd.DataFrame(records)
    _, fdr, _, _ = multipletests(df['pvalue'], method='fdr_bh')
    df['fdr']         = fdr
    df['significant'] = df['fdr'] < 0.05
    df = df.sort_values('fold_enrichment', ascending=False)
    df.to_csv(f'{OUTDIR}/{label}_enrichment_table.csv', index=False)
    return df

df_vmh = compute_enrichment(vmh_hq, ref_meta, 'VMH_antero')
df_mm  = compute_enrichment(mm_hq,  ref_meta, 'MM_antero')


# In[29]:


# %% Cell 8 — Print enrichment results
glial_types = {
    'Astrocytes', 'Oligodendrocytes', 'OPCs', 'Microglia',
    'VLMCs', 'Ependymal_cells', 'Tanycytes', 'Pericytes',
    'Endothelial_cells',
}

# Neurons are anything starting with GABA_, Glu_, or Chat_
# so we can also define it positively as a safety net
def is_neuronal(ct):
    return (ct.startswith('GABA_') or 
            ct.startswith('Glu_') or 
            ct.startswith('Chat_'))

for label, df in [('VMH anterograde', df_vmh), ('MM anterograde', df_mm)]:
    sig = df[
        df['significant'] &
        (df['n_query'] >= 5) &
        (~df['cell_type'].isin(glial_types))
    ]
    print(f"\n{'='*70}")
    print(f"Significantly enriched neurons — {label} (FDR<0.05, n≥5):")
    print(f"{'Cell type':30s} {'n':>6s} {'%_query':>8s} {'%_atlas':>8s} "
          f"{'fold':>6s} {'FDR':>10s}")
    print("-" * 70)
    for _, r in sig.iterrows():
        print(f"{r['cell_type']:30s} {r['n_query']:>6.0f} {r['pct_query']:>7.1f}% "
              f"{r['pct_atlas']:>7.1f}% {r['fold_enrichment']:>6.2f} "
              f"{r['fdr']:>10.2e}")


# In[33]:


# %% Cell 9 — Enrichment bar chart (neurons only, both groups side by side)

fig, axes = plt.subplots(1, 2, figsize=(14, 9), sharey=False)

for ax, (label, df, data) in zip(axes, [
    ('VMH anterograde', df_vmh, vmh_hq),
    ('MM anterograde',  df_mm,  mm_hq),
]):
    exclude = glial_types | {'Undefined', '01_Undefined'}

    df_plot = (df[~df['cell_type'].isin(exclude)]
               .nlargest(20, 'n_query')
               .sort_values('pct_query', ascending=True))

    colors = ['#E8A87C' if r['significant'] else '#C0C0C0' for _, r in df_plot.iterrows()]

    y = range(len(df_plot))
    ax.barh(y, df_plot['pct_query'].values, color=colors,
            edgecolor='white', linewidth=0.5, height=0.7)

    for i, (_, r) in enumerate(df_plot.iterrows()):
        ax.plot(r['pct_atlas'], i, 'k|', markersize=8, markeredgewidth=1.5, zorder=5)

    for i, (_, r) in enumerate(df_plot.iterrows()):
        if r['significant'] and r['fold_enrichment'] > 1:
            ax.text(r['pct_query'] + 0.3, i, f'{r["fold_enrichment"]:.1f}×',
                    va='center', fontsize=7, fontweight='bold',
                    color='#8B0000')

    ax.set_yticks(y)
    ax.set_yticklabels(df_plot['cell_type'].values, fontsize=7)
    ax.set_xlabel('% of sorted GFP+ nuclei')
    n_animals = data.obs['animal_id'].nunique()
    ax.set_title(f'Cell type composition\n{label}\n'
                 f'(n={data.n_obs} cells, {n_animals} animals)',
                 fontweight='bold')



    legend_el = [
        Patch(fc='#E8A87C', label='FDR<0.05'),
        Patch(fc='#C0C0C0', label='Not enriched'),
        plt.Line2D([0], [0], color='k', marker='|', ls='None',
                   markersize=8, markeredgewidth=1.5, label='Atlas proportion'),
    ]
    ax.legend(handles=legend_el, fontsize=7, loc='lower right')
    ax.grid(axis='x', alpha=0.3)

plt.suptitle('DVC cell types receiving hypothalamic input\n'
             'VMH (anterior) vs MM (mammillary) anterograde tracing',
             fontweight='bold', fontsize=12)
plt.tight_layout()
plt.savefig(f'{OUTDIR}/Fig_VMH_MM_enrichment_bar.pdf', bbox_inches='tight')
plt.show()


# In[40]:


# %% Cell 10b — Enrichment bar chart including glia
fig, axes = plt.subplots(1, 2, figsize=(14, 12), sharey=False)

for ax, (label, df, data) in zip(axes, [
    ('VMH anterograde', df_vmh, vmh_hq),
    ('MM anterograde',  df_mm,  mm_hq),
]):
    # No glia exclusion — only drop Undefined populations
    df_plot = (df[~df['cell_type'].isin({'Undefined', '01_Undefined'})]
               .nlargest(30, 'n_query')
               .sort_values('pct_query', ascending=True))

    colors = []
    for _, r in df_plot.iterrows():
        if r['cell_type'] in glial_types and r['significant']:
            colors.append('#4878CF')   # blue = significant glia
        elif r['cell_type'] in glial_types:
            colors.append('#A8C4E8')   # light blue = non-significant glia
        elif r['significant']:
            colors.append('#E8A87C')   # orange = significant neuron
        else:
            colors.append('#C0C0C0')   # grey = non-significant neuron

    y = range(len(df_plot))
    ax.barh(y, df_plot['pct_query'].values, color=colors,
            edgecolor='white', linewidth=0.5, height=0.7)

    for i, (_, r) in enumerate(df_plot.iterrows()):
        ax.plot(r['pct_atlas'], i, 'k|', markersize=8, markeredgewidth=1.5, zorder=5)

    for i, (_, r) in enumerate(df_plot.iterrows()):
        if r['significant'] and r['fold_enrichment'] > 1:
            ax.text(r['pct_query'] + 0.3, i, f'{r["fold_enrichment"]:.1f}×',
                    va='center', fontsize=7, fontweight='bold', color='#444')

    ax.set_yticks(y)
    ax.set_yticklabels(df_plot['cell_type'].values, fontsize=7)
    ax.set_xlabel('% of sorted GFP+ nuclei')
    n_animals = data.obs['animal_id'].nunique()
    ax.set_title(f'Cell type composition\n{label}\n'
                 f'(n={data.n_obs} cells, {n_animals} animals)',
                 fontweight='bold')

    legend_el = [
        Patch(fc='#4878CF', label='Glia + FDR<0.05'),
        Patch(fc='#A8C4E8', label='Glia (n.s.)'),
        Patch(fc='#E8A87C', label='Neuron + FDR<0.05'),
        Patch(fc='#C0C0C0', label='Neuron (n.s.)'),
        plt.Line2D([0], [0], color='k', marker='|', ls='None',
                   markersize=8, markeredgewidth=1.5, label='Atlas proportion'),
    ]
    ax.legend(handles=legend_el, fontsize=7, loc='lower right')
    ax.grid(axis='x', alpha=0.3)

plt.suptitle('DVC cell types receiving hypothalamic input (all cell classes)\n'
             'VMH (anterior) vs MM (mammillary) anterograde tracing',
             fontweight='bold', fontsize=12)
plt.tight_layout()
plt.savefig(f'{OUTDIR}/Fig_VMH_MM_enrichment_bar_withglia.pdf', bbox_inches='tight')
plt.show()


# In[45]:


# %% Cell 3 (revised) — Assign tracing groups
vmh_ids      = [166, 169, 177]
mm_ids       = [170, 178, 153]
excluded_ids = [180, 162, 174]   # failed injections → negative control

adata_2p.obs['animal_id'] = (adata_2p.obs['sample_id']
                               .str.extract(r'(\d+)')
                               .astype(int).values)

adata_2p.obs['tracing_group'] = 'other'
adata_2p.obs.loc[adata_2p.obs['animal_id'].isin(vmh_ids),      'tracing_group'] = 'VMH_antero'
adata_2p.obs.loc[adata_2p.obs['animal_id'].isin(mm_ids),       'tracing_group'] = 'MM_antero'
adata_2p.obs.loc[adata_2p.obs['animal_id'].isin(excluded_ids), 'tracing_group'] = 'failed_ctrl'

print(adata_2p.obs['tracing_group'].value_counts())
print()
print(pd.crosstab(adata_2p.obs['tracing_group'], adata_2p.obs['sample_id']))

# %% Cell 4 (revised) — Subset
vmh      = adata_2p[adata_2p.obs['tracing_group'] == 'VMH_antero'].copy()
mm       = adata_2p[adata_2p.obs['tracing_group'] == 'MM_antero'].copy()
ctrl     = adata_2p[adata_2p.obs['tracing_group'] == 'failed_ctrl'].copy()

for name, data in [('VMH', vmh), ('MM', mm), ('Failed ctrl', ctrl)]:
    print(f"{name}: {data.n_obs} cells from "
          f"{data.obs['animal_id'].nunique()} animals "
          f"({sorted(data.obs['animal_id'].unique())})")

# %% Cell 7 (revised) — Enrichment vs failed controls
def compute_enrichment(query_adata, ref_adata_or_meta, label,
                       min_n=5, ct_col='detailed_celltype_annotation'):
    """Fisher's exact enrichment — ref can be AnnData or plain DataFrame."""
    if hasattr(ref_adata_or_meta, 'obs'):
        ref_obs = ref_adata_or_meta.obs
    else:
        ref_obs = ref_adata_or_meta

    ref_cts   = ref_obs[ct_col].value_counts()
    query_cts = query_adata.obs[ct_col].value_counts()
    all_types = sorted(set(ref_cts.index) | set(query_cts.index))

    n_query = len(query_adata)
    n_ref   = len(ref_obs)

    records = []
    for ct in all_types:
        n_ct_q = query_cts.get(ct, 0)
        n_ct_r = ref_cts.get(ct, 0)
        pct_q  = n_ct_q / n_query * 100
        pct_r  = n_ct_r / n_ref   * 100
        fold   = (pct_q / pct_r) if pct_r > 0 else np.inf

        table      = np.array([[n_ct_q, n_query - n_ct_q],
                               [n_ct_r, n_ref   - n_ct_r]])
        odds, pval = fisher_exact(table, alternative='greater')

        records.append({
            'cell_type': ct, 'n_query': n_ct_q, 'pct_query': pct_q,
            'n_ctrl':  n_ct_r,  'pct_ctrl':  pct_r,
            'fold_enrichment': fold, 'odds_ratio': odds, 'pvalue': pval,
        })

    df = pd.DataFrame(records)
    _, fdr, _, _ = multipletests(df['pvalue'], method='fdr_bh')
    df['fdr']         = fdr
    df['significant'] = df['fdr'] < 0.05
    df = df.sort_values('fold_enrichment', ascending=False)
    df.to_csv(f'{OUTDIR}/{label}_enrichment_vs_ctrl.csv', index=False)
    print(f"Saved {label}_enrichment_vs_ctrl.csv")
    return df

# Use failed controls as reference instead of atlas
ctrl_hq = ctrl[ctrl.obs[score_col] > 0.1].copy()
vmh_hq  = vmh[vmh.obs[score_col]   > 0.1].copy()
mm_hq   = mm[mm.obs[score_col]    > 0.1].copy()

print(f"\nFailed ctrl: {ctrl_hq.n_obs}/{ctrl.n_obs} pass QC")
print(f"VMH:         {vmh_hq.n_obs}/{vmh.n_obs} pass QC")
print(f"MM:          {mm_hq.n_obs}/{mm.n_obs} pass QC")

df_vmh = compute_enrichment(vmh_hq, ctrl_hq, 'VMH_antero_vs_ctrl')
df_mm  = compute_enrichment(mm_hq,  ctrl_hq, 'MM_antero_vs_ctrl')

# %% Cell 8 — Print results
for label, df in [('VMH vs failed ctrl', df_vmh), ('MM vs failed ctrl', df_mm)]:
    sig = df[
        df['significant'] &
        (df['n_query'] >= 5) &
        (~df['cell_type'].isin(glial_types))
    ]
    print(f"\n{'='*70}")
    print(f"Enriched neurons — {label} (FDR<0.05, n≥5):")
    print(f"{'Cell type':30s} {'n_query':>8s} {'%_query':>8s} "
          f"{'%_ctrl':>8s} {'fold':>6s} {'FDR':>10s}")
    print("-" * 70)
    for _, r in sig.iterrows():
        print(f"{r['cell_type']:30s} {r['n_query']:>8.0f} {r['pct_query']:>7.1f}% "
              f"{r['pct_ctrl']:>7.1f}% {r['fold_enrichment']:>6.2f} "
              f"{r['fdr']:>10.2e}")

# %% Cell 9 — Bar chart vs failed controls (neurons + glia, two versions)
for include_glia, suffix in [(False, 'neurons'), (True, 'all')]:
    fig, axes = plt.subplots(1, 2, figsize=(14, 10 if include_glia else 9),
                              sharey=False)

    for ax, (label, df, data) in zip(axes, [
        ('VMH anterograde', df_vmh, vmh_hq),
        ('MM anterograde',  df_mm,  mm_hq),
    ]):
        exclude = {'Undefined', '01_Undefined'}
        if not include_glia:
            exclude |= glial_types

        df_plot = (df[~df['cell_type'].isin(exclude)]
                   .nlargest(25, 'n_query')
                   .sort_values('pct_query', ascending=True))

        colors = []
        for _, r in df_plot.iterrows():
            is_glia = r['cell_type'] in glial_types
            if is_glia and r['significant']:      colors.append('#4878CF')
            elif is_glia:                          colors.append('#A8C4E8')
            elif r['significant']:                 colors.append('#E8A87C')
            else:                                  colors.append('#C0C0C0')

        y = range(len(df_plot))
        ax.barh(y, df_plot['pct_query'].values, color=colors,
                edgecolor='white', linewidth=0.5, height=0.7)

        # Failed control proportion as black ticks
        for i, (_, r) in enumerate(df_plot.iterrows()):
            ax.plot(r['pct_ctrl'], i, 'k|',
                    markersize=8, markeredgewidth=1.5, zorder=5)

        for i, (_, r) in enumerate(df_plot.iterrows()):
            if r['significant'] and r['fold_enrichment'] > 1:
                ax.text(r['pct_query'] + 0.3, i,
                        f'{r["fold_enrichment"]:.1f}×',
                        va='center', fontsize=7, fontweight='bold', color='#444')

        ax.set_yticks(y)
        ax.set_yticklabels(df_plot['cell_type'].values, fontsize=7)
        ax.set_xlabel('% of sorted GFP+ nuclei')
        ax.set_title(f'{label}\nvs failed injection controls\n'
                     f'(n={data.n_obs} cells, '
                     f'{data.obs["animal_id"].nunique()} animals)',
                     fontweight='bold')

        legend_el = [
            Patch(fc='#E8A87C', label='Neuron FDR<0.05'),
            Patch(fc='#C0C0C0', label='Neuron n.s.'),
            *([ Patch(fc='#4878CF', label='Glia FDR<0.05'),
                Patch(fc='#A8C4E8', label='Glia n.s.')] if include_glia else []),
            plt.Line2D([0], [0], color='k', marker='|', ls='None',
                       markersize=8, markeredgewidth=1.5,
                       label='Failed ctrl proportion'),
        ]
        ax.legend(handles=legend_el, fontsize=7, loc='lower right')
        ax.grid(axis='x', alpha=0.3)

    plt.suptitle(
        f'Cell type enrichment vs failed injection controls\n'
        f'VMH vs MM anterograde tracing '
        f'({"neurons only" if not include_glia else "all cell types"})\n'
        f'Controls: animals {excluded_ids}',
        fontweight='bold', fontsize=11)
    plt.tight_layout()
    plt.savefig(f'{OUTDIR}/Fig_VMH_MM_vs_ctrl_{suffix}.pdf', bbox_inches='tight')
    plt.show()


# In[46]:


# %% Cell — Pooled anterograde (VMH + MM) vs failed controls

import anndata
pooled_hq = anndata.concat([vmh_hq, mm_hq], label='source',
                            keys=['VMH', 'MM'], join='outer')

print(f"Pooled anterograde: {pooled_hq.n_obs} cells from "
      f"{pooled_hq.obs['animal_id'].nunique()} animals "
      f"({sorted(pooled_hq.obs['animal_id'].unique())})")
print(f"Failed ctrl:        {ctrl_hq.n_obs} cells from "
      f"{ctrl_hq.obs['animal_id'].nunique()} animals "
      f"({sorted(ctrl_hq.obs['animal_id'].unique())})")
print(f"\nPooled source breakdown:")
print(pooled_hq.obs['source'].value_counts())

df_pooled = compute_enrichment(pooled_hq, ctrl_hq, 'pooled_antero_vs_ctrl')

# ── Print results ─────────────────────────────────────────────────────────────
for include_glia, label in [(False, 'neurons only'), (True, 'all cell types')]:
    exclude = {'Undefined', '01_Undefined'}
    if not include_glia:
        exclude |= glial_types
    sig = df_pooled[
        df_pooled['significant'] &
        (df_pooled['n_query'] >= 5) &
        (~df_pooled['cell_type'].isin(exclude))
    ]
    print(f"\n{'='*70}")
    print(f"Enriched in pooled anterograde vs failed ctrl ({label}, FDR<0.05, n≥5):")
    print(f"{'Cell type':30s} {'n_query':>8s} {'%_query':>8s} "
          f"{'%_ctrl':>8s} {'fold':>6s} {'FDR':>10s}")
    print("-" * 70)
    for _, r in sig.iterrows():
        print(f"{r['cell_type']:30s} {r['n_query']:>8.0f} {r['pct_query']:>7.1f}% "
              f"{r['pct_ctrl']:>7.1f}% {r['fold_enrichment']:>6.2f} "
              f"{r['fdr']:>10.2e}")

# ── Bar chart ─────────────────────────────────────────────────────────────────
for include_glia, suffix in [(False, 'neurons'), (True, 'all')]:
    fig, ax = plt.subplots(figsize=(7, 10))

    exclude = {'Undefined', '01_Undefined'}
    if not include_glia:
        exclude |= glial_types

    df_plot = (df_pooled[~df_pooled['cell_type'].isin(exclude)]
               .nlargest(25, 'n_query')
               .sort_values('pct_query', ascending=True))

    colors = []
    for _, r in df_plot.iterrows():
        is_glia = r['cell_type'] in glial_types
        if is_glia and r['significant']:  colors.append('#4878CF')
        elif is_glia:                      colors.append('#A8C4E8')
        elif r['significant']:             colors.append('#E8A87C')
        else:                              colors.append('#C0C0C0')

    y = range(len(df_plot))
    ax.barh(y, df_plot['pct_query'].values, color=colors,
            edgecolor='white', linewidth=0.5, height=0.7)

    for i, (_, r) in enumerate(df_plot.iterrows()):
        ax.plot(r['pct_ctrl'], i, 'k|',
                markersize=8, markeredgewidth=1.5, zorder=5)

    for i, (_, r) in enumerate(df_plot.iterrows()):
        if r['significant'] and r['fold_enrichment'] > 1:
            ax.text(r['pct_query'] + 0.3, i,
                    f'{r["fold_enrichment"]:.1f}×',
                    va='center', fontsize=7, fontweight='bold', color='#444')

    ax.set_yticks(y)
    ax.set_yticklabels(df_plot['cell_type'].values, fontsize=7)
    ax.set_xlabel('% of sorted GFP+ nuclei')
    ax.set_title(f'Pooled VMH + MM anterograde\nvs failed injection controls\n'
                 f'(n={pooled_hq.n_obs} cells, '
                 f'{pooled_hq.obs["animal_id"].nunique()} animals)',
                 fontweight='bold')

    legend_el = [
        Patch(fc='#E8A87C', label='Neuron FDR<0.05'),
        Patch(fc='#C0C0C0', label='Neuron n.s.'),
        *([ Patch(fc='#4878CF', label='Glia FDR<0.05'),
            Patch(fc='#A8C4E8', label='Glia n.s.')] if include_glia else []),
        plt.Line2D([0], [0], color='k', marker='|', ls='None',
                   markersize=8, markeredgewidth=1.5,
                   label='Failed ctrl proportion'),
    ]
    ax.legend(handles=legend_el, fontsize=7, loc='lower right')
    ax.grid(axis='x', alpha=0.3)

    plt.suptitle(
        f'Cell type enrichment: pooled HYP anterograde vs failed controls\n'
        f'VMH animals {vmh_ids} + MM animals {mm_ids}\n'
        f'Controls: {excluded_ids}',
        fontweight='bold', fontsize=11)
    plt.tight_layout()
    plt.savefig(f'{OUTDIR}/Fig_pooled_antero_vs_ctrl_{suffix}.pdf', bbox_inches='tight')
    plt.show()


# In[ ]:




