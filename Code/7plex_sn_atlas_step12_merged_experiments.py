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


# In[4]:


#adata = sc.read("allgenes_bipca_hypomap.h5ad")
adata_atlas = sc.read("allgenes_bipca_annotated.h5ad", cache=True)
adata_atlas


# In[3]:


ax = sc.pl.embedding(
    adata_atlas,
    basis="tsne",
    color="detailed_celltype_annotation",
    legend_loc="on data",
    legend_fontsize="x-small",
    return_fig=True)
ax=plt.gca()
texts = ax.texts
adjustText.adjust_text(texts, 
                       ax=ax,
                       arrowprops=dict(arrowstyle="-", color='black', lw=0.5))


# In[5]:


import os
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
import anndata as ad
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from pathlib import Path


# ── 1. Load all three objects ─────────────────────────────────────────────────
adata_2p    = sc.read_h5ad("HYP_integrated_2plex_scran.h5ad")
adata_lpbn  = sc.read_h5ad("HYP_integrated_LPBN_scran.h5ad")

print(f"Atlas:      {adata_atlas.n_obs} cells")
print(f"2-plex:     {adata_2p.n_obs} cells")
print(f"LPBN:       {adata_lpbn.n_obs} cells")

# ── 2. Filter 2-plex to NTS retrograde animals only ───────────────────────────
nts_retro_ids = [168, 176, 149, 159, 165, 188]

adata_2p.obs['animal_id'] = (adata_2p.obs['sample_id']
                               .str.extract(r'(\d+)')
                               .astype(int).values)

nts_retro = adata_2p[adata_2p.obs['animal_id'].isin(nts_retro_ids)].copy()
print(f"NTS retrograde cells: {nts_retro.n_obs} "
      f"from animals {sorted(nts_retro.obs['animal_id'].unique())}")

# All LPBN animals are retrograde
lpbn_retro = adata_lpbn.copy()
n_lpbn_animals = lpbn_retro.obs['sample_id'].str.extract(r'(\d+)')[0].nunique()
print(f"LPBN retrograde cells: {lpbn_retro.n_obs} from {n_lpbn_animals} animals")

# ── 3. Tag dataset source ─────────────────────────────────────────────────────
adata_atlas.obs['dataset'] = 'Atlas (7-plex)'
nts_retro.obs['dataset']   = 'NTS retrograde (v147)'
lpbn_retro.obs['dataset']  = 'LPBN retrograde (v36)'


# In[6]:


# ══════════════════════════════════════════════════════════════════════════════
# STEP 4 — Unify tracer genes and build clean source objects
# ══════════════════════════════════════════════════════════════════════════════

shared_genes = adata_atlas.var_names

def unify_tracer(adata, src_gene='v36', dst_gene='v147'):
    """
    Merge src_gene into dst_gene across all layers and .X, then drop src_gene.
    Handles three cases:
      - both present → sum into dst, drop src
      - only src present → rename src to dst
      - src absent → no-op
    """
    has_src = src_gene in adata.var_names
    has_dst = dst_gene in adata.var_names

    if not has_src:
        print(f"  {src_gene} not found — no action needed")
        return adata

    if has_src and has_dst:
        print(f"  Both {src_gene} and {dst_gene} found — summing into {dst_gene}")
        idx_src = list(adata.var_names).index(src_gene)
        idx_dst = list(adata.var_names).index(dst_gene)

        def _sum_cols(mat, idx_dst, idx_src):
            if sparse.issparse(mat):
                mat = mat.tolil()
                mat[:, idx_dst] = mat[:, idx_dst] + mat[:, idx_src]
                return mat.tocsr()
            else:
                mat[:, idx_dst] += mat[:, idx_src]
                return mat

        for layer_name in list(adata.layers.keys()):
            adata.layers[layer_name] = _sum_cols(
                adata.layers[layer_name], idx_dst, idx_src)

        adata.X = _sum_cols(adata.X, idx_dst, idx_src)
        adata    = adata[:, adata.var_names != src_gene].copy()
        print(f"  Dropped {src_gene}")

    elif has_src and not has_dst:
        print(f"  Only {src_gene} found — renaming to {dst_gene}")
        adata.var_names = adata.var_names.str.replace(
            f'^{src_gene}$', dst_gene, regex=True)

    # Verify no duplicates remain
    assert not adata.var_names.duplicated().any(), \
        "Duplicate var_names after tracer unification!"

    n_pos = (adata[:, dst_gene].layers['counts'] > 0).sum() \
            if 'counts' in adata.layers \
            else (adata[:, dst_gene].X > 0).sum()
    print(f"  {dst_gene} positive cells: {n_pos}/{adata.n_obs} "
          f"({n_pos/adata.n_obs*100:.1f}%)")
    return adata

# ── Apply to each source object ───────────────────────────────────────────────
print("2-plex (NTS dataset):")
adata_2p   = unify_tracer(adata_2p)

print("LPBN retrograde:")
lpbn_retro = unify_tracer(lpbn_retro)

# ── Subset NTS AFTER tracer unification so v147 is already merged ─────────────
nts_retro_ids = [168, 176, 149, 159, 165, 188]
adata_2p.obs['animal_id'] = (adata_2p.obs['sample_id']
                               .str.extract(r'(\d+)')[0]
                               .astype(int).values)
nts_retro = adata_2p[adata_2p.obs['animal_id'].isin(nts_retro_ids)].copy()
nts_retro.obs['dataset'] = 'NTS retrograde (v147)'
lpbn_retro.obs['dataset'] = 'LPBN retrograde (v36)'
adata_atlas.obs['dataset'] = 'Atlas (7-plex)'

print(f"\nNTS retrograde: {nts_retro.n_obs} cells from "
      f"animals {sorted(nts_retro.obs['animal_id'].unique())}")

# ── Final gene state check ────────────────────────────────────────────────────
print("\nFinal tracer gene status:")
for name, adata in [('Atlas', adata_atlas), ('NTS', nts_retro), ('LPBN', lpbn_retro)]:
    has_v36  = 'v36'  in adata.var_names
    has_v147 = 'v147' in adata.var_names
    dupes    = adata.var_names.duplicated().any()
    print(f"  {name:6s}: v36={has_v36}, v147={has_v147}, duplicates={dupes}")

# ── Align to atlas gene space ─────────────────────────────────────────────────
def align_to_atlas(adata, shared_genes, layer='counts'):
    common   = adata.var_names.isin(shared_genes)
    sub      = adata[:, common].copy()
    x        = sub.layers[layer] if layer in sub.layers else sub.X
    if not sparse.issparse(x):
        x = sparse.csr_matrix(x)
    full     = sparse.lil_matrix((sub.n_obs, len(shared_genes)), dtype=np.float32)
    gene_idx = [list(shared_genes).index(g) for g in sub.var_names]
    full[:, gene_idx] = x
    return sparse.csr_matrix(full)

print("\nAligning to atlas gene space...")
counts_atlas = adata_atlas.layers['counts']
if not sparse.issparse(counts_atlas):
    counts_atlas = sparse.csr_matrix(counts_atlas)
counts_nts   = align_to_atlas(nts_retro,  shared_genes)
counts_lpbn  = align_to_atlas(lpbn_retro, shared_genes)

merged_counts = sparse.vstack([counts_atlas, counts_nts, counts_lpbn]).tocsr()
print(f"Merged counts: {merged_counts.shape}")


# In[15]:


# ── 5. Concatenate into one AnnData ──────────────────────────────────────────
# Keep minimal obs columns that exist across all three
keep_obs = ['detailed_celltype_annotation', 'dataset', 'sex']

def trim_obs(adata, keep_cols):
    cols = [c for c in keep_cols if c in adata.obs.columns]
    return adata.obs[cols].copy()

obs_atlas = trim_obs(adata_atlas, keep_obs)
obs_nts   = trim_obs(nts_retro,   keep_obs)
obs_lpbn  = trim_obs(lpbn_retro,  keep_obs)

# Make unique cell barcodes
obs_atlas.index = 'atlas_'  + obs_atlas.index
obs_nts.index   = 'nts_'    + obs_nts.index
obs_lpbn.index  = 'lpbn_'   + obs_lpbn.index

merged_counts = sparse.vstack([counts_atlas, counts_nts, counts_lpbn]).tocsr()
merged_obs    = pd.concat([obs_atlas, obs_nts, obs_lpbn])

adata_merged = ad.AnnData(
    X   = merged_counts,
    obs = merged_obs,
    var = adata_atlas.var[[]].copy()   # keep var index, drop all var columns
)

print(f"\nMerged object: {adata_merged.n_obs} cells × {adata_merged.n_vars} genes")
print(adata_merged.obs['dataset'].value_counts())


# In[16]:


# ── Project via cell type centroid + jitter ───────────────────────────────────
# For each new cell, sample a t-SNE position from atlas cells of the same
# detailed_celltype_annotation label. This is honest because:
#   (1) we already trust the label transfer
#   (2) it places cells exactly where their cell type lives in the atlas
#   (3) jitter prevents overplotting while preserving cluster structure
# Extract atlas t-SNE coordinates
tsne_atlas = adata_atlas.obsm['tsne']
print(f"Atlas t-SNE shape: {tsne_atlas.shape}")
print(f"t-SNE range: x=[{tsne_atlas[:,0].min():.1f}, {tsne_atlas[:,0].max():.1f}], "
      f"y=[{tsne_atlas[:,1].min():.1f}, {tsne_atlas[:,1].max():.1f}]")

# Build per-cell-type t-SNE coordinate pools from atlas
print("Building atlas cell type t-SNE coordinate pools...")
ct_tsne_pool = {}
for ct, grp in adata_atlas.obs.groupby(ct_col, observed=True):
    idx = adata_atlas.obs_names.isin(grp.index)
    ct_tsne_pool[ct] = tsne_atlas[idx]
ct_col = 'detailed_celltype_annotation'

# Build per-cell-type t-SNE coordinate pool from atlas
print("Building atlas cell type t-SNE coordinate pools...")
ct_tsne_pool = {}
for ct, grp in adata_atlas.obs.groupby(ct_col, observed=True):
    idx = adata_atlas.obs_names.isin(grp.index)
    ct_tsne_pool[ct] = tsne_atlas[idx]   # all atlas t-SNE coords for this ct

def assign_tsne_by_celltype(obs_df, ct_tsne_pool, tsne_atlas, rng_seed=42):
    """
    For each cell, sample a t-SNE coordinate from its cell type pool.
    Falls back to global random atlas cell if label not in pool.
    """
    rng    = np.random.default_rng(rng_seed)
    coords = np.zeros((len(obs_df), 2), dtype=np.float32)

    # Compute per-cell-type jitter scale (proportional to cluster spread)
    for i, (_, row) in enumerate(obs_df.iterrows()):
        ct   = row[ct_col]
        pool = ct_tsne_pool.get(ct, None)

        if pool is None or len(pool) == 0:
            # Fallback: random atlas cell
            j = rng.integers(len(tsne_atlas))
            coords[i] = tsne_atlas[j]
        else:
            # Sample from pool + small jitter (5% of cluster std)
            j      = rng.integers(len(pool))
            spread = pool.std(axis=0).mean() * 0.15
            jitter = rng.normal(0, spread, size=2)
            coords[i] = pool[j] + jitter

    return coords

print("Assigning t-SNE coordinates to NTS retrograde cells...")
tsne_nts  = assign_tsne_by_celltype(nts_retro.obs,  ct_tsne_pool, tsne_atlas)

print("Assigning t-SNE coordinates to LPBN retrograde cells...")
tsne_lpbn = assign_tsne_by_celltype(lpbn_retro.obs, ct_tsne_pool, tsne_atlas,
                                     rng_seed=43)

# Stack into merged obsm
tsne_merged = np.vstack([tsne_atlas, tsne_nts, tsne_lpbn])
adata_merged.obsm['X_tsne'] = tsne_merged

# Sanity check — do centroids match?
print("\nSanity check — mean t-SNE coords per dataset:")
for name, coords in [('Atlas',   tsne_atlas),
                     ('NTS',     tsne_nts),
                     ('LPBN',    tsne_lpbn)]:
    print(f"  {name:8s}: x={coords[:,0].mean():.1f}, y={coords[:,1].mean():.1f}, "
          f"range x=[{coords[:,0].min():.1f}, {coords[:,0].max():.1f}]")


# In[17]:


# ── Replot ────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(18, 6))

dataset_colors = {
    'Atlas (7-plex)':        '#D0D0D0',
    'NTS retrograde (v147)': '#C0392B',
    'LPBN retrograde (v36)': '#2980B9',
}

for ax, highlight in zip(axes, [None,
                                  'NTS retrograde (v147)',
                                  'LPBN retrograde (v36)']):
    mask_atlas = adata_merged.obs['dataset'] == 'Atlas (7-plex)'
    ax.scatter(adata_merged.obsm['X_tsne'][mask_atlas, 0],
               adata_merged.obsm['X_tsne'][mask_atlas, 1],
               c='#D0D0D0', s=0.3, alpha=0.3, rasterized=True, zorder=1)

    for ds, color in dataset_colors.items():
        if ds == 'Atlas (7-plex)':
            continue
        if highlight is not None and ds != highlight:
            continue
        mask = adata_merged.obs['dataset'].values == ds
        ax.scatter(adata_merged.obsm['X_tsne'][mask, 0],
                   adata_merged.obsm['X_tsne'][mask, 1],
                   c=color, s=1.5, alpha=0.6, rasterized=True,
                   zorder=2, label=ds)

    title = 'All datasets' if highlight is None else highlight
    ax.set_title(title, fontweight='bold', fontsize=10)
    ax.set_xlabel('t-SNE 1'); ax.set_ylabel('t-SNE 2')
    ax.legend(fontsize=7, markerscale=5, loc='lower right')
    ax.set_xticks([]); ax.set_yticks([])
    ax.spines[['top', 'right', 'left', 'bottom']].set_visible(False)

plt.suptitle('Atlas t-SNE: NTS and LPBN retrograde overlaid\n'
             '(positioned by label-transferred cell type)',
             fontweight='bold', fontsize=11)
plt.tight_layout()
plt.savefig("Fig_atlas_overlay_NTS_LPBN_celltype.pdf", bbox_inches='tight', dpi=300)
plt.show()


# In[18]:


# ── 7. Save merged object for CELLEX ─────────────────────────────────────────
adata_merged.write_h5ad("atlas_merged_NTS_LPBN_retro.h5ad")
adata_merged


# In[19]:


# ══════════════════════════════════════════════════════════════════════════════
# BUILD MERGED OBS — single clean pipeline
# ══════════════════════════════════════════════════════════════════════════════

ct_col = 'detailed_celltype_annotation'

# ── 1. Patch missing obs onto source objects before building ──────────────────

# LPBN: sex from surgery sheet
lpbn_surgery_sex = {
    50: 'female', 52: 'female', 54: 'female', 55: 'female',
    56: 'male',   58: 'male',   59: 'female', 61: 'female',
    62: 'male',   63: 'female', 67: 'female', 68: 'male',
    71: 'female', 73: 'female'
}
lpbn_retro.obs['animal_id_int'] = (lpbn_retro.obs['sample_id']
                                    .str.extract(r'(\d+)')[0].astype(int))
lpbn_retro.obs['sex']  = lpbn_retro.obs['animal_id_int'].map(lpbn_surgery_sex)
lpbn_retro.obs['area'] = 'PBN'

# ── 2. AAV and AAV_ratio (computed from counts layer) ─────────────────────────

def compute_aav(adata, gene='v147', ct_col='detailed_celltype_annotation'):
    if gene not in adata.var_names:
        print(f"  WARNING: {gene} not found — AAV=False, ratio=0")
        return (pd.Series(False, index=adata.obs_names, name='AAV'),
                pd.Series(0.0,   index=adata.obs_names, name='AAV_ratio'))
    x    = adata[:, gene].layers['counts'] if 'counts' in adata.layers \
           else adata[:, gene].X
    expr = x.toarray().flatten() if sparse.issparse(x) else np.array(x).flatten()
    aav  = pd.Series(expr > 0, index=adata.obs_names, name='AAV')
    ratio_map = (adata.obs.assign(AAV=aav.values)
                 .groupby(ct_col, observed=True)['AAV']
                 .mean().rename('AAV_ratio'))
    ratio = adata.obs[ct_col].map(ratio_map).rename('AAV_ratio')
    print(f"  {gene}: {aav.sum()}/{len(aav)} positive ({aav.mean()*100:.1f}%)")
    return aav, ratio

print("Computing AAV positivity:")
print("  Atlas:")
aav_atlas,  ratio_atlas  = compute_aav(adata_atlas)
print("  NTS retrograde:")
aav_nts,    ratio_nts    = compute_aav(nts_retro)
print("  LPBN retrograde:")
aav_lpbn,   ratio_lpbn   = compute_aav(lpbn_retro)

for adata, aav, ratio in [(adata_atlas, aav_atlas,  ratio_atlas),
                           (nts_retro,   aav_nts,    ratio_nts),
                           (lpbn_retro,  aav_lpbn,   ratio_lpbn)]:
    adata.obs['AAV']       = aav.values
    adata.obs['AAV_ratio'] = ratio.values

# ── 3. Cell-type-linked lookup table from atlas ───────────────────────────────
# Columns that are a property of the cell type label, not the individual cell

ct_linked_cols = [
    'neurotransmitter', 'predicted_spatial_region',
    'predicted_spatial_region_norm', 'predicted_spatial_region_filtered',
    'predicted_spatial_region_filtered_norm',
    'cellect_bmi_detailed', 'cellect_bmi_detailed_BH', 'cellect_bmi_log10',
    'detailed_significant_BH',
    # Non-1:1 but useful to carry over via majority vote
    'manual_celltype_annotation', 'predicted_brain_region',
    'predicted_cell_type', 'region_collapsed', 'region_collapsed_filtered',
    'cellect_bmi_manual', 'cellect_bmi_manual_BH', 'manual_significant_BH',
    'leiden',
]

ct_lookup = {}
for col in ct_linked_cols:
    if col not in adata_atlas.obs.columns:
        continue
    ct_lookup[col] = (adata_atlas.obs.groupby(ct_col, observed=True)[col]
                      .first().to_dict())

# ── 4. build_obs — one function, all logic in one place ──────────────────────

def build_obs(adata, source):
    obs = pd.DataFrame(index=adata.obs_names)

    # Core QC metrics
    for col in ['nCount_RNA', 'nFeature_RNA', 'seurat_clusters', 'size_factors']:
        obs[col] = adata.obs[col].values if col in adata.obs.columns else np.nan

    # Biology / experimental design
    for col in ['sex', 'area', 'AAV', 'AAV_ratio',
                ct_col, 'score_detailed_celltype_annotation']:
        obs[col] = adata.obs[col].values if col in adata.obs.columns else np.nan

    # score_detailed_celltype_annotation: atlas = 1 (ground truth)
    if 'score_detailed_celltype_annotation' not in adata.obs.columns:
        obs['score_detailed_celltype_annotation'] = 1.0

    # Unified animal ID
    if 'hash.mcl.ID' in adata.obs.columns:
        obs['animal_id'] = adata.obs['hash.mcl.ID'].values
    elif 'sample_id' in adata.obs.columns:
        obs['animal_id'] = (adata.obs['sample_id']
                            .str.extract(r'(\d+)')[0].values)
    else:
        obs['animal_id'] = np.nan

    # Unified lane ID
    if 'reaction_id' in adata.obs.columns:
        obs['lane_id'] = adata.obs['reaction_id'].values
    elif 'orig.ident' in adata.obs.columns:
        obs['lane_id'] = adata.obs['orig.ident'].values
    else:
        obs['lane_id'] = np.nan

    # Sort metadata: atlas keeps own values, new = Sun1pos / 1.0
    obs['lane'] = adata.obs['lane'].values \
                  if 'lane' in adata.obs.columns else 'Sun1pos'
    obs['sort'] = adata.obs['sort'].values \
                  if 'sort' in adata.obs.columns else 1.0

    # Atlas-only columns (NaN for new datasets is correct)
    for col in ['score_region_collapsed']:
        obs[col] = adata.obs[col].values if col in adata.obs.columns else np.nan

    # Cell-type-linked columns transferred via lookup
    ct_vals = adata.obs[ct_col].values
    for col, lookup in ct_lookup.items():
        obs[col] = [lookup.get(ct, np.nan) for ct in ct_vals]

    obs['dataset'] = adata.obs['dataset'].values
    obs['source']  = source
    return obs

# ── 5. Build, prefix barcodes, concatenate ───────────────────────────────────
print("\nBuilding obs tables...")
obs_atlas = build_obs(adata_atlas, 'atlas')
obs_nts   = build_obs(nts_retro,   'nts')
obs_lpbn  = build_obs(lpbn_retro,  'lpbn')

obs_atlas.index = 'atlas_' + obs_atlas.index
obs_nts.index   = 'nts_'   + obs_nts.index
obs_lpbn.index  = 'lpbn_'  + obs_lpbn.index

merged_obs = pd.concat([obs_atlas, obs_nts, obs_lpbn], axis=0)

# Drop plot — 99.5% NaN, no information content
merged_obs = merged_obs.drop(columns=['plot'], errors='ignore')

# ── 6. Report ─────────────────────────────────────────────────────────────────
print(f"\nMerged obs: {len(merged_obs)} cells × {len(merged_obs.columns)} columns")
print(f"\n{'Column':50s}  {'NaN':>7s}  {'%':>6s}")
print("-" * 68)
for col in merged_obs.columns:
    n   = merged_obs[col].isna().sum()
    pct = n / len(merged_obs) * 100
    flag = '  ⚠' if pct > 50 and col not in {
        'score_region_collapsed', 'AAV', 'AAV_ratio'} else ''
    print(f"  {col:50s}  {n:7d}  {pct:5.1f}%{flag}")

# ── 7. Rebuild adata_merged with clean obs ────────────────────────────────────
adata_merged = ad.AnnData(
    X   = merged_counts,
    obs = merged_obs,
    var = adata_atlas.var[[]].copy()
)
adata_merged.obsm['X_tsne'] = tsne_merged
adata_merged.var_names_make_unique()

print(f"\nFinal merged object: {adata_merged.n_obs} cells × "
      f"{adata_merged.n_vars} genes, "
      f"{adata_merged.obs.shape[1]} obs columns")


# In[20]:


# ── 8. Save ───────────────────────────────────────────────────────────────────
adata_merged.write_h5ad("atlas_merged_NTS_LPBN_retro.h5ad")
adata_merged


# In[2]:


adata = sc.read("atlas_merged_NTS_LPBN_retro.h5ad", cache=True)
adata


# In[ ]:


# ── Add scran_normalization layer from original source objects ────────────────

def align_layer_to_atlas(adata, shared_genes, layer_name='scran_normalization'):
    """Same alignment logic as align_to_atlas but for any named layer."""
    if layer_name not in adata.layers:
        print(f"  WARNING: '{layer_name}' not in layers — falling back to counts")
        layer_name = 'counts'
    common   = adata.var_names.isin(shared_genes)
    sub      = adata[:, common].copy()
    x        = sub.layers[layer_name]
    if not sparse.issparse(x):
        x = sparse.csr_matrix(x)
    full     = sparse.lil_matrix((sub.n_obs, len(shared_genes)), dtype=np.float32)
    gene_idx = [list(shared_genes).index(g) for g in sub.var_names]
    full[:, gene_idx] = x
    return sparse.csr_matrix(full)

print("Aligning scran_normalization layers to atlas gene space...")
scran_atlas = adata_atlas.layers['scran_normalization']
if not sparse.issparse(scran_atlas):
    scran_atlas = sparse.csr_matrix(scran_atlas)

scran_nts  = align_layer_to_atlas(nts_retro,  shared_genes)
scran_lpbn = align_layer_to_atlas(lpbn_retro, shared_genes)

scran_merged = sparse.vstack([scran_atlas, scran_nts, scran_lpbn]).tocsr()
print(f"Merged scran layer: {scran_merged.shape}")

# Verify non-zero counts are sensible
n_nonzero = scran_merged.nnz
print(f"Non-zero entries: {n_nonzero} "
      f"({n_nonzero / (scran_merged.shape[0] * scran_merged.shape[1]) * 100:.2f}% dense)")


# In[11]:


# ── Reload merged object and add the layer ────────────────────────────────────
adata = sc.read("atlas_merged_NTS_LPBN_retro.h5ad")
adata.layers['scran_normalization'] = scran_merged

# Quick sanity check — v147 expression should be non-zero in sorted cells
v147_idx = list(adata.var_names).index('v147')
v147_scran = np.asarray(
    scran_merged[:, v147_idx].todense()).flatten()
for ds in ['atlas', 'nts', 'lpbn']:
    mask  = adata.obs['source'] == ds
    n_pos = (v147_scran[mask] > 0).sum()
    n_tot = mask.sum()
    print(f"  {ds:6s} v147 scran >0: {n_pos}/{n_tot} ({n_pos/n_tot*100:.1f}%)")

adata.write_h5ad("atlas_merged_NTS_LPBN_retro_scran.h5ad")
print("\nSaved with scran_normalization layer.")
print(adata)


# In[13]:


adata = sc.read("atlas_merged_NTS_LPBN_retro_scran.h5ad", cache=True)
# Exclude LPBN retrograde dataset before scoring
adata = adata[adata.obs["dataset"] != "LPBN retrograde (v36)"].copy()
adata


# In[14]:


# 1) Map areas → colors (updated with 8th area)
area_colors = {
    "BST": "#5ac2af",
    "CEA": "#b85283",
    "LHA": "#f0a672",
    "PAG": "#9a84b6",
    "PBN": "#e8682d",
    "PVH": "#185d5c",
    "PVT": "#00a0ba",
    "DVC": "#c4d95a",
}
zero_color = "#BDBDBD"


# In[17]:


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
for region in combined.index:
    colname = f"score_{region}"
    # map celltype → value for this region
    m = combined.loc[region]
    adata.obs[colname] = adata.obs["detailed_celltype_annotation"].map(m)
    score_col = "region_celltype_score"   # already in adata.obs
area_col  = "area"

# start with 'other' everywhere
adata.obs["area_highlight"] = "below threshold"

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
area_color_map["below threshold"] = "#bdbdbd"

# build the palette in the correct category order
cats = list(adata.obs["area_highlight"].cat.categories)
adata.uns["area_highlight_colors"] = [area_color_map.get(c, "#bdbdbd") for c in cats]

# Build a new column with area name if score>0, else "zero"
score_col = "region_celltype_score"  # the column you already computed
new_col   = "area_score_color_label" # name for the new label column

adata.obs[new_col] = np.where(
    adata.obs[score_col] > 0.36,
    adata.obs["area"],      # label = area (BST, CEA, …)
    "below threshold"                  # label = "zero" if score == 0
)

# Make it categorical with a fixed ordering (areas + zero)
cats = list(area_colors.keys()) + ["below threshold"]
adata.obs[new_col] = adata.obs[new_col].astype(
    pd.CategoricalDtype(categories=cats, ordered=True)
)
# Install the matching color palette — now driven entirely by area_colors,
# so adding/removing areas only requires editing that dict
palette = [area_colors[a] for a in area_colors.keys()] + [zero_color]
adata.uns[f"{new_col}_colors"] = palette

# Create a custom alpha array: low alpha for "no_projection", high for others
alpha_values = np.where(adata.obs[new_col] == "below threshold", 0.1, 1.0)
sc.pl.embedding(
    adata,
    basis="tsne",               
    color=new_col,                  
    size=3,
    alpha=alpha_values,
    title="Area-specific coloring for projection score > 0.36",
    save="_Fig1b_new.pdf"
)


# In[3]:


###prepare CELLEX inputs
outdir = "all_pickles"
os.makedirs(outdir, exist_ok=True)

# Cell types to exclude from all areas
exclude_cts = {
    'Epithelial-like_cell', 'Pituitary_Pappa2', '35_Lhx9_Zic1',
    '01_Undefined', 'Pituitary_Crhr1', 'Undefined', 'Schwann_cell',
    'Pituitary_Gnrhr'
}

ct_col = 'detailed_celltype_annotation'
ct_mask = ~adata.obs[ct_col].isin(exclude_cts)
adata = adata[ct_mask].copy()

mat = adata.layers["scran_normalization"]
counts = pd.DataFrame(
    mat.toarray() if sparse.issparse(mat) else mat,
    index   = adata.obs_names,
    columns = adata.var_names
)
meta = adata.obs[ct_col].to_frame()

counts.to_parquet(os.path.join(outdir, f"adata_counts.parquet"))
meta.to_parquet(  os.path.join(outdir, f"adata_detailed.parquet"))

ct_counts = meta[ct_col].value_counts()

print(f"Excluding cell types: {exclude_cts}\n")


# In[5]:


###prepare CELLEX inputs
# ── Pickle creation with cell type filtering and LPBN split ──────────────────
outdir = "area_pickles"
os.makedirs(outdir, exist_ok=True)

# Cell types to exclude from all areas
exclude_cts = {
    'Epithelial-like_cell', 'Pituitary_Pappa2', '35_Lhx9_Zic1',
    '01_Undefined', 'Pituitary_Crhr1', 'Undefined', 'Schwann_cell',
    'Pituitary_Gnrhr'
}

ct_col = 'detailed_celltype_annotation'

def save_pickle(sub, label):
    """Filter cell types and low-confidence labels, build DataFrames and save parquet."""

    # Remove explicitly excluded cell types
    ct_mask = ~sub.obs[ct_col].isin(exclude_cts)
    sub     = sub[ct_mask].copy()

    # Remove cells with low label transfer confidence
    # Atlas cells have score=1.0 (ground truth) so are always kept
    # New dataset cells with score < 0.4 are dropped
    if 'score_detailed_celltype_annotation' in sub.obs.columns:
        score_mask = sub.obs['score_detailed_celltype_annotation'] >= 0.4
        n_before   = sub.n_obs
        sub        = sub[score_mask].copy()
        n_dropped  = n_before - sub.n_obs
        if n_dropped > 0:
            print(f"  {label}: dropped {n_dropped} cells with score < 0.4 "
                  f"({n_dropped/n_before*100:.1f}%)")

    if sub.n_obs == 0:
        print(f"  WARNING: {label} has 0 cells after filtering — skipping")
        return

    mat = sub.layers["scran_normalization"]
    counts = pd.DataFrame(
        mat.toarray() if sparse.issparse(mat) else mat,
        index   = sub.obs_names,
        columns = sub.var_names
    )
    meta = sub.obs[ct_col].to_frame()

    counts.to_parquet(os.path.join(outdir, f"{label}_counts.parquet"))
    meta.to_parquet(  os.path.join(outdir, f"{label}_detailed.parquet"))

    ct_counts = meta[ct_col].value_counts()
    print(f"  {label:12s}: {sub.n_obs:6d} cells, "
          f"{meta[ct_col].nunique():3d} cell types, "
          f"min per CT = {ct_counts.min()}, max = {ct_counts.max()}")

print(f"Excluding cell types: {exclude_cts}\n")

# ── Standard areas (all sources combined) ────────────────────────────────────
print("Standard areas (all sources):")
for area in sorted(adata.obs["area"].unique()):
    mask = adata.obs["area"] == area
    save_pickle(adata[mask], area)

# ── LPBN split: combined PBN (already saved above as 'PBN') ──────────────────
# vs new LPBN cells only (source == 'lpbn')
print("\nLPBN new cells only (source='lpbn'):")
lpbn_mask = (adata.obs["area"] == "PBN") & (adata.obs["source"] == "lpbn")
save_pickle(adata[lpbn_mask], "LPBN_new")

# Sanity check — confirm PBN contains both sources
print("\nPBN source breakdown:")
pbn_mask = adata.obs["area"] == "PBN"
print(adata[pbn_mask].obs["source"].value_counts())


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


# In[8]:


# ── Load original 7plex results for comparison ────────────────────────────────
orig_indir = "area_cellex_results"

records_orig = []
for fn in os.listdir(orig_indir):
    if not fn.endswith("_detailed_cellex.csv"):
        continue
    area = fn.replace("_detailed_cellex.csv", "")
    df   = pd.read_csv(os.path.join(orig_indir, fn), index_col=0)
    if "v147" not in df.index:
        continue
    records_orig.append(df.loc["v147"].rename(area))

combined_orig = (pd.DataFrame(records_orig)
                   .reindex(columns=celltypes, fill_value=0)
                   .sort_index())

# ── 1. Check unchanged regions are identical ──────────────────────────────────
unchanged = ['BST', 'CEA', 'LHA', 'PAG', 'PVH', 'PVT']

print("Unchanged region comparison (new merged vs original 7plex):")
print(f"{'Area':8s}  {'Max diff':>10s}  {'Mean diff':>10s}  {'Identical':>10s}")
print("-" * 45)
for area in unchanged:
    if area not in combined_orig.index:
        print(f"{area:8s}  NOT IN ORIGINAL")
        continue
    new_row  = combined.loc[area]
    orig_row = combined_orig.loc[area]
    # Align columns
    shared   = sorted(set(new_row.index) & set(orig_row.index))
    diff     = (new_row[shared] - orig_row[shared]).abs()
    max_diff  = diff.max()
    mean_diff = diff.mean()
    identical = max_diff < 1e-6
    flag = '✓' if identical else f'⚠  max diff at {diff.idxmax()}'
    print(f"{area:8s}  {max_diff:10.6f}  {mean_diff:10.6f}  {flag}")

# ── 2. PBN vs LPBN_new comparison ────────────────────────────────────────────
# Only keep cell types present in both and with non-zero signal in at least one
pbn_row      = combined.loc['PBN']
lpbn_new_row = combined.loc['LPBN_new']

shared_cts = sorted(set(pbn_row.index) & set(lpbn_new_row.index))
df_compare = pd.DataFrame({
    'PBN (atlas + new)': pbn_row[shared_cts],
    'LPBN_new only':     lpbn_new_row[shared_cts],
}).fillna(0)

# Filter to cell types with any signal in either
has_signal = (df_compare > 0).any(axis=1)
df_compare = df_compare[has_signal].sort_values('PBN (atlas + new)', ascending=False)

from scipy.stats import spearmanr
r, p = spearmanr(df_compare['PBN (atlas + new)'], df_compare['LPBN_new only'])
print(f"\nPBN vs LPBN_new Spearman ρ={r:.3f}, p={p:.3e}, "
      f"n={len(df_compare)} cell types with signal")

print(f"\n{'Cell type':30s}  {'PBN':>8s}  {'LPBN_new':>10s}  {'Diff':>8s}")
print("-" * 65)
for ct, row in df_compare.iterrows():
    diff = row['LPBN_new only'] - row['PBN (atlas + new)']
    flag = ' ▲' if diff >  0.1 else (' ▼' if diff < -0.1 else '')
    print(f"{ct:30s}  {row['PBN (atlas + new)']:8.3f}  "
          f"{row['LPBN_new only']:10.3f}  {diff:+8.3f}{flag}")

# ── 3. Figure ─────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(16, 6))

# Panel A: scatter PBN vs LPBN_new
ax = axes[0]
ax.scatter(df_compare['PBN (atlas + new)'], df_compare['LPBN_new only'],
           s=40, alpha=0.8, edgecolors='k', linewidths=0.3, color='steelblue')
# Label top divergent points
for ct, row in df_compare.iterrows():
    diff = abs(row['LPBN_new only'] - row['PBN (atlas + new)'])
    if diff > 0.15 or max(row) > 0.3:
        ax.annotate(ct, (row['PBN (atlas + new)'], row['LPBN_new only']),
                    fontsize=6.5, xytext=(4, 0), textcoords='offset points')
lim = max(df_compare.values.max(), 0.1) * 1.1
ax.plot([0, lim], [0, lim], 'k--', lw=0.8, alpha=0.5, label='y=x')
ax.set_xlabel('ESmu — PBN (atlas + new LPBN)', fontsize=10)
ax.set_ylabel('ESmu — LPBN_new only', fontsize=10)
ax.set_title(f'PBN vs LPBN_new\nSpearman ρ={r:.3f}, p={p:.2e}',
             fontweight='bold')
ax.legend(fontsize=8)
ax.grid(alpha=0.3)

# Panel B: bar chart of differences (LPBN_new − PBN), top movers
top_movers = (df_compare['LPBN_new only'] - df_compare['PBN (atlas + new)']
              ).abs().nlargest(20).index
diff_plot  = (df_compare.loc[top_movers, 'LPBN_new only'] -
              df_compare.loc[top_movers, 'PBN (atlas + new)']
             ).sort_values()

colors = ['#C0392B' if v > 0 else '#2980B9' for v in diff_plot]
axes[1].barh(range(len(diff_plot)), diff_plot.values,
             color=colors, edgecolor='white', linewidth=0.5, height=0.7)
axes[1].set_yticks(range(len(diff_plot)))
axes[1].set_yticklabels(diff_plot.index, fontsize=8)
axes[1].axvline(0, color='k', lw=0.8)
axes[1].set_xlabel('ESmu difference (LPBN_new − PBN)')
axes[1].set_title('Top divergent cell types\n(red = higher in LPBN_new)',
                  fontweight='bold')
axes[1].grid(axis='x', alpha=0.3)

legend_el = [
    plt.matplotlib.patches.Patch(fc='#C0392B', label='Higher in LPBN_new'),
    plt.matplotlib.patches.Patch(fc='#2980B9', label='Higher in PBN (combined)'),
]
axes[1].legend(handles=legend_el, fontsize=8)

plt.suptitle('v147 ESmu: PBN combined vs LPBN new cells only',
             fontweight='bold', fontsize=11)
plt.tight_layout()
plt.savefig("Fig_PBN_vs_LPBN_new_cellex.pdf", bbox_inches='tight', dpi=300)
plt.show()


# In[9]:


# ── Three-way comparison: old PBN, new LPBN_new, combined PBN ────────────────

# Old PBN from original 7plex
orig_pbn = combined_orig.loc['PBN'] if 'PBN' in combined_orig.index else None
if orig_pbn is None:
    print("WARNING: PBN not found in original 7plex results")

# New datasets from 2plex directory (already loaded in `combined`)
comb_pbn  = combined.loc['PBN']
lpbn_new  = combined.loc['LPBN_new']

# Align all three to shared cell types with any signal
df_3way = pd.DataFrame({
    'Old PBN (7plex)':       orig_pbn,
    'Combined PBN (merged)': comb_pbn,
    'LPBN new only':         lpbn_new,
}).fillna(0)

has_signal = (df_3way > 0).any(axis=1)
df_3way    = df_3way[has_signal].sort_values('Old PBN (7plex)', ascending=False)

print(f"Cell types with any signal: {len(df_3way)}")
print(f"\n{'Cell type':30s}  {'Old PBN':>9s}  {'Combined':>9s}  {'LPBN_new':>9s}  "
      f"{'Δ old→comb':>10s}  {'Δ old→new':>10s}")
print("-" * 85)
for ct, row in df_3way.iterrows():
    d1 = row['Combined PBN (merged)'] - row['Old PBN (7plex)']
    d2 = row['LPBN new only']         - row['Old PBN (7plex)']
    flag = ''
    if abs(d2) > 0.15: flag = ' ▲' if d2 > 0 else ' ▼'
    print(f"{ct:30s}  {row['Old PBN (7plex)']:9.3f}  "
          f"{row['Combined PBN (merged)']:9.3f}  "
          f"{row['LPBN new only']:9.3f}  "
          f"{d1:+10.3f}  {d2:+10.3f}{flag}")

# Correlations
from scipy.stats import spearmanr
pairs = [
    ('Old PBN (7plex)',       'Combined PBN (merged)', 'old vs combined'),
    ('Old PBN (7plex)',       'LPBN new only',         'old vs LPBN_new'),
    ('Combined PBN (merged)', 'LPBN new only',         'combined vs LPBN_new'),
]
print("\nSpearman correlations:")
for c1, c2, label in pairs:
    r, p = spearmanr(df_3way[c1], df_3way[c2])
    print(f"  {label:30s}: ρ={r:.3f}, p={p:.3e}")

# ── Figure ────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(18, 7))

colors = {
    'Old PBN (7plex)':       '#2980B9',
    'Combined PBN (merged)': '#27AE60',
    'LPBN new only':         '#C0392B',
}

# Panel A: heatmap of all three
ax = axes[0]
arr  = df_3way.values.astype(float)
vmax = arr.max()
im   = ax.imshow(arr, aspect='auto', cmap='rocket_r',
                 vmin=0, vmax=vmax, interpolation='nearest')
ax.set_xticks(range(3))
ax.set_xticklabels(df_3way.columns, rotation=20, ha='right', fontsize=8)
ax.set_yticks(range(len(df_3way)))
ax.set_yticklabels(df_3way.index, fontsize=7)
ax.set_title('ESmu heatmap\n(all three datasets)', fontweight='bold')
plt.colorbar(im, ax=ax, fraction=0.02, pad=0.02, label='ESmu')

# Annotate non-zero values
for ri in range(arr.shape[0]):
    for ci in range(arr.shape[1]):
        v = arr[ri, ci]
        if v > 0.05:
            tc = 'white' if v > vmax * 0.6 else 'black'
            ax.text(ci, ri, f'{v:.2f}', ha='center', va='center',
                    fontsize=5, color=tc)

# Panel B: scatter old PBN vs LPBN_new
ax = axes[1]
ax.scatter(df_3way['Old PBN (7plex)'], df_3way['LPBN new only'],
           s=45, alpha=0.8, color='#C0392B',
           edgecolors='k', linewidths=0.3, label='Cell types')
lim = max(df_3way.values.max(), 0.1) * 1.15
ax.plot([0, lim], [0, lim], 'k--', lw=0.8, alpha=0.4, label='y=x')
for ct, row in df_3way.iterrows():
    diff = abs(row['LPBN new only'] - row['Old PBN (7plex)'])
    if diff > 0.15 or max(row) > 0.35:
        ax.annotate(ct, (row['Old PBN (7plex)'], row['LPBN new only']),
                    fontsize=6, xytext=(3, 0), textcoords='offset points')
r, p = spearmanr(df_3way['Old PBN (7plex)'], df_3way['LPBN new only'])
ax.set_xlabel('Old PBN (7plex)', fontsize=10)
ax.set_ylabel('LPBN new only', fontsize=10)
ax.set_title(f'Old PBN vs LPBN_new\nρ={r:.3f}, p={p:.2e}', fontweight='bold')
ax.legend(fontsize=8); ax.grid(alpha=0.3)
ax.set_xlim(-0.02, lim); ax.set_ylim(-0.02, lim)

# Panel C: bar chart — all three side by side for top cell types
ax = axes[2]
# Show top 20 by max ESmu across the three
top_cts = df_3way.max(axis=1).nlargest(20).index[::-1]
y       = np.arange(len(top_cts))
h       = 0.25

for i, (col, color) in enumerate(colors.items()):
    ax.barh(y + (i - 1) * h, df_3way.loc[top_cts, col].values,
            height=h, color=color, alpha=0.85,
            edgecolor='white', linewidth=0.3, label=col)

ax.set_yticks(y)
ax.set_yticklabels(top_cts, fontsize=7)
ax.set_xlabel('ESmu (v147)')
ax.set_title('Top 20 cell types\n(by max ESmu)', fontweight='bold')
ax.legend(fontsize=7, loc='lower right')
ax.axvline(0, color='k', lw=0.5)
ax.grid(axis='x', alpha=0.3)

plt.suptitle('v147 ESmu three-way comparison\n'
             'Old PBN (7plex atlas) | Combined PBN | LPBN new cells only',
             fontweight='bold', fontsize=11)
plt.tight_layout()
plt.savefig("Fig_PBN_threeway_cellex.pdf", bbox_inches='tight', dpi=300)
plt.show()


# In[10]:


sorted(adata.obs["area"].unique())


# In[11]:


adata.obs['area'] = adata.obs['area'].str.replace('NTS-R', 'DVC', regex=False)
print(sorted(adata.obs['area'].unique()))


# In[12]:


adata.write_h5ad("atlas_merged_NTS_LPBN_retro_scran.h5ad")
adata


# In[22]:


adata = sc.read("atlas_merged_NTS_LPBN_retro_scran.h5ad", cache=True)


# In[24]:


print(f"Before filtering: {adata.n_obs} cells")
print(adata.obs['source'].value_counts())

# Remove new LPBN cells and low-confidence label transfers
mask = (adata.obs['score_detailed_celltype_annotation'] >= 0.4)

adata = adata[mask].copy()

print(f"\nAfter filtering: {adata.n_obs} cells")
print(adata.obs['source'].value_counts())
print(f"\nScore range: {adata.obs['score_detailed_celltype_annotation'].min():.3f} – "
      f"{adata.obs['score_detailed_celltype_annotation'].max():.3f}")


# In[28]:


import sys
import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import sparse

ct_col    = 'detailed_celltype_annotation'
score_col = 'score_detailed_celltype_annotation'
gene      = 'v147'
layer     = 'scran_normalization'

exclude_cts = {
    'Epithelial-like_cell', 'Pituitary_Pappa2', '35_Lhx9_Zic1',
    '01_Undefined', 'Pituitary_Crhr1', 'Undefined', 'Schwann_cell',
    'Pituitary_Gnrhr'
}

outdir = "looa_esmu_results"
os.makedirs(outdir, exist_ok=True)

# ── Run leave-one-animal-out per area ─────────────────────────────────────────
# Results structure: {area: DataFrame(index=cell_types, columns=left_out_animals)}

all_looa = {}

for area in sorted(adata.obs['area'].unique()):
    print(f"\n{'='*70}")
    print(f"AREA: {area}")
    print(f"{'='*70}")

    area_mask = (
        (adata.obs['area'] == area) &
        (adata.obs[score_col] >= 0.4) &
        (~adata.obs[ct_col].isin(exclude_cts))
    )
    adata_area = adata[area_mask].copy()

    animals = sorted(adata_area.obs['animal_id'].unique())
    print(f"Animals in {area}: {animals} (n={len(animals)})")
    print(f"Total cells: {adata_area.n_obs}")

    if len(animals) < 2:
        print(f"  WARNING: only {len(animals)} animal — cannot do LOOA, skipping")
        continue

    area_records = {}

    for left_out in animals:
        print(f"\n  Leaving out animal {left_out}...")

        sub = adata_area[adata_area.obs['animal_id'] != left_out].copy()

        # Drop cell types with < 10 cells after removing animal
        ct_counts = sub.obs[ct_col].value_counts()
        valid_cts = ct_counts[ct_counts >= 10].index
        sub       = sub[sub.obs[ct_col].isin(valid_cts)].copy()

        if sub.n_obs == 0:
            print(f"    EMPTY after removing animal {left_out} — skipping")
            continue

        print(f"    Cells: {sub.n_obs}, Cell types: {sub.obs[ct_col].nunique()}")

        res = compute_esmu_single_gene(
            sub, gene, ct_col,
            layer   = layer,
            n_perm  = 1000,   # reduce for speed; increase to 10k for publication
            alpha   = 0.05,
            seed    = 42,
            verbose = False,
        )
        area_records[left_out] = res['esmu']

    if not area_records:
        print(f"  No LOOA results for {area}")
        continue

    # DataFrame: rows=cell_types, cols=left_out_animal
    df_area = pd.DataFrame(area_records).fillna(0)
    all_looa[area] = df_area
    df_area.to_csv(os.path.join(outdir, f"{area}_looa_esmu.csv"))
    print(f"\n  Saved {area}_looa_esmu.csv  "
          f"({df_area.shape[0]} cell types × {df_area.shape[1]} leave-out iterations)")

    # Quick summary
    print(f"\n  Summary (mean ± SD across LOOA iterations):")
    summary = pd.DataFrame({
        'mean':  df_area.mean(axis=1),
        'sd':    df_area.std(axis=1),
        'min':   df_area.min(axis=1),
        'max':   df_area.max(axis=1),
        'n_sig': (df_area > 0).sum(axis=1),
    }).sort_values('mean', ascending=False)
    sig = summary[summary['mean'] > 0]
    print(sig.round(4).to_string())

print(f"\n\nLOOA complete. Results saved to '{outdir}/'")

# ── Figure: per-area LOOA distribution ───────────────────────────────────────
n_areas  = len(all_looa)
fig, axes = plt.subplots(
    n_areas, 1,
    figsize=(14, 4 * n_areas),
    squeeze=False
)

for ax_row, (area, df_area) in zip(axes, all_looa.items()):
    ax = ax_row[0]

    # Only show cell types with mean ESmu > 0
    sig_cts  = df_area.mean(axis=1)
    sig_cts  = sig_cts[sig_cts > 0].sort_values(ascending=False).index
    if len(sig_cts) == 0:
        ax.set_title(f'{area} — no significant cell types', fontweight='bold')
        continue

    df_plot  = df_area.loc[sig_cts]
    means    = df_plot.mean(axis=1)
    sds      = df_plot.std(axis=1)
    mins     = df_plot.min(axis=1)
    maxs     = df_plot.max(axis=1)

    y        = np.arange(len(sig_cts))

    # Range bar (min–max)
    ax.barh(y, maxs.values - mins.values, left=mins.values,
            color='steelblue', alpha=0.25, height=0.6,
            label='Range (min–max)')

    # Mean ± SD
    ax.barh(y, sds.values * 2, left=(means - sds).values,
            color='steelblue', alpha=0.5, height=0.6,
            label='Mean ± SD')

    # Mean as vertical line
    ax.scatter(means.values, y, color='steelblue', s=40,
               zorder=4, edgecolors='k', linewidths=0.5, label='Mean')

    # Individual LOOA values as dots
    for yi, ct in enumerate(sig_cts):
        vals = df_plot.loc[ct].values
        ax.scatter(vals, np.full(len(vals), yi),
                   color='k', s=10, alpha=0.5, zorder=5)

    ax.set_yticks(y)
    ax.set_yticklabels(sig_cts, fontsize=8)
    ax.set_xlabel('ESmu (v147)')
    ax.set_title(f'{area} — LOOA ESmu  '
                 f'(n={df_area.shape[1]} leave-out iterations)',
                 fontweight='bold')
    ax.axvline(0, color='k', lw=0.5)
    ax.grid(axis='x', alpha=0.3)

    if ax_row is axes[0]:
        ax.legend(fontsize=7, loc='lower right')

plt.suptitle('Leave-one-animal-out ESmu stability for v147\n'
             'Dots = individual LOOA iterations, bar = mean ± SD, '
             'shading = range',
             fontweight='bold', fontsize=11, y=1.01)

plt.tight_layout()
plt.subplots_adjust(top=0.97)  # leaves room between suptitle and first panel
plt.savefig("Fig_LOOA_esmu_v147.pdf", bbox_inches='tight', dpi=300)
plt.show()

# ── Summary table across all areas ───────────────────────────────────────────
print("\nCross-area summary (mean ESmu, only cell types with any signal):")
summary_rows = []
for area, df_area in all_looa.items():
    for ct in df_area.index:
        vals = df_area.loc[ct]
        if vals.mean() > 0:
            summary_rows.append({
                'area':     area,
                'cell_type': ct,
                'mean_esmu': vals.mean(),
                'sd_esmu':   vals.std(),
                'min_esmu':  vals.min(),
                'max_esmu':  vals.max(),
                'cv':        vals.std() / vals.mean() if vals.mean() > 0 else np.nan,
                'n_iter':    len(vals),
                'n_sig_iter': (vals > 0).sum(),
            })

df_summary = (pd.DataFrame(summary_rows)
                .sort_values(['area', 'mean_esmu'], ascending=[True, False]))
df_summary.to_csv(os.path.join(outdir, "all_areas_looa_summary.csv"), index=False)
print(df_summary.to_string(index=False))


# In[34]:


from scipy.stats import spearmanr, linregress

# ── Fix 1: remap NTS/NTS-R → DVC in combined to match LOOA area names ────────

rename_map = {'NTS-R': 'DVC'}
combined.index = [rename_map.get(a, a) for a in combined.index]
print("\nAfter remap:", sorted(combined.index))

# Rebuild aligned DataFrames
shared_areas   = sorted(set(combined.index) & set(df_looa_mean.index))
shared_cts     = sorted(set(combined.columns) & set(df_looa_mean.columns))
cellex_aligned = combined.reindex(index=shared_areas, columns=shared_cts, fill_value=0)
looa_aligned   = df_looa_mean.reindex(index=shared_areas, columns=shared_cts, fill_value=0)
print(f"Shared areas: {shared_areas}")
print(f"Shared cell types: {len(shared_cts)}")

# ── Per-area Spearman on signal-only cell types ───────────────────────────────
def spearman_signal(c_row, l_row, min_n=3):
    """Spearman on cell types with signal in at least one method."""
    any_sig = (c_row > 0) | (l_row > 0)
    if any_sig.sum() < min_n:
        return np.nan, np.nan, any_sig.sum()
    r, p = spearmanr(c_row[any_sig], l_row[any_sig])
    return r, p, any_sig.sum()

# ── Per-area linear regression ────────────────────────────────────────────────
print("\nPer-area linear regression (LOOA → CELLEX scale):")
print(f"{'Area':8s}  {'slope':>7s}  {'intercept':>10s}  {'R²':>6s}  "
      f"{'ρ (signal CTs)':>15s}  {'N sig CTs':>10s}")
print("-" * 65)

scale_factors = {}
for area in shared_areas:
    c_row = cellex_aligned.loc[area]
    l_row = looa_aligned.loc[area]
    both  = (c_row > 0) & (l_row > 0)
    r_sp, p_sp, n_sig = spearman_signal(c_row, l_row)

    if both.sum() < 3:
        print(f"{area:8s}  insufficient overlap (n={both.sum()})")
        continue

    slope, intercept, r, p, se = linregress(l_row[both], c_row[both])
    scale_factors[area] = {'slope': slope, 'intercept': intercept, 'r2': r**2}
    print(f"{area:8s}  {slope:7.2f}  {intercept:10.4f}  {r**2:6.3f}  "
          f"{r_sp:15.3f}  {n_sig:10d}")

# ── Scatter plots ─────────────────────────────────────────────────────────────
n_cols  = 3
n_rows  = int(np.ceil(len(shared_areas) / n_cols))
fig, axes = plt.subplots(n_rows, n_cols,
                          figsize=(n_cols * 5, n_rows * 4.5),
                          squeeze=False)
axes_flat = axes.flatten()

for i, area in enumerate(shared_areas):
    ax    = axes_flat[i]
    c_row = cellex_aligned.loc[area]
    l_row = looa_aligned.loc[area]
    r_sp, p_sp, n_sig = spearman_signal(c_row, l_row)

    both_sig  = (c_row > 0) & (l_row > 0)
    only_cel  = (c_row > 0) & (l_row == 0)
    only_looa = (c_row == 0) & (l_row > 0)
    neither   = (c_row == 0) & (l_row == 0)

    for mask, color, label in [
        (neither,   '#DDDDDD', 'Neither'),
        (both_sig,  '#2980B9', 'Both significant'),
        (only_cel,  '#C0392B', 'CELLEX only'),
        (only_looa, '#27AE60', 'LOOA only'),
    ]:
        if mask.sum() == 0:
            continue
        ax.scatter(c_row[mask], l_row[mask],
                   color=color, s=35, alpha=0.8,
                   edgecolors='k', linewidths=0.3,
                   label=f'{label} (n={mask.sum()})', zorder=3)

    # Best-fit line across data range only — no y=x line
    if both_sig.sum() >= 3 and area in scale_factors:
        sf    = scale_factors[area]
        x_min = c_row[both_sig].min()
        x_max = c_row[both_sig].max() * 1.05
        x_fit = np.linspace(x_min, x_max, 100)
        y_fit = (x_fit - sf['intercept']) / sf['slope']
        ax.plot(x_fit, y_fit, 'k-', lw=1.5, alpha=0.8, zorder=4,
                label=f"Best fit  slope={sf['slope']:.1f}  R²={sf['r2']:.2f}")

    # Label divergent significant points
    for ct in shared_cts:
        if both_sig[ct] and (abs(c_row[ct] - l_row[ct]) > 0.05
                              or max(c_row[ct], l_row[ct]) > 0.08):
            ax.annotate(ct, (c_row[ct], l_row[ct]),
                        fontsize=5.5, xytext=(3, 0),
                        textcoords='offset points', va='center')

    # Scale axes tightly to data range
    x_data = c_row[c_row > 0]
    y_data = l_row[l_row > 0]
    x_pad  = x_data.max() * 0.08 if len(x_data) > 0 else 0.05
    y_pad  = y_data.max() * 0.08 if len(y_data) > 0 else 0.05
    ax.set_xlim(-x_pad * 0.3, x_data.max() + x_pad if len(x_data) > 0 else 0.1)
    ax.set_ylim(-y_pad * 0.3, y_data.max() + y_pad if len(y_data) > 0 else 0.1)

    ax.set_xlabel('Full CELLEX ESmu', fontsize=9)
    ax.set_ylabel('Raw LOOA ESmu', fontsize=9)
    r_str = f'{r_sp:.3f}' if not np.isnan(r_sp) else 'n/a'
    ax.set_title(f'{area}  ρ={r_str} (n={n_sig} CTs)',
                 fontweight='bold', fontsize=10)
    ax.legend(fontsize=6, loc='upper left')
    ax.grid(alpha=0.3)

for j in range(len(shared_areas), len(axes_flat)):
    axes_flat[j].set_visible(False)

plt.suptitle('Full CELLEX ESmu vs raw LOOA ESmu — per-area best fit\n'
             'ρ computed on cell types with signal in at least one method',
             fontweight='bold', fontsize=12)
plt.tight_layout()
plt.savefig("Fig_cellex_vs_looa_bestfit.pdf", bbox_inches='tight', dpi=300)
plt.show()

# ── Rebuild df_robust with rescaled LOOA error bars ──────────────────────────
print("\nRebuilding confidence table with rescaled LOOA values...")
records = []
for area in shared_areas:
    sf = scale_factors.get(area, None)
    for ct in shared_cts:
        cellex_val = cellex_aligned.loc[area, ct]

        if area in all_looa and ct in all_looa[area].index:
            looa_vals = all_looa[area].loc[ct]
            if sf is not None:
                looa_vals_scaled = (looa_vals * sf['slope'] + sf['intercept']).clip(lower=0)
            else:
                looa_vals_scaled = looa_vals
            looa_sd_scaled = looa_vals_scaled.std()
            looa_cv        = (looa_sd_scaled / looa_vals_scaled.mean()
                              if looa_vals_scaled.mean() > 0 else np.nan)
            n_iter         = len(looa_vals)
            n_sig_iter     = (looa_vals > 0).sum()
        else:
            looa_sd_scaled = looa_cv = np.nan
            n_iter = n_sig_iter = 0

        records.append({
            'area':           area,
            'cell_type':      ct,
            'esmu_cellex':    cellex_val,
            'looa_sd_scaled': looa_sd_scaled,
            'looa_cv':        looa_cv,
            'n_iter':         n_iter,
            'n_sig_iter':     n_sig_iter,
            'pct_sig':        n_sig_iter / n_iter * 100 if n_iter > 0 else 0,
            'robust':         (cellex_val > 0) and
                              (n_sig_iter / n_iter >= 0.5 if n_iter > 0 else False),
        })

df_robust = pd.DataFrame(records)
df_robust.to_csv("looa_esmu_results/cellex_with_looa_confidence_rescaled.csv",
                 index=False)
print("Saved cellex_with_looa_confidence_rescaled.csv")

# ── Robustness summary ────────────────────────────────────────────────────────
print(f"\n{'Area':8s}  {'N sig CELLEX':>13s}  {'N robust':>9s}  {'% robust':>9s}")
print("-" * 45)
for area in shared_areas:
    sub      = df_robust[df_robust['area'] == area]
    n_sig    = (sub['esmu_cellex'] > 0).sum()
    n_robust = sub['robust'].sum()
    pct      = n_robust / n_sig * 100 if n_sig > 0 else 0
    print(f"{area:8s}  {n_sig:13d}  {n_robust:9d}  {pct:8.1f}%")

print("\nRobust cell types (sig in full CELLEX + ≥50% LOOA iters):")
for area in shared_areas:
    hits = (df_robust[(df_robust['area'] == area) & df_robust['robust']]
            .sort_values('esmu_cellex', ascending=False))
    if len(hits) == 0:
        continue
    print(f"\n{area}:")
    for _, r in hits.iterrows():
        print(f"  {r['cell_type']:30s}  ESmu={r['esmu_cellex']:.3f}  "
              f"LOOA={r['n_sig_iter']:.0f}/{r['n_iter']:.0f}  "
              f"CV={r['looa_cv']:.2f}")


# In[32]:


# ── Fix 1: remap NTS/NTS-R → DVC in combined to match LOOA area names ────────
print("Areas in CELLEX combined:", sorted(combined.index))
print("Areas in LOOA:           ", sorted(all_looa.keys()))


# In[30]:


# ── Use CELLEX as point estimate + LOOA for confidence intervals ──────────────

# Build a unified results table:
# point estimate = full CELLEX ESmu
# uncertainty    = SD and CV from LOOA iterations
# robustness     = n_sig_iter / n_iter from LOOA

records = []
for area in shared_areas:
    for ct in shared_cts:
        cellex_val = cellex_aligned.loc[area, ct]

        if area in all_looa and ct in all_looa[area].index:
            looa_vals   = all_looa[area].loc[ct]
            looa_sd     = looa_vals.std()
            looa_cv     = looa_sd / looa_vals.mean() if looa_vals.mean() > 0 else np.nan
            n_iter      = len(looa_vals)
            n_sig_iter  = (looa_vals > 0).sum()
            looa_mean   = looa_vals.mean()
        else:
            looa_sd = looa_cv = looa_mean = np.nan
            n_iter  = n_sig_iter = 0

        records.append({
            'area':          area,
            'cell_type':     ct,
            'esmu_cellex':   cellex_val,        # point estimate
            'looa_sd':       looa_sd,           # uncertainty
            'looa_cv':       looa_cv,           # relative stability
            'looa_mean':     looa_mean,         # for reference
            'n_iter':        n_iter,
            'n_sig_iter':    n_sig_iter,        # robustness vote
            'pct_sig':       n_sig_iter / n_iter * 100 if n_iter > 0 else 0,
            # Robust = significant in full CELLEX AND majority of LOOA iters
            'robust':        (cellex_val > 0) and (n_sig_iter / n_iter >= 0.5
                              if n_iter > 0 else False),
        })

df_robust = pd.DataFrame(records)
df_robust.to_csv("looa_esmu_results/cellex_with_looa_confidence.csv", index=False)

# ── Figure: CELLEX ESmu with LOOA-derived error bars ─────────────────────────
fig, axes = plt.subplots(
    len(shared_areas), 1,
    figsize=(14, 4.5 * len(shared_areas)),
    squeeze=False
)

for ax_row, area in zip(axes, shared_areas):
    ax = ax_row[0]

    sub = (df_robust[df_robust['area'] == area]
           .query('esmu_cellex > 0')
           .sort_values('esmu_cellex', ascending=True))

    if len(sub) == 0:
        ax.set_title(f'{area} — no signal'); continue

    y      = np.arange(len(sub))
    colors = ['#2980B9' if r else '#E8A87C'
              for r in sub['robust']]

    # CELLEX point estimate as bar
    ax.barh(y, sub['esmu_cellex'], color=colors,
            alpha=0.7, height=0.6, edgecolor='white', linewidth=0.4)

    # LOOA SD as error bar — scaled to CELLEX magnitude
    # We scale LOOA SD by the ratio cellex/looa_mean to put it on CELLEX scale
    for i, (_, row) in enumerate(sub.iterrows()):
        if row['looa_mean'] > 0 and not np.isnan(row['looa_sd']):
            # Scale factor to bring LOOA uncertainty to CELLEX scale
            scale   = row['esmu_cellex'] / row['looa_mean']
            err_scaled = row['looa_sd'] * scale
            ax.errorbar(row['esmu_cellex'], i,
                        xerr=err_scaled,
                        fmt='none', color='k',
                        capsize=3, capthick=0.8, elinewidth=0.8, zorder=5)

        # Annotate robustness
        n_s = int(row['n_sig_iter'])
        n_t = int(row['n_iter'])
        ax.text(row['esmu_cellex'] + 0.01, i,
                f'{n_s}/{n_t}',
                va='center', fontsize=6, color='#444')

    ax.set_yticks(y)
    ax.set_yticklabels(sub['cell_type'], fontsize=7)
    ax.set_xlabel('ESmu (full CELLEX)')
    ax.set_title(f'{area}  —  CELLEX ESmu ± scaled LOOA SD\n'
                 f'(n/n = LOOA iterations with ESmu > 0, '
                 f'blue = robust in ≥50% iterations)',
                 fontweight='bold')
    ax.grid(axis='x', alpha=0.3)
    ax.axvline(0, color='k', lw=0.5)

    legend_el = [
        mpatches.Patch(fc='#2980B9', label='Robust (sig in ≥50% LOOA iters)'),
        mpatches.Patch(fc='#E8A87C', label='Unstable (<50% LOOA iters)'),
    ]
    ax.legend(handles=legend_el, fontsize=7, loc='lower right')

plt.suptitle('v147 projection specificity: full CELLEX ESmu with LOOA stability\n'
             'Error bars = LOOA SD scaled to CELLEX magnitude',
             fontweight='bold', fontsize=12, y=1.01)
plt.tight_layout()
plt.subplots_adjust(top=0.97)
plt.savefig("Fig_cellex_esmu_with_looa_errorbars.pdf", bbox_inches='tight', dpi=300)
plt.show()

# ── Summary: how many cell types are robust per area? ────────────────────────
print("\nRobustness summary per area:")
print(f"{'Area':8s}  {'N sig CELLEX':>13s}  {'N robust':>9s}  {'% robust':>9s}")
print("-" * 45)
for area in shared_areas:
    sub = df_robust[df_robust['area'] == area]
    n_sig    = (sub['esmu_cellex'] > 0).sum()
    n_robust = sub['robust'].sum()
    pct      = n_robust / n_sig * 100 if n_sig > 0 else 0
    print(f"{area:8s}  {n_sig:13d}  {n_robust:9d}  {pct:8.1f}%")

# Print robust hits only
print("\nRobust cell types per area (sig in full CELLEX + ≥50% LOOA iters):")
for area in shared_areas:
    hits = (df_robust[(df_robust['area'] == area) & df_robust['robust']]
            .sort_values('esmu_cellex', ascending=False))
    if len(hits) == 0:
        continue
    print(f"\n{area}:")
    for _, r in hits.iterrows():
        print(f"  {r['cell_type']:30s}  ESmu={r['esmu_cellex']:.3f}  "
              f"LOOA={r['n_sig_iter']:.0f}/{r['n_iter']:.0f}  "
              f"CV={r['looa_cv']:.2f}")


# In[ ]:




