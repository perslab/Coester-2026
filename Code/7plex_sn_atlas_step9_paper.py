#!/usr/bin/env python
# coding: utf-8

# In[2]:


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
os.chdir("/projects/perslab/people/qbx224/7plex_paper/atlas/")
from pathlib import Path

outdir = Path("/projects/perslab/people/qbx224/7plex_paper/atlas/paper/")
sc.settings.figdir = outdir  # <- controls where scanpy writes


# In[3]:


#adata = sc.read("allgenes_bipca_hypomap.h5ad")
adata = sc.read("allgenes_bipca_annotated.h5ad", cache=True)
adata


# In[3]:


ax = sc.pl.embedding(
    adata,
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
plt.savefig("paper/Fig_1a.pdf")


# In[10]:


import re
import matplotlib as mpl
def transform_label(label):
    if re.match(r'^\d', label):              # neuron: starts with a digit
        return label.split('_')[0]           # "06_Tbx3_Npy1r" → "06"
    elif 'Pituitary' in label:               # merge all Pituitary subtypes
        return 'Pituitary'
    else:                                    # glial cell
        return label.replace('_', ' ')      # "Endothelial_cell" → "Endothelial cell"

adata.obs['label_short'] = (
    adata.obs['detailed_celltype_annotation']
    .astype(str)
    .map(transform_label)
)

ax = sc.pl.embedding(
    adata,
    basis="tsne",
    color="label_short",
    legend_loc="on data",
    legend_fontsize="small",
    return_fig=True,
    title="Annotated hypothalamic nuclei \n (Neurons are numbered)")

ax = plt.gca()

# Add white outline to all text labels for contrast
for text in ax.texts:
    text.set_path_effects([
        mpl.patheffects.withStroke(linewidth=2, foreground='white')
    ])

adjustText.adjust_text(
    ax.texts,
    ax=ax,
    arrowprops=dict(arrowstyle="-", color='black', lw=1.5))

# Re-apply path effects after adjustText (it may reset them)
for text in ax.texts:
    text.set_path_effects([
        mpl.patheffects.withStroke(linewidth=2, foreground='white')
    ])
plt.savefig("paper/Fig_1a_new.pdf")


# In[15]:


celltypes = ['04_Crnde_Lmx1a', '27_Lmx1a_Ttc6', '32_Gucy2c_Ttc6', "36_Ebf2_Crhbp"]

ax = sc.pl.embedding(
    adata,
    basis="tsne",
    color="Lmx1a",
    layer="scran_normalization", 
    color_map="YlOrBr",
    title="Marker gene expression for Lmx1a-Cre driver gene selection",
    return_fig=True,
    vmax=1.5)

ax = plt.gca()

# Compute centroid of each celltype in tSNE space and annotate
tsne_coords = adata.obsm["tsne"]
for ct in celltypes:
    mask = adata.obs["detailed_celltype_annotation"] == ct
    if mask.sum() == 0:
        continue
    cx, cy = tsne_coords[mask].mean(axis=0)
    short = ct.split("_")[0]   # just the number, e.g. "04"
    ax.text(
        cx, cy, short,
        fontsize=7,
        fontweight="bold",
        ha="center", va="center",
        color="black",
        path_effects=[mpl.patheffects.withStroke(linewidth=2, foreground="white")]
    )

plt.savefig("paper/_Lmx1a.pdf")


# In[17]:


# --- 1. Pull expression + cell type into a DataFrame ---
gene = "Lmx1a"
expr = pd.Series(
    np.array(adata[:, gene].layers["counts"].todense()).flatten(),
    index=adata.obs_names,
    name=gene
)
df = adata.obs[["detailed_celltype_annotation"]].copy()
df[gene] = expr.values

# --- 2. Rank cell types by mean expression ---
ct_order = (
    df.groupby("detailed_celltype_annotation")[gene]
    .mean()
    .sort_values(ascending=False)
    .index.tolist()
)

# --- 3. Pull the existing color palette from adata.uns ---
cats   = list(adata.obs["detailed_celltype_annotation"].cat.categories)
colors = adata.uns["detailed_celltype_annotation_colors"]
color_map = dict(zip(cats, colors))

# --- 4. Plot ---
fig, ax = plt.subplots(figsize=(len(ct_order) * 0.4, 5))

sns.violinplot(
    data=df,
    x="detailed_celltype_annotation",
    y=gene,
    order=ct_order,
    palette=color_map,
    inner="box",       # small box-and-whisker inside each violin
    linewidth=0.5,
    cut=0,             # don't extend beyond data range
    ax=ax
)

# --- 5. Tidy up axes ---
ax.set_xticklabels(
    [ct.split("_")[0] for ct in ct_order],   # short numeric labels
    rotation=90, fontsize=6
)
ax.set_xlabel("Cell type (ranked by mean expression)")
ax.set_ylabel(f"{gene} (counts)")
ax.set_title(f"{gene} expression across cell types")
sns.despine()

plt.tight_layout()
plt.savefig("paper/Fig_Lmx1a_violin.pdf")


# In[13]:


adata.X = adata.layers["scran_normalization"]
sc.tl.rank_genes_groups(adata, groupby="detailed_celltype_annotation", method="wilcoxon", key_added="dea_leiden")


# In[14]:


sc.tl.dendrogram(adata, groupby="detailed_celltype_annotation", use_rep = "pcs", n_pcs=adata.uns["bipca"]['rank'])


# In[15]:


sc.tl.filter_rank_genes_groups(
    adata,
    min_in_group_fraction=0.2,
    max_out_group_fraction=0.2,
    key="dea_leiden",
    key_added="dea_leiden_filtered"
)
#To preserve the original structure of adata.uns[‘rank_genes_groups’], filtered genes are set to NaN.
#This needs to be removed because it raises an error while saving to h5ad.


# In[16]:


sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby="detailed_celltype_annotation",
    standard_scale="var",
    n_genes=4,
    key="dea_leiden_filtered",
    save="dotplot_filtered.pdf"
)


# In[5]:


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
adata.obs["area_highlight"] = "no_projection"

# cells that should be colored by their area
mask = score_col > 0.1

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
    adata.obs[score_col] > 0.36,
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
# Create a custom alpha array: low alpha for "no_projection", high for others
alpha_values = np.where(adata.obs[new_col] == "no_projection", 0.1, 1.0)
sc.pl.embedding(
    adata,
    basis="tsne",               
    color=new_col,                  
    size=3,
    alpha=alpha_values,
    title="Area-specific coloring for score > 0.36",
    save="_Fig1b.svg"
)


# In[2]:


#del adata
adata = sc.read("filtered_bipca_full_annotation.h5ad", cache=True)
adata


# In[9]:


sc.pl.correlation_matrix(
    adata, 
    "detailed_celltype_annotation",
    save="_corr_mat_hvg_glia.pdf"
)


# In[3]:


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


# In[11]:


import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch
from matplotlib.colors import to_rgb

# Get the dendrogram ordering
dendro_info = adata.uns['dendrogram_detailed_celltype_annotation']
cell_type_order = dendro_info['categories_ordered']

# Get neurotransmitter for each cell type (taking the most common per cell type)
neurotransmitter_map = {}
for ct in cell_type_order:
    mask = adata.obs['detailed_celltype_annotation'] == ct
    nt_values = adata.obs.loc[mask, 'neurotransmitter'].value_counts()
    neurotransmitter_map[ct] = nt_values.index[0] if len(nt_values) > 0 else 'Unknown'

# Your specific color scheme
cols = ('#FF4A46', '#1CE6FF', '#FFFF00', '#FF34FF', "#008941", "#006FA6")
nt = ("Glut", "GABA", "Chat", "Hist", "Pituitary", "Glia")
nt_color_map = dict(zip(nt, cols))

# Convert hex colors to RGB for imshow
nt_color_map_rgb = {n: to_rgb(c) for n, c in nt_color_map.items()}
nt_color_map_rgb['Unknown'] = to_rgb('#CCCCCC')

# Create annotation array with RGB values
nt_array = np.array([[nt_color_map_rgb.get(neurotransmitter_map[ct], nt_color_map_rgb['Unknown']) 
                      for ct in cell_type_order]])

# Calculate figure height to make squares (1 row, so height = width / n_columns)
n_cols = len(cell_type_order)
fig_width = 18
fig_height = fig_width / n_cols  # This makes each cell square

# Create figure
fig, ax = plt.subplots(figsize=(fig_width, fig_height))

# Plot annotation with equal aspect to maintain squares
ax.imshow(nt_array, aspect='equal', interpolation='nearest')

# Add faint grid
ax.set_xticks(np.arange(-.5, len(cell_type_order), 1), minor=True)
ax.set_yticks(np.arange(-.5, 1, 1), minor=True)
ax.grid(which='minor', color='gray', linestyle='-', linewidth=0.5, alpha=0.3)

# Set labels
ax.set_yticks([0])
ax.set_yticklabels(['Neurotransmitter'], fontsize=12)
ax.set_xticks(range(len(cell_type_order)))
ax.set_xticklabels(cell_type_order, rotation=90, ha='right', fontsize=8)

# Add legend
legend_elements = [Patch(facecolor=nt_color_map[n], label=n) for n in nt]
ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.02, 1), 
          fontsize=10, title='Neurotransmitter')

plt.tight_layout()
plt.savefig('neurotransmitter_annotation_glia.pdf', bbox_inches='tight', dpi=300)
plt.show()


# In[14]:


###This is used to run area-wise CELLEX with leaving out 1 animal at a time
area_to_animals = (
    adata.obs.groupby("area")["hash.mcl.ID"]
    .unique()
    .apply(list)
)
meta = adata.obs[["hash.mcl.ID", "detailed_celltype_annotation", "area"]].copy()
meta.to_parquet("cellex_metadata_area_glia.parquet")
area_to_animals


# In[16]:


#-------------------------------------
#Singlecell CELLEX 1 animal left out
#-------------------------------------
indir = "LOO_cellex_results_glia"
records = []

# Loop over each area folder
for area in os.listdir(indir):
    area_dir = os.path.join(indir, area)
    if not os.path.isdir(area_dir):
        continue

    print(f"Scanning area directory: {area}")

    # Loop over all CELLEX result files inside this area folder
    for fn in os.listdir(area_dir):
        if not fn.startswith("CELLEX_ESmu_"):
            continue
        if not fn.endswith(".csv"):
            continue

        path = os.path.join(area_dir, fn)

        # Example filename:
        # "CELLEX_ESmu_LHA_LOO_animal_124.csv"
        # Extract left-out animal ID
        try:
            left_out_animal = fn.split("_animal_")[1].replace(".csv", "")
        except:
            print(f"Could not extract animal ID from {fn}, skipping")
            continue

        # Load CELLEX ESmu table
        df = pd.read_csv(path, index_col=0)

        if "v147" not in df.index:
            print(f"Skipping {fn}: v147 not found")
            continue

        s = df.loc["v147"]  # ESmu values per cluster

        # Build long format records
        for cluster, esmu in s.items():
            records.append({
                "area": area,
                "left_out_animal": left_out_animal,
                "cluster": cluster,
                "esmu": esmu
            })

# Combine into one DataFrame
combined_LOO = pd.DataFrame(records)
combined_LOO.to_csv('df_combined_LOO_new.csv', index=False)
#-------------------------------------
#Singlecell CELLEX
#-------------------------------------
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
combined.to_csv('df_combined_new.csv', index=False)
#-------------------------------------
#Spatial CELLEX
#-------------------------------------
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
combined_spatial.to_csv('df_combined_spatial_new.csv', index=False)


# In[5]:


df = combined_LOO.copy()

# group by area+cluster to compute stability across LOO folds
stab = (
    df.groupby(["area", "cluster"])
      .agg(
        mean_esmu = ("esmu", "mean"),
        sd_esmu   = ("esmu", "std"),
        min_esmu  = ("esmu", "min"),
        max_esmu  = ("esmu", "max"),
        n_folds   = ("esmu", "size"),
        zero_frac = ( "esmu", lambda x: np.mean(x == 0) )
      )
)

stab["range_esmu"] = stab["max_esmu"] - stab["min_esmu"]
stab["cv_esmu"] = stab["sd_esmu"] / stab["mean_esmu"].abs().replace(0, np.nan)

stab = stab.reset_index()
area_summary = (
    stab.groupby("area")
        .agg(
            mean_cv=("cv_esmu", "mean"),
            mean_sd=("sd_esmu", "mean"),
            mean_range=("range_esmu", "mean"),
            frac_unstable=("cv_esmu", lambda x: np.mean(x > 1.0))
        )
)
area_summary


# In[8]:


import statsmodels.formula.api as smf
import statsmodels.stats.multitest as smm
import statsmodels.api as sm

df = combined_LOO.copy()
df['left_out_animal'] = df['left_out_animal'].astype(str)
df['cluster'] = df['cluster'].astype(str)
df['area'] = df['area'].astype(str)
anova_results = []

for area in sorted(df['area'].unique()):
    d = df[df['area'] == area].copy()

    # linear model (fixed effects only)
    lm = smf.ols("esmu ~ C(cluster) + C(left_out_animal)", data=d).fit()

    aov = sm.stats.anova_lm(lm, typ=2)

    ss_cluster = aov.loc["C(cluster)", "sum_sq"]
    ss_animal  = aov.loc["C(left_out_animal)", "sum_sq"]
    ss_resid   = aov.loc["Residual", "sum_sq"]

    tot = ss_cluster + ss_animal + ss_resid

    anova_results.append({
        "area": area,
        "ss_cluster": ss_cluster,
        "ss_animal": ss_animal,
        "ss_resid": ss_resid,
        "prop_cluster": ss_cluster / tot,
        "prop_animal": ss_animal / tot,
        "prop_resid": ss_resid / tot,
    })

anova_df = pd.DataFrame(anova_results)
display(anova_df)


# In[9]:


stab_list = []

for area in sorted(df.area.unique()):
    d = df[df.area==area]

    stats = (
        d.groupby("cluster")["esmu"]
        .agg(sd="std", iqr=lambda x: x.quantile(0.75)-x.quantile(0.25),
             mad=lambda x: (x - x.median()).abs().median(),
             min="min", max="max")
    )

    stats["range"] = stats["max"] - stats["min"]

    stab_list.append({
        "area": area,
        "mean_sd": stats["sd"].mean(),
        "max_sd": stats["sd"].max(),
        "mean_range": stats["range"].mean(),
        "max_range": stats["range"].max(),
        "clusters_high_instability": (stats["sd"] > 0.1).sum()
    })

stab_df = pd.DataFrame(stab_list)
display(stab_df)


# In[10]:


# --- 1. compute instability per cluster per area ---
stab = (
    df.groupby(["area", "cluster"])["esmu"]
    .agg(
        mad=lambda x: (x - x.median()).abs().median(),
        sd="std",
        iqr=lambda x: x.quantile(0.75) - x.quantile(0.25),
        range=lambda x: x.max() - x.min()
    )
    .reset_index()
)

# use MAD as final metric
stab_matrix = stab.pivot(index="cluster", columns="area", values="mad")

# --- 2. sort clusters by overall instability (descending) ---
stab_matrix["instability_mean"] = stab_matrix.mean(axis=1)
stab_matrix = stab_matrix.sort_values("instability_mean", ascending=False)
stab_matrix = stab_matrix.drop(columns=["instability_mean"])

# --- 3. plot heatmap ---
plt.figure(figsize=(12, 14))
sns.heatmap(
    stab_matrix,
    cmap="viridis",
    linewidths=0.2,
    linecolor="gray",
    cbar_kws={"label": "ESμ instability (MAD across animals)"}
)
plt.title("Cluster ESμ Instability (MAD across Leave-One-Animal-Out runs)")
plt.xlabel("Area")
plt.ylabel("Cluster")
plt.tight_layout()
plt.savefig("Cluster ESμ Instability.png", dpi=300)
plt.show()


# In[11]:


# compute mean MAD across all areas for ranking
cluster_instab = (
    stab.groupby("cluster")["mad"]
    .mean()
    .sort_values(ascending=False)
)

top4_clusters = cluster_instab.head(10).index.tolist()
print("Top 10 unstable clusters:", top4_clusters)


# In[15]:


area_to_animals = {
    "BST": [60, 19, 20, 37, 72, 12, 40, 48, 47],
    "CEA": [126, 57, 81, 95, 100, 90, 67, 120, 73],
    "LHA": [124, 61, 93, 98, 122, 65],
    "PAG": [15, 10, 18, 46, 45, 38, 33, 58, 66],
    "PBN": [31, 34, 64, 51, 44, 56, 43, 76, 16],
    "PVH": [63, 75, 138, 97, 85, 88, 118, 69],
    "PVT": [53, 62, 41, 42, 25, 74, 54, 49],
}
region_order = ["BST", "CEA", "LHA", "PAG", "PBN", "PVH", "PVT"] 
colors = ["#4D67C8","#A7EA51","#5DCEAF","#702FA0","#FFC001","#FF8021","#5DCBF3"] 
region_to_color = dict(zip(region_order, colors))
# top4_clusters should be your list of clusters to inspect
for cl in top4_clusters:
    df_cl = df[df["cluster"] == cl].copy()
    df_cl["left_out_animal"] = df_cl["left_out_animal"].astype(int)

    # Set up figure: one subplot per area
    n_areas = len(df_cl["area"].unique())
    ncols = 10
    nrows = int(np.ceil(n_areas / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4*ncols, 3*nrows), sharey=True)
    axes = axes.flatten()

    for i, area in enumerate(sorted(df_cl["area"].unique(), key=lambda x: region_order.index(x))):
        ax = axes[i]
        # Filter to only the animals that belong to this area
        animals_in_area = area_to_animals[area]
        df_area = df_cl[(df_cl["area"] == area) & 
                        (df_cl["left_out_animal"].isin(animals_in_area))].copy()
        df_area["left_out_animal"] = df_area["left_out_animal"].astype(str)

        # Get area-specific color
        color = region_to_color.get(area, "gray")

        # Plot points
        sns.stripplot(
            data=df_area,
            x="left_out_animal",
            y="esmu",
            size=6,
            color=color,
            ax=ax
        )
        # Connect points with lines
        sns.lineplot(
            data=df_area,
            x="left_out_animal",
            y="esmu",
            estimator=None,
            sort=True,
            color=color,
            ax=ax
        )

        # Add full-data ESμ reference line (black, dotted)
        try:
            ref_esmu = combined.loc[area, cl]
            ax.axhline(ref_esmu, color='black', linestyle='--', linewidth=1.5, label='Full-data ESμ')
        except KeyError:
            print(f"Cluster {cl} not found in combined ESμ for area {area}")

        ax.set_title(area)
        ax.set_xlabel("Left-out animal")
        ax.set_ylabel("ESμ value")
        ax.tick_params(axis='x', rotation=45)
        ax.legend()

    # Hide any empty subplots
    for j in range(i+1, len(axes)):
        axes[j].axis('off')

    fig.suptitle(f"Instability Diagnostic: {cl}", fontsize=14, y=1.03)
    plt.tight_layout()
    plt.savefig(f"Instability Diagnostic: {cl}.png", dpi=300)
    plt.show()


# In[17]:


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch
from matplotlib.colors import to_rgb

# Get the dendrogram ordering
dendro_info = adata.uns['dendrogram_detailed_celltype_annotation']
cell_type_order = dendro_info['categories_ordered']

# Assuming your dataframes are named 'combined' and 'combined_spatial'
# Create the flag dataframe
flag_df = pd.DataFrame(index=combined.index, columns=combined.columns)

# Determine flags based on threshold > 0.36
for region in combined.index:
    for celltype in combined.columns:
        snrna_pos = combined.loc[region, celltype] > 0.36
        spatial_pos = combined_spatial.loc[region, celltype] > 0.36 if region in combined_spatial.index else False

        if snrna_pos and spatial_pos:
            flag_df.loc[region, celltype] = 'both'
        elif snrna_pos:
            flag_df.loc[region, celltype] = 'snRNA'
        elif spatial_pos:
            flag_df.loc[region, celltype] = 'spatial'
        else:
            flag_df.loc[region, celltype] = 'negative'

# Reorder columns to match dendrogram order
flag_df = flag_df[cell_type_order]

# Define colors for each category
category_colors = {
    'negative': '#FFFFFF',  # white
    'snRNA': '#CC0744', 
    'spatial': '#A079BF',  
    'both': '#000035'     
}

# Convert to RGB
category_colors_rgb = {cat: to_rgb(col) for cat, col in category_colors.items()}

# Create annotation array with RGB values
annotation_array = np.array([[category_colors_rgb[flag_df.loc[region, ct]] 
                              for ct in cell_type_order] 
                             for region in flag_df.index])

# Create figure with more height
fig, ax = plt.subplots(figsize=(18, 4))

# Plot annotation
ax.imshow(annotation_array, aspect='auto')

# Add faint grid
ax.set_xticks(np.arange(-.5, len(cell_type_order), 1), minor=True)
ax.set_yticks(np.arange(-.5, len(flag_df.index), 1), minor=True)
ax.grid(which='minor', color='gray', linestyle='-', linewidth=0.5, alpha=0.3)

# Set labels with larger font
ax.set_yticks(range(len(flag_df.index)))
ax.set_yticklabels(flag_df.index, fontsize=12)
ax.set_xticks(range(len(cell_type_order)))
ax.set_xticklabels(cell_type_order, rotation=90, ha='right', fontsize=8)

# Add legend
legend_elements = [
    Patch(facecolor=category_colors['snRNA'], label='snRNA only'),
    Patch(facecolor=category_colors['spatial'], label='Spatial only'),
    Patch(facecolor=category_colors['both'], label='Both')
]
ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.02, 1), 
          fontsize=10, title='Detection')

plt.tight_layout()
plt.savefig('region_celltype_annotation_glia.pdf', bbox_inches='tight', dpi=300)
plt.show()


# In[18]:


df_cellect = pd.read_csv('df_cellect.csv')
df_cellect


# In[19]:


# Get the dendrogram ordering
dendro_info = adata.uns['dendrogram_detailed_celltype_annotation']
cell_type_order = dendro_info['categories_ordered']

# Pivot the dataframe to get GWAS x cell_type matrix
pivot_df = df_cellect.pivot(index='gwas', columns='annotation', values='mlog10p')

# Also pivot p-values for significance testing
pivot_pval = df_cellect.pivot(index='gwas', columns='annotation', values='pvalue')

# Reorder columns to match dendrogram order
common_celltypes = [ct for ct in cell_type_order if ct in pivot_df.columns]
pivot_df = pivot_df[common_celltypes]
pivot_pval = pivot_pval[common_celltypes]

# Bonferroni correction threshold
n_tests = 60
bonferroni_threshold = 0.05 / n_tests

# Calculate figure size to maintain square cells
n_rows = len(pivot_df.index)
n_cols = len(common_celltypes)
fig_width = 18
fig_height = fig_width * (n_rows / n_cols)  # Scale height based on number of rows

# Create figure
fig, ax = plt.subplots(figsize=(fig_width, fig_height))

# Plot heatmap with mako_r colormap
im = ax.imshow(pivot_df.values, aspect='equal', cmap='mako_r', interpolation='nearest')

# Add asterisks for significant p-values
for i in range(len(pivot_df.index)):
    for j in range(len(common_celltypes)):
        if pivot_pval.iloc[i, j] < bonferroni_threshold:
            ax.text(j, i, '*', ha='center', va='center', 
                   color='white', fontsize=16, fontweight='bold',
                   clip_on=False)

# Add faint grid
ax.set_xticks(np.arange(-.5, len(common_celltypes), 1), minor=True)
ax.set_yticks(np.arange(-.5, len(pivot_df.index), 1), minor=True)
ax.grid(which='minor', color='gray', linestyle='-', linewidth=0.5, alpha=0.3)

# Set labels
ax.set_yticks(range(len(pivot_df.index)))
ax.set_yticklabels(pivot_df.index, fontsize=10)
ax.set_xticks(range(len(common_celltypes)))
ax.set_xticklabels(common_celltypes, rotation=90, ha='right', fontsize=8)

# Add colorbar
cbar = plt.colorbar(im, ax=ax, label='-log10(p-value)', fraction=0.046, pad=0.04)

plt.tight_layout()
plt.savefig('gwas_celltype_heatmap_glia.pdf', bbox_inches='tight', dpi=300)
plt.show()


# In[21]:


# Dictionary to rename GWAS to human-readable labels
gwas_labels = {
    'BMI_Loh2018_UKBB_munged_txt': 'Body Mass Index (BMI)',
    'BMI_Yengo2018': 'Body Mass Index (BMI) - duplicate',
    'BMI_UPDATE_Yengo2018': 'Body Mass Index (BMI) - duplicate',
    'LEPR_P48357_OID20128_v1_Cardiometabolic': 'Leptin Receptor (LEPR)',
    'FLT4_P35916_OID21491_v1_Oncology': 'VEGF Receptor 3 (FLT4)',
    'DiastolicBP_GCST90310295_munged_txt': 'Diastolic Blood Pressure',
    'Pulse_GCST90310296_munged_txt': 'Pulse Rate',
    'forced_expiratory_volume_GCST90244092_munged_txt': 'Forced Expiratory Volume (FEV)'
}

# Function to get readable label
def get_gwas_label(gwas_name):
    return gwas_labels.get(gwas_name, gwas_name)

# Bonferroni correction threshold
n_tests = df_cellect['annotation'].nunique()
bonferroni_threshold = 0.05 / n_tests

print(f"\nBonferroni threshold: p < {bonferroni_threshold:.6f} (0.05 / {n_tests} cell types)")

# Identify GWAS with at least one significant cell type
df_cellect['significant'] = df_cellect['pvalue'] < bonferroni_threshold
gwas_with_sig = df_cellect[df_cellect['significant']].groupby('gwas').size().sort_values(ascending=False)

print(f"\nGWAS with ≥1 significant cell type: {len(gwas_with_sig)}")
for gwas, count in gwas_with_sig.items():
    print(f"  • {get_gwas_label(gwas)}: {count} cell types")

# Remove duplicates (keep Loh2018 as it's the largest BMI study)
gwas_to_remove = ['BMI_Yengo2018', 'BMI_UPDATE_Yengo2018']
print(f"\nRemoving duplicate BMI GWAS (keeping Loh2018):")
for g in gwas_to_remove:
    if g in gwas_with_sig.index:
        print(f"  • {get_gwas_label(g)}")

# Filter to GWAS with significance AND not in removal list
gwas_to_keep = [g for g in gwas_with_sig.index if g not in gwas_to_remove]

print(f"\nFinal GWAS to plot: {len(gwas_to_keep)}")
for g in gwas_to_keep:
    n_sig = df_cellect[(df_cellect['gwas'] == g) & (df_cellect['significant'])].shape[0]
    sig_cts = df_cellect[(df_cellect['gwas'] == g) & (df_cellect['significant'])]['annotation'].tolist()
    print(f"  • {get_gwas_label(g)} ({n_sig} significant: {', '.join(sig_cts)})")

# Filter dataframe to selected GWAS
df_filtered = df_cellect[df_cellect['gwas'].isin(gwas_to_keep)].copy()

# Pivot the dataframe
pivot_df = df_filtered.pivot(index='gwas', columns='annotation', values='mlog10p')
pivot_pval = df_filtered.pivot(index='gwas', columns='annotation', values='pvalue')

# Get dendrogram ordering from your AnnData object
dendro_info = adata.uns['dendrogram_detailed_celltype_annotation']
cell_type_order = dendro_info['categories_ordered']

# Reorder columns to match dendrogram order
common_celltypes = [ct for ct in cell_type_order if ct in pivot_df.columns]
pivot_df = pivot_df[common_celltypes]
pivot_pval = pivot_pval[common_celltypes]

# Sort GWAS by number of significant cell types (descending)
gwas_sig_counts = (pivot_pval < bonferroni_threshold).sum(axis=1).sort_values(ascending=False)
pivot_df = pivot_df.loc[gwas_sig_counts.index]
pivot_pval = pivot_pval.loc[gwas_sig_counts.index]

print(f"\nUsing dendrogram ordering ({len(common_celltypes)} cell types)")

# Calculate figure size to maintain aspect ratio
n_rows = len(pivot_df.index)
n_cols = len(common_celltypes)
fig_width = 20
fig_height = fig_width * (n_rows / n_cols)

# Ensure reasonable figure height
fig_height = max(5, min(fig_height, 10))

print(f"\nFigure dimensions: {fig_width} × {fig_height:.1f}")
print(f"Matrix: {n_rows} GWAS × {n_cols} cell types")

# Create figure
fig, ax = plt.subplots(figsize=(fig_width, fig_height))

# Plot heatmap with mako_r colormap (requires seaborn)
try:
    import seaborn as sns
    cmap = sns.color_palette("mako_r", as_cmap=True)
except:
    # Fallback if seaborn not available
    cmap = 'viridis_r'
    print("\nNote: Using viridis_r instead of mako_r (install seaborn for mako_r)")

im = ax.imshow(pivot_df.values, aspect='equal', cmap=cmap, interpolation='nearest')

# Add asterisks for significant p-values
for i in range(len(pivot_df.index)):
    for j in range(len(common_celltypes)):
        if pivot_pval.iloc[i, j] < bonferroni_threshold:
            ax.text(j, i, '*', ha='center', va='center', 
                   color='white', fontsize=20, fontweight='bold',
                   clip_on=False)

# Add faint grid
ax.set_xticks(np.arange(-.5, len(common_celltypes), 1), minor=True)
ax.set_yticks(np.arange(-.5, len(pivot_df.index), 1), minor=True)
ax.grid(which='minor', color='gray', linestyle='-', linewidth=0.5, alpha=0.3)

# Set labels with readable GWAS names
ax.set_yticks(range(len(pivot_df.index)))
readable_gwas_labels = [get_gwas_label(g) for g in pivot_df.index]
ax.set_yticklabels(readable_gwas_labels, fontsize=12, fontweight='bold')

ax.set_xticks(range(len(common_celltypes)))
ax.set_xticklabels(common_celltypes, rotation=90, ha='right', fontsize=9)

# Add colorbar
cbar = plt.colorbar(im, ax=ax, label='-log₁₀(p-value)', fraction=0.046, pad=0.04)
cbar.ax.tick_params(labelsize=10)

# Add title
ax.set_title(f'GWAS × Cell Type Enrichment (CELLECT)\n'
             f'{n_rows} GWAS Traits with Bonferroni-Significant Cell Type Associations\n'
             f'* = p < {bonferroni_threshold:.1e}', 
             fontsize=14, fontweight='bold', pad=15)

plt.tight_layout()
plt.savefig('gwas_celltype_heatmap_filtered.pdf', bbox_inches='tight', dpi=300)
plt.show()

# Print summary with readable labels
print("\n" + "="*80)
print("SUMMARY")
print("="*80)

total_sig = (pivot_pval < bonferroni_threshold).sum().sum()
total_tests = pivot_df.size

print(f"\nFiltered heatmap:")
print(f"  GWAS plotted: {n_rows}")
print(f"  Cell types: {n_cols}")
print(f"  Total cells: {total_tests}")
print(f"  Significant: {total_sig} ({100*total_sig/total_tests:.1f}%)")

print(f"\nSignificant cell types per GWAS:")
for gwas in pivot_df.index:
    n_sig = (pivot_pval.loc[gwas] < bonferroni_threshold).sum()
    sig_celltypes = pivot_pval.loc[gwas][pivot_pval.loc[gwas] < bonferroni_threshold].index.tolist()
    print(f"\n  {get_gwas_label(gwas)}: {n_sig} cell types")
    for ct in sig_celltypes:
        mlog10p = pivot_df.loc[gwas, ct]
        pval = pivot_pval.loc[gwas, ct]
        print(f"    - {ct}: -log₁₀(p) = {mlog10p:.2f} (p = {pval:.2e})")

# Save the pivot tables with readable labels
pivot_df_readable = pivot_df.copy()
pivot_df_readable.index = [get_gwas_label(g) for g in pivot_df_readable.index]
pivot_pval_readable = pivot_pval.copy()
pivot_pval_readable.index = [get_gwas_label(g) for g in pivot_pval_readable.index]

pivot_df_readable.to_csv('gwas_celltype_mlog10p_filtered.csv')
pivot_pval_readable.to_csv('gwas_celltype_pval_filtered.csv')

print("\n✓ Data matrices saved with readable GWAS labels:")
print("  • gwas_celltype_mlog10p_filtered.csv")
print("  • gwas_celltype_pval_filtered.csv")


# In[22]:


import pandas as pd
import numpy as np
from statsmodels.regression.mixed_linear_model import MixedLM
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm

# Load data
df_cellect = pd.read_csv('df_cellect.csv')
combined_LOO = pd.read_csv('df_combined_LOO_new.csv')

print("="*80)
print("STAGE 1: Identify GWAS with significant cell type enrichments")
print("="*80)

# Bonferroni correction for cell types (41 tests per GWAS)
n_celltypes = df_cellect['annotation'].nunique()
celltype_threshold = 0.05 / n_celltypes

print(f"\nNumber of cell types: {n_celltypes}")
print(f"Bonferroni threshold: 0.05 / {n_celltypes} = {celltype_threshold:.6f}")

# Find GWAS with at least one significant cell type
df_cellect['significant_celltype'] = df_cellect['pvalue'] < celltype_threshold

gwas_with_enrichment = df_cellect[df_cellect['significant_celltype']]['gwas'].unique()

print(f"\nTotal GWAS tested: {df_cellect['gwas'].nunique()}")
print(f"GWAS with ≥1 significant cell type: {len(gwas_with_enrichment)}")
print(f"\nGWAS passing Stage 1:")
for gwas in sorted(gwas_with_enrichment):
    n_sig = df_cellect[(df_cellect['gwas'] == gwas) & 
                       (df_cellect['significant_celltype'])]['annotation'].nunique()
    print(f"  {gwas}: {n_sig} significant cell type(s)")

print("\n" + "="*80)
print("STAGE 2: Test projection associations for qualified GWAS only")
print("="*80)

# Prepare the results list
results_list = []

# Define areas
area_order = ['BST', 'CEA', 'LHA', 'PAG', 'PBN', 'PVH', 'PVT']

# For each qualified GWAS and each area
for gwas in tqdm(gwas_with_enrichment, desc="Processing qualified GWAS"):
    # Get GWAS data for this trait
    gwas_data = df_cellect[df_cellect['gwas'] == gwas][['annotation', 'z']].copy()
    gwas_data.columns = ['cluster', 'z_score']

    for area in area_order:
        # Filter combined_LOO for this specific area
        area_data = combined_LOO[combined_LOO['area'] == area].copy()

        # Merge with GWAS z-scores
        merged_data = area_data.merge(gwas_data, on='cluster', how='inner')

        if len(merged_data) == 0:
            continue

        # Add intercept
        merged_data['intercept'] = 1.0

        # Run mixed model: esmu ~ z_score + (1|animal)
        try:
            model = MixedLM(
                endog=merged_data['esmu'],
                exog=merged_data[['intercept', 'z_score']],
                groups=merged_data['left_out_animal']
            )
            result = model.fit(method='lbfgs', maxiter=1000)

            # Extract results
            results_list.append({
                'gwas': gwas,
                'area': area,
                'beta_intercept': result.params['intercept'],
                'beta_z': result.params['z_score'],
                'se_intercept': result.bse['intercept'],
                'se_z': result.bse['z_score'],
                't_stat': result.tvalues['z_score'],
                'pvalue': result.pvalues['z_score'],
                'ci_lower': result.conf_int().loc['z_score', 0],
                'ci_upper': result.conf_int().loc['z_score', 1],
                'converged': result.converged,
                'n_obs': len(merged_data),
                'n_animals': merged_data['left_out_animal'].nunique(),
                'n_celltypes': merged_data['cluster'].nunique()
            })
        except Exception as e:
            print(f"Error for {gwas} - {area}: {str(e)}")

# Create DataFrame
df_stage2 = pd.DataFrame(results_list)
df_stage2['mlog10p'] = -np.log10(df_stage2['pvalue'])

print(f"\nStage 2 complete: {len(df_stage2)} tests run")
print(f"({len(gwas_with_enrichment)} GWAS × 7 areas)")

# Apply FDR correction across all Stage 2 tests
reject_fdr, pvals_fdr, _, _ = multipletests(
    df_stage2['pvalue'], 
    alpha=0.05, 
    method='fdr_bh'
)
df_stage2['significant_fdr'] = reject_fdr
df_stage2['pvalue_fdr_corrected'] = pvals_fdr

# Also add global Bonferroni for comparison
global_bonf_threshold = 0.05 / len(df_stage2)
df_stage2['significant_global_bonf'] = df_stage2['pvalue'] < global_bonf_threshold

print("\n" + "="*80)
print("MULTIPLE TESTING CORRECTION (Stage 2 tests only)")
print("="*80)
print(f"Total Stage 2 tests: {len(df_stage2)}")
print(f"Global Bonferroni (p<{global_bonf_threshold:.6f}): {df_stage2['significant_global_bonf'].sum()} significant")
print(f"Global FDR (q<0.05): {df_stage2['significant_fdr'].sum()} significant")
print(f"Percentage significant (FDR): {100 * df_stage2['significant_fdr'].sum() / len(df_stage2):.1f}%")

# Save results
df_stage2.to_csv('mixedmodel_twostage_results_glia.csv', index=False)
print(f"\nResults saved to: mixedmodel_twostage_results_glia.csv")

# Summary by area
print("\n" + "="*80)
print("FDR-SIGNIFICANT RESULTS BY AREA")
print("="*80)
pivot_sig = df_stage2.pivot(index='gwas', columns='area', values='significant_fdr')
pivot_sig = pivot_sig[area_order]
print(pivot_sig.sum().sort_values(ascending=False))

# Summary by GWAS
print("\n" + "="*80)
print("TOP GWAS BY NUMBER OF SIGNIFICANT AREAS (FDR)")
print("="*80)
print(pivot_sig.sum(axis=1).sort_values(ascending=False).head(15))

print("\n" + "="*80)
print("COMPARISON WITH ORIGINAL ANALYSIS")
print("="*80)
df_original = pd.read_csv('mixedmodel_results_by_area.csv')
print(f"Original analysis: {df_original['gwas'].nunique()} GWAS tested, {len(df_original)} total tests")
print(f"Two-stage analysis: {len(gwas_with_enrichment)} GWAS tested, {len(df_stage2)} total tests")
print(f"Reduction: {100 * (1 - len(df_stage2)/len(df_original)):.1f}% fewer tests")


# In[24]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load the two-stage results
df_stage2 = pd.read_csv('mixedmodel_twostage_results_glia.csv')

# Dictionary to rename GWAS to human-readable labels
gwas_labels = {
    'BMI_Loh2018_UKBB_munged_txt': 'Body Mass Index (BMI)',
    'BMI_Yengo2018': 'Body Mass Index (BMI) - duplicate',
    'BMI_UPDATE_Yengo2018': 'Body Mass Index (BMI) - duplicate',
    'LEPR_P48357_OID20128_v1_Cardiometabolic': 'Leptin Receptor (LEPR)',
    'FLT4_P35916_OID21491_v1_Oncology': 'VEGF Receptor 3 (FLT4)',
    'DiastolicBP_GCST90310295_munged_txt': 'Diastolic Blood Pressure',
    'Pulse_GCST90310296_munged_txt': 'Pulse Rate',
    'forced_expiratory_volume_GCST90244092_munged_txt': 'Forced Expiratory Volume (FEV)'
}

# Function to get readable label
def get_gwas_label(gwas_name):
    return gwas_labels.get(gwas_name, gwas_name)

# Filter out duplicate Yengo2018 studies (keep Loh2018)
gwas_to_remove = ['BMI_Yengo2018', 'BMI_UPDATE_Yengo2018']
df_stage2_filtered = df_stage2[~df_stage2['gwas'].isin(gwas_to_remove)].copy()

print("="*80)
print("FILTERING DUPLICATES")
print("="*80)
print(f"Original GWAS count: {df_stage2['gwas'].nunique()}")
print(f"After removing duplicates: {df_stage2_filtered['gwas'].nunique()}")
print(f"\nRemoved:")
for g in gwas_to_remove:
    if g in df_stage2['gwas'].values:
        print(f"  • {get_gwas_label(g)}")

# Apply readable labels
df_stage2_filtered['gwas_readable'] = df_stage2_filtered['gwas'].apply(get_gwas_label)

# Define area order
area_order = ['BST', 'CEA', 'LHA', 'PAG', 'PBN', 'PVH', 'PVT']

# Pivot to create GWAS x Area matrix (using readable labels)
pivot_beta = df_stage2_filtered.pivot(index='gwas_readable', columns='area', values='beta_z')
pivot_sig_fdr = df_stage2_filtered.pivot(index='gwas_readable', columns='area', values='significant_fdr')

# Reorder columns by area
pivot_beta = pivot_beta[area_order]
pivot_sig_fdr = pivot_sig_fdr[area_order]

# Sort rows by number of significant areas (descending)
row_order = pivot_sig_fdr.sum(axis=1).sort_values(ascending=False).index
pivot_beta = pivot_beta.loc[row_order]
pivot_sig_fdr = pivot_sig_fdr.loc[row_order]

# Calculate figure size to maintain square cells
n_rows = len(pivot_beta.index)
n_cols = len(area_order)
fig_width = 10
fig_height = fig_width * (n_rows / n_cols)

# Ensure reasonable figure height
fig_height = max(4, min(fig_height, 12))

# Create figure
fig, ax = plt.subplots(figsize=(fig_width, fig_height))

# Determine color scale limits (symmetric around 0)
vmax = np.abs(pivot_beta.values).max()
vmin = -vmax

# Plot heatmap with RdBu_r colormap (red=positive, blue=negative)
im = ax.imshow(pivot_beta.values, aspect='equal', cmap='RdBu_r', 
               vmin=vmin, vmax=vmax, interpolation='nearest')

# Add asterisks for FDR-significant results
for i in range(len(pivot_beta.index)):
    for j in range(len(area_order)):
        if pivot_sig_fdr.iloc[i, j]:
            ax.text(j, i, '*', ha='center', va='center', 
                   color='black', fontsize=16, fontweight='bold',
                   clip_on=False)

# Add faint grid
ax.set_xticks(np.arange(-.5, len(area_order), 1), minor=True)
ax.set_yticks(np.arange(-.5, len(pivot_beta.index), 1), minor=True)
ax.grid(which='minor', color='gray', linestyle='-', linewidth=0.5, alpha=0.3)

# Set labels
ax.set_yticks(range(len(pivot_beta.index)))
ax.set_yticklabels(pivot_beta.index, fontsize=11, fontweight='bold')
ax.set_xticks(range(len(area_order)))
ax.set_xticklabels(area_order, fontsize=12, fontweight='bold')

# Position x-axis labels at top
ax.xaxis.tick_top()
ax.xaxis.set_label_position('top')

# Add colorbar
cbar = plt.colorbar(im, ax=ax, label='Beta coefficient (effect size)', 
                    fraction=0.046, pad=0.04)

# Add title with two-stage info
title_text = ('Two-Stage Analysis: ESμ ~ GWAS z-score + (1|animal)\n'
              'Stage 1: GWAS with ≥1 significant cell type (Bonferroni)\n'
              'Stage 2: * = FDR significant projection (q < 0.05)')
ax.set_title(title_text, fontsize=11, fontweight='bold', pad=25, y=1.12)

plt.tight_layout()
plt.savefig('gwas_area_heatmap_twostage_filtered.pdf', bbox_inches='tight', dpi=300)
plt.show()

# Print summary statistics
print("\n" + "="*80)
print("TWO-STAGE HEATMAP SUMMARY (FILTERED)")
print("="*80)
print(f"\nGWAS included (passed Stage 1, duplicates removed): {len(pivot_beta.index)}")
print(f"Brain areas: {len(area_order)}")
print(f"Total tests: {len(df_stage2_filtered)}")
print(f"FDR-significant results: {pivot_sig_fdr.sum().sum()}")
print(f"Percentage significant: {100 * pivot_sig_fdr.sum().sum() / len(df_stage2_filtered):.1f}%")

print(f"\n{'Area':<6} {'Significant':>11}")
print("-"*20)
for area in area_order:
    n_sig = pivot_sig_fdr[area].sum()
    print(f"{area:<6} {n_sig:>11}")

print(f"\n{'GWAS (Readable Name)':<50} {'N Areas':>8}")
print("-"*60)
gwas_summary = pivot_sig_fdr.sum(axis=1).sort_values(ascending=False)
for gwas, n_areas in gwas_summary.items():
    gwas_short = gwas[:47] + '...' if len(gwas) > 50 else gwas
    print(f"{gwas_short:<50} {int(n_areas):>8}")

# Create a comparison plot showing before/after filtering
df_original = pd.read_csv('mixedmodel_results_by_area.csv')

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))

# Original analysis p-value distribution
ax1.hist(df_original['pvalue'], bins=50, edgecolor='black', alpha=0.7)
ax1.axvline(0.05, color='red', linestyle='--', label='p=0.05')
ax1.set_xlabel('P-value', fontsize=11)
ax1.set_ylabel('Frequency', fontsize=11)
ax1.set_title(f'Original Analysis\n({len(df_original)} tests, all GWAS)', fontsize=11)
ax1.legend()

# Two-stage analysis p-value distribution (filtered)
ax2.hist(df_stage2_filtered['pvalue'], bins=50, edgecolor='black', alpha=0.7, color='green')
ax2.axvline(0.05, color='red', linestyle='--', label='p=0.05')
ax2.set_xlabel('P-value', fontsize=11)
ax2.set_ylabel('Frequency', fontsize=11)
ax2.set_title(f'Two-Stage Analysis (Filtered)\n({len(df_stage2_filtered)} tests, no duplicates)', fontsize=11)
ax2.legend()

plt.tight_layout()
plt.savefig('pvalue_distribution_comparison_filtered.pdf', dpi=300, bbox_inches='tight')
plt.show()

# Print GWAS that were filtered
print("\n" + "="*80)
print("FILTERING DETAILS")
print("="*80)
print(f"\nOriginal tests: {len(df_stage2)}")
print(f"After filtering: {len(df_stage2_filtered)}")
print(f"Tests removed: {len(df_stage2) - len(df_stage2_filtered)}")

# Save filtered results with readable labels
df_stage2_filtered.to_csv('mixedmodel_twostage_results_filtered.csv', index=False)
print("\n✓ Filtered results saved: mixedmodel_twostage_results_filtered.csv")


# In[27]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import chi2

# Load the two-stage results
df_stage2 = pd.read_csv('mixedmodel_twostage_results_glia.csv')

# Filter to BMI_Loh2018 only (removing Yengo2018 duplicates)
bmi_gwas = 'BMI_Loh2018_UKBB_munged_txt'
df_bmi = df_stage2[df_stage2['gwas'] == bmi_gwas].copy()

print("="*80)
print(f"BMI ENRICHMENT BY AREA: {bmi_gwas}")
print("="*80)

# Define area colors
area_colors = {
    "BST": "#4D67C8",
    "CEA": "#A7EA51",
    "LHA": "#5DCEAF",
    "PAG": "#702FA0",
    "PBN": "#FFC001",
    "PVH": "#FF8021",
    "PVT": "#5DCBF3",
}

# Define area order (sorted by beta coefficient)
df_bmi_sorted = df_bmi.sort_values('beta_z')

areas = df_bmi_sorted['area'].values
betas = df_bmi_sorted['beta_z'].values
ses = df_bmi_sorted['se_z'].values
pvals = df_bmi_sorted['pvalue'].values
significant = df_bmi_sorted['significant_fdr'].values

print(f"\n{'Area':<6} {'Beta':<8} {'SE':<8} {'P-value':<12} {'FDR Sig':<10}")
print("-"*50)
for area, beta, se, pval, sig in zip(areas, betas, ses, pvals, significant):
    sig_mark = '*' if sig else ''
    print(f"{area:<6} {beta:>7.4f} {se:>7.4f} {pval:>11.2e} {sig_mark:>10}")

# Create figure
fig, ax = plt.subplots(figsize=(10, 6))

# Plot points with error bars (95% CI = ±1.96*SE)
for i, (area, beta, se, sig) in enumerate(zip(areas, betas, ses, significant)):
    # Calculate 95% confidence interval
    ci_lower = beta - 1.96 * se
    ci_upper = beta + 1.96 * se

    # Get area-specific color
    color = area_colors.get(area, '#888888')

    # Use hollow circle if not significant, filled if significant
    markerfacecolor = color if sig else 'white'

    ax.errorbar(beta, i, xerr=[[beta - ci_lower], [ci_upper - beta]], 
               fmt='o', color=color, markerfacecolor=markerfacecolor,
               markeredgecolor=color, markeredgewidth=2.5,
               markersize=14, capsize=5, capthick=2, linewidth=2)

# Add vertical line at x=0
ax.axvline(0, color='gray', linestyle='--', linewidth=1.5, alpha=0.7, zorder=0)

# Set labels and title
ax.set_xlabel('β₁ (Effect of ESμ on BMI -log₁₀(p))', fontsize=12, fontweight='bold')
ax.set_yticks(range(len(areas)))
ax.set_yticklabels(areas, fontsize=11, fontweight='bold')
ax.set_title('ESμ → BMI Enrichment: Slope by Area (Body Mass Index - Loh2018)', 
            fontsize=13, fontweight='bold', pad=15)

# Add grid
ax.grid(axis='x', alpha=0.3, linestyle=':', linewidth=0.5)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Add legend
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', 
           markeredgecolor='gray', markersize=10, label='FDR significant', 
           markeredgewidth=2),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='white', 
           markeredgecolor='gray', markersize=10, label='Not significant',
           markeredgewidth=2)
]
ax.legend(handles=legend_elements, loc='lower right', fontsize=10)

plt.tight_layout()
plt.savefig('bmi_loh2018_area_slopes.pdf', dpi=300, bbox_inches='tight')
plt.show()

# Summary statistics
print("\n" + "="*80)
print("SUMMARY")
print("="*80)
print(f"\nTotal areas: {len(areas)}")
print(f"FDR-significant: {significant.sum()}")
print(f"Range of betas: [{betas.min():.4f}, {betas.max():.4f}]")
print(f"\nAreas with positive effects: {(betas > 0).sum()}")
print(f"Areas with negative effects: {(betas < 0).sum()}")

print(f"\nStrongest positive effect: {areas[np.argmax(betas)]} (β = {betas.max():.4f})")
print(f"Strongest negative effect: {areas[np.argmin(betas)]} (β = {betas.min():.4f})")

# Test for heterogeneity across areas
# Cochran's Q test for heterogeneity
# Q = Σ(wi * (betai - beta_mean)^2) where wi = 1/se^2
weights = 1 / (ses ** 2)
weighted_mean = np.sum(weights * betas) / np.sum(weights)
Q = np.sum(weights * (betas - weighted_mean) ** 2)
df = len(betas) - 1
p_heterogeneity = 1 - chi2.cdf(Q, df)

# I² statistic (percentage of variation due to heterogeneity)
I2 = max(0, ((Q - df) / Q) * 100)

print("\n" + "="*80)
print("HETEROGENEITY TEST")
print("="*80)
print(f"Cochran's Q: {Q:.2f}")
print(f"df: {df}")
print(f"p-value: {p_heterogeneity:.2e}")
print(f"I² statistic: {I2:.1f}%")

if I2 > 75:
    print("\n→ SUBSTANTIAL heterogeneity across areas")
elif I2 > 50:
    print("\n→ MODERATE heterogeneity across areas")
else:
    print("\n→ LOW heterogeneity across areas")

print("\n✓ Figure saved: bmi_loh2018_area_slopes.pdf")


# In[9]:


# Cell 1: Imports and load data
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import gridspec

# Load the CSV files
df_single = pd.read_csv('threshold_optimization_single.csv')
df_grid = pd.read_csv('threshold_optimization_grid.csv')

# Display first few rows
print("Single threshold optimization:")
print(df_single.head())
print("\nGrid search:")
print(df_grid.head())


# In[10]:


# Cell 2: Define constants and colors
# Define colors
color_kappa = '#1f77b4'  # Blue
color_mcc = '#2ca02c'    # Green
color_f1 = '#d62728'     # Red
color_current = '#bcbcbc'  # Gray
color_optimal = 'black'

# Current and optimal thresholds
current_threshold = 0.1
optimal_threshold = 0.36

# Find closest matches (to handle floating point precision)
idx_current = (df_single['threshold'] - current_threshold).abs().idxmin()
idx_optimal = (df_single['threshold'] - optimal_threshold).abs().idxmin()

current_kappa = df_single.loc[idx_current, 'kappa']
current_agreement = df_single.loc[idx_current, 'agreement'] * 100

optimal_kappa_single = df_single.loc[idx_optimal, 'kappa']
optimal_agreement_single = df_single.loc[idx_optimal, 'agreement'] * 100

# Find optimal from grid (max kappa)
optimal_grid_idx = df_grid['kappa'].idxmax()
optimal_sc = df_grid.loc[optimal_grid_idx, 'sc_threshold']
optimal_spatial = df_grid.loc[optimal_grid_idx, 'spatial_threshold']
optimal_kappa_grid = df_grid.loc[optimal_grid_idx, 'kappa']
optimal_agreement_grid = df_grid.loc[optimal_grid_idx, 'agreement'] * 100

print(f"Current threshold: {current_threshold}")
print(f"  Kappa: {current_kappa:.3f}")
print(f"  Agreement: {current_agreement:.1f}%")
print(f"\nOptimal single threshold: {optimal_threshold}")
print(f"  Kappa: {optimal_kappa_single:.3f}")
print(f"  Agreement: {optimal_agreement_single:.1f}%")
print(f"\nOptimal grid: SC={optimal_sc:.3f}, Spatial={optimal_spatial:.3f}")
print(f"  Kappa: {optimal_kappa_grid:.3f}")


# In[17]:


# Cell 3: Panel A - Single Shared Threshold Optimization
fig, ax_a = plt.subplots(figsize=(8, 6))

# Plot metrics
ax_a.plot(df_single['threshold'], df_single['kappa'], 
         color=color_kappa, linewidth=2.5, label="Cohen's Kappa")
ax_a.plot(df_single['threshold'], df_single['mcc'], 
         color=color_mcc, linewidth=2, linestyle='--', label='MCC')
ax_a.plot(df_single['threshold'], df_single['f1'], 
         color=color_f1, linewidth=2, linestyle=':', label='F1 Score')

# Mark optimal threshold
ax_a.axvline(optimal_threshold, color=color_optimal, linestyle='--', 
            linewidth=2, alpha=0.5)
ax_a.plot(optimal_threshold, optimal_kappa_single, marker='*', 
         markersize=25, color=color_optimal, zorder=10, 
         markeredgecolor='white', markeredgewidth=2)

# Annotation for optimal
ax_a.annotate(f'Optimal\n{optimal_threshold:.3f}', 
             xy=(optimal_threshold, optimal_kappa_single),
             xytext=(optimal_threshold + 0.05, optimal_kappa_single + 0.02),
             fontsize=11, fontweight='bold',
             bbox=dict(boxstyle='round,pad=0.5', facecolor='white', 
                      edgecolor='black', linewidth=2))

ax_a.set_xlabel('Threshold (same for both datasets)', fontsize=12, fontweight='bold')
ax_a.set_ylabel('Agreement Metric', fontsize=12, fontweight='bold')
ax_a.set_title('A. Single Shared Threshold Optimization', fontsize=14, fontweight='bold', pad=10)
ax_a.legend(loc='upper left', frameon=True, fontsize=10)
ax_a.grid(alpha=0.3, linestyle='--')
ax_a.set_xlim(0, 0.5)
ax_a.spines['top'].set_visible(False)
ax_a.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig('panel_a_single_threshold.pdf', dpi=300, bbox_inches='tight')
plt.show()


# In[12]:


# Cell 4: Panel B - Independent Threshold Optimization (Heatmap)
fig, ax_b = plt.subplots(figsize=(8, 6))

# Pivot for heatmap
heatmap_data = df_grid.pivot(index='spatial_threshold', 
                             columns='sc_threshold', 
                             values='kappa')

# Create heatmap
im = ax_b.imshow(heatmap_data, aspect='auto', cmap='RdYlGn', 
                vmin=0.25, vmax=0.52, origin='lower', interpolation='nearest')

# Set ticks
n_ticks_x = 6
n_ticks_y = 6
tick_positions_x = np.linspace(0, len(heatmap_data.columns)-1, n_ticks_x).astype(int)
tick_labels_x = [f'{heatmap_data.columns[i]:.2f}' for i in tick_positions_x]
ax_b.set_xticks(tick_positions_x)
ax_b.set_xticklabels(tick_labels_x, fontsize=10)

tick_positions_y = np.linspace(0, len(heatmap_data.index)-1, n_ticks_y).astype(int)
tick_labels_y = [f'{heatmap_data.index[i]:.2f}' for i in tick_positions_y]
ax_b.set_yticks(tick_positions_y)
ax_b.set_yticklabels(tick_labels_y, fontsize=10)

ax_b.set_xlabel('Single-cell threshold', fontsize=12, fontweight='bold')
ax_b.set_ylabel('Spatial threshold', fontsize=12, fontweight='bold')
ax_b.set_title("B. Independent Threshold Optimization\n(Heatmap of Cohen's Kappa)", 
              fontsize=14, fontweight='bold', pad=10)

# Find pixel coordinates for optimal point
sc_idx = np.argmin(np.abs(heatmap_data.columns - optimal_sc))
spatial_idx = np.argmin(np.abs(heatmap_data.index - optimal_spatial))

ax_b.plot(sc_idx, spatial_idx, marker='*', markersize=30, 
         color='red', markeredgecolor='darkred', markeredgewidth=2.5, zorder=10)

# Add text box
textstr = f'★ Optimal:\nSC: {optimal_sc:.3f}\nSpatial: {optimal_spatial:.3f}\nKappa: {optimal_kappa_grid:.3f}'
ax_b.text(0.98, 0.02, textstr, transform=ax_b.transAxes,
         fontsize=10, verticalalignment='bottom', horizontalalignment='right',
         bbox=dict(boxstyle='round,pad=0.8', facecolor='white', 
                  edgecolor='orange', linewidth=3))

# Colorbar
cbar = plt.colorbar(im, ax=ax_b, fraction=0.046, pad=0.04)
cbar.set_label("Cohen's Kappa", rotation=270, labelpad=20, fontweight='bold', fontsize=11)

# Add "90th" label
ax_b.text(1.15, 0.98, '90th', transform=ax_b.transAxes, fontsize=12, 
         fontweight='bold', ha='left', va='top')

plt.tight_layout()
plt.savefig('panel_b_heatmap.pdf', dpi=300, bbox_inches='tight')
plt.show()


# In[13]:


# Cell 5: Panel C - Percentile-Based Threshold
fig, ax_c = plt.subplots(figsize=(8, 6))

# Calculate percentile-based metrics
# Use quantiles of the single-cell dataset as thresholds
percentiles = np.arange(50, 96, 5)
percentile_kappa = []
percentile_mcc = []

for p in percentiles:
    # Find threshold closest to this percentile
    threshold_at_p = df_single.iloc[(df_single['threshold'] - np.percentile(df_single['threshold'], p)).abs().argsort()[0]]
    percentile_kappa.append(threshold_at_p['kappa'])
    percentile_mcc.append(threshold_at_p['mcc'])

# Plot
ax_c.plot(percentiles, percentile_kappa, 
         color=color_kappa, linewidth=2.5, marker='o', markersize=8, label="Cohen's Kappa")
ax_c.plot(percentiles, percentile_mcc, 
         color=color_mcc, linewidth=2, marker='s', markersize=7, 
         linestyle='--', label='MCC')

# Mark 90th percentile
idx_90 = np.argmin(np.abs(percentiles - 90))
ax_c.plot(90, percentile_kappa[idx_90], marker='*', markersize=25, 
         color=color_optimal, zorder=10, markeredgecolor='white', markeredgewidth=2)

ax_c.set_xlabel('Percentile cutoff', fontsize=12, fontweight='bold')
ax_c.set_ylabel('Agreement Metric', fontsize=12, fontweight='bold')
ax_c.set_title('C. Percentile-Based Threshold', fontsize=14, fontweight='bold', pad=10)
ax_c.legend(loc='upper right', frameon=True, fontsize=10)
ax_c.grid(alpha=0.3, linestyle='--')
ax_c.set_xlim(95, 45)
ax_c.invert_xaxis()
ax_c.spines['top'].set_visible(False)
ax_c.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig('panel_c_percentile.pdf', dpi=300, bbox_inches='tight')
plt.show()


# In[18]:


# Cell 6: Panel D - Overall Agreement vs Threshold
fig, ax_d = plt.subplots(figsize=(8, 6))

ax_d.plot(df_single['threshold'], df_single['agreement'] * 100, 
         color=color_kappa, linewidth=3)

# Mark optimal
ax_d.axvline(optimal_threshold, color=color_optimal, linestyle='--', linewidth=2, alpha=0.5)
ax_d.plot(optimal_threshold, optimal_agreement_single, marker='*', 
         markersize=25, color=color_optimal, zorder=10,
         markeredgecolor='white', markeredgewidth=2)

ax_d.set_xlabel('Threshold', fontsize=12, fontweight='bold')
ax_d.set_ylabel('Agreement (%)', fontsize=12, fontweight='bold')
ax_d.set_title('D. Overall Agreement vs Threshold', fontsize=14, fontweight='bold', pad=10)
ax_d.legend(loc='lower right', frameon=True, fontsize=10)
ax_d.grid(alpha=0.3, linestyle='--')
ax_d.set_xlim(0, 0.5)
ax_d.spines['top'].set_visible(False)
ax_d.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig('panel_d_agreement.pdf', dpi=300, bbox_inches='tight')
plt.show()


# In[19]:


# Cell 7: Panel E - Dataset Balance vs Threshold
fig, ax_e = plt.subplots(figsize=(8, 6))

ax_e.plot(df_single['threshold'], df_single['n_sc_high'], 
         color=color_kappa, linewidth=2.5, label='Single-cell', marker='o', 
         markersize=4, markevery=5)
ax_e.plot(df_single['threshold'], df_single['n_spatial_high'], 
         color=color_mcc, linewidth=2.5, linestyle='--', label='Spatial',
         marker='s', markersize=4, markevery=5)

ax_e.axvline(optimal_threshold, color='red', linestyle='--', 
            linewidth=2, alpha=0.7, label='Optimal')

ax_e.set_xlabel('Threshold', fontsize=12, fontweight='bold')
ax_e.set_ylabel('Number of "high" projectors', fontsize=12, fontweight='bold')
ax_e.set_title('E. Dataset Balance vs Threshold', fontsize=14, fontweight='bold', pad=10)
ax_e.legend(loc='upper right', frameon=True, fontsize=10)
ax_e.grid(alpha=0.3, linestyle='--')
ax_e.set_xlim(0, 0.5)
ax_e.spines['top'].set_visible(False)
ax_e.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig('panel_e_balance.pdf', dpi=300, bbox_inches='tight')
plt.show()


# In[16]:


# Cell 8: Summary Statistics
# Get values at optimal threshold
n_sc_opt = df_single.loc[idx_optimal, 'n_sc_high']
n_spatial_opt = df_single.loc[idx_optimal, 'n_spatial_high']

print("="*60)
print("THRESHOLD OPTIMIZATION SUMMARY")
print("="*60)
print(f"\nCurrent (0.1 for both):")
print(f"  Kappa:      {current_kappa:.3f}")
print(f"  Agreement:  {current_agreement:.1f}%")
print(f"\nOption 1 - Single Threshold ({optimal_threshold}):")
print(f"  Kappa:      {optimal_kappa_single:.3f} (+{optimal_kappa_single - current_kappa:.3f})")
print(f"  Agreement:  {optimal_agreement_single:.1f}% (+{optimal_agreement_single - current_agreement:.1f}%)")
print(f"  ✓ Simple, interpretable")
print(f"\nOption 2 - Independent Thresholds:")
print(f"  SC:         {optimal_sc:.3f}")
print(f"  Spatial:    {optimal_spatial:.3f}")
print(f"  Kappa:      {optimal_kappa_grid:.3f} (+{optimal_kappa_grid - current_kappa:.3f})")
print(f"  Agreement:  {optimal_agreement_grid:.1f}% (+{optimal_agreement_grid - current_agreement:.1f}%)")
print(f"  ✓ Maximum agreement")
print(f"\nRECOMMENDATION:")
print(f"Use threshold = {optimal_threshold} for both")
print(f"(or SC={optimal_sc:.2f}, Spatial={optimal_spatial:.2f})")
print(f"\nThis achieves:")
print(f"• {optimal_agreement_single:.0f}% agreement (vs {current_agreement:.0f}%)")
print(f"• Kappa {optimal_kappa_single:.2f} (vs {current_kappa:.2f})")
print(f"• ~{n_sc_opt}-{n_spatial_opt} 'high' per dataset (balanced)")
print("="*60)


# In[5]:


del adata
adata = sc.read("allgenes_bipca_annotated.h5ad", cache=True)
adata


# In[6]:


# Cell 1: Load and explore the data
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.stats import fisher_exact, chi2_contingency
import scanpy as sc

# Check the data
print("Available regions:")
print(adata.obs['region_collapsed'].value_counts())
print("\nScore distribution:")
print(adata.obs['score_region_collapsed'].describe())
print(f"\nNumber of cell types: {adata.obs['detailed_celltype_annotation'].nunique()}")

# Visualize score distribution
fig, axes = plt.subplots(1, 2, figsize=(12, 4))
axes[0].hist(adata.obs['score_region_collapsed'], bins=50, edgecolor='black')
axes[0].axvline(0.4, color='red', linestyle='--', label='Current cutoff (0.4)')
axes[0].set_xlabel('Region assignment score')
axes[0].set_ylabel('Frequency')
axes[0].set_title('Distribution of region assignment scores')
axes[0].legend()

# Score by region
adata.obs.boxplot(column='score_region_collapsed', by='region_collapsed', ax=axes[1], rot=45)
axes[1].set_title('Score distribution by region')
axes[1].set_xlabel('Region')
plt.suptitle('')
plt.tight_layout()
plt.show()


# In[7]:


# Cell 2: Benchmark different score cutoffs
def compute_enrichment_at_cutoff(adata, score_cutoff, exclude_regions=['undefined', 'NA']):
    """
    Compute cell type enrichment in regions at a given score cutoff.
    Returns contingency tables and enrichment statistics.
    """
    # Filter high-confidence assignments
    mask_highconf = adata.obs['score_region_collapsed'] > score_cutoff
    mask_valid_region = ~adata.obs['region_collapsed'].isin(exclude_regions)
    mask = mask_highconf & mask_valid_region

    # Get filtered data
    df = adata.obs[mask][['detailed_celltype_annotation', 'region_collapsed']].copy()

    # Compute contingency table
    contingency = pd.crosstab(
        df['detailed_celltype_annotation'],
        df['region_collapsed']
    )

    # Compute proportions (normalize by region size)
    proportions = contingency.div(contingency.sum(axis=0), axis=1)

    # Compute expected proportions (if no enrichment)
    expected_prop = contingency.sum(axis=1) / contingency.sum().sum()

    # Compute fold-enrichment
    fold_enrichment = proportions.div(expected_prop, axis=0)

    # Compute statistics for each cell type - region pair
    results = []
    for celltype in contingency.index:
        for region in contingency.columns:
            # 2x2 contingency table for Fisher's exact test
            a = contingency.loc[celltype, region]  # celltype in region
            b = contingency.loc[celltype, :].sum() - a  # celltype in other regions
            c = contingency.loc[:, region].sum() - a  # other celltypes in region
            d = contingency.sum().sum() - a - b - c  # other celltypes in other regions

            if a + b > 0:  # Only test if celltype exists
                oddsratio, pval = fisher_exact([[a, b], [c, d]], alternative='greater')

                results.append({
                    'celltype': celltype,
                    'region': region,
                    'n_observed': a,
                    'n_celltype_total': a + b,
                    'n_region_total': a + c,
                    'proportion_in_region': a / (a + c) if (a + c) > 0 else 0,
                    'proportion_expected': (a + b) / contingency.sum().sum(),
                    'fold_enrichment': fold_enrichment.loc[celltype, region],
                    'odds_ratio': oddsratio,
                    'pvalue': pval,
                    'score_cutoff': score_cutoff
                })

    df_results = pd.DataFrame(results)

    # FDR correction
    from statsmodels.stats.multitest import multipletests
    if len(df_results) > 0:
        _, df_results['pvalue_adj'], _, _ = multipletests(
            df_results['pvalue'], 
            method='fdr_bh'
        )

    return df_results, contingency, proportions

# Test multiple cutoffs
cutoffs = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
cutoff_results = {}

for cutoff in cutoffs:
    df_enrich, contingency, proportions = compute_enrichment_at_cutoff(adata, cutoff)
    cutoff_results[cutoff] = {
        'enrichment': df_enrich,
        'contingency': contingency,
        'proportions': proportions,
        'n_cells': len(adata.obs[adata.obs['score_region_collapsed'] > cutoff]),
        'n_significant': (df_enrich['pvalue_adj'] < 0.05).sum()
    }
    print(f"Cutoff {cutoff}: {cutoff_results[cutoff]['n_cells']} cells, "
          f"{cutoff_results[cutoff]['n_significant']} significant enrichments")


# In[8]:


# Cell 3: Visualize cutoff benchmark
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Panel A: Number of cells retained
ax = axes[0, 0]
n_cells = [cutoff_results[c]['n_cells'] for c in cutoffs]
ax.plot(cutoffs, n_cells, marker='o', linewidth=2, markersize=8)
ax.axvline(0.4, color='red', linestyle='--', alpha=0.5, label='Current (0.4)')
ax.set_xlabel('Score cutoff', fontweight='bold')
ax.set_ylabel('Number of cells retained', fontweight='bold')
ax.set_title('A. Cell retention vs cutoff', fontweight='bold')
ax.grid(alpha=0.3)
ax.legend()

# Panel B: Number of significant enrichments
ax = axes[0, 1]
n_sig = [cutoff_results[c]['n_significant'] for c in cutoffs]
ax.plot(cutoffs, n_sig, marker='o', linewidth=2, markersize=8, color='green')
ax.axvline(0.4, color='red', linestyle='--', alpha=0.5, label='Current (0.4)')
ax.set_xlabel('Score cutoff', fontweight='bold')
ax.set_ylabel('Number of significant enrichments (FDR<0.05)', fontweight='bold')
ax.set_title('B. Statistical power vs cutoff', fontweight='bold')
ax.grid(alpha=0.3)
ax.legend()

# Panel C: Mean fold-enrichment of significant hits
ax = axes[1, 0]
mean_fold = []
for c in cutoffs:
    sig_hits = cutoff_results[c]['enrichment'][cutoff_results[c]['enrichment']['pvalue_adj'] < 0.05]
    if len(sig_hits) > 0:
        mean_fold.append(sig_hits['fold_enrichment'].mean())
    else:
        mean_fold.append(0)
ax.plot(cutoffs, mean_fold, marker='o', linewidth=2, markersize=8, color='purple')
ax.axvline(0.4, color='red', linestyle='--', alpha=0.5, label='Current (0.4)')
ax.set_xlabel('Score cutoff', fontweight='bold')
ax.set_ylabel('Mean fold-enrichment', fontweight='bold')
ax.set_title('C. Enrichment strength vs cutoff', fontweight='bold')
ax.grid(alpha=0.3)
ax.legend()

# Panel D: Tradeoff plot
ax = axes[1, 1]
ax.scatter(n_cells, n_sig, s=100, c=cutoffs, cmap='viridis', edgecolors='black', linewidth=1.5)
for i, c in enumerate(cutoffs):
    ax.annotate(f'{c}', (n_cells[i], n_sig[i]), 
               xytext=(5, 5), textcoords='offset points', fontsize=9)
# Mark 0.4
idx_04 = cutoffs.index(0.4)
ax.scatter(n_cells[idx_04], n_sig[idx_04], s=200, marker='*', 
          color='red', edgecolors='darkred', linewidth=2, zorder=10)
ax.set_xlabel('Number of cells retained', fontweight='bold')
ax.set_ylabel('Number of significant enrichments', fontweight='bold')
ax.set_title('D. Power-Sensitivity tradeoff', fontweight='bold')
ax.grid(alpha=0.3)
cbar = plt.colorbar(ax.collections[0], ax=ax)
cbar.set_label('Score cutoff', rotation=270, labelpad=15)

plt.tight_layout()
plt.savefig('cutoff_benchmark.pdf', dpi=300, bbox_inches='tight')
plt.show()

print("\nRecommended cutoff based on balance:")
print(f"Current (0.4): {cutoff_results[0.4]['n_cells']} cells, {cutoff_results[0.4]['n_significant']} enrichments")


# In[23]:


# Cell 4: Use optimal cutoff for detailed analysis
optimal_cutoff = 0.4  # Adjust based on Cell 3 results

df_enrichment, contingency, proportions = compute_enrichment_at_cutoff(adata, optimal_cutoff)

# Filter significant results
df_sig = df_enrichment[df_enrichment['pvalue_adj'] < 0.05].copy()
df_sig = df_sig.sort_values('fold_enrichment', ascending=False)

print(f"Significant enrichments (FDR < 0.05): {len(df_sig)}")
print("\nTop 20 enrichments:")
print(df_sig[['celltype', 'region', 'n_observed', 'fold_enrichment', 'pvalue_adj']].head(20))

# Check glia/pituitary as negative controls
glia_pituitary = ['Glia', 'Pituitary']  # Adjust names as needed
glia_results = df_enrichment[df_enrichment['celltype'].str.contains('|'.join(glia_pituitary), case=False)]
print(f"\n\nGlia/Pituitary enrichments (should be minimal):")
print(glia_results[glia_results['pvalue_adj'] < 0.05])


# In[24]:


# Cell 5: Heatmap of fold-enrichment (rotated 90°, regions on Y-axis)
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Get the dendrogram ordering from the neurotransmitter plot
# (This is the order from top to bottom in your uploaded PDF)
dendrogram_order = [
    '17_Dkk2_Ebf3', '42_Foxb1_Myo5b', '09_Nfia_Ccbe1', '18_Ust_Bcl11b', 
    '23_Ust_Nkd2', '05_Qrfpr_Pdzrn3', '02_Fezf1_Lama1', '00_Hmcn2_Gpr149',
    '03_Rai14_Acvr1c', '34_Pla2r1_Sim1', '04_Crnde_Lmx1a', '32_Gucy2c_Ttc6',
    '35_Lhx9_Zic1', '41_Otx1_Zic1', '27_Lmx1a_Ttc6', '36_Ebf2_Crhbp',
    '13_Pth2r_Otp', '07_Nfib_Six3', '29_Dnah11_Pth2r', '25_Nxph4_Ebf3',
    '30_Trh_Sim1', '12_Ntng2_Onecut1', '19_Bace2_Tspan18', '39_Nfib_Nfia',
    '16_Gal_Sox6', '24_Nts_Ecel1', '01_Undefined', '15_Onecut1_Slc32a1',
    '20_Tmem114_Npsr1', '26_Dach2_Npsr1', '08_Lef1_Rspo2', '38_Hdc_Prph',
    '31_Lhx9_Npvf', '37_Qrfp_Qrfprl', '11_Tafa4_Kank4', '06_Tbx3_Npy1r',
    '28_Kiss1_Slc38a3', '14_Lhx1_Th', '40_Gpr50_Tac2', '22_Crabp1_Vgll3',
    '10_Npy_Agrp', '21_Bnc2_Glp1r', '33_Chat_Ghrh',
    'Pituitary_Mtnr1a', 'Pituitary_Crhr1', 'Pituitary_Pax7', 
    'Pituitary_Gnrhr', 'Pituitary_Pappa2', 'Pituitary_Rxfp1',
    'OPC', 'Microglia', 'Oligodendrocyte', 'Undefined', 'Astrocyte',
    'VLMC', 'Epithelial-like_cell', 'Schwann_cell', 'Endothelial_cell',
    'Pericyte', 'Pituicyte', 'Ependymal_cell', 'Tanycyte'
]

# Pivot to wide format for heatmap (TRANSPOSED: regions as rows, celltypes as columns)
pivot_fold = df_enrichment.pivot(index='region', columns='celltype', values='fold_enrichment')
pivot_pval = df_enrichment.pivot(index='region', columns='celltype', values='pvalue_adj')

# Reorder columns to match dendrogram order
common_celltypes = [ct for ct in dendrogram_order if ct in pivot_fold.columns]
pivot_fold = pivot_fold[common_celltypes]
pivot_pval = pivot_pval[common_celltypes]

# Create significance mask (FDR < 0.05)
sig_mask = pivot_pval < 0.05

# Calculate figure size to maintain square cells (now transposed)
n_rows = len(pivot_fold.index)  # regions
n_cols = len(pivot_fold.columns)  # celltypes
fig_height = 8  # Base height for regions
fig_width = fig_height * (n_cols / n_rows)  # Scale width for square cells

# Create figure
fig, ax = plt.subplots(figsize=(fig_width, fig_height))

# Plot heatmap with colormap centered at 1
im = ax.imshow(pivot_fold.values, aspect='equal', cmap='RdPu', 
               vmin=0, vmax=5, interpolation='nearest')

# Add asterisks for significant enrichments
for i in range(len(pivot_fold.index)):
    for j in range(len(pivot_fold.columns)):
        if sig_mask.iloc[i, j]:
            ax.text(j, i, '*', ha='center', va='center', 
                   color='white', fontsize=10, fontweight='bold',
                   clip_on=False)

# Add faint grid
ax.set_xticks(np.arange(-.5, len(pivot_fold.columns), 1), minor=True)
ax.set_yticks(np.arange(-.5, len(pivot_fold.index), 1), minor=True)
ax.grid(which='minor', color='gray', linestyle='-', linewidth=0.5, alpha=0.3)

# Set labels (X-axis: cell types, Y-axis: regions)
ax.set_xticks(range(len(pivot_fold.columns)))
ax.set_xticklabels(pivot_fold.columns, rotation=90, ha='right', fontsize=7)
ax.set_yticks(range(len(pivot_fold.index)))
ax.set_yticklabels(pivot_fold.index, fontsize=9)

# Labels
ax.set_xlabel('Cell Type (dendrogram ordered)', fontweight='bold', fontsize=11)
ax.set_ylabel('Region', fontweight='bold', fontsize=11)
ax.set_title(f'Cell Type Enrichment by Region (score > {optimal_cutoff})\n* = FDR < 0.05', 
            fontweight='bold', fontsize=13, pad=10)

# Add colorbar
cbar = plt.colorbar(im, ax=ax, label='Fold enrichment', fraction=0.046, pad=0.04)

plt.tight_layout()
plt.savefig('enrichment_heatmap_all_rotated.pdf', dpi=300, bbox_inches='tight')
plt.show()


# In[25]:


# Cell 6: Focused heatmap - only significant enrichments (rotated 90°)
# Filter to cell types with at least one significant enrichment
celltypes_with_sig = df_sig['celltype'].unique()

# Filter pivot tables to only significant celltypes
common_celltypes_sig = [ct for ct in dendrogram_order if ct in celltypes_with_sig]
pivot_fold_filtered = pivot_fold[common_celltypes_sig]
pivot_pval_filtered = pivot_pval[common_celltypes_sig]
sig_mask_filtered = pivot_pval_filtered < 0.05

# Calculate figure size for square cells
n_rows = len(pivot_fold_filtered.index)  # regions
n_cols = len(pivot_fold_filtered.columns)  # celltypes
fig_height = 8
fig_width = fig_height * (n_cols / n_rows)

# Create figure
fig, ax = plt.subplots(figsize=(fig_width, fig_height))

# Plot heatmap
im = ax.imshow(pivot_fold_filtered.values, aspect='equal', cmap='RdBu_r', 
               vmin=0, vmax=5, interpolation='nearest')

# Add asterisks for significant enrichments
for i in range(len(pivot_fold_filtered.index)):
    for j in range(len(pivot_fold_filtered.columns)):
        if sig_mask_filtered.iloc[i, j]:
            ax.text(j, i, '*', ha='center', va='center', 
                   color='black', fontsize=10, fontweight='bold',
                   clip_on=False)

# Add faint grid
ax.set_xticks(np.arange(-.5, len(pivot_fold_filtered.columns), 1), minor=True)
ax.set_yticks(np.arange(-.5, len(pivot_fold_filtered.index), 1), minor=True)
ax.grid(which='minor', color='gray', linestyle='-', linewidth=0.5, alpha=0.3)

# Set labels
ax.set_xticks(range(len(pivot_fold_filtered.columns)))
ax.set_xticklabels(pivot_fold_filtered.columns, rotation=90, ha='right', fontsize=7)
ax.set_yticks(range(len(pivot_fold_filtered.index)))
ax.set_yticklabels(pivot_fold_filtered.index, fontsize=9)

# Labels
ax.set_xlabel('Cell Type (dendrogram ordered)', fontweight='bold', fontsize=11)
ax.set_ylabel('Region', fontweight='bold', fontsize=11)
ax.set_title(f'Significant Cell Type Enrichments (score > {optimal_cutoff})\n* = FDR < 0.05', 
            fontweight='bold', fontsize=13, pad=10)

# Add colorbar
cbar = plt.colorbar(im, ax=ax, label='Fold enrichment', fraction=0.046, pad=0.04)

plt.tight_layout()
plt.savefig('enrichment_heatmap_significant_rotated.pdf', dpi=300, bbox_inches='tight')
plt.show()


# In[12]:


# Cell 7: Region-specific cell type composition
# For each region, show top enriched cell types
fig, axes = plt.subplots(1, len(proportions.columns), figsize=(4*len(proportions.columns), 6))
if len(proportions.columns) == 1:
    axes = [axes]

for idx, region in enumerate(proportions.columns):
    # Get top cell types in this region
    top_celltypes = proportions[region].nlargest(15)

    # Get enrichment significance
    region_enrichments = df_enrichment[df_enrichment['region'] == region].set_index('celltype')
    colors = ['red' if region_enrichments.loc[ct, 'pvalue_adj'] < 0.05 else 'gray' 
             for ct in top_celltypes.index]

    axes[idx].barh(range(len(top_celltypes)), top_celltypes.values, color=colors)
    axes[idx].set_yticks(range(len(top_celltypes)))
    axes[idx].set_yticklabels(top_celltypes.index, fontsize=9)
    axes[idx].set_xlabel('Proportion', fontweight='bold')
    axes[idx].set_title(f'{region}\n(red = enriched)', fontweight='bold')
    axes[idx].invert_yaxis()
    axes[idx].spines['top'].set_visible(False)
    axes[idx].spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig('region_composition.pdf', dpi=300, bbox_inches='tight')
plt.show()


# In[13]:


# Cell 8: Export results to CSV
# Save full enrichment table
df_enrichment.to_csv('celltype_region_enrichment_all.csv', index=False)

# Save significant enrichments only
df_sig.to_csv('celltype_region_enrichment_significant.csv', index=False)

# Save contingency tables
contingency.to_csv('celltype_region_contingency.csv')
proportions.to_csv('celltype_region_proportions.csv')

print("Results saved:")
print("  - celltype_region_enrichment_all.csv")
print("  - celltype_region_enrichment_significant.csv")
print("  - celltype_region_contingency.csv")
print("  - celltype_region_proportions.csv")

# Summary statistics
print(f"\n\nSummary at cutoff {optimal_cutoff}:")
print(f"Total cells analyzed: {len(adata.obs[adata.obs['score_region_collapsed'] > optimal_cutoff])}")
print(f"Total cell types: {df_enrichment['celltype'].nunique()}")
print(f"Total regions: {df_enrichment['region'].nunique()}")
print(f"Total enrichments tested: {len(df_enrichment)}")
print(f"Significant enrichments (FDR<0.05): {len(df_sig)} ({100*len(df_sig)/len(df_enrichment):.1f}%)")
print(f"\nCell types with strongest enrichments:")
for _, row in df_sig.head(10).iterrows():
    print(f"  {row['celltype']:<30} → {row['region']:<15} (fold={row['fold_enrichment']:.2f}, p={row['pvalue_adj']:.2e})")


# In[13]:


# Cell 1: Explore the crossed experimental design
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.stats import spearmanr, pearsonr
import scanpy as sc

print("="*80)
print("CROSSED EXPERIMENTAL DESIGN")
print("="*80)

# Get v147 expression from normalized layer
v147_expr = adata[:, 'v147'].layers['scran_normalization'].toarray().flatten() \
    if hasattr(adata[:, 'v147'].layers['scran_normalization'], 'toarray') \
    else adata[:, 'v147'].layers['scran_normalization'].flatten()

adata.obs['v147_expr'] = v147_expr
adata.obs['v147_positive'] = (v147_expr > 0).astype(int)

# Build a mapping of animal-reaction combinations
animal_reaction_data = []

for idx, row in adata.obs.iterrows():
    animal_reaction_data.append({
        'animal': row['hash.mcl.ID'],
        'reaction': row['orig.ident'],
        'sort_prob': row['sort'],
        'v147_expr': row['v147_expr'],
        'v147_positive': row['v147_positive']
    })

df_cells = pd.DataFrame(animal_reaction_data)

# Check the design
print("\n1. Reactions per animal:")
reactions_per_animal = df_cells.groupby('animal')['reaction'].nunique()
print(f"   Mean: {reactions_per_animal.mean():.1f}")
print(f"   Mode: {reactions_per_animal.mode()[0]}")
print(f"   Range: {reactions_per_animal.min()}-{reactions_per_animal.max()}")

print("\n2. Animals per reaction:")
animals_per_reaction = df_cells.groupby('reaction')['animal'].nunique()
print(f"   Mean: {animals_per_reaction.mean():.1f}")
print(f"   Range: {animals_per_reaction.min()}-{animals_per_reaction.max()}")

print("\n3. Total structure:")
print(f"   Animals: {df_cells['animal'].nunique()}")
print(f"   Reactions: {df_cells['reaction'].nunique()}")
print(f"   Animal-Reaction combinations: {df_cells.groupby(['animal', 'reaction']).ngroups}")

# Check if each animal sees different sort probabilities
print("\n4. Sort probability variation per animal:")
sort_by_animal = df_cells.groupby('animal')['sort_prob'].agg(['nunique', 'min', 'max', 'std'])
print(sort_by_animal.describe())

print("\n5. Example animals across reactions:")
example_animals = df_cells['animal'].value_counts().head(5).index
for animal in example_animals:
    animal_data = df_cells[df_cells['animal'] == animal].groupby('reaction')['sort_prob'].first()
    print(f"   Animal {animal}: {len(animal_data)} reactions, sort probs: {animal_data.values}")


# In[14]:


# Cell 2: Aggregate to animal-reaction level
# This is now our unit of analysis: each animal in each reaction

animal_reaction_summary = []

for (animal, reaction), group in df_cells.groupby(['animal', 'reaction']):
    animal_reaction_summary.append({
        'animal': animal,
        'reaction': reaction,
        'sort_prob': group['sort_prob'].iloc[0],  # Constant within reaction
        'n_cells': len(group),
        'n_v147_pos': (group['v147_positive'] == 1).sum(),
        'pct_v147_pos': 100 * (group['v147_positive'] == 1).sum() / len(group),
        'mean_v147_all': group['v147_expr'].mean(),
        'mean_v147_pos': group[group['v147_positive'] == 1]['v147_expr'].mean() 
                         if (group['v147_positive'] == 1).sum() > 0 else 0,
    })

df_animal_reaction = pd.DataFrame(animal_reaction_summary)

print("\n" + "="*80)
print("ANIMAL-REACTION LEVEL SUMMARY")
print("="*80)
print(f"\nTotal animal-reaction combinations: {len(df_animal_reaction)}")
print(f"Expected if complete design: {df_cells['animal'].nunique() * 3} (58 animals × 3 reactions)")

print("\nAnimal-reaction level statistics:")
print(df_animal_reaction.describe())

# Save
df_animal_reaction.to_csv('sun1_v147_animal_reaction_level.csv', index=False)
print("\n✓ Animal-reaction level data saved")


# In[15]:


# Cell 3: Within-animal analysis
# For each animal, look at correlation between sort_prob and v147 expression

within_animal_results = []

for animal in df_animal_reaction['animal'].unique():
    animal_data = df_animal_reaction[df_animal_reaction['animal'] == animal]

    if len(animal_data) < 2:
        continue

    # Check if there's variation in sort probability for this animal
    if animal_data['sort_prob'].std() < 0.001:
        continue  # Skip if no variation in sort

    # Correlations
    r_pct, p_pct = spearmanr(animal_data['sort_prob'], animal_data['pct_v147_pos'])
    r_mean, p_mean = spearmanr(animal_data['sort_prob'], animal_data['mean_v147_all'])

    within_animal_results.append({
        'animal': animal,
        'n_reactions': len(animal_data),
        'sort_range': animal_data['sort_prob'].max() - animal_data['sort_prob'].min(),
        'r_pct_v147': r_pct,
        'p_pct_v147': p_pct,
        'r_mean_v147': r_mean,
        'p_mean_v147': p_mean,
        'mean_pct_v147': animal_data['pct_v147_pos'].mean(),
    })

df_within_animal = pd.DataFrame(within_animal_results)

print("\n" + "="*80)
print("WITHIN-ANIMAL CORRELATIONS")
print("="*80)
print(f"\nAnimals with sufficient variation: {len(df_within_animal)}")

if len(df_within_animal) > 0:
    print("\nWithin-animal correlation distribution (% v147+):")
    print(df_within_animal['r_pct_v147'].describe())

    # Meta-analysis: average correlation across animals
    mean_r = df_within_animal['r_pct_v147'].mean()
    print(f"\nMean within-animal correlation: {mean_r:.3f}")

    # Fisher's z-transformation for proper averaging
    z_scores = np.arctanh(df_within_animal['r_pct_v147'].values)
    mean_z = z_scores.mean()
    mean_r_fisher = np.tanh(mean_z)
    print(f"Mean correlation (Fisher's z): {mean_r_fisher:.3f}")

    # One-sample t-test: is mean correlation > 0?
    from scipy.stats import ttest_1samp
    t_stat, p_val = ttest_1samp(z_scores, 0)
    print(f"One-sample t-test (H0: r=0): t={t_stat:.3f}, p={p_val:.4e}")

    df_within_animal.to_csv('sun1_v147_within_animal_correlations.csv', index=False)
    print("\n✓ Within-animal correlations saved")


# In[16]:


# Cell 4: Visualization of within-animal effects
if len(df_within_animal) > 0:
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    # Panel A: Distribution of within-animal correlations
    ax = axes[0, 0]
    ax.hist(df_within_animal['r_pct_v147'], bins=20, edgecolor='black', alpha=0.7, color='steelblue')
    ax.axvline(0, color='red', linestyle='--', linewidth=2, label='No correlation')
    ax.axvline(df_within_animal['r_pct_v147'].mean(), color='green', linestyle='--', 
              linewidth=2, label=f'Mean r={df_within_animal["r_pct_v147"].mean():.3f}')
    ax.set_xlabel('Within-animal Spearman ρ', fontweight='bold', fontsize=12)
    ax.set_ylabel('Number of animals', fontweight='bold', fontsize=12)
    ax.set_title('A. Distribution of within-animal correlations\n(Sun1 prob vs % v147+)',
                fontweight='bold', fontsize=12)
    ax.legend()
    ax.grid(alpha=0.3, axis='y')

    # Panel B: Example animals with strongest positive correlations
    ax = axes[0, 1]
    top_animals = df_within_animal.nlargest(6, 'r_pct_v147')['animal'].values

    for i, animal in enumerate(top_animals[:6]):
        animal_data = df_animal_reaction[df_animal_reaction['animal'] == animal]
        ax.plot(animal_data['sort_prob'], animal_data['pct_v147_pos'], 
               marker='o', linewidth=2, markersize=8, label=f'Animal {animal}', alpha=0.7)

    ax.set_xlabel('Sun1 sorting probability', fontweight='bold', fontsize=12)
    ax.set_ylabel('% v147+ cells', fontweight='bold', fontsize=12)
    ax.set_title('B. Top 6 animals by positive correlation',
                fontweight='bold', fontsize=12)
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)

    # Panel C: Correlation vs sort range
    ax = axes[1, 0]
    scatter = ax.scatter(df_within_animal['sort_range'], df_within_animal['r_pct_v147'],
                        s=100, alpha=0.6, edgecolors='black')
    ax.set_xlabel('Sort probability range (per animal)', fontweight='bold', fontsize=12)
    ax.set_ylabel('Within-animal correlation', fontweight='bold', fontsize=12)
    ax.set_title('C. Correlation strength vs sort probability range',
                fontweight='bold', fontsize=12)
    ax.axhline(0, color='red', linestyle='--', alpha=0.5)
    ax.grid(alpha=0.3)

    # Panel D: Correlation vs mean v147+ percentage
    ax = axes[1, 1]
    scatter = ax.scatter(df_within_animal['mean_pct_v147'], df_within_animal['r_pct_v147'],
                        s=100, alpha=0.6, edgecolors='black', c=df_within_animal['n_reactions'],
                        cmap='viridis')
    ax.set_xlabel('Mean % v147+ (per animal)', fontweight='bold', fontsize=12)
    ax.set_ylabel('Within-animal correlation', fontweight='bold', fontsize=12)
    ax.set_title('D. Correlation vs overall v147+ rate',
                fontweight='bold', fontsize=12)
    ax.axhline(0, color='red', linestyle='--', alpha=0.5)
    ax.grid(alpha=0.3)
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('Number of reactions', rotation=270, labelpad=15)

    plt.tight_layout()
    plt.savefig('sun1_v147_within_animal_analysis.pdf', dpi=300, bbox_inches='tight')
    plt.show()


# In[17]:


# Cell 5: Mixed-effects model with crossed random effects
# Animals and Reactions are both random effects

try:
    import statsmodels.api as sm
    from statsmodels.regression.mixed_linear_model import MixedLM

    # Center predictors
    df_animal_reaction['sort_centered'] = df_animal_reaction['sort_prob'] - df_animal_reaction['sort_prob'].mean()

    print("\n" + "="*80)
    print("CROSSED RANDOM EFFECTS MODEL")
    print("="*80)
    print("\nModel: pct_v147_pos ~ sort_prob + (1|animal) + (1|reaction)")
    print("Note: Reaction random effect captures technical variation")
    print("      Animal random effect captures biological variation")

    # Model 1: % v147+ with animal random intercept only (simpler)
    print("\n" + "-"*80)
    print("MODEL 1: Random intercept for animals only")
    print("-"*80)

    md1 = MixedLM.from_formula('pct_v147_pos ~ sort_centered',
                               data=df_animal_reaction,
                               groups=df_animal_reaction['animal'])
    mdf1 = md1.fit(method='lbfgs')
    print(mdf1.summary())

    # Model 2: Mean v147 expression (all cells)
    print("\n" + "-"*80)
    print("MODEL 2: Mean v147 expression (all cells)")
    print("-"*80)

    md2 = MixedLM.from_formula('mean_v147_all ~ sort_centered',
                               data=df_animal_reaction,
                               groups=df_animal_reaction['animal'])
    mdf2 = md2.fit(method='lbfgs')
    print(mdf2.summary())

    # Extract results
    print("\n" + "="*80)
    print("MIXED-EFFECTS RESULTS SUMMARY")
    print("="*80)

    print("\n1. % v147+ cells:")
    print(f"   Fixed effect (sort): β = {mdf1.params['sort_centered']:.4f}")
    print(f"   Standard error: {mdf1.bse['sort_centered']:.4f}")
    print(f"   p-value: {mdf1.pvalues['sort_centered']:.4e}")
    print(f"   95% CI: [{mdf1.conf_int().loc['sort_centered', 0]:.4f}, {mdf1.conf_int().loc['sort_centered', 1]:.4f}]")
    print(f"   Interpretation: {mdf1.params['sort_centered']:.2f}% change in v147+ per unit Sun1 probability")

    print("\n2. Mean v147 expression (all cells):")
    print(f"   Fixed effect (sort): β = {mdf2.params['sort_centered']:.6f}")
    print(f"   Standard error: {mdf2.bse['sort_centered']:.6f}")
    print(f"   p-value: {mdf2.pvalues['sort_centered']:.4e}")
    print(f"   95% CI: [{mdf2.conf_int().loc['sort_centered', 0]:.6f}, {mdf2.conf_int().loc['sort_centered', 1]:.6f}]")

    # Variance components
    print("\n3. Variance components:")
    print(f"   Between-animal variance (% v147+): {mdf1.cov_re.values[0][0]:.4f}")
    print(f"   Residual variance: {mdf1.scale:.4f}")
    print(f"   ICC (animal): {mdf1.cov_re.values[0][0] / (mdf1.cov_re.values[0][0] + mdf1.scale):.3f}")

except ImportError:
    print("\nstatsmodels not available")
except Exception as e:
    print(f"\nMixed model error: {e}")


# In[18]:


# Cell 6: Summary figure showing the crossed design
fig, axes = plt.subplots(2, 2, figsize=(16, 12))

# Panel A: Raw data - all animal-reaction combinations
ax = axes[0, 0]
for animal in df_animal_reaction['animal'].unique()[:20]:  # Show first 20 animals
    animal_data = df_animal_reaction[df_animal_reaction['animal'] == animal]
    ax.plot(animal_data['sort_prob'], animal_data['pct_v147_pos'],
           marker='o', alpha=0.3, linewidth=1, markersize=6)

ax.set_xlabel('Sun1 sorting probability', fontweight='bold', fontsize=12)
ax.set_ylabel('% v147+ cells', fontweight='bold', fontsize=12)
ax.set_title('A. Individual animal trajectories\n(first 20 animals shown)',
            fontweight='bold', fontsize=12)
ax.grid(alpha=0.3)

# Panel B: Average by reaction (pooling animals)
ax = axes[0, 1]
reaction_avg = df_animal_reaction.groupby('reaction').agg({
    'sort_prob': 'first',
    'pct_v147_pos': 'mean',
    'n_cells': 'sum'
})
scatter = ax.scatter(reaction_avg['sort_prob'], reaction_avg['pct_v147_pos'],
                    s=reaction_avg['n_cells']/50, alpha=0.7, edgecolors='black')
ax.set_xlabel('Sun1 sorting probability (reaction)', fontweight='bold', fontsize=12)
ax.set_ylabel('Mean % v147+ cells (across animals)', fontweight='bold', fontsize=12)
ax.set_title('B. Reaction-level averages',
            fontweight='bold', fontsize=12)

r_rxn, p_rxn = spearmanr(reaction_avg['sort_prob'], reaction_avg['pct_v147_pos'])
ax.text(0.05, 0.95, f"Spearman ρ = {r_rxn:.3f}\np = {p_rxn:.3e}",
       transform=ax.transAxes, verticalalignment='top',
       bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))

z = np.polyfit(reaction_avg['sort_prob'], reaction_avg['pct_v147_pos'], 1)
p_fit = np.poly1d(z)
x_line = np.linspace(reaction_avg['sort_prob'].min(), reaction_avg['sort_prob'].max(), 100)
ax.plot(x_line, p_fit(x_line), "r--", linewidth=2)
ax.grid(alpha=0.3)

# Panel C: Average by animal (across reactions)
ax = axes[1, 0]
animal_avg = df_animal_reaction.groupby('animal').agg({
    'sort_prob': 'mean',
    'pct_v147_pos': 'mean',
    'n_cells': 'sum'
})
scatter = ax.scatter(animal_avg['sort_prob'], animal_avg['pct_v147_pos'],
                    s=animal_avg['n_cells']/50, alpha=0.7, edgecolors='black')
ax.set_xlabel('Mean Sun1 sorting probability (per animal)', fontweight='bold', fontsize=12)
ax.set_ylabel('Mean % v147+ cells (per animal)', fontweight='bold', fontsize=12)
ax.set_title('C. Animal-level averages',
            fontweight='bold', fontsize=12)

r_animal, p_animal = spearmanr(animal_avg['sort_prob'], animal_avg['pct_v147_pos'])
ax.text(0.05, 0.95, f"Spearman ρ = {r_animal:.3f}\np = {p_animal:.3e}",
       transform=ax.transAxes, verticalalignment='top',
       bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))

z2 = np.polyfit(animal_avg['sort_prob'], animal_avg['pct_v147_pos'], 1)
p_fit2 = np.poly1d(z2)
x_line2 = np.linspace(animal_avg['sort_prob'].min(), animal_avg['sort_prob'].max(), 100)
ax.plot(x_line2, p_fit2(x_line2), "r--", linewidth=2)
ax.grid(alpha=0.3)

# Panel D: Summary text
ax = axes[1, 1]
ax.axis('off')

summary_text = f"""EXPERIMENTAL DESIGN SUMMARY

Structure:
  • {df_cells['animal'].nunique()} animals
  • {df_cells['reaction'].nunique()} reactions  
  • {len(df_animal_reaction)} animal-reaction combinations
  • Animals per reaction: {animals_per_reaction.mean():.1f}
  • Reactions per animal: {reactions_per_animal.mean():.1f}

Analysis Levels:

1. REACTION-LEVEL (pooled animals):
   Spearman ρ = {r_rxn:.3f}, p = {p_rxn:.2e}

2. ANIMAL-LEVEL (across reactions):
   Spearman ρ = {r_animal:.3f}, p = {p_animal:.2e}

3. WITHIN-ANIMAL (preferred):
   Mean correlation: {df_within_animal['r_pct_v147'].mean():.3f}
   Animals tested: {len(df_within_animal)}

4. MIXED-EFFECTS MODEL:
   Fixed effect (sort): β = {mdf1.params['sort_centered']:.2f}%
   p-value: {mdf1.pvalues['sort_centered']:.2e}

CONCLUSION:
The mixed-effects model accounts for:
  ✓ Animal-level clustering (biological variation)
  ✓ Reaction-level effects (technical variation)
  ✓ Crossed design (animals × reactions)
"""

ax.text(0.1, 0.9, summary_text, transform=ax.transAxes,
       fontsize=10, verticalalignment='top', family='monospace',
       bbox=dict(boxstyle='round,pad=1', facecolor='lightyellow', 
                edgecolor='black', linewidth=2))

plt.tight_layout()
plt.savefig('sun1_v147_crossed_design_summary.pdf', dpi=300, bbox_inches='tight')
plt.show()


# In[2]:


adata = sc.read("neurons_bipca_reclustered.h5ad")
adata


# In[3]:


ax = sc.pl.embedding(
    adata,
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


# In[3]:


# Use only highly variable genes
adata_hvg = adata[:, adata.var['highly_deviant']].copy()

sc.tl.dendrogram(
    adata_hvg, 
    groupby="detailed_celltype_annotation"
)

sc.pl.correlation_matrix(
    adata_hvg, 
    "detailed_celltype_annotation",
    save="_dendrogram_neurons.pdf"
)


# In[4]:


import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch
from matplotlib.colors import to_rgb

# Get the dendrogram ordering
dendro_info = adata_hvg.uns['dendrogram_detailed_celltype_annotation']
cell_type_order = dendro_info['categories_ordered']

# Get neurotransmitter for each cell type (taking the most common per cell type)
neurotransmitter_map = {}
for ct in cell_type_order:
    mask = adata.obs['detailed_celltype_annotation'] == ct
    nt_values = adata.obs.loc[mask, 'neurotransmitter'].value_counts()
    neurotransmitter_map[ct] = nt_values.index[0] if len(nt_values) > 0 else 'Unknown'

# Your specific color scheme
cols = ('#FF4A46', '#1CE6FF', '#FFFF00', '#FF34FF')
nt = ("Glut", "GABA", "Chat", "Hist")
nt_color_map = dict(zip(nt, cols))

# Convert hex colors to RGB for imshow
nt_color_map_rgb = {n: to_rgb(c) for n, c in nt_color_map.items()}
nt_color_map_rgb['Unknown'] = to_rgb('#CCCCCC')

# Create annotation array with RGB values
nt_array = np.array([[nt_color_map_rgb.get(neurotransmitter_map[ct], nt_color_map_rgb['Unknown']) 
                      for ct in cell_type_order]])

# Calculate figure height to make squares (1 row, so height = width / n_columns)
n_cols = len(cell_type_order)
fig_width = 18
fig_height = fig_width / n_cols  # This makes each cell square

# Create figure
fig, ax = plt.subplots(figsize=(fig_width, fig_height))

# Plot annotation with equal aspect to maintain squares
ax.imshow(nt_array, aspect='equal', interpolation='nearest')

# Add faint grid
ax.set_xticks(np.arange(-.5, len(cell_type_order), 1), minor=True)
ax.set_yticks(np.arange(-.5, 1, 1), minor=True)
ax.grid(which='minor', color='gray', linestyle='-', linewidth=0.5, alpha=0.3)

# Set labels
ax.set_yticks([0])
ax.set_yticklabels(['Neurotransmitter'], fontsize=12)
ax.set_xticks(range(len(cell_type_order)))
ax.set_xticklabels(cell_type_order, rotation=90, ha='right', fontsize=8)

# Add legend
legend_elements = [Patch(facecolor=nt_color_map[n], label=n) for n in nt]
ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.02, 1), 
          fontsize=10, title='Neurotransmitter')

plt.tight_layout()
plt.savefig('neurotransmitter_annotation_neurons.pdf', bbox_inches='tight', dpi=300)
plt.show()


# In[5]:


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
combined


# In[6]:


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
combined_spatial.to_csv('df_spatial.csv', index=False)
combined_spatial


# In[7]:


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch
from matplotlib.colors import to_rgb

# Get the dendrogram ordering
dendro_info = adata_hvg.uns['dendrogram_detailed_celltype_annotation']
cell_type_order = dendro_info['categories_ordered']

# Assuming your dataframes are named 'combined' and 'combined_spatial'
# Create the flag dataframe
flag_df = pd.DataFrame(index=combined.index, columns=combined.columns)

# Determine flags based on threshold > 0.1
for region in combined.index:
    for celltype in combined.columns:
        snrna_pos = combined.loc[region, celltype] > 0.1
        spatial_pos = combined_spatial.loc[region, celltype] > 0.1 if region in combined_spatial.index else False

        if snrna_pos and spatial_pos:
            flag_df.loc[region, celltype] = 'both'
        elif snrna_pos:
            flag_df.loc[region, celltype] = 'snRNA'
        elif spatial_pos:
            flag_df.loc[region, celltype] = 'spatial'
        else:
            flag_df.loc[region, celltype] = 'negative'

# Reorder columns to match dendrogram order
flag_df = flag_df[cell_type_order]

# Define colors for each category
category_colors = {
    'negative': '#FFFFFF',  # white
    'snRNA': '#CC0744', 
    'spatial': '#A079BF',  
    'both': '#000035'     
}

# Convert to RGB
category_colors_rgb = {cat: to_rgb(col) for cat, col in category_colors.items()}

# Create annotation array with RGB values
annotation_array = np.array([[category_colors_rgb[flag_df.loc[region, ct]] 
                              for ct in cell_type_order] 
                             for region in flag_df.index])

# Create figure with more height
fig, ax = plt.subplots(figsize=(18, 4))

# Plot annotation
ax.imshow(annotation_array, aspect='auto')

# Add faint grid
ax.set_xticks(np.arange(-.5, len(cell_type_order), 1), minor=True)
ax.set_yticks(np.arange(-.5, len(flag_df.index), 1), minor=True)
ax.grid(which='minor', color='gray', linestyle='-', linewidth=0.5, alpha=0.3)

# Set labels with larger font
ax.set_yticks(range(len(flag_df.index)))
ax.set_yticklabels(flag_df.index, fontsize=12)
ax.set_xticks(range(len(cell_type_order)))
ax.set_xticklabels(cell_type_order, rotation=90, ha='right', fontsize=8)

# Add legend
legend_elements = [
    Patch(facecolor=category_colors['snRNA'], label='snRNA only'),
    Patch(facecolor=category_colors['spatial'], label='Spatial only'),
    Patch(facecolor=category_colors['both'], label='Both')
]
ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.02, 1), 
          fontsize=10, title='Detection')

plt.tight_layout()
plt.savefig('region_celltype_annotation_neurons.pdf', bbox_inches='tight', dpi=300)
plt.show()


# In[8]:


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

# Get the dendrogram ordering
dendro_info     = adata_hvg.uns['dendrogram_detailed_celltype_annotation']
cell_type_order = dendro_info['categories_ordered']

# Build the combined score dataframe (0-1 scale)
# PVH has only snRNA data → use raw ESµ directly
# All other regions have both datasets → average the two ESµ values
score_df = pd.DataFrame(index=combined.index, columns=cell_type_order, dtype=float)

single_dataset_regions = {'PVH', "DVC"}  # only snRNA available, do not divide by 2

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
plt.savefig('region_celltype_annotation_combined_score.pdf', bbox_inches='tight', dpi=300)
plt.show()


# In[11]:


df_cellect = pd.read_csv('/projects/perslab/people/qbx224/CELLECT/CELLECT-ALLPLEX/CELLECT-LDSC/results/prioritization.csv')
df_cellect["mlog10p"] = -np.log10(df_cellect["pvalue"])
df_cellect["z"] = df_cellect["beta"] / df_cellect["beta_se"]
df_cellect.to_csv('df_cellect.csv', index=False)
df_cellect


# In[12]:


# Get the dendrogram ordering
dendro_info = adata_hvg.uns['dendrogram_detailed_celltype_annotation']
cell_type_order = dendro_info['categories_ordered']

# Pivot the dataframe to get GWAS x cell_type matrix
pivot_df = df_cellect.pivot(index='gwas', columns='annotation', values='mlog10p')

# Also pivot p-values for significance testing
pivot_pval = df_cellect.pivot(index='gwas', columns='annotation', values='pvalue')

# Reorder columns to match dendrogram order
common_celltypes = [ct for ct in cell_type_order if ct in pivot_df.columns]
pivot_df = pivot_df[common_celltypes]
pivot_pval = pivot_pval[common_celltypes]

# Bonferroni correction threshold
n_tests = 41
bonferroni_threshold = 0.05 / n_tests

# Calculate figure size to maintain square cells
n_rows = len(pivot_df.index)
n_cols = len(common_celltypes)
fig_width = 18
fig_height = fig_width * (n_rows / n_cols)  # Scale height based on number of rows

# Create figure
fig, ax = plt.subplots(figsize=(fig_width, fig_height))

# Plot heatmap with mako_r colormap
im = ax.imshow(pivot_df.values, aspect='equal', cmap='mako_r', interpolation='nearest')

# Add asterisks for significant p-values
for i in range(len(pivot_df.index)):
    for j in range(len(common_celltypes)):
        if pivot_pval.iloc[i, j] < bonferroni_threshold:
            ax.text(j, i, '*', ha='center', va='center', 
                   color='white', fontsize=16, fontweight='bold',
                   clip_on=False)

# Add faint grid
ax.set_xticks(np.arange(-.5, len(common_celltypes), 1), minor=True)
ax.set_yticks(np.arange(-.5, len(pivot_df.index), 1), minor=True)
ax.grid(which='minor', color='gray', linestyle='-', linewidth=0.5, alpha=0.3)

# Set labels
ax.set_yticks(range(len(pivot_df.index)))
ax.set_yticklabels(pivot_df.index, fontsize=10)
ax.set_xticks(range(len(common_celltypes)))
ax.set_xticklabels(common_celltypes, rotation=90, ha='right', fontsize=8)

# Add colorbar
cbar = plt.colorbar(im, ax=ax, label='-log10(p-value)', fraction=0.046, pad=0.04)

plt.tight_layout()
plt.savefig('gwas_celltype_heatmap.pdf', bbox_inches='tight', dpi=300)
plt.show()


# In[13]:


df_cellect = pd.read_csv('df_cellect.csv')
df_cellect = df_cellect[df_cellect['gwas'] == 'BMI_Loh2018_UKBB_munged_txt'].copy()
df_cellect['gwas'] = 'BMI'  # clean y-axis label
# Get the dendrogram ordering
dendro_info = adata_hvg.uns['dendrogram_detailed_celltype_annotation']
cell_type_order = dendro_info['categories_ordered']
# Filter to single GWAS trait

# Pivot the dataframe to get GWAS x cell_type matrix
pivot_df = df_cellect.pivot(index='gwas', columns='annotation', values='mlog10p')

# Also pivot p-values for significance testing
pivot_pval = df_cellect.pivot(index='gwas', columns='annotation', values='pvalue')

# Reorder columns to match dendrogram order
common_celltypes = [ct for ct in cell_type_order if ct in pivot_df.columns]
pivot_df = pivot_df[common_celltypes]
pivot_pval = pivot_pval[common_celltypes]

# Bonferroni correction threshold
n_tests = 41
bonferroni_threshold = 0.05 / n_tests

# Calculate figure size to maintain square cells
n_rows = len(pivot_df.index)
n_cols = len(common_celltypes)
fig_width = 18
fig_height = fig_width * (n_rows / n_cols)  # Scale height based on number of rows

# Create figure
fig, ax = plt.subplots(figsize=(fig_width, fig_height))

# Plot heatmap with mako_r colormap
im = ax.imshow(pivot_df.values, aspect='equal', cmap='mako_r', interpolation='nearest')

# Add asterisks for significant p-values
for i in range(len(pivot_df.index)):
    for j in range(len(common_celltypes)):
        if pivot_pval.iloc[i, j] < bonferroni_threshold:
            ax.text(j, i, '*', ha='center', va='center', 
                   color='white', fontsize=16, fontweight='bold',
                   clip_on=False)

# Add faint grid
ax.set_xticks(np.arange(-.5, len(common_celltypes), 1), minor=True)
ax.set_yticks(np.arange(-.5, len(pivot_df.index), 1), minor=True)
ax.grid(which='minor', color='gray', linestyle='-', linewidth=0.5, alpha=0.3)

# Set labels
ax.set_yticks(range(len(pivot_df.index)))
ax.set_yticklabels(pivot_df.index, fontsize=10)
ax.set_xticks(range(len(common_celltypes)))
ax.set_xticklabels(common_celltypes, rotation=90, ha='right', fontsize=8)

# Add colorbar
cbar = plt.colorbar(im, ax=ax, label='-log10(p-value)', fraction=0.046, pad=0.04)

plt.tight_layout()
plt.savefig('gwas_celltype_heatmap_BMI_neurons.pdf', bbox_inches='tight', dpi=300)
plt.show()


# In[14]:


# Cell 1: Load and explore the data
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.stats import fisher_exact, chi2_contingency
import scanpy as sc

# Check the data
print("Available regions:")
print(adata.obs['region_collapsed'].value_counts())
print("\nScore distribution:")
print(adata.obs['score_region_collapsed'].describe())
print(f"\nNumber of cell types: {adata.obs['detailed_celltype_annotation'].nunique()}")

# Visualize score distribution
fig, axes = plt.subplots(1, 2, figsize=(12, 4))
axes[0].hist(adata.obs['score_region_collapsed'], bins=50, edgecolor='black')
axes[0].axvline(0.4, color='red', linestyle='--', label='Current cutoff (0.4)')
axes[0].set_xlabel('Region assignment score')
axes[0].set_ylabel('Frequency')
axes[0].set_title('Distribution of region assignment scores')
axes[0].legend()

# Score by region
adata.obs.boxplot(column='score_region_collapsed', by='region_collapsed', ax=axes[1], rot=45)
axes[1].set_title('Score distribution by region')
axes[1].set_xlabel('Region')
plt.suptitle('')
plt.tight_layout()
plt.show()


# In[15]:


# Cell 2: Benchmark different score cutoffs
def compute_enrichment_at_cutoff(adata, score_cutoff, exclude_regions=['undefined', 'NA']):
    """
    Compute cell type enrichment in regions at a given score cutoff.
    Returns contingency tables and enrichment statistics.
    """
    # Filter high-confidence assignments
    mask_highconf = adata.obs['score_region_collapsed'] > score_cutoff
    mask_valid_region = ~adata.obs['region_collapsed'].isin(exclude_regions)
    mask = mask_highconf & mask_valid_region

    # Get filtered data
    df = adata.obs[mask][['detailed_celltype_annotation', 'region_collapsed']].copy()

    # Compute contingency table
    contingency = pd.crosstab(
        df['detailed_celltype_annotation'],
        df['region_collapsed']
    )

    # Compute proportions (normalize by region size)
    proportions = contingency.div(contingency.sum(axis=0), axis=1)

    # Compute expected proportions (if no enrichment)
    expected_prop = contingency.sum(axis=1) / contingency.sum().sum()

    # Compute fold-enrichment
    fold_enrichment = proportions.div(expected_prop, axis=0)

    # Compute statistics for each cell type - region pair
    results = []
    for celltype in contingency.index:
        for region in contingency.columns:
            # 2x2 contingency table for Fisher's exact test
            a = contingency.loc[celltype, region]  # celltype in region
            b = contingency.loc[celltype, :].sum() - a  # celltype in other regions
            c = contingency.loc[:, region].sum() - a  # other celltypes in region
            d = contingency.sum().sum() - a - b - c  # other celltypes in other regions

            if a + b > 0:  # Only test if celltype exists
                oddsratio, pval = fisher_exact([[a, b], [c, d]], alternative='greater')

                results.append({
                    'celltype': celltype,
                    'region': region,
                    'n_observed': a,
                    'n_celltype_total': a + b,
                    'n_region_total': a + c,
                    'proportion_in_region': a / (a + c) if (a + c) > 0 else 0,
                    'proportion_expected': (a + b) / contingency.sum().sum(),
                    'fold_enrichment': fold_enrichment.loc[celltype, region],
                    'odds_ratio': oddsratio,
                    'pvalue': pval,
                    'score_cutoff': score_cutoff
                })

    df_results = pd.DataFrame(results)

    # FDR correction
    from statsmodels.stats.multitest import multipletests
    if len(df_results) > 0:
        _, df_results['pvalue_adj'], _, _ = multipletests(
            df_results['pvalue'], 
            method='fdr_bh'
        )

    return df_results, contingency, proportions

# Test multiple cutoffs
cutoffs = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
cutoff_results = {}

for cutoff in cutoffs:
    df_enrich, contingency, proportions = compute_enrichment_at_cutoff(adata, cutoff)
    cutoff_results[cutoff] = {
        'enrichment': df_enrich,
        'contingency': contingency,
        'proportions': proportions,
        'n_cells': len(adata.obs[adata.obs['score_region_collapsed'] > cutoff]),
        'n_significant': (df_enrich['pvalue_adj'] < 0.05).sum()
    }
    print(f"Cutoff {cutoff}: {cutoff_results[cutoff]['n_cells']} cells, "
          f"{cutoff_results[cutoff]['n_significant']} significant enrichments")


# In[16]:


# Cell 3: Visualize cutoff benchmark
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Panel A: Number of cells retained
ax = axes[0, 0]
n_cells = [cutoff_results[c]['n_cells'] for c in cutoffs]
ax.plot(cutoffs, n_cells, marker='o', linewidth=2, markersize=8)
ax.axvline(0.4, color='red', linestyle='--', alpha=0.5, label='Current (0.4)')
ax.set_xlabel('Score cutoff', fontweight='bold')
ax.set_ylabel('Number of cells retained', fontweight='bold')
ax.set_title('A. Cell retention vs cutoff', fontweight='bold')
ax.grid(alpha=0.3)
ax.legend()

# Panel B: Number of significant enrichments
ax = axes[0, 1]
n_sig = [cutoff_results[c]['n_significant'] for c in cutoffs]
ax.plot(cutoffs, n_sig, marker='o', linewidth=2, markersize=8, color='green')
ax.axvline(0.4, color='red', linestyle='--', alpha=0.5, label='Current (0.4)')
ax.set_xlabel('Score cutoff', fontweight='bold')
ax.set_ylabel('Number of significant enrichments (FDR<0.05)', fontweight='bold')
ax.set_title('B. Statistical power vs cutoff', fontweight='bold')
ax.grid(alpha=0.3)
ax.legend()

# Panel C: Mean fold-enrichment of significant hits
ax = axes[1, 0]
mean_fold = []
for c in cutoffs:
    sig_hits = cutoff_results[c]['enrichment'][cutoff_results[c]['enrichment']['pvalue_adj'] < 0.05]
    if len(sig_hits) > 0:
        mean_fold.append(sig_hits['fold_enrichment'].mean())
    else:
        mean_fold.append(0)
ax.plot(cutoffs, mean_fold, marker='o', linewidth=2, markersize=8, color='purple')
ax.axvline(0.4, color='red', linestyle='--', alpha=0.5, label='Current (0.4)')
ax.set_xlabel('Score cutoff', fontweight='bold')
ax.set_ylabel('Mean fold-enrichment', fontweight='bold')
ax.set_title('C. Enrichment strength vs cutoff', fontweight='bold')
ax.grid(alpha=0.3)
ax.legend()

# Panel D: Tradeoff plot
ax = axes[1, 1]
ax.scatter(n_cells, n_sig, s=100, c=cutoffs, cmap='viridis', edgecolors='black', linewidth=1.5)
for i, c in enumerate(cutoffs):
    ax.annotate(f'{c}', (n_cells[i], n_sig[i]), 
               xytext=(5, 5), textcoords='offset points', fontsize=9)
# Mark 0.4
idx_04 = cutoffs.index(0.4)
ax.scatter(n_cells[idx_04], n_sig[idx_04], s=200, marker='*', 
          color='red', edgecolors='darkred', linewidth=2, zorder=10)
ax.set_xlabel('Number of cells retained', fontweight='bold')
ax.set_ylabel('Number of significant enrichments', fontweight='bold')
ax.set_title('D. Power-Sensitivity tradeoff', fontweight='bold')
ax.grid(alpha=0.3)
cbar = plt.colorbar(ax.collections[0], ax=ax)
cbar.set_label('Score cutoff', rotation=270, labelpad=15)

plt.tight_layout()
plt.savefig('cutoff_benchmark.pdf', dpi=300, bbox_inches='tight')
plt.show()

print("\nRecommended cutoff based on balance:")
print(f"Current (0.4): {cutoff_results[0.4]['n_cells']} cells, {cutoff_results[0.4]['n_significant']} enrichments")


# In[17]:


# Cell 4: Use optimal cutoff for detailed analysis
optimal_cutoff = 0.4  # Adjust based on Cell 3 results

df_enrichment, contingency, proportions = compute_enrichment_at_cutoff(adata, optimal_cutoff)

# Filter significant results
df_sig = df_enrichment[df_enrichment['pvalue_adj'] < 0.05].copy()
df_sig = df_sig.sort_values('fold_enrichment', ascending=False)

print(f"Significant enrichments (FDR < 0.05): {len(df_sig)}")
print("\nTop 20 enrichments:")
print(df_sig[['celltype', 'region', 'n_observed', 'fold_enrichment', 'pvalue_adj']].head(20))


# In[18]:


# Cell 5: Heatmap of fold-enrichment (rotated 90°, regions on Y-axis)
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Get the dendrogram ordering from the neurotransmitter plot
dendro_info = adata_hvg.uns['dendrogram_detailed_celltype_annotation']
cell_type_order = dendro_info['categories_ordered']
dendrogram_order = list(cell_type_order)

# Pivot to wide format for heatmap (TRANSPOSED: regions as rows, celltypes as columns)
pivot_fold = df_enrichment.pivot(index='region', columns='celltype', values='fold_enrichment')
pivot_pval = df_enrichment.pivot(index='region', columns='celltype', values='pvalue_adj')

# Reorder columns to match dendrogram order
common_celltypes = [ct for ct in dendrogram_order if ct in pivot_fold.columns]
pivot_fold = pivot_fold[common_celltypes]
pivot_pval = pivot_pval[common_celltypes]

# Create significance mask (FDR < 0.05)
sig_mask = pivot_pval < 0.05

# Calculate figure size to maintain square cells (now transposed)
n_rows = len(pivot_fold.index)  # regions
n_cols = len(pivot_fold.columns)  # celltypes
fig_height = 8  # Base height for regions
fig_width = fig_height * (n_cols / n_rows)  # Scale width for square cells

# Create figure
fig, ax = plt.subplots(figsize=(fig_width, fig_height))

# Plot heatmap with colormap centered at 1
im = ax.imshow(pivot_fold.values, aspect='equal', cmap='RdPu', 
               vmin=0, vmax=5, interpolation='nearest')

# Add asterisks for significant enrichments
for i in range(len(pivot_fold.index)):
    for j in range(len(pivot_fold.columns)):
        if sig_mask.iloc[i, j]:
            ax.text(j, i, '*', ha='center', va='center', 
                   color='white', fontsize=10, fontweight='bold',
                   clip_on=False)

# Add faint grid
ax.set_xticks(np.arange(-.5, len(pivot_fold.columns), 1), minor=True)
ax.set_yticks(np.arange(-.5, len(pivot_fold.index), 1), minor=True)
ax.grid(which='minor', color='gray', linestyle='-', linewidth=0.5, alpha=0.3)

# Set labels (X-axis: cell types, Y-axis: regions)
ax.set_xticks(range(len(pivot_fold.columns)))
ax.set_xticklabels(pivot_fold.columns, rotation=90, ha='right', fontsize=7)
ax.set_yticks(range(len(pivot_fold.index)))
ax.set_yticklabels(pivot_fold.index, fontsize=9)

# Labels
ax.set_xlabel('Cell Type (dendrogram ordered)', fontweight='bold', fontsize=11)
ax.set_ylabel('Region', fontweight='bold', fontsize=11)
ax.set_title(f'Cell Type Enrichment by Region (score > {optimal_cutoff})\n* = FDR < 0.05', 
            fontweight='bold', fontsize=13, pad=10)

# Add colorbar
cbar = plt.colorbar(im, ax=ax, label='Fold enrichment', fraction=0.046, pad=0.04)

plt.tight_layout()
plt.savefig('enrichment_heatmap_all_neurons.pdf', dpi=300, bbox_inches='tight')
plt.show()


# In[19]:


import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl
import numpy as np
import pandas as pd
from matplotlib.patches import Patch
from matplotlib.colors import to_rgb, Normalize
from matplotlib.cm import ScalarMappable

# ── 0. Global font sizes ──────────────────────────────────────────────────────
mpl.rcParams.update({
    'font.size':             15,
    'ytick.labelsize':       15,
    'xtick.labelsize':       12,
    'legend.fontsize':       15,
    'legend.title_fontsize': 15,
})

# ── 1. Shared setup ───────────────────────────────────────────────────────────
dendro_info     = adata_hvg.uns['dendrogram_detailed_celltype_annotation']
cell_type_order = dendro_info['categories_ordered']
n_cols          = len(cell_type_order)

# ── 2. Prepare all data up front ──────────────────────────────────────────────

# --- Neurotransmitter ---
neurotransmitter_map = {}
for ct in cell_type_order:
    mask    = adata.obs['detailed_celltype_annotation'] == ct
    nt_vals = adata.obs.loc[mask, 'neurotransmitter'].value_counts()
    neurotransmitter_map[ct] = nt_vals.index[0] if len(nt_vals) > 0 else 'Unknown'

nt_cols      = ('#FF4A46', '#1CE6FF', '#FFFF00', '#FF34FF')
nt_labels    = ('Glut', 'GABA', 'Chat', 'Hist')
nt_color_map = dict(zip(nt_labels, nt_cols))
nt_color_rgb = {n: to_rgb(c) for n, c in nt_color_map.items()}
nt_color_rgb['Unknown'] = to_rgb('#CCCCCC')
nt_array = np.array([[nt_color_rgb.get(neurotransmitter_map[ct], nt_color_rgb['Unknown'])
                       for ct in cell_type_order]])

# --- Combined projection score ---
score_df = pd.DataFrame(index=combined.index, columns=cell_type_order, dtype=float)
for region in combined.index:
    for celltype in cell_type_order:
        snrna_val   = combined.loc[region, celltype] \
                      if celltype in combined.columns else np.nan
        spatial_val = combined_spatial.loc[region, celltype] \
                      if (region in combined_spatial.index
                          and celltype in combined_spatial.columns) else np.nan
        if np.isnan(snrna_val) and np.isnan(spatial_val):
            score_df.loc[region, celltype] = np.nan
        else:
            score_df.loc[region, celltype] = np.nansum([snrna_val, spatial_val])
score_array = score_df.values.astype(float)
score_norm  = Normalize(vmin=0, vmax=2)
score_cmap  = plt.get_cmap('rocket_r')
score_cmap.set_bad(color='lightgray')

# --- Enrichment ---
pivot_fold  = df_enrichment.pivot(index='region', columns='celltype', values='fold_enrichment')
pivot_epval = df_enrichment.pivot(index='region', columns='celltype', values='pvalue_adj')
common_enrich = [ct for ct in cell_type_order if ct in pivot_fold.columns]
pivot_fold  = pivot_fold[common_enrich]
pivot_epval = pivot_epval[common_enrich]
sig_mask    = pivot_epval < 0.05

# --- GWAS BMI ---
df_cellect = pd.read_csv('df_cellect.csv')
df_cellect = df_cellect[df_cellect['gwas'] == 'BMI_Loh2018_UKBB_munged_txt'].copy()
df_cellect['gwas'] = 'BMI'
pivot_gwas  = df_cellect.pivot(index='gwas', columns='annotation', values='mlog10p')
pivot_gpval = df_cellect.pivot(index='gwas', columns='annotation', values='pvalue')
common_gwas = [ct for ct in cell_type_order if ct in pivot_gwas.columns]
pivot_gwas  = pivot_gwas[common_gwas]
pivot_gpval = pivot_gpval[common_gwas]
bonferroni  = 0.05 / 41

# ── 3. Figure layout ──────────────────────────────────────────────────────────
fig_width  = 20
cell_size  = fig_width / n_cols

n_rows_nt     = 1
n_rows_score  = len(score_df.index)
n_rows_enrich = len(pivot_fold.index)
n_rows_gwas   = len(pivot_gwas.index)

dendro_h = 1.2
height_ratios = [
    dendro_h,
    n_rows_nt     * cell_size,
    n_rows_score  * cell_size,
    n_rows_enrich * cell_size,
    n_rows_gwas   * cell_size,
]
fig_height = sum(height_ratios)

fig = plt.figure(figsize=(fig_width, fig_height))
gs  = gridspec.GridSpec(
    5, 1,
    height_ratios=height_ratios,
    hspace=0.04,
    left=0.08, right=0.82,
    top=0.99, bottom=0.04
)

ax_dendro = fig.add_subplot(gs[0])
ax_nt     = fig.add_subplot(gs[1])
ax_score  = fig.add_subplot(gs[2])
ax_enrich = fig.add_subplot(gs[3])
ax_gwas   = fig.add_subplot(gs[4])

xlim = (-0.5, n_cols - 0.5)

# ── 4. Dendrogram ─────────────────────────────────────────────────────────────
dendro_data   = dendro_info['dendrogram_info']
icoord        = np.array(dendro_data['icoord'])
dcoord        = np.array(dendro_data['dcoord'])
icoord_scaled = (icoord - 5) / 10

for xs, ys in zip(icoord_scaled, dcoord):
    ax_dendro.plot(xs, ys, 'k-', linewidth=0.8)

ax_dendro.set_xlim(xlim)
ax_dendro.set_ylim(bottom=0)
ax_dendro.axis('off')

# ── 5. Neurotransmitter ───────────────────────────────────────────────────────
ax_nt.imshow(nt_array, aspect='auto', interpolation='nearest')
ax_nt.set_xlim(xlim)
ax_nt.set_yticks([0])
ax_nt.set_yticklabels(['Neurotransmitter'], fontsize=13)
ax_nt.set_xticks([])
ax_nt.tick_params(left=False)
ax_nt.set_xticks(np.arange(-.5, n_cols, 1), minor=True)
ax_nt.set_yticks(np.arange(-.5, 1, 1), minor=True)
ax_nt.grid(which='minor', color='gray', linestyle='-', linewidth=0.5, alpha=0.3)

legend_nt = [Patch(facecolor=nt_color_map[n], label=n) for n in nt_labels]
ax_nt.legend(handles=legend_nt, loc='lower left',
             bbox_to_anchor=(1.01, 0.8),
             fontsize=11, title='Neurotransmitter',
             title_fontsize=11, frameon=True,
             bbox_transform=ax_nt.transAxes)

# ── 6. Combined projection score ──────────────────────────────────────────────
single_dataset_regions = {'PVH', "DVC"}

score_df = pd.DataFrame(index=combined.index, columns=cell_type_order, dtype=float)
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
            score_df.loc[region, celltype] = np.nansum([snrna_val, spatial_val])
        else:
            score_df.loc[region, celltype] = np.nansum([snrna_val, spatial_val]) / 2

score_array  = score_df.values.astype(float)
score_norm   = Normalize(vmin=0, vmax=1)
score_cmap   = plt.get_cmap('rocket_r')
score_cmap.set_bad(color='lightgray')
n_rows_score = len(score_df.index)

im_score = ax_score.imshow(np.ma.masked_invalid(score_array),
                            aspect='auto', cmap=score_cmap, norm=score_norm,
                            interpolation='nearest')
ax_score.set_xlim(xlim)
ax_score.set_yticks(range(n_rows_score))
ax_score.set_yticklabels(score_df.index, fontsize=13)
ax_score.set_xticks([])
ax_score.set_xticks(np.arange(-.5, n_cols, 1), minor=True)
ax_score.set_yticks(np.arange(-.5, n_rows_score, 1), minor=True)
ax_score.grid(which='minor', color='gray', linestyle='-', linewidth=0.5, alpha=0.3)

# Asterisks: uniform threshold > 0.36 across all regions
for row_i, region in enumerate(score_df.index):
    for col_i in range(n_cols):
        val = score_array[row_i, col_i]
        if np.isnan(val):
            continue
        if val > 0.36:
            brightness = score_norm(val)
            text_color = 'white' if brightness > 0.5 else 'black'
            ax_score.text(col_i, row_i, '*', ha='center', va='center',
                          fontsize=15, color=text_color, fontweight='bold')

# ── 7. Enrichment heatmap ─────────────────────────────────────────────────────
col_idx_map = {ct: i for i, ct in enumerate(cell_type_order)}

enrich_full = np.full((n_rows_enrich, n_cols), np.nan)
for j, ct in enumerate(common_enrich):
    if ct in col_idx_map:
        enrich_full[:, col_idx_map[ct]] = pivot_fold.values[:, j]

enrich_cmap = plt.get_cmap('RdPu')
enrich_cmap.set_bad(color='lightgray')
im_enrich = ax_enrich.imshow(np.ma.masked_invalid(enrich_full),
                              aspect='auto', cmap=enrich_cmap,
                              vmin=0, vmax=5, interpolation='nearest')
ax_enrich.set_xlim(xlim)
ax_enrich.set_yticks(range(n_rows_enrich))
ax_enrich.set_yticklabels(pivot_fold.index, fontsize=13)
ax_enrich.set_xticks([])
ax_enrich.set_xticks(np.arange(-.5, n_cols, 1), minor=True)
ax_enrich.set_yticks(np.arange(-.5, n_rows_enrich, 1), minor=True)
ax_enrich.grid(which='minor', color='gray', linestyle='-', linewidth=0.5, alpha=0.3)

for i in range(n_rows_enrich):
    for j, ct in enumerate(common_enrich):
        if sig_mask.iloc[i, j] and ct in col_idx_map:
            col_pos    = col_idx_map[ct]
            brightness = enrich_full[i, col_pos] / 5
            ax_enrich.text(col_pos, i, '*', ha='center', va='center',
                           color='white' if brightness > 0.5 else 'black',
                           fontsize=15, fontweight='bold')

# ── 8. GWAS BMI ───────────────────────────────────────────────────────────────
gwas_full = np.full((n_rows_gwas, n_cols), np.nan)
for j, ct in enumerate(common_gwas):
    if ct in col_idx_map:
        gwas_full[:, col_idx_map[ct]] = pivot_gwas.values[:, j]

gwas_cmap = plt.get_cmap('mako_r')
gwas_cmap.set_bad(color='lightgray')
im_gwas = ax_gwas.imshow(np.ma.masked_invalid(gwas_full),
                          aspect='auto', cmap=gwas_cmap, interpolation='nearest')
ax_gwas.set_xlim(xlim)
ax_gwas.set_yticks(range(n_rows_gwas))
ax_gwas.set_yticklabels(pivot_gwas.index, fontsize=13)
ax_gwas.set_xticks(range(n_cols))
ax_gwas.set_xticklabels(cell_type_order, rotation=90, ha='right', fontsize=11)
ax_gwas.set_xticks(np.arange(-.5, n_cols, 1), minor=True)
ax_gwas.set_yticks(np.arange(-.5, n_rows_gwas, 1), minor=True)
ax_gwas.grid(which='minor', color='gray', linestyle='-', linewidth=0.5, alpha=0.3)

for i in range(n_rows_gwas):
    for j, ct in enumerate(common_gwas):
        if pivot_gpval.iloc[i, j] < bonferroni and ct in col_idx_map:
            ax_gwas.text(col_idx_map[ct], i, '*', ha='center', va='center',
                         color='white', fontsize=15, fontweight='bold')

# ── 9. Colorbars via fig.add_axes (fixes PDF export rendering bug) ────────────
fig.canvas.draw()  # finalise layout so .get_position() is accurate

cbar_padding = 0.005
cbar_width   = 0.008

# Combined projection score
pos = ax_score.get_position()
cax_score  = fig.add_axes([pos.x1 + cbar_padding, pos.y0 + pos.height * 0.1,
                            cbar_width, pos.height * 0.8])
cbar_score = fig.colorbar(ScalarMappable(norm=score_norm, cmap=score_cmap), cax=cax_score)
cbar_score.set_ticks([0, 0.25, 0.5, 0.75, 1.0])
cbar_score.set_label('Combined score\n(* mean ESµ > 0.36)', fontsize=11)
cbar_score.ax.tick_params(labelsize=10)

# Enrichment
pos = ax_enrich.get_position()
cax_enrich  = fig.add_axes([pos.x1 + cbar_padding, pos.y0 + pos.height * 0.1,
                             cbar_width, pos.height * 0.8])
cbar_enrich = fig.colorbar(im_enrich, cax=cax_enrich)
cbar_enrich.set_ticks([0, 1, 2, 3, 4, 5])
cbar_enrich.set_label('Fold enrichment\n(* FDR < 0.05)', fontsize=11)
cbar_enrich.ax.tick_params(labelsize=10)

# GWAS BMI
pos = ax_gwas.get_position()
cax_gwas  = fig.add_axes([pos.x1 + cbar_padding, pos.y0 + pos.height * 0.1,
                           cbar_width, pos.height * 0.8])
cbar_gwas = fig.colorbar(im_gwas, cax=cax_gwas)
cbar_gwas.set_ticks([1, 3])
cbar_gwas.set_label('-log₁₀(p)\n(* Bonferroni < 0.05)', fontsize=11)
cbar_gwas.ax.tick_params(labelsize=10)

# ── 10. Save ──────────────────────────────────────────────────────────────────
plt.savefig('combined_panel_figure.pdf', bbox_inches='tight', dpi=300)
plt.show()


# In[33]:


###Here we read in LOO v147 ESµ
indir = "LOO_cellex_results"
records = []

# Loop over each area folder
for area in os.listdir(indir):
    area_dir = os.path.join(indir, area)
    if not os.path.isdir(area_dir):
        continue

    print(f"Scanning area directory: {area}")

    # Loop over all CELLEX result files inside this area folder
    for fn in os.listdir(area_dir):
        if not fn.startswith("CELLEX_ESmu_"):
            continue
        if not fn.endswith(".csv"):
            continue

        path = os.path.join(area_dir, fn)

        # Example filename:
        # "CELLEX_ESmu_LHA_LOO_animal_124.csv"
        # Extract left-out animal ID
        try:
            left_out_animal = fn.split("_animal_")[1].replace(".csv", "")
        except:
            print(f"Could not extract animal ID from {fn}, skipping")
            continue

        # Load CELLEX ESmu table
        df = pd.read_csv(path, index_col=0)

        if "v147" not in df.index:
            print(f"Skipping {fn}: v147 not found")
            continue

        s = df.loc["v147"]  # ESmu values per cluster

        # Build long format records
        for cluster, esmu in s.items():
            records.append({
                "area": area,
                "left_out_animal": left_out_animal,
                "cluster": cluster,
                "esmu": esmu
            })

# Combine into one DataFrame
combined_LOO = pd.DataFrame(records)
combined_LOO.to_csv('combined_LOO.csv', index=False)
combined_LOO


# In[30]:


import pandas as pd
import numpy as np
from statsmodels.regression.mixed_linear_model import MixedLM
from tqdm import tqdm

# Prepare the results list
results_list = []

# Get unique GWAS traits and areas
gwas_traits = df_cellect['gwas'].unique()
areas = combined_LOO['area'].unique()

# For each GWAS and each area separately
for gwas in tqdm(gwas_traits, desc="Processing GWAS"):
    # Get GWAS data for this trait
    gwas_data = df_cellect[df_cellect['gwas'] == gwas][['annotation', 'z']].copy()
    gwas_data.columns = ['cluster', 'z_score']

    for area in areas:
        # Filter combined_LOO for this specific area
        area_data = combined_LOO[combined_LOO['area'] == area].copy()

        # Merge with GWAS z-scores
        merged_data = area_data.merge(gwas_data, on='cluster', how='inner')

        if len(merged_data) == 0:
            continue

        # Add intercept
        merged_data['intercept'] = 1.0

        # Run mixed model: esmu ~ z_score + (1|animal)
        try:
            model = MixedLM(
                endog=merged_data['esmu'],
                exog=merged_data[['intercept', 'z_score']],
                groups=merged_data['left_out_animal']
            )
            result = model.fit(method='lbfgs', maxiter=1000)

            # Extract results
            results_list.append({
                'gwas': gwas,
                'area': area,
                'beta_intercept': result.params['intercept'],
                'beta_z': result.params['z_score'],
                'se_intercept': result.bse['intercept'],
                'se_z': result.bse['z_score'],
                't_stat': result.tvalues['z_score'],
                'pvalue': result.pvalues['z_score'],
                'ci_lower': result.conf_int().loc['z_score', 0],
                'ci_upper': result.conf_int().loc['z_score', 1],
                'converged': result.converged,
                'n_obs': len(merged_data),
                'n_animals': merged_data['left_out_animal'].nunique(),
                'n_celltypes': merged_data['cluster'].nunique()
            })
        except Exception as e:
            print(f"Error for {gwas} - {area}: {str(e)}")
            results_list.append({
                'gwas': gwas,
                'area': area,
                'beta_z': np.nan,
                'se_z': np.nan,
                'pvalue': np.nan,
                'converged': False,
                'error': str(e)
            })

# Create DataFrame
df_mixedmodel_area = pd.DataFrame(results_list)
df_mixedmodel_area['mlog10p'] = -np.log10(df_mixedmodel_area['pvalue'])

# Add significance stars for Bonferroni correction (7 tests per GWAS)
bonferroni_threshold = 0.05 / 7
df_mixedmodel_area['significant'] = df_mixedmodel_area['pvalue'] < bonferroni_threshold

print(f"\nTotal models run: {len(df_mixedmodel_area)}")
print(f"Converged models: {df_mixedmodel_area['converged'].sum()}")
print(f"Significant results (Bonferroni): {df_mixedmodel_area['significant'].sum()}")

# Save results
df_mixedmodel_area.to_csv('mixedmodel_results_by_area.csv', index=False)
print("\nResults saved to mixedmodel_results_by_area.csv")

# Display example for BMI
bmi_results = df_mixedmodel_area[df_mixedmodel_area['gwas'].str.contains('BMI', case=False)]
print("\n=== BMI Results by Area ===")
print(bmi_results[['area', 'beta_z', 'se_z', 'pvalue', 'significant']].sort_values('area'))


# In[31]:


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Load the results
df_area = pd.read_csv('mixedmodel_results_by_area.csv')

# Define area order to match your example
area_order = ['BST', 'CEA', 'LHA', 'PAG', 'PBN', 'PVH', 'PVT']

# Define colors for each area
area_colors = {
    'BST': '#1f77b4',
    'CEA': '#ff7f0e', 
    'LHA': '#2ca02c',
    'PAG': '#9467bd',
    'PBN': '#8c564b',
    'PVH': '#e377c2',
    'PVT': '#17becf'
}

def plot_gwas_forest(gwas_name, df, save_path=None):
    """
    Create forest plot for a single GWAS showing all 7 areas
    """
    # Filter for this GWAS
    gwas_data = df[df['gwas'] == gwas_name].copy()

    # Reorder by area
    gwas_data['area'] = pd.Categorical(gwas_data['area'], categories=area_order, ordered=True)
    gwas_data = gwas_data.sort_values('area')

    # Create figure
    fig, ax = plt.subplots(figsize=(8, 6))

    # Plot each area
    for idx, row in gwas_data.iterrows():
        area = row['area']
        beta = row['beta_z']
        ci_lower = row['ci_lower']
        ci_upper = row['ci_upper']
        significant = row['significant']

        # Get position
        y_pos = area_order.index(area)

        # Plot confidence interval
        ax.plot([ci_lower, ci_upper], [y_pos, y_pos], 
               color=area_colors[area], linewidth=2, alpha=0.8)

        # Plot point estimate
        marker = '*' if significant else 'o'
        markersize = 12 if significant else 10
        ax.scatter(beta, y_pos, color=area_colors[area], 
                  marker=marker, s=markersize**2, zorder=3,
                  edgecolors='black' if significant else 'none', linewidth=1.5)

    # Add reference line at 0
    ax.axvline(0, color='gray', linestyle='--', linewidth=1, alpha=0.5)

    # Formatting
    ax.set_yticks(range(len(area_order)))
    ax.set_yticklabels(area_order)
    ax.set_xlabel('Fixed-effect slope (coeff) ± 95% CI', fontsize=11)
    ax.set_ylabel('Region', fontsize=11)
    ax.set_title(f'Mixed model: esmu ~ z + (1|left_out_animal)\n{gwas_name}', 
                fontsize=11, pad=15)
    ax.grid(axis='x', alpha=0.3, linestyle=':')

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')

    plt.show()

# Example: Plot for BMI
bmi_gwas = [g for g in df_area['gwas'].unique() if 'BMI' in g][0]
plot_gwas_forest(bmi_gwas, df_area, save_path='forest_plot_BMI.pdf')

# You can also create plots for all GWAS
# for gwas in df_area['gwas'].unique():
#     plot_gwas_forest(gwas, df_area, save_path=f'forest_plot_{gwas}.pdf')


# In[36]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from statsmodels.stats.multitest import multipletests

# Load results
df_area = pd.read_csv('mixedmodel_results_by_area.csv')

# Apply Global FDR correction (RECOMMENDED)
reject_fdr, pvals_fdr, _, _ = multipletests(
    df_area['pvalue'], 
    alpha=0.05, 
    method='fdr_bh'
)
df_area['significant_fdr'] = reject_fdr
df_area['pvalue_fdr_corrected'] = pvals_fdr

# Also add other correction methods for comparison
global_bonf_threshold = 0.05 / len(df_area)
df_area['significant_global_bonf'] = df_area['pvalue'] < global_bonf_threshold

# Print summary
print("CORRECTION METHOD COMPARISON:")
print(f"Original (per-GWAS Bonf, p<0.05/7): {df_area['significant'].sum()} / {len(df_area)}")
print(f"Global Bonferroni (p<{global_bonf_threshold:.6f}): {df_area['significant_global_bonf'].sum()} / {len(df_area)}")
print(f"Global FDR (q<0.05): {df_area['significant_fdr'].sum()} / {len(df_area)}")

# Define area order
area_order = ['BST', 'CEA', 'LHA', 'PAG', 'PBN', 'PVH', 'PVT']

# Pivot to create GWAS x Area matrix
pivot_beta = df_area.pivot(index='gwas', columns='area', values='beta_z')
pivot_sig_fdr = df_area.pivot(index='gwas', columns='area', values='significant_fdr')

# Reorder columns by area
pivot_beta = pivot_beta[area_order]
pivot_sig_fdr = pivot_sig_fdr[area_order]

# Calculate figure size to maintain square cells
n_rows = len(pivot_beta.index)
n_cols = len(area_order)
fig_width = 10
fig_height = fig_width * (n_rows / n_cols)

# Create figure
fig, ax = plt.subplots(figsize=(fig_width, fig_height))

# Determine color scale limits (symmetric around 0)
vmax = np.abs(pivot_beta.values).max()
vmin = -vmax

# Plot heatmap with RdBu_r colormap (red=positive, blue=negative)
im = ax.imshow(pivot_beta.values, aspect='equal', cmap='RdBu_r', 
               vmin=vmin, vmax=vmax, interpolation='nearest')

# Add asterisks for FDR-significant results
for i in range(len(pivot_beta.index)):
    for j in range(len(area_order)):
        if pivot_sig_fdr.iloc[i, j]:
            ax.text(j, i, '*', ha='center', va='center', 
                   color='black', fontsize=16, fontweight='bold',
                   clip_on=False)

# Add faint grid
ax.set_xticks(np.arange(-.5, len(area_order), 1), minor=True)
ax.set_yticks(np.arange(-.5, len(pivot_beta.index), 1), minor=True)
ax.grid(which='minor', color='gray', linestyle='-', linewidth=0.5, alpha=0.3)

# Set labels
ax.set_yticks(range(len(pivot_beta.index)))
ax.set_yticklabels(pivot_beta.index, fontsize=9)
ax.set_xticks(range(len(area_order)))
ax.set_xticklabels(area_order, fontsize=12, fontweight='bold')

# Position x-axis labels at top
ax.xaxis.tick_top()
ax.xaxis.set_label_position('top')

# Add colorbar
cbar = plt.colorbar(im, ax=ax, label='Beta coefficient (effect size)', 
                    fraction=0.046, pad=0.04)

# Add title with FDR correction info
ax.set_title('Mixed Model Results: ESμ ~ GWAS z-score + (1|animal)\n* = FDR significant (q < 0.05, Benjamini-Hochberg)', 
             fontsize=12, fontweight='bold', pad=20, y=1.08)

plt.tight_layout()
plt.savefig('gwas_area_heatmap_fdr.pdf', bbox_inches='tight', dpi=300)
plt.show()

# Print detailed summary statistics
print(f"\n=== SUMMARY STATISTICS ===")
print(f"Total GWAS traits: {len(pivot_beta.index)}")
print(f"Total areas: {len(area_order)}")
print(f"Total tests: {len(df_area)}")
print(f"FDR-significant results: {pivot_sig_fdr.sum().sum()}")
print(f"Percentage significant: {100 * pivot_sig_fdr.sum().sum() / len(df_area):.1f}%")

print(f"\n=== FDR-SIGNIFICANT RESULTS PER AREA ===")
print(pivot_sig_fdr.sum().sort_values(ascending=False))

print(f"\n=== TOP 15 GWAS BY NUMBER OF FDR-SIGNIFICANT AREAS ===")
print(pivot_sig_fdr.sum(axis=1).sort_values(ascending=False).head(15))

# Save corrected results
df_area.to_csv('mixedmodel_results_with_fdr.csv', index=False)
print(f"\nSaved results with FDR correction to: mixedmodel_results_with_fdr.csv")


# In[38]:


import pandas as pd
import numpy as np
from statsmodels.regression.mixed_linear_model import MixedLM
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm

# Load data
df_cellect = pd.read_csv('df_cellect.csv')
combined_LOO = pd.read_csv('combined_LOO.csv')

print("="*80)
print("STAGE 1: Identify GWAS with significant cell type enrichments")
print("="*80)

# Bonferroni correction for cell types (41 tests per GWAS)
n_celltypes = df_cellect['annotation'].nunique()
celltype_threshold = 0.05 / n_celltypes

print(f"\nNumber of cell types: {n_celltypes}")
print(f"Bonferroni threshold: 0.05 / {n_celltypes} = {celltype_threshold:.6f}")

# Find GWAS with at least one significant cell type
df_cellect['significant_celltype'] = df_cellect['pvalue'] < celltype_threshold

gwas_with_enrichment = df_cellect[df_cellect['significant_celltype']]['gwas'].unique()

print(f"\nTotal GWAS tested: {df_cellect['gwas'].nunique()}")
print(f"GWAS with ≥1 significant cell type: {len(gwas_with_enrichment)}")
print(f"\nGWAS passing Stage 1:")
for gwas in sorted(gwas_with_enrichment):
    n_sig = df_cellect[(df_cellect['gwas'] == gwas) & 
                       (df_cellect['significant_celltype'])]['annotation'].nunique()
    print(f"  {gwas}: {n_sig} significant cell type(s)")

print("\n" + "="*80)
print("STAGE 2: Test projection associations for qualified GWAS only")
print("="*80)

# Prepare the results list
results_list = []

# Define areas
area_order = ['BST', 'CEA', 'LHA', 'PAG', 'PBN', 'PVH', 'PVT']

# For each qualified GWAS and each area
for gwas in tqdm(gwas_with_enrichment, desc="Processing qualified GWAS"):
    # Get GWAS data for this trait
    gwas_data = df_cellect[df_cellect['gwas'] == gwas][['annotation', 'z']].copy()
    gwas_data.columns = ['cluster', 'z_score']

    for area in area_order:
        # Filter combined_LOO for this specific area
        area_data = combined_LOO[combined_LOO['area'] == area].copy()

        # Merge with GWAS z-scores
        merged_data = area_data.merge(gwas_data, on='cluster', how='inner')

        if len(merged_data) == 0:
            continue

        # Add intercept
        merged_data['intercept'] = 1.0

        # Run mixed model: esmu ~ z_score + (1|animal)
        try:
            model = MixedLM(
                endog=merged_data['esmu'],
                exog=merged_data[['intercept', 'z_score']],
                groups=merged_data['left_out_animal']
            )
            result = model.fit(method='lbfgs', maxiter=1000)

            # Extract results
            results_list.append({
                'gwas': gwas,
                'area': area,
                'beta_intercept': result.params['intercept'],
                'beta_z': result.params['z_score'],
                'se_intercept': result.bse['intercept'],
                'se_z': result.bse['z_score'],
                't_stat': result.tvalues['z_score'],
                'pvalue': result.pvalues['z_score'],
                'ci_lower': result.conf_int().loc['z_score', 0],
                'ci_upper': result.conf_int().loc['z_score', 1],
                'converged': result.converged,
                'n_obs': len(merged_data),
                'n_animals': merged_data['left_out_animal'].nunique(),
                'n_celltypes': merged_data['cluster'].nunique()
            })
        except Exception as e:
            print(f"Error for {gwas} - {area}: {str(e)}")

# Create DataFrame
df_stage2 = pd.DataFrame(results_list)
df_stage2['mlog10p'] = -np.log10(df_stage2['pvalue'])

print(f"\nStage 2 complete: {len(df_stage2)} tests run")
print(f"({len(gwas_with_enrichment)} GWAS × 7 areas)")

# Apply FDR correction across all Stage 2 tests
reject_fdr, pvals_fdr, _, _ = multipletests(
    df_stage2['pvalue'], 
    alpha=0.05, 
    method='fdr_bh'
)
df_stage2['significant_fdr'] = reject_fdr
df_stage2['pvalue_fdr_corrected'] = pvals_fdr

# Also add global Bonferroni for comparison
global_bonf_threshold = 0.05 / len(df_stage2)
df_stage2['significant_global_bonf'] = df_stage2['pvalue'] < global_bonf_threshold

print("\n" + "="*80)
print("MULTIPLE TESTING CORRECTION (Stage 2 tests only)")
print("="*80)
print(f"Total Stage 2 tests: {len(df_stage2)}")
print(f"Global Bonferroni (p<{global_bonf_threshold:.6f}): {df_stage2['significant_global_bonf'].sum()} significant")
print(f"Global FDR (q<0.05): {df_stage2['significant_fdr'].sum()} significant")
print(f"Percentage significant (FDR): {100 * df_stage2['significant_fdr'].sum() / len(df_stage2):.1f}%")

# Save results
df_stage2.to_csv('mixedmodel_twostage_results.csv', index=False)
print(f"\nResults saved to: mixedmodel_twostage_results.csv")

# Summary by area
print("\n" + "="*80)
print("FDR-SIGNIFICANT RESULTS BY AREA")
print("="*80)
pivot_sig = df_stage2.pivot(index='gwas', columns='area', values='significant_fdr')
pivot_sig = pivot_sig[area_order]
print(pivot_sig.sum().sort_values(ascending=False))

# Summary by GWAS
print("\n" + "="*80)
print("TOP GWAS BY NUMBER OF SIGNIFICANT AREAS (FDR)")
print("="*80)
print(pivot_sig.sum(axis=1).sort_values(ascending=False).head(15))

print("\n" + "="*80)
print("COMPARISON WITH ORIGINAL ANALYSIS")
print("="*80)
df_original = pd.read_csv('mixedmodel_results_by_area.csv')
print(f"Original analysis: {df_original['gwas'].nunique()} GWAS tested, {len(df_original)} total tests")
print(f"Two-stage analysis: {len(gwas_with_enrichment)} GWAS tested, {len(df_stage2)} total tests")
print(f"Reduction: {100 * (1 - len(df_stage2)/len(df_original)):.1f}% fewer tests")


# In[39]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load the two-stage results
df_stage2 = pd.read_csv('mixedmodel_twostage_results.csv')

# Define area order
area_order = ['BST', 'CEA', 'LHA', 'PAG', 'PBN', 'PVH', 'PVT']

# Pivot to create GWAS x Area matrix
pivot_beta = df_stage2.pivot(index='gwas', columns='area', values='beta_z')
pivot_sig_fdr = df_stage2.pivot(index='gwas', columns='area', values='significant_fdr')

# Reorder columns by area
pivot_beta = pivot_beta[area_order]
pivot_sig_fdr = pivot_sig_fdr[area_order]

# Calculate figure size to maintain square cells
n_rows = len(pivot_beta.index)
n_cols = len(area_order)
fig_width = 10
fig_height = fig_width * (n_rows / n_cols)

# Create figure
fig, ax = plt.subplots(figsize=(fig_width, fig_height))

# Determine color scale limits (symmetric around 0)
vmax = np.abs(pivot_beta.values).max()
vmin = -vmax

# Plot heatmap with RdBu_r colormap (red=positive, blue=negative)
im = ax.imshow(pivot_beta.values, aspect='equal', cmap='RdBu_r', 
               vmin=vmin, vmax=vmax, interpolation='nearest')

# Add asterisks for FDR-significant results
for i in range(len(pivot_beta.index)):
    for j in range(len(area_order)):
        if pivot_sig_fdr.iloc[i, j]:
            ax.text(j, i, '*', ha='center', va='center', 
                   color='black', fontsize=16, fontweight='bold',
                   clip_on=False)

# Add faint grid
ax.set_xticks(np.arange(-.5, len(area_order), 1), minor=True)
ax.set_yticks(np.arange(-.5, len(pivot_beta.index), 1), minor=True)
ax.grid(which='minor', color='gray', linestyle='-', linewidth=0.5, alpha=0.3)

# Set labels
ax.set_yticks(range(len(pivot_beta.index)))
ax.set_yticklabels(pivot_beta.index, fontsize=10)
ax.set_xticks(range(len(area_order)))
ax.set_xticklabels(area_order, fontsize=12, fontweight='bold')

# Position x-axis labels at top
ax.xaxis.tick_top()
ax.xaxis.set_label_position('top')

# Add colorbar
cbar = plt.colorbar(im, ax=ax, label='Beta coefficient (effect size)', 
                    fraction=0.046, pad=0.04)

# Add title with two-stage info
title_text = ('Two-Stage Analysis: ESμ ~ GWAS z-score + (1|animal)\n'
              'Stage 1: GWAS with ≥1 significant cell type (Bonferroni)\n'
              'Stage 2: * = FDR significant projection (q < 0.05)')
ax.set_title(title_text, fontsize=11, fontweight='bold', pad=25, y=1.12)

plt.tight_layout()
plt.savefig('gwas_area_heatmap_twostage.pdf', bbox_inches='tight', dpi=300)
plt.show()

# Print summary statistics
print("="*80)
print("TWO-STAGE HEATMAP SUMMARY")
print("="*80)
print(f"\nGWAS included (passed Stage 1): {len(pivot_beta.index)}")
print(f"Brain areas: {len(area_order)}")
print(f"Total tests: {len(df_stage2)}")
print(f"FDR-significant results: {pivot_sig_fdr.sum().sum()}")
print(f"Percentage significant: {100 * pivot_sig_fdr.sum().sum() / len(df_stage2):.1f}%")

print(f"\n{'Area':<6} {'Significant':>11}")
print("-"*20)
for area in area_order:
    n_sig = pivot_sig_fdr[area].sum()
    print(f"{area:<6} {n_sig:>11}")

print(f"\n{'GWAS':<50} {'N Areas':>8}")
print("-"*60)
gwas_summary = pivot_sig_fdr.sum(axis=1).sort_values(ascending=False)
for gwas, n_areas in gwas_summary.head(15).items():
    gwas_short = gwas[:47] + '...' if len(gwas) > 50 else gwas
    print(f"{gwas_short:<50} {int(n_areas):>8}")

# Create a comparison plot showing before/after filtering
df_original = pd.read_csv('mixedmodel_results_by_area.csv')

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))

# Original analysis p-value distribution
ax1.hist(df_original['pvalue'], bins=50, edgecolor='black', alpha=0.7)
ax1.axvline(0.05, color='red', linestyle='--', label='p=0.05')
ax1.set_xlabel('P-value', fontsize=11)
ax1.set_ylabel('Frequency', fontsize=11)
ax1.set_title(f'Original Analysis\n({len(df_original)} tests, all GWAS)', fontsize=11)
ax1.legend()

# Two-stage analysis p-value distribution
ax2.hist(df_stage2['pvalue'], bins=50, edgecolor='black', alpha=0.7, color='green')
ax2.axvline(0.05, color='red', linestyle='--', label='p=0.05')
ax2.set_xlabel('P-value', fontsize=11)
ax2.set_ylabel('Frequency', fontsize=11)
ax2.set_title(f'Two-Stage Analysis\n({len(df_stage2)} tests, filtered GWAS)', fontsize=11)
ax2.legend()

plt.tight_layout()
plt.savefig('pvalue_distribution_comparison.pdf', dpi=300, bbox_inches='tight')
plt.show()

print(f"\nPlots saved:")
print(f"  - gwas_area_heatmap_twostage.pdf")
print(f"  - pvalue_distribution_comparison.pdf")


# In[41]:


import pandas as pd
import numpy as np
from statsmodels.regression.mixed_linear_model import MixedLM
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm

# Load data
df_cellect = pd.read_csv('df_cellect.csv')
combined_LOO = pd.read_csv('combined_LOO.csv')

print("="*80)
print("WEIGHTED MIXED MODEL ANALYSIS (via data replication)")
print("="*80)
print("\nWeighting strategy: -log10(p-value) from CELLECT")
print("Cell types with stronger GWAS associations get replicated more")

# Prepare the results list
results_list = []

# Get unique GWAS traits and areas
gwas_traits = df_cellect['gwas'].unique()
area_order = ['BST', 'CEA', 'LHA', 'PAG', 'PBN', 'PVH', 'PVT']

# For each GWAS and each area
for gwas in tqdm(gwas_traits, desc="Processing GWAS"):
    # Get GWAS data for this trait
    gwas_data = df_cellect[df_cellect['gwas'] == gwas][['annotation', 'z', 'pvalue']].copy()
    gwas_data.columns = ['cluster', 'z_score', 'cellect_pvalue']

    # Calculate weights from CELLECT p-values
    # Use -log10(p) as weights, with a floor to avoid extreme weights
    gwas_data['weight_raw'] = -np.log10(gwas_data['cellect_pvalue'] + 1e-300)

    # Normalize weights to sum to number of cell types (so average weight = 1)
    n_celltypes = len(gwas_data)
    gwas_data['weight'] = gwas_data['weight_raw'] / gwas_data['weight_raw'].mean()

    # Cap extreme weights to avoid single cell types dominating
    weight_cap = gwas_data['weight'].quantile(0.90)
    gwas_data['weight_capped'] = np.minimum(gwas_data['weight'], weight_cap)

    for area in area_order:
        # Filter combined_LOO for this specific area
        area_data = combined_LOO[combined_LOO['area'] == area].copy()

        # Merge with GWAS z-scores and weights
        merged_data = area_data.merge(gwas_data, on='cluster', how='inner')

        if len(merged_data) == 0:
            continue

        # Apply weights by replicating observations
        # Convert weights to integer counts (multiply by 10 and round)
        merged_data['weight_count'] = np.maximum(1, np.round(merged_data['weight_capped'] * 10).astype(int))

        # Replicate rows based on weights
        weighted_data = merged_data.loc[merged_data.index.repeat(merged_data['weight_count'])].copy()

        # Add intercept
        weighted_data['intercept'] = 1.0

        # Run mixed model on weighted (replicated) data
        try:
            model = MixedLM(
                endog=weighted_data['esmu'],
                exog=weighted_data[['intercept', 'z_score']],
                groups=weighted_data['left_out_animal']
            )
            result = model.fit(method='lbfgs', maxiter=1000)

            # Extract results
            results_list.append({
                'gwas': gwas,
                'area': area,
                'beta_intercept': result.params['intercept'],
                'beta_z': result.params['z_score'],
                'se_intercept': result.bse['intercept'],
                'se_z': result.bse['z_score'],
                't_stat': result.tvalues['z_score'],
                'pvalue': result.pvalues['z_score'],
                'ci_lower': result.conf_int().loc['z_score', 0],
                'ci_upper': result.conf_int().loc['z_score', 1],
                'converged': result.converged,
                'n_obs_original': len(merged_data),
                'n_obs_weighted': len(weighted_data),
                'n_animals': merged_data['left_out_animal'].nunique(),
                'n_celltypes': merged_data['cluster'].nunique(),
                'mean_weight': merged_data['weight_capped'].mean(),
                'max_weight': merged_data['weight_capped'].max(),
                'min_weight': merged_data['weight_capped'].min()
            })
        except Exception as e:
            print(f"Error for {gwas} - {area}: {str(e)}")

# Create DataFrame
df_weighted = pd.DataFrame(results_list)

if len(df_weighted) == 0:
    print("ERROR: No results generated!")
else:
    df_weighted['mlog10p'] = -np.log10(df_weighted['pvalue'])

    print(f"\nAnalysis complete: {len(df_weighted)} tests run")

    # Apply FDR correction
    reject_fdr, pvals_fdr, _, _ = multipletests(
        df_weighted['pvalue'], 
        alpha=0.05, 
        method='fdr_bh'
    )
    df_weighted['significant_fdr'] = reject_fdr
    df_weighted['pvalue_fdr_corrected'] = pvals_fdr

    # Also add global Bonferroni
    global_bonf_threshold = 0.05 / len(df_weighted)
    df_weighted['significant_global_bonf'] = df_weighted['pvalue'] < global_bonf_threshold

    print("\n" + "="*80)
    print("MULTIPLE TESTING CORRECTION")
    print("="*80)
    print(f"Total tests: {len(df_weighted)}")
    print(f"Global Bonferroni (p<{global_bonf_threshold:.6f}): {df_weighted['significant_global_bonf'].sum()} significant")
    print(f"Global FDR (q<0.05): {df_weighted['significant_fdr'].sum()} significant")
    print(f"Percentage significant (FDR): {100 * df_weighted['significant_fdr'].sum() / len(df_weighted):.1f}%")

    # Save results
    df_weighted.to_csv('mixedmodel_weighted_results.csv', index=False)
    print(f"\nResults saved to: mixedmodel_weighted_results.csv")

    # Summary by area
    print("\n" + "="*80)
    print("FDR-SIGNIFICANT RESULTS BY AREA")
    print("="*80)
    pivot_sig = df_weighted.pivot(index='gwas', columns='area', values='significant_fdr')
    pivot_sig = pivot_sig[area_order]
    print(pivot_sig.sum().sort_values(ascending=False))

    # Summary by GWAS
    print("\n" + "="*80)
    print("TOP 15 GWAS BY NUMBER OF SIGNIFICANT AREAS (FDR)")
    print("="*80)
    print(pivot_sig.sum(axis=1).sort_values(ascending=False).head(15))

    # Weight statistics
    print("\n" + "="*80)
    print("WEIGHT DISTRIBUTION STATISTICS")
    print("="*80)
    print(f"Mean weight: {df_weighted['mean_weight'].mean():.2f}")
    print(f"Mean observations after weighting: {df_weighted['n_obs_weighted'].mean():.1f}")
    print(f"Mean observations before weighting: {df_weighted['n_obs_original'].mean():.1f}")
    print(f"Average replication factor: {(df_weighted['n_obs_weighted'] / df_weighted['n_obs_original']).mean():.2f}x")


# In[42]:


import pandas as pd
import numpy as np
from statsmodels.formula.api import mixedlm
from statsmodels.stats.multitest import multipletests

# Load data
df_cellect = pd.read_csv('df_cellect.csv')
combined_LOO = pd.read_csv('combined_LOO.csv')

print("="*80)
print("UNIFIED MIXED MODEL: All GWAS and Areas Together")
print("="*80)

# Stage 1: Filter to GWAS with significant cell types
n_celltypes = df_cellect['annotation'].nunique()
celltype_threshold = 0.05 / n_celltypes

df_cellect['significant_celltype'] = df_cellect['pvalue'] < celltype_threshold
gwas_with_enrichment = df_cellect[df_cellect['significant_celltype']]['gwas'].unique()

print(f"\nStage 1: {len(gwas_with_enrichment)} GWAS with significant cell type enrichment")

# Prepare unified dataset
gwas_data = df_cellect[df_cellect['gwas'].isin(gwas_with_enrichment)][['gwas', 'annotation', 'z']].copy()
gwas_data.columns = ['gwas', 'cluster', 'z_score']

# Merge with projection data
merged_data = combined_LOO.merge(gwas_data, on='cluster', how='inner')

print(f"\nUnified dataset: {len(merged_data)} observations")
print(f"  {merged_data['gwas'].nunique()} GWAS")
print(f"  {merged_data['area'].nunique()} areas")
print(f"  {merged_data['cluster'].nunique()} cell types")
print(f"  {merged_data['left_out_animal'].nunique()} animals")

# Create dummy variables for area (reference = BST)
merged_data = pd.get_dummies(merged_data, columns=['area'], drop_first=False, prefix='area')
area_cols = [col for col in merged_data.columns if col.startswith('area_')]

# Remove reference category (BST)
area_cols.remove('area_BST')

print("\n" + "="*80)
print("MODEL 1: Main Effect Only")
print("="*80)
print("ESμ ~ GWAS_z-score + (1|animal) + (1|GWAS)")

# Model 1: Main effect only
formula1 = 'esmu ~ z_score'
model1 = mixedlm(formula1, merged_data, groups=merged_data['left_out_animal'],
                 re_formula='1', vc_formula={'gwas': '0 + C(gwas)'})
result1 = model1.fit(method='lbfgs', maxiter=1000)

print("\n### RESULTS ###")
print(result1.summary().tables[1])
print(f"\nMain effect beta: {result1.params['z_score']:.6f}")
print(f"P-value: {result1.pvalues['z_score']:.2e}")
print(f"95% CI: [{result1.conf_int().loc['z_score', 0]:.6f}, {result1.conf_int().loc['z_score', 1]:.6f}]")

print("\n" + "="*80)
print("MODEL 2: With Area Interactions")
print("="*80)
print("ESμ ~ GWAS_z-score × Area + (1|animal) + (1|GWAS)")

# Model 2: With interactions
interaction_terms = ' + '.join([f'z_score:{col}' for col in area_cols])
formula2 = f'esmu ~ z_score + {" + ".join(area_cols)} + {interaction_terms}'

model2 = mixedlm(formula2, merged_data, groups=merged_data['left_out_animal'],
                 re_formula='1', vc_formula={'gwas': '0 + C(gwas)'})
result2 = model2.fit(method='lbfgs', maxiter=1000)

print("\n### RESULTS ###")
print(result2.summary().tables[1])

# Extract interaction effects
print("\n### AREA-SPECIFIC EFFECTS ###")
print(f"{'Area':<15} {'Effect (β)':<12} {'P-value':<12} {'95% CI':<25}")
print("-"*65)

# BST (reference)
bst_effect = result2.params['z_score']
bst_pval = result2.pvalues['z_score']
bst_ci = result2.conf_int().loc['z_score']
print(f"{'BST (ref)':<15} {bst_effect:>11.6f} {bst_pval:>11.2e} [{bst_ci[0]:.6f}, {bst_ci[1]:.6f}]")

# Other areas (reference + interaction)
area_names = {'area_CEA': 'CEA', 'area_LHA': 'LHA', 'area_PAG': 'PAG', 
              'area_PBN': 'PBN', 'area_PVH': 'PVH', 'area_PVT': 'PVT'}

for col in area_cols:
    interaction_param = f'z_score:{col}'
    if interaction_param in result2.params:
        total_effect = bst_effect + result2.params[interaction_param]
        pval = result2.pvalues[interaction_param]
        ci = result2.conf_int().loc[interaction_param]
        print(f"{area_names[col]:<15} {total_effect:>11.6f} {pval:>11.2e} [interaction: {result2.params[interaction_param]:.6f}]")

# Model comparison
print("\n" + "="*80)
print("MODEL COMPARISON")
print("="*80)
print(f"Model 1 (main effect) AIC: {result1.aic:.1f}")
print(f"Model 2 (with interactions) AIC: {result2.aic:.1f}")
print(f"AIC difference: {result1.aic - result2.aic:.1f}")

if result2.aic < result1.aic:
    print("\n✓ Model 2 (with interactions) is preferred")
else:
    print("\n✓ Model 1 (main effect only) is preferred")

# Save results
results_summary = {
    'model': ['Main Effect', 'With Interactions'],
    'aic': [result1.aic, result2.aic],
    'main_effect_beta': [result1.params['z_score'], result2.params['z_score']],
    'main_effect_pval': [result1.pvalues['z_score'], result2.pvalues['z_score']]
}

df_model_comparison = pd.DataFrame(results_summary)
df_model_comparison.to_csv('unified_model_results.csv', index=False)

print("\n✓ Results saved to: unified_model_results.csv")


# In[43]:


import pandas as pd
import numpy as np
from scipy import stats

# Load two-stage results
df_twostage = pd.read_csv('mixedmodel_twostage_results.csv')

# For each GWAS, test if beta coefficients differ across areas
print("="*80)
print("TESTING HETEROGENEITY ACROSS AREAS")
print("="*80)

for gwas in df_twostage['gwas'].unique():
    gwas_data = df_twostage[df_twostage['gwas'] == gwas]

    # Chi-square test for heterogeneity
    # Under null: all areas have same effect
    betas = gwas_data['beta_z'].values
    ses = gwas_data['se_z'].values

    # Weighted mean
    weights = 1 / (ses ** 2)
    weighted_mean = np.sum(betas * weights) / np.sum(weights)

    # Q statistic (heterogeneity test)
    Q = np.sum(weights * (betas - weighted_mean) ** 2)
    df = len(betas) - 1
    p_heterogeneity = 1 - stats.chi2.cdf(Q, df)

    # I-squared (proportion of variance due to heterogeneity)
    I2 = max(0, (Q - df) / Q * 100)

    print(f"\n{gwas}")
    print(f"  Q statistic: {Q:.2f} (df={df})")
    print(f"  P-heterogeneity: {p_heterogeneity:.4f}")
    print(f"  I²: {I2:.1f}%")

    if p_heterogeneity < 0.05:
        print(f"  ✓ Significant heterogeneity - area-specific effects are real")
        print(f"  Beta range: [{betas.min():.4f}, {betas.max():.4f}]")


# In[ ]:


#Neuronal Cell Types (3):

#04_Crnde_Lmx1a
#Significant for all 3 BMI GWAS (Yengo2018, UPDATE_Yengo2018, Loh2018)
#z-scores: 3.24-3.51

#24_Nts_Ecel1
#Significant for LEPR (Leptin Receptor)
#z = 3.19
#Neurotensin (Nts) expressing neuron


#02_Fezf1_Lama1

#Significant for BMI_Loh2018 (secondary)
#z = 3.23



#Vascular Cell Types (3):

#VLMC (Vascular Leptomeningeal Cells)
#Significant for Forced Expiratory Volume (z = 4.82, strongest signal!)
#Significant for Pulse (z = 3.39)


#Endothelial_cell
#Significant for FLT4/VEGFR3 (z = 3.16) - makes perfect sense! (VEGFR3 is on endothelial cells)
#Significant for Diastolic BP (z = 3.53)


#Pericyte
#Significant for Diastolic BP (z = 3.70)
#Vascular support cells


# In[2]:


#adata = sc.read("allgenes_bipca_hypomap.h5ad")
adata = sc.read("allgenes_bipca_annotated.h5ad", cache=True)
adata


# In[4]:


celltypes = ['04_Crnde_Lmx1a', '27_Lmx1a_Ttc6', '32_Gucy2c_Ttc6', '36_Ebf2_Crhbp']

# Assign as plain strings first (breaks the inherited category dtype)
adata.obs['plot'] = adata.obs['detailed_celltype_annotation'].astype(str).where(
    adata.obs['detailed_celltype_annotation'].isin(celltypes)
)
# Re-cast as categorical with ONLY the 4 categories
adata.obs['plot'] = pd.Categorical(
    adata.obs['plot'],
    categories=celltypes
)

# Pull original colors
all_cats   = list(adata.obs['detailed_celltype_annotation'].cat.categories)
all_colors = adata.uns['detailed_celltype_annotation_colors']
color_map  = dict(zip(all_cats, all_colors))
palette    = [color_map[ct] for ct in celltypes]

ax = sc.pl.embedding(
    adata,
    basis='tsne',
    color='plot',
    groups=celltypes,
    palette=palette,
    na_color='#E0E0E0',
    na_in_legend=False,
    legend_loc='on data',
    legend_fontsize='x-small',
    return_fig=True,
)

ax = plt.gca()
adjustText.adjust_text(
    ax.texts,
    ax=ax,
    arrowprops=dict(arrowstyle='-', color='black', lw=0.5)
)

plt.savefig('paper/Fig_4a.pdf', bbox_inches='tight')
plt.show()


# In[3]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import spearmanr, pearsonr, kruskal
import statsmodels.api as sm

# -----------------------------
# Helpers
# -----------------------------
def extract_stat_p(res):
    """Compatibility helper for scipy correlation outputs across versions."""
    if hasattr(res, "statistic") and hasattr(res, "pvalue"):
        return res.statistic, res.pvalue
    return res[0], res[1]

def make_connectivity_strata(series):
    """
    Stratify connectivity as:
      - Zero
      - tertiles among nonzero values
    Handles tied values by dropping duplicate cut points if needed.
    """
    out = pd.Series(index=series.index, dtype="object")

    zero_mask = series == 0
    nonzero_mask = series > 0

    out.loc[zero_mask] = "Zero"

    nonzero_vals = series.loc[nonzero_mask]

    if nonzero_vals.empty:
        return pd.Categorical(out, categories=["Zero"], ordered=True)

    # Try tertiles among nonzero values; drop duplicate bin edges if needed
    qcut_result, bins = pd.qcut(
        nonzero_vals,
        q=3,
        retbins=True,
        duplicates="drop"
    )

    n_bins = len(bins) - 1

    if n_bins == 3:
        labels = ["Low nonzero", "Mid nonzero", "High nonzero"]
    elif n_bins == 2:
        labels = ["Low nonzero", "High nonzero"]
    elif n_bins == 1:
        labels = ["Nonzero"]
    else:
        raise ValueError("Unexpected number of bins created for nonzero connectivity.")

    out.loc[nonzero_mask] = pd.qcut(
        nonzero_vals,
        q=n_bins,
        labels=labels,
        duplicates="drop"
    ).astype(str)

    categories = ["Zero"] + labels
    return pd.Categorical(out, categories=categories, ordered=True)

# -----------------------------
# 1) Load the two source tables
# -----------------------------
df_cellect = pd.read_csv("df_cellect.csv")

loo_candidates = ["combined_LOO.csv", "df_combined_LOO_new.csv"]
loo_path = next((p for p in loo_candidates if Path(p).exists()), None)
if loo_path is None:
    raise FileNotFoundError(f"None of these files were found: {loo_candidates}")

combined_LOO = pd.read_csv(loo_path)
print("Using projection file:", loo_path)

# ---------------------------------------------
# 2) Aggregate AAV ESmu to a cell-type measure
# ---------------------------------------------
# Mean across leave-one-out animals within each area x cell type
proj_mean = (
    combined_LOO
    .groupby(["cluster", "area"], as_index=False)["esmu"]
    .mean()
)

# Wide matrix: rows = cell type, cols = area
proj_wide = (
    proj_mean
    .pivot(index="cluster", columns="area", values="esmu")
    .fillna(0)
)

# Save area columns before adding summary metrics
area_cols = proj_wide.columns.tolist()

# Primary connectivity metric: total ESmu across all projection targets
proj_wide["connectivity_sum"] = proj_wide[area_cols].sum(axis=1)

# Additional summary metrics
proj_wide["connectivity_mean"] = proj_wide[area_cols].mean(axis=1)
proj_wide["connectivity_max"] = proj_wide[area_cols].max(axis=1)
proj_wide["n_proj_0p10"] = (proj_wide[area_cols] > 0.10).sum(axis=1)
proj_wide["n_proj_0p36"] = (proj_wide[area_cols] > 0.36).sum(axis=1)

proj_summary = proj_wide.reset_index()

# ----------------------------------------------------
# 3) Pull the BMI GWAS effect for each cell type
# ----------------------------------------------------
# Primary analysis: use the non-duplicate BMI study emphasized in the notebook
bmi_gwas = "BMI_Loh2018_UKBB_munged_txt"

bmi = (
    df_cellect.loc[
        df_cellect["gwas"] == bmi_gwas,
        ["annotation", "beta", "beta_se", "pvalue", "mlog10p", "z"]
    ]
    .rename(columns={
        "annotation": "cluster",
        "beta": "bmi_beta",
        "beta_se": "bmi_beta_se",
        "pvalue": "bmi_p",
        "mlog10p": "bmi_mlog10p",
        "z": "bmi_z",
    })
    .copy()
)

# ----------------------------------------------------
# 4) Merge projection connectivity with BMI enrichment
# ----------------------------------------------------
dat = proj_summary.merge(bmi, on="cluster", how="inner").copy()
print(f"Overlapping cell types in merged analysis: {dat['cluster'].nunique()}")

# Bonferroni threshold on the merged set
bonf = 0.05 / dat["cluster"].nunique()
dat["bmi_sig_bonf"] = dat["bmi_p"] < bonf

# Connectivity strata that handle zero inflation / tied values
dat["connectivity_strata"] = make_connectivity_strata(dat["connectivity_sum"])
print("\nConnectivity strata counts:")
print(dat["connectivity_strata"].value_counts(dropna=False).sort_index())

# ----------------------------------------------------
# 5) Correlation tests
# ----------------------------------------------------
# Use z as primary continuous BMI measure; beta as requested; mlog10p as robustness
spearman_z_res = spearmanr(dat["connectivity_sum"], dat["bmi_z"])
pearson_z_res = pearsonr(dat["connectivity_sum"], dat["bmi_z"])

spearman_beta_res = spearmanr(dat["connectivity_sum"], dat["bmi_beta"])
pearson_beta_res = pearsonr(dat["connectivity_sum"], dat["bmi_beta"])

spearman_mlog10p_res = spearmanr(dat["connectivity_sum"], dat["bmi_mlog10p"])

spearman_z_stat, spearman_z_p = extract_stat_p(spearman_z_res)
pearson_z_stat, pearson_z_p = extract_stat_p(pearson_z_res)
spearman_beta_stat, spearman_beta_p = extract_stat_p(spearman_beta_res)
pearson_beta_stat, pearson_beta_p = extract_stat_p(pearson_beta_res)
spearman_mlog10p_stat, spearman_mlog10p_p = extract_stat_p(spearman_mlog10p_res)

print("\n=== CORRELATIONS ===")
print(f"Spearman(connectivity_sum, BMI z):     rho={spearman_z_stat:.3f}, p={spearman_z_p:.3e}")
print(f"Pearson(connectivity_sum, BMI z):      r={pearson_z_stat:.3f}, p={pearson_z_p:.3e}")
print(f"Spearman(connectivity_sum, BMI beta):  rho={spearman_beta_stat:.3f}, p={spearman_beta_p:.3e}")
print(f"Pearson(connectivity_sum, BMI beta):   r={pearson_beta_stat:.3f}, p={pearson_beta_p:.3e}")
print(f"Spearman(connectivity_sum, -log10 p):  rho={spearman_mlog10p_stat:.3f}, p={spearman_mlog10p_p:.3e}")

# Robust linear models
X = sm.add_constant(dat["connectivity_sum"])
fit_z = sm.OLS(dat["bmi_z"], X).fit(cov_type="HC3")
fit_beta = sm.OLS(dat["bmi_beta"], X).fit(cov_type="HC3")

print("\n=== ROBUST OLS ===")
print("\nBMI z ~ connectivity_sum")
print(fit_z.summary())

print("\nBMI beta ~ connectivity_sum")
print(fit_beta.summary())

# ----------------------------------------------------
# 6) Stratified summary by connectivity strata
# ----------------------------------------------------
strat_summary = (
    dat.groupby("connectivity_strata", observed=True)
    .agg(
        n_celltypes=("cluster", "size"),
        mean_connectivity=("connectivity_sum", "mean"),
        mean_bmi_beta=("bmi_beta", "mean"),
        mean_bmi_z=("bmi_z", "mean"),
        mean_bmi_mlog10p=("bmi_mlog10p", "mean"),
        n_bonf_sig=("bmi_sig_bonf", "sum"),
    )
)

print("\n=== STRATIFIED SUMMARY ===")
print(strat_summary)

# Nonparametric test across strata
groups = [
    g["bmi_z"].values
    for _, g in dat.groupby("connectivity_strata", observed=True)
    if len(g) > 0
]
kw = kruskal(*groups)
print(f"\nKruskal-Wallis across connectivity strata for BMI z: H={kw.statistic:.3f}, p={kw.pvalue:.3e}")

# ----------------------------------------------------
# 7) Plot
# ----------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(13, 5))

# Scatter: connectivity vs BMI z
axes[0].scatter(dat["connectivity_sum"], dat["bmi_z"], alpha=0.8, s=60)

# Label Bonferroni-significant BMI cell types and top-connectivity cell types
label_mask = dat["bmi_sig_bonf"] | (
    dat["connectivity_sum"] >= dat["connectivity_sum"].quantile(0.90)
)
for _, row in dat.loc[label_mask].iterrows():
    axes[0].text(row["connectivity_sum"], row["bmi_z"], row["cluster"], fontsize=8)

# Add fitted line from robust OLS coefficients
xline = np.linspace(dat["connectivity_sum"].min(), dat["connectivity_sum"].max(), 200)
yline = fit_z.params["const"] + fit_z.params["connectivity_sum"] * xline
axes[0].plot(xline, yline, linewidth=2)

axes[0].set_xlabel("Aggregate connectivity (sum of mean AAV ESμ across areas)")
axes[0].set_ylabel("BMI cell-type enrichment (z-score)")
axes[0].set_title(f"Cell-type connectivity vs BMI enrichment\n{bmi_gwas}")

# Boxplot by strata
strata_order = [cat for cat in dat["connectivity_strata"].cat.categories if cat in dat["connectivity_strata"].astype(str).unique()]
strata_data = [
    dat.loc[dat["connectivity_strata"] == s, "bmi_z"].values
    for s in strata_order
]

axes[1].boxplot(strata_data, labels=strata_order)
axes[1].set_xlabel("Connectivity strata")
axes[1].set_ylabel("BMI cell-type enrichment (z-score)")
axes[1].set_title("BMI enrichment stratified by connectivity")

plt.tight_layout()
plt.show()

# ----------------------------------------------------
# 8) Helpful tables to inspect
# ----------------------------------------------------
print("\n=== TOP CONNECTED CELL TYPES ===")
print(
    dat.sort_values("connectivity_sum", ascending=False)[
        ["cluster", "connectivity_sum", "bmi_beta", "bmi_z", "bmi_mlog10p", "bmi_sig_bonf", "connectivity_strata"]
    ].head(15)
)

print("\n=== TOP BMI-ENRICHED CELL TYPES ===")
print(
    dat.sort_values("bmi_z", ascending=False)[
        ["cluster", "connectivity_sum", "bmi_beta", "bmi_z", "bmi_mlog10p", "bmi_sig_bonf", "connectivity_strata"]
    ].head(15)
)


# In[4]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import spearmanr, pearsonr
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

# -----------------------------
# Helpers
# -----------------------------
def extract_stat_p(res):
    if hasattr(res, "statistic") and hasattr(res, "pvalue"):
        return res.statistic, res.pvalue
    return res[0], res[1]

def make_area_strata(series):
    """
    Stratify one area's ESmu as:
      - Zero
      - tertiles among nonzero values
    """
    out = pd.Series(index=series.index, dtype="object")

    zero_mask = series == 0
    nonzero_mask = series > 0

    out.loc[zero_mask] = "Zero"
    nonzero_vals = series.loc[nonzero_mask]

    if nonzero_vals.empty:
        return pd.Categorical(out, categories=["Zero"], ordered=True)

    qcut_result, bins = pd.qcut(
        nonzero_vals,
        q=3,
        retbins=True,
        duplicates="drop"
    )
    n_bins = len(bins) - 1

    if n_bins == 3:
        labels = ["Low nonzero", "Mid nonzero", "High nonzero"]
    elif n_bins == 2:
        labels = ["Low nonzero", "High nonzero"]
    elif n_bins == 1:
        labels = ["Nonzero"]
    else:
        raise ValueError("Unexpected number of bins created.")

    out.loc[nonzero_mask] = pd.qcut(
        nonzero_vals,
        q=n_bins,
        labels=labels,
        duplicates="drop"
    ).astype(str)

    categories = ["Zero"] + labels
    return pd.Categorical(out, categories=categories, ordered=True)

# -----------------------------
# 1) Load tables
# -----------------------------
df_cellect = pd.read_csv("df_cellect.csv")

loo_candidates = ["combined_LOO.csv", "df_combined_LOO_new.csv"]
loo_path = next((p for p in loo_candidates if Path(p).exists()), None)
if loo_path is None:
    raise FileNotFoundError(f"None of these files were found: {loo_candidates}")

combined_LOO = pd.read_csv(loo_path)
print("Using projection file:", loo_path)

# -----------------------------
# 2) Aggregate ESmu by cluster x area
# -----------------------------
proj_mean = (
    combined_LOO
    .groupby(["cluster", "area"], as_index=False)["esmu"]
    .mean()
)

proj_wide = (
    proj_mean
    .pivot(index="cluster", columns="area", values="esmu")
    .fillna(0)
)

area_cols = proj_wide.columns.tolist()

# optional overall summaries
proj_wide["connectivity_sum"] = proj_wide[area_cols].sum(axis=1)
proj_wide["connectivity_mean"] = proj_wide[area_cols].mean(axis=1)

proj_summary = proj_wide.reset_index()

# -----------------------------
# 3) Pull BMI GWAS per-cell-type effect
# -----------------------------
bmi_gwas = "BMI_Loh2018_UKBB_munged_txt"

bmi = (
    df_cellect.loc[
        df_cellect["gwas"] == bmi_gwas,
        ["annotation", "beta", "beta_se", "pvalue", "mlog10p", "z"]
    ]
    .rename(columns={
        "annotation": "cluster",
        "beta": "bmi_beta",
        "beta_se": "bmi_beta_se",
        "pvalue": "bmi_p",
        "mlog10p": "bmi_mlog10p",
        "z": "bmi_z",
    })
    .copy()
)

# -----------------------------
# 4) Merge
# -----------------------------
dat = proj_summary.merge(bmi, on="cluster", how="inner").copy()
print(f"Overlapping cell types in merged analysis: {dat['cluster'].nunique()}")

bonf = 0.05 / dat["cluster"].nunique()
dat["bmi_sig_bonf"] = dat["bmi_p"] < bonf

# -----------------------------
# 5) Per-area association tests
# -----------------------------
results = []

for area in area_cols:
    x = dat[area].values
    y_z = dat["bmi_z"].values
    y_beta = dat["bmi_beta"].values

    # Correlations
    sp_z = spearmanr(x, y_z)
    pe_z = pearsonr(x, y_z)
    sp_beta = spearmanr(x, y_beta)
    pe_beta = pearsonr(x, y_beta)

    sp_z_stat, sp_z_p = extract_stat_p(sp_z)
    pe_z_stat, pe_z_p = extract_stat_p(pe_z)
    sp_beta_stat, sp_beta_p = extract_stat_p(sp_beta)
    pe_beta_stat, pe_beta_p = extract_stat_p(pe_beta)

    # Robust OLS: BMI z ~ area ESmu
    X = sm.add_constant(dat[area])
    fit_z = sm.OLS(dat["bmi_z"], X).fit(cov_type="HC3")
    fit_beta = sm.OLS(dat["bmi_beta"], X).fit(cov_type="HC3")

    ols_z_coef = fit_z.params[area]
    ols_z_p = fit_z.pvalues[area]

    ols_beta_coef = fit_beta.params[area]
    ols_beta_p = fit_beta.pvalues[area]

    results.append({
        "area": area,
        "n_nonzero": int((dat[area] > 0).sum()),
        "mean_esmu": dat[area].mean(),
        "spearman_rho_bmi_z": sp_z_stat,
        "spearman_p_bmi_z": sp_z_p,
        "pearson_r_bmi_z": pe_z_stat,
        "pearson_p_bmi_z": pe_z_p,
        "ols_coef_bmi_z": ols_z_coef,
        "ols_p_bmi_z": ols_z_p,
        "spearman_rho_bmi_beta": sp_beta_stat,
        "spearman_p_bmi_beta": sp_beta_p,
        "pearson_r_bmi_beta": pe_beta_stat,
        "pearson_p_bmi_beta": pe_beta_p,
        "ols_coef_bmi_beta": ols_beta_coef,
        "ols_p_bmi_beta": ols_beta_p,
    })

results_df = pd.DataFrame(results)

# FDR correction across the 7 areas
results_df["fdr_spearman_p_bmi_z"] = multipletests(
    results_df["spearman_p_bmi_z"], method="fdr_bh"
)[1]

results_df["fdr_pearson_p_bmi_z"] = multipletests(
    results_df["pearson_p_bmi_z"], method="fdr_bh"
)[1]

results_df["fdr_ols_p_bmi_z"] = multipletests(
    results_df["ols_p_bmi_z"], method="fdr_bh"
)[1]

results_df["fdr_spearman_p_bmi_beta"] = multipletests(
    results_df["spearman_p_bmi_beta"], method="fdr_bh"
)[1]

results_df["fdr_pearson_p_bmi_beta"] = multipletests(
    results_df["pearson_p_bmi_beta"], method="fdr_bh"
)[1]

results_df["fdr_ols_p_bmi_beta"] = multipletests(
    results_df["ols_p_bmi_beta"], method="fdr_bh"
)[1]

# Sort by strongest BMI z association
results_df = results_df.sort_values(
    ["fdr_spearman_p_bmi_z", "spearman_rho_bmi_z"],
    ascending=[True, False]
).reset_index(drop=True)

print("\n=== PER-AREA ASSOCIATION RESULTS ===")
print(results_df[[
    "area",
    "n_nonzero",
    "mean_esmu",
    "spearman_rho_bmi_z",
    "spearman_p_bmi_z",
    "fdr_spearman_p_bmi_z",
    "pearson_r_bmi_z",
    "pearson_p_bmi_z",
    "fdr_pearson_p_bmi_z",
    "ols_coef_bmi_z",
    "ols_p_bmi_z",
    "fdr_ols_p_bmi_z",
    "spearman_rho_bmi_beta",
    "spearman_p_bmi_beta",
    "fdr_spearman_p_bmi_beta",
]].to_string(index=False))

# -----------------------------
# 6) Stratified summaries within each area
# -----------------------------
stratified_tables = {}

for area in area_cols:
    tmp = dat.copy()
    tmp["area_strata"] = make_area_strata(tmp[area])

    summary = (
        tmp.groupby("area_strata", observed=True)
        .agg(
            n_celltypes=("cluster", "size"),
            mean_area_esmu=(area, "mean"),
            mean_bmi_z=("bmi_z", "mean"),
            mean_bmi_beta=("bmi_beta", "mean"),
            mean_bmi_mlog10p=("bmi_mlog10p", "mean"),
            n_bonf_sig=("bmi_sig_bonf", "sum"),
        )
    )
    stratified_tables[area] = summary

print("\n=== STRATIFIED SUMMARY BY AREA ===")
for area in area_cols:
    print(f"\n--- {area} ---")
    print(stratified_tables[area])

# -----------------------------
# 7) Scatter plots for all 7 areas
# -----------------------------
n_areas = len(area_cols)
ncols = 3
nrows = int(np.ceil(n_areas / ncols))

fig, axes = plt.subplots(nrows, ncols, figsize=(5*ncols, 4*nrows))
axes = np.array(axes).reshape(-1)

results_lookup = results_df.set_index("area")

for i, area in enumerate(area_cols):
    ax = axes[i]
    ax.scatter(dat[area], dat["bmi_z"], alpha=0.8, s=55)

    # Fit line
    X = sm.add_constant(dat[area])
    fit = sm.OLS(dat["bmi_z"], X).fit(cov_type="HC3")
    xline = np.linspace(dat[area].min(), dat[area].max(), 200)
    yline = fit.params["const"] + fit.params[area] * xline
    ax.plot(xline, yline, linewidth=2)

    # Label Bonferroni BMI cell types
    for _, row in dat.loc[dat["bmi_sig_bonf"]].iterrows():
        ax.text(row[area], row["bmi_z"], row["cluster"], fontsize=7)

    rho = results_lookup.loc[area, "spearman_rho_bmi_z"]
    p = results_lookup.loc[area, "fdr_spearman_p_bmi_z"]

    ax.set_title(f"{area}\nSpearman rho={rho:.2f}, FDR p={p:.3g}")
    ax.set_xlabel(f"{area} ESmu")
    ax.set_ylabel("BMI z-score")

# Hide unused axes
for j in range(i + 1, len(axes)):
    axes[j].axis("off")

plt.tight_layout()
plt.show()

# -----------------------------
# 8) Top cell types for each area
# -----------------------------
print("\n=== TOP CONNECTED CELL TYPES WITHIN EACH AREA ===")
for area in area_cols:
    print(f"\n--- {area} ---")
    cols = ["cluster", area, "bmi_z", "bmi_beta", "bmi_mlog10p", "bmi_sig_bonf"]
    print(
        dat.sort_values(area, ascending=False)[cols].head(10).to_string(index=False)
    )


# In[5]:


import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.api as sm

from scipy.stats import spearmanr, pearsonr
from statsmodels.stats.multitest import multipletests


# =========================================================
# Settings
# =========================================================
projection_file_candidates = ["combined_LOO.csv", "df_combined_LOO_new.csv"]
primary_bmi_trait = "BMI_Loh2018_UKBB_munged_txt"
duplicate_bmi_traits = {"BMI_Yengo2018", "BMI_UPDATE_Yengo2018"}
drop_duplicate_bmi = True

# Rank and plot settings
top_n_to_print = 15
top_n_heatmap = 20
selected_metric_for_barplot = "connectivity_sum"   # or one of: BST, CEA, LHA, PAG, PBN, PVH, PVT


# =========================================================
# Helpers
# =========================================================
def safe_spearman(x, y):
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = spearmanr(x, y)
        if hasattr(res, "statistic"):
            return float(res.statistic), float(res.pvalue)
        return float(res[0]), float(res[1])
    except Exception:
        return np.nan, np.nan


def safe_pearson(x, y):
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = pearsonr(x, y)
        if hasattr(res, "statistic"):
            return float(res.statistic), float(res.pvalue)
        return float(res[0]), float(res[1])
    except Exception:
        return np.nan, np.nan


def safe_hc3_ols(x, y):
    """
    Robust OLS: y ~ x
    Returns slope, p-value, R^2
    """
    try:
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)

        if np.all(np.isclose(x, x[0])):
            return np.nan, np.nan, np.nan

        X = sm.add_constant(x)
        fit = sm.OLS(y, X).fit(cov_type="HC3")
        slope = float(fit.params[1])
        pval = float(fit.pvalues[1])
        r2 = float(fit.rsquared)
        return slope, pval, r2
    except Exception:
        return np.nan, np.nan, np.nan


def get_metric_label(metric):
    if metric == "connectivity_sum":
        return "Overall connectivity"
    return metric


gwas_labels = {
    "BMI_Loh2018_UKBB_munged_txt": "Body Mass Index (BMI)",
    "BMI_Yengo2018": "Body Mass Index (BMI) - duplicate",
    "BMI_UPDATE_Yengo2018": "Body Mass Index (BMI) - duplicate",
    "LEPR_P48357_OID20128_v1_Cardiometabolic": "Leptin Receptor (LEPR)",
    "FLT4_P35916_OID21491_v1_Oncology": "VEGF Receptor 3 (FLT4)",
    "DiastolicBP_GCST90310295_munged_txt": "Diastolic Blood Pressure",
    "Pulse_GCST90310296_munged_txt": "Pulse Rate",
    "forced_expiratory_volume_GCST90244092_munged_txt": "Forced Expiratory Volume (FEV)",
}

def get_gwas_label(gwas_name):
    return gwas_labels.get(gwas_name, gwas_name)


# =========================================================
# 1) Load data
# =========================================================
df_cellect = pd.read_csv("df_cellect.csv")

loo_path = next((p for p in projection_file_candidates if Path(p).exists()), None)
if loo_path is None:
    raise FileNotFoundError(f"Could not find any of: {projection_file_candidates}")

combined_LOO = pd.read_csv(loo_path)

print("Using projection file:", loo_path)
print("Total GWAS traits in df_cellect:", df_cellect["gwas"].nunique())
print("Total cell types in df_cellect:", df_cellect["annotation"].nunique())


# =========================================================
# 2) Aggregate AAV ESmu by cell type x area
# =========================================================
proj_mean = (
    combined_LOO
    .groupby(["cluster", "area"], as_index=False)["esmu"]
    .mean()
)

proj_wide = (
    proj_mean
    .pivot(index="cluster", columns="area", values="esmu")
    .fillna(0)
)

# Area order from notebook when available
preferred_area_order = ["BST", "CEA", "LHA", "PAG", "PBN", "PVH", "PVT"]
area_cols = [a for a in preferred_area_order if a in proj_wide.columns]
area_cols += [a for a in proj_wide.columns if a not in area_cols]

proj_wide["connectivity_sum"] = proj_wide[area_cols].sum(axis=1)

proj_summary = proj_wide.reset_index()

print("\nProjection metrics available:")
print(["connectivity_sum"] + area_cols)
print("Overlapping projection cell types:", proj_summary["cluster"].nunique())


# =========================================================
# 3) Optionally remove duplicate BMI traits
# =========================================================
all_traits = list(df_cellect["gwas"].unique())

if drop_duplicate_bmi:
    traits_to_keep = [g for g in all_traits if g not in duplicate_bmi_traits]
else:
    traits_to_keep = all_traits

df_cellect_use = df_cellect[df_cellect["gwas"].isin(traits_to_keep)].copy()

print("\nTraits used in comparative analysis:", df_cellect_use["gwas"].nunique())
if drop_duplicate_bmi:
    print("Dropped duplicate BMI traits:", sorted(duplicate_bmi_traits))


# =========================================================
# 4) Comparative analysis across all traits
#    Primary outcome: GWAS z-score by cell type
#    Secondary: beta
# =========================================================
metrics = ["connectivity_sum"] + area_cols
results = []

for gwas in sorted(df_cellect_use["gwas"].unique()):
    gwas_dat = (
        df_cellect_use.loc[
            df_cellect_use["gwas"] == gwas,
            ["annotation", "z", "beta", "pvalue", "mlog10p"]
        ]
        .rename(columns={
            "annotation": "cluster",
            "z": "gwas_z",
            "beta": "gwas_beta",
            "pvalue": "gwas_p",
            "mlog10p": "gwas_mlog10p",
        })
        .copy()
    )

    dat = proj_summary.merge(gwas_dat, on="cluster", how="inner").copy()

    if dat.empty:
        continue

    for metric in metrics:
        x = dat[metric].astype(float).values
        y_z = dat["gwas_z"].astype(float).values
        y_beta = dat["gwas_beta"].astype(float).values

        sp_z_rho, sp_z_p = safe_spearman(x, y_z)
        pe_z_r, pe_z_p = safe_pearson(x, y_z)
        ols_z_beta, ols_z_p, ols_z_r2 = safe_hc3_ols(x, y_z)

        sp_beta_rho, sp_beta_p = safe_spearman(x, y_beta)
        pe_beta_r, pe_beta_p = safe_pearson(x, y_beta)
        ols_beta_coef, ols_beta_p, ols_beta_r2 = safe_hc3_ols(x, y_beta)

        results.append({
            "gwas": gwas,
            "gwas_label": get_gwas_label(gwas),
            "metric": metric,
            "metric_label": get_metric_label(metric),
            "n_celltypes": int(dat["cluster"].nunique()),
            "n_nonzero_metric": int((dat[metric] > 0).sum()),
            "mean_metric": float(dat[metric].mean()),
            "spearman_rho_z": sp_z_rho,
            "spearman_p_z": sp_z_p,
            "pearson_r_z": pe_z_r,
            "pearson_p_z": pe_z_p,
            "ols_beta_z": ols_z_beta,
            "ols_p_z": ols_z_p,
            "ols_r2_z": ols_z_r2,
            "spearman_rho_beta": sp_beta_rho,
            "spearman_p_beta": sp_beta_p,
            "pearson_r_beta": pe_beta_r,
            "pearson_p_beta": pe_beta_p,
            "ols_beta_beta": ols_beta_coef,
            "ols_p_beta": ols_beta_p,
            "ols_r2_beta": ols_beta_r2,
        })

results_df = pd.DataFrame(results)

if results_df.empty:
    raise ValueError("No comparative results were produced. Check that df_cellect and projection data overlap on cell types.")

print("\nComparative result rows:", len(results_df))
print("Metrics tested:", results_df["metric"].nunique())
print("Traits compared:", results_df["gwas"].nunique())


# =========================================================
# 5) Multiple-testing correction
#    A) within each connectivity metric across all traits
#    B) global across all trait x metric tests
# =========================================================
results_df["spearman_fdr_metric_z"] = np.nan
results_df["pearson_fdr_metric_z"] = np.nan
results_df["ols_fdr_metric_z"] = np.nan

results_df["spearman_fdr_metric_beta"] = np.nan
results_df["pearson_fdr_metric_beta"] = np.nan
results_df["ols_fdr_metric_beta"] = np.nan

for metric in metrics:
    idx = results_df["metric"] == metric

    results_df.loc[idx, "spearman_fdr_metric_z"] = multipletests(
        results_df.loc[idx, "spearman_p_z"].fillna(1.0),
        method="fdr_bh"
    )[1]

    results_df.loc[idx, "pearson_fdr_metric_z"] = multipletests(
        results_df.loc[idx, "pearson_p_z"].fillna(1.0),
        method="fdr_bh"
    )[1]

    results_df.loc[idx, "ols_fdr_metric_z"] = multipletests(
        results_df.loc[idx, "ols_p_z"].fillna(1.0),
        method="fdr_bh"
    )[1]

    results_df.loc[idx, "spearman_fdr_metric_beta"] = multipletests(
        results_df.loc[idx, "spearman_p_beta"].fillna(1.0),
        method="fdr_bh"
    )[1]

    results_df.loc[idx, "pearson_fdr_metric_beta"] = multipletests(
        results_df.loc[idx, "pearson_p_beta"].fillna(1.0),
        method="fdr_bh"
    )[1]

    results_df.loc[idx, "ols_fdr_metric_beta"] = multipletests(
        results_df.loc[idx, "ols_p_beta"].fillna(1.0),
        method="fdr_bh"
    )[1]

results_df["spearman_fdr_global_z"] = multipletests(
    results_df["spearman_p_z"].fillna(1.0),
    method="fdr_bh"
)[1]

results_df["pearson_fdr_global_z"] = multipletests(
    results_df["pearson_p_z"].fillna(1.0),
    method="fdr_bh"
)[1]

results_df["ols_fdr_global_z"] = multipletests(
    results_df["ols_p_z"].fillna(1.0),
    method="fdr_bh"
)[1]


# =========================================================
# 6) Add rankings
#    - rank by positive Spearman rho (descending)
#    - rank by absolute Spearman rho
# =========================================================
results_df["rank_spearman_z_within_metric"] = np.nan
results_df["rank_abs_spearman_z_within_metric"] = np.nan

for metric in metrics:
    idx = results_df["metric"] == metric

    metric_df = results_df.loc[idx].copy()

    metric_df = metric_df.sort_values(
        ["spearman_rho_z", "spearman_p_z"],
        ascending=[False, True]
    ).reset_index()

    metric_df["rank_spearman_z_within_metric"] = np.arange(1, len(metric_df) + 1)

    abs_df = results_df.loc[idx].copy()
    abs_df["abs_spearman_rho_z"] = abs_df["spearman_rho_z"].abs()
    abs_df = abs_df.sort_values(
        ["abs_spearman_rho_z", "spearman_p_z"],
        ascending=[False, True]
    ).reset_index()
    abs_df["rank_abs_spearman_z_within_metric"] = np.arange(1, len(abs_df) + 1)

    results_df.loc[metric_df["index"], "rank_spearman_z_within_metric"] = metric_df["rank_spearman_z_within_metric"].values
    results_df.loc[abs_df["index"], "rank_abs_spearman_z_within_metric"] = abs_df["rank_abs_spearman_z_within_metric"].values


# =========================================================
# 7) Print top traits for each connectivity metric
# =========================================================
cols_to_show = [
    "gwas_label",
    "gwas",
    "n_celltypes",
    "n_nonzero_metric",
    "spearman_rho_z",
    "spearman_p_z",
    "spearman_fdr_metric_z",
    "spearman_fdr_global_z",
    "pearson_r_z",
    "pearson_p_z",
    "ols_beta_z",
    "ols_p_z",
    "rank_spearman_z_within_metric",
]

print("\n" + "=" * 100)
print("TOP POSITIVE TRAITS FOR EACH CONNECTIVITY METRIC (ranked by Spearman rho with GWAS z-score)")
print("=" * 100)

for metric in metrics:
    print(f"\n--- {get_metric_label(metric)} ---")
    metric_df = (
        results_df[results_df["metric"] == metric]
        .sort_values(["spearman_rho_z", "spearman_p_z"], ascending=[False, True])
        .head(top_n_to_print)
    )
    print(metric_df[cols_to_show].to_string(index=False))


# =========================================================
# 8) BMI-specific comparison summary
# =========================================================
bmi_rows = results_df[results_df["gwas"] == primary_bmi_trait].copy()

if bmi_rows.empty:
    print(f"\nWARNING: {primary_bmi_trait} not found in results.")
else:
    bmi_summary = (
        bmi_rows[[
            "metric",
            "metric_label",
            "spearman_rho_z",
            "spearman_p_z",
            "spearman_fdr_metric_z",
            "spearman_fdr_global_z",
            "pearson_r_z",
            "pearson_p_z",
            "ols_beta_z",
            "ols_p_z",
            "rank_spearman_z_within_metric",
            "rank_abs_spearman_z_within_metric",
        ]]
        .sort_values("rank_spearman_z_within_metric")
        .reset_index(drop=True)
    )

    print("\n" + "=" * 100)
    print(f"BMI RANK SUMMARY: {get_gwas_label(primary_bmi_trait)}")
    print("=" * 100)
    print(bmi_summary.to_string(index=False))


# =========================================================
# 9) Does BMI stand out overall?
#    Summarize BMI relative to all traits by percentile rank
# =========================================================
if not bmi_rows.empty:
    bmi_compare = results_df.copy()

    metric_counts = bmi_compare.groupby("metric")["gwas"].nunique().to_dict()

    bmi_percentiles = []
    for _, row in bmi_rows.iterrows():
        n_traits_here = metric_counts[row["metric"]]
        rank_here = int(row["rank_spearman_z_within_metric"])
        percentile_here = 100 * (1 - (rank_here - 1) / n_traits_here)

        bmi_percentiles.append({
            "metric": row["metric"],
            "metric_label": row["metric_label"],
            "bmi_rank_positive": rank_here,
            "n_traits": n_traits_here,
            "bmi_percentile_positive": percentile_here,
            "bmi_spearman_rho_z": row["spearman_rho_z"],
            "bmi_spearman_p_z": row["spearman_p_z"],
            "bmi_spearman_fdr_metric_z": row["spearman_fdr_metric_z"],
        })

    bmi_percentiles_df = pd.DataFrame(bmi_percentiles).sort_values("bmi_rank_positive")
    print("\n" + "=" * 100)
    print("BMI POSITION RELATIVE TO ALL OTHER TRAITS")
    print("=" * 100)
    print(bmi_percentiles_df.to_string(index=False))


# =========================================================
# 10) Traits with the strongest connectivity relationship
#     anywhere across metrics
# =========================================================
best_any_metric = (
    results_df.sort_values(
        ["gwas", "spearman_rho_z", "spearman_p_z"],
        ascending=[True, False, True]
    )
    .groupby("gwas", as_index=False)
    .first()
)

best_any_metric["abs_spearman_rho_z"] = best_any_metric["spearman_rho_z"].abs()

print("\n" + "=" * 100)
print("TRAITS WITH THE STRONGEST BEST-CASE CONNECTIVITY RELATIONSHIP ACROSS ANY METRIC")
print("=" * 100)

print(
    best_any_metric.sort_values(
        ["spearman_rho_z", "spearman_p_z"],
        ascending=[False, True]
    )[[
        "gwas_label",
        "gwas",
        "metric_label",
        "spearman_rho_z",
        "spearman_p_z",
        "spearman_fdr_metric_z",
        "spearman_fdr_global_z",
        "pearson_r_z",
        "ols_beta_z",
        "ols_p_z",
    ]].head(top_n_to_print).to_string(index=False)
)


# =========================================================
# 11) Save full results
# =========================================================
results_df.to_csv("all_gwas_connectivity_comparative_results.csv", index=False)
print("\nSaved full results to: all_gwas_connectivity_comparative_results.csv")


# =========================================================
# 12) Plot A: bar plot of top traits for one selected metric
# =========================================================
plot_df = (
    results_df[results_df["metric"] == selected_metric_for_barplot]
    .sort_values(["spearman_rho_z", "spearman_p_z"], ascending=[False, True])
    .head(top_n_to_print)
    .copy()
)

if not plot_df.empty:
    labels = plot_df["gwas_label"].tolist()
    vals = plot_df["spearman_rho_z"].tolist()

    # mark BMI in label
    labels = [
        f"{lab}  [BMI]" if g == primary_bmi_trait else lab
        for lab, g in zip(labels, plot_df["gwas"])
    ]

    plt.figure(figsize=(10, max(5, 0.4 * len(plot_df))))
    y = np.arange(len(plot_df))
    plt.barh(y, vals)
    plt.yticks(y, labels)
    plt.gca().invert_yaxis()
    plt.xlabel("Spearman rho with GWAS z-score")
    plt.title(f"Top traits for {get_metric_label(selected_metric_for_barplot)}")
    plt.tight_layout()
    plt.show()


# =========================================================
# 13) Plot B: heatmap of top traits across metrics
#     top traits chosen by best positive rho in any metric
# =========================================================
top_traits_for_heatmap = (
    best_any_metric.sort_values(["spearman_rho_z", "spearman_p_z"], ascending=[False, True])
    .head(top_n_heatmap)["gwas"]
    .tolist()
)

heatmap_df = (
    results_df[results_df["gwas"].isin(top_traits_for_heatmap)]
    .pivot(index="gwas_label", columns="metric_label", values="spearman_rho_z")
)

# Preserve top-trait order
top_trait_labels = [
    get_gwas_label(g) for g in top_traits_for_heatmap
]
heatmap_df = heatmap_df.reindex(top_trait_labels)

# Preserve metric order
metric_label_order = [get_metric_label(m) for m in metrics]
heatmap_df = heatmap_df[metric_label_order]

fig, ax = plt.subplots(figsize=(2 + 1.2 * len(metric_label_order), 0.45 * len(heatmap_df) + 2))
im = ax.imshow(heatmap_df.values, aspect="auto", interpolation="nearest")

ax.set_xticks(np.arange(len(heatmap_df.columns)))
ax.set_xticklabels(heatmap_df.columns, rotation=45, ha="right")
ax.set_yticks(np.arange(len(heatmap_df.index)))
ax.set_yticklabels(heatmap_df.index)

# annotate BMI row
for i, trait_name in enumerate(heatmap_df.index):
    if trait_name == get_gwas_label(primary_bmi_trait):
        ax.get_yticklabels()[i].set_fontweight("bold")

cbar = plt.colorbar(im, ax=ax)
cbar.set_label("Spearman rho (GWAS z vs connectivity)")

ax.set_title("Top connectivity-associated GWAS traits across metrics")
plt.tight_layout()
plt.show()


# =========================================================
# 14) Optional: direct BMI vs all traits table for overall connectivity
# =========================================================
overall_metric = "connectivity_sum"
overall_table = (
    results_df[results_df["metric"] == overall_metric]
    .sort_values(["spearman_rho_z", "spearman_p_z"], ascending=[False, True])
    .reset_index(drop=True)
)

print("\n" + "=" * 100)
print("FULL RANKING FOR OVERALL CONNECTIVITY")
print("=" * 100)
print(
    overall_table[[
        "gwas_label",
        "gwas",
        "spearman_rho_z",
        "spearman_p_z",
        "spearman_fdr_metric_z",
        "pearson_r_z",
        "ols_beta_z",
        "ols_p_z",
        "rank_spearman_z_within_metric",
    ]].to_string(index=False)
)


# In[6]:


import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests


# =========================================================
# Settings
# =========================================================
projection_file_candidates = ["combined_LOO.csv", "df_combined_LOO_new.csv"]

# Primary BMI trait to highlight
primary_bmi_trait = "BMI_Loh2018_UKBB_munged_txt"

# Duplicate BMI traits to optionally drop from the comparative panel
duplicate_bmi_traits = {"BMI_Yengo2018", "BMI_UPDATE_Yengo2018"}
drop_duplicate_bmi = True

# Which outcome to treat as primary
# Use z-score as primary for cross-trait comparisons because raw betas are not
# directly comparable across different GWAS traits.
primary_outcome = "gwas_z"      # options: "gwas_z" or "gwas_beta"
secondary_outcome = "gwas_beta" # fitted as a supplement

top_n_to_print = 15
selected_metric_for_plot = "PVH"   # e.g. "connectivity_sum", "PVH", "CEA", ...


# =========================================================
# Labels (optional nicer names)
# =========================================================
gwas_labels = {
    "BMI_Loh2018_UKBB_munged_txt": "Body Mass Index (BMI)",
    "BMI_Yengo2018": "Body Mass Index (BMI) - duplicate",
    "BMI_UPDATE_Yengo2018": "Body Mass Index (BMI) - duplicate",
    "LEPR_P48357_OID20128_v1_Cardiometabolic": "Leptin Receptor (LEPR)",
    "FLT4_P35916_OID21491_v1_Oncology": "VEGF Receptor 3 (FLT4)",
    "DiastolicBP_GCST90310295_munged_txt": "Diastolic Blood Pressure",
    "Pulse_GCST90310296_munged_txt": "Pulse Rate",
    "forced_expiratory_volume_GCST90244092_munged_txt": "Forced Expiratory Volume (FEV)",
}

def get_gwas_label(gwas_name):
    return gwas_labels.get(gwas_name, gwas_name)

def get_metric_label(metric):
    if metric == "connectivity_sum":
        return "Overall connectivity"
    return metric


# =========================================================
# 1) Load data
# =========================================================
df_cellect = pd.read_csv("df_cellect.csv")

loo_path = next((p for p in projection_file_candidates if Path(p).exists()), None)
if loo_path is None:
    raise FileNotFoundError(f"Could not find any of: {projection_file_candidates}")

combined_LOO = pd.read_csv(loo_path)

print("Using projection file:", loo_path)
print("Total GWAS traits in df_cellect:", df_cellect["gwas"].nunique())
print("Total cell types in df_cellect:", df_cellect["annotation"].nunique())


# =========================================================
# 2) Aggregate connectivity by cell type x area
# =========================================================
proj_mean = (
    combined_LOO
    .groupby(["cluster", "area"], as_index=False)["esmu"]
    .mean()
)

proj_wide = (
    proj_mean
    .pivot(index="cluster", columns="area", values="esmu")
    .fillna(0)
)

preferred_area_order = ["BST", "CEA", "LHA", "PAG", "PBN", "PVH", "PVT"]
area_cols = [a for a in preferred_area_order if a in proj_wide.columns]
area_cols += [a for a in proj_wide.columns if a not in area_cols]

proj_wide["connectivity_sum"] = proj_wide[area_cols].sum(axis=1)

proj_summary = proj_wide.reset_index()

print("\nProjection metrics available:")
print(["connectivity_sum"] + area_cols)
print("Overlapping projection cell types:", proj_summary["cluster"].nunique())


# =========================================================
# 3) Optionally remove duplicate BMI traits
# =========================================================
all_traits = list(df_cellect["gwas"].unique())

if drop_duplicate_bmi:
    traits_to_keep = [g for g in all_traits if g not in duplicate_bmi_traits]
else:
    traits_to_keep = all_traits

df_cellect_use = df_cellect[df_cellect["gwas"].isin(traits_to_keep)].copy()

print("\nTraits used in comparative analysis:", df_cellect_use["gwas"].nunique())
if drop_duplicate_bmi:
    print("Dropped duplicate BMI traits:", sorted(duplicate_bmi_traits))


# =========================================================
# 4) Merge to long format: one row per cell type x trait
# =========================================================
long_df = (
    df_cellect_use.rename(columns={
        "annotation": "cluster",
        "z": "gwas_z",
        "beta": "gwas_beta",
        "pvalue": "gwas_p",
        "mlog10p": "gwas_mlog10p",
    })
    .merge(proj_summary, on="cluster", how="inner")
    .copy()
)

print("\nLong table shape (cell type x trait):", long_df.shape)
print("Traits in long table:", long_df["gwas"].nunique())
print("Cell types in long table:", long_df["cluster"].nunique())

metrics = ["connectivity_sum"] + area_cols


# =========================================================
# 5) Standardize each connectivity metric across cell types
#    so slopes are comparable across metrics
# =========================================================
for metric in metrics:
    mean_val = long_df[metric].mean()
    std_val = long_df[metric].std(ddof=0)

    if np.isclose(std_val, 0):
        raise ValueError(f"Metric {metric} has near-zero variance and cannot be standardized.")

    long_df[f"{metric}_std"] = (long_df[metric] - mean_val) / std_val


# =========================================================
# 6) Helper to fit cell-type-fixed-effect model
#
# Model:
#   outcome ~ 0 + C(cluster) + C(gwas) + 0 + C(gwas):metric_std
#
# Interpretation:
# - C(cluster): keeps track of cell type explicitly (cell-type fixed effects)
# - C(gwas): trait intercepts
# - C(gwas):metric_std gives a separate connectivity slope for each trait
#
# Because connectivity is constant within cell type, its main effect is absorbed
# by C(cluster), so we estimate trait-specific slopes rather than a global slope.
# =========================================================
def fit_trait_specific_fe_model(df, outcome, metric_std_col):
    formula = f"{outcome} ~ 0 + C(cluster) + C(gwas) + 0 + C(gwas):{metric_std_col}"

    try:
        fit = smf.ols(formula, data=df).fit(
            cov_type="cluster",
            cov_kwds={"groups": df["cluster"]}
        )
        cov_label = "cluster-robust by cell type"
    except Exception:
        fit = smf.ols(formula, data=df).fit(cov_type="HC3")
        cov_label = "HC3"

    return fit, cov_label


def extract_trait_slopes(fit, trait_list, metric_std_col):
    """
    Extract one trait-specific slope and p-value for every trait.
    """
    out = []
    params = fit.params
    pvals = fit.pvalues
    conf = fit.conf_int()

    for g in trait_list:
        candidate_names = [
            f"C(gwas)[{g}]:{metric_std_col}",
            f"C(gwas)[T.{g}]:{metric_std_col}",
        ]

        found_name = None
        for nm in candidate_names:
            if nm in params.index:
                found_name = nm
                break

        # fallback regex match if needed
        if found_name is None:
            pattern = re.compile(rf"^C\(gwas\)\[(T\.)?{re.escape(g)}\]:{re.escape(metric_std_col)}$")
            matches = [nm for nm in params.index if pattern.match(nm)]
            if len(matches) == 1:
                found_name = matches[0]

        if found_name is None:
            out.append({
                "gwas": g,
                "slope": np.nan,
                "pvalue": np.nan,
                "ci_low": np.nan,
                "ci_high": np.nan,
            })
        else:
            out.append({
                "gwas": g,
                "slope": float(params[found_name]),
                "pvalue": float(pvals[found_name]),
                "ci_low": float(conf.loc[found_name, 0]),
                "ci_high": float(conf.loc[found_name, 1]),
            })

    return pd.DataFrame(out)


# =========================================================
# 7) Fit models for every connectivity metric
#    Primary: GWAS z
#    Secondary: GWAS beta
# =========================================================
traits = sorted(long_df["gwas"].unique())
all_results = []

for metric in metrics:
    metric_std_col = f"{metric}_std"

    # Primary model on GWAS z
    fit_z, cov_label_z = fit_trait_specific_fe_model(long_df, "gwas_z", metric_std_col)
    slopes_z = extract_trait_slopes(fit_z, traits, metric_std_col)
    slopes_z = slopes_z.rename(columns={
        "slope": "fe_slope_z",
        "pvalue": "fe_p_z",
        "ci_low": "fe_ci_low_z",
        "ci_high": "fe_ci_high_z",
    })

    # Secondary model on raw GWAS beta
    fit_beta, cov_label_beta = fit_trait_specific_fe_model(long_df, "gwas_beta", metric_std_col)
    slopes_beta = extract_trait_slopes(fit_beta, traits, metric_std_col)
    slopes_beta = slopes_beta.rename(columns={
        "slope": "fe_slope_beta",
        "pvalue": "fe_p_beta",
        "ci_low": "fe_ci_low_beta",
        "ci_high": "fe_ci_high_beta",
    })

    res = slopes_z.merge(slopes_beta, on="gwas", how="outer")
    res["gwas_label"] = res["gwas"].map(get_gwas_label)
    res["metric"] = metric
    res["metric_label"] = get_metric_label(metric)
    res["n_celltypes"] = long_df["cluster"].nunique()
    res["covariance_z"] = cov_label_z
    res["covariance_beta"] = cov_label_beta

    # FDR within this metric across traits
    res["fe_fdr_metric_z"] = multipletests(res["fe_p_z"].fillna(1.0), method="fdr_bh")[1]
    res["fe_fdr_metric_beta"] = multipletests(res["fe_p_beta"].fillna(1.0), method="fdr_bh")[1]

    # ranks within metric
    res["rank_fe_slope_z_within_metric"] = (
        res["fe_slope_z"]
        .rank(method="min", ascending=False)
    )
    res["rank_abs_fe_slope_z_within_metric"] = (
        res["fe_slope_z"].abs()
        .rank(method="min", ascending=False)
    )

    all_results.append(res)

results_df = pd.concat(all_results, ignore_index=True)

# Global FDR across all trait x metric tests
results_df["fe_fdr_global_z"] = multipletests(
    results_df["fe_p_z"].fillna(1.0),
    method="fdr_bh"
)[1]

results_df["fe_fdr_global_beta"] = multipletests(
    results_df["fe_p_beta"].fillna(1.0),
    method="fdr_bh"
)[1]


# =========================================================
# 8) Print top traits for each metric
# =========================================================
print("\n" + "=" * 110)
print("CELL-TYPE-FIXED-EFFECT ANALYSIS")
print("Trait-specific connectivity slopes after accounting for cell-type identity")
print("=" * 110)

cols_show = [
    "gwas_label",
    "gwas",
    "fe_slope_z",
    "fe_p_z",
    "fe_fdr_metric_z",
    "fe_fdr_global_z",
    "fe_slope_beta",
    "fe_p_beta",
    "rank_fe_slope_z_within_metric",
]

for metric in metrics:
    print(f"\n--- {get_metric_label(metric)} ---")
    tmp = (
        results_df[results_df["metric"] == metric]
        .sort_values(["fe_slope_z", "fe_p_z"], ascending=[False, True])
        .head(top_n_to_print)
    )
    print(tmp[cols_show].to_string(index=False))


# =========================================================
# 9) BMI summary
# =========================================================
bmi_rows = results_df[results_df["gwas"] == primary_bmi_trait].copy()

if bmi_rows.empty:
    print(f"\nWARNING: {primary_bmi_trait} not found in fixed-effect results.")
else:
    bmi_summary = bmi_rows[[
        "metric",
        "metric_label",
        "fe_slope_z",
        "fe_p_z",
        "fe_fdr_metric_z",
        "fe_fdr_global_z",
        "fe_slope_beta",
        "fe_p_beta",
        "rank_fe_slope_z_within_metric",
        "rank_abs_fe_slope_z_within_metric",
    ]].sort_values("rank_fe_slope_z_within_metric")

    print("\n" + "=" * 110)
    print(f"BMI SUMMARY IN CELL-TYPE-FIXED-EFFECT MODEL: {get_gwas_label(primary_bmi_trait)}")
    print("=" * 110)
    print(bmi_summary.to_string(index=False))


# =========================================================
# 10) BMI percentile position within each metric
# =========================================================
if not bmi_rows.empty:
    metric_counts = results_df.groupby("metric")["gwas"].nunique().to_dict()

    bmi_percentiles = []
    for _, row in bmi_rows.iterrows():
        n_traits_here = metric_counts[row["metric"]]
        rank_here = int(row["rank_fe_slope_z_within_metric"])
        percentile_here = 100 * (1 - (rank_here - 1) / n_traits_here)

        bmi_percentiles.append({
            "metric": row["metric"],
            "metric_label": row["metric_label"],
            "bmi_rank_positive": rank_here,
            "n_traits": n_traits_here,
            "bmi_percentile_positive": percentile_here,
            "bmi_fe_slope_z": row["fe_slope_z"],
            "bmi_fe_p_z": row["fe_p_z"],
            "bmi_fe_fdr_metric_z": row["fe_fdr_metric_z"],
            "bmi_fe_fdr_global_z": row["fe_fdr_global_z"],
        })

    bmi_percentiles_df = pd.DataFrame(bmi_percentiles).sort_values("bmi_rank_positive")

    print("\n" + "=" * 110)
    print("BMI POSITION RELATIVE TO ALL OTHER TRAITS (cell-type-fixed-effect model)")
    print("=" * 110)
    print(bmi_percentiles_df.to_string(index=False))


# =========================================================
# 11) Save full results
# =========================================================
results_df.to_csv("all_gwas_connectivity_celltype_FE_results.csv", index=False)
print("\nSaved full results to: all_gwas_connectivity_celltype_FE_results.csv")


# =========================================================
# 12) Optional plot for one selected metric
# =========================================================
plot_df = (
    results_df[results_df["metric"] == selected_metric_for_plot]
    .sort_values(["fe_slope_z", "fe_p_z"], ascending=[False, True])
    .head(top_n_to_print)
    .copy()
)

if not plot_df.empty:
    labels = [
        f"{lab}  [BMI]" if g == primary_bmi_trait else lab
        for lab, g in zip(plot_df["gwas_label"], plot_df["gwas"])
    ]

    plt.figure(figsize=(10, max(5, 0.4 * len(plot_df))))
    y = np.arange(len(plot_df))
    plt.barh(y, plot_df["fe_slope_z"].values)
    plt.yticks(y, labels)
    plt.gca().invert_yaxis()
    plt.xlabel("Trait-specific connectivity slope (GWAS z, cell-type FE model)")
    plt.title(f"Top traits for {get_metric_label(selected_metric_for_plot)}\ncell-type-fixed-effect model")
    plt.tight_layout()
    plt.show()


# =========================================================
# 13) Optional: inspect raw BMI cell-type table
#     This is useful for seeing which cell types drive the BMI slope.
# =========================================================
if primary_bmi_trait in long_df["gwas"].unique():
    bmi_raw = (
        long_df[long_df["gwas"] == primary_bmi_trait][
            ["cluster", "gwas_z", "gwas_beta", "gwas_mlog10p"] + metrics
        ]
        .sort_values("gwas_z", ascending=False)
        .reset_index(drop=True)
    )

    print("\n" + "=" * 110)
    print("RAW BMI CELL-TYPE TABLE")
    print("=" * 110)
    print(bmi_raw.to_string(index=False))

    bmi_raw.to_csv("BMI_celltype_connectivity_table.csv", index=False)
    print("\nSaved BMI cell-type table to: BMI_celltype_connectivity_table.csv")


# In[7]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

# =========================================================
# Settings
# =========================================================
bmi_gwas = "BMI_Loh2018_UKBB_munged_txt"
projection_file_candidates = ["combined_LOO.csv", "df_combined_LOO_new.csv"]

# Define "top BMI cell types"
top_k_bmi = 10

# Define "high-projecting" within each area:
# use top tertile among NONZERO-projecting cell types
high_proj_quantile = 0.67

# Preferred area order
preferred_area_order = ["BST", "CEA", "LHA", "PAG", "PBN", "PVH", "PVT"]


# =========================================================
# 1) Load data
# =========================================================
df_cellect = pd.read_csv("df_cellect.csv")

loo_path = next((p for p in projection_file_candidates if Path(p).exists()), None)
if loo_path is None:
    raise FileNotFoundError(f"Could not find any of: {projection_file_candidates}")

combined_LOO = pd.read_csv(loo_path)

print("Using projection file:", loo_path)


# =========================================================
# 2) Aggregate ESmu by cluster x area
# =========================================================
proj_mean = (
    combined_LOO
    .groupby(["cluster", "area"], as_index=False)["esmu"]
    .mean()
)

proj_wide = (
    proj_mean
    .pivot(index="cluster", columns="area", values="esmu")
    .fillna(0)
)

area_cols = [a for a in preferred_area_order if a in proj_wide.columns]
area_cols += [a for a in proj_wide.columns if a not in area_cols]

proj_summary = proj_wide.reset_index()


# =========================================================
# 3) Pull BMI GWAS per-cell-type values
# =========================================================
bmi = (
    df_cellect.loc[
        df_cellect["gwas"] == bmi_gwas,
        ["annotation", "z", "beta", "pvalue", "mlog10p"]
    ]
    .rename(columns={
        "annotation": "cluster",
        "z": "bmi_z",
        "beta": "bmi_beta",
        "pvalue": "bmi_p",
        "mlog10p": "bmi_mlog10p",
    })
    .copy()
)

dat = proj_summary.merge(bmi, on="cluster", how="inner").copy()
print("Overlapping cell types:", dat["cluster"].nunique())


# =========================================================
# 4) Define BMI-enriched cell-type sets
# =========================================================
bonf = 0.05 / dat["cluster"].nunique()
dat["bmi_sig_bonf"] = dat["bmi_p"] < bonf

dat = dat.sort_values("bmi_z", ascending=False).reset_index(drop=True)
dat["bmi_topk"] = False
dat.loc[:top_k_bmi - 1, "bmi_topk"] = True

print("\nBMI set sizes:")
print("Bonferroni-significant BMI cell types:", int(dat["bmi_sig_bonf"].sum()))
print(f"Top {top_k_bmi} BMI cell types:", int(dat["bmi_topk"].sum()))

print("\nTop BMI cell types:")
print(dat.loc[dat["bmi_topk"], ["cluster", "bmi_z", "bmi_p"]].to_string(index=False))


# =========================================================
# 5) Helper: define high-projecting cell types within an area
# =========================================================
def get_high_projecting_mask(series, quantile=0.67):
    """
    High-projecting = top quantile among NONZERO-projecting cell types.
    Returns:
      high_mask, cutoff, n_nonzero
    """
    nonzero = series[series > 0]

    if len(nonzero) == 0:
        return pd.Series(False, index=series.index), np.nan, 0

    cutoff = nonzero.quantile(quantile)
    high_mask = (series > 0) & (series >= cutoff)

    return high_mask, float(cutoff), int((series > 0).sum())


# =========================================================
# 6) Enrichment test for each area
#    Fisher exact test, alternative='greater'
# =========================================================
results = []
overlap_tables = {}

bmi_sets = {
    f"BMI_top{top_k_bmi}": "bmi_topk",
    "BMI_bonf": "bmi_sig_bonf",
}

for set_name, set_col in bmi_sets.items():
    bmi_mask = dat[set_col].astype(bool)
    n_bmi = int(bmi_mask.sum())

    if n_bmi == 0:
        continue

    for area in area_cols:
        high_mask, cutoff, n_nonzero = get_high_projecting_mask(
            dat[area],
            quantile=high_proj_quantile
        )

        n_high = int(high_mask.sum())
        N = len(dat)

        # Contingency table:
        #                high   not_high
        # BMI set          a       b
        # not BMI          c       d
        a = int((bmi_mask & high_mask).sum())
        b = int((bmi_mask & ~high_mask).sum())
        c = int((~bmi_mask & high_mask).sum())
        d = int((~bmi_mask & ~high_mask).sum())

        oddsratio, pval = fisher_exact([[a, b], [c, d]], alternative="greater")

        # expected overlap under random sampling
        expected_overlap = n_bmi * (n_high / N) if N > 0 else np.nan

        # fold enrichment
        if n_bmi > 0 and n_high > 0:
            frac_high_in_bmi = a / n_bmi
            frac_high_overall = n_high / N
            fold_enrichment = frac_high_in_bmi / frac_high_overall if frac_high_overall > 0 else np.nan
        else:
            frac_high_in_bmi = np.nan
            frac_high_overall = np.nan
            fold_enrichment = np.nan

        overlap_clusters = dat.loc[bmi_mask & high_mask, "cluster"].tolist()

        results.append({
            "bmi_set": set_name,
            "area": area,
            "n_total": N,
            "n_bmi_set": n_bmi,
            "n_nonzero_area": n_nonzero,
            "n_high_area": n_high,
            "high_cutoff_esmu": cutoff,
            "overlap_n": a,
            "expected_overlap": expected_overlap,
            "fold_enrichment": fold_enrichment,
            "frac_high_in_bmi_set": frac_high_in_bmi,
            "frac_high_overall": frac_high_overall,
            "odds_ratio": oddsratio,
            "fisher_p": pval,
            "overlap_celltypes": ", ".join(overlap_clusters) if overlap_clusters else "",
        })

        overlap_tables[(set_name, area)] = dat.loc[
            bmi_mask | high_mask,
            ["cluster", "bmi_z", "bmi_p", area, set_col]
        ].copy()

results_df = pd.DataFrame(results)

# FDR within each BMI set across 7 areas
results_df["fdr_by_set"] = np.nan
for set_name in results_df["bmi_set"].unique():
    idx = results_df["bmi_set"] == set_name
    results_df.loc[idx, "fdr_by_set"] = multipletests(
        results_df.loc[idx, "fisher_p"],
        method="fdr_bh"
    )[1]

# Sort
results_df = results_df.sort_values(
    ["bmi_set", "fisher_p", "fold_enrichment"],
    ascending=[True, True, False]
).reset_index(drop=True)


# =========================================================
# 7) Print results
# =========================================================
print("\n" + "=" * 100)
print("ENRICHMENT OF BMI CELL TYPES AMONG HIGH-PROJECTING CELL TYPES")
print("=" * 100)

for set_name in results_df["bmi_set"].unique():
    print(f"\n--- {set_name} ---")
    tmp = results_df[results_df["bmi_set"] == set_name].copy()
    print(
        tmp[[
            "area",
            "n_bmi_set",
            "n_high_area",
            "overlap_n",
            "expected_overlap",
            "fold_enrichment",
            "odds_ratio",
            "fisher_p",
            "fdr_by_set",
            "overlap_celltypes",
        ]].to_string(index=False)
    )


# =========================================================
# 8) Save results
# =========================================================
results_df.to_csv("BMI_high_projection_enrichment_by_area.csv", index=False)
print("\nSaved enrichment table to: BMI_high_projection_enrichment_by_area.csv")


# =========================================================
# 9) Plot fold-enrichment
# =========================================================
fig, axes = plt.subplots(1, len(results_df["bmi_set"].unique()), figsize=(6 * len(results_df["bmi_set"].unique()), 4), sharey=True)

if len(results_df["bmi_set"].unique()) == 1:
    axes = [axes]

for ax, set_name in zip(axes, results_df["bmi_set"].unique()):
    tmp = results_df[results_df["bmi_set"] == set_name].copy()

    # preserve area order
    tmp["area"] = pd.Categorical(tmp["area"], categories=area_cols, ordered=True)
    tmp = tmp.sort_values("area")

    ax.bar(tmp["area"], tmp["fold_enrichment"])
    ax.axhline(1.0, linestyle="--", linewidth=1)

    for i, row in enumerate(tmp.itertuples(index=False)):
        ax.text(i, row.fold_enrichment, f"p={row.fisher_p:.3g}", ha="center", va="bottom", fontsize=8)

    ax.set_title(set_name)
    ax.set_xlabel("Projection target")
    ax.set_ylabel("Fold-enrichment of BMI set in high-projecting cells")

plt.tight_layout()
plt.show()


# =========================================================
# 10) Optional: inspect exact overlaps area by area
# =========================================================
print("\n" + "=" * 100)
print("EXACT OVERLAPS")
print("=" * 100)

for set_name in results_df["bmi_set"].unique():
    for area in area_cols:
        row = results_df[(results_df["bmi_set"] == set_name) & (results_df["area"] == area)].iloc[0]
        print(f"\n[{set_name}] {area}")
        print(
            f"overlap={row['overlap_n']} | expected={row['expected_overlap']:.2f} | "
            f"fold={row['fold_enrichment']:.2f} | p={row['fisher_p']:.3g} | fdr={row['fdr_by_set']:.3g}"
        )
        if row["overlap_celltypes"]:
            print("Overlapping cell types:", row["overlap_celltypes"])
        else:
            print("Overlapping cell types: none")


# In[10]:


import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests


# =========================================================
# Settings
# =========================================================
df_cellect_path = "df_cellect.csv"
proj_candidates = ["combined_LOO.csv", "df_combined_LOO_new.csv"]

area_order = ["BST", "CEA", "LHA", "PAG", "PBN", "PVH", "PVT"]

# Need at least 2 BH-significant cell types for a GWAS to estimate a slope
min_sig_celltypes_per_gwas = 2

# BH mode:
#   "within_gwas" = adjust p-values within each GWAS across cell types
#   "global"      = adjust across all celltype x GWAS rows at once
bh_mode = "within_gwas"

# Output file names
out_results = "mixedmodel_twostage_BHsigOnly_results.csv"
out_heatmap = "gwas_area_heatmap_twostage_BHsigOnly.pdf"
out_compare = "heatmap_comparison_original_vs_BHsigOnly.pdf"
out_pval_compare = "pvalue_distribution_comparison_BHsigOnly.pdf"


# =========================================================
# Helpers
# =========================================================
def pick_first_existing(columns, candidates, label):
    for c in candidates:
        if c in columns:
            return c
    raise ValueError(
        f"Could not find {label}. Tried: {candidates}\n"
        f"Available columns: {list(columns)}"
    )

def fit_model_with_fallback(tmp, esmu_col, z_col, animal_col):
    """
    Primary: mixed model with random intercept for animal
    Fallback: OLS with cluster-robust SE by animal
    """
    formula = f"{esmu_col} ~ {z_col}"

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            model = smf.mixedlm(formula, data=tmp, groups=tmp[animal_col])
            fit = model.fit(reml=False, method="lbfgs", maxiter=500, disp=False)

        beta = float(fit.params[z_col])
        pval = float(fit.pvalues[z_col])
        converged = bool(getattr(fit, "converged", True))
        method = "MixedLM"
        return beta, pval, converged, method

    except Exception as e:
        fit = smf.ols(formula, data=tmp).fit(
            cov_type="cluster",
            cov_kwds={"groups": tmp[animal_col]}
        )
        beta = float(fit.params[z_col])
        pval = float(fit.pvalues[z_col])
        converged = np.nan
        method = f"OLS_clustered_fallback[{type(e).__name__}]"
        return beta, pval, converged, method

def add_bh_significance(df, gwas_col, pval_col, mode="within_gwas"):
    """
    Add BH-adjusted p-values and a boolean BH-significant mask.

    mode="within_gwas": BH correction separately within each GWAS across cell types
    mode="global":      BH correction across all rows together
    """
    df = df.copy()

    if mode == "within_gwas":
        df["pvalue_bh"] = np.nan
        for g, idx in df.groupby(gwas_col).groups.items():
            p = df.loc[idx, pval_col].astype(float).values
            q = multipletests(p, method="fdr_bh")[1]
            df.loc[idx, "pvalue_bh"] = q

    elif mode == "global":
        df["pvalue_bh"] = multipletests(
            df[pval_col].astype(float).values,
            method="fdr_bh"
        )[1]

    else:
        raise ValueError("bh_mode must be one of: 'within_gwas', 'global'")

    df["significant_bh"] = df["pvalue_bh"] < 0.05
    return df


# =========================================================
# 1) Load inputs
# =========================================================
df_cellect = pd.read_csv(df_cellect_path)

proj_path = next((p for p in proj_candidates if Path(p).exists()), None)
if proj_path is None:
    raise FileNotFoundError(f"Could not find any of: {proj_candidates}")

df_proj = pd.read_csv(proj_path)

print("Using projection file:", proj_path)
print("df_cellect shape:", df_cellect.shape)
print("df_proj shape:", df_proj.shape)

print("\ndf_cellect columns:")
print(df_cellect.columns.tolist())

print("\ndf_proj columns:")
print(df_proj.columns.tolist())


# =========================================================
# 2) Resolve columns
# =========================================================
cluster_cellect_col = pick_first_existing(
    df_cellect.columns,
    ["annotation", "cluster", "celltype", "detailed_celltype_annotation"],
    "cell-type column in df_cellect"
)

gwas_col = pick_first_existing(
    df_cellect.columns,
    ["gwas"],
    "GWAS column"
)

z_col = pick_first_existing(
    df_cellect.columns,
    ["z", "gwas_z"],
    "GWAS z-score column"
)

pval_col = pick_first_existing(
    df_cellect.columns,
    ["pvalue"],
    "p-value column in df_cellect"
)

cluster_proj_col = pick_first_existing(
    df_proj.columns,
    ["cluster", "annotation", "celltype", "detailed_celltype_annotation"],
    "cell-type column in projection table"
)

area_col = pick_first_existing(
    df_proj.columns,
    ["area"],
    "area column in projection table"
)

animal_col = pick_first_existing(
    df_proj.columns,
    ["left_out_animal", "animal", "hash.mcl.ID", "orig.ident", "mouse", "sample"],
    "animal column in projection table"
)

esmu_col = pick_first_existing(
    df_proj.columns,
    ["esmu", "ESmu", "ESMU"],
    "ESmu column in projection table"
)

print("\nResolved columns")
print("-" * 60)
print("cluster_cellect_col:", cluster_cellect_col)
print("gwas_col:", gwas_col)
print("z_col:", z_col)
print("pval_col:", pval_col)
print("cluster_proj_col:", cluster_proj_col)
print("area_col:", area_col)
print("animal_col:", animal_col)
print("esmu_col:", esmu_col)


# =========================================================
# 3) Compute BH significance from p-values
# =========================================================
df_cellect = add_bh_significance(
    df=df_cellect,
    gwas_col=gwas_col,
    pval_col=pval_col,
    mode=bh_mode
)

print("\nBH computation")
print("-" * 60)
print("BH mode:", bh_mode)
print("Total BH-significant celltype x GWAS rows:", int(df_cellect["significant_bh"].sum()))

sig_per_gwas = (
    df_cellect.groupby(gwas_col)["significant_bh"]
    .sum()
    .sort_values(ascending=False)
)

print("\nTop GWAS by number of BH-significant cell types:")
print(sig_per_gwas.head(15).to_string())


# =========================================================
# 4) Keep only BH-significant celltype x GWAS rows
# =========================================================
df_sig = df_cellect.loc[
    df_cellect["significant_bh"],
    [cluster_cellect_col, gwas_col, z_col, pval_col, "pvalue_bh"]
].copy()

df_sig = df_sig.dropna(subset=[cluster_cellect_col, gwas_col, z_col])

# Restrict to cell types present in projection table
proj_clusters = set(df_proj[cluster_proj_col].dropna().unique())
df_sig = df_sig[df_sig[cluster_cellect_col].isin(proj_clusters)].copy()

sig_counts = (
    df_sig.groupby(gwas_col)[cluster_cellect_col]
    .nunique()
    .sort_values(ascending=False)
)

gwas_keep = sig_counts[sig_counts >= min_sig_celltypes_per_gwas].index.tolist()
df_sig = df_sig[df_sig[gwas_col].isin(gwas_keep)].copy()

print("\nBH filtering summary")
print("-" * 60)
print("BH-significant projected celltype x GWAS rows kept:", len(df_sig))
print("GWAS with >= 2 BH-significant projected cell types:", len(gwas_keep))
print("Top GWAS by # BH-significant projected cell types:")
print(sig_counts.head(15).to_string())


# =========================================================
# 5) Re-fit area models using ONLY BH-significant cell types
# =========================================================
areas_present = [a for a in area_order if a in df_proj[area_col].unique()]
areas_present += [a for a in df_proj[area_col].unique() if a not in areas_present]

results = []

for gwas in gwas_keep:
    sub_sig = (
        df_sig[df_sig[gwas_col] == gwas][[cluster_cellect_col, z_col]]
        .drop_duplicates()
        .copy()
    )

    n_sig_ct_total = sub_sig[cluster_cellect_col].nunique()

    for area in areas_present:
        proj_area = df_proj[df_proj[area_col] == area].copy()

        tmp = proj_area.merge(
            sub_sig,
            left_on=cluster_proj_col,
            right_on=cluster_cellect_col,
            how="inner"
        )

        n_rows = len(tmp)
        n_animals = tmp[animal_col].nunique()
        n_sig_ct_in_area = tmp[cluster_proj_col].nunique()
        n_unique_z = tmp[z_col].nunique()

        # Need at least 2 unique z values and at least 2 animals
        if n_sig_ct_in_area < 2 or n_animals < 2 or n_unique_z < 2:
            results.append({
                "gwas": gwas,
                "area": area,
                "beta_z": np.nan,
                "pvalue": np.nan,
                "pvalue_fdr": np.nan,
                "significant_fdr": False,
                "n_sig_celltypes_total": n_sig_ct_total,
                "n_sig_celltypes_in_area": n_sig_ct_in_area,
                "n_animals": n_animals,
                "n_rows": n_rows,
                "model": "skipped_insufficient_data",
                "converged": np.nan,
            })
            continue

        beta, pval, converged, method = fit_model_with_fallback(
            tmp=tmp,
            esmu_col=esmu_col,
            z_col=z_col,
            animal_col=animal_col
        )

        results.append({
            "gwas": gwas,
            "area": area,
            "beta_z": beta,
            "pvalue": pval,
            "pvalue_fdr": np.nan,
            "significant_fdr": False,
            "n_sig_celltypes_total": n_sig_ct_total,
            "n_sig_celltypes_in_area": n_sig_ct_in_area,
            "n_animals": n_animals,
            "n_rows": n_rows,
            "model": method,
            "converged": converged,
        })

df_bh_only = pd.DataFrame(results)

# FDR across all non-null tests
valid_mask = df_bh_only["pvalue"].notna()
if valid_mask.sum() > 0:
    df_bh_only.loc[valid_mask, "pvalue_fdr"] = multipletests(
        df_bh_only.loc[valid_mask, "pvalue"],
        method="fdr_bh"
    )[1]

df_bh_only["significant_fdr"] = df_bh_only["pvalue_fdr"] < 0.05
df_bh_only.to_csv(out_results, index=False)

print("\nSaved:", out_results)


# =========================================================
# 6) Heatmap for BH-significant-celltype-only rerun
# =========================================================
pivot_beta = df_bh_only.pivot(index="gwas", columns="area", values="beta_z")
pivot_sig_fdr = df_bh_only.pivot(index="gwas", columns="area", values="significant_fdr")

pivot_beta = pivot_beta.reindex(columns=areas_present)
pivot_sig_fdr = pivot_sig_fdr.reindex(columns=areas_present)

# Drop GWAS where all betas are NaN after filtering
keep_rows = ~pivot_beta.isna().all(axis=1)
pivot_beta = pivot_beta.loc[keep_rows]
pivot_sig_fdr = pivot_sig_fdr.loc[keep_rows]

# Order GWAS by number of FDR-significant areas, then max |beta|
row_order = (
    pivot_sig_fdr.sum(axis=1).sort_values(ascending=False).index.tolist()
)
pivot_beta = pivot_beta.loc[row_order]
pivot_sig_fdr = pivot_sig_fdr.loc[row_order]

n_rows = len(pivot_beta.index)
n_cols = len(areas_present)
fig_width = 10
fig_height = max(4, fig_width * (n_rows / max(n_cols, 1)))

fig, ax = plt.subplots(figsize=(fig_width, fig_height))

vals = pivot_beta.values.astype(float)
if np.isfinite(vals).any():
    vmax = np.nanmax(np.abs(vals))
else:
    vmax = 1.0
vmin = -vmax

im = ax.imshow(vals, aspect="equal", cmap="RdBu_r", vmin=vmin, vmax=vmax, interpolation="nearest")

for i in range(len(pivot_beta.index)):
    for j in range(len(areas_present)):
        val = pivot_beta.iloc[i, j]
        if pd.notna(val) and bool(pivot_sig_fdr.iloc[i, j]):
            ax.text(
                j, i, "*",
                ha="center", va="center",
                color="black", fontsize=16, fontweight="bold",
                clip_on=False
            )

ax.set_xticks(np.arange(-0.5, len(areas_present), 1), minor=True)
ax.set_yticks(np.arange(-0.5, len(pivot_beta.index), 1), minor=True)
ax.grid(which="minor", color="gray", linestyle="-", linewidth=0.5, alpha=0.3)

ax.set_yticks(range(len(pivot_beta.index)))
ax.set_yticklabels(pivot_beta.index, fontsize=10)
ax.set_xticks(range(len(areas_present)))
ax.set_xticklabels(areas_present, fontsize=12, fontweight="bold")

ax.xaxis.tick_top()
ax.xaxis.set_label_position("top")

plt.colorbar(im, ax=ax, label="Beta coefficient (effect size)", fraction=0.046, pad=0.04)

title_text = (
    f"Two-Stage Analysis (BH-significant cell types only)\n"
    f"Stage 1: celltype x GWAS BH-significant only ({bh_mode})\n"
    f"Stage 2: ESmu ~ GWAS z-score + (1|{animal_col}), re-fit within each area\n"
    f"* = FDR significant projection (q < 0.05)"
)
ax.set_title(title_text, fontsize=11, fontweight="bold", pad=25, y=1.10)

plt.tight_layout()
plt.savefig(out_heatmap, bbox_inches="tight", dpi=300)
plt.show()


# =========================================================
# 7) Summary
# =========================================================
print("\n" + "=" * 80)
print("BH-SIGNIFICANT-CELLTYPE-ONLY TWO-STAGE SUMMARY")
print("=" * 80)
print(f"\nGWAS included: {len(pivot_beta.index)}")
print(f"Brain areas: {len(areas_present)}")
print(f"Total tests attempted: {len(df_bh_only)}")
print(f"Valid tests fit: {int(valid_mask.sum())}")
print(f"FDR-significant results: {int(df_bh_only['significant_fdr'].sum())}")

if valid_mask.sum() > 0:
    pct = 100 * df_bh_only["significant_fdr"].sum() / valid_mask.sum()
    print(f"Percentage significant among valid tests: {pct:.1f}%")

print(f"\n{'Area':<6} {'Significant':>11}")
print("-" * 20)
for area in areas_present:
    n_sig = int(df_bh_only.loc[df_bh_only["area"] == area, "significant_fdr"].sum())
    print(f"{area:<6} {n_sig:>11}")

print(f"\n{'GWAS':<50} {'N Areas':>8}")
print("-" * 60)
gwas_summary = pivot_sig_fdr.sum(axis=1).sort_values(ascending=False)
for gwas, n_areas in gwas_summary.head(15).items():
    gwas_short = gwas[:47] + "..." if len(gwas) > 50 else gwas
    print(f"{gwas_short:<50} {int(n_areas):>8}")


# =========================================================
# 8) Compare against original two-stage results if present
# =========================================================
shared_gwas = []

if Path("mixedmodel_twostage_results.csv").exists():
    df_original_stage2 = pd.read_csv("mixedmodel_twostage_results.csv")

    # P-value distribution comparison
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))

    ax1.hist(df_original_stage2["pvalue"].dropna(), bins=50, edgecolor="black", alpha=0.7)
    ax1.axvline(0.05, color="red", linestyle="--", label="p=0.05")
    ax1.set_xlabel("P-value", fontsize=11)
    ax1.set_ylabel("Frequency", fontsize=11)
    ax1.set_title(f"Original Two-Stage\n({df_original_stage2['pvalue'].notna().sum()} valid tests)", fontsize=11)
    ax1.legend()

    ax2.hist(df_bh_only["pvalue"].dropna(), bins=50, edgecolor="black", alpha=0.7)
    ax2.axvline(0.05, color="red", linestyle="--", label="p=0.05")
    ax2.set_xlabel("P-value", fontsize=11)
    ax2.set_ylabel("Frequency", fontsize=11)
    ax2.set_title(f"BH-significant only rerun\n({df_bh_only['pvalue'].notna().sum()} valid tests)", fontsize=11)
    ax2.legend()

    plt.tight_layout()
    plt.savefig(out_pval_compare, dpi=300, bbox_inches="tight")
    plt.show()

    # Side-by-side heatmap comparison on shared GWAS
    orig_beta = df_original_stage2.pivot(index="gwas", columns="area", values="beta_z").reindex(columns=areas_present)
    orig_sig = df_original_stage2.pivot(index="gwas", columns="area", values="significant_fdr").reindex(columns=areas_present)

    bh_beta = pivot_beta.copy()
    bh_sig = pivot_sig_fdr.copy()

    shared_gwas = [g for g in orig_beta.index if g in bh_beta.index]

    if len(shared_gwas) > 0:
        orig_beta2 = orig_beta.loc[shared_gwas]
        orig_sig2 = orig_sig.loc[shared_gwas]
        bh_beta2 = bh_beta.loc[shared_gwas]
        bh_sig2 = bh_sig.loc[shared_gwas]

        concat_vals = np.concatenate([
            orig_beta2.values.flatten(),
            bh_beta2.values.flatten()
        ])
        concat_vals = concat_vals[np.isfinite(concat_vals)]
        vmax_compare = np.max(np.abs(concat_vals)) if len(concat_vals) > 0 else 1.0

        n_rows_compare = len(shared_gwas)
        fig_width = 14
        fig_height = max(4, fig_width * (n_rows_compare / (2 * max(len(areas_present), 1))))

        fig, axes = plt.subplots(1, 2, figsize=(fig_width, fig_height), sharey=True)

        for ax, beta_mat, sig_mat, title in [
            (axes[0], orig_beta2, orig_sig2, "Original two-stage"),
            (axes[1], bh_beta2, bh_sig2, "BH-significant celltypes only"),
        ]:
            im = ax.imshow(
                beta_mat.values,
                aspect="equal",
                cmap="RdBu_r",
                vmin=-vmax_compare,
                vmax=vmax_compare,
                interpolation="nearest"
            )

            for i in range(len(beta_mat.index)):
                for j in range(len(areas_present)):
                    val = beta_mat.iloc[i, j]
                    if pd.notna(val) and bool(sig_mat.iloc[i, j]):
                        ax.text(
                            j, i, "*",
                            ha="center", va="center",
                            color="black", fontsize=14, fontweight="bold"
                        )

            ax.set_xticks(np.arange(-0.5, len(areas_present), 1), minor=True)
            ax.set_yticks(np.arange(-0.5, len(beta_mat.index), 1), minor=True)
            ax.grid(which="minor", color="gray", linestyle="-", linewidth=0.5, alpha=0.3)

            ax.set_xticks(range(len(areas_present)))
            ax.set_xticklabels(areas_present, fontsize=11, fontweight="bold")
            ax.xaxis.tick_top()
            ax.xaxis.set_label_position("top")
            ax.set_title(title, fontsize=12, fontweight="bold", pad=20)

        axes[0].set_yticks(range(len(shared_gwas)))
        axes[0].set_yticklabels(shared_gwas, fontsize=9)

        plt.colorbar(im, ax=axes.ravel().tolist(), fraction=0.025, pad=0.02, label="Beta coefficient (effect size)")
        plt.tight_layout()
        plt.savefig(out_compare, dpi=300, bbox_inches="tight")
        plt.show()

    print("\nSaved comparison plots:")
    print(" -", out_pval_compare)
    if len(shared_gwas) > 0:
        print(" -", out_compare)

print("\nSaved:")
print(" -", out_results)
print(" -", out_heatmap)


# In[11]:


# Side-by-side heatmap comparison on shared GWAS
orig_beta = df_original_stage2.pivot(index="gwas", columns="area", values="beta_z").reindex(columns=areas_present)
orig_sig = df_original_stage2.pivot(index="gwas", columns="area", values="significant_fdr").reindex(columns=areas_present)

bh_beta = pivot_beta.copy()
bh_sig = pivot_sig_fdr.copy()

shared_gwas = [g for g in orig_beta.index if g in bh_beta.index]

if len(shared_gwas) > 0:
    orig_beta2 = orig_beta.loc[shared_gwas]
    orig_sig2 = orig_sig.loc[shared_gwas]
    bh_beta2 = bh_beta.loc[shared_gwas]
    bh_sig2 = bh_sig.loc[shared_gwas]

    # Shared symmetric color scale
    concat_vals = np.concatenate([
        orig_beta2.values.flatten(),
        bh_beta2.values.flatten()
    ])
    concat_vals = concat_vals[np.isfinite(concat_vals)]
    vmax_compare = np.max(np.abs(concat_vals)) if len(concat_vals) > 0 else 1.0

    n_rows_compare = len(shared_gwas)
    n_cols_compare = len(areas_present)

    # Make figure a bit taller and give colorbar its own axis
    fig_height = max(3.5, 0.5 * n_rows_compare + 1.5)
    fig = plt.figure(figsize=(14, fig_height), constrained_layout=True)
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 0.06])

    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1], sharey=ax1)
    cax = fig.add_subplot(gs[0, 2])

    for ax, beta_mat, sig_mat, title in [
        (ax1, orig_beta2, orig_sig2, "Original two-stage"),
        (ax2, bh_beta2, bh_sig2, "BH-significant celltypes only"),
    ]:
        im = ax.imshow(
            beta_mat.values,
            aspect="auto",          # less cramped than "equal"
            cmap="RdBu_r",
            vmin=-vmax_compare,
            vmax=vmax_compare,
            interpolation="nearest"
        )

        # Significance stars
        for i in range(len(beta_mat.index)):
            for j in range(len(areas_present)):
                val = beta_mat.iloc[i, j]
                if pd.notna(val) and bool(sig_mat.iloc[i, j]):
                    ax.text(
                        j, i, "*",
                        ha="center", va="center",
                        color="black", fontsize=13, fontweight="bold"
                    )

        # Grid
        ax.set_xticks(np.arange(-0.5, len(areas_present), 1), minor=True)
        ax.set_yticks(np.arange(-0.5, len(beta_mat.index), 1), minor=True)
        ax.grid(which="minor", color="gray", linestyle="-", linewidth=0.5, alpha=0.3)

        # Ticks
        ax.set_xticks(range(len(areas_present)))
        ax.set_xticklabels(areas_present, fontsize=11, fontweight="bold")
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position("top")
        ax.tick_params(axis="x", pad=6, length=0)
        ax.tick_params(axis="y", length=0)

        ax.set_title(title, fontsize=12, fontweight="bold", pad=18)

    # Y labels only on left panel
    ax1.set_yticks(range(len(shared_gwas)))
    ax1.set_yticklabels(shared_gwas, fontsize=9)
    ax2.tick_params(axis="y", labelleft=False)

    # Dedicated colorbar axis
    cbar = fig.colorbar(im, cax=cax)
    cbar.set_label("Beta coefficient (effect size)", rotation=90, labelpad=10)
    cbar.ax.tick_params(labelsize=9)

    plt.savefig(out_compare, dpi=300, bbox_inches="tight")
    plt.show()


# In[ ]:




