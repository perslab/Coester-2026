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
import statsmodels.stats.multitest as smm
import statsmodels.api as sm
import statsmodels.formula.api as smf
import adjustText
from scipy.stats import pearsonr, spearmanr
from scipy.stats import chi2_contingency
from scipy import stats
os.chdir("/projects/perslab/people/qbx224/7plex_paper/atlas/")
import bipca
from bipca import plotting
from bipca.utils import stabilize_matrix
from bipca.experiments import new_svd as SVD
from openTSNE import TSNE
from multiprocessing import Pool
import matplotlib.pyplot as plt
from tqdm.contrib.concurrent import thread_map,process_map
from functools import partial
import torch
from threadpoolctl import threadpool_limits


# In[2]:


adata = sc.read("allgenes_bipca_hypomap.h5ad")
adata


# In[3]:


ax = sc.pl.embedding(
    adata,
    basis="tsne",
    color="neurotransmitter",
    legend_loc="on data",
    legend_fontsize="xx-small",
    return_fig=True)
plt.show()


# In[4]:


celltypes = ["GABA", "Glut", "Hist", "Chat"]
adata.obs['plot'] = adata.obs['neurotransmitter'].where(
    adata.obs['neurotransmitter'].isin(celltypes)
)

sc.pl.embedding(
    adata,
    basis="tsne",
    color="plot",
    legend_loc="on data",
    legend_fontsize="xx-small"
)


# In[5]:


adata = adata[adata.obs['plot'].notna()].copy()
adata


# In[6]:


print(adata.X[:5, :5].toarray())


# In[7]:


X_sub, (cell_idx, gene_idx) = stabilize_matrix(adata.X, threshold=100)
X_sub


# In[8]:


adata = adata[cell_idx, gene_idx].copy()
adata


# In[9]:


adata.X = X_sub
adata.obs.head()


# In[10]:


n_threads = 16
torch.set_num_threads(n_threads)
with threadpool_limits(limits=n_threads):
    op = bipca.BiPCA(n_components=-1,seed=42) # get the BiPCA operator, here n_components=-1 is to do full SVD
    Z = op.fit_transform(X_sub.toarray()) # take the count data as input and run bipca
    op.get_plotting_spectrum() # get the fitting parameters


# In[11]:


op.write_to_adata(adata) # save the results back to anndata


# In[12]:


plotting.MP_histograms_from_bipca(adata,median=False,bins=100) 


# In[13]:


adata.obsm['pcs'] = SVD(adata.layers['Z_biwhite'],adata.uns["bipca"]['rank'])
sc.pp.neighbors(adata, use_rep = "pcs",method="gauss", n_pcs=adata.uns["bipca"]['rank'])
adata.obsm['tsne'] = np.array(TSNE(n_jobs=32).fit(adata.obsm['pcs']))


# In[37]:


sc.pl.embedding(adata,
                basis="tsne",
                color="detailed_celltype_annotation", 
                legend_loc="on data", 
                legend_fontsize="xx-small",
                save="_detailed_celltype_annotation_recluster.pdf", 
               )


# In[17]:


resolutions = [0.2, 0.4, 0.6, 0.8, 1.0]

for res in resolutions:
    sc.tl.leiden(
        adata,
        resolution=res,
        random_state=42,
        flavor="leidenalg",
        n_iterations=-1,
        key_added=f"leiden_{res:.1f}")
    print(f"leiden_{res:.1f}")


# In[19]:


resolutions = [1.2, 1.4, 1.6, 1.8, 2.0]

for res in resolutions:
    sc.tl.leiden(
        adata,
        resolution=res,
        random_state=42,
        flavor="leidenalg",
        n_iterations=-1,
        key_added=f"leiden_{res:.1f}")
    print(f"leiden_{res:.1f}")


# In[20]:


leiden = ['leiden_0.2', 
          'leiden_0.4', 
          'leiden_0.6', 
          'leiden_0.8', 
          'leiden_1.0', 
          'leiden_1.2', 
          'leiden_1.4', 
          'leiden_1.6', 
          'leiden_1.8', 
          'leiden_2.0']
sc.pl.embedding(adata,
                basis="tsne",
                color=leiden, 
                legend_loc="on data", 
                legend_fontsize="x-small",
                save="_leiden_neurons_bipca.pdf"
               )


# In[ ]:


#adata.write("neurons_bipca_reclustered.h5ad")
adata = sc.read("neurons_bipca_reclustered.h5ad")
adata


# In[ ]:


adata.obs['n_counts'] = adata.X.sum(axis=1)
adata.obs['n_genes'] = (adata.X > 0).sum(axis=1)
scales_counts = sc.pp.normalize_total(adata, target_sum=None, inplace=False)
# log1p transform
adata.layers["log1p_norm"] = sc.pp.log1p(scales_counts["X"], copy=True)
fig, axes = plt.subplots(1, 2, figsize=(10, 5))
p1 = sns.histplot(adata.obs["n_counts"], bins=100, kde=False, ax=axes[0])
axes[0].set_title("Total counts")
p2 = sns.histplot(adata.layers["log1p_norm"].sum(1), bins=100, kde=False, ax=axes[1])
axes[1].set_title("Shifted logarithm")
plt.show()


# In[5]:


sc.pl.embedding(
    adata,
    basis="tsne",                       
    color="Yengo_BMI",
    cmap="RdBu_r"
)


# In[6]:


# Use normalized expression (use_raw=False)
df = sc.get.obs_df(
    adata,
    keys=["Yengo_BMI", "v147"],
    use_raw=False
)

# Drop cells with missing values
df = df.dropna(subset=["Yengo_BMI", "v147"])

x = df["Yengo_BMI"].values
y = df["v147"].values

pear_r, pear_p = pearsonr(x, y)
spear_r, spear_p = spearmanr(x, y)


# In[7]:


plt.figure(figsize=(5, 4))

sns.regplot(
    data=df,
    x="Yengo_BMI",
    y="v147",
    scatter_kws={"s": 10, "alpha": 0.4},
    line_kws={"linewidth": 2}
)

plt.xlabel("Yengo_BMI")
plt.ylabel("Normalized expression of v147")
plt.title(
    f"v147 vs Yengo_BMI\n"
    f"Pearson r={pear_r:.2f}, p={pear_p:.1e}; "
    f"Spearman ρ={spear_r:.2f}, p={spear_p:.1e}"
)

plt.tight_layout()
plt.show()


# In[8]:


leiden_cols = [
    'leiden_0.2',
    'leiden_0.4',
    'leiden_0.6',
    'leiden_0.8',
    'leiden_1.0',
    'leiden_1.2',
    'leiden_1.4',
    'leiden_1.6',
    'leiden_1.8',
    'leiden_2.0',
]

for col in leiden_cols:
    # extract resolution part, e.g. "0.2" from "leiden_0.2"
    res_str = col.split('_')[1]     # "0.2"
    res_tag = res_str.replace('.', '')  # "02", "04", ..., "20"
    # make sure they’re strings, append suffix, and (optionally) back to category
    adata.obs[col] = (
        adata.obs[col].astype(str) + f"_res{res_tag}"
    ).astype('category')
adata.obs["leiden_2.0"]


# In[38]:


sc.pl.embedding(adata,
                basis="tsne",
                color=leiden_cols, 
                legend_loc="on data", 
                legend_fontsize="x-small",
                save="_leiden_neurons_recluster.pdf"
               )


# In[11]:


counts = pd.DataFrame(adata.layers["log1p_norm"].toarray(),
                  index=adata.obs_names,
                  columns=adata.var_names)
counts.to_parquet("neurons_counts.parquet")


# In[12]:


leiden = ['leiden_0.2', 
          'leiden_0.4', 
          'leiden_0.6', 
          'leiden_0.8', 
          'leiden_1.0', 
          'leiden_1.2', 
          'leiden_1.4', 
          'leiden_1.6', 
          'leiden_1.8', 
          'leiden_2.0']
for res in leiden:
    meta = adata.obs[res].to_frame()
    meta.to_parquet(f"{res}_clusters.parquet")


# In[16]:


#Run CELLECT with config_leiden.yaml


# In[10]:


df = pd.read_csv('/projects/perslab/people/qbx224/CELLECT/CELLECT-leiden/CELLECT-LDSC/results/prioritization.csv')
df


# In[11]:


leiden_cols = [
    'leiden_0.2',
    'leiden_0.4',
    'leiden_0.6',
    'leiden_0.8',
    'leiden_1.0',
    'leiden_1.2',
    'leiden_1.4',
    'leiden_1.6',
    'leiden_1.8',
    'leiden_2.0',]
gwas = "BMI_Yengo2018"

for col in leiden_cols:
    # col example: "leiden_0.2"
    res_str = col.split('_')[1]      # "0.2"
    res_tag = res_str.replace('.', '')  # "02"
    spec_id = f"leiden_{res_tag}"    # "leiden_02"

    # subset df to this GWAS + resolution
    sub = df[(df["gwas"] == gwas) & (df["specificity_id"] == spec_id)]

    # build mapping: annotation ("3_res02") -> pvalue
    ann_to_p = sub.set_index("annotation")["pvalue"].to_dict()

    # new column name in adata.obs, e.g. "BMI_leiden_0.2"
    new_col = f"BMI_{col}"

    # map cluster labels in adata.obs[col] to p-values
    adata.obs[new_col] = adata.obs[col].map(ann_to_p).astype(float)


# In[39]:


BMI_leiden = ['BMI_leiden_0.2', 'BMI_leiden_0.4', 'BMI_leiden_0.6', 'BMI_leiden_0.8', 'BMI_leiden_1.0', 'BMI_leiden_1.2', 'BMI_leiden_1.4', 'BMI_leiden_1.6', 'BMI_leiden_1.8', 'BMI_leiden_2.0']
sc.pl.embedding(adata,
                basis="tsne",
                color=BMI_leiden,
                color_map="flare",
                save="_leiden_BMI_pvalue.pdf")


# In[13]:


df = pd.read_csv('/projects/perslab/people/qbx224/CELLECT/CELLECT-leiden/CELLECT-LDSC/results/prioritization.csv')
df["pvalue_BH"] = np.nan

for col in leiden_cols:
    # convert "leiden_0.2" → "leiden_02"
    res_str = col.split("_")[1]          # "0.2"
    res_tag = res_str.replace(".", "")   # "02"
    spec_id = f"leiden_{res_tag}"       # "leiden_02"

    # subset to this GWAS + resolution
    mask = (df["gwas"] == gwas) & (df["specificity_id"] == spec_id)
    sub = df.loc[mask].copy()
    if sub.empty:
        print(f"No rows for {spec_id}, skipping")
        continue

    # BH-correct pvalues for this resolution
    _, pvals_corr, _, _ = smm.multipletests(sub["pvalue"].values, method="fdr_bh")
    sub["pvalue_BH"] = pvals_corr

    # write back into the main df (optional but nice to keep)
    df.loc[sub.index, "pvalue_BH"] = sub["pvalue_BH"]

    # build mapping annotation → BH pvalue
    metadata_dict = dict(zip(sub["annotation"], sub["pvalue_BH"]))

    # new column name in adata.obs, e.g. "BMI_leiden_0.2_BH"
    new_col = f"BMI_{col}_BH"

    # map cluster labels in adata.obs[col] (e.g. "3_res02") to BH pvalue
    adata.obs[new_col] = adata.obs[col].map(metadata_dict)
    adata.obs[new_col] = pd.to_numeric(adata.obs[new_col], errors="coerce")


# In[40]:


BMI_leiden_BH = ['BMI_leiden_0.2_BH', 'BMI_leiden_0.4_BH', 'BMI_leiden_0.6_BH', 'BMI_leiden_0.8_BH', 'BMI_leiden_1.0_BH', 'BMI_leiden_1.2_BH', 'BMI_leiden_1.4_BH', 'BMI_leiden_1.6_BH', 'BMI_leiden_1.8_BH', 'BMI_leiden_2.0_BH'] 
sc.pl.embedding(adata,
                basis="tsne",
                color=BMI_leiden_BH,
                color_map="flare",
                save="_leiden_BMI_pvalue_BH.pdf")


# In[15]:


sc.pl.embedding(adata,
                basis="tsne",
                color="detailed_celltype_annotation", 
                legend_loc="on data", 
                legend_fontsize="xx-small"
               )


# In[16]:


df["neglog10_p"] = -np.log10(df["pvalue"])
leiden_cols = [
    "leiden_0.2",
    "leiden_0.4",
    "leiden_0.6",
    "leiden_0.8",
    "leiden_1.0",
    "leiden_1.2",
    "leiden_1.4",
    "leiden_1.6",
    "leiden_1.8",
    "leiden_2.0",
]
gwas = "BMI_Yengo2018"

for col in leiden_cols:
    # convert "leiden_0.2" → "leiden_02" to match df["specificity_id"]
    res_str = col.split("_")[1]        # "0.2"
    res_tag = res_str.replace(".", "") # "02"
    spec_id = f"leiden_{res_tag}"      # "leiden_02"

    # subset df to this GWAS + resolution
    mask = (df["gwas"] == gwas) & (df["specificity_id"] == spec_id)
    sub = df.loc[mask].copy()
    if sub.empty:
        print(f"No rows for {spec_id}, skipping")
        continue

    # mapping: annotation ("3_res02") → −log10(BH p)
    ann_to_logp = dict(zip(sub["annotation"], sub["neglog10_p"]))

    # create a new column in adata.obs, e.g. "BMI_leiden_0.2_log10"
    new_col = f"BMI_{col}_log10"

    adata.obs[new_col] = adata.obs[col].map(ann_to_logp)
    adata.obs[new_col] = pd.to_numeric(adata.obs[new_col], errors="coerce")


# In[41]:


for col in leiden_cols:
    log_col = f"BMI_{col}_log10"   # e.g. "BMI_leiden_0.2_log10"
    if log_col not in adata.obs.columns:
        print(f"{log_col} not found in adata.obs, skipping")
        continue

    df_obs = adata.obs

    # Explicitly set observed=False to match old behavior and silence the FutureWarning
    order = (
        df_obs.groupby(col, observed=False)[log_col]
              .median()
              .sort_values()
              .index
              .tolist()
    )

    plt.figure(figsize=(12, 6))
    sns.set_style("whitegrid")

    sns.pointplot(
        data=df_obs,
        x=col,
        y=log_col,
        order=order,
        # new API bits:
        errorbar="sd",             # replaces ci='sd'
        capsize=0.1,
        err_kws={"linewidth": 1},  # replaces errwidth=1
        linestyle="none",          # replaces join=False
    )

    plt.axhline(3.1, color="red", linestyle="--", label="y = 3.1")

    plt.xticks(rotation=45, ha="right")
    plt.xlabel(f"Leiden clusters ({col})")
    plt.ylabel("-log10(BH p-value)")
    plt.title(f"Cellect BMI by Leiden clusters ({col})")
    plt.tight_layout()

    plt.savefig(f"cellect_BMI_{col}_clusters_log10.png", dpi=300)
    plt.show()
    # plt.close()


# In[29]:


outdir = "area_parquets"
os.makedirs(outdir, exist_ok=True)

leiden_cols = [
    "leiden_0.2",
    "leiden_0.4",
    "leiden_0.6",
    "leiden_0.8",
    "leiden_1.0",
    "leiden_1.2",
    "leiden_1.4",
    "leiden_1.6",
    "leiden_1.8",
    "leiden_2.0",
]

for area in adata.obs["area"].unique():
    # subset AnnData by area
    mask = adata.obs["area"] == area
    sub = adata[mask].copy()

    # build counts DataFrame (one per area)
    counts = pd.DataFrame(
        sub.layers["scran_normalization"].toarray(),
        index=sub.obs_names,
        columns=sub.var_names,
    )

    # save counts to parquet
    area_safe = str(area).replace(" ", "_")
    counts.to_parquet(os.path.join(outdir, f"{area_safe}_counts.parquet"))

    # now save *one meta parquet per Leiden resolution* for this area
    for col in leiden_cols:
        if col not in sub.obs.columns:
            print(f"Warning: {col} not in sub.obs for area={area!r}, skipping")
            continue

        meta = sub.obs[[col]].copy()  # single column DataFrame
        # make filename filesystem-friendly (no dots in the col name)
        col_safe = col.replace(".", "_")  # e.g. "leiden_2.0" → "leiden_2_0"

        meta.to_parquet(
            os.path.join(outdir, f"{area_safe}_{col_safe}_clusters.parquet")
        )

    print(f"Saved parquet files for area={area!r}")


# In[2]:


###Run CELLEX in environment
#adata.write("neurons_bipca_leiden.h5ad")
adata = sc.read("neurons_bipca_leiden.h5ad")
adata


# In[3]:


indir = "leiden_cellex_results"
records = []

for fn in os.listdir(indir):
    if not fn.endswith("_cellex.csv"):
        continue
    path = os.path.join(indir, fn)

    # strip suffix → 'AreaX_leiden_0_2'
    base = fn[:-len("_cellex.csv")]

    # split into area and leiden part
    if "_leiden_" not in base:
        print(f"Skipping {fn}: no '_leiden_' in name")
        continue
    area, leiden_part = base.split("_leiden_", 1)  # e.g. area='AreaX', leiden_part='0_2'

    # reconstruct the resolution label if you like (optional)
    # '0_2' → 'leiden_0.2'
    res_str = leiden_part.replace("_", ".")
    resolution = f"leiden_{res_str}"

    # read CELLEX results
    df = pd.read_csv(path, index_col=0)
    if "v147" not in df.index:
        print(f"Skipping {fn}: v147 not found")
        continue

    # v147 row: index = clusters, values = ESMU
    s = df.loc["v147"]

    # make long-format records: one per cluster
    for cluster, esmu in s.items():
        records.append({
            "area": area,
            "resolution": resolution,
            "cluster": cluster,
            "esmu": esmu,
        })

# one big long DataFrame
combined_long = pd.DataFrame(records)
print(combined_long)


# In[4]:


# resolutions and matching obs fields
resolutions = [
    0.2, 0.4, 0.6, 0.8, 1.0, 
    1.2, 1.4, 1.6, 1.8, 2.0
]

records = []

for res in resolutions:
    leiden_col = f"leiden_{res}"
    bmi_col    = f"BMI_leiden_{res}_log10"

    # dataframe: cell → cluster label at this resolution
    clusters = adata.obs[leiden_col]

    # dataframe: cell → BMI score for this resolution
    bmi_vals = adata.obs[bmi_col]

    # group BMI per cluster
    bmi_per_cluster = bmi_vals.groupby(clusters).mean()

    # store long format
    for cluster, bmi in bmi_per_cluster.items():
        records.append({
            "resolution": f"leiden_{res}",
            "cluster": cluster,        
            "bmi_neglog10": bmi,      
        })

bmi_long = pd.DataFrame(records)
bmi_long


# In[19]:


df_plot = combined_long.merge(
    bmi_long,
    on=["resolution", "cluster"],
    how="inner"
)
df_plot


# In[6]:


plt.figure(figsize=(8, 6))

sns.scatterplot(
    data=df_plot,
    x="bmi_neglog10",
    y="esmu",
    hue="area",
    alpha=0.8,
    s=40
)

plt.xlabel("-log10(BMI)")
plt.ylabel("ESMU")
plt.tight_layout()
plt.savefig("all_clusters_log10.png", dpi=300)
plt.show()


# In[7]:


sns.relplot(
    data=df_plot,
    x="bmi_neglog10",
    y="esmu",
    hue="area",
    col="resolution",
    col_wrap=4,
    kind="scatter",
    alpha=0.7
)


# In[8]:


corr_per_area = (
    df_plot
    .groupby(["resolution", "area"])
    .apply(lambda g: g["esmu"].corr(g["bmi_neglog10"]))
    .reset_index(name="corr")
)

corr_per_area


# In[9]:


corr_per_resolution = (
    df_plot
    .groupby("resolution")
    .apply(lambda g: g["esmu"].corr(g["bmi_neglog10"]))
    .reset_index(name="corr")
)

corr_per_resolution


# In[10]:


final_corr_table = (
    df_plot
    .groupby(["resolution", "area"])
    .apply(lambda g: g["esmu"].corr(g["bmi_neglog10"]))
    .unstack("area")                         # columns = areas
)

# Add overall collective correlation per resolution
final_corr_table["ALL"] = (
    df_plot
    .groupby("resolution")
    .apply(lambda g: g["esmu"].corr(g["bmi_neglog10"]))
)

final_corr_table


# In[11]:


plt.figure(figsize=(10,5))
sns.heatmap(final_corr_table, annot=True, cmap="coolwarm", center=0)
plt.title("Correlation(ESMU vs -log10(BMI)) per area and per resolution")
plt.savefig("ESmu_clusters_log10.png", dpi=300)
plt.show()


# In[24]:


df_plot


# In[21]:


# --------- 1. Build cluster-level covariates for each resolution ----------
# Resolutions list inferred from df_plot/resolution unique values:
resolutions = sorted(df_plot["resolution"].unique())

cluster_stats_records = []
total_counts = adata.obs['n_counts'].values
print(total_counts)
n_genes_per_cell = adata.obs['n_genes'].values
print(n_genes_per_cell)


# In[23]:


# We'll also capture number of nuclei per cluster (cluster_size)
for res in resolutions:
    # map resolution name to the obs leiden column (they should match)
    leiden_col = res  # e.g., 'leiden_1.0'
    if leiden_col not in adata.obs.columns:
        # try replacing '.' with '_' if necessary
        alt = leiden_col.replace('.', '_')
        if alt in adata.obs.columns:
            leiden_col = alt
        else:
            raise KeyError(f"Leiden column {res} or {alt} not found in adata.obs")

    cluster_labels = adata.obs[leiden_col].astype(str).values  # ensure str
    # Build a temporary DataFrame per cell to group easily
    tmp = pd.DataFrame({
        "cluster": cluster_labels,
        "n_counts": total_counts,
        "n_genes": n_genes_per_cell
    }, index=adata.obs_names)

    # group to cluster-level stats
    g = tmp.groupby("cluster").agg(
        cluster_size = ("n_counts", "size"),
        mean_total_counts = ("n_counts", "mean"),
        mean_n_genes = ("n_genes", "mean")
    ).reset_index()

    # attach resolution label and push to records
    g["resolution"] = res
    cluster_stats_records.append(g)

cluster_stats = pd.concat(cluster_stats_records, ignore_index=True)
# normalize types
cluster_stats["cluster"] = cluster_stats["cluster"].astype(str)
cluster_stats


# In[25]:


# --------- 2. Merge cluster stats into df_plot ----------
# Ensure the cluster strings in df_plot are str
df_plot["cluster"] = df_plot["cluster"].astype(str)

merged = df_plot.merge(
    cluster_stats,
    on=["resolution","cluster"],
    how="left", validate="m:1"
)

# Report any missing merges
n_missing = merged["cluster_size"].isna().sum()
if n_missing:
    print(f"[WARN] {n_missing} rows in df_plot had no matching cluster stats. Check naming conventions.")

# Add log cluster size
merged["log_cluster_size"] = np.log(merged["cluster_size"].astype(float) + 1.0)  # +1 to be safe
merged


# In[26]:


# --------- 3. Function to fit model and collect summary ----------
from collections import defaultdict

def fit_model_and_summarize(df, outcome="bmi_neglog10", predictor="esmu", covariate="log_cluster_size"):
    """
    Fit outcome ~ predictor + covariate using OLS (with intercept).
    Returns dict with coef, se, pvalue, 95% CI for predictor term, plus n_obs.
    """
    formula = f"{outcome} ~ {predictor} + {covariate}"
    model = smf.ols(formula=formula, data=df).fit()
    # Extract predictor row (predictor)
    if predictor not in model.params.index:
        # unexpected
        return None
    coef = model.params[predictor]
    se = model.bse[predictor]
    pval = model.pvalues[predictor]
    ci_low, ci_high = model.conf_int().loc[predictor].values
    return {
        "coef": float(coef),
        "se": float(se),
        "pval": float(pval),
        "ci_lower": float(ci_low),
        "ci_upper": float(ci_high),
        "n": int(model.nobs),
        "r2": float(model.rsquared)
    }


# In[27]:


# --------- 4. Per-area, per-resolution models ----------
rows = []
for res in resolutions:
    df_res = merged[merged["resolution"] == res].copy()
    for area in sorted(df_res["area"].unique()):
        df_sub = df_res[df_res["area"] == area]
        # require minimal number of points to fit
        if len(df_sub) < 5:
            rows.append({
                "resolution": res, "area": area,
                "note": "too_few_obs"
            })
            continue
        res_summary = fit_model_and_summarize(df_sub)
        if res_summary is None:
            rows.append({
                "resolution": res, "area": area,
                "note": "fit_failed"
            })
        else:
            rows.append({
                "resolution": res,
                "area": area,
                "coef": res_summary["coef"],
                "se": res_summary["se"],
                "pval": res_summary["pval"],
                "ci_lower": res_summary["ci_lower"],
                "ci_upper": res_summary["ci_upper"],
                "n": res_summary["n"],
                "r2": res_summary["r2"],
                "note": ""
            })

corr_area_res_df = pd.DataFrame(rows)
corr_area_res_df


# In[28]:


# --------- 5. Pooled model per resolution with area fixed effects ----------
pooled_rows = []
for res in resolutions:
    df_res = merged[merged["resolution"] == res].copy()
    if len(df_res) < 10:
        pooled_rows.append({"resolution": res, "note": "too_few_obs"})
        continue
    # formula with area fixed effects
    formula = "bmi_neglog10 ~ esmu + log_cluster_size + C(area)"
    model = smf.ols(formula=formula, data=df_res).fit()
    # extract esmu coefficient
    coef = model.params.get("esmu", np.nan)
    se = model.bse.get("esmu", np.nan)
    pval = model.pvalues.get("esmu", np.nan)
    ci_low, ci_high = model.conf_int().loc["esmu"].values
    pooled_rows.append({
        "resolution": res,
        "coef": float(coef),
        "se": float(se),
        "pval": float(pval),
        "ci_lower": float(ci_low),
        "ci_upper": float(ci_high),
        "n": int(model.nobs),
        "r2": float(model.rsquared),
        "note": ""
    })

pooled_corr_df = pd.DataFrame(pooled_rows)
pooled_corr_df


# In[34]:


# --------- 7. Visualization: coef plot across areas, faceted by resolution ----------

# Keep only fitted rows
plot_df = corr_area_res_df[corr_area_res_df["note"].eq("")].copy()
plot_df["area"] = plot_df["area"].astype(str)

# Custom color palette for 7 regions
region_order = ["BST", "CEA", "LHA", "PAG", "PBN", "PVH", "PVT"]
colors = ["#4D67C8","#A7EA51","#5DCEAF","#702FA0","#FFC001","#FF8021","#5DCBF3"]
region_to_color = dict(zip(region_order, colors))

# Grid layout: 2 × 5
nrows, ncols = 2, 5
fig, axes = plt.subplots(nrows=nrows, ncols=ncols,
                         figsize=(4*ncols, 3*nrows),
                         sharey=False)
axes = axes.flatten()

for i, res in enumerate(resolutions):
    ax = axes[i]
    d = plot_df[plot_df["resolution"] == res]

    if d.empty:
        ax.set_title(res + "\n(no fits)")
        ax.axis("off")
        continue

    # Sort areas by coef for vertical order
    d = d.sort_values("coef").reset_index(drop=True)
    y_positions = np.arange(len(d))

    # Plot effect sizes
    for j, row in d.iterrows():
        area = row["area"]
        color = region_to_color.get(area, "gray")  # default if unexpected region appears
        significant = (row["pval"] < 0.05)

        ax.errorbar(
            x=row["coef"],
            y=j,
            xerr=row["se"] * 1.96,
            fmt='o',
            color=color,
            ecolor=color,
            markersize=8 if significant else 6,
            markeredgecolor='black' if significant else 'none',
            markeredgewidth=1.2 if significant else 0.0,
            zorder=3,
        )

    # Y labels = areas
    ax.set_yticks(y_positions)
    ax.set_yticklabels(d["area"])

    # Reference vertical line at β = 0
    ax.axvline(0, linestyle='--', color='gray', linewidth=1.0)

    ax.set_title(res)
    ax.set_xlabel("β1 (ESMU effect)")
    ax.grid(axis='x', linestyle=':', alpha=0.3)

# Turn off unused axes (if resolutions < 10)
for k in range(len(resolutions), len(axes)):
    axes[k].axis("off")
plt.tight_layout()
plt.savefig("ESmu_clustersize_regression.png", dpi=300)
plt.show()


# In[37]:


# --------- 8. Additional models: mean_total_counts and mean_n_genes outcomes ----------
# We'll fit the same predictors but with different outcomes:
outcomes = ["mean_total_counts", "mean_n_genes"]
outcome_results = []
for outcome in outcomes:
    for res in resolutions:
        df_res = merged[merged["resolution"] == res].copy()
        if df_res[outcome].isna().all() or len(df_res) < 5:
            outcome_results.append({
                "outcome": outcome,
                "resolution": res,
                "note": "no_data_or_too_few"
            })
            continue
        # fit: outcome ~ esmu + log_cluster_size
        summary = fit_model_and_summarize(df_res, outcome=outcome, predictor="esmu", covariate="log_cluster_size")
        if summary is None:
            outcome_results.append({
                "outcome": outcome,
                "resolution": res,
                "note": "fit_failed"
            })
        else:
            outcome_results.append({
                "outcome": outcome,
                "resolution": res,
                "coef": summary["coef"],
                "se": summary["se"],
                "pval": summary["pval"],
                "ci_lower": summary["ci_lower"],
                "ci_upper": summary["ci_upper"],
                "n": summary["n"],
                "r2": summary["r2"],
                "note": ""
            })

outcome_results_df = pd.DataFrame(outcome_results)
outcome_results_df


# In[44]:


# Filter to valid fits
plot_df = outcome_results_df[outcome_results_df["note"] == ""].copy()

outcomes = ["mean_total_counts", "mean_n_genes"]

fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=False)

for ax, outcome in zip(axes, outcomes):
    d = plot_df[plot_df["outcome"] == outcome].copy()
    d = d.sort_values("resolution").reset_index(drop=True)

    # Create gradient colormap
    cmap = plt.cm.viridis
    colors = cmap(np.linspace(0, 1, len(d)))

    # --- Plot each point individually ---
    for j, row in d.iterrows():
        is_sig = (row["pval"] < 0.05)

        ax.errorbar(
            x=row["coef"],
            y=j,
            xerr=1.96 * row["se"],
            fmt="o",
            color=colors[j],
            ecolor=colors[j],
            markersize=8 if is_sig else 6,
            markeredgecolor="black" if is_sig else "none",
            markeredgewidth=1.2 if is_sig else 0.0,
            zorder=3,
        )

    # Add vertical dotted line at 0
    ax.axvline(0, linestyle="--", color="gray", linewidth=1.0)

    # Y-axis resolution labels
    ax.set_yticks(np.arange(len(d)))
    ax.set_yticklabels(d["resolution"])

    ax.set_xlabel("β₁ (ESμ) ± 95% CI")
    ax.set_title(f"Outcome: {outcome}")

plt.tight_layout()
plt.savefig("ESmu_gene_regression.png", dpi=300)
plt.show()


# In[4]:


###Save a new metadata file to generate CELLECT input
# Columns to include in CELLEX meta
meta = adata.obs[["hash.mcl.ID", "detailed_celltype_annotation"]].copy()
meta.to_parquet("cellex_metadata.parquet")
meta


# In[4]:


###This is used to run area-wise CELLEX with leaving out 1 animal at a time
area_to_animals = (
    adata.obs.groupby("area")["hash.mcl.ID"]
    .unique()
    .apply(list)
)
meta = adata.obs[["hash.mcl.ID", "detailed_celltype_annotation", "area"]].copy()
meta.to_parquet("cellex_metadata_area.parquet")
area_to_animals


# In[20]:


###Here we read in those results again
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
combined_LOO


# In[3]:


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


# In[4]:


# Assume combined_LOO is already loaded in the namespace
# columns: ['area','left_out_animal','cluster','esmu']

df = combined_LOO.copy()
df['left_out_animal'] = df['left_out_animal'].astype(str)
df['cluster'] = df['cluster'].astype(str)
df['area'] = df['area'].astype(str)


# In[5]:


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


# In[6]:


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


# In[7]:


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


# In[8]:


# compute mean MAD across all areas for ranking
cluster_instab = (
    stab.groupby("cluster")["mad"]
    .mean()
    .sort_values(ascending=False)
)

top4_clusters = cluster_instab.head(4).index.tolist()
print("Top 4 unstable clusters:", top4_clusters)


# In[9]:


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
    ncols = 4
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

        ax.set_title(area)
        ax.set_xlabel("Left-out animal")
        ax.set_ylabel("ESμ value")
        ax.tick_params(axis='x', rotation=45)

    # Hide any empty subplots
    for j in range(i+1, len(axes)):
        axes[j].axis('off')

    fig.suptitle(f"Instability Diagnostic: {cl}", fontsize=14, y=1.03)
    plt.tight_layout()
    plt.show()


# In[10]:


indir = "LOO_cellex_results"
#Loop through each CSV, pull out the v147 row, and build a per‐area DataFrame
records = []
for fn in os.listdir(indir):
    if not fn.endswith("_ALL.csv"):
        continue
    area = fn.replace("_LOO_animal_ALL.csv","")
    area = area.replace("CELLEX_ESmu_","")
    dfb   = pd.read_csv(os.path.join(indir, fn), index_col=0)
    if "v147" not in dfb.index:
        continue
    # extract the row as a Series of cell‐type → ESMU
    s = dfb.loc["v147"].rename(area)
    records.append(s)

#Combine into one DataFrame: rows=areas, cols=cell‐types
combined = pd.DataFrame(records)
#Make sure columns match your actual adata.obs categories
#celltypes = sorted(adata.obs["detailed_celltype_annotation"].unique())
#combined = combined.reindex(columns=celltypes, fill_value=0)
combined = combined.sort_index()
print(combined)


# In[11]:


# top4_clusters should be your list of clusters to inspect
for cl in top4_clusters:
    df_cl = df[df["cluster"] == cl].copy()
    df_cl["left_out_animal"] = df_cl["left_out_animal"].astype(int)

    # Set up figure: one subplot per area
    n_areas = len(df_cl["area"].unique())
    ncols = 4
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


# In[12]:


df_loo = pd.read_csv('/projects/perslab/people/qbx224/CELLECT/CELLECT-LOO/CELLECT-LDSC/results/prioritization.csv')
df_loo


# In[13]:


# extract the left-out animal from specificity_id
df_loo["left_out_animal"] = df_loo["specificity_id"].str.replace("loo_", "").astype(int)
# compute -log10(p)
df_loo["mlog10p"] = -np.log10(df_loo["pvalue"])
# Define MAD function
def mad(x):
    return np.median(np.abs(x - np.median(x)))

cluster_stats = (
    df_loo.groupby("annotation")
      .agg(
          median_p = ("pvalue", "median"),
          median_mlog10p = ("mlog10p", "median"),
          mad_mlog10p = ("mlog10p", mad),
          n_animals = ("left_out_animal", "nunique"),
          frac_sig_005 = ("pvalue", lambda x: np.mean(x < 0.05)),
          frac_sig_001 = ("pvalue", lambda x: np.mean(x < 0.01))
      )
      .reset_index()
)

cluster_stats


# In[14]:


df_long = df_loo[["annotation", "left_out_animal", "pvalue", "mlog10p"]].copy()
# Sort clusters by instability
plot_df = cluster_stats.sort_values("mad_mlog10p", ascending=False)

plt.figure(figsize=(6, 12))
sns.barplot(
    data=plot_df,
    y="annotation",
    x="mad_mlog10p",
    palette="viridis"
)
plt.xlabel("Instability (MAD of –log10(p))")
plt.ylabel("Cluster")
plt.title("Cluster CELLECT Instability Across Leave-One-Animal-Out Runs")
plt.tight_layout()
plt.savefig("Cluster CELLECT Instability Across Leave-One-Animal-Out Runs.png", dpi=300)
plt.show()


# In[15]:


# get top 7 unstable clusters
top7 = (
    cluster_stats
    .sort_values("mad_mlog10p", ascending=False)
    .head(7)["annotation"]
    .tolist()
)

top7


# In[16]:


# OPTIONAL: dict mapping cluster → full-data -log10(p)
# Example:
# full_data_mlog10p = {"04_Crnde_Lmx1a": 4.2, ...}
full_data_mlog10p = {}

for cl in top7:
    df_cl = df_loo[df_loo["annotation"] == cl].copy()
    df_cl["left_out_animal"] = df_cl["left_out_animal"].astype(str)

    fig, ax = plt.subplots(figsize=(10, 5))

    # jittered points
    sns.stripplot(
        data=df_cl,
        x="left_out_animal",
        y="mlog10p",
        size=6,
        color="#4D67C8",   # a clean distinct color
        ax=ax
    )

    # connecting line
    sns.lineplot(
        data=df_cl,
        x="left_out_animal",
        y="mlog10p",
        estimator=None,
        sort=True,
        color="#4D67C8",
        linewidth=1.2,
        ax=ax
    )

    # reference line (full-data CELLECT)
    if cl in full_data_mlog10p:
        ax.axhline(full_data_mlog10p[cl],
                   color="black",
                   linestyle="--",
                   linewidth=1,
                   label="Full-data –log10(p)")

        ax.legend()

    ax.set_title(f"CELLECT Instability Diagnostic for Cluster: {cl}")
    ax.set_xlabel("Left-out Animal ID")
    ax.set_ylabel("−log10(p)")
    ax.tick_params(axis='x', rotation=90)

    plt.tight_layout()
    plt.show()


# In[17]:


df = pd.read_csv('/projects/perslab/people/qbx224/CELLECT/CELLECT-LOO/CELLECT-LDSC/results/prioritization.csv')
df["mlog10p"] = -np.log10(df["pvalue"])

summary = (
    df.groupby("annotation")["mlog10p"]
      .agg(["mean", "std", "count"])
      .rename(columns={"mean": "mlog10p_mean", "std": "mlog10p_sd"})
)

summary["sem"] = summary["mlog10p_sd"] / np.sqrt(summary["count"])
summary["ci95_low"] = summary["mlog10p_mean"] - 1.96 * summary["sem"]
summary["ci95_high"] = summary["mlog10p_mean"] + 1.96 * summary["sem"]

summary = summary.reset_index()
print(summary.head())


# In[18]:


summary_sorted = summary.sort_values("mlog10p_mean", ascending=False)

fig, ax = plt.subplots(figsize=(14, 8))

ax.errorbar(
    summary_sorted["annotation"],
    summary_sorted["mlog10p_mean"],
    yerr=1.96 * summary_sorted["sem"],
    fmt="o",
    ecolor="gray",
    capsize=3,
    markersize=5
)

ax.set_ylabel("−log10(p) mean (± 95% CI)")
ax.set_xlabel("Cluster")
ax.set_title("Cluster BMI Enrichment (LOO Mean ± 95% CI)")

plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig("Cluster BMI Enrichment (LOO Mean ± 95% CI).png", dpi=300)
plt.show()


# In[19]:


df_cellect = pd.read_csv('/projects/perslab/people/qbx224/CELLECT/CELLECT-LOO/CELLECT-LDSC/results/prioritization.csv')
df_cellect["mlog10p"] = -np.log10(df_cellect["pvalue"])
df_cellect


# In[22]:


df_cellect = df_cellect.copy()
df_cellect["left_out_animal"] = df_cellect["specificity_id"].str.replace("loo_", "").astype(str)
df_cellect = df_cellect.rename(columns={"annotation": "cluster"})
df_merged = df_cellect.merge(
    combined_LOO[["left_out_animal", "cluster", "esmu"]],
    on=["left_out_animal", "cluster"],
    how="inner"
)
df_merged


# In[23]:


model = smf.mixedlm(
    "mlog10p ~ esmu",
    data=df_merged,
    groups=df_merged["left_out_animal"]
).fit()
print(model.summary())


# In[29]:


df_merged = df_cellect.merge(
    combined_LOO[["area", "left_out_animal", "cluster", "esmu"]],
    on=["left_out_animal", "cluster"],
    how="inner"
)
df_merged["left_out_animal"] = df_merged["left_out_animal"].astype(int)
df_merged


# In[30]:


area_to_animals = { "BST": [60, 19, 20, 37, 72, 12, 40, 48, 47], 
                   "CEA": [126, 57, 81, 95, 100, 90, 67, 120, 73], 
                   "LHA": [124, 61, 93, 98, 122, 65], 
                   "PAG": [15, 10, 18, 46, 45, 38, 33, 58, 66], 
                   "PBN": [31, 34, 64, 51, 44, 56, 43, 76, 16], 
                   "PVH": [63, 75, 138, 97, 85, 88, 118, 69], 
                   "PVT": [53, 62, 41, 42, 25, 74, 54, 49], }
area_results = []

for area, animals in area_to_animals.items():

    df_area = df_merged[df_merged["left_out_animal"].isin(animals)].copy()

    # Safety: skip if not enough observations
    if df_area["esmu"].nunique() < 2:
        print(f"Skipping {area}: no variability in ESmu.")
        continue

    # Fit mixed model
    model = smf.mixedlm(
        "mlog10p ~ esmu",
        data=df_area,
        groups=df_area["left_out_animal"]
    ).fit(reml=True)

    area_results.append({
        "area": area,
        "n_obs": len(df_area),
        "n_animals": df_area["left_out_animal"].nunique(),
        "beta_esmu": model.params["esmu"],
        "pvalue_esmu": model.pvalues["esmu"],
        "intercept": model.params["Intercept"],
        "random_effect_var": model.cov_re.iloc[0,0],
        "residual_var": model.scale
    })

    print(f"\n=== {area} ===")
    print(model.summary())


# In[31]:


df_area_results = pd.DataFrame(area_results)
df_area_results


# In[32]:


plt.figure(figsize=(7,5))
sns.barplot(
    data=df_area_results.sort_values("beta_esmu", ascending=False),
    x="beta_esmu", y="area", palette="viridis"
)
plt.axvline(0, color="black", linestyle="--")
plt.xlabel("Effect of ESμ on -log10(p)")
plt.ylabel("Area")
plt.title("Mixed Model Effect Sizes per Area")
plt.tight_layout()
plt.show()


# In[35]:


# ---------------------------------------------------------
# area_betas_df must contain: area, beta1, se, pval
# ---------------------------------------------------------
area_betas_df = df_area_results.copy()
area_betas_df["se"] = np.sqrt(area_betas_df["residual_var"] / area_betas_df["n_obs"])
area_betas_df["beta1"] = area_betas_df["beta_esmu"]
area_betas_df["pval"] = area_betas_df["pvalue_esmu"]
# Area ordering and colors
region_order = ["BST", "CEA", "LHA", "PAG", "PBN", "PVH", "PVT"]
colors = ["#4D67C8","#A7EA51","#5DCEAF","#702FA0","#FFC001","#FF8021","#5DCBF3"]
region_to_color = dict(zip(region_order, colors))

# Filter only these 7 areas, ensure ordering
df = area_betas_df.copy()
df = df[df["area"].isin(region_order)]
df["area"] = pd.Categorical(df["area"], categories=region_order, ordered=True)

# Sort vertically by beta1
df = df.sort_values("beta_esmu").reset_index(drop=True)

# Prepare positions
y_positions = np.arange(len(df))

# ---------------------------------------------------------
# PLOT
# ---------------------------------------------------------
fig, ax = plt.subplots(figsize=(7, 4))

for j, row in df.iterrows():
    area = row["area"]
    color = region_to_color[area]
    sig = (row["pval"] < 0.05)

    ax.errorbar(
        x=row["beta_esmu"],
        y=j,
        xerr=row["se"] * 1.96,
        fmt='o',
        color=color,
        ecolor=color,
        markersize=9 if sig else 7,
        markeredgecolor='black' if sig else 'none',
        markeredgewidth=1.2 if sig else 0.0,
        capsize=3,
        zorder=3,
    )

# Horizontal positions
ax.set_yticks(y_positions)
ax.set_yticklabels(df["area"])

# Reference line at 0
ax.axvline(0, linestyle='--', color='gray', linewidth=1)

# Labels
ax.set_xlabel("β₁ (Effect of ESμ on BMI −log10(p))")
ax.set_title("ESμ → BMI Enrichment: Slope by Area (Single Annotation)")
ax.grid(axis='x', linestyle=':', alpha=0.4)

plt.tight_layout()
plt.savefig("ESμ_BMI Enrichment: Slope by Area.png", dpi=300)
plt.show()


# In[9]:


df_cellect = pd.read_csv('/projects/perslab/people/qbx224/CELLECT/CELLECT-UKBB/CELLECT-LDSC/results/prioritization.csv')
df_cellect["mlog10p"] = -np.log10(df_cellect["pvalue"])
df_cellect


# In[15]:


df_plot= df_cellect[df_cellect["gwas"]=="BMI"]
order = (
    df_plot.groupby('annotation')['mlog10p']
    .median()
    .sort_values()
    .index
    .tolist()
)
plt.figure(figsize=(12, 6))
sns.set_style("whitegrid")
sns.pointplot(
    data=df_plot,
    x='annotation',
    y='mlog10p',
    order=order,
    linestyle="none"
)

# 3) add the intersect line at y=3.1
#plt.axhline(3.1, color='red', linestyle='--', label='y = 3.1')
# polish up
plt.xticks(rotation=45, ha='right')
plt.xlabel('Cell Type')
plt.ylabel('log10(Cellect BMI)')
plt.title('Cellect BMI (log10) by Cell Type')
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
plt.show()


# In[17]:


# -------------------------
# 1. Convert ESmu (combined) to long format
# -------------------------
df_esmu = combined.copy()
df_esmu["area"] = df_esmu.index
df_esmu = df_esmu.reset_index(drop=True)

df_long = df_esmu.melt(id_vars="area",
                       var_name="annotation",
                       value_name="esmu")

# -------------------------
# 2. Merge with CELLECT output
# -------------------------
dfm = df_cellect.merge(df_long, on="annotation", how="inner")
dfm


# In[23]:


# -------------------------
# 3. Function: OLS per GWAS × area
# -------------------------
results = []

for gwas_name in dfm["gwas"].unique():
    df_gwas = dfm[dfm["gwas"] == gwas_name]

    for area in df_gwas["area"].unique():
        d = df_gwas[df_gwas["area"] == area]
        d = d.replace([np.inf, -np.inf], np.nan).dropna(subset=["esmu", "mlog10p"])
        # require at least 3 points to fit
        if len(d) < 3:
            continue

        X = sm.add_constant(d["esmu"])
        y = d["mlog10p"]

        model = sm.OLS(y, X).fit()

        results.append({
            "gwas": gwas_name,
            "area": area,
            "n_annotations": len(d),
            "beta_esmu": model.params["esmu"],
            "beta_esmu_se": model.bse["esmu"],
            "pvalue_esmu": model.pvalues["esmu"],
            "intercept": model.params["const"],
            "r_squared": model.rsquared
        })
# Convert to dataframe
df_regression = pd.DataFrame(results)
df_regression


# In[29]:


df = df_regression.copy()
df["mlog10p"] = -np.log10(df["pvalue_esmu"])
mat = df.pivot(index="gwas", columns="area", values="beta_esmu")

plt.figure(figsize=(10, 8))
sns.heatmap(
    mat,
    cmap="coolwarm",
    center=0,
    annot=True,           # write numbers inside cells
    fmt=".2f",            # formatting
    linewidths=0.5,
    cbar_kws={"label": "β (EμS)"}
)

plt.title("GWAS × Area – Effect Size β")
plt.xlabel("Area")
plt.ylabel("GWAS")
plt.tight_layout()
plt.show()


# In[26]:


df = dfm.copy()
final_corr_table = (
    df
    .groupby(["gwas", "area"])
    .apply(lambda g: g["esmu"].corr(g["mlog10p"]))
    .unstack("area")                         # columns = areas
)

# Add overall collective correlation per resolution
final_corr_table["ALL"] = (
    df
    .groupby("gwas")
    .apply(lambda g: g["esmu"].corr(g["mlog10p"]))
)

final_corr_table


# In[27]:


plt.figure(figsize=(10,5))
sns.heatmap(final_corr_table, annot=True, cmap="coolwarm", center=0)
plt.title("Correlation(ESMU vs -log10(BMI)) per area and per GWAS trait")
plt.savefig("ESmu_GWAS_log10.png", dpi=300)
plt.show()


# In[ ]:




