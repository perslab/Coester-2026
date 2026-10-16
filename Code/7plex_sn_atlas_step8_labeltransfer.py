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


# In[ ]:


#adata.write("neurons_bipca_leiden.h5ad")
adata = sc.read("neurons_bipca_leiden.h5ad")
adata


# In[5]:


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
sc.pl.embedding(adata,
                basis="tsne",
                color=leiden_cols, 
                legend_loc="on data", 
                legend_fontsize="x-small")


# In[69]:


###Load in the virus esmu for the leiden clusters
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


# In[70]:


###Read in the CELLECT results for all leiden resolutions
df = pd.read_csv('CELLECT/CELLECT-leiden/CELLECT-LDSC/results/prioritization.csv')
# rename the column
df = df.rename(columns={"specificity_id": "resolution"})
df = df.rename(columns={"annotation": "cluster"})
# convert values like "leiden_02" -> "leiden_0.2"
df["resolution"] = (
    df["resolution"]
      .str.replace("leiden_", "", regex=False)      # "02"
      .astype(int)                                 # 2
      .map(lambda x: f"leiden_{x/10:.1f}")          # "leiden_0.2", ..., "leiden_2.0"
)
df


# In[71]:


merged = combined_long.merge(
    df,
    on=["resolution", "cluster"],
    how="left",          # keep all area×cluster rows; add beta/pvalue where available
    validate="many_to_one"  # many areas per (resolution,cluster), but only one CELLECT row
)
merged["neglog10_p"] = -np.log10(merged["pvalue"])
merged


# In[72]:


corr_per_area = (
    merged
    .groupby(["resolution", "area"])
    .apply(lambda g: g["esmu"].corr(g["beta"]))
    .reset_index(name="corr")
)
corr_per_resolution = (
    merged
    .groupby("resolution")
    .apply(lambda g: g["esmu"].corr(g["beta"]))
    .reset_index(name="corr")
)
final_corr_table = (
    merged
    .groupby(["resolution", "area"])
    .apply(lambda g: g["esmu"].corr(g["beta"]))
    .unstack("area")                         # columns = areas
)

# Add overall collective correlation per resolution
final_corr_table["ALL"] = (
    merged
    .groupby("resolution")
    .apply(lambda g: g["esmu"].corr(g["beta"]))
)
final_corr_table


# In[73]:


plt.figure(figsize=(10,5))
sns.heatmap(final_corr_table, annot=True, cmap="coolwarm", center=0)
plt.title("Correlation(ESMU vs beta(BMI)) per area and per resolution")
plt.savefig("ESmu_clusters_beta.png", dpi=300)
plt.show()


# In[74]:


corr_per_area = (
    merged
    .groupby(["resolution", "area"])
    .apply(lambda g: g["esmu"].corr(g["neglog10_p"]))
    .reset_index(name="corr")
)
corr_per_resolution = (
    merged
    .groupby("resolution")
    .apply(lambda g: g["esmu"].corr(g["neglog10_p"]))
    .reset_index(name="corr")
)
final_corr_table = (
    merged
    .groupby(["resolution", "area"])
    .apply(lambda g: g["esmu"].corr(g["neglog10_p"]))
    .unstack("area")                         # columns = areas
)

# Add overall collective correlation per resolution
final_corr_table["ALL"] = (
    merged
    .groupby("resolution")
    .apply(lambda g: g["esmu"].corr(g["neglog10_p"]))
)
plt.figure(figsize=(10,5))
sns.heatmap(final_corr_table, annot=True, cmap="coolwarm", center=0)
plt.title("Correlation(ESMU vs neglog10_p(BMI)) per area and per resolution")
plt.savefig("ESmu_clusters_pvalue.png", dpi=300)
plt.show()


# In[75]:


df = merged.copy()

# Z-score
df["z"] = df["beta"] / df["beta_se"]

# BH correction within each (area, resolution)
df["pvalue_BH"] = (
    df.groupby(["area", "resolution"])["pvalue"]
      .transform(lambda p: smm.multipletests(p.values, method="fdr_bh")[1])
)

# safe -log10 transforms
eps = np.nextafter(0, 1)
df["neglog10_p"] = -np.log10(df["pvalue"].clip(lower=eps))
df["neglog10_fdr"] = -np.log10(df["pvalue_BH"].clip(lower=eps))
df


# In[76]:


metric_col = "esmu"
records = []
for (area, res), sub in df.groupby(["area", "resolution"]):
    sub = sub.dropna(subset=[metric_col, "neglog10_p", "z", "beta"])
    if len(sub) < 3:
        continue

    rho_fdr, p_fdr = spearmanr(sub[metric_col], sub["neglog10_p"])
    rho_z,   p_z   = spearmanr(sub[metric_col], sub["z"])
    rho_b,   p_b   = spearmanr(sub[metric_col], sub["beta"])

    records.append({
        "area": area,
        "resolution": res,
        "n_clusters": len(sub),
        "rho_esmu_vs_neglog10_p": rho_fdr,
        "p_esmu_vs_neglog10_p": p_fdr,
        "rho_esmu_vs_z": rho_z,
        "p_esmu_vs_z": p_z,
        "rho_esmu_vs_beta": rho_b,
        "p_esmu_vs_beta": p_b,
    })

corr_area_res = pd.DataFrame(records).sort_values(["area", "resolution"])
corr_area_res


# In[77]:


corr_long = corr_area_res.melt(
    id_vars=["area", "resolution"],
    value_vars=[
        "rho_esmu_vs_neglog10_p",
        "rho_esmu_vs_z",
        "rho_esmu_vs_beta",
    ],
    var_name="comparison",
    value_name="rho",
)

# nicer labels
corr_long["comparison"] = corr_long["comparison"].map({
    "rho_esmu_vs_neglog10_p": "−log10(p)",
    "rho_esmu_vs_z": "Z-score",
    "rho_esmu_vs_beta": "Beta",
})
heatmaps = {
    comp: sub.pivot(index="area", columns="resolution", values="rho")
    for comp, sub in corr_long.groupby("comparison")
}
sns.set(style="white", font_scale=1.0)

fig, axes = plt.subplots(
    1, 3,
    figsize=(18, 6),
    sharey=True,
    constrained_layout=True
)

vlim = 1  # Spearman rho is in [-1, 1]

for ax, (comp, mat) in zip(axes, heatmaps.items()):
    sns.heatmap(
        mat,
        ax=ax,
        cmap="coolwarm",
        vmin=-vlim,
        vmax=vlim,
        center=0,
        annot=True,
        fmt=".2f",
        linewidths=0.5,
        cbar=ax is axes[-1],  # only show colorbar once
        cbar_kws={"label": "Spearman ρ"}
    )
    ax.set_title(f"esmu vs {comp}")
    ax.set_xlabel("Resolution")
    ax.set_ylabel("Area")
plt.savefig("ESmu_clusters_spearman.png", dpi=300)
plt.show()


# In[78]:


comparisons = {
    "−log10(p)": "neglog10_p",
    "Z-score": "z",
    "Beta": "beta",
}

records = []
for (area, res), sub in df.groupby(["area", "resolution"]):
    for label, ycol in comparisons.items():
        tmp = sub.dropna(subset=[metric_col, ycol])
        if len(tmp) < 3:
            rho, p = np.nan, np.nan
        else:
            rho, p = spearmanr(tmp[metric_col], tmp[ycol])

        records.append({
            "area": area,
            "resolution": res,
            "comparison": label,
            "rho": rho,
            "p": p,
            "n": len(tmp),
        })

corr_long = pd.DataFrame(records)
def p_to_stars(p):
    if pd.isna(p):
        return ""
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return ""

corr_long["stars"] = corr_long["p"].apply(p_to_stars)

sns.set(style="white", font_scale=1.0)

# optional: order resolutions numerically
res_order = sorted(corr_long["resolution"].unique(), key=lambda x: float(x.split("_")[1]))
area_order = sorted(corr_long["area"].unique())

fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharey=True, constrained_layout=True)

for ax, comp in zip(axes, comparisons.keys()):
    sub = corr_long[corr_long["comparison"] == comp]

    mat_rho = sub.pivot(index="area", columns="resolution", values="rho").reindex(index=area_order, columns=res_order)
    mat_p   = sub.pivot(index="area", columns="resolution", values="p").reindex(index=area_order, columns=res_order)
    mat_st  = sub.pivot(index="area", columns="resolution", values="stars").reindex(index=area_order, columns=res_order)

    # Make annotation strings like "0.42\n**"
    annot = mat_rho.copy().astype(object)
    for i in range(mat_rho.shape[0]):
        for j in range(mat_rho.shape[1]):
            r = mat_rho.iat[i, j]
            st = mat_st.iat[i, j]
            annot.iat[i, j] = "" if pd.isna(r) else f"{r:.2f}\n{st}"

    sns.heatmap(
        mat_rho,
        ax=ax,
        cmap="coolwarm",
        vmin=-1, vmax=1, center=0,
        linewidths=0.5,
        annot=annot, fmt="",
        cbar=(ax is axes[-1]),
        cbar_kws={"label": "Spearman ρ"} if ax is axes[-1] else None,
    )

    ax.set_title(f"{metric_col} vs {comp}")
    ax.set_xlabel("Resolution")
    ax.set_ylabel("Area")
plt.savefig("ESmu_clusters_spearman.png", dpi=300)
plt.show()


# In[79]:


def corr_pair(x, y):
    return {
        "pearson": pearsonr(x, y)[0],
        "spearman": spearmanr(x, y)[0],
    }

# example for one area × resolution
g = df.query("area == 'CEA' and resolution == 'leiden_1.0'")

{
    "raw_p": corr_pair(g["esmu"], g["neglog10_p"]),
    "fdr": corr_pair(g["esmu"], g["neglog10_fdr"]),
}


# In[88]:


###Start reading in CELLECT esmu values with left-out animals
area_to_animals = {
    "BST": [60, 19, 20, 37, 72, 12, 40, 48, 47],
    "CEA": [126, 57, 81, 95, 100, 90, 67, 120, 73],
    "LHA": [124, 61, 93, 98, 122, 65],
    "PAG": [15, 10, 18, 46, 45, 38, 33, 58, 66],
    "PBN": [31, 34, 64, 51, 44, 56, 43, 76, 16],
    "PVH": [63, 75, 138, 97, 85, 88, 118, 69],
    "PVT": [53, 62, 41, 42, 25, 74, 54, 49],
}
cellect_df = pd.read_csv('/CELLECT/CELLECT-7plex/CELLECT-LDSC/results/prioritization.csv')
cellect_df = cellect_df[cellect_df["specificity_id"]=="detailed"]
# Z-score
cellect_df["z"] = cellect_df["beta"] / cellect_df["beta_se"]

# safe -log10 transforms
eps = np.nextafter(0, 1)
cellect_df["neglog10_p"] = -np.log10(cellect_df["pvalue"].clip(lower=eps))
cellect_df = cellect_df.rename(columns={"annotation": "cluster"})
cellect_df


# In[82]:


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


# In[90]:


df = combined_LOO.merge(
    cellect_df[["cluster", "z", "beta"]],
    on=["cluster"],
    how="left",
    validate="many_to_one",
)
df


# In[91]:


records = []

for area, sub in df.groupby("area"):
    sub = sub.dropna(subset=["esmu", "z", "beta", "left_out_animal"])

    # Need enough data to fit a mixed model
    if sub["left_out_animal"].nunique() < 3:
        continue

    # ---- esmu ~ z ----
    try:
        m_z = smf.mixedlm(
            "esmu ~ z",
            sub,
            groups=sub["left_out_animal"],
            re_formula="1",
        ).fit(reml=False)

        records.append({
            "area": area,
            "predictor": "z",
            "coef": m_z.params["z"],
            "se": m_z.bse["z"],
            "z_stat": m_z.tvalues["z"],
            "pvalue": m_z.pvalues["z"],
            "n_obs": len(sub),
            "n_animals": sub["left_out_animal"].nunique(),
        })
    except Exception as e:
        print(f"[WARN] z-model failed for area={area}: {e}")

    # ---- esmu ~ beta ----
    try:
        m_b = smf.mixedlm(
            "esmu ~ beta",
            sub,
            groups=sub["left_out_animal"],
            re_formula="1",
        ).fit(reml=False)

        records.append({
            "area": area,
            "predictor": "beta",
            "coef": m_b.params["beta"],
            "se": m_b.bse["beta"],
            "z_stat": m_b.tvalues["beta"],
            "pvalue": m_b.pvalues["beta"],
            "n_obs": len(sub),
            "n_animals": sub["left_out_animal"].nunique(),
        })
    except Exception as e:
        print(f"[WARN] beta-model failed for area={area}: {e}")


# In[92]:


mixed_results = pd.DataFrame(records)

# BH correction across areas per predictor
from statsmodels.stats.multitest import multipletests

mixed_results["pvalue_BH"] = (
    mixed_results
    .groupby("predictor")["pvalue"]
    .transform(lambda p: multipletests(p, method="fdr_bh")[1])
)

mixed_results["neglog10_fdr"] = -np.log10(mixed_results["pvalue_BH"].clip(lower=1e-300))

mixed_results


# In[94]:


region_order = ["BST", "CEA", "LHA", "PAG", "PBN", "PVH", "PVT"]
colors = ["#4D67C8","#A7EA51","#5DCEAF","#702FA0","#FFC001","#FF8021","#5DCBF3"]
region_to_color = dict(zip(region_order, colors))

plot_df = mixed_results.copy()
plot_df["area"] = pd.Categorical(plot_df["area"], categories=region_order, ordered=True)

plot_df["ci_low"]  = plot_df["coef"] - 1.96 * plot_df["se"]
plot_df["ci_high"] = plot_df["coef"] + 1.96 * plot_df["se"]

pcol = "pvalue_BH" if "pvalue_BH" in plot_df.columns else "pvalue"

def stars(p):
    if pd.isna(p): return ""
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    if p < 0.05:  return "*"
    return ""

plot_df["stars"] = plot_df[pcol].apply(stars)

def plot_model(pred, title):
    sub = plot_df[plot_df["predictor"] == pred].copy()
    sub = sub.set_index("area").reindex(region_order).reset_index()

    x = np.arange(len(region_order))
    c = [region_to_color[r] for r in sub["area"]]

    plt.figure(figsize=(7, 4))

    # error bars (CI)
    plt.errorbar(
        x=x,
        y=sub["coef"],
        yerr=[sub["coef"] - sub["ci_low"], sub["ci_high"] - sub["coef"]],
        fmt="o",
        capsize=4,
        linewidth=1,
        color="black",
        zorder=2,
    )

    # recolor points by region
    for xi, yi, colr in zip(x, sub["coef"], c):
        plt.scatter(xi, yi, s=70, color=colr, zorder=3)

    # stars
    for xi, yi, st in zip(x, sub["coef"], sub["stars"]):
        if st:
            plt.text(xi, yi, st, ha="center", va="bottom", fontsize=11)

    plt.axhline(0, linestyle="--", linewidth=1)
    plt.xticks(x, region_order, rotation=45, ha="right")
    plt.xlabel("Region")
    plt.ylabel("Fixed-effect slope (coef) ± 95% CI")
    plt.title(title)
    plt.grid(True, axis="y", alpha=0.3)
    plt.tight_layout()
    plt.show()

# separate plots (separate y-axes)
plot_model("z",    "Mixed model: esmu ~ z + (1|left_out_animal)")
plot_model("beta", "Mixed model: esmu ~ beta + (1|left_out_animal)")


# In[95]:


region_order = ["BST", "CEA", "LHA", "PAG", "PBN", "PVH", "PVT"]
colors = ["#4D67C8","#A7EA51","#5DCEAF","#702FA0","#FFC001","#FF8021","#5DCBF3"]
region_to_color = dict(zip(region_order, colors))

plot_df = mixed_results.copy()
plot_df["area"] = pd.Categorical(plot_df["area"], categories=region_order, ordered=True)

plot_df["ci_low"]  = plot_df["coef"] - 1.96 * plot_df["se"]
plot_df["ci_high"] = plot_df["coef"] + 1.96 * plot_df["se"]

pcol = "pvalue_BH" if "pvalue_BH" in plot_df.columns else "pvalue"

def stars(p):
    if pd.isna(p): return ""
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    if p < 0.05:  return "*"
    return ""

plot_df["stars"] = plot_df[pcol].apply(stars)

outdir = "mixedlm_plots"
os.makedirs(outdir, exist_ok=True)

def plot_model(pred, title, outpath):
    sub = plot_df[plot_df["predictor"] == pred].copy()
    sub = sub.set_index("area").reindex(region_order).reset_index()

    x = np.arange(len(region_order))
    c = [region_to_color[r] for r in sub["area"]]

    plt.figure(figsize=(7, 4))

    plt.errorbar(
        x=x,
        y=sub["coef"],
        yerr=[sub["coef"] - sub["ci_low"], sub["ci_high"] - sub["coef"]],
        fmt="o",
        capsize=4,
        linewidth=1,
        color="black",
        zorder=2,
    )

    for xi, yi, colr in zip(x, sub["coef"], c):
        plt.scatter(xi, yi, s=70, color=colr, zorder=3)

    for xi, yi, st in zip(x, sub["coef"], sub["stars"]):
        if st:
            plt.text(xi, yi, st, ha="center", va="bottom", fontsize=11)

    plt.axhline(0, linestyle="--", linewidth=1)
    plt.xticks(x, region_order, rotation=45, ha="right")
    plt.xlabel("Region")
    plt.ylabel("Fixed-effect slope (coef) ± 95% CI")
    plt.title(title)
    plt.grid(True, axis="y", alpha=0.3)
    plt.tight_layout()

    plt.savefig(outpath, format="pdf", bbox_inches="tight")
    plt.close()

# Save separate PDFs
plot_model("z",    "Mixed model: esmu ~ z + (1|left_out_animal)",
           os.path.join(outdir, "mixed_model_esmu_vs_z.pdf"))

plot_model("beta", "Mixed model: esmu ~ beta + (1|left_out_animal)",
           os.path.join(outdir, "mixed_model_esmu_vs_beta.pdf"))

print("Saved:")
print(" -", os.path.join(outdir, "mixed_model_esmu_vs_z.pdf"))
print(" -", os.path.join(outdir, "mixed_model_esmu_vs_beta.pdf"))


# In[ ]:




