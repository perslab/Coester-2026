# Analysis code: Coester et al. (2026)

**Hypothalamic projection architecture links genetic risk for common obesity to *Lmx1a/b* neurons**

Bernd Coester, Charlotte Høy Kruse, Oliver Christensen, Mathias Mathiasen, Dylan Belmont-Rausch, Tune H. Pers

This folder contains the analysis code for the manuscript. This document maps each section of the Methods and each figure panel to the script that produced it. Not everything in this code has been used for publication, and some analyses contain conceptual errors that were not carried forward but also not deleted. Upstream steps performed by core facilities or commercial partners are listed in Section 8. Code version: 2026-10-16. This overview was generated with the help of AI and manually reviewed. Please reach out if you have specific questions.

------------------------------------------------------------------------

## Contents

1.  [How to read this folder](#1-how-to-read-this-folder)
2.  [Computational environments](#2-computational-environments)
3.  [Workflow overview](#3-workflow-overview)
4.  [Methods to code](#4-methods-to-code)
5.  [Figures to code](#5-figures-to-code)
6.  [Supplementary tables to code](#6-supplementary-tables-to-code)
7.  [Execution order](#7-execution-order)
8.  [Steps not contained in this folder](#8-steps-not-contained-in-this-folder)
9.  [Files not used for the manuscript figures](#9-files-not-used-for-the-manuscript-figures)
10. [Data availability](#10-data-availability)
11. [Use of generative AI](#11-use-of-generative-ai)

------------------------------------------------------------------------

## 1. How to read this folder

All analyses were performed interactively in Jupyter notebooks on a compute cluster. Each `.py` file in this folder is an export of one notebook (`jupyter nbconvert --to script`). Cell boundaries are marked as `# In[n]:`, where `n` is the execution count of that cell, so cell numbers are not sequential. Notebook outputs are not contained in the exports as they take up too much space.

Three file types occur:

| Type | Files | Language |
|------------------------|------------------------|------------------------|
| Python notebooks | `7plex_sn_atlas_step*.py`, `Cellex_7plex.py`, `Imagej_*.py`, `7plex_SABLE_*.py` and others | Python |
| Python notebooks with R cells | `7plex_sn_atlas_step0_preprocessing.py`, `step2`, `step6`, `step10`, `step11`, `2plex_step1_anterograde.py`, `7plex_atlas_spatial_squidpy.py` | Python; scran is called via rpy2 (`%%R` cells) |
| R notebooks which were ran on an R kernel in ipynb | `7plex_spatial_atlas_step1` to `step5`, `2plex_Seurat_ImportToAnnData.R` | R |
| R scripts (batch jobs) | `*.R`, `*.Rmd` | R |

The two notebooks that render the final figure panels are located in the repository root:

| Notebook | Kernel | Content |
|------------------------|------------------------|------------------------|
| `8plex_sn_atlas_figures.ipynb` | `bipca_env` (Python) | Figures 1, 3, 4, 5a and Supplementary Figures 1 to 6; Supplementary Tables |
| `8plex_spatial_atlas_figures.ipynb` | R 4.4.0 | Figures 2, 4c, 4h, 4i, 6b and Supplementary Figure 2b |

The prefix "7plex" refers to the seven retrograde tracing experiments that form the atlas. The prefix "8plex" in the figure notebooks refers to the same atlas together with the dorsal vagal complex (DVC) experiment (see Section 4.6). Please note that this last tracing (DVC retro/2 v36) was added on from a larger pilot experiment where we tested if anterograde tracing with AAV1 would also be possible with this method. It is not. As a workaround, we split off the DVC part and renamed the AAV gene to v147 just for convenience when running CELLEX.

To summarize, the sequencing experiments were done in three batches: 4-plex was the initial batch of retrograde tracings from PAG, PBN, PVT and BST. 3-plex was the second batch with CEA, LHA, and PVH. The 3-plex sequencing run was modified by spiking in AAV-specific primers into the library prep. This resulted in a higher signal for v147, but it has amplified the ambient RNA too. CELLEX can deal with this quite elegantly while other methods have resulted in a clear batch effect between 4-plex and 3-plex. If you do your own analysis, be acutely aware of this. Finally, 2plex was an experiment where we switched to a single-complementary AAV construct and combined an v36-retro/2 injection into the NTS/DVC with anterograde pilots using AAV-1 injections in the NTS/DVC, the VMH and the Mammillary area. The latter AAV-1 injections have yielded almost no signal at all, only the retrograde experiment was carried forward.

Hence, the bulk of the analysis was done on the combined 4-plex plus 3-plex data, while the last DVC piece was just added on but not fully integrated. The reason is that the AAVs we used are technically different and might not be 100% congruent in the results they produce.

As a takeaway, the targeted amplification of captured transcripts is somewhat useless, both in the library prep and also post hoc through amplification of the leftover cDNA. If you aim to refine this approach, implement the targeted capture of transcripts instead or use a probe-based approach.

------------------------------------------------------------------------

## 2. Computational environments

Environment files (YAML and pip requirements) are provided in `Conda_envs/`.

| Environment file | Conda name | Used for | Main packages |
|------------------|------------------|------------------|------------------|
| `bipca.yml`, `bipca.txt` | `pca_env` | snRNA-seq processing, BiPCA, clustering, figures | Python 3.10, scanpy 1.11.0, anndata 0.11.3, bipca 0.1.dev948, openTSNE 1.0.2, leidenalg 0.10.2, rpy2 3.5.17, R 4.3.3 with scran |
| `R_Seurat.yml`, `R_Seurat.txt` | `env_R4.4.0` | Seurat integration, spatial data, label transfer, HYPOMAP mapping | R 4.4, Seurat 5, scran 1.34 |
| `CELLEX.yml`, `CELLEX.txt` | `cellex_env` | Expression specificity (ESµ) | Python 3.10, CELLEX 1.2.2 |
| `CELLECT.yml`, `CELLECT.txt`, `cellect_env.yaml` | `CELLECT_main` | CELLECT-LDSC and CELLECT-MAGMA | Python 3.8, snakemake 6.0.5 |
| `Imaging.yml`, `Imaging.txt` | `imagej` | Fluorescence microscopy processing | Python 3.13, bioio 3.0, scikit-image 0.25.2 |

Batch R jobs were submitted with `rscript.sh` (SLURM, 1 CPU, 1200 GB memory). This is necessary because the Label transfer utility in Seurat is not optimized for large objects and single-threaded.

------------------------------------------------------------------------

## 3. Workflow overview

Note: Several notebooks write an intermediate object once and the corresponding line was afterwards commented out (for example `# adata.write("neurons_bipca_leiden.h5ad")`).

``` mermaid
flowchart TD
    subgraph SN["Single-nucleus RNA-seq"]
        A0["Cell Ranger + hashtag demultiplexing<br/>(core facility)"] --> A1["Tapseq_import.Rmd<br/>Seurat import, rPCA integration"]
        A1 --> A2["step0_preprocessing<br/>BiPCA, Leiden, scran"]
        A2 --> A3["step1_annotation<br/>coarse annotation"]
        A3 --> A4n["step0_neurons"] & A4g["step0_glia"]
        A4n --> A5["step1_annotation<br/>merge, 43 neuronal + 13 glial types"]
        A4g --> A5
        A5 --> A6["step2_reclustering_all_genes<br/>full-genome object with v147"]
        A6 --> A7["step3_validation_labeling<br/>batch QC, exports for Seurat"]
        A7 --> A8["step4_downstream_analysis<br/>per-area count tables"]
    end
    subgraph SP["Spatial transcriptomics"]
        B0["Molecular Cartography, Cellpose<br/>(segmentation upstream)"] --> B1["resolve_preprocessing.R<br/>SCTransform per section"]
        B1 --> B2["resolve_bipca_labeltransfer.R<br/>CCA label transfer"]
        B2 --> B3["spatial step1 to step3<br/>merge, scran, viral probes"]
        B3 --> B4["spatial step4 to step5<br/>plots, benchmarks"]
    end
    A6 --> R1["7plex_Bipca_Seurat.py<br/>98-gene reference"]
    R1 --> B2
    B3 --> R2["resolve_bipca_labeltransfer_reverse.R<br/>regions to nuclei"]
    R2 --> A7
    A8 --> C1["Cellex_7plex.py<br/>ESµ per area and per spatial run"]
    B3 --> C1
    A6 --> C1
    C1 --> C2["CELLECT-LDSC / MAGMA<br/>(snakemake)"]
    C2 --> F1["step9_paper, step12<br/>combined projection scores"]
    A6 --> H1["7plex_Seurat_hypomap.py<br/>HYPOMAP label transfer"]
    H1 --> H2["step5_hypomap, step6_subset<br/>cross-species analyses"]
    D0["2plex DVC experiment<br/>step10_DVC"] --> F1
    F1 --> FIG["8plex_*_figures.ipynb"]
    H2 --> FIG
    B4 --> FIG
```

------------------------------------------------------------------------

## 4. Methods to code

Each subsection corresponds to a Methods heading in the manuscript. Parameters are given as set in the code.

### 4.1 Single-nucleus data import and quality control

The core facility has provided objects that were demultiplexed based on hashed nuclei. These cell barcode identities were carried forward while we re-ran Cell Ranger on the FASTQ with a manually modified mouse genome that includes the AAV reference. The alignment of AAV reads was checked on the Interactive Genome Viewer with the BAM files and they nicely pile up at the 3' end and also correctly detect the splicing site.

| Methods section | Script | Parameters in code | Output |
|------------------|------------------|------------------|------------------|
| Read alignment, hashtag demultiplexing | Not contained (Section 8) | Cell Ranger with mm10 v2.0 reference plus AAV sequence; MCL demultiplexing | Per-lane Seurat objects (`*_rna-seurat.rds`) |
| Import of lanes, integration | `Tapseq_import.Rmd` | Barcodes intersected with core facility objects; `CreateSeuratObject(min.cells = 10)`; rPCA anchors, 50 dims, reference lanes 1, 2, 7, 10, 14, 18, 19 | `7plex_integrated.qs` |
| Gene filtering, re-integration | `7plex_raw.R` | Restriction to 16,934 genes of `filtered_gene_names.txt` | `7plex_raw_filtered.qs` |
| Transfer of original cluster labels | `7plex_rerun.py` (R) |  | `7plex_original_seurat_clusters.csv`, `7plex_metadata.csv` |
| Preprocessing and QC | `7plex_sn_atlas_step0_preprocessing.py` | `stabilize_matrix(threshold = 100)` on cells and genes | `raw_filtered_bipca.h5ad` |

### 4.2 Dimensionality reduction, clustering and annotation

For clustering in the iterative labelling approach, we deliberately removed the pseudogenes to avoid induced variation that is biologically not interpretable. Additionally, the BIPCA pipeline removes very sparsely expressed genes `stabilize_matrix(threshold = 100)`. The full matrix is later added back to allow looking up any gene or psudogene of interest, just note that they are not used for cluster annotation. Cluster annotations are strictly based on the most specific marker genes and not by literature, this means that for example "Pomc" neurons are denoted as "Tbx3_Npy1r". This is by design to allow a first pass analysis of cluster-based projections without being biased by prior literature, as this dataset was not created to define neuronal subtypes but rather to stratify their projection signal across several clustering resolutions.

| Methods section | Script | Parameters in code | Output |
|------------------|------------------|------------------|------------------|
| Dimensionality reduction with BiPCA | `step0_preprocessing`, `step0_neurons`, `step0_glia` | `bipca.BiPCA(n_components = -1, seed = 42)`; rank from Marchenko-Pastur fit | `Z_biwhite` layer, `pcs` |
| Neighborhood graph and clustering | same | `sc.pp.neighbors(use_rep = "pcs", method = "gauss", n_pcs = rank)`; `sc.tl.leiden(resolution = 1, random_state = 42, flavor = "leidenalg", n_iterations = -1)`; `openTSNE.TSNE(n_jobs = 32)`, default perplexity | `leiden`, `tsne` |
| Normalization | same | Preliminary clustering: `normalize_total`, `log1p`, PCA (15 PCs), Leiden (igraph); `scran::computeSumFactors(min.mean = 0.1)`; counts divided by size factors, `log1p` | `scran_normalization` layer |
| Differential expression and markers | `7plex_sn_altas_step1_annotation_celltypes.py`, `step0_neurons`, `step0_glia` | Wilcoxon on `scran_normalization` (`key_added = "dea_leiden"`); `filter_rank_genes_groups(min_in_group_fraction = 0.2, max_out_group_fraction = 0.2)`; log2FC ≥ 2 (all nuclei), ≥ 1.5 (neurons), ≥ 3 (glia) | `dea_*_filtered.csv` |
| Iterative annotation | same | Coarse round: clusters labelled Neuron, Glia or dropped (mitochondrial, immune, histone, unknown); fine round: `Glia-like` and `Neuron-like` clusters removed | `neurons_bipca_annotated.h5ad`, `glia_bipca_annotated.h5ad`, `filtered_bipca_full_annotation.h5ad` |
| Full-genome reclustering | `7plex_sn_atlas_step2_reclustering_all_genes.py` | Complete transcriptome including `v147`, restricted to annotated barcodes; same BiPCA, Leiden and t-SNE settings; scran recomputed; labels mapped by barcode | `allgenes_bipca_annotated.h5ad`, `allgenes_counts.pkl` |

The t-SNE embeddings shown in Figure 1b, Figure 4a and Supplementary Figure 2a are those of the full-genome object (`allgenes_bipca_annotated.h5ad`). Cell type labels in these panels were assigned on the filtered gene set (Section 4.1) and transferred by barcode.

### 4.3 Batch assessment and cell type quality control

Nothing was excluded by statistical results here but rather inspected to see whether some clusters are exclusively populated by a few or single animals. Based on the quality assessment, the lack of specific marker genes and the imbalance in animal x cluster representation, "01_Unassigned" was considered an uninformative cluster. It is largely retained in the dataset but excluded from most high-level analyses.

| Methods section | Script | Parameters in code |
|------------------------|------------------------|------------------------|
| Batch assessment | `7plex_sn_atlas_step3_validation_labeling.py` | `scipy.stats.chi2_contingency` on animal × cell type tables; standardized residuals (O − E)/√E; two-tailed normal p-values; Benjamini-Hochberg; flag at FDR \< 0.001. Repeated by injection area and by sex. |
| Per-animal composition | `7plex_sn_atlas_step4_downstream_analysis.py` | Proportions per animal and z-scores per cell type |

### 4.4 Spatial transcriptomics

Spatial injections are unilateral. This information is retained and special care was taken to exclude possible spillover from injection sites into the analysis of projections by manually inspecting and excluding these regions.

| Methods section | Script | Parameters in code | Output |
|------------------|------------------|------------------|------------------|
| Cell segmentation, transcript assignment | Not contained (Section 8) | Cellpose on merged DAPI and poly(T); transcripts assigned with `sf` | `*object_projection_call.qs` |
| Per-section preprocessing | `resolve_preprocessing.R` | Viral probes removed from the expression assay; `SCTransform`; 30 PCs | `*_transformed.qs` |
| 98-gene single-nucleus reference | `7plex_Bipca_Seurat.py` (R) | Single-nucleus object restricted to the panel genes; `SCTransform` | `7plex_bipca_n98.qs` |
| Label transfer (nuclei to spatial cells) | `resolve_bipca_labeltransfer.R` | `FindTransferAnchors(reduction = "cca", normalization.method = "SCT")`; `TransferData(weight.reduction = "cca")` | `*_bipca.qs` |
| Merging and scran normalization | `7plex_spatial_atlas_step1_seurat.py`, `step2_scran.py` (R) | Sections merged; size factors with transferred labels as clusters (`min.mean = 0.1`); cells with size factor \> 20 removed | `spatial_24_bipca_annotated_scran.qs`, per-run objects |
| Prediction score threshold | `8plex_spatial_atlas_figures.ipynb`, `step3_virus.py` | Label kept if `score_detailed_celltype_annotation > 0.4`, otherwise `no_label` |  |
| Region assignment (spatial to nuclei) | `resolve_bipca_labeltransfer_reverse.R` | CCA transfer with spatial cells as reference; transferred field `region.collapsed` | `7plex_bipca_n98_region_labeltransfer.qs`, `7plex_region_labeltransfer_meta.csv` |
| Viral probe detection | `7plex_spatial_atlas_step3_virus.py` (R) | Probes `p147-mCherry-2A-iCre` and `p36-NLS-Cre`; scran counts and labels exported per run | `area_pickles/spatial{1,2,3}_*.pkl` |

Assignment of viral probes to target regions: run 1, v36 = PAG and v147 = BST; run 2, v36 = PBN and v147 = LHA; run 3, v36 = CEA and v147 = PVT.

### 4.5 Projection identity (ESµ of viral transcripts)

CELLEX works best on the different noise levels in this dataset. The expression specificity scores of AAVs is used as a proxy for projection strength. This approach is semi-quantitative and will inevitable exclude weak signals even if they are biologically relevant.

| Methods section | Script | Parameters in code | Output |
|------------------|------------------|------------------|------------------|
| Per-area count tables (nuclei) | `step4_downstream_analysis.py` | `scran_normalization` layer split by `area` | `area_pickles/{area}_*.pkl` |
| CELLEX per area and per spatial run | `Cellex_7plex.py` | `cellex.ESObject(normalize = False)` on scran counts; annotation `detailed_celltype_annotation` | `area_cellex_results/*_detailed_cellex.csv` |
| Agreement between datasets, threshold 0.36 | `step9_paper.py`, figure notebook cell 19 | Cohen's kappa, MCC and F1 over thresholds 0 to 0.5; maximum kappa at 0.36 | Supplementary Figure 1 |
| Combined projection score | `step9_paper.py`, `step12_merged_experiments.py`, figure notebook cell 10 | Mean of snRNA-seq and spatial ESµ; PVH and DVC from snRNA-seq only; call at ESµ \> 0.36 | Figure 3b, Figure 5a |
| Leave-one-animal-out stability | `Cellex_7plex.py`, `step12_merged_experiments.py` | ESµ recomputed with each animal removed | `looa_esmu_results/` |

### 4.6 Standalone tracing experiments

| Experiment | Script | Content |
|------------------------|------------------------|------------------------|
| DVC retrograde (scAAV-retro2-v36), May 2023 | `Tapseq_import_new.Rmd`, `2plex_Seurat_ImportToAnnData.R`, `7plex_sn_atlas_step10_DVC.py` | Import, CCA label transfer from the atlas, scran normalization, enrichment of labelled nuclei |
| Merge of DVC nuclei with the atlas | `7plex_sn_atlas_step12_merged_experiments.py` | `v36` counts summed into `v147`; nuclei concatenated with the atlas; scran layers aligned to the atlas gene space |
| Anterograde AAV1 test cohort | `2plex_step1_anterograde.py` | Not used in the manuscript |
| LPBN retrograde | `7plex_sn_atlas_step11_LPBN.py` | Excluded from all figures (`dataset != "LPBN retrograde (v36)"`) |

### 4.7 GWAS enrichment

| Methods section | Script | Parameters in code |
|------------------------|------------------------|------------------------|
| Orthologue mapping | `Cellex_7plex.py` | Ensembl BioMart one-to-one orthologues (`mart_homologues.txt`); JAX `HOM_MouseHumanSequence.rpt.txt`; Ensembl IDs from the Pool lab GRCh38 v2.0 annotation |
| CELLEX for CELLECT | `Cellex_7plex.py` | Scran counts, detailed annotation |
| CELLECT-LDSC and MAGMA | Not contained (Section 8) | Results read from `CELLECT-ALLPLEX/CELLECT-LDSC/results/prioritization.csv` |
| Bonferroni correction (mouse) | Figure notebook cell 10 | GWAS `BMI_Loh2018_UKBB`; threshold 0.05/41 neuronal cell types |
| Resolution sensitivity | `7plex_sn_atlas_step7_recluster.py`, `step8_labeltransfer.py` | Leiden 0.2 to 2.0 in steps of 0.2 on neurons; CELLEX and CELLECT per resolution; Pearson correlation between cluster ESµ (`v147`) and mean BMI −log10(p) |

### 4.8 Cross-species analyses

| Methods section | Script | Parameters in code |
|------------------------|------------------------|------------------------|
| HYPOMAP preprocessing | `7plex_sn_atlas_step5_hypomap.py` | Raw counts from `human_HYPOMAP_snRNASeq.h5ad`; slim object and normalized layer |
| Orthologue subsetting, SCTransform | `7plex_Seurat_hypomap.py` (R) | One-to-one orthologues with non-empty symbols; `SCTransform` per species |
| Label transfer | `7plex_Seurat_hypomap.py`, `hypomap_bipca_labeltransfer.R` | HYPOMAP as reference, mouse as query; CCA; levels C1 to C4, region and cell type |
| Filtering of transferred labels | Figure notebook cell 14 | Two-component Gaussian mixture on `score_C4_HM` (`random_state = 42`); cells retained at P(signal) ≥ 0.5 (threshold 0.421) |
| Genuine and contaminant C4 clusters | `step5_hypomap.py` | ≥ 50 % of cells with a given C4 label from the four mouse *Lmx1a* types; human *LMX1A* mean normalized expression \> 0.5 and \> 10 % of cells expressing |
| Human CELLEX, CELLECT | `Cellex_7plex.py` (`hypomap_c3_cellex.csv`, `hypomap_c4_cellex.csv`); CELLECT not contained | GWAS `BMI_Yengo2018`; Bonferroni 0.05/415 (C4) and 0.05/156 (C3) |
| Effector genes | Figure notebook cells 18 and 30 | CELLECT-GENES output; ESµ ≥ 0.3 and MAGMA p \< 0.05 |
| *Lmx1a* subset reclustering | `7plex_sn_atlas_step6_subset.py` | Nuclei of clusters 04, 27, 32 and 36; BiPCA and Leiden as in 4.2; sex bias as the mean difference of male and female Pearson residuals |

### 4.9 Fluorescence microscopy

| Methods section | Script | Parameters in code |
|------------------------|------------------------|------------------------|
| Per-channel thresholding and normalization | `Imagej_LMX_Cohort2.py`, `Imagej_LMX_Cohort2_tracing.py`, `Imagej_LMX_Cohort3.py` | Otsu lower bound; upper bound at the 99.8th to 99.999th percentile of foreground pixels; linear rescaling to 8 bit |
| Batch-consistent thresholding | `Imagej_LMX_thresholding.py`, `Imagej_LMX_Cohort2_tracing.py` | Cohort median of per-image thresholds; outliers \> 10 × MAD flagged |
| Multi-channel compositing | `Imagej_LMX_Cohort2.py` | Screen blending 1 − ∏(1 − I_c · color_c); logistic soft edge with k = 10 |

### 4.10 Steps without code in this folder

Light-sheet segmentation and registration (TrailMap, Perens atlas), body weight statistics and DEXA (GraphPad Prism) are described in Section 8.

------------------------------------------------------------------------

## 5. Figures to code

"sn" denotes `8plex_sn_atlas_figures.ipynb`, "sp" denotes `8plex_spatial_atlas_figures.ipynb`. Cell indices count code cells from 0. Panels that are schematics were drawn in Inkscape. Layout and labels of all panels were finalized in Inkscape.

| Panel | Notebook cell | File written | Input objects |
|------------------|------------------|------------------|------------------|
| 1a | Schematic |  |  |
| 1b | sn 3, 4 | `Figure_1b.pdf`, `Figure_1b_legend.pdf` | `allgenes_bipca_annotated.h5ad` |
| 1c | sn 7, 8, 9 | `Figure_1c.pdf`, `Figure_1c_legend.pdf` | `atlas_merged_NTS_LPBN_retro_scran.h5ad`, `area_cellex_results/` |
| 2a | Schematic |  |  |
| 2b, left | sp 3 | `Figure_2b.pdf` | `spatial_24_bipca_downstream.qs`, section S3B2 |
| 2b, right | sp 4 | `Figure_2c.pdf` | same |
| 3a to 3c | sn 10 | `Figure_3.pdf` | `area_cellex_results/`, `df_cellect.csv`, `celltype_region_enrichment_all.csv` |
| 4a | sn 5 | `Figure_4a.pdf` | `allgenes_bipca_annotated.h5ad` |
| 4b | sn 6 | `Figure_4b.pdf` | `Cellex_7plex_detailed_annotation.csv` |
| 4c | sp 5 | `Figure_4c.pdf` | section S3C2 |
| 4d | sn 14, 15 | `LMX_C4_sankey_GMM_filtered.pdf` | `LMX_bipca.h5ad` |
| 4e | sn 16 | `Figure_4f.pdf` | `hypomap_slim_norm.h5ad` |
| 4f | sn 17 | `Figure_C4_248_411_bars.pdf` | `hypomap_cellect_unified.csv` |
| 4g | sn 18 | `effector_crossspecies_04_C4248_shared.pdf` | CELLECT-GENES effector tables |
| 4h | sp 6 | `Figure_4h.pdf` | section S3C2 |
| 4i | sp 7 | `coexpr_Lmx1b_Bdnf_top10_normalized.pdf` | raw counts, \> 2 transcripts per gene |
| 5a | sn 11 | `Figure_4d.pdf` | output of sn 10 |
| 5b, 5c | Light-sheet pipeline (Section 8) |  |  |
| 6a | Schematic |  |  |
| 6b | sp 8 | `Figure_6b.pdf` | section S2C2 |
| 6c to 6h | `Imagej_LMX_Cohort2.py`, `Imagej_LMX_Cohort2_tracing.py` | PNG composites | CZI slide scans |
| 6i, 6j | GraphPad Prism (Section 8) |  |  |
| S1 | sn 19 | `Supplementary_Figure_1.pdf` | `threshold_optimization_single.csv`, `threshold_optimization_grid.csv` |
| S2a | sn 21 | `Supplementary_Figure_2a.pdf` | `allgenes_bipca_annotated.h5ad` |
| S2b | sp 9 | `Supplementary_Figure_2b.pdf` | section S3B2 |
| S3a | sn 23 | `Supplementary_Figure_3a.pdf` | `neurons_bipca_leiden.h5ad` |
| S3b | sn 24 | `Supplementary_Figure_3b.pdf` | `leiden_cellex_results/`, per-cell BMI −log10(p) |
| S4a | sn 25 | `Supplementary_Figure_4a.pdf` | `allgenes_bipca_annotated.h5ad` |
| S4b | sn 26 | `Supplementary_Figure_4b.pdf` | linear mixed model `frac ~ C(region) + (1 | animal)`; ≥ 10 nuclei per animal and region, ≥ 3 animals per region |
| S4c | sn 28 | `Supplementary_Figure_4c.pdf` | `LMX_bipca.h5ad` |
| S4d | sn 29 | `Supplementary_Figure_4d.pdf` | `LMX_bipca.h5ad` |
| S5 | sn 30 | `Supplementary_Figure_5.pdf` | CELLECT-GENES effector tables |
| S6 | sn 31 | `Supplementary_Figure_6.pdf` | `per_sample_region_density.csv`, `annotation_info.csv` |
| S7 | GraphPad Prism (Section 8) |  |  |

------------------------------------------------------------------------

## 6. Supplementary tables to code

| Table | Notebook cell | File written |
|------------------------|------------------------|------------------------|
| Dataset overview (nuclei, cell types, spatial cells) | sn 32 to 37 | `dataset_summary_stats.xlsx` |
| Overview, per sample, per cell type, *Lmx1a* subset | sn 38 | `suppl_table1_overview.csv` to `suppl_table5_lmx1a_per_celltype.csv` |

The per-sample table lists, for each animal and lane, the injection target, sex, number of nuclei, median UMI and gene counts and the percentage of AAV-positive nuclei.

------------------------------------------------------------------------

## 7. Execution order

Steps are listed in the order in which they consume each other's outputs.

**snRNA-seq atlas** (`pca_env`, R steps in `env_R4.4.0`) 1. `Tapseq_import.Rmd` 2. `7plex_raw.R`, `7plex_rerun.py` 3. `7plex_sn_atlas_step0_preprocessing.py` 4. `7plex_sn_altas_step1_annotation_celltypes.py` (coarse annotation, first part) 5. `7plex_sn_atlas_step0_neurons.py`, `7plex_sn_atlas_step0_glia.py` 6. `7plex_sn_altas_step1_annotation_celltypes.py` (merge, second part) 7. `7plex_sn_atlas_step2_reclustering_all_genes.py` 8. `7plex_sn_atlas_step3_validation_labeling.py` 9. `7plex_sn_atlas_step4_downstream_analysis.py`

**Spatial atlas** (`env_R4.4.0`) 10. `7plex_Bipca_Seurat.py` 11. `resolve_preprocessing.R`, `resolve_bipca_labeltransfer.R` 12. `7plex_spatial_atlas_step1_seurat.py` to `step3_virus.py` 13. `resolve_bipca_labeltransfer_reverse.R`, then `step3_validation_labeling.py` (import of region labels) 14. `7plex_spatial_atlas_step4_plotting.py`, `step5_paper.py`

**Projection identity and GWAS** (`cellex_env`, `CELLECT_main`) 15. `Cellex_7plex.py` 16. CELLECT (Section 8) 17. `7plex_sn_atlas_step7_recluster.py`, `step8_labeltransfer.py`

**Cross-species** (`env_R4.4.0`, `pca_env`) 18. `7plex_sn_atlas_step5_hypomap.py` (preprocessing part), `7plex_Seurat_hypomap.py`, `hypomap_bipca_labeltransfer.R` 19. `7plex_sn_atlas_step5_hypomap.py` (analysis part), `7plex_sn_atlas_step6_subset.py`

**Standalone experiments and merge** (`pca_env`) 20. `Tapseq_import_new.Rmd`, `2plex_Seurat_ImportToAnnData.py`, `7plex_sn_atlas_step10_DVC.py` 21. `7plex_sn_atlas_step12_merged_experiments.py` 22. `7plex_sn_atlas_step9_paper.py`

**Figures** 23. `8plex_sn_atlas_figures.ipynb`, `8plex_spatial_atlas_figures.ipynb`

------------------------------------------------------------------------

## 8. Steps not contained in this folder

| Step | Performed by | Output used here |
|------------------------|------------------------|------------------------|
| Cell Ranger alignment, cell calling, hashtag demultiplexing (MCL) | Single-cell Omics Platform, CBMR | `*_filtered_feature_bc_matrix.h5`, `*_rna-seurat.rds` |
| Molecular Cartography imaging, transcript calling | Resolve Biosciences | Transcript coordinates per section |
| Cellpose segmentation, transcript assignment, ROI import | Separate pipeline | `*object_projection_call.qs` |
| CELLECT-LDSC, CELLECT-MAGMA, CELLECT-GENES | CELLECT snakemake workflow | `prioritization.csv`, `*.effector_genes.txt` |
| Light-sheet segmentation (TrailMap), registration (Perens atlas), density quantification | Commercial partner | `per_sample_region_density.csv`, `annotation_info.csv`, density maps |
| Body weight and DEXA statistics | GraphPad Prism | Figure 6i, 6j, Supplementary Figure 7 |

------------------------------------------------------------------------

## 9. Files not used for the manuscript figures

These files are retained for completeness. They document earlier analysis routes or experiments outside the manuscript.

| File | Content |
|------------------------------------|------------------------------------|
| `7plex.Rmd`, `7plex_copy.Rmd`, `7plex_preprocessing.R`, `sn_object_integration.R` | Earlier Seurat-based annotation of neurons (Louvain clustering, `nCount_RNA > 1500`), superseded by the BiPCA workflow |
| `resolve_labels.R` | Earlier label transfer with a previous reference |
| `7plex_filtered.py` | Export of embeddings from an earlier BiPCA run |
| `optimal_cellex_cutoffs.rds` | Per-area ESµ cutoffs from Youden's J (`spatial_step5_paper.py`), superseded by the shared threshold of 0.36 |
| `7plex_atlas_SpatialData_integration.py`, `7plex_atlas_spatial_squidpy.py`, `resolve_to_xenium.py` | Conversion of spatial data to SpatialData/Xenium format for the browser |
| `2plex_step1_anterograde.py`, `7plex_sn_atlas_step11_LPBN.py` | Anterograde AAV1 test cohort and LPBN retrograde data, not included in the manuscript |
| `Imagej_env.py`, `ImageJ_parallel.py`, `ImageJ_panel.py`, `ImageJ_rename.py`, `ImageJ_macro.py`, `Imagej_LMX_Cohort1_FOS.py` | Earlier injection-site panels with fixed display ranges; c-Fos test cohort |
| `7plex_calR_step0_preprocessing.py`, `7plex_SABLE_step0_preprocessing.py`, `7plex_SABLE_step1_analysis.py`, `7plex_sable.py`, `7plex_sable_reassignment.py` | Indirect calorimetry of a further TeTx cohort (Cohort 3) |
| `7plex_SABLE_step2_DREADD_preprocessing.py`, `7plex_SABLE_step3_DREADD_analysis.py` | Chemogenetic activation study with CNO and vehicle (only 6 pilot animals) |

------------------------------------------------------------------------

## 10. Data availability

Raw and processed snRNA-seq and spatial transcriptomics data will be deposited in the Gene Expression Omnibus upon publication. Annotated AnnData and Seurat objects will be provided together with this code. The spatial dataset can be browsed at <https://cbmr-rmpp.shinyapps.io/spatial_hypothalamus_app/>.

------------------------------------------------------------------------

## 11. Use of generative AI

Generative AI was used to write code for downstream analysis and figure generation, to draft method summaries of the notebooks and to compile this document. All code and this document were reviewed by the authors.
