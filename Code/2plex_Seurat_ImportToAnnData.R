#!/usr/bin/env python
# coding: utf-8

# In[1]:


library(dplyr)
library(Seurat)
library(ggplot2)
library(stats)
library(patchwork)
library(scran)
library(scater)
library(BiocParallel)
library(Matrix)
'%!in%' <- function(x,y)!('%in%'(x,y))
options(repr.plot.width = 12, repr.plot.height = 12)
setwd("/maps/projects/perslab/people/qbx224/7plex_paper/2plex")
options(future.globals.maxSize= 50*1024*1024^2)


# In[2]:


godsnot_102 <- c(
  "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
  "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF",
  "#997D87", "#5A0007", "#809693", "#FEFFE6", "#1B4400", "#4FC601", "#3B5DFF",
  "#4A3B53", "#FF2F80", "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92",
  "#FF90C9", "#B903AA", "#D16100", "#DDEFFF", "#000035", "#7B4F4B", "#A1C299",
  "#300018", "#0AA6D8", "#013349", "#00846F", "#372101", "#FFB500", "#C2FFED",
  "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09", "#00489C", "#6F0062",
  "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66", "#885578",
  "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
  "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F",
  "#938A81", "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757",
  "#C8A1A1", "#1E6E00", "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C",
  "#772600", "#D790FF", "#9B9700", "#549E79", "#FFF69F", "#201625", "#72418F",
  "#BC23FF", "#99ADC0", "#3A2465", "#922329", "#5B4534", "#FDE8DC", "#404E55",
  "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C", "#83AB58", "#001C1E",
  "#D1F7CE", "#004B28", "#C8D0F6", "#A3A489", "#806C66", "#222800", "#BF5650",
  "#E83000", "#66796D", "#DA007C", "#FF1A59", "#8ADBB4", "#1E0200", "#5B4E51",
  "#C895C5", "#320033", "#FF6832", "#66E1D3", "#CFCDAC", "#D0AC94", "#7ED379",
  "#012C58"
)
glial_cells_detailed <- c(
  "Astrocyte", "Endothelial_cell", "Ependymal_cell", "Epithelial-like_cell", 
  "Microglia", "OPC", "Oligodendrocyte", "Pericyte", "Pituicyte", 
  "Pituitary_Crhr1", "Pituitary_Gnrhr", "Pituitary_Mtnr1a", "Pituitary_Pappa2", 
  "Pituitary_Pax7", "Pituitary_Rxfp1", "Schwann_cell", "Tanycyte", "Undefined", "VLMC"
)

cells_detailed <- c('00_Hmcn2_Gpr149', '01_Undefined', '02_Fezf1_Lama1', '03_Rai14_Acvr1c',
       '04_Crnde_Lmx1a', '05_Qrfpr_Pdzrn3', '06_Tbx3_Npy1r', '07_Nfib_Six3',
       '08_Lef1_Rspo2', '09_Nfia_Ccbe1', '10_Npy_Agrp', '11_Tafa4_Kank4',
       '12_Ntng2_Onecut1', '13_Pth2r_Otp', '14_Lhx1_Th', '15_Onecut1_Slc32a1',
       '16_Gal_Sox6', '17_Dkk2_Ebf3', '18_Ust_Bcl11b', '19_Bace2_Tspan18',
       '20_Tmem114_Npsr1', '21_Bnc2_Glp1r', '22_Crabp1_Vgll3', '23_Ust_Nkd2',
       '24_Nts_Ecel1', '25_Nxph4_Ebf3', '26_Dach2_Npsr1', '27_Lmx1a_Ttc6',
       '28_Kiss1_Slc38a3', '29_Dnah11_Pth2r', '30_Trh_Sim1', '31_Lhx9_Npvf',
       '32_Gucy2c_Ttc6', '33_Chat_Ghrh', '34_Pla2r1_Sim1', '35_Lhx9_Zic1',
       '36_Ebf2_Crhbp', '37_Qrfp_Qrfprl', '38_Hdc_Prph', '39_Nfib_Nfia',
       '40_Gpr50_Tac2', '41_Otx1_Zic1', '42_Foxb1_Myo5b', 'Astrocyte',
       'Endothelial_cell', 'Ependymal_cell', 'Epithelial-like_cell',
       'Microglia', 'OPC', 'Oligodendrocyte', 'Pericyte', 'Pituicyte',
       'Pituitary_Crhr1', 'Pituitary_Gnrhr', 'Pituitary_Mtnr1a',
       'Pituitary_Pappa2', 'Pituitary_Pax7', 'Pituitary_Rxfp1', 'Schwann_cell',
       'Tanycyte', 'Undefined', 'VLMC', "no_label")

colors_detailed <- c('#ffff00', '#1ce6ff', '#ff34ff', '#ff4a46', 
                      '#008941', '#006fa6', '#a30059', '#ffdbe5', 
                      '#7a4900', '#0000a6', '#63ffac', '#b79762', 
                      '#004d43', '#8fb0ff', '#997d87', '#5a0007', 
                      '#809693', '#6a3a4c', '#1b4400', '#4fc601', 
                      '#3b5dff', '#4a3b53', '#ff2f80', '#61615a', 
                      '#ba0900', '#6b7900', '#00c2a0', '#ffaa92', 
                      '#ff90c9', '#b903aa', '#d16100', '#ddefff', 
                      '#000035', '#7b4f4b', '#a1c299', '#300018', 
                      '#0aa6d8', '#013349', '#00846f', '#372101', 
                      '#ffb500', '#c2ffed', '#a079bf', '#cc0744', 
                      '#c0b9b2', '#c2ff99', '#001e09', '#00489c', 
                      '#6f0062', '#0cbd66', '#eec3ff', '#456d75', 
                      '#b77b68', '#7a87a1', '#788d66', '#885578', 
                      '#fad09f', '#ff8a9a', '#d157a0', '#bec459', 
                      '#456648', '#0086ed', "#242424")

names(colors_detailed) <- cells_detailed


# In[18]:


raw <- qs::qread("HYP_integrated_2plex.qs")


# In[19]:


DimPlot(raw, raster = F, label = T, repel = F, cols = godsnot_102) + NoLegend()


# In[5]:


raw[[]]


# In[6]:


colnames(raw[[]])


# In[7]:


table(raw$sample_id, raw$old.ident)


# In[8]:


raw@assays$RNA


# In[22]:


DefaultAssay(raw) <- "RNA"


# In[23]:


FeaturePlot(raw, features = "v147")


# In[24]:


FeaturePlot(raw, features = "v36")


# In[11]:


ref <- qs::qread("/maps/projects/perslab/people/qbx224/7plex_paper/atlas/7plex_Seurat_bipca_annotated.qs")
ref


# In[12]:


colnames(ref[[]])


# In[13]:


probes <- rownames(raw)
miss <- probes %in% rownames(ref[["RNA"]])
probes <- probes[miss]
print(paste("Keeping", paste0(length(probes)), "genes.."))
write.csv(probes, "HYP_2plex_probes.csv", row.names = FALSE)


# In[14]:


#Subsetting the probes in the sn data
#Pull the raw expression matrix from the original Seurat object containing only the genes of interest
subset_matrix <- ref@assays$RNA$counts[probes,]
#Create a new Seurat object with just the genes of interest
sn_object <- CreateSeuratObject(subset_matrix)
#Pull the identities from the original Seurat object as a data.frame
meta_data <- ref@meta.data[,colnames(ref@meta.data)]
sn_object <- AddMetaData(object = sn_object, metadata = meta_data) # Add the idents to the meta.data slot
sn_object
rm(ref)


# In[15]:


sn_object$nFeature_RNA <- Matrix::colSums(sn_object@assays$RNA$counts > 0)
sn_object$nCount_RNA <- Matrix::colSums(sn_object@assays$RNA$counts)
sn_object<- subset(sn_object, subset = nCount_RNA > 0)
sn_object <- SCTransform(sn_object, verbose = F)
## Transfer neuron and glia labels
raw <- SCTransform(raw, assay = "RNA", verbose = F)


# In[17]:


###Load and transfer
TransferAnchors <- FindTransferAnchors(reference = sn_object, 
                                       query = raw, 
                                       reduction = "cca",
                                       normalization.method = "SCT",
                                       recompute.residuals = F)
#Data to transfer
predictions <- TransferData(anchorset = TransferAnchors, weight.reduction = "cca",
                             refdata = sn_object$detailed_celltype_annotation)
raw$detailed_celltype_annotation <- predictions$predicted.id
raw$score_detailed_celltype_annotation <- predictions$prediction.score.max
####
qs::qsave(raw, "HYP_integrated_2plex_labelled.qs")


# In[3]:


seurat_object <- qs::qread("HYP_integrated_2plex_labelled.qs")
seurat_object


# In[4]:


DimPlot(seurat_object, 
        raster = F, 
        label = T, 
        repel = F, 
        group.by = "detailed_celltype_annotation", 
        cols = colors_detailed) + NoLegend()


# In[5]:


table(seurat_object$sample_id, seurat_object$area)


# In[6]:


library(sceasy)
library(anndata)
library(reticulate)
use_condaenv('pca_env')


# In[7]:


DefaultAssay(seurat_object) <- "RNA"
seurat_object[["SCT"]] <- NULL
seurat_object <- JoinLayers(seurat_object)


# In[11]:


# Convert Assay5 -> Assay3 (legacy)
seurat_object[["RNA"]] <- as(object = seurat_object[["RNA"]], Class = "Assay")

# Verify it's now a legacy assay
class(seurat_object[["RNA"]])  # Should print "Assay", not "Assay5"


# In[12]:


seurat_object@assays$RNA@data <- seurat_object@assays$RNA@counts # Sceasy uses the data layer which causes the normalized count matrix to be exported instead of the raw counts


# In[13]:


sceasy::convertFormat(seurat_object, 
                      from="seurat", 
                      to="anndata", 
                      assay="RNA", 
                      outFile='HYP_integrated_2plex_labelled.h5ad')


# In[3]:


raw <- qs::qread("DVC_integrated_2plex.qs")


# In[4]:


DimPlot(raw, raster = F, label = T, repel = F, cols = godsnot_102) + NoLegend()


# In[5]:


colnames(raw[[]])


# In[6]:


table(raw$sample_id, raw$old.ident)


# In[7]:


raw <- subset(raw, subset = sample_id %!in% c("1_GDF", "2_GDF", "3_GDF"))
table(raw$sample_id, raw$old.ident)


# In[9]:


DefaultAssay(raw) <- "RNA"
FeaturePlot(raw, features = "v63")


# In[10]:


ref1 <- qs::qread("/maps/projects/perslab/people/qbx224/cagrisema_paper/dvc_neurons_mouse.qs")
ref2 <- qs::qread("/maps/projects/perslab/people/qbx224/cagrisema_paper/dvc_glia_mouse.qs")
ref <- merge(ref1, ref2)
rm(ref1, ref2)
ref


# In[13]:


DefaultAssay(ref) <- "mouse_RNA"
ref[["SCT"]] <- NULL


# In[15]:


ref <- SCTransform(ref, assay = "mouse_RNA", verbose = F)
## Transfer neuron and glia labels
raw <- SCTransform(raw, assay = "RNA", verbose = F)


# In[16]:


table(ref$cell.type)


# In[17]:


TransferAnchors <- FindTransferAnchors(reference = ref, 
                                       query = raw, 
                                       reduction = "cca",
                                       normalization.method = "SCT",
                                       recompute.residuals = F)
#Data to transfer
predictions <- TransferData(anchorset = TransferAnchors, weight.reduction = "cca",
                             refdata = ref$cell.type)
raw$detailed_celltype_annotation <- predictions$predicted.id
raw$score_detailed_celltype_annotation <- predictions$prediction.score.max
####
qs::qsave(raw, "DVC_integrated_2plex_labelled.qs")


# In[18]:


DimPlot(raw, 
        raster = F, 
        label = T, 
        repel = F, 
        group.by = "detailed_celltype_annotation", 
        cols = godsnot_102) + NoLegend()


# In[19]:


library(sceasy)
library(anndata)
library(reticulate)
use_condaenv('pca_env')


# In[20]:


DefaultAssay(raw) <- "RNA"
raw[["SCT"]] <- NULL
raw <- JoinLayers(raw)


# In[21]:


# Convert Assay5 -> Assay3 (legacy)
raw[["RNA"]] <- as(object = raw[["RNA"]], Class = "Assay")

# Verify it's now a legacy assay
class(raw[["RNA"]])  # Should print "Assay", not "Assay5"


# In[22]:


raw@assays$RNA@data <- raw@assays$RNA@counts
# Sceasy uses the data layer which causes the normalized count matrix to be exported instead of the raw counts


# In[23]:


sceasy::convertFormat(raw, 
                      from="seurat", 
                      to="anndata", 
                      assay="RNA", 
                      outFile='DVC_integrated_2plex_labelled.h5ad')


# In[24]:


ref


# In[27]:


meta <- ref@meta.data
write.csv(meta, "DVC_mouse_ref_metadata.csv", row.names = TRUE)


# In[ ]:


DefaultAssay(ref) <- "mouse_RNA"
# Keep only the one assay you need, drop everything else
ref_slim <- DietSeurat(
    ref,
    assays  = "mouse_RNA",   # drop SCT, integrated, rat_RNA, etc.
    dimreducs = NULL,         # drop PCA/UMAP
    graphs  = NULL
)

# Now convert — should be ~6x faster with 1 assay instead of 7
ref_slim[["mouse_RNA"]] <- as(ref_slim[["mouse_RNA"]], "Assay")
DefaultAssay(ref_slim) <- "mouse_RNA"
ref_slim@assays$mouse_RNA@data <- ref_slim@assays$mouse_RNA@counts
sceasy::convertFormat(ref_slim, from="seurat", to="anndata",
                      assay="mouse_RNA", outFile='DVC_mouse_ref_2026.h5ad')


# In[24]:


setwd("/maps/projects/perslab/people/qbx224/cellranger_ouputs_Aug2024/LPBN")


# In[37]:


SCOP_rename <- function(h5_cb, ind, sun) {
  newind <- paste0("_", ind)
  newobj <- substr(deparse(substitute(h5_cb)), 1, 2)
  colnames(h5_cb) <- gsub('-1',newind,colnames(h5_cb)) #SCOP pipeline changes dashes to underscores
  print(paste0(newobj, " appended _", ind)) #check object name and appended underscore

  keep.cells <- intersect(colnames(h5_cb), colnames(sun)) #intersect cellbender and SCOP file to only keep called cells
  h5_cb <- h5_cb[, keep.cells] #subset cellbender matrix
  print(paste0("Cells from cellranger.h5 present in SCOP object: ", (length(keep.cells))))

  h5_cb <- CreateSeuratObject(counts = h5_cb, min.cells = 10) #Seurat object

  lane <- paste0("SI-TT-", newobj) #get lane ID by pulling object name
  #print(paste0("Cells of corresponding lane present in SCOP object: ", (sum(sun$orig.ident %in% c(lane)))))#Check number of cells with corresponding lane metadata

  #lane <- paste0(lane, "_", deparse(substitute(sun))) #create unique lane identities
  h5_cb$old.ident <- lane #add to metadata
  hashes <- sun@meta.data #take SCOP metadata
  h5_cb <- AddMetaData(h5_cb, hashes) #add SCOP metadata

  DefaultAssay(h5_cb) <- "RNA"
  h5_cb <- NormalizeData(h5_cb)
  h5_cb <- FindVariableFeatures(h5_cb, selection.method = "vst", verbose = F)
  return(h5_cb)
}

Seurat_normalize_integrated <- function(seurat_obj, res = 0.8, dims = 1:50, ...) {
  #seurat_obj <- JoinLayers(seurat_obj)
  seurat_obj <- ScaleData(object = seurat_obj)
  seurat_obj <- RunPCA(object = seurat_obj)
  seurat_obj <- FindNeighbors(object = seurat_obj, dims = dims)
  seurat_obj <- FindClusters(object = seurat_obj, resolution = res)
  seurat_obj <- RunUMAP(object = seurat_obj, dims = dims)
  return(seurat_obj)
}


# In[36]:


raw <- readRDS("singlets.RDS")
raw


# In[38]:


E3 <- Read10X_h5("E3_filtered_feature_bc_matrix.h5")
E3 <- SCOP_rename(E3, "1", raw)


# In[39]:


E7 <- Read10X_h5("E7_filtered_feature_bc_matrix.h5")
E7 <- SCOP_rename(E7, "2", raw)


# In[40]:


E10 <- Read10X_h5("E10_filtered_feature_bc_matrix.h5")
E10 <- SCOP_rename(E10, "3", raw)


# In[41]:


E12 <- Read10X_h5("E12_filtered_feature_bc_matrix.h5")
E12 <- SCOP_rename(E12, "4", raw)


# In[42]:


variable_names <- c(E3, E7, E10, E12)
features <- SelectIntegrationFeatures(variable_names) #intregration features across all Seurat objects
variable_names <- lapply(X = variable_names, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
}) #apply these features for scaling and PCA on all Seurat objects separately
anchors <- FindIntegrationAnchors(object.list = variable_names, 
                                   anchor.features = features, 
                                   reduction = "rpca", 
                                   dims = 1:50) #create anchorset with all Seurat objects and their variable features


# In[43]:


HYP <- IntegrateData(anchorset = anchors, dims = 1:50) #integrate objects
HYP <- Seurat_normalize_integrated(HYP)
qs::qsave(HYP, file = "HYP_integrated_LPBN.qs")


# In[47]:


DefaultAssay(HYP) <- "RNA"
FeaturePlot(HYP, features = "v147")


# In[48]:


# Key regional markers
# Hypothalamus: Agrp, Npy, Pomc, Hcrt, Avp, Oxt, Sim1
# Hindbrain/DVC: Phox2b, Dbx1, Chat (NTS/DMV neurons), Th (area postrema)

markers <- c(
    "Agrp",  # Agrp  — ARC hypothalamus
    "Npy",  # Npy   — ARC hypothalamus  
    "Pomc",  # Pomc  — ARC hypothalamus
    "Hcrt",  # Hcrt  — LH hypothalamus
    "Avp",  # Avp   — PVH/SON hypothalamus
    "Phox2b",  # Phox2b — hindbrain/DVC
    "Chat",  # Chat  — hindbrain motor/DVC
    "Th"   # Th    — area postrema/LC
)

# Quick expression check
DotPlot(HYP, features = markers) + RotatedAxis()


# In[50]:


HYP <- JoinLayers(HYP)
# Sum of counts across all cells — a zero means the region doesn't express it
rowSums(HYP@assays$RNA$counts[markers, ])


# In[3]:


raw <- qs::qread("HYP_integrated_LPBN.qs")
DefaultAssay(raw) <- "RNA"
raw <- JoinLayers(raw)
raw


# In[4]:


ref <- qs::qread("/maps/projects/perslab/people/qbx224/7plex_paper/atlas/7plex_Seurat_bipca_annotated.qs")
ref


# In[5]:


ref <- SCTransform(ref, verbose = F)
## Transfer neuron and glia labels
raw <- SCTransform(raw, assay = "RNA", verbose = F)


# In[6]:


###Load and transfer
TransferAnchors <- FindTransferAnchors(reference = ref, 
                                       query = raw, 
                                       reduction = "cca",
                                       normalization.method = "SCT",
                                       recompute.residuals = F)
#Data to transfer
predictions <- TransferData(anchorset = TransferAnchors, weight.reduction = "cca",
                             refdata = ref$detailed_celltype_annotation)
raw$detailed_celltype_annotation <- predictions$predicted.id
raw$score_detailed_celltype_annotation <- predictions$prediction.score.max
####
qs::qsave(raw, "HYP_integrated_LPBN_labelled.qs")


# In[7]:


DimPlot(raw, 
        raster = F, 
        label = T, 
        repel = T, 
        group.by = "detailed_celltype_annotation", 
        cols = colors_detailed) + NoLegend()


# In[3]:


seurat_object <- qs::qread("HYP_integrated_LPBN_labelled.qs")
seurat_object


# In[5]:


table(seurat_object$sample_id, seurat_object$old.ident)


# In[6]:


library(sceasy)
library(anndata)
library(reticulate)
use_condaenv('pca_env')


# In[7]:


DefaultAssay(seurat_object) <- "RNA"
seurat_object[["SCT"]] <- NULL
seurat_object <- JoinLayers(seurat_object)


# In[8]:


# Convert Assay5 -> Assay3 (legacy)
seurat_object[["RNA"]] <- as(object = seurat_object[["RNA"]], Class = "Assay")

# Verify it's now a legacy assay
class(seurat_object[["RNA"]])  # Should print "Assay", not "Assay5"


# In[9]:


seurat_object@assays$RNA@data <- seurat_object@assays$RNA@counts # Sceasy uses the data layer which causes the normalized count matrix to be exported instead of the raw counts


# In[10]:


sceasy::convertFormat(seurat_object, 
                      from="seurat", 
                      to="anndata", 
                      assay="RNA", 
                      outFile='HYP_integrated_LPBN_labelled.h5ad')


# In[ ]:




