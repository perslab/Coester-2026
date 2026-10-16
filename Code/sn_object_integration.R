library(dplyr)
library(Seurat)
options(future.globals.maxSize= 50*1024*1024^2)

Seurat_large <- function(seurat_obj, split = "orig.ident", sct = F, res = 0.8, dims = 1:50, reference = c(1, 2, 7, 10, 14, 18, 19), ...) {
  DefaultAssay(seurat_obj) <- "RNA"
  #seurat_obj[['SCT']] <- NULL
  seurat.list <- SplitObject(seurat_obj, split.by = split)
  seurat.list <- lapply(X = seurat.list, FUN = function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, selection.method = "vst", verbose = FALSE)
  })
  features <- SelectIntegrationFeatures(seurat.list)
  seurat.list <- lapply(X = seurat.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
  })
  anchors <- FindIntegrationAnchors(object.list = seurat.list, reference = reference, anchor.features = features, reduction = "rpca", dims = dims)
  seurat.integrated <- IntegrateData(anchorset = anchors, dims = dims, )
  seurat.integrated <- ScaleData(seurat.integrated, verbose = FALSE)
  seurat.integrated <- RunPCA(seurat.integrated, verbose = FALSE)
  seurat.integrated <- FindNeighbors(seurat.integrated, dims = dims)
  seurat.integrated <- FindClusters(seurat.integrated, resolution = res)
  seurat.integrated <- RunUMAP(seurat.integrated, dims = dims)
  DefaultAssay(seurat.integrated) <- "RNA"
  seurat.integrated <- JoinLayers(seurat.integrated)
  return(seurat.integrated)
}
filepath <- "/7plex_paper/"
sn_object <- qs::qread(paste0(filepath, "7plex_reference_combined.qs"))
sn_object <- Seurat_large(sn_object, split = "old.ident")
qs::qsave(sn_object, paste0(filepath, "7plex_reference_combined_integrated.qs"))
