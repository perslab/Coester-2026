library(dplyr)
library(Seurat)
library(ggplot2)
'%!in%' <- function(x,y)!('%in%'(x,y))
options(future.globals.maxSize= 50*1024*1024^2)
atlaspath <- "/7plex_paper/atlas/"
###Load and transfer spatial data MOUSE HYPOTHALAMUS
spatial_object <- qs::qread(paste0(atlaspath, "spatial_24_bipca_annotated.qs"))
spatial_object <- JoinLayers(spatial_object)
spatial_object <- SCTransform(spatial_object, assay = "RNA")
sn_object <- qs::qread(paste0(atlaspath,"7plex_bipca_n98.qs"))
DefaultAssay(sn_object) <- "SCT"

TransferAnchors <- FindTransferAnchors(reference = spatial_object, query = sn_object, reduction = "cca",normalization.method = "SCT",recompute.residuals = F)
#Data to transfer
predictions1 <- TransferData(anchorset = TransferAnchors, weight.reduction = "cca", refdata = spatial_object$region.collapsed)
sn_object$region_collapsed <- predictions1$predicted.id
sn_object$score_region_collapsed <- predictions1$prediction.score.max

qs::qsave(sn_object, (paste0(atlaspath,"7plex_bipca_n98_region_labeltransfer.qs")))
