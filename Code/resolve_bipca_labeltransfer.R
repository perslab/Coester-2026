library(dplyr)
library(Seurat)
library(ggplot2)
'%!in%' <- function(x,y)!('%in%'(x,y))
options(future.globals.maxSize= 5*1024*1024^2)
###Load and subset snSEQ data MOUSE HYPOTHALAMUS
filepath <- "/processed_data"
sn_object <- qs::qread("/atlas/7plex_bipca_n98.qs")
files <- list.files(path=filepath, recursive=F, pattern = "*transformed.qs")
for (f in files) {
  if (file.exists(paste0(filepath,"/",substr(f,0,11),"_bipca.qs"))){#check if already processed previously
    print(paste("Processed file already found for", f))
  }
  else{
    print(paste("Loading spatial object", f))
    spatial_object <- qs::qread(paste0(filepath,"/",f))
    DefaultAssay(spatial_object) <- "SCT"
    TransferAnchors <- FindTransferAnchors(reference = sn_object, query = spatial_object, reduction = "cca",normalization.method = "SCT",recompute.residuals = F)
    #Data to transfer. Multiple rows
    predictions1 <- TransferData(anchorset = TransferAnchors, weight.reduction = "cca", refdata = sn_object$detailed_celltype_annotation)
    predictions2 <- TransferData(anchorset = TransferAnchors, weight.reduction = "cca", refdata = sn_object$manual_celltype_annotation)
    predictions3 <- TransferData(anchorset = TransferAnchors, weight.reduction = "cca", refdata = sn_object$predicted_brain_region)
    spatial_object$detailed_celltype_annotation <- predictions1$predicted.id
    spatial_object$manual_celltype_annotation <- predictions2$predicted.id
    spatial_object$predicted_brain_region <- predictions3$predicted.id
    spatial_object$score_detailed_celltype_annotation <- predictions1$prediction.score.max
    spatial_object$score_manual_celltype_annotation <- predictions2$prediction.score.max
    spatial_object$score_predicted_brain_region <- predictions3$prediction.score.max
    #Save labelled spatial objects
    print(paste("Saving label transfer spatial object", f))
    qs::qsave(spatial_object, (paste0(filepath,"/",substr(f,0,11),"_bipca.qs")))
  }
}


