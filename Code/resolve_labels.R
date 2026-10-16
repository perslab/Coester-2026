library(dplyr)
library(Seurat)
library(ggplot2)
'%!in%' <- function(x,y)!('%in%'(x,y))
options(future.globals.maxSize= 5*1024*1024^2)
###Load and subset snSEQ data MOUSE HYPOTHALAMUS
filepath <- "/processed_data"
sn_object <- qs::qread("7plex_reference_n98_resolve.qs")
files <- list.files(path=filepath, recursive=F, pattern = "*transformed.qs")
for (f in files) {
  if (file.exists(paste0(filepath,"/",substr(f,0,11),"_annotated.qs"))){#check if already processed previously
    print(paste("Processed file already found for", f))
  }
  else{
    print(paste("Loading spatial object", f))
    spatial_object <- qs::qread(paste0(filepath,"/",f))
    DefaultAssay(spatial_object) <- "SCT"
    TransferAnchors <- FindTransferAnchors(reference = sn_object, query = spatial_object, reduction = "cca",normalization.method = "SCT",recompute.residuals = F)
    #Data to transfer. Multiple rows
    predictions1 <- TransferData(anchorset = TransferAnchors, weight.reduction = "cca", refdata = sn_object$minor.cell.type)
    predictions2 <- TransferData(anchorset = TransferAnchors, weight.reduction = "cca", refdata = sn_object$original.cell.type)
    predictions3 <- TransferData(anchorset = TransferAnchors, weight.reduction = "cca", refdata = sn_object$anatomical.area)
    spatial_object$predicted.minor.cell.type <- predictions1$predicted.id
    spatial_object$predicted.original.cell.type <- predictions2$predicted.id
    spatial_object$predicted.anatomical.area <- predictions3$predicted.id
    spatial_object$prediction.score.minor.cell.type <- predictions1$prediction.score.max
    spatial_object$prediction.score.original.cell.type <- predictions2$prediction.score.max
    spatial_object$prediction.score.anatomical.area <- predictions3$prediction.score.max
    #Save labelled spatial objects
    print("Saving label transfer spatial object..")
    qs::qsave(spatial_object, (paste0(filepath,"/",substr(f,0,11),"_annotated.qs")))
  }
}


