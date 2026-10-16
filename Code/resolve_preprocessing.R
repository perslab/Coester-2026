library(dplyr)
library(Seurat)
library(ggplot2)
'%!in%' <- function(x,y)!('%in%'(x,y))
options(future.globals.maxSize= 5*1024*1024^2)
###Preprocessing RESOLVE HYPOTHALAMUS
filepath <- "/resolve"
savepath <- "/processed_data"
files <- list.files(path=filepath, recursive=T, pattern = "*object_projection_call.qs")

for (f in files) {
  if (file.exists(paste0(savepath,"/",substr(f,0,8),"_",substr(f,10,11),"_transformed.qs"))){#check if already processed previously
    print(paste("Processed file already found for", f))
  }
  else{
    print(paste("Loading spatial object", f))
    spatial_object <- qs::qread(paste0(filepath,"/",f))
    spatial_object <- subset(spatial_object, subset = nCount_raw > 0)#Make sure there are no empty cells
    probes <- rownames(spatial_object)
    drop <- c("p147-mCherry-2A-iCre", "p36-NLS-Cre")
    miss <- probes %!in% drop
    probes <- probes[miss]
    #Subsetting the probes in the RNA assay
    #Pull the assay the original Seurat object containing only the genes of interest
    subset_matrix <- spatial_object[["raw"]]
    subset_matrix <- subset(subset_matrix, features = probes)
    #Add the virus-free assay to the Seurat object
    spatial_object[["RNA"]] <- subset_matrix
    #Standard processing
    DefaultAssay(spatial_object) <- "RNA"
    spatial_object$nFeature_RNA <- Matrix::colSums(spatial_object@assays$RNA$counts > 0)
    spatial_object$nCount_RNA <- Matrix::colSums(spatial_object@assays$RNA$counts)
    spatial_object <- subset(spatial_object, subset = nCount_RNA > 0)#Make sure there are no empty cells
    spatial_object <- SCTransform(spatial_object, assay = "RNA")
    spatial_object <- RunPCA(spatial_object, npcs = 30, features = probes)
    spatial_object <- RunUMAP(spatial_object, dims = 1:30)
    spatial_object <- FindNeighbors(spatial_object, reduction = "pca", dims = 1:30)
    spatial_object <- FindClusters(spatial_object, resolution = 0.3)
    DefaultAssay(spatial_object) <- "raw"
    ImageDimPlot(spatial_object, cols = "polychrome", border.color = NA, boundaries = "segmentation")+ ggtitle(substr(f,0,11))
    ggsave(paste0(savepath,"/",substr(f,0,8),"_",substr(f,10,11),"_DimPlot.pdf"))
    print("Saving preprocessed spatial object..")
    qs::qsave(spatial_object, (paste0(savepath,"/",substr(f,0,8),"_",substr(f,10,11),"_transformed.qs")))
  }
}
