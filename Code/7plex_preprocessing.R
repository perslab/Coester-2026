library(Seurat)
library(data.table)
filepath <- ""
all_genes <- data.table::fread(paste0(filepath, "gene_names.txt"), sep = "\t")
all_genes <- as.matrix(all_genes) #32141
filtered_genes <- data.table::fread(paste0(filepath, "filtered_gene_names.txt"), sep = "\t")
filtered_genes <- as.matrix(filtered_genes) #20833
neurons <- qs::qread(file = paste0(filepath,"7plex_neurons_labelled_annotated.qs"))
print(length(rownames(neurons))) #23957
match <- intersect(rownames(neurons), all_genes) #23946
print(length(match))
filtered_match <- intersect(rownames(neurons), filtered_genes) #16934
print(length(filtered_match))
subset_matrix <- neurons[["RNA"]]
subset_matrix <- subset(subset_matrix, features = filtered_match)
#Add the filtered assay to the Seurat object
neurons[["RNA_filtered"]] <- subset_matrix
#Standard processing
DefaultAssay(neurons) <- "RNA"
rm(subset_matrix)
dims = 1:50
res = 0.8
DefaultAssay(neurons) <- "RNA_filtered"

seurat.list <- SplitObject(neurons, split.by = "old.ident")
seurat.list <- lapply(X = seurat.list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, selection.method = "vst", verbose = FALSE, assay = "RNA_filtered")
})
features <- SelectIntegrationFeatures(seurat.list, assay = c("RNA_filtered","RNA_filtered","RNA_filtered","RNA_filtered","RNA_filtered","RNA_filtered","RNA_filtered","RNA_filtered","RNA_filtered","RNA_filtered","RNA_filtered","RNA_filtered","RNA_filtered","RNA_filtered","RNA_filtered","RNA_filtered","RNA_filtered","RNA_filtered","RNA_filtered","RNA_filtered","RNA_filtered"))
seurat.list <- lapply(X = seurat.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})
anchors <- FindIntegrationAnchors(object.list = seurat.list, reference = c(1, 2, 7, 10, 14, 18, 19), anchor.features = features, reduction = "rpca", dims = dims)
seurat.integrated <- IntegrateData(anchorset = anchors, dims = dims)
seurat.integrated <- ScaleData(seurat.integrated, verbose = FALSE)
seurat.integrated <- RunPCA(seurat.integrated, verbose = FALSE)
seurat.integrated <- FindNeighbors(seurat.integrated, dims = dims)
seurat.integrated <- FindClusters(seurat.integrated, resolution = res)
seurat.integrated <- RunUMAP(seurat.integrated, dims = dims)
DefaultAssay(seurat.integrated) <- "RNA"
seurat.integrated <- JoinLayers(seurat.integrated)
qs::qsave(seurat.integrated, paste0(filepath, "7plex_neurons_reprocessed.qs"))
