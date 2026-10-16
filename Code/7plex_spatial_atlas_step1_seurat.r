library(dplyr)
library(Seurat)
library(ggplot2)
library(stats)
library(scran)
library(scater)
library(scran)
library(BiocParallel)
'%!in%' <- function(x,y)!('%in%'(x,y))
options(repr.plot.width = 12, repr.plot.height = 12)
spatialpath <- "/processed_data/"
atlaspath <- "/atlas/"
options(future.globals.maxSize= 50*1024*1024^2)

resolve_save <- function(plot, path){
    ggsave(path, 
           plot = plot, 
           device = cairo_pdf, 
           width = 12, 
           height = 12, 
           units = "in", 
           dpi = 300, 
           bg = "transparent"
          )
}

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

cells_manual <- c('00_Hmcn2_Rai14', '01_Astrocyte', '02_Sim1_Chrm2', '03_Oligodendrocyte',
       '04_Ttc6_Lmx1a', '05_Fezf1_Lama1', '06_Trhde_Ecel1', '07_Crabp1_Tbx3',
       '08_Qrfpr_Alk', '09_Tanycyte', '10_Lef1_Rspo2', '11_Slc17a8_Lgr5',
       '12_Tafa4_Kank4', '14_OPC', '15_Bcl11b_Cdh7', '16_Ghrh_Ebf1',
       '17_Npy_Agrp', '18_Pou1f1_Rxfp1', '19_Ntng2_Onecut1',
       '20_Tmem114_Npsr1', '21_Ebf3_Nxph4', '22_Dkk2_Lhx1', '23_Microglia',
       '24_Pth2r_Tacr1', '25_Pomc_Npy1r', '26_Lhx1_Th', '27_Pmfbp1_Onecut1',
       '28_Nts', '29_Lhx9_Qrfp', '30_Bace2_Ros1', '31_Hdc_Wif1', '32_VLMC',
       '33_Kiss1_Tac2', '34_Fibroblast', '35_Dnah11_Pth2r', '38_Sim1_Pla2r1',
       '39_Lhx9_Zic1', '40_Corticotropic_cell', '41_Endothelial_cell',
       '42_Gnrhr_Fshb', '43_Gpr50', '44_VLMC', '45_OPC', '47_Ependymal_cell',
       '49_Otx1_Igf1', '50_Schwann_cell', '51_Oligodendrocyte',
       '52_Fibroblast', '53_Astrocyte')

colors_manual <- c('#ffff00', '#1ce6ff', '#ff34ff', '#ff4a46', 
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
                   '#6f0062')

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

colors_region <-  c(
  "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
  "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF",
  "#997D87", "#5A0007", "#809693", "#FEFFE6", "#1B4400", "#3B5DFF", "#4A3B53", "#242424"
)

cells_region <- c('anterior_hypothalamus','arcuate_nucleus','bed_nucleus','dorsal_hypothalamus',
                  'dorsomedial_nucleus','lateral_amygdaloid_nucleus','lateral_hypothalamus',
                  'mammillary_nucleus','medial_preoptic_nucleus','midbrain',
                  'paraventricular_nucleus','posterior_hypothalamus','premammillary_nucleus',
                  'retrochiasmatic_area','retromammillary_area','suprachiasmatic_nucleus',
                  'supraoptic_nucleus','thalamus','tuberal_nucleus',
                  'ventromedial_nucleus','zona_incerta','undefined')

names(colors_manual) <- cells_manual
names(colors_detailed) <- cells_detailed
names(colors_region) <- cells_region

files <- list.files(path=spatialpath, recursive=F, pattern = "*bipca.qs")
files

###MERGE spatial objects
options(future.globals.maxSize= 50*1024*1024^2)
filepath <- spatialpath
S1A1 <- qs::qread(paste0(filepath, files[1]))
S1A1$old.ident <- "S1A1"
S1A2 <- qs::qread(paste0(filepath, files[2]))
S1A2$old.ident <- "S1A2"
S1B1 <- qs::qread(paste0(filepath, files[3]))
S1B1$old.ident <- "S1B1"
S1B2 <- qs::qread(paste0(filepath, files[4]))
S1B2$old.ident <- "S1B2"
S1C1 <- qs::qread(paste0(filepath, files[5]))
S1C1$old.ident <- "S1C1"
S1C2 <- qs::qread(paste0(filepath, files[6]))
S1C2$old.ident <- "S1C2"
S1D1 <- qs::qread(paste0(filepath, files[7]))
S1D1$old.ident <- "S1D1"
S1D2 <- qs::qread(paste0(filepath, files[8]))
S1D2$old.ident <- "S1D2"
S2A1 <- qs::qread(paste0(filepath, files[9]))
S2A1$old.ident <- "S2A1"
S2A2 <- qs::qread(paste0(filepath, files[10]))
S2A2$old.ident <- "S2A2"
S2B1 <- qs::qread(paste0(filepath, files[11]))
S2B1$old.ident <- "S2B1"
S2B2 <- qs::qread(paste0(filepath, files[12]))
S2B2$old.ident <- "S2B2"
S2C1 <- qs::qread(paste0(filepath, files[13]))
S2C1$old.ident <- "S2C1"
S2C2 <- qs::qread(paste0(filepath, files[14]))
S2C2$old.ident <- "S2C2"
S2D1 <- qs::qread(paste0(filepath, files[15]))
S2D1$old.ident <- "S2D1"
S2D2 <- qs::qread(paste0(filepath, files[16]))
S2D2$old.ident <- "S2D2"
S3A1 <- qs::qread(paste0(filepath, files[17]))
S3A1$old.ident <- "S3A1"
S3A2 <- qs::qread(paste0(filepath, files[18]))
S3A2$old.ident <- "S3A2"
S3B1 <- qs::qread(paste0(filepath, files[19]))
S3B1$old.ident <- "S3B1"
S3B2 <- qs::qread(paste0(filepath, files[20]))
S3B2$old.ident <- "S3B2"
S3C1 <- qs::qread(paste0(filepath, files[21]))
S3C1$old.ident <- "S3C1"
S3C2 <- qs::qread(paste0(filepath, files[22]))
S3C2$old.ident <- "S3C2"
S3D1 <- qs::qread(paste0(filepath, files[23]))
S3D1$old.ident <- "S3D1"
S3D2 <- qs::qread(paste0(filepath, files[24]))
S3D2$old.ident <- "S3D2"

object_names <- c("S1A1", "S1A2", "S1B1", "S1B2", 
                  "S1C1", "S1C2", "S1D1", "S1D2", 
                  "S2A1", "S2A2", "S2B1", "S2B2", 
                  "S2C1", "S2C2", "S2D1", "S2D2", 
                  "S3A1", "S3A2", "S3B1", "S3B2", 
                  "S3C1", "S3C2", "S3D1", "S3D2")
object_list <- list(S1A1, S1A2, S1B1, S1B2, S1C1, S1C2, S1D1, S1D2, 
                    S2A1, S2A2, S2B1, S2B2, S2C1, S2C2, S2D1, S2D2, 
                    S3A1, S3A2, S3B1, S3B2, S3C1, S3C2, S3D1, S3D2)

# Rename the FOVs
for (i in seq_along(object_list)) {
  current_object <- object_list[[i]]
  new_fov_name <- object_names[i]
  
  names(current_object@images) <- new_fov_name
  
  object_list[[i]] <- current_object
}
object_list[[1]]

# Reassign unique keys to each Seurat object
for (i in seq_along(object_list)) {
  object_list[[i]] <- RenameCells(object_list[[i]], 
                                 new.names = paste0(object_names[[i]], "_", colnames(object_list[[i]])))
}
head(colnames(object_list[[1]]))

# Merge all objects
seurat_merged <- object_list[[1]]
object_list <- object_list[-1]
seurat_merged <- merge(seurat_merged, y = object_list)

DefaultAssay(seurat_merged) <- "raw"
seurat_merged@assays$SCT <- NULL
qs::qsave(seurat_merged, paste0(atlaspath, "spatial_24_bipca_labels.qs"))

seurat_merged <- qs::qread(paste0(atlaspath, "spatial_24_bipca_labels.qs"))
seurat_merged

resolve_plot_dark <- function(spatial_merged, fov){
    plot <- ImageDimPlot(spatial_merged, border.color = NA, 
                         boundaries = "segmentation", 
                         group.by = "plot", 
                         cols = colors_detailed,
                         fov = fov
            ) + theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      legend.key.size = unit(0.5, "lines"),
                      legend.text = element_text(size = 6), 
                      plot.background = element_rect(fill = "black"),
                      panel.background = element_rect(fill = "black"), 
                      plot.margin = unit(c(0, 0, 0, 0), "cm") 
                     ) + scale_fill_manual(values = colors_detailed, 
                                           breaks = cells_detailed, name = "Detailed Cell Type Annotation")
}

resolve_plot_white <- function(spatial_merged, fov){
    plot <- ImageDimPlot(spatial_merged, 
                         border.color = "black",
                         border.size = 0.01,
                         boundaries = "segmentation", 
                         group.by = "plot", 
                         cols = colors_detailed,
                         fov = fov,
                         dark.background = FALSE
            ) + theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      legend.key.size = unit(0.5, "lines"),
                      legend.text = element_text(size = 9), 
                      plot.background = element_rect(fill = "white"),
                      panel.background = element_rect(fill = "white"), 
                      plot.margin = unit(c(0, 0, 0, 0), "cm") 
                     ) + scale_fill_manual(values = colors_detailed, 
                                           breaks = cells_detailed, name = "Detailed Cell Type Annotation")
}

ggplot(seurat_merged@meta.data, aes(x = score_detailed_celltype_annotation)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "white") +
  theme_minimal() +
  labs(
    title = "Distribution of Cell Type Annotation Scores",
    x = "score_detailed_celltype_annotation",
    y = "Number of cells"
  )

seurat_merged$plot <- seurat_merged$detailed_celltype_annotation
seurat_merged@meta.data[["plot"]] <- ifelse(seurat_merged@meta.data[["score_detailed_celltype_annotation"]] > 0.4, 
                                            seurat_merged@meta.data[["plot"]], "no_label")

plot <- resolve_plot_dark(seurat_merged, fov = "S3B2")
show(plot)

resolve_save(plot, paste0(spatialpath, "plots/testplot.pdf"))

fov_list <- c("S1A1", "S1A2", "S1B1", "S1B2", 
                  "S1C1", "S1C2", "S1D1", "S1D2", 
                  "S2A1", "S2A2", "S2B1", "S2B2", 
                  "S2C1", "S2C2", "S2D1", "S2D2", 
                  "S3A1", "S3A2", "S3B1", "S3B2", 
                  "S3C1", "S3C2", "S3D1", "S3D2")
for (f in fov_list) {
    plot <- resolve_plot_dark(seurat_merged, fov = f)
    resolve_save(plot, paste0(atlaspath, "plots/",f,"_detailed_celltype.pdf"))
}

fov_list <- c("S1A1", "S1A2", "S1B1", "S1B2", 
                  "S1C1", "S1C2", "S1D1", "S1D2", 
                  "S2A1", "S2A2", "S2B1", "S2B2", 
                  "S2C1", "S2C2", "S2D1", "S2D2", 
                  "S3A1", "S3A2", "S3B1", "S3B2", 
                  "S3C1", "S3C2", "S3D1", "S3D2")
for (f in fov_list) {
    plot <- resolve_plot_white(seurat_merged, fov = f)
    resolve_save(plot, paste0(atlaspath, "plots/",f,"_detailed_celltype_bright.pdf"))
}

metadata_sn <- read.csv(paste0(atlaspath,"allgenes_bipca_annotated_obs.csv"))
metadata_sn <- metadata_sn[,-1] #remove index column
#metadata_dict = dict(zip(meta["detailed_celltype_annotation"], meta["cellect_bmi_detailed"]))
#metadata_dict
metadata_sn

celltype_map <- setNames(metadata_sn$cellect_bmi_detailed,
                        metadata_sn$detailed_celltype_annotation)

metadata <- seurat_merged@meta.data
metadata$cellect_bmi_detailed <- celltype_map[as.character(metadata$detailed_celltype_annotation)]
seurat_merged <- AddMetaData(seurat_merged, 
                           metadata = metadata)

seurat_merged$cellect_bmi_log10 <- -log10(seurat_merged$cellect_bmi_detailed)

plot <- ImageFeaturePlot(seurat_merged, border.color = NA, 
                         boundaries = "segmentation", 
                         features = "cellect_bmi_log10",
                         fov = "S3B2") + scale_fill_distiller(palette = "Spectral")
show(plot)

resolve_plot_BMI <- function(spatial_merged, fov){
    plot <- ImageFeaturePlot(spatial_merged, 
                         border.color = NA,
                         boundaries = "segmentation", 
                         features = "cellect_bmi_log10",
                         fov = fov
            ) + theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(), 
                      plot.background = element_rect(fill = "black"),
                      panel.background = element_rect(fill = "black"), 
                      plot.margin = unit(c(0, 0, 0, 0), "cm") 
                     ) + scale_fill_distiller(palette = "Spectral")
}
fov_list <- c("S1A1", "S1A2", "S1B1", "S1B2", 
                  "S1C1", "S1C2", "S1D1", "S1D2", 
                  "S2A1", "S2A2", "S2B1", "S2B2", 
                  "S2C1", "S2C2", "S2D1", "S2D2", 
                  "S3A1", "S3A2", "S3B1", "S3B2", 
                  "S3C1", "S3C2", "S3D1", "S3D2")

for (f in fov_list) {
    plot <- resolve_plot_BMI(seurat_merged, fov = f)
    resolve_save(plot, paste0(atlaspath, "plots/",f,"_BMI.pdf"))
}

plot <- ImageDimPlot(seurat_merged,
                     border.color = NA,
                     boundaries = "segmentation",
                     group.by = "sun1", 
                     fov = "S3B2",
                     cols = "polychrome",
                     dark.background = TRUE
                    )
show(plot)

seurat_merged$plot <- seurat_merged$predicted_brain_region
seurat_merged@meta.data[["plot"]] <- ifelse(seurat_merged@meta.data[["score_predicted_brain_region"]] > 0.4, 
                                            seurat_merged@meta.data[["plot"]], NA)
ImageDimPlot(seurat_merged, 
             border.color = NA, 
             boundaries = "segmentation", 
             group.by = "plot", 
             fov = "S3B2",
             cols = "polychrome")

seurat_merged$virus1_log10 <- -log10(seurat_merged$virus1)
seurat_merged$virus2_log10 <- -log10(seurat_merged$virus2)

plot <- ImageFeaturePlot(seurat_merged, border.color = NA, 
                         boundaries = "segmentation", 
                         features = "virus1_log10",
                         max.cutoff = "q75",
                         fov = "S3B2") + scale_fill_distiller(palette = "Spectral")
show(plot)

qs::qsave(seurat_merged, paste0(atlaspath, "spatial_24_bipca_annotated.qs"))

spatial_object <- qs::qread(paste0(atlaspath, "spatial_24_bipca_annotated.qs"))
spatial_object

table(spatial_object$detailed_celltype_annotation, spatial_object$run_index)

unique(spatial_object$region) %>% sort()


ant <- Cells(spatial_object)[spatial_object@meta.data$region %in% c("anterior_hyopthalamus_left","anterior_hypothalamus","anterior_hypothalamus_left","anterior_hypothalamus_right")]
arc <- Cells(spatial_object)[spatial_object@meta.data$region %in% c("arcuate","arcuate_left","arcuate_right"  )]
bed <- Cells(spatial_object)[spatial_object@meta.data$region %in% c( "bed_nucleus","bed_nucleus_left","Bed_nucleus_left","bed_nucleus_right","Bed_nucleus_right")]
dor <- Cells(spatial_object)[spatial_object@meta.data$region %in% c("dorsal_hyothalamus","dorsal_hypothalamus")]
dmh <- Cells(spatial_object)[spatial_object@meta.data$region %in% c("dorsomedial","dorsomedial_hypothalamus","dorsomedial_left","dorsomedial_right")]
mid <- Cells(spatial_object)[spatial_object@meta.data$region %in% c("interfascicular","interpeduncular","interpreduncular","pontine_gray","pontine_grey","ventral_tegmental_area" )]
amy <- Cells(spatial_object)[spatial_object@meta.data$region %in% c("lateral_amygdaloid_left","lateral_amygdaloid_right")]
lat <- Cells(spatial_object)[spatial_object@meta.data$region %in% c("lateral_hypothalamus","lateral_hypothalamus_left","lateral_hypothalamus_right","lateral_hypothalamus_right-1")]
mam <- Cells(spatial_object)[spatial_object@meta.data$region %in% c("mammillary","mammillary_left","mammillary_right")]
mpn <- Cells(spatial_object)[spatial_object@meta.data$region %in% c( "medial_preoptic","medial_preoptic_left","medial_preoptic_right")]
pvn <- Cells(spatial_object)[spatial_object@meta.data$region %in% c("paraventricular","paraventricular_left","Paraventricular_left","paraventricular_right","Paraventricular_right")]
pos <- Cells(spatial_object)[spatial_object@meta.data$region %in% c("posterior_hypothalamus","posterior_hypothalamus_left","posterior_hypothalamus_right")]
prm <- Cells(spatial_object)[spatial_object@meta.data$region %in% c("premammillary","premammillary_dorsal_left","premammillary_dorsal_right","premammillary_left","premammillary_right","premammillary_ventral_left","premammillary_ventral_right","dorsal_tuberomammillary")]
rch <- Cells(spatial_object)[spatial_object@meta.data$region %in% c("retrochiasmaric_left","retrochiasmatic","retrochiasmatic_left","retrochiasmatic_right")]
rmm <- Cells(spatial_object)[spatial_object@meta.data$region %in% c("retromammillary","retromammillary_left","retromammillary_right" )]
sch <- Cells(spatial_object)[spatial_object@meta.data$region %in% c("suprachiamsatic","suprachiasmatic","suprachiasmatic_left","Suprachiasmatic_left","suprachiasmatic_right","Suprachiasmatic_right")]
son <- Cells(spatial_object)[spatial_object@meta.data$region %in% c("supraoptic_left","supraoptic_right")]
thm <- Cells(spatial_object)[spatial_object@meta.data$region %in% c("reuniens_thalamus","thalamus","thalamus_reuniens" )]
tub <- Cells(spatial_object)[spatial_object@meta.data$region %in% c("tuberal_left","tuberal_right")]
vmh <- Cells(spatial_object)[spatial_object@meta.data$region %in% c("ventromedial","ventromedial_left","ventromedial_right")]
zoi <- Cells(spatial_object)[spatial_object@meta.data$region %in% c( "zona_incerta","zona_incerta_left","zona_incerta_right" )]

spatial_object$region.collapsed <- case_when(rownames(spatial_object@meta.data) %in% ant ~ "anterior_hypothalamus",
                                             rownames(spatial_object@meta.data) %in% arc ~ "arcuate_nucleus",
                                             rownames(spatial_object@meta.data) %in% bed ~ "bed_nucleus",
                                             rownames(spatial_object@meta.data) %in% dor ~ "dorsal_hypothalamus",
                                             rownames(spatial_object@meta.data) %in% dmh ~ "dorsomedial_nucleus",
                                             rownames(spatial_object@meta.data) %in% mid ~ "midbrain",
                                             rownames(spatial_object@meta.data) %in% amy ~ "lateral_amygdaloid_nucleus",
                                             rownames(spatial_object@meta.data) %in% lat ~ "lateral_hypothalamus",
                                             rownames(spatial_object@meta.data) %in% mam ~ "mammillary_nucleus",
                                             rownames(spatial_object@meta.data) %in% mpn ~ "medial_preoptic_nucleus",
                                             rownames(spatial_object@meta.data) %in% pvn ~ "paraventricular_nucleus",
                                             rownames(spatial_object@meta.data) %in% pos ~ "posterior_hypothalamus",
                                             rownames(spatial_object@meta.data) %in% prm ~ "premammillary_nucleus",
                                             rownames(spatial_object@meta.data) %in% rch ~ "retrochiasmatic_area",
                                             rownames(spatial_object@meta.data) %in% rmm ~ "retromammillary_area",
                                             rownames(spatial_object@meta.data) %in% sch ~ "suprachiasmatic_nucleus",
                                             rownames(spatial_object@meta.data) %in% son ~ "supraoptic_nucleus",
                                             rownames(spatial_object@meta.data) %in% thm ~ "thalamus",
                                             rownames(spatial_object@meta.data) %in% tub ~ "tuberal_nucleus",
                                             rownames(spatial_object@meta.data) %in% vmh ~ "ventromedial_nucleus",
                                             rownames(spatial_object@meta.data) %in% zoi ~ "zona_incerta",
                                             T  ~ "undefined")

unique(spatial_object$region.collapsed) %>% sort()

qs::qsave(spatial_object, paste0(atlaspath, "spatial_24_bipca_annotated.qs"))

seurat_merged <- qs::qread(paste0(atlaspath, "spatial_24_bipca_annotated.qs"))
seurat_merged

seurat_merged$virus1_pos <- case_when(seurat_merged$virus1 <= 0.05 & seurat_merged$run_index == "spatial1" ~ "BST",
                                      seurat_merged$virus1 <= 0.05 & seurat_merged$run_index == "spatial2" ~ "LHA",
                                      seurat_merged$virus1 <= 0.05 & seurat_merged$run_index == "spatial3" ~ "PVT",
                                      T  ~ "no_projection")
seurat_merged$virus2_pos <- case_when(seurat_merged$virus2 <= 0.05 & seurat_merged$run_index == "spatial1" ~ "PAG",
                                      seurat_merged$virus2 <= 0.05 & seurat_merged$run_index == "spatial2" ~ "PBN",
                                      seurat_merged$virus2 <= 0.05 & seurat_merged$run_index == "spatial3" ~ "CeA",
                                      T  ~ "no_projection")

seurat_merged$projection <- case_when(seurat_merged$virus1_pos == "BST" & seurat_merged$virus2_pos == "PAG" ~ "BST+PAG",
                                      seurat_merged$virus1_pos == "LHA" & seurat_merged$virus2_pos == "PBN" ~ "LHA+PBN",
                                      seurat_merged$virus1_pos == "PVT" & seurat_merged$virus2_pos == "CeA" ~ "PVT+CeA",
                                      seurat_merged$virus1_pos == "BST" & seurat_merged$virus2_pos == "no_projection" ~ "BST",
                                      seurat_merged$virus1_pos == "LHA" & seurat_merged$virus2_pos == "no_projection" ~ "LHA",
                                      seurat_merged$virus1_pos == "PVT" & seurat_merged$virus2_pos == "no_projection" ~ "PVT",
                                      seurat_merged$virus1_pos == "no_projection" & seurat_merged$virus2_pos == "PAG" ~ "PAG",
                                      seurat_merged$virus1_pos == "no_projection" & seurat_merged$virus2_pos == "PBN" ~ "PBN",
                                      seurat_merged$virus1_pos == "no_projection" & seurat_merged$virus2_pos == "CeA" ~ "CeA",
                                      T  ~ "negative")

table(seurat_merged$virus1_pos, seurat_merged$run_index)
table(seurat_merged$virus2_pos, seurat_merged$run_index)
table(seurat_merged$projection, seurat_merged$run_index)

plot <- ImageDimPlot(seurat_merged, 
             border.color = NA, 
             boundaries = "segmentation", 
             group.by = "projection",
             molecules = c("p147_mCherry_2A_iCre", "p36_NLS_Cre"),
             mols.cols = c("cyan", "green"),
             mols.size = 0.01,
             fov = "S3B2",
             cols = "polychrome")+ theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      legend.key.size = unit(0.5, "lines"),
                      legend.text = element_text(size = 9), 
                      plot.background = element_rect(fill = "black"),
                      panel.background = element_rect(fill = "black"), 
                      plot.margin = unit(c(0, 0, 0, 0), "cm"))
show(plot)
resolve_save(plot, paste0(atlaspath, "plots/testplot.pdf"))

seurat_merged <- JoinLayers(seurat_merged)
seurat_merged

#qs::qsave(seurat_merged, paste0(atlaspath, "spatial_24_bipca_annotated.qs"))
seurat_merged <- qs::qread(paste0(atlaspath, "spatial_24_bipca_annotated.qs"))
seurat_merged

seurat_merged$comp1 <- seurat_merged$detailed_celltype_annotation
seurat_merged@meta.data[["comp1"]] <- ifelse(seurat_merged@meta.data[["score_detailed_celltype_annotation"]] > 0.4, 
                                            seurat_merged@meta.data[["comp1"]], NA)
table(seurat_merged$comp1)

metadata_sn <- read.csv(paste0(atlaspath,"allgenes_bipca_annotated_obs.csv"))
metadata_sn <- metadata_sn[,-1] #remove index column
celltype_map <- setNames(metadata_sn$AAV_ratio,
                        metadata_sn$detailed_celltype_annotation)
metadata <- seurat_merged@meta.data
metadata$AAV_ratio <- celltype_map[as.character(metadata$detailed_celltype_annotation)]
seurat_merged <- AddMetaData(seurat_merged, 
                           metadata = metadata)

plot <- ImageFeaturePlot(seurat_merged, border.color = NA, 
                         boundaries = "segmentation",
                         split.by = "area",
                         features = "AAV_ratio",
                         fov = "S3B2") 
#resolve_save(plot, paste0(atlaspath, "plots/testplot.pdf"))
show(plot)

library(sceasy)
library(anndata)
library(reticulate)
use_condaenv('pca_env')

###downgrade assays from assay5
seurat_object <- seurat_merged
seurat_object[["RNA"]] <- NULL
seurat_object@assays$raw$data <- seurat_object@assays$raw$counts
seurat_object

seurat_object@assays$raw

seurat_object@assays$raw$data <- seurat_object@assays$raw$counts
seurat_object <- NormalizeData(seurat_object)
seurat_object <- FindVariableFeatures(seurat_object)
seurat_object <- ScaleData(seurat_object)
seurat_object

seurat_object[["raw"]] <- as(seurat_object[["raw"]], "Assay")
seurat_object@assays$raw$scale.data <- NULL
seurat_object@assays$raw$data <- NULL
seurat_object

sceasy::convertFormat(seurat_object, from="seurat", to="anndata", assay="raw", outFile=paste0(atlaspath,'spatial_24_bipca_annotated.h5ad'))

hist(seurat_merged$nCount_raw, breaks = 1000, xlim=c(0,500))

g <- ggplot(na.omit(seurat_merged[[]][, c("region.collapsed", "nCount_raw")]), 
            aes(x = region.collapsed, y = nCount_raw, fill = region.collapsed))
g <- g + geom_violin() + xlab("region.collapsed") + ylab("nCount")
g <- g + scale_fill_manual(values = colors_region, breaks = cells_region)
g


ImageDimPlot(seurat_merged, 
             border.color = NA, 
             boundaries = "segmentation", 
             group.by = "region.collapsed",
             fov = "S3B2")+ scale_fill_manual(values = colors_region, 
                                           breaks = cells_region, name = "Region Annotation")

plot <- ImageFeaturePlot(seurat_merged, border.color = NA, 
                         boundaries = "segmentation", 
                         features = "nCount_raw",
                         max.cutoff = "q75",
                         fov = "S3B2") + scale_fill_distiller(palette = "Spectral")
show(plot)

seurat_merged@images$S1A1

glia <- c("10_Npy_Agrp",
"13_Pth2r_Otp",
"14_Lhx1_Th",
"25_Nxph4_Ebf3",
"32_Gucy2c_Ttc6",
"40_Gpr50_Tac2")

seurat_merged$plot <- seurat_merged$detailed_celltype_annotation
seurat_merged@meta.data[["plot"]] <- ifelse(seurat_merged@meta.data[["plot"]] %in% glia, seurat_merged@meta.data[["plot"]], "no_label")
seurat_merged@meta.data[["plot"]] <- ifelse(seurat_merged@meta.data[["score_detailed_celltype_annotation"]] > 0.4, seurat_merged@meta.data[["plot"]], "no_label")

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
                      '#456648', '#0086ed', "#121212")
names(colors_detailed) <- cells_detailed

ImageDimPlot(seurat_merged, 
             border.color = NA, 
             boundaries = "segmentation", 
             group.by = "plot", 
             cols = colors_detailed,
             fov = "S3B2",
             dark.background = TRUE,
            )
ggsave(paste0(atlaspath,"testplot.png"))

glia <- c("01_Undefined", "31_Lhx9_Npvf")
seurat_merged$plot <- seurat_merged$detailed_celltype_annotation
seurat_merged@meta.data[["plot"]] <- ifelse(seurat_merged@meta.data[["plot"]] %in% glia, seurat_merged@meta.data[["plot"]], "no_label")
seurat_merged@meta.data[["plot"]] <- ifelse(seurat_merged@meta.data[["score_detailed_celltype_annotation"]] > 0.4, seurat_merged@meta.data[["plot"]], "no_label")
ImageDimPlot(seurat_merged, 
             border.color = NA, 
             boundaries = "segmentation", 
             group.by = "plot", 
             cols = colors_detailed,
             fov = "S3D1",
             dark.background = TRUE,
            )
#ggsave(paste0(atlaspath,"testplot.png"))

df <- table(seurat_merged$projection, seurat_merged$comp1)
write.csv(df, "projection_calls_spatial.csv", row.names=TRUE)

df <- table(seurat_merged$comp1)
write.csv(df, "cells_calls_spatial.csv", row.names=TRUE)

sn_object <- qs::qread(paste0(atlaspath, "7plex_bipca_n98_region_labeltransfer.qs"))
sn_object

DimPlot(sn_object, 
        raster = F, 
        label = T, 
        repel = F, 
        group.by = "detailed_celltype_annotation", 
        reduction = "tsne", 
        cols = colors_detailed) + NoLegend()

DimPlot(sn_object, 
        raster = F, 
        group.by = "region_collapsed", 
        reduction = "tsne",
        cols = godsnot_102
       )+NoLegend()
ggsave("region_collapsed_tsne.pdf")

regions <- c("anterior_hypothalamus","arcuate_nucleus","bed_nucleus",
             "dorsal_hypothalamus","dorsomedial_nucleus",
             "lateral_amygdaloid_nucleus","lateral_hypothalamus",
             "mammillary_nucleus","medial_preoptic_nucleus","midbrain",
             "paraventricular_nucleus","posterior_hypothalamus",
             "premammillary_nucleus","retrochiasmatic_area",
             "retromammillary_area","suprachiasmatic_nucleus",
             "supraoptic_nucleus","thalamus","tuberal_nucleus",
             "undefined","ventromedial_nucleus","zona_incerta", "no_label")

colors_regions <- c(
  "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
  "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF",
  "#997D87", "#5A0007", "#809693", "#FF2F80", "#1B4400", "#4FC601", "#3B5DFF",
  "#4A3B53", "#FEFFE6")
names(colors_regions) <- regions

length(regions)

sn_object$plot <- sn_object$region_collapsed
sn_object@meta.data[["plot"]] <- ifelse(sn_object@meta.data[["score_region_collapsed"]] > 0.4, 
                                            sn_object@meta.data[["plot"]], "no_label")
DimPlot(sn_object, 
        raster = F, 
        group.by = "plot", 
        reduction = "tsne",
        cols = colors_regions
       ) + NoLegend()
#ggsave("region_collapsed_tsne_filtered.pdf")

sn_subset <- subset(sn_object, cells = Cells(sn_object)[sn_object$plot %in% c("no_label")], invert = T)

DimPlot(sn_subset, 
        raster = F, 
        group.by = "plot", 
        reduction = "tsne",
        cols = colors_regions
       ) + NoLegend()
ggsave("region_collapsed_tsne_filtered.pdf")

DimPlot(sn_subset, 
        raster = F, 
        group.by = "plot", 
        reduction = "tsne",
        cols = colors_regions
       )
ggsave("region_collapsed_tsne_filtered_legend.pdf")

df <- table(sn_subset$detailed_celltype_annotation,sn_subset$region_collapsed)
write.csv(df, "regions_filtered.csv", row.names=TRUE)

df <- table(sn_object$detailed_celltype_annotation,sn_object$region_collapsed)
write.csv(df, "regions.csv", row.names=TRUE)


write.csv(sn_object[[]], "7plex_region_labeltransfer_meta.csv", row.names=TRUE)

spatial <- readRDS("spatial_object_preprocessed+labelled.rds")
spatial

resolve_plot_dark <- function(spatial_merged, fov, cols){
    plot <- ImageDimPlot(spatial_merged, border.color = NA, 
                         boundaries = "segmentation", 
                         group.by = "plot", 
                         cols = cols,
                         fov = fov
            ) + theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      legend.key.size = unit(0.5, "lines"),
                      legend.text = element_text(size = 6), 
                      plot.background = element_rect(fill = "black"),
                      panel.background = element_rect(fill = "black"), 
                      plot.margin = unit(c(0, 0, 0, 0), "cm") 
                     )
}


spatial$plot <- spatial$predicted.cell.type
spatial@meta.data[["plot"]] <- ifelse(spatial@meta.data[["prediction.score.cell.type"]] > 0.6, 
                                            spatial@meta.data[["plot"]], "no_label")

plot <- resolve_plot_dark(spatial, fov = "B1.2", cols = godsnot_102)
show(plot)

glia <- c("Sero_0_0", "Glu_4_2")
spatial@meta.data[["plot"]] <- ifelse(spatial@meta.data[["plot"]] %in% glia, spatial@meta.data[["plot"]], "no_label")
plot <- resolve_plot_dark(spatial, fov = "B1.2", cols = "polychrome")
show(plot)
#ggsave(paste0(atlaspath,"testplot.pdf"))

plot <- ImageDimPlot(spatial, 
             border.color = NA, 
             boundaries = "segmentation", 
             group.by = "plot",
             molecules = c("Adcyap1", "Glp1r"),
             mols.cols = c("white", "cyan"),
             mols.size = 0.01,
             mols.alpha = 0.1,
             nmols=1000000,
             fov = "B2.2",
             cols = godsnot_102)+ theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      legend.key.size = unit(0.5, "lines"),
                      legend.text = element_text(size = 9), 
                      plot.background = element_rect(fill = "black"),
                      panel.background = element_rect(fill = "black"), 
                      plot.margin = unit(c(0, 0, 0, 0), "cm")) + NoLegend()
show(plot)
ggsave(paste0(atlaspath,"adcyap_nolegend.pdf"))

DefaultAssay(spatial) <- "raw"
spatial <- JoinLayers(spatial)
# Access raw counts
counts <- GetAssayData(spatial, layer = "counts")
# Find cells with Tph2 count > 1
cells_to_keep <- colnames(spatial)[counts["Adcyap1", ] > 0]
# Subset the Seurat object
spatial_subset <- subset(spatial, cells = cells_to_keep)
table(spatial_subset$plot, spatial_subset$region)

table(spatial$plot, spatial$region)

spatial


