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
spatialpath <- "/projects/perslab/people/qbx224/7plex_paper/processed_data/"
atlaspath <- "/projects/perslab/people/qbx224/7plex_paper/atlas/"
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

seurat_merged <- qs::qread(paste0(atlaspath, "spatial_24_bipca_annotated.qs"))
seurat_merged

seurat_merged$region[is.na(seurat_merged$region)] <- "no_annotation"

id6 <- Cells(seurat_merged)[seurat_merged@meta.data$old.ident %in% c("S1A1", "S1A2", "S1B1", "S1B2", "S1C1", "S1C2")]
id7 <- Cells(seurat_merged)[seurat_merged@meta.data$old.ident %in% c("S1D1", "S1D2")]
id43 <- Cells(seurat_merged)[seurat_merged@meta.data$old.ident %in% c("S2A1", "S2A2", "S2B1", "S2B2", "S2C1", "S2C2")]
id39 <- Cells(seurat_merged)[seurat_merged@meta.data$old.ident %in% c("S2D1", "S2D2")]
id45 <- Cells(seurat_merged)[seurat_merged@meta.data$old.ident %in% c("S3A1", "S3A2", "S3B1", "S3B2", "S3C1", "S3C2")]
id41 <- Cells(seurat_merged)[seurat_merged@meta.data$old.ident %in% c("S3D1", "S3D2")]


seurat_merged$animal_id <- case_when(rownames(seurat_merged@meta.data) %in% id6 ~ "6",
                                             rownames(seurat_merged@meta.data) %in% id7 ~ "7",
                                             rownames(seurat_merged@meta.data) %in% id43 ~ "43",
                                             rownames(seurat_merged@meta.data) %in% id39 ~ "39",
                                             rownames(seurat_merged@meta.data) %in% id45 ~ "45",
                                             rownames(seurat_merged@meta.data) %in% id41 ~ "41",
                                             T  ~ NA)
table(seurat_merged$animal_id, useNA = "always")

seurat_merged$injection_side <- case_when(rownames(seurat_merged@meta.data) %in% id6 ~ "left",
                                             rownames(seurat_merged@meta.data) %in% id7 ~ "right",
                                             rownames(seurat_merged@meta.data) %in% id43 ~ "right",
                                             rownames(seurat_merged@meta.data) %in% id39 ~ "left",
                                             rownames(seurat_merged@meta.data) %in% id45 ~ "right",
                                             rownames(seurat_merged@meta.data) %in% id41 ~ "left",
                                             T  ~ NA)
table(seurat_merged$injection_side, useNA = "always")

seurat_merged$plot <- seurat_merged$detailed_celltype_annotation
seurat_merged@meta.data[["plot"]] <- ifelse(seurat_merged@meta.data[["score_detailed_celltype_annotation"]] > 0.4, 
                                            seurat_merged@meta.data[["plot"]], "no_label")

plot <- resolve_plot_dark(seurat_merged, fov = "S3B2")
show(plot)

ImageDimPlot(seurat_merged, 
             border.color = NA, 
             boundaries = "segmentation", 
             group.by = "region", 
             cols = godsnot_102,
             fov = "S3B2",
             dark.background = TRUE,
            )

#the "arcuate" annotation is in S3A2 in a very anterior slice of the ARC
seurat_merged$region <- tolower(seurat_merged$region)
corrections <- c(
  "anterior_hyopthalamus_left" = "anterior_hypothalamus_left",
  "retrochiasmaric_left" = "retrochiasmatic_left",
  "lateral_hypothalamus_right-1" = "lateral_hypothalamus_left",
  "dorsal_hyothalamus" = "dorsal_hypothalamus",
  "thalamus_reuniens" = "reuniens_thalamus",
  "arcuate" = "anterior_arcuate")
seurat_merged$region <- ifelse(
  seurat_merged$region %in% names(corrections),
  corrections[seurat_merged$region],
  seurat_merged$region
)
unique(seurat_merged$region) %>% sort()

ImageDimPlot(seurat_merged, 
             border.color = NA, 
             boundaries = "segmentation", 
             group.by = "region", 
             cols = godsnot_102,
             fov = "S3B2",
             dark.background = TRUE,
            )

#tabulating stats for comparison with snRNAseq
seurat_merged$comp1 <- seurat_merged$detailed_celltype_annotation
seurat_merged@meta.data[["comp1"]] <- ifelse(seurat_merged@meta.data[["score_detailed_celltype_annotation"]] > 0.4, 
                                            seurat_merged@meta.data[["comp1"]], NA)
#projection calls
df <- table(seurat_merged$projection, seurat_merged$comp1)
write.csv(df, "projection_calls_spatial.csv", row.names=TRUE)
#total counts
df <- table(seurat_merged$comp1)
write.csv(df, "cells_calls_spatial.csv", row.names=TRUE)

#############NEUROTRANSMITTER ANNOTATION
glut <- c(
  "00_Hmcn2_Gpr149", "02_Fezf1_Lama1", "03_Rai14_Acvr1c",
  "04_Crnde_Lmx1a", "05_Qrfpr_Pdzrn3", "09_Nfia_Ccbe1",
  "11_Tafa4_Kank4", "12_Ntng2_Onecut1", "17_Dkk2_Ebf3",
  "18_Ust_Bcl11b", "19_Bace2_Tspan18", "23_Ust_Nkd2",
  "25_Nxph4_Ebf3", "27_Lmx1a_Ttc6", "28_Kiss1_Slc38a3",
  "30_Trh_Sim1", "31_Lhx9_Npvf", "32_Gucy2c_Ttc6",
  "34_Pla2r1_Sim1", "35_Lhx9_Zic1", "36_Ebf2_Crhbp",
  "37_Qrfp_Qrfprl", "41_Otx1_Zic1", "42_Foxb1_Myo5b"
)
gaba_glut <- c(
  "01_Undefined"
)
gaba <- c(
  "06_Tbx3_Npy1r", "07_Nfib_Six3", "08_Lef1_Rspo2",
  "10_Npy_Agrp", "13_Pth2r_Otp", "14_Lhx1_Th",
  "15_Onecut1_Slc32a1", "16_Gal_Sox6", "20_Tmem114_Npsr1",
  "21_Bnc2_Glp1r", "22_Crabp1_Vgll3", "24_Nts_Ecel1",
  "26_Dach2_Npsr1", "29_Dnah11_Pth2r", "39_Nfib_Nfia",
  "40_Gpr50_Tac2"
)
chat <- c(
  "33_Chat_Ghrh"
)
hist <- c(
  "38_Hdc_Prph"
)
glia <- c(
  "Astrocyte", "Ependymal_cell", "Endothelial_cell",
  "Epithelial-like_cell", "Microglia", "Oligodendrocyte",
  "OPC", "Pericyte", "Pituicyte", "Schwann_cell",
  "Tanycyte", "Undefined", "VLMC"
)
pituitary <- c(
  "Pituitary_Crhr1", "Pituitary_Gnrhr", "Pituitary_Mtnr1a",
  "Pituitary_Pappa2", "Pituitary_Pax7", "Pituitary_Rxfp1"
)
seurat_merged$neurotransmitter <- case_when(
  seurat_merged$detailed_celltype_annotation %in% glut       ~ "Glut",
  seurat_merged$detailed_celltype_annotation %in% gaba      ~ "GABA",
  seurat_merged$detailed_celltype_annotation %in% gaba_glut ~ "GABA_Glut",
  seurat_merged$detailed_celltype_annotation %in% chat      ~ "Chat",
  seurat_merged$detailed_celltype_annotation %in% hist      ~ "Hist",
  seurat_merged$detailed_celltype_annotation %in% pituitary ~ "Pituitary",
  seurat_merged$detailed_celltype_annotation %in% glia      ~ "Glia",
  T ~ NA
)
table(seurat_merged$neurotransmitter, useNA = "always")

###EXCLUDING INJECTION REGION TO DISTINGUISH TRACING FROM DIFFUSION
inj_mask <- (
   (seurat_merged$old.ident == "S1A1" & seurat_merged$region %in% c(
     "bed_nucleus_left", "reuniens_thalamus", "paraventricular_left",
     "anterior_hypothalamus_left", "medial_preoptic_left",
     "lateral_amygdaloid_left", "suprachiasmatic_left"
   )) |
   (seurat_merged$old.ident == "S1D1" & seurat_merged$region %in% c(
     "bed_nucleus", "medial_preoptic", "zona_incerta", "paraventricular"
   )) |
   (seurat_merged$old.ident == "S1D2" & seurat_merged$region %in% c(
     "bed_nucleus", "medial_preoptic", "paraventricular"
   )) |
   (seurat_merged$old.ident == "S2A1" & seurat_merged$region == "lateral_hypothalamus_right") |
   (seurat_merged$old.ident == "S2A2" & seurat_merged$region %in% c(
     "lateral_hypothalamus_right", "anterior_hypothalamus_right"
   )) |
   (seurat_merged$old.ident == "S2B1" & seurat_merged$region %in% c(
     "lateral_hypothalamus_right", "anterior_hypothalamus_right"
   )) |
   (seurat_merged$old.ident == "S2B2" & seurat_merged$region %in% c(
     "lateral_hypothalamus_right", "dorsomedial_right", "zona_incerta_right"
   )) |
   (seurat_merged$old.ident == "S2C1" & seurat_merged$region %in% c(
     "lateral_hypothalamus_right", "dorsomedial_right", "zona_incerta_right",
     "dorsal_hypothalamus", "reuniens_thalamus", "tuberal_right"
   )) |
   (seurat_merged$old.ident == "S2C2" & seurat_merged$region %in% c(
     "lateral_hypothalamus_right", "premammillary_dorsal_right",
     "premammillary_ventral_right", "retromammillary_right",
     "posterior_hypothalamus_right"
   )) |
   (seurat_merged$old.ident == "S2D1" & seurat_merged$region %in% c(
     "posterior_hypothalamus", "retromammillary"
   )) |
   (seurat_merged$old.ident == "S2D2" & seurat_merged$region %in% c(
     "posterior_hypothalamus", "retromammillary"
   ))
)
seurat_merged$injection_filter <- ifelse(inj_mask, "injection", "projection")
table(seurat_merged$injection_filter)

##################SCRAN NORMALIZATION######################
library(SingleCellExperiment)
library(scran)
library(BiocParallel)

data_mat <- GetAssayData(seurat_merged, assay = "raw", layer = "counts")
input_groups = seurat_merged$detailed_celltype_annotation
size_factors = sizeFactors(
    computeSumFactors(
        SingleCellExperiment(
            list(counts=data_mat)), 
            clusters = input_groups,
            min.mean = 0.1,
            BPPARAM = MulticoreParam()
    )
)
# Add the size factors as metadata
seurat_merged$size_factors <- size_factors

plot(
  seurat_merged$size_factors,
  seurat_merged$area,
  pch   = 16,
  cex   = 0.6,
  xlab  = "Size factors",
  ylab  = "Area",
  main  = "Size Factors vs. Area in all cells"
)

# Add linear fit
fit <- lm(area ~ size_factors, data = seurat_merged@meta.data)
abline(fit, col = "red", lwd = 2)

# Optionally display R² and p‐value
s <- summary(fit)
mtext(
  bquote(R^2 == .(format(s$r.squared, digits = 3))),
  side = 3, line = 0.5, adj = 1, col = "blue"
)

# Pull out the raw‐count matrix (genes × cells)
counts_mat <- GetAssayData(seurat_merged, assay="raw", layer="counts")
# Divide each column (cell) by its size factor
# equivalent to: scran = adata.X / adata.obs["size_factors"][:,None]
scran_norm <- sweep(counts_mat,     # matrix
                    2,              # margin = 2 (columns/cells)
                    seurat_merged$size_factors,
                    FUN = "/")
# log1p‐transform and coerce back to a dgCMatrix (CSC)
scran_log1p <- log1p(scran_norm)
scran_log1p <- as(scran_log1p, "dgCMatrix")
# Store the result in your Seurat object as a new “layer” or assay
# as a new assay called “scran”
seurat_merged <- SetAssayData(
  object   = seurat_merged,
  assay    = "raw",
  layer    = "scran_normalization",
  new.data = scran_log1p
)
seurat_merged@assays$raw

total_counts <- colSums(
  GetAssayData(seurat_merged, assay = "raw", layer = "counts")
)
scran_sums <- colSums(
  seurat_merged@assays$raw@layers[["scran_normalization"]]
)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
hist(total_counts,
     breaks = 100,
     main   = "Total counts",
     xlab   = "Counts per cell",
     ylab   = "Number of cells")
hist(scran_sums,
     breaks = 100,
     main   = "log1p w/ Scran size factors",
     xlab   = "Sum of log1p normalized counts",
     ylab   = "Number of cells")

qs::qsave(seurat_merged, paste0(atlaspath, "spatial_24_bipca_annotated_scran.qs"))
seurat_merged

seurat_merged <- qs::qread(paste0(atlaspath, "spatial_24_bipca_annotated_scran.qs"))
seurat_merged

par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
hist(seurat_merged$size_factors,
     breaks = 10000,
     xlab   = "size factors",
     ylab   = "Number of cells")
hist(seurat_merged$area,
     breaks = 10000,
     xlab   = "area",
     ylab   = "Number of cells")

# Define threshold
thr <- 20

# Identify cells where size_factors > 20
outlier_cells <- colnames(seurat_merged)[
  which(!is.na(seurat_merged$size_factors) & seurat_merged$size_factors > thr)
]
length(outlier_cells)
seurat_merged$sf_outlier <- factor(
  ifelse(!is.na(seurat_merged$size_factors) & seurat_merged$size_factors > thr, "outlier", "other"),
  levels = c("other", "outlier")
)
table(seurat_merged$sf_outlier, seurat_merged$old.ident)

hist(seurat_merged$area,
     breaks = 10000,
     xlim   = c(0, 1000),
     xlab   = "area",
     ylab   = "Number of cells")

# Define threshold
thr_upper <- 20000
thr_lower <- 250
# Identify cells
outlier_cells <- colnames(seurat_merged)[
  which(seurat_merged$area < thr_lower | seurat_merged$area > thr_upper)
]
length(outlier_cells)
seurat_merged$ar_outlier <- factor(
  ifelse(!is.na(seurat_merged$size_factors) & seurat_merged$size_factors > thr, "outlier", "other"),
  levels = c("other", "outlier")
)
table(seurat_merged$ar_outlier, seurat_merged$old.ident)

plot <- ImageDimPlot(seurat_merged, 
             border.color = NA, 
             boundaries = "segmentation", 
             group.by = "sf_outlier", 
             cols = c("grey", "firebrick"),
             fov = "S3C1",
             dark.background = FALSE,
            )+ theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      legend.key.size = unit(0.5, "lines"),
                      legend.text = element_text(size = 9), 
                      plot.background = element_rect(fill = "white"),
                      panel.background = element_rect(fill = "white"), 
                      plot.margin = unit(c(0, 0, 0, 0), "cm"))
show(plot)

plot <- ImageDimPlot(seurat_merged, 
             border.color = NA, 
             boundaries = "segmentation", 
             group.by = "ar_outlier", 
             cols = c("grey", "darkgreen"),
             fov = "S3C1",
             dark.background = FALSE,
            )+ theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      legend.key.size = unit(0.5, "lines"),
                      legend.text = element_text(size = 9), 
                      plot.background = element_rect(fill = "white"),
                      panel.background = element_rect(fill = "white"), 
                      plot.margin = unit(c(0, 0, 0, 0), "cm"))
show(plot)

seurat_merged <- subset(seurat_merged, subset = sf_outlier == "outlier" | ar_outlier == "outlier", invert = TRUE)
seurat_merged

plot(
  seurat_merged$size_factors,
  seurat_merged$area,
  pch   = 16,
  cex   = 0.6,
  xlab  = "Size factors",
  ylab  = "Area",
  main  = "Size Factors vs. Area in all cells"
)

# Add linear fit
fit <- lm(area ~ size_factors, data = seurat_merged@meta.data)
abline(fit, col = "red", lwd = 2)

# Optionally display R² and p‐value
s <- summary(fit)
mtext(
  bquote(R^2 == .(format(s$r.squared, digits = 3))),
  side = 3, line = 0.5, adj = 1, col = "blue"
)

qs::qsave(seurat_merged, paste0(atlaspath, "spatial_24_bipca_annotated_scran_qc.qs"))
seurat_merged

###############INJECTION EXCLUDED SUBSET
seurat_projecting <- subset(seurat_merged, subset = injection_filter %in% c("projection"))
qs::qsave(seurat_projecting, paste0(atlaspath, "spatial_24_bipca_annotated_projecting.qs"))
seurat_projecting

###############CORONAL SUBSET
fov_list <- c("A1", "A2", "B1", "B2", "C1", "C2")
seurat_coronal <- subset(seurat_merged, subset = section_index %in% fov_list)
# Create a new metadata column 'region_side' based on suffix
seurat_coronal$region_side <- ifelse(
  grepl("_left$", seurat_coronal$region), "left",
  ifelse(grepl("_right$", seurat_coronal$region), "right", NA)
)
# Remove "_left" or "_right" to extract base region name
seurat_coronal$region_base <- gsub("_(left|right)$", "", seurat_coronal$region)
qs::qsave(seurat_coronal, paste0(atlaspath, "spatial_24_bipca_annotated_coronal.qs"))
seurat_coronal

###############RUN INDEX SUBSET
seurat_index1 <- subset(seurat_merged, subset = run_index == "spatial1")
seurat_index2 <- subset(seurat_merged, subset = run_index == "spatial2")
seurat_index3 <- subset(seurat_merged, subset = run_index == "spatial3")

qs::qsave(seurat_index1, paste0(atlaspath, "spatial_24_bipca_annotated_spatial1.qs"))
qs::qsave(seurat_index2, paste0(atlaspath, "spatial_24_bipca_annotated_spatial2.qs"))
qs::qsave(seurat_index3, paste0(atlaspath, "spatial_24_bipca_annotated_spatial3.qs"))
seurat_index1

#######VARIOUS PLOTTING FUNCTIONS#############

fov_list <- c("S1A1", "S1A2", "S1B1", "S1B2", 
                  "S1C1", "S1C2", "S1D1", "S1D2", 
                  "S2A1", "S2A2", "S2B1", "S2B2", 
                  "S2C1", "S2C2", "S2D1", "S2D2", 
                  "S3A1", "S3A2", "S3B1", "S3B2", 
                  "S3C1", "S3C2", "S3D1", "S3D2")
stack_plot <- function(spatial_object, fov){
plot <- ImageDimPlot(spatial_object, 
             border.color = NA, 
             boundaries = "segmentation", 
             group.by = "injection_filter",
             molecules = c("p147_mCherry_2A_iCre", "p36_NLS_Cre"),
             mols.cols = c("white", "green"),
             mols.size = 0.01,
             mols.alpha = 0.1,
             nmols = 1000000,
             fov = fov,
             cols = godsnot_102)+ theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      legend.key.size = unit(0.5, "lines"),
                      legend.text = element_text(size = 9), 
                      plot.background = element_rect(fill = "black"),
                      panel.background = element_rect(fill = "black"), 
                      plot.margin = unit(c(0, 0, 0, 0), "cm"))
}
for (f in fov_list) {
    plot <- stack_plot(seurat_merged, fov = f)
    resolve_save(plot, paste0(atlaspath, "plots/",f,"_projections_filtered.pdf"))
}

stack_plot <- function(spatial_object, fov){
    plot <- ImageFeaturePlot(spatial_object, 
                         border.color = NA,
                         boundaries = "segmentation", 
                         features = "region_side",
                         fov = fov
            ) + theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(), 
                      plot.background = element_rect(fill = "black"),
                      panel.background = element_rect(fill = "black"), 
                      plot.margin = unit(c(0, 0, 0, 0), "cm") 
                     ) + scale_fill_manual(values = c("red", "blue", "gray"))
}
fov_list <- c("S1A1", "S1A2", "S1B1", "S1B2", 
                  "S1C1", "S1C2", 
                  "S2A1", "S2A2", "S2B1", "S2B2", 
                  "S2C1", "S2C2", 
                  "S3A1", "S3A2", "S3B1", "S3B2", 
                  "S3C1", "S3C2")

for (f in fov_list) {
    plot <- stack_plot(spatial_object, fov = f)
    resolve_save(plot, paste0(atlaspath, "plots/",f,"_leftright.pdf"))
}

fov_list <- c("S1A1", "S1A2", "S1B1", "S1B2", 
                  "S1C1", "S1C2", "S1D1", "S1D2", 
                  "S2A1", "S2A2", "S2B1", "S2B2", 
                  "S2C1", "S2C2", "S2D1", "S2D2", 
                  "S3A1", "S3A2", "S3B1", "S3B2", 
                  "S3C1", "S3C2", "S3D1", "S3D2")
stack_plot <- function(spatial_object, fov){
plot <- ImageDimPlot(spatial_object, 
             border.color = NA, 
             boundaries = "segmentation", 
             group.by = "projection",
             molecules = c("p147_mCherry_2A_iCre", "p36_NLS_Cre"),
             mols.cols = c("cyan", "green"),
             mols.size = 0.01,
             mols.alpha = 0.5,
             nmols = 1000000,
             fov = fov,
             cols = "polychrome")+ theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      legend.key.size = unit(0.5, "lines"),
                      legend.text = element_text(size = 9), 
                      plot.background = element_rect(fill = "black"),
                      panel.background = element_rect(fill = "black"), 
                      plot.margin = unit(c(0, 0, 0, 0), "cm"))
}
#resolve_save(plot, paste0(atlaspath, "plots/testprojection.pdf"))
for (f in fov_list) {
    plot <- stack_plot(seurat_merged, fov = f)
    resolve_save(plot, paste0(atlaspath, "plots/",f,"_projections.pdf"))
}
###S1 PAG v63, BST v147
###S2 PBN v63, LHA v147
###S3 CEA v63, PVT v147

plot <- ImageDimPlot(spatial_object, 
             border.color = NA, 
             boundaries = "segmentation", 
             group.by = "region",
             molecules = c("p147_mCherry_2A_iCre", "p36_NLS_Cre"),
             mols.cols = c("white", "green"),
             mols.size = 0.01,
             mols.alpha = 0.1,
             nmols = 1000000,
             fov = "S1A1",
             cols = godsnot_102)+ theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      legend.key.size = unit(0.5, "lines"),
                      legend.text = element_text(size = 9), 
                      plot.background = element_rect(fill = "black"),
                      panel.background = element_rect(fill = "black"), 
                      plot.margin = unit(c(0, 0, 0, 0), "cm"))
show(plot)

plot <- ImageDimPlot(spatial_object, 
             border.color = NA, 
             boundaries = "segmentation", 
             group.by = "region",
             molecules = c("p147_mCherry_2A_iCre", "p36_NLS_Cre"),
             mols.cols = c("white", "green"),
             mols.size = 0.01,
             mols.alpha = 0.1,
             nmols = 1000000,
             fov = "S2C2",
             cols = godsnot_102)+ theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      legend.key.size = unit(0.5, "lines"),
                      legend.text = element_text(size = 9), 
                      plot.background = element_rect(fill = "black"),
                      panel.background = element_rect(fill = "black"), 
                      plot.margin = unit(c(0, 0, 0, 0), "cm"))
show(plot)
