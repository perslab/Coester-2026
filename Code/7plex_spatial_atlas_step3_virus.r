library(dplyr)
library(Seurat)
library(ggplot2)
library(reshape2)
library(stats)
library(scran)
library(scater)
library(BiocParallel)
library(tidyr)
library(pheatmap)
library(tibble)
library(irr)
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
downstream <- c("01_Undefined", "35_Lhx9_Zic1", "Undefined", 
                "Schwann_cell", "Epithelial-like_cell", 
                "Pituitary_Crhr1", "Pituitary_Gnrhr", "Pituitary_Pappa2")

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

colors_detailed_light <- c('#ffff00', '#1ce6ff', '#ff34ff', '#ff4a46', 
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
                      '#456648', '#0086ed', "#d4d4d4")

names(colors_manual) <- cells_manual
names(colors_detailed) <- cells_detailed
names(colors_region) <- cells_region
names(colors_detailed_light) <- cells_detailed

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

###Loading projecting dataset for v147 calls, this one has injection sites excluded
seurat_merged <- qs::qread(paste0(atlaspath, "spatial_24_bipca_annotated_projecting.qs"))
seurat_merged

#tabulating stats for comparison with snRNAseq
seurat_merged$comp1 <- seurat_merged$detailed_celltype_annotation
seurat_merged@meta.data[["comp1"]] <- ifelse(seurat_merged@meta.data[["score_detailed_celltype_annotation"]] > 0.4, 
                                            seurat_merged@meta.data[["comp1"]], NA)
#projection calls
df <- table(seurat_merged$projection, seurat_merged$comp1)
write.csv(df, "projection_calls_spatial_o.csv", row.names=TRUE)
#total counts
df <- table(seurat_merged$comp1)
write.csv(df, "cells_calls_spatial_o.csv", row.names=TRUE)

###Loading projecting dataset for v147 calls CORONAL
seurat_merged <- qs::qread(paste0(atlaspath, "spatial_24_bipca_annotated_coronal.qs"))
seurat_merged$projection[seurat_merged$projection == "CeA"] <- "CEA"
seurat_merged$projection[seurat_merged$projection == "PVT+CeA"] <- "PVT+CEA"
table(seurat_merged$projection, seurat_merged$run_index)
#tabulating stats for comparison with snRNAseq
seurat_merged$comp1 <- seurat_merged$detailed_celltype_annotation
seurat_merged@meta.data[["comp1"]] <- ifelse(seurat_merged@meta.data[["score_detailed_celltype_annotation"]] > 0.4, 
                                            seurat_merged@meta.data[["comp1"]], NA)
#projection calls
df <- table(seurat_merged$projection, seurat_merged$comp1)
write.csv(df, "projection_calls_spatial_c.csv", row.names=TRUE)
#total counts
df <- table(seurat_merged$comp1)
write.csv(df, "cells_calls_spatial_c.csv", row.names=TRUE)

###Loading projecting dataset for v147 calls, this one has injection sites excluded
seurat_merged <- qs::qread(paste0(atlaspath, "spatial_24_bipca_annotated_projecting.qs"))
viral_gene <- "p147-mCherry-2A-iCre"
glial_label <- "Glia"
expr_all <- seurat_merged@assays$raw$scran_normalization[viral_gene, ]
animal_ids <- unique(seurat_merged$old.ident)
z_scores <- rep(NA, length(expr_all))
names(z_scores) <- names(expr_all)
expr_hist <- expr_all[expr_all > 0]
hist(expr_hist,
     breaks = 100,
     main   = "log1p w/ Scran size factors",
     xlab   = "Sum of log1p normalized counts",
     ylab   = "Number of cells")

stats <- data.frame(
  animal   = character(0),
  n_glia   = integer(0),
  mean_bg  = numeric(0),
  std_bg   = numeric(0),
  stringsAsFactors = FALSE
)
for(animal in animal_ids) {
  mask_animal <- seurat_merged$old.ident == animal
  expr_animal <- expr_all[mask_animal]
  mask_glia <- seurat_merged$neurotransmitter[mask_animal] == glial_label
  n_glia <- sum(mask_glia)
  glia_expr <- expr_animal[mask_glia]
  mean_bg   <- mean(glia_expr)
  std_bg    <- sd(glia_expr)
  stats <- rbind(stats, data.frame(animal, n_glia, mean_bg, std_bg,
                                   stringsAsFactors=FALSE))
  z_scores[mask_animal] <- (expr_animal - mean_bg)/std_bg
}
stats

seurat_merged$viral_zscore_147 <- z_scores
df_stats <- stats
rownames(df_stats) <- df_stats$animal
table(is.na(seurat_merged$viral_zscore_147))

p <- ggplot(df_stats, aes(x = animal, y = mean_bg)) +
  geom_col(fill = "steelblue") +
  geom_errorbar(aes(
      ymin = mean_bg - std_bg,
      ymax = mean_bg + std_bg
    ),
    width = 0.2,
    size = 0.7
  ) +
  labs(
    x = "Animal ID",
    y = "Glial background expression (mean ± SD)",
    title = paste0(viral_gene, " background per Section")
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title  = element_text(hjust = 0.5)
  )
plot(p)

plot <- ImageFeaturePlot(seurat_merged, border.color = NA, 
                         boundaries = "segmentation", 
                         features = "viral_zscore_147",
                         max.cutoff = "q75",
                         fov = "S3B2") + scale_fill_distiller(palette = "Spectral")
show(plot)

# Get viral_zscore values from the subset
zscore_vec <- seurat_merged$viral_zscore_147
zscore_df <- data.frame(viral_zscore_147 = zscore_vec)
rownames(zscore_df) <- names(zscore_vec)

###Loading projecting dataset for v36 calls, this has no regions excluded
seurat_merged <- qs::qread(paste0(atlaspath, "spatial_24_bipca_annotated_scran_qc.qs"))
seurat_merged

# Add v147 zscores to the original object
seurat_merged <- AddMetaData(seurat_merged, metadata = zscore_df)

viral_gene <- "p36-NLS-Cre"
glial_label <- "Glia"
expr_all <- seurat_merged@assays$raw$scran_normalization[viral_gene, ]
animal_ids <- unique(seurat_merged$old.ident)
stats <- data.frame(
  animal   = character(0),
  n_glia   = integer(0),
  mean_bg  = numeric(0),
  std_bg   = numeric(0),
  stringsAsFactors = FALSE
)
z_scores <- rep(NA, length(expr_all))
names(z_scores) <- names(expr_all)

expr_hist <- expr_all[expr_all > 0]
hist(expr_hist,
     breaks = 100,
     main   = "log1p w/ Scran size factors",
     xlab   = "Sum of log1p normalized counts",
     ylab   = "Number of cells")

for(animal in animal_ids) {
  mask_animal <- seurat_merged$old.ident == animal
  expr_animal <- expr_all[mask_animal]
  
  mask_glia <- seurat_merged$neurotransmitter[mask_animal] == glial_label
  n_glia <- sum(mask_glia)
  if(n_glia == 0) { 
    warning("No glia for ", animal); next 
  }
  
  glia_expr <- expr_animal[mask_glia]
  mean_bg   <- mean(glia_expr)
  std_bg    <- sd(glia_expr)
  if(std_bg == 0) {
    warning("Zero‐variance for ", animal); next
  }
  
  stats <- rbind(stats, data.frame(animal, n_glia, mean_bg, std_bg,
                                   stringsAsFactors=FALSE))
  
  z_scores[mask_animal] <- (expr_animal - mean_bg)/std_bg
}
stats

seurat_merged$viral_zscore_36 <- z_scores
df_stats <- stats
rownames(df_stats) <- df_stats$animal
table(is.na(seurat_merged$viral_zscore_36))

p <- ggplot(df_stats, aes(x = animal, y = mean_bg)) +
  geom_col(fill = "steelblue") +
  geom_errorbar(aes(
      ymin = mean_bg - std_bg,
      ymax = mean_bg + std_bg
    ),
    width = 0.2,
    size = 0.7
  ) +
  labs(
    x = "Animal ID",
    y = "Glial background expression (mean ± SD)",
    title = paste0(viral_gene, " background per Section")
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title  = element_text(hjust = 0.5)
  )
plot(p)

seurat_merged@meta.data <- seurat_merged@meta.data %>%
  mutate(
    projection_zscore = case_when(

      # spatial1: both positive → PAG+BST
      run_index == "spatial1" & as.numeric(viral_zscore_36)  >= 2 &
                                as.numeric(viral_zscore_147) >= 2             ~ "BST+PAG",
      run_index == "spatial1" & as.numeric(viral_zscore_36)  >= 2             ~ "PAG",
      run_index == "spatial1" & as.numeric(viral_zscore_147) >= 2             ~ "BST",
      
      # spatial2: both positive → PBN+LHA
      run_index == "spatial2" & as.numeric(viral_zscore_36)  >= 2 &
                                as.numeric(viral_zscore_147) >= 2             ~ "LHA+PBN",
      run_index == "spatial2" & as.numeric(viral_zscore_36)  >= 2             ~ "PBN",
      run_index == "spatial2" & as.numeric(viral_zscore_147) >= 2             ~ "LHA",
      
      # spatial3: both positive → CEA+PVT
      run_index == "spatial3" & as.numeric(viral_zscore_36)  >= 2 &
                                as.numeric(viral_zscore_147) >= 2             ~ "PVT+CEA",
      run_index == "spatial3" & as.numeric(viral_zscore_36)  >= 2             ~ "CEA",
      run_index == "spatial3" & as.numeric(viral_zscore_147) >= 2             ~ "PVT",

      is.na(viral_zscore_147)                                                ~ "excluded",
      
      # default
      TRUE                                                                     ~ "negative"
    )
  )
# Quick check
table(seurat_merged$projection_zscore, seurat_merged$run_index, useNA = "always")

seurat_merged$projection[seurat_merged$projection == "CeA"] <- "CEA"
seurat_merged$projection[seurat_merged$projection == "PVT+CeA"] <- "PVT+CEA"
table(seurat_merged$projection, seurat_merged$run_index)

ct <- table(
  seurat_merged$projection,
  seurat_merged$projection_zscore)
ct <- ct[rownames(ct) != "excluded", colnames(ct) != "excluded"]
ct <- ct[rownames(ct) != "negative", colnames(ct) != "negative"]
mat <- as.matrix(ct)
df_long <- melt(mat, 
                varnames = c("projection", "zscore"), 
                value.name = "count")

# 4b. Plot
ggplot(df_long, aes(x = zscore, y = projection, fill = count)) +
  geom_tile(color = "grey70") +
  geom_text(aes(label = count), size = 3) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(
    x    = "Z-score category",
    y    = "Projection",
    fill = "Count",
    title= "Projection vs. Z-score"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

#tabulating stats for comparison with snRNAseq
seurat_merged$comp1 <- seurat_merged$detailed_celltype_annotation
seurat_merged@meta.data[["comp1"]] <- ifelse(seurat_merged@meta.data[["score_detailed_celltype_annotation"]] > 0.4, 
                                            seurat_merged@meta.data[["comp1"]], NA)
#projection calls
df <- table(seurat_merged$projection_zscore, seurat_merged$comp1)
write.csv(df, "projection_calls_spatial_z.csv", row.names=TRUE)
#total counts
df <- table(seurat_merged$comp1)
write.csv(df, "cells_calls_spatial_z.csv", row.names=TRUE)

plot <- ImageDimPlot(seurat_merged, 
             border.color = NA, 
             boundaries = "segmentation", 
             group.by = "projection",
             molecules = c("p147_mCherry_2A_iCre", "p36_NLS_Cre"),
             mols.cols = c("cyan", "green"),
             mols.size = 0.01,
             mols.alpha = 0.1,
             nmols=1000000,
             fov = "S3B2",
             cols = "polychrome")+ theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      legend.key.size = unit(0.5, "lines"),
                      legend.text = element_text(size = 9), 
                      plot.background = element_rect(fill = "black"),
                      panel.background = element_rect(fill = "black"), 
                      plot.margin = unit(c(0, 0, 0, 0), "cm"))
show(plot)
resolve_save(plot, paste0(atlaspath, "plots/test_projection.pdf"))

plot <- ImageDimPlot(seurat_merged, 
             border.color = NA, 
             boundaries = "segmentation", 
             group.by = "projection_zscore",
             molecules = c("p147_mCherry_2A_iCre", "p36_NLS_Cre"),
             mols.cols = c("cyan", "green"),
             mols.size = 0.01,
             mols.alpha = 0.1,
             nmols=1000000,
             fov = "S3B2",
             cols = "polychrome")+ theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      legend.key.size = unit(0.5, "lines"),
                      legend.text = element_text(size = 9), 
                      plot.background = element_rect(fill = "black"),
                      panel.background = element_rect(fill = "black"), 
                      plot.margin = unit(c(0, 0, 0, 0), "cm"))
show(plot)
resolve_save(plot, paste0(atlaspath, "plots/test_projection_z.pdf"))

qs::qsave(seurat_merged, paste0(atlaspath, "spatial_24_bipca_viruscalls.qs"))
#seurat_merged <- qs::qread(paste0(atlaspath, "spatial_24_bipca_viruscalls.qs"))
seurat_merged

library(viridisLite)
rocket_pal  <- viridisLite::rocket(100)
rocket_r_pal <- rev(rocket_pal)
# 1. Filter out excluded cells
df <- seurat_merged@meta.data %>%
  filter(projection_zscore != "excluded") %>% filter (detailed_celltype_annotation %!in% downstream)

# 2. Tabulate counts of projection_zscore vs. celltype
count_table <- df %>%
  group_by(detailed_celltype_annotation, projection_zscore) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = projection_zscore, values_from = count, values_fill = 0)

# 3. Calculate enrichment ratios: each column vs. "negative"
projection_labels <- setdiff(names(count_table)[-1], "negative")

# Add ratio columns
ratio_table <- count_table %>%
  mutate(across(all_of(projection_labels), ~ .x / negative, .names = "ratio_{.col}"))

# 4. Extract ratio matrix and prepare for heatmap
ratio_mat <- ratio_table %>%
  select(detailed_celltype_annotation, starts_with("ratio_")) %>%
  column_to_rownames("detailed_celltype_annotation") %>%
  as.matrix()

# 5. Clean column names
colnames(ratio_mat) <- gsub("^ratio_", "", colnames(ratio_mat))

# 6. Plot heatmap
pheatmap(ratio_mat,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = rocket_r_pal,
         fontsize_row = 10,
         fontsize_col = 10,
         display_numbers = TRUE,
         number_format = "%.2f",
         cellwidth = 30,
         cellheight = 30,
         main = "Projection Z-Score Enrichment vs. Negative (by Cell Type)",
         na_col = "gray90",
         filename = paste0(atlaspath,"heatmap_projection_enrichment.png"))

#tabulating stats for comparison with snRNAseq
seurat_merged$comp1 <- seurat_merged$detailed_celltype_annotation
seurat_merged@meta.data[["comp1"]] <- ifelse(seurat_merged@meta.data[["score_detailed_celltype_annotation"]] > 0.4, 
                                            seurat_merged@meta.data[["comp1"]], "excluded")
# 1. Filter out excluded cells
df <- seurat_merged@meta.data %>%
  filter(projection_zscore != "excluded") %>% filter (comp1 %!in% downstream) %>% filter(comp1 != "excluded")

# 2. Tabulate counts of projection_zscore vs. celltype
count_table <- df %>%
  group_by(comp1, projection_zscore) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = projection_zscore, values_from = count, values_fill = 0)

# 3. Calculate enrichment ratios: each column vs. "negative"
projection_labels <- setdiff(names(count_table)[-1], "negative")

# Add ratio columns
ratio_table <- count_table %>%
  mutate(across(all_of(projection_labels), ~ .x / negative, .names = "ratio_{.col}"))

# 4. Extract ratio matrix and prepare for heatmap
ratio_mat <- ratio_table %>%
  select(comp1, starts_with("ratio_")) %>%
  column_to_rownames("comp1") %>%
  as.matrix()

# 5. Clean column names
colnames(ratio_mat) <- gsub("^ratio_", "", colnames(ratio_mat))

# 6. Plot heatmap
pheatmap(ratio_mat,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = rocket_r_pal,
         fontsize_row = 10,
         fontsize_col = 10,
         display_numbers = TRUE,
         number_format = "%.2f",
         cellwidth = 30,
         cellheight = 30,
         main = "Projection Z-Score Enrichment vs. Negative (by Cell Type)",
         na_col = "gray90",
         filename = paste0(atlaspath,"heatmap_projection_enrichment_04.png"))

library(reticulate)
library(data.table)
pd <- import("pandas", convert = FALSE)
outdir <- paste0(atlaspath,"area_pickles")

for (run in unique(seurat_merged$run_index)) {
  #Subset by run_index
  subpath <- paste0(atlaspath, "spatial_24_bipca_annotated_", run, ".qs")
  sub <- qs::qread(subpath)
  sub$comp1 <- sub$detailed_celltype_annotation
  sub@meta.data[["comp1"]] <- ifelse(
      sub@meta.data[["score_detailed_celltype_annotation"]] > 0.4,
      sub@meta.data[["comp1"]], "no_label")
  #Extract the scran‐normalized counts matrix
  counts_mat <- sub@assays$raw@layers[["scran_normalization"]]
  counts_mat <- t(counts_mat)
  #make sure it’s a dense matrix
  counts_mat <- as.matrix(counts_mat)
  
  #Wrap in pandas.DataFrame (cells × genes)
  df_counts <- pd$DataFrame(
    counts_mat,
    index   = rownames(counts_mat),
    columns = colnames(counts_mat)
  )
  
  #Extract metadata columns
  meta1 <- pd$DataFrame(
    sub@meta.data[ , "comp1", drop = FALSE],
    index = rownames(sub@meta.data)
  )
  
  #Write out pickles
  df_counts$to_pickle(file.path(outdir, paste0(run, "_counts.pkl")))
  meta1     $to_pickle(file.path(outdir, paste0(run, "_detailed.pkl")))
  
  message("Saved pickles for run_index = ", run)
}

counts_mat <- seurat_merged@assays$raw@layers[["scran_normalization"]]
counts_mat <- t(counts_mat)

#tabulating stats for comparison with snRNAseq
seurat_merged$comp1 <- seurat_merged$detailed_celltype_annotation
seurat_merged@meta.data[["comp1"]] <- ifelse(seurat_merged@meta.data[["score_detailed_celltype_annotation"]] > 0.4, 
                                            seurat_merged@meta.data[["comp1"]], NA)
seurat_merged$sun1_num <- as.integer(seurat_merged$sun1)
#projection calls
df <- table(seurat_merged$sun1_num, seurat_merged$comp1)
write.csv(df, "projection_calls_spatial_sun.csv", row.names=TRUE)
#total counts
df <- table(seurat_merged$comp1)
write.csv(df, "cells_calls_spatial_sun.csv", row.names=TRUE)

table(seurat_merged$sun1_num, seurat_merged$projection_zscore)

sun1_num <- as.integer(seurat_merged$sun1)
test_147 <- cor.test(
  sun1_num,
  seurat_merged$viral_zscore_147,
  use = "complete.obs",
  method = "pearson"
)
test_36  <- cor.test(
  sun1_num,
  seurat_merged$viral_zscore_36,
  use = "complete.obs",
  method = "pearson"
)

# Inspect results
test_147
test_36

# 1) Create two new metadata columns, TRUE when z ≥ 2
seurat_merged$z36_pos  <- seurat_merged$viral_zscore_36  >= 2
seurat_merged$z147_pos <- seurat_merged$viral_zscore_147 >= 2

# 2) View the 2×2 contingency with sun1
table(sun1 = seurat_merged$sun1, z36_pos  = seurat_merged$z36_pos)
table(sun1 = seurat_merged$sun1, z147_pos = seurat_merged$z147_pos)

# 3) (Optional) Compute a “phi” correlation between two binaries
phi_36  <- cor(as.integer(seurat_merged$sun1),
               as.integer(seurat_merged$z36_pos),
               use = "complete.obs")
phi_147 <- cor(as.integer(seurat_merged$sun1),
               as.integer(seurat_merged$z147_pos),
               use = "complete.obs")
phi_36; phi_147

#Cohen’s kappa for agreement

kappa36  <- kappa2(data.frame(sun1 = seurat_merged$sun1,
                              z36  = seurat_merged$z36_pos))
kappa147 <- kappa2(data.frame(sun1 = seurat_merged$sun1,
                              z147 = seurat_merged$z147_pos))
kappa36; kappa147

