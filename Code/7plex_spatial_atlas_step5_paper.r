library(dplyr)
library(Seurat)
library(ggplot2)
library(reshape2)
library(tidyr)
library(pheatmap)
library(tibble)
library(irr)
'%!in%' <- function(x,y)!('%in%'(x,y))
options(repr.plot.width = 12, repr.plot.height = 12)
spatialpath <- "/projects/perslab/people/qbx224/7plex_paper/atlas/"
atlaspath <- "/projects/perslab/people/qbx224/7plex_paper/atlas/paper/"
options(future.globals.maxSize= 50*1024*1024^2)
setwd(spatialpath)
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

names(colors_detailed) <- cells_detailed
names(colors_detailed_light) <- cells_detailed

plot_spatial <- function(seurat_obj, 
                        fov,
                        group.by = "plot",
                        colors = colors_detailed,
                        save_pdf = FALSE,
                        output_dir = NULL,
                        show_legend = TRUE,
                        rotation = 0) {
  
  # Validate rotation
  if (!rotation %in% c(0, 90, 180, 270)) {
    stop("rotation must be 0, 90, 180, or 270")
  }
  
  # Get bregma value for cells in this FOV
  fov_mask <- seurat_obj$old.ident == fov
  bregma_val <- unique(seurat_obj$bregma[fov_mask])
  
  # Handle case where there might be multiple bregma values (take first)
  if (length(bregma_val) > 1) {
    warning("Multiple bregma values found for FOV ", fov, ". Using first value.")
    bregma_val <- bregma_val[1]
  }
  
  # Get levels and sort alphabetically for legend ordering
  current_levels <- levels(seurat_obj@meta.data[[group.by]])
  if (is.null(current_levels)) {
    current_levels <- unique(seurat_obj@meta.data[[group.by]])
  }
  sorted_levels <- sort(current_levels)
  
  # Reorder the factor in the metadata
  seurat_obj@meta.data[[group.by]] <- factor(seurat_obj@meta.data[[group.by]], 
                                              levels = sorted_levels)
  
  # Create base plot
  plot <- ImageDimPlot(seurat_obj, 
                       border.color = NA, 
                       boundaries = "segmentation", 
                       group.by = group.by, 
                       cols = colors,
                       fov = fov,
                       dark.background = TRUE) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.key.size = unit(0.5, "lines"),
          legend.text = element_text(size = 9), 
          plot.background = element_rect(fill = "black"),
          panel.background = element_rect(fill = "black"), 
          plot.margin = unit(c(0, 0, 0, 0), "cm"))
  
  # Hide legend if requested
  if (!show_legend) {
    plot <- plot + theme(legend.position = "none")
  }
  
  # Display plot
  show(plot)
  
  # Save if requested
  if (save_pdf) {
    if (is.null(output_dir)) {
      stop("output_dir must be provided when save_pdf = TRUE")
    }
    # Convert group.by to string for filename
    groupby_str <- as.character(group.by)
    # Create filename with fov and bregma value
    bregma_str <- sprintf("%.2f", bregma_val)
    filename <- paste0(fov, "_bregma_", bregma_str, "_", groupby_str, ".pdf")
    filepath <- file.path(output_dir, filename)
    
    # Save with rotation applied via viewport if needed
    cairo_pdf(filepath, width = 12, height = 12, bg = "transparent")
    grid::grid.newpage()
    
    if (rotation != 0) {
      vp <- grid::viewport(angle = rotation)
      grid::pushViewport(vp)
      print(plot, newpage = FALSE)
      grid::popViewport()
    } else {
      print(plot, newpage = FALSE)
    }
    
    dev.off()
    
    message("Plot saved to: ", filepath)
  }
  
  return(invisible(plot))
}

# Plot continuous numeric features on spatial data
# Similar to plot_spatial but for numeric/continuous variables

plot_spatial_continuous <- function(seurat_obj, 
                                   fov,
                                   feature,
                                   palette = "Spectral",
                                   palette_direction = -1,
                                   save_pdf = FALSE,
                                   output_dir = NULL,
                                   show_legend = TRUE,
                                   rotation = 0,
                                   limits = NULL,
                                   na_value = "grey50") {
  
  # Validate rotation
  if (!rotation %in% c(0, 90, 180, 270)) {
    stop("rotation must be 0, 90, 180, or 270")
  }
  
  # Get bregma value for cells in this FOV
  fov_mask <- seurat_obj$old.ident == fov
  bregma_val <- unique(seurat_obj$bregma[fov_mask])
  
  # Handle case where there might be multiple bregma values (take first)
  if (length(bregma_val) > 1) {
    warning("Multiple bregma values found for FOV ", fov, ". Using first value.")
    bregma_val <- bregma_val[1]
  }
  
  # Create base plot
  plot <- ImageFeaturePlot(seurat_obj, 
                          border.color = NA, 
                          boundaries = "segmentation", 
                          features = feature,
                          fov = fov,
                          dark.background = TRUE) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.key.size = unit(0.5, "lines"),
          legend.text = element_text(size = 9), 
          legend.title = element_text(size = 10),
          plot.background = element_rect(fill = "black"),
          panel.background = element_rect(fill = "black"), 
          plot.margin = unit(c(0, 0, 0, 0), "cm"))
  
  # Apply color palette
  if (!is.null(limits)) {
    plot <- plot + scale_fill_distiller(
      palette = palette,
      direction = palette_direction,
      limits = limits,
      na.value = na_value
    )
  } else {
    plot <- plot + scale_fill_distiller(
      palette = palette,
      direction = palette_direction,
      na.value = na_value
    )
  }
  
  # Hide legend if requested
  if (!show_legend) {
    plot <- plot + theme(legend.position = "none")
  }
  
  # Display plot
  show(plot)
  
  # Save if requested
  if (save_pdf) {
    if (is.null(output_dir)) {
      stop("output_dir must be provided when save_pdf = TRUE")
    }
    
    # Convert feature to string for filename
    feature_str <- as.character(feature)
    
    # Create filename with fov, bregma value, and feature
    bregma_str <- sprintf("%.2f", bregma_val)
    filename <- paste0(fov, "_bregma_", bregma_str, "_", feature_str, ".pdf")
    filepath <- file.path(output_dir, filename)
    
    # Save with rotation applied via viewport if needed
    cairo_pdf(filepath, width = 12, height = 12, bg = "transparent")
    grid::grid.newpage()
    
    if (rotation != 0) {
      vp <- grid::viewport(angle = rotation)
      grid::pushViewport(vp)
      print(plot, newpage = FALSE)
      grid::popViewport()
    } else {
      print(plot, newpage = FALSE)
    }
    
    dev.off()
    
    message("Plot saved to: ", filepath)
  }
  
  return(invisible(plot))
}

# Usage examples:

# Basic usage with default Spectral palette (reversed)
# plot_spatial_continuous(seurat_merged, 
#                         fov = "S3C2", 
#                         feature = "cellect_bmi_score")

# Save to PDF
# plot_spatial_continuous(seurat_merged, 
#                         fov = "S3C2", 
#                         feature = "cellect_bmi_score",
#                         save_pdf = TRUE,
#                         output_dir = atlaspath)

# Use different palette (viridis-style)
# plot_spatial_continuous(seurat_merged, 
#                         fov = "S3C2", 
#                         feature = "cellect_bmi_score",
#                         palette = "RdYlBu",
#                         palette_direction = 1)

# Available RColorBrewer palettes for scale_fill_distiller:
# Sequential: Blues, BuGn, BuPu, GnBu, Greens, Greys, Oranges, OrRd, 
#             PuBu, PuBuGn, PuRd, Purples, RdPu, Reds, YlGn, YlGnBu, YlOrBr, YlOrRd
# Diverging: BrBG, PiYG, PRGn, PuOr, RdBu, RdGy, RdYlBu, RdYlGn, Spectral

# Set specific limits
# plot_spatial_continuous(seurat_merged, 
#                         fov = "S3C2", 
#                         feature = "cellect_bmi_score",
#                         palette = "Spectral",
#                         limits = c(0, 5))

# With rotation
# plot_spatial_continuous(seurat_merged, 
#                         fov = "S3C2", 
#                         feature = "cellect_bmi_score",
#                         rotation = 90,
#                         save_pdf = TRUE,
#                         output_dir = atlaspath)

# Hide legend
# plot_spatial_continuous(seurat_merged, 
#                         fov = "S3C2", 
#                         feature = "cellect_bmi_score",
#                         show_legend = FALSE)

seurat_merged <- qs::qread(paste0(spatialpath, "spatial_24_bipca_downstream.qs"))
seurat_merged

unique(rownames(seurat_merged@assays$raw))

seurat_merged$detailed_celltype <- seurat_merged$detailed_celltype_annotation
seurat_merged@meta.data[["detailed_celltype"]] <- ifelse(seurat_merged@meta.data[["score_detailed_celltype_annotation"]] > 0.4, 
                                                         seurat_merged@meta.data[["detailed_celltype"]], "no_label")
plot_spatial(seurat_merged, 
             fov = "S3C2",
             group.by = "detailed_celltype",
             save_pdf = TRUE,
             output_dir = atlaspath,
             show_legend = TRUE,
             rotation = 0)

plot_spatial(seurat_merged, 
             fov = "S3B2",
             group.by = "detailed_celltype",
             save_pdf = TRUE,
             output_dir = atlaspath,
             show_legend = FALSE,
             rotation = 270)

plot_legend_only <- function(seurat_obj,
                             group.by = "detailed_celltype_annotation",
                             colors = colors_detailed,
                             output_path,
                             ncol = 1) {
  
  # Get levels and sort alphabetically
  current_levels <- levels(seurat_obj@meta.data[[group.by]])
  if (is.null(current_levels)) {
    current_levels <- unique(seurat_obj@meta.data[[group.by]])
  }
  sorted_levels <- sort(current_levels)
  
  # Move "no_label" to the end if it exists
  if ("no_label" %in% sorted_levels) {
    sorted_levels <- c(sorted_levels[sorted_levels != "no_label"], "no_label")
  }
  
  # Create a minimal dummy plot just to extract the legend
  dummy_data <- data.frame(
    x = 1:length(sorted_levels),
    y = 1:length(sorted_levels),
    category = factor(sorted_levels, levels = sorted_levels)
  )
  
  p <- ggplot(dummy_data, aes(x = x, y = y, color = category)) +
    geom_point(size = 4) +
    scale_color_manual(
      name = NULL,
      values = colors[sorted_levels],
      breaks = sorted_levels
    ) +
    theme_void() +
    theme(
      legend.key.size = unit(0.5, "lines"),
      legend.text = element_text(size = 9, color = "white"),
      legend.background = element_rect(fill = "black", color = NA),
      plot.background = element_rect(fill = "black", color = NA),
      panel.background = element_rect(fill = "black", color = NA)
    ) +
    guides(color = guide_legend(ncol = ncol, override.aes = list(size = 3)))
  
  # Extract just the legend
  legend <- cowplot::get_legend(p)
  
  # Create a plot with only the legend
  legend_plot <- cowplot::ggdraw() + 
    cowplot::draw_plot(legend) +
    theme(plot.background = element_rect(fill = "black", color = NA))
  
  # Save the legend
  ggsave(output_path,
         plot = legend_plot,
         device = cairo_pdf,
         width = 4,
         height = 8,
         units = "in",
         dpi = 300,
         bg = "black")
  
  message("Legend saved to: ", output_path)
  
  return(invisible(legend_plot))
}

# Usage:
plot_legend_only(seurat_merged,
                 group.by = "plot",
                 colors = colors_detailed,
                 output_path = paste0(atlaspath, "legend_only.pdf"),
                 ncol = 1)

glia <- c("GABA_Glut")
seurat_merged$NT_filtered <- seurat_merged$neurotransmitter
seurat_merged@meta.data[["NT_filtered"]] <- ifelse(seurat_merged@meta.data[["NT_filtered"]] %!in% glia, seurat_merged@meta.data[["NT_filtered"]], "no_label")
cols <- c('#FF4A46', '#1CE6FF', '#FFFF00', '#FF34FF', '#006FA6', '#008941', "#242424")
cells <- c("Glut", "GABA", "Chat", "Hist", "Glia", "Pituitary", "no_label")
names(cols) <- cells
plot_spatial(seurat_merged, 
             fov = "S3B2",
             group.by = "NT_filtered",
             colors = cols,
             save_pdf = TRUE,
             output_dir = atlaspath,
             show_legend = TRUE,
             rotation = 0)

plot_spatial(seurat_merged, 
             fov = "S3B2",
             group.by = "NT_filtered",
             colors = cols,
             save_pdf = TRUE,
             output_dir = atlaspath,
             show_legend = FALSE,
             rotation = 270)

plot_spatial_continuous(seurat_merged, 
                        fov = "S3B2", 
                        feature = "cellect_bmi_score",
                        palette = "Purples",
                        palette_direction = -1,
                        save_pdf = TRUE,
                        output_dir = atlaspath)

plot_spatial_continuous(seurat_merged, 
                        fov = "S3B2", 
                        feature = "cellect_bmi_score",
                        palette = "Purples",
                        palette_direction = -1,
                        save_pdf = TRUE,
                        output_dir = atlaspath,
                        rotation = 270)

cols <- c('#efa471', '#5bc2af', '#FF2F80', "#666666")
cells <- c("CEA", "PVT", "PVT+CEA", "negative")
names(cols) <- cells
plot_spatial(seurat_merged, 
             fov = "S3B2",
             group.by = "projection",
             colors = cols,
             save_pdf = TRUE,
             output_dir = atlaspath,
             show_legend = TRUE,
             rotation = 270)

seurat_merged$plot <- "Cell boundaries"
plot <- ImageDimPlot(seurat_merged, 
             border.color = NA, 
             boundaries = "segmentation",
             group.by = "plot",
             cols = "#666666",
             molecules = c("p36_NLS_Cre", "p147_mCherry_2A_iCre"),
             mols.cols = c("#efa471", "#5bc2af"),
             mols.size = 1,
             mols.alpha = 0.5,
             nmols=1000000,
             fov = "S3B2",
             dark.background = TRUE,
            )+ theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      legend.key.size = unit(0.5, "lines"),
                      legend.text = element_text(size = 9), 
                      plot.background = element_rect(fill = "black"),
                      panel.background = element_rect(fill = "black"), 
                      plot.margin = unit(c(0, 0, 0, 0), "cm"))
show(plot)
resolve_save(plot, paste0(atlaspath, "AAV_plot.pdf"))

# Spatial CELLEX - R translation
# Build mapping of projection areas to cell type enrichment scores

indir <- "area_cellex_results"
runs <- c("spatial1", "spatial2", "spatial3")
genes <- c("p36-NLS-Cre", "p147-mCherry-2A-iCre")

# Map each (run → gene) to your desired area label
gene_area_map <- list(
  spatial1 = list(
    "p36-NLS-Cre" = "PAG",
    "p147-mCherry-2A-iCre" = "BST"
  ),
  spatial2 = list(
    "p36-NLS-Cre" = "PBN",
    "p147-mCherry-2A-iCre" = "LHA"
  ),
  spatial3 = list(
    "p36-NLS-Cre" = "CEA",
    "p147-mCherry-2A-iCre" = "PVT"
  )
)

records <- list()

for (run in runs) {
  fn <- paste0(run, "_detailed_cellex.csv")
  path <- file.path(indir, fn)
  
  # Read CSV with first column as row names (rows=genes, cols=cell-types)
  df <- read.csv(path, row.names = 1, check.names = FALSE)
  
  for (gene in genes) {
    # Check if gene exists in the data
    if (!gene %in% rownames(df)) {
      stop(sprintf("Gene '%s' not found in %s", gene, fn))
    }
    
    # Get the area name for this run/gene combination
    area <- gene_area_map[[run]][[gene]]
    
    # Extract that row (cell-type → enrichment) as a data frame
    gene_row <- df[gene, , drop = FALSE]
    rownames(gene_row) <- area
    
    # Add to records list
    records[[area]] <- gene_row
  }
}

# Get all unique column names across all data frames
all_columns <- unique(unlist(lapply(records, colnames)))

# Make sure each record has all columns (fill missing with 0)
records_filled <- lapply(records, function(df) {
  missing_cols <- setdiff(all_columns, colnames(df))
  for (col in missing_cols) {
    df[[col]] <- 0
  }
  # Reorder to match all_columns
  df[, all_columns, drop = FALSE]
})

# Build combined data frame: rows=areas, cols=cell-types
combined_spatial <- do.call(rbind, records_filled)

# Ensure you have all your cell-types in the right order
celltypes <- sort(unique(seurat_merged$detailed_celltype_annotation))

# Reindex columns to match cell types, fill missing with 0
missing_cols <- setdiff(celltypes, colnames(combined_spatial))
for (col in missing_cols) {
  combined_spatial[[col]] <- 0
}

# Reorder columns to match celltypes order
combined_spatial <- combined_spatial[, celltypes, drop = FALSE]

cat("\n=== Spatial CELLEX Combined Matrix ===\n")
cat(sprintf("Dimensions: %d areas × %d cell types\n", nrow(combined_spatial), ncol(combined_spatial)))
cat("\nAreas (rows):\n")
print(rownames(combined_spatial))
cat("\nFirst few columns:\n")
print(head(combined_spatial[, 1:min(5, ncol(combined_spatial))]))

# Now add these scores to seurat_merged metadata
# For each area, create a column: seurat_merged$cellex_spatial_{area}

areas <- rownames(combined_spatial)

cat("\n=== Adding spatial CELLEX scores to Seurat metadata ===\n")

for (area in areas) {
  # Create column name
  col_name <- paste0("cellex_spatial_", area)
  
  # Get the scores for this area (one value per cell type)
  area_scores <- combined_spatial[area, ]
  
  # Get cell type annotation for each cell
  cell_types <- seurat_merged$detailed_celltype_annotation
  
  # Create a named vector for mapping (cell type -> score)
  score_lookup <- as.numeric(area_scores)
  names(score_lookup) <- colnames(combined_spatial)
  
  # Map each cell's type to its score
  scores <- score_lookup[as.character(cell_types)]
  
  # Replace NAs with 0
  scores[is.na(scores)] <- 0
  
  # Remove names before adding to metadata
  scores <- unname(scores)
  
  # Add to Seurat metadata
  seurat_merged[[col_name]] <- scores
  
  cat(sprintf("✓ Added %s (range: %.4f to %.4f)\n", 
              col_name, 
              min(scores, na.rm = TRUE), 
              max(scores, na.rm = TRUE)))
}

cat("\n=== Summary ===\n")
cat(sprintf("Added %d new metadata columns:\n", length(areas)))
cat(paste("  -", paste0("cellex_spatial_", areas), collapse = "\n"))

# Show example of what was added
cat("\n=== Example: First few cells ===\n")
example_cols <- paste0("cellex_spatial_", areas[1:min(3, length(areas))])
print(head(seurat_merged@meta.data[, c("detailed_celltype_annotation", example_cols)]))

# Save the combined spatial matrix for reference
write.csv(combined_spatial, "combined_spatial_cellex.csv", row.names = TRUE)

cat("\n✓ Spatial CELLEX scores added to Seurat metadata!\n")
cat("✓ Combined matrix saved to: combined_spatial_cellex.csv\n")

# Benchmark CELLEX spatial scores against ground truth projections
# Find optimal cutoff values for each projection area

library(dplyr)
library(ggplot2)
library(tidyr)
library(pROC)

# Define projection areas and their corresponding CELLEX columns
areas <- c("PAG", "BST", "PBN", "LHA", "CEA", "PVT")
cellex_cols <- paste0("cellex_spatial_", areas)

# Get ground truth projections
ground_truth <- seurat_merged$projection

# Create binary ground truth for each area
# For each area, cells are positive if they have that projection (alone or in combination)
create_binary_truth <- function(projection_col, area) {
  grepl(area, projection_col, fixed = TRUE)
}

cat("=== Benchmark CELLEX scores vs Ground Truth Projections ===\n\n")

# Store results for each area
results_list <- list()
optimal_cutoffs <- list()

for (i in seq_along(areas)) {
  area <- areas[i]
  cellex_col <- cellex_cols[i]
  
  cat(sprintf("\n--- Analyzing %s ---\n", area))
  
  # Get CELLEX scores for this area
  scores <- seurat_merged@meta.data[[cellex_col]]
  
  # Create binary ground truth (TRUE if cell projects to this area)
  truth <- create_binary_truth(ground_truth, area)
  
  # Remove cells with "negative" projection for cleaner analysis
  valid_cells <- ground_truth != "negative"
  scores_valid <- scores[valid_cells]
  truth_valid <- truth[valid_cells]
  
  cat(sprintf("Total cells: %d\n", length(truth_valid)))
  cat(sprintf("Positive cells (project to %s): %d (%.1f%%)\n", 
              area, sum(truth_valid), 100 * mean(truth_valid)))
  
  # Calculate ROC curve
  roc_obj <- roc(truth_valid, scores_valid, quiet = TRUE)
  auc_value <- auc(roc_obj)
  
  cat(sprintf("AUC: %.4f\n", auc_value))
  
  # Find optimal cutoff using Youden's J statistic
  coords_obj <- coords(roc_obj, "best", ret = c("threshold", "sensitivity", "specificity", 
                                                  "accuracy", "precision", "recall"))
  optimal_threshold <- coords_obj$threshold
  
  cat(sprintf("Optimal threshold: %.4f\n", optimal_threshold))
  cat(sprintf("  Sensitivity: %.4f\n", coords_obj$sensitivity))
  cat(sprintf("  Specificity: %.4f\n", coords_obj$specificity))
  cat(sprintf("  Accuracy: %.4f\n", coords_obj$accuracy))
  
  # Test multiple threshold values
  thresholds <- seq(min(scores_valid), max(scores_valid), length.out = 100)
  
  metrics <- data.frame(
    threshold = thresholds,
    accuracy = NA,
    sensitivity = NA,
    specificity = NA,
    precision = NA,
    f1_score = NA
  )
  
  for (j in seq_along(thresholds)) {
    thresh <- thresholds[j]
    predicted <- scores_valid >= thresh
    
    tp <- sum(predicted & truth_valid)
    tn <- sum(!predicted & !truth_valid)
    fp <- sum(predicted & !truth_valid)
    fn <- sum(!predicted & truth_valid)
    
    metrics$accuracy[j] <- (tp + tn) / length(truth_valid)
    metrics$sensitivity[j] <- if (tp + fn > 0) tp / (tp + fn) else 0
    metrics$specificity[j] <- if (tn + fp > 0) tn / (tn + fp) else 0
    metrics$precision[j] <- if (tp + fp > 0) tp / (tp + fp) else 0
    
    # F1 score
    if (metrics$precision[j] + metrics$sensitivity[j] > 0) {
      metrics$f1_score[j] <- 2 * (metrics$precision[j] * metrics$sensitivity[j]) / 
                                  (metrics$precision[j] + metrics$sensitivity[j])
    } else {
      metrics$f1_score[j] <- 0
    }
  }
  
  # Find threshold with best accuracy
  best_acc_idx <- which.max(metrics$accuracy)
  best_acc_threshold <- metrics$threshold[best_acc_idx]
  
  # Find threshold with best F1 score
  best_f1_idx <- which.max(metrics$f1_score)
  best_f1_threshold <- metrics$threshold[best_f1_idx]
  
  cat(sprintf("Best accuracy threshold: %.4f (accuracy: %.4f)\n", 
              best_acc_threshold, metrics$accuracy[best_acc_idx]))
  cat(sprintf("Best F1 threshold: %.4f (F1: %.4f)\n", 
              best_f1_threshold, metrics$f1_score[best_f1_idx]))
  
  # Store results
  results_list[[area]] <- list(
    area = area,
    auc = auc_value,
    optimal_threshold = optimal_threshold,
    best_acc_threshold = best_acc_threshold,
    best_f1_threshold = best_f1_threshold,
    metrics = metrics,
    roc = roc_obj,
    n_positive = sum(truth_valid),
    n_total = length(truth_valid)
  )
  
  optimal_cutoffs[[area]] <- optimal_threshold
}

# Create summary table
summary_df <- data.frame(
  area = areas,
  auc = sapply(results_list, function(x) x$auc),
  optimal_threshold = sapply(results_list, function(x) x$optimal_threshold),
  best_acc_threshold = sapply(results_list, function(x) x$best_acc_threshold),
  best_f1_threshold = sapply(results_list, function(x) x$best_f1_threshold),
  n_positive = sapply(results_list, function(x) x$n_positive),
  n_total = sapply(results_list, function(x) x$n_total),
  percent_positive = sapply(results_list, function(x) 100 * x$n_positive / x$n_total)
)

cat("\n\n=== SUMMARY TABLE ===\n")
print(summary_df)

# Save summary table
write.csv(summary_df, paste0(atlaspath, "cellex_benchmark_summary.csv"), row.names = FALSE)

# Create plots for each area
pdf(paste0(atlaspath, "cellex_benchmark_plots.pdf"), width = 14, height = 10)

for (area in areas) {
  res <- results_list[[area]]
  
  # Create 2x2 layout
  par(mfrow = c(2, 2))
  
  # 1. ROC curve
  plot(res$roc, main = sprintf("%s - ROC Curve (AUC = %.3f)", area, res$auc),
       col = "blue", lwd = 2)
  
  # 2. Metrics vs threshold
  plot(res$metrics$threshold, res$metrics$accuracy, type = "l", col = "black", lwd = 2,
       xlab = "Threshold", ylab = "Metric Value", 
       main = sprintf("%s - Performance Metrics", area),
       ylim = c(0, 1))
  lines(res$metrics$threshold, res$metrics$sensitivity, col = "red", lwd = 2)
  lines(res$metrics$threshold, res$metrics$specificity, col = "blue", lwd = 2)
  lines(res$metrics$threshold, res$metrics$f1_score, col = "green", lwd = 2)
  abline(v = res$optimal_threshold, lty = 2, col = "purple", lwd = 2)
  legend("topright", legend = c("Accuracy", "Sensitivity", "Specificity", "F1", "Optimal"),
         col = c("black", "red", "blue", "green", "purple"), lty = c(1, 1, 1, 1, 2), lwd = 2)
  
  # 3. Distribution of scores by ground truth
  truth <- create_binary_truth(ground_truth, area)
  valid_cells <- ground_truth != "negative"
  scores <- seurat_merged@meta.data[[paste0("cellex_spatial_", area)]][valid_cells]
  truth_valid <- truth[valid_cells]
  
  boxplot(scores ~ truth_valid, 
          names = c("Negative", "Positive"),
          xlab = "Ground Truth", ylab = "CELLEX Score",
          main = sprintf("%s - Score Distribution", area),
          col = c("lightgray", "lightblue"))
  abline(h = res$optimal_threshold, col = "red", lty = 2, lwd = 2)
  
  # 4. Confusion matrix at optimal threshold
  predicted <- scores >= res$optimal_threshold
  cm <- table(Predicted = predicted, Truth = truth_valid)
  
  # Plot confusion matrix
  image(1:2, 1:2, t(cm)[2:1, ], col = c("white", "lightblue", "lightcoral", "darkblue"),
        xlab = "Truth", ylab = "Predicted", main = sprintf("%s - Confusion Matrix", area),
        axes = FALSE)
  axis(1, at = 1:2, labels = c("Negative", "Positive"))
  axis(2, at = 1:2, labels = c("Positive", "Negative"))
  text(expand.grid(1:2, 1:2), labels = c(cm[2, 1], cm[2, 2], cm[1, 1], cm[1, 2]), cex = 2)
}

dev.off()

cat("\n\n=== Files saved ===\n")
cat("1. cellex_benchmark_summary.csv - Summary table with AUC and optimal thresholds\n")
cat("2. cellex_benchmark_plots.pdf - Diagnostic plots for each area\n")

# Store optimal cutoffs in a named vector for easy access
optimal_cutoffs_vec <- unlist(optimal_cutoffs)
cat("\n=== Optimal Cutoffs (for easy reference) ===\n")
print(optimal_cutoffs_vec)

# Save optimal cutoffs
saveRDS(optimal_cutoffs_vec, paste0(atlaspath, "optimal_cellex_cutoffs.rds"))
cat("\n3. optimal_cellex_cutoffs.rds - Optimal cutoffs saved as R object\n")

cat("\n✓ Benchmark complete!\n")

# Spatial statistics approach for projection prediction
# Uses: 1) Extracellular AAV molecules, 2) Glia exclusion, 3) Run-specific area mapping

library(dplyr)
library(ggplot2)
library(FNN)  # for k-nearest neighbors

cat("=== Spatial Statistics Projection Prediction ===\n\n")

# Define run-specific area mappings
run_area_map <- list(
  spatial1 = c("p36_NLS_Cre" = "PAG", "p147_mCherry_2A_iCre" = "BST"),
  spatial2 = c("p36_NLS_Cre" = "PBN", "p147_mCherry_2A_iCre" = "LHA"),
  spatial3 = c("p36_NLS_Cre" = "CEA", "p147_mCherry_2A_iCre" = "PVT")
)

# Get all FOVs and their run indices
fov_run_map <- unique(seurat_merged@meta.data[, c("old.ident", "run_index")])
rownames(fov_run_map) <- NULL
colnames(fov_run_map) <- c("fov", "run_index")

cat("FOV to run mapping:\n")
print(fov_run_map)

# Initialize prediction column
seurat_merged$projection_spatial_pred <- "negative"

# Process each FOV
for (fov_idx in 1:nrow(fov_run_map)) {
  fov <- fov_run_map$fov[fov_idx]
  run <- fov_run_map$run_index[fov_idx]
  
  cat(sprintf("\n--- Processing %s (run: %s) [%d/%d] ---\n", fov, run, fov_idx, nrow(fov_run_map)))
  
  # Get the two areas for this run
  area_map <- run_area_map[[run]]
  gene1 <- names(area_map)[1]  # p36-NLS-Cre
  gene2 <- names(area_map)[2]  # p147-mCherry-2A-iCre
  area1 <- area_map[1]  # e.g., PAG
  area2 <- area_map[2]  # e.g., BST
  
  cat(sprintf("  %s -> %s\n", gene1, area1))
  cat(sprintf("  %s -> %s\n", gene2, area2))
  
  # Get molecule coordinates for this FOV
  cat("  Getting molecule coordinates...\n")
  flush.console()
  
  tryCatch({
    molecules <- GetTissueCoordinates(seurat_merged, image = fov, which = "molecules")
    cat(sprintf("  ✓ Loaded %d total molecules\n", nrow(molecules)))
  }, error = function(e) {
    cat(sprintf("  ✗ Error loading molecules: %s\n", e$message))
    cat("  Skipping this FOV...\n")
    return(NULL)
  })
  
  if (!exists("molecules") || is.null(molecules)) {
    next
  }
  
  # Check column name for gene/feature
  gene_col <- if("molecule" %in% colnames(molecules)) {
    "molecule"
  } else if("gene" %in% colnames(molecules)) {
    "gene"
  } else {
    "feature"
  }
  
  # DIAGNOSTIC: Show what genes/features are actually present (first FOV only)
  if (fov_idx == 1) {
    cat("\n  DIAGNOSTIC INFO (first FOV only):\n")
    cat(sprintf("  Molecule columns: %s\n", paste(colnames(molecules), collapse = ", ")))
    cat(sprintf("  Molecule data structure:\n"))
    print(str(molecules))
    cat(sprintf("  First few rows of molecule data:\n"))
    print(head(molecules, 10))
    
    cat(sprintf("  Gene/Feature column: %s\n", gene_col))
    
    # Check if it's actually called "molecule" instead
    if ("molecule" %in% colnames(molecules)) {
      cat("  Using 'molecule' column instead\n")
      gene_col <- "molecule"
      unique_genes <- unique(molecules[[gene_col]])
    } else {
      unique_genes <- unique(molecules[[gene_col]])
    }
    
    cat("  Unique genes/features (first 20):\n")
    print(head(unique_genes, 20))
    cat(sprintf("  Total unique genes: %d\n", length(unique_genes)))
    
    # Check for AAV-like names
    if (length(unique_genes) > 0) {
      aav_genes <- grep("p36|p147|AAV|Cre|mCherry|iCre|NLS|2A", unique_genes, value = TRUE, ignore.case = TRUE)
      if (length(aav_genes) > 0) {
        cat("  Possible AAV-related genes found:\n")
        print(aav_genes)
      }
    }
    cat("\n")
  }
  
  # Filter for the two AAV genes
  mol_gene1 <- molecules[molecules[[gene_col]] == gene1, ]
  mol_gene2 <- molecules[molecules[[gene_col]] == gene2, ]
  
  cat(sprintf("  Molecules - %s: %d, %s: %d\n", 
              gene1, nrow(mol_gene1), gene2, nrow(mol_gene2)))
  
  # Get cell coordinates for this FOV
  cells <- GetTissueCoordinates(seurat_merged, image = fov, which = "centroids")
  
  # Get cell metadata - filter from full Seurat object for this FOV
  fov_cell_mask <- seurat_merged$old.ident == fov
  cell_meta <- seurat_merged@meta.data[fov_cell_mask, ]
  cell_ids_in_fov <- rownames(seurat_merged@meta.data)[fov_cell_mask]
  
  # The centroids have numeric rownames (1, 2, 3...) but we need full cell IDs
  # The order should match, so we can assign the full IDs to the centroids
  if (nrow(cells) == length(cell_ids_in_fov)) {
    rownames(cells) <- cell_ids_in_fov
    if (fov_idx == 1) {
      cat("  ✓ Successfully mapped cell coordinates to cell IDs\n")
      cat(sprintf("  Remapped centroid rownames: %s\n", 
                  paste(head(rownames(cells), 3), collapse = ", ")))
    }
  } else {
    cat(sprintf("  Warning: Mismatch in cell counts - centroids: %d, metadata: %d\n",
                nrow(cells), length(cell_ids_in_fov)))
    next
  }
  
  # DIAGNOSTIC: Check neurotransmitter and injection_filter (first FOV only)
  if (fov_idx == 1) {
    cat("\n  DIAGNOSTIC - Cell metadata:\n")
    cat(sprintf("  Total cells in FOV: %d\n", nrow(cell_meta)))
    cat(sprintf("  neurotransmitter column exists: %s\n", "neurotransmitter" %in% colnames(cell_meta)))
    cat(sprintf("  injection_filter column exists: %s\n", "injection_filter" %in% colnames(cell_meta)))
    
    if ("neurotransmitter" %in% colnames(cell_meta)) {
      cat("  neurotransmitter values in this FOV:\n")
      print(table(cell_meta$neurotransmitter, useNA = "always"))
    }
    
    if ("injection_filter" %in% colnames(cell_meta)) {
      cat("  injection_filter values in this FOV:\n")
      print(table(cell_meta$injection_filter, useNA = "always"))
    }
    cat("\n")
  }
  
  # Filter out glia (they don't project) - handle NAs
  is_neuron <- !is.na(cell_meta$neurotransmitter) & cell_meta$neurotransmitter != "Glia"
  
  # Also check injection filter status - handle NAs
  is_projection_site <- !is.na(cell_meta$injection_filter) & cell_meta$injection_filter == "projection"
  
  cat(sprintf("  Cells - Total: %d, Neurons: %d, Glia/NA: %d\n", 
              nrow(cell_meta), sum(is_neuron), sum(!is_neuron)))
  cat(sprintf("  Neurons - Projection sites: %d, Injection sites/NA: %d\n",
              sum(is_neuron & is_projection_site), 
              sum(is_neuron & !is_projection_site)))
  
  # Only predict for neurons at projection sites (not injection sites)
  is_valid <- is_neuron & is_projection_site
  
  if (sum(is_valid, na.rm = TRUE) == 0) {
    cat("  No valid neurons (projection sites) in this FOV, skipping...\n")
    next
  }
  
  # Get valid cell IDs
  valid_cell_ids <- cell_ids_in_fov[is_valid]
  
  # Get coordinates for valid cells - match by rownames
  # The centroids might have all cells, we need to subset
  neuron_cells <- cells[rownames(cells) %in% valid_cell_ids, ]
  
  if (nrow(neuron_cells) == 0) {
    cat("  Warning: Valid neurons found in metadata but no coordinates available\n")
    cat("  This suggests cell ID mismatch between metadata and spatial coordinates\n")
    next
  }
  
  neuron_ids <- rownames(neuron_cells)
  
  # Calculate spatial enrichment for each neuron
  # Method: Count molecules within a radius around each cell
  # Use both inside-cell and outside-cell molecules
  
  # DIAGNOSTIC: Check coordinate ranges (first FOV only)
  if (fov_idx == 1) {
    cat("\n  DIAGNOSTIC - Coordinate ranges:\n")
    cat(sprintf("  Neuron X range: %.0f to %.0f\n", min(neuron_cells$x), max(neuron_cells$x)))
    cat(sprintf("  Neuron Y range: %.0f to %.0f\n", min(neuron_cells$y), max(neuron_cells$y)))
    if (nrow(mol_gene2) > 0) {
      cat(sprintf("  BST molecule X range: %.0f to %.0f\n", min(mol_gene2$x), max(mol_gene2$x)))
      cat(sprintf("  BST molecule Y range: %.0f to %.0f\n", min(mol_gene2$y), max(mol_gene2$y)))
      
      # Calculate sample distances
      sample_neuron <- neuron_cells[1, c("x", "y")]
      sample_dists <- sqrt((sample_neuron$x - mol_gene2$x[1:min(10, nrow(mol_gene2))])^2 + 
                          (sample_neuron$y - mol_gene2$y[1:min(10, nrow(mol_gene2))])^2)
      cat(sprintf("  Sample distances from first neuron to first 10 BST molecules:\n"))
      cat(sprintf("    Min: %.1f, Max: %.1f, Mean: %.1f\n", 
                  min(sample_dists), max(sample_dists), mean(sample_dists)))
    }
    cat("\n")
  }
  
  search_radius <- 40  # pixels
  
  # For each neuron, count nearby molecules of each AAV type
  if (nrow(mol_gene1) > 0 && nrow(neuron_cells) > 0) {
    # Calculate distances from each neuron to all gene1 molecules
    neuron_coords <- as.matrix(neuron_cells[, c("x", "y")])
    mol1_coords <- as.matrix(mol_gene1[, c("x", "y")])
    
    # Calculate pairwise distances
    count_gene1 <- sapply(1:nrow(neuron_coords), function(i) {
      dists <- sqrt((neuron_coords[i, 1] - mol1_coords[, 1])^2 + 
                    (neuron_coords[i, 2] - mol1_coords[, 2])^2)
      sum(dists <= search_radius)
    })
  } else {
    count_gene1 <- rep(0, nrow(neuron_cells))
  }
  
  if (nrow(mol_gene2) > 0 && nrow(neuron_cells) > 0) {
    neuron_coords <- as.matrix(neuron_cells[, c("x", "y")])
    mol2_coords <- as.matrix(mol_gene2[, c("x", "y")])
    
    count_gene2 <- sapply(1:nrow(neuron_coords), function(i) {
      dists <- sqrt((neuron_coords[i, 1] - mol2_coords[, 1])^2 + 
                    (neuron_coords[i, 2] - mol2_coords[, 2])^2)
      sum(dists <= search_radius)
    })
  } else {
    count_gene2 <- rep(0, nrow(neuron_cells))
  }
  
  # DIAGNOSTIC: Show distribution of molecule counts (first FOV only)
  if (fov_idx == 1) {
    cat(sprintf("  Molecule counts per neuron - %s: min=%d, max=%d, mean=%.1f\n",
                gene1, min(count_gene1), max(count_gene1), mean(count_gene1)))
    cat(sprintf("  Molecule counts per neuron - %s: min=%d, max=%d, mean=%.1f\n",
                gene2, min(count_gene2), max(count_gene2), mean(count_gene2)))
  }
  
  # Determine projection based on molecule counts
  # Use a threshold approach
  threshold <- 1  # lowered from 3 to be more permissive
  
  has_proj1 <- count_gene1 >= threshold
  has_proj2 <- count_gene2 >= threshold
  
  # Assign projections
  predictions <- rep("negative", length(neuron_ids))
  predictions[has_proj1 & !has_proj2] <- area1
  predictions[!has_proj1 & has_proj2] <- area2
  predictions[has_proj1 & has_proj2] <- paste0(area1, "+", area2)
  
  # Update seurat object
  seurat_merged$projection_spatial_pred[neuron_ids] <- predictions
  
  cat(sprintf("  Predictions - %s: %d, %s: %d, %s+%s: %d, negative: %d\n",
              area1, sum(predictions == area1),
              area2, sum(predictions == area2),
              area1, area2, sum(predictions == paste0(area1, "+", area2)),
              sum(predictions == "negative")))
}

# Compare predictions to ground truth
cat("\n\n=== Comparison to Ground Truth ===\n")

# Create comparison dataframe
# Exclude: glia, negative predictions, and injection sites (not long-range projections)
# Handle NAs properly
valid_cells <- !is.na(seurat_merged$neurotransmitter) &
               seurat_merged$neurotransmitter != "Glia" & 
               !is.na(seurat_merged$projection) &
               seurat_merged$projection != "negative" &
               !is.na(seurat_merged$projection_spatial_pred) &
               seurat_merged$projection_spatial_pred != "negative" &
               !is.na(seurat_merged$injection_filter) &
               seurat_merged$injection_filter == "projection"  # Exclude injection sites

comparison_df <- data.frame(
  ground_truth = seurat_merged$projection[valid_cells],
  predicted = seurat_merged$projection_spatial_pred[valid_cells],
  run_index = seurat_merged$run_index[valid_cells],
  stringsAsFactors = FALSE
)

cat(sprintf("Cells for comparison (projection sites only, excluding injections): %d\n", nrow(comparison_df)))
cat(sprintf("Cells excluded (injection sites): %d\n", 
            sum(!is.na(seurat_merged$injection_filter) & 
                seurat_merged$injection_filter == "injection" & 
                !is.na(seurat_merged$neurotransmitter) &
                seurat_merged$neurotransmitter != "Glia", na.rm = TRUE)))
cat(sprintf("Cells with NA metadata: %d\n",
            sum(is.na(seurat_merged$neurotransmitter) | 
                is.na(seurat_merged$injection_filter))))

# Overall accuracy
accuracy <- mean(comparison_df$ground_truth == comparison_df$predicted)
cat(sprintf("Overall accuracy: %.2f%%\n", 100 * accuracy))

# Confusion matrix
cat("\nConfusion Matrix:\n")
conf_matrix <- table(Predicted = comparison_df$predicted, 
                     Truth = comparison_df$ground_truth)
print(conf_matrix)

# Calculate metrics per area
all_areas <- unique(c(comparison_df$ground_truth, comparison_df$predicted))
all_areas <- all_areas[all_areas != "negative"]

metrics_per_area <- data.frame(
  area = character(),
  precision = numeric(),
  recall = numeric(),
  f1_score = numeric(),
  n_true = integer(),
  n_predicted = integer(),
  stringsAsFactors = FALSE
)

for (area in all_areas) {
  tp <- sum(comparison_df$ground_truth == area & comparison_df$predicted == area)
  fp <- sum(comparison_df$ground_truth != area & comparison_df$predicted == area)
  fn <- sum(comparison_df$ground_truth == area & comparison_df$predicted != area)
  
  precision <- if (tp + fp > 0) tp / (tp + fp) else 0
  recall <- if (tp + fn > 0) tp / (tp + fn) else 0
  f1 <- if (precision + recall > 0) 2 * (precision * recall) / (precision + recall) else 0
  
  n_true <- sum(comparison_df$ground_truth == area)
  n_pred <- sum(comparison_df$predicted == area)
  
  metrics_per_area <- rbind(metrics_per_area, data.frame(
    area = area,
    precision = precision,
    recall = recall,
    f1_score = f1,
    n_true = n_true,
    n_predicted = n_pred
  ))
}

cat("\nMetrics per area:\n")
print(metrics_per_area)

# Accuracy by run
cat("\nAccuracy by run:\n")
accuracy_by_run <- comparison_df %>%
  group_by(run_index) %>%
  summarize(
    n_cells = n(),
    accuracy = mean(ground_truth == predicted),
    .groups = "drop"
  )
print(accuracy_by_run)

# Visualizations
pdf(paste0(atlaspath, "spatial_projection_prediction_analysis.pdf"), width = 14, height = 10)

# 1. Confusion matrix heatmap
library(reshape2)
conf_matrix_melted <- melt(conf_matrix)
colnames(conf_matrix_melted) <- c("Predicted", "Truth", "Count")

p1 <- ggplot(conf_matrix_melted, aes(x = Truth, y = Predicted, fill = Count)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Count), size = 4) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Confusion Matrix: Spatial Prediction vs Ground Truth",
       x = "Ground Truth", y = "Predicted")
print(p1)

# 2. Metrics per area
p2 <- metrics_per_area %>%
  pivot_longer(cols = c(precision, recall, f1_score), 
               names_to = "metric", values_to = "value") %>%
  ggplot(aes(x = area, y = value, fill = metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Performance Metrics by Projection Area",
       x = "Area", y = "Score", fill = "Metric") +
  scale_fill_brewer(palette = "Set2")
print(p2)

# 3. Accuracy by run
p3 <- ggplot(accuracy_by_run, aes(x = run_index, y = accuracy, fill = run_index)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = sprintf("%.1f%%", 100*accuracy)), vjust = -0.5) +
  theme_minimal() +
  ylim(0, 1) +
  labs(title = "Prediction Accuracy by Run",
       x = "Run Index", y = "Accuracy") +
  theme(legend.position = "none")
print(p3)

dev.off()

# Save results
write.csv(comparison_df, paste0(atlaspath, "spatial_projection_predictions.csv"), row.names = FALSE)
write.csv(metrics_per_area, paste0(atlaspath, "spatial_prediction_metrics.csv"), row.names = FALSE)

cat("\n\n=== Files saved ===\n")
cat("1. spatial_projection_predictions.csv - Cell-level predictions vs ground truth\n")
cat("2. spatial_prediction_metrics.csv - Performance metrics per area\n")
cat("3. spatial_projection_prediction_analysis.pdf - Visualization plots\n")

cat("\n✓ Spatial prediction complete!\n")
cat("New column added: seurat_merged$projection_spatial_pred\n")

# Optimize Spatial Statistics Parameters for Best Agreement with Spatial CELLEX
# Iterate: search_radius, threshold, and potentially other parameters
# Goal: Find parameters that maximize correlation between:
#   - Spatial stats aggregated by cell type
#   - Spatial CELLEX scores (cell type x area)

library(dplyr)
library(ggplot2)
library(tidyr)

cat("=== Optimizing Spatial Statistics Parameters ===\n\n")

# Read spatial CELLEX for comparison
cellex_spatial <- read.csv("combined_spatial_cellex.csv", 
                           row.names = 1, check.names = FALSE)

areas_common <- c("PAG", "BST", "PBN", "LHA", "CEA", "PVT")
cellex_spatial <- cellex_spatial[areas_common, ]

cat("Reference: Spatial CELLEX matrix\n")
cat(sprintf("Dimensions: %d areas × %d cell types\n\n", 
            nrow(cellex_spatial), ncol(cellex_spatial)))

# Define run-specific area mappings
run_area_map <- list(
  spatial1 = c("p36_NLS_Cre" = "PAG", "p147_mCherry_2A_iCre" = "BST"),
  spatial2 = c("p36_NLS_Cre" = "PBN", "p147_mCherry_2A_iCre" = "LHA"),
  spatial3 = c("p36_NLS_Cre" = "CEA", "p147_mCherry_2A_iCre" = "PVT")
)

# Get FOV mapping
fov_run_map <- unique(seurat_merged@meta.data[, c("old.ident", "run_index")])
rownames(fov_run_map) <- NULL
colnames(fov_run_map) <- c("fov", "run_index")

cat(sprintf("Processing %d FOVs\n\n", nrow(fov_run_map)))

# PARAMETER GRID
search_radius_values <- c(20, 30, 40, 50, 60, 75, 100)  # pixels
threshold_values <- c(1, 2, 3, 5, 10, 15, 20)  # minimum molecules

cat("Testing parameter combinations:\n")
cat(sprintf("  Search radius: %s\n", paste(search_radius_values, collapse = ", ")))
cat(sprintf("  Thresholds: %s\n\n", paste(threshold_values, collapse = ", ")))

# Store results
optimization_results <- expand.grid(
  search_radius = search_radius_values,
  threshold = threshold_values,
  correlation_with_cellex = NA,
  mean_max_score = NA,
  prop_predictions = NA,
  n_predictions = NA,
  stringsAsFactors = FALSE
)

total_combinations <- nrow(optimization_results)
cat(sprintf("Total combinations to test: %d\n\n", total_combinations))

# Run optimization
for (param_idx in 1:nrow(optimization_results)) {
  search_radius <- optimization_results$search_radius[param_idx]
  threshold <- optimization_results$threshold[param_idx]
  
  cat(sprintf("[%d/%d] Testing radius=%d, threshold=%d... ", 
              param_idx, total_combinations, search_radius, threshold))
  flush.console()
  
  # Initialize temporary prediction column
  temp_predictions <- rep("negative", nrow(seurat_merged@meta.data))
  names(temp_predictions) <- rownames(seurat_merged@meta.data)
  
  # Process each FOV with current parameters
  for (fov_idx in 1:nrow(fov_run_map)) {
    fov <- fov_run_map$fov[fov_idx]
    run <- fov_run_map$run_index[fov_idx]
    
    area_map <- run_area_map[[run]]
    gene1 <- names(area_map)[1]
    gene2 <- names(area_map)[2]
    area1 <- area_map[1]
    area2 <- area_map[2]
    
    # Get molecules and cells
    tryCatch({
      molecules <- GetTissueCoordinates(seurat_merged, image = fov, which = "molecules")
      cells <- GetTissueCoordinates(seurat_merged, image = fov, which = "centroids")
      
      gene_col <- if("molecule" %in% colnames(molecules)) "molecule" else "feature"
      
      mol_gene1 <- molecules[molecules[[gene_col]] == gene1, ]
      mol_gene2 <- molecules[molecules[[gene_col]] == gene2, ]
      
      # Get cell metadata
      fov_cell_mask <- seurat_merged$old.ident == fov
      cell_meta <- seurat_merged@meta.data[fov_cell_mask, ]
      cell_ids_in_fov <- rownames(seurat_merged@meta.data)[fov_cell_mask]
      
      # Map cell IDs to centroids
      if (nrow(cells) == length(cell_ids_in_fov)) {
        rownames(cells) <- cell_ids_in_fov
      } else {
        next
      }
      
      # Filter for valid neurons at projection sites
      is_neuron <- !is.na(cell_meta$neurotransmitter) & 
                   cell_meta$neurotransmitter != "Glia"
      is_projection_site <- !is.na(cell_meta$injection_filter) & 
                           cell_meta$injection_filter == "projection"
      is_valid <- is_neuron & is_projection_site
      
      if (sum(is_valid, na.rm = TRUE) == 0) next
      
      valid_cell_ids <- cell_ids_in_fov[is_valid]
      neuron_cells <- cells[rownames(cells) %in% valid_cell_ids, ]
      
      if (nrow(neuron_cells) == 0) next
      
      neuron_ids <- rownames(neuron_cells)
      
      # Count molecules within radius
      if (nrow(mol_gene1) > 0 && nrow(neuron_cells) > 0) {
        neuron_coords <- as.matrix(neuron_cells[, c("x", "y")])
        mol1_coords <- as.matrix(mol_gene1[, c("x", "y")])
        
        count_gene1 <- sapply(1:nrow(neuron_coords), function(i) {
          dists <- sqrt((neuron_coords[i, 1] - mol1_coords[, 1])^2 + 
                       (neuron_coords[i, 2] - mol1_coords[, 2])^2)
          sum(dists <= search_radius)
        })
      } else {
        count_gene1 <- rep(0, nrow(neuron_cells))
      }
      
      if (nrow(mol_gene2) > 0 && nrow(neuron_cells) > 0) {
        neuron_coords <- as.matrix(neuron_cells[, c("x", "y")])
        mol2_coords <- as.matrix(mol_gene2[, c("x", "y")])
        
        count_gene2 <- sapply(1:nrow(neuron_coords), function(i) {
          dists <- sqrt((neuron_coords[i, 1] - mol2_coords[, 1])^2 + 
                       (neuron_coords[i, 2] - mol2_coords[, 2])^2)
          sum(dists <= search_radius)
        })
      } else {
        count_gene2 <- rep(0, nrow(neuron_cells))
      }
      
      # Apply threshold and assign predictions
      has_proj1 <- count_gene1 >= threshold
      has_proj2 <- count_gene2 >= threshold
      
      predictions <- rep("negative", length(neuron_ids))
      predictions[has_proj1 & !has_proj2] <- area1
      predictions[!has_proj1 & has_proj2] <- area2
      predictions[has_proj1 & has_proj2] <- paste0(area1, "+", area2)
      
      temp_predictions[neuron_ids] <- predictions
      
    }, error = function(e) {
      # Skip FOV on error
    })
  }
  
  # Aggregate predictions by cell type
  prediction_data <- data.frame(
    cell_type = seurat_merged$plot,
    prediction = temp_predictions,
    neurotransmitter = seurat_merged$neurotransmitter,
    injection_filter = seurat_merged$injection_filter,
    stringsAsFactors = FALSE
  ) %>%
    filter(!is.na(neurotransmitter) & neurotransmitter != "Glia",
           !is.na(injection_filter) & injection_filter == "projection",
           prediction != "negative")
  
  if (nrow(prediction_data) == 0) {
    cat("No predictions\n")
    next
  }
  
  # Calculate proportion projecting to each area by cell type
  spatial_stats_summary <- prediction_data %>%
    group_by(cell_type) %>%
    summarize(
      n_cells = n(),
      prop_PAG = sum(grepl("PAG", prediction)) / n(),
      prop_BST = sum(grepl("BST", prediction)) / n(),
      prop_PBN = sum(grepl("PBN", prediction)) / n(),
      prop_LHA = sum(grepl("LHA", prediction)) / n(),
      prop_CEA = sum(grepl("CEA", prediction)) / n(),
      prop_PVT = sum(grepl("PVT", prediction)) / n(),
      .groups = "drop"
    )
  
  # Find common cell types
  common_cell_types <- intersect(colnames(cellex_spatial), 
                                 spatial_stats_summary$cell_type)
  
  if (length(common_cell_types) < 10) {
    cat("Too few cell types\n")
    next
  }
  
  # Create spatial stats matrix
  spatial_stats_mat <- matrix(NA, nrow = length(areas_common), 
                              ncol = length(common_cell_types),
                              dimnames = list(areas_common, common_cell_types))
  
  for (area in areas_common) {
    col_name <- paste0("prop_", area)
    for (ct in common_cell_types) {
      val <- spatial_stats_summary %>% 
        filter(cell_type == ct) %>% 
        pull(!!col_name)
      spatial_stats_mat[area, ct] <- ifelse(length(val) > 0, val, 0)
    }
  }
  
  # Get CELLEX matrix for same cell types
  cellex_mat <- as.matrix(cellex_spatial[areas_common, common_cell_types])
  cellex_mat[is.na(cellex_mat)] <- 0
  
  # Calculate correlation
  correlations <- numeric(length(areas_common))
  for (i in seq_along(areas_common)) {
    area <- areas_common[i]
    has_data <- !is.na(cellex_mat[area, ]) & !is.na(spatial_stats_mat[area, ])
    
    if (sum(has_data) > 3) {
      correlations[i] <- cor(cellex_mat[area, has_data], 
                            spatial_stats_mat[area, has_data], 
                            method = "spearman", use = "complete.obs")
    } else {
      correlations[i] <- NA
    }
  }
  
  mean_cor <- mean(correlations, na.rm = TRUE)
  
  # Calculate other metrics
  max_scores <- apply(spatial_stats_mat, 2, max, na.rm = TRUE)
  mean_max <- mean(max_scores[is.finite(max_scores)])
  
  n_pred <- nrow(prediction_data)
  prop_pred <- n_pred / sum(!is.na(seurat_merged$neurotransmitter) & 
                            seurat_merged$neurotransmitter != "Glia" &
                            !is.na(seurat_merged$injection_filter) &
                            seurat_merged$injection_filter == "projection")
  
  # Store results
  optimization_results$correlation_with_cellex[param_idx] <- mean_cor
  optimization_results$mean_max_score[param_idx] <- mean_max
  optimization_results$prop_predictions[param_idx] <- prop_pred
  optimization_results$n_predictions[param_idx] <- n_pred
  
  cat(sprintf("cor=%.3f, n=%d\n", mean_cor, n_pred))
}

cat("\n=== Optimization Complete ===\n\n")

# Find best parameters
best_idx <- which.max(optimization_results$correlation_with_cellex)
best_params <- optimization_results[best_idx, ]

cat("=== OPTIMAL PARAMETERS ===\n")
cat(sprintf("Search radius: %d pixels\n", best_params$search_radius))
cat(sprintf("Threshold: %d molecules\n", best_params$threshold))
cat(sprintf("Correlation with spatial CELLEX: %.3f\n", best_params$correlation_with_cellex))
cat(sprintf("Mean max score: %.3f\n", best_params$mean_max_score))
cat(sprintf("Proportion with predictions: %.2f%%\n", 100 * best_params$prop_predictions))
cat(sprintf("Number of predictions: %d\n", best_params$n_predictions))

cat("\n=== Top 10 Parameter Combinations ===\n")
top_10 <- optimization_results %>%
  arrange(desc(correlation_with_cellex)) %>%
  head(10) %>%
  mutate(rank = 1:10) %>%
  select(rank, search_radius, threshold, correlation_with_cellex, 
         mean_max_score, prop_predictions, n_predictions)
print(top_10)

# VISUALIZATIONS
pdf(paste0(atlaspath, "spatial_stats_optimization.pdf"), width = 14, height = 12)

# 1. Heatmap of correlation vs parameters
cor_matrix <- optimization_results %>%
  select(search_radius, threshold, correlation_with_cellex) %>%
  pivot_wider(names_from = threshold, values_from = correlation_with_cellex) %>%
  select(-search_radius) %>%
  as.matrix()

rownames(cor_matrix) <- sprintf("r=%d", search_radius_values)

image(1:ncol(cor_matrix), 1:nrow(cor_matrix), t(cor_matrix),
      col = colorRampPalette(c("blue", "white", "red"))(100),
      xlab = "Threshold (molecules)", ylab = "Search Radius (pixels)",
      main = "Correlation with Spatial CELLEX\n(Red = Better)",
      axes = FALSE)
axis(1, at = 1:ncol(cor_matrix), labels = colnames(cor_matrix))
axis(2, at = 1:nrow(cor_matrix), labels = rownames(cor_matrix), las = 1)
contour(1:ncol(cor_matrix), 1:nrow(cor_matrix), t(cor_matrix), add = TRUE)

# Mark optimum
opt_thresh_idx <- which(threshold_values == best_params$threshold)
opt_radius_idx <- which(search_radius_values == best_params$search_radius)
points(opt_thresh_idx, opt_radius_idx, pch = 8, cex = 3, col = "red", lwd = 3)

# 2. Effect of search radius
par(mfrow = c(2, 2))

for (thresh in c(1, 3, 5, 10)) {
  thresh_data <- optimization_results %>%
    filter(threshold == thresh)
  
  plot(thresh_data$search_radius, thresh_data$correlation_with_cellex,
       type = "b", pch = 19, col = "blue",
       xlab = "Search Radius (pixels)", ylab = "Correlation",
       main = sprintf("Effect of Radius (threshold=%d)", thresh),
       ylim = range(optimization_results$correlation_with_cellex, na.rm = TRUE))
  grid()
  
  if (thresh == best_params$threshold) {
    abline(v = best_params$search_radius, col = "red", lty = 2, lwd = 2)
  }
}

# 3. Effect of threshold
par(mfrow = c(2, 2))

for (rad in c(30, 50, 75, 100)) {
  rad_data <- optimization_results %>%
    filter(search_radius == rad)
  
  plot(rad_data$threshold, rad_data$correlation_with_cellex,
       type = "b", pch = 19, col = "darkgreen",
       xlab = "Threshold (molecules)", ylab = "Correlation",
       main = sprintf("Effect of Threshold (radius=%d)", rad),
       ylim = range(optimization_results$correlation_with_cellex, na.rm = TRUE))
  grid()
  
  if (rad == best_params$search_radius) {
    abline(v = best_params$threshold, col = "red", lty = 2, lwd = 2)
  }
}

# 4. Correlation vs number of predictions
par(mfrow = c(1, 1))
plot(optimization_results$n_predictions, optimization_results$correlation_with_cellex,
     pch = 19, col = rgb(0, 0, 1, 0.5),
     xlab = "Number of Predictions", ylab = "Correlation with CELLEX",
     main = "Trade-off: Prediction Count vs Quality")
points(best_params$n_predictions, best_params$correlation_with_cellex, 
       pch = 8, cex = 2, col = "red", lwd = 2)
text(best_params$n_predictions, best_params$correlation_with_cellex,
     sprintf("r=%d, t=%d", best_params$search_radius, best_params$threshold),
     pos = 4, col = "red")
grid()

dev.off()

# Save results
write.csv(optimization_results, 
          paste0(atlaspath, "spatial_stats_parameter_optimization.csv"), 
          row.names = FALSE)

optimal_params_summary <- data.frame(
  parameter = c("search_radius", "threshold", "correlation_with_cellex",
                "mean_max_score", "prop_predictions", "n_predictions"),
  value = c(best_params$search_radius, best_params$threshold, 
            best_params$correlation_with_cellex, best_params$mean_max_score,
            best_params$prop_predictions, best_params$n_predictions)
)
write.csv(optimal_params_summary, 
          paste0(atlaspath, "optimal_spatial_stats_parameters.csv"), 
          row.names = FALSE)

cat("\n=== Files Saved ===\n")
cat("1. spatial_stats_optimization.pdf - Optimization visualizations\n")
cat("2. spatial_stats_parameter_optimization.csv - All tested combinations\n")
cat("3. optimal_spatial_stats_parameters.csv - Best parameters summary\n")

cat("\n✓ Optimization complete!\n")
cat(sprintf("\nRecommendation: Use radius=%d pixels, threshold=%d molecules\n",
            best_params$search_radius, best_params$threshold))

# Optimize CELLEX Prediction Strategy for Agreement with Spatial Stats
# Goal: Find best way to convert CELLEX ESμ scores into area predictions
# Compare different strategies against seurat_merged$projection_spatial_pred

library(dplyr)
library(ggplot2)
library(tidyr)

cat("=== Optimizing CELLEX Prediction Strategy ===\n\n")

# Read spatial CELLEX
cellex_spatial <- read.csv("combined_spatial_cellex.csv", 
                           row.names = 1, check.names = FALSE)

areas_common <- c("PAG", "BST", "PBN", "LHA", "CEA", "PVT")
cellex_spatial <- cellex_spatial[areas_common, ]

cat("Spatial CELLEX matrix:\n")
cat(sprintf("  %d areas × %d cell types\n\n", nrow(cellex_spatial), ncol(cellex_spatial)))

# Get spatial stats predictions - aggregate to cell type level
spatial_stats_data <- data.frame(
  cell_type = seurat_merged$plot,
  prediction = seurat_merged$projection_spatial_pred,
  neurotransmitter = seurat_merged$neurotransmitter,
  injection_filter = seurat_merged$injection_filter,
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(neurotransmitter) & neurotransmitter != "Glia",
         !is.na(injection_filter) & injection_filter == "projection",
         prediction != "negative")

cat(sprintf("Spatial stats: %d cells with projections\n", nrow(spatial_stats_data)))

# Aggregate to cell type × area proportions
spatial_stats_summary <- spatial_stats_data %>%
  group_by(cell_type) %>%
  summarize(
    n_cells = n(),
    prop_PAG = sum(grepl("PAG", prediction)) / n(),
    prop_BST = sum(grepl("BST", prediction)) / n(),
    prop_PBN = sum(grepl("PBN", prediction)) / n(),
    prop_LHA = sum(grepl("LHA", prediction)) / n(),
    prop_CEA = sum(grepl("CEA", prediction)) / n(),
    prop_PVT = sum(grepl("PVT", prediction)) / n(),
    .groups = "drop"
  )

cat(sprintf("Cell types: %d\n\n", nrow(spatial_stats_summary)))

# STRATEGY 1: Absolute threshold on ESμ
# Test different thresholds
abs_threshold_values <- seq(0, 0.5, by = 0.05)

# STRATEGY 2: Relative threshold (fraction of max ESμ for that cell type)
rel_threshold_values <- seq(0.3, 0.9, by = 0.1)

# STRATEGY 3: Top N areas (regardless of score)
top_n_values <- 1:3

# Store all results
all_results <- data.frame(
  strategy = character(),
  parameter = numeric(),
  correlation = numeric(),
  mean_agreement = numeric(),
  n_predictions = integer(),
  stringsAsFactors = FALSE
)

cat("Testing prediction strategies...\n\n")

# STRATEGY 1: Absolute threshold
cat("Strategy 1: Absolute ESμ threshold\n")
for (thresh in abs_threshold_values) {
  cat(sprintf("  Testing threshold %.2f... ", thresh))
  
  # For each cell type, which areas are above threshold?
  predictions <- data.frame(
    cell_type = colnames(cellex_spatial),
    stringsAsFactors = FALSE
  )
  
  for (i in 1:ncol(cellex_spatial)) {
    ct <- colnames(cellex_spatial)[i]
    scores <- cellex_spatial[, ct]
    scores[is.na(scores)] <- 0
    
    above_thresh <- areas_common[scores >= thresh]
    predictions$prediction[i] <- ifelse(length(above_thresh) == 0, 
                                       "negative",
                                       paste(sort(above_thresh), collapse = "+"))
  }
  
  # Compare to spatial stats at proportion level
  common_cts <- intersect(spatial_stats_summary$cell_type, predictions$cell_type)
  
  # Create binary matrices
  cellex_binary <- matrix(0, nrow = length(areas_common), ncol = length(common_cts),
                         dimnames = list(areas_common, common_cts))
  
  for (ct in common_cts) {
    pred <- predictions$prediction[predictions$cell_type == ct]
    if (pred != "negative") {
      pred_areas <- strsplit(pred, "\\+")[[1]]
      cellex_binary[pred_areas, ct] <- 1
    }
  }
  
  # Spatial stats as continuous proportions
  spatial_cont <- matrix(0, nrow = length(areas_common), ncol = length(common_cts),
                        dimnames = list(areas_common, common_cts))
  
  for (area in areas_common) {
    col_name <- paste0("prop_", area)
    for (ct in common_cts) {
      val <- spatial_stats_summary %>% filter(cell_type == ct) %>% pull(!!col_name)
      spatial_cont[area, ct] <- ifelse(length(val) > 0, val, 0)
    }
  }
  
  # Calculate correlation per area
  cors <- sapply(areas_common, function(area) {
    if (sum(cellex_binary[area, ]) > 0) {
      cor(cellex_binary[area, ], spatial_cont[area, ], method = "spearman")
    } else {
      NA
    }
  })
  
  mean_cor <- mean(cors, na.rm = TRUE)
  
  # Agreement: do they agree on which areas are "top" for each cell type?
  agreements <- sapply(common_cts, function(ct) {
    cellex_areas <- areas_common[cellex_binary[, ct] > 0]
    spatial_areas <- areas_common[spatial_cont[, ct] > 0.15]  # 15% threshold
    
    if (length(cellex_areas) == 0 & length(spatial_areas) == 0) return(TRUE)
    if (length(cellex_areas) == 0 | length(spatial_areas) == 0) return(FALSE)
    
    length(intersect(cellex_areas, spatial_areas)) > 0
  })
  
  mean_agree <- mean(agreements)
  n_pred <- sum(predictions$prediction != "negative")
  
  all_results <- rbind(all_results, data.frame(
    strategy = "absolute_threshold",
    parameter = thresh,
    correlation = mean_cor,
    mean_agreement = mean_agree,
    n_predictions = n_pred
  ))
  
  cat(sprintf("cor=%.3f, agree=%.2f, n=%d\n", mean_cor, mean_agree, n_pred))
}

# STRATEGY 2: Relative threshold (as fraction of max)
cat("\nStrategy 2: Relative threshold (fraction of max ESμ)\n")
for (rel_thresh in rel_threshold_values) {
  cat(sprintf("  Testing %.1f of max... ", rel_thresh))
  
  predictions <- data.frame(
    cell_type = colnames(cellex_spatial),
    stringsAsFactors = FALSE
  )
  
  for (i in 1:ncol(cellex_spatial)) {
    ct <- colnames(cellex_spatial)[i]
    scores <- cellex_spatial[, ct]
    scores[is.na(scores)] <- 0
    
    max_score <- max(scores)
    if (max_score == 0) {
      predictions$prediction[i] <- "negative"
    } else {
      above_thresh <- areas_common[scores >= rel_thresh * max_score]
      predictions$prediction[i] <- ifelse(length(above_thresh) == 0,
                                         "negative",
                                         paste(sort(above_thresh), collapse = "+"))
    }
  }
  
  # Same comparison logic
  common_cts <- intersect(spatial_stats_summary$cell_type, predictions$cell_type)
  
  cellex_binary <- matrix(0, nrow = length(areas_common), ncol = length(common_cts),
                         dimnames = list(areas_common, common_cts))
  
  for (ct in common_cts) {
    pred <- predictions$prediction[predictions$cell_type == ct]
    if (!is.na(pred) && pred != "negative") {
      pred_areas <- strsplit(pred, "\\+")[[1]]
      cellex_binary[pred_areas, ct] <- 1
    }
  }
  
  spatial_cont <- matrix(0, nrow = length(areas_common), ncol = length(common_cts),
                        dimnames = list(areas_common, common_cts))
  
  for (area in areas_common) {
    col_name <- paste0("prop_", area)
    for (ct in common_cts) {
      val <- spatial_stats_summary %>% filter(cell_type == ct) %>% pull(!!col_name)
      spatial_cont[area, ct] <- ifelse(length(val) > 0, val, 0)
    }
  }
  
  cors <- sapply(areas_common, function(area) {
    if (sum(cellex_binary[area, ]) > 0) {
      cor(cellex_binary[area, ], spatial_cont[area, ], method = "spearman")
    } else {
      NA
    }
  })
  
  mean_cor <- mean(cors, na.rm = TRUE)
  
  agreements <- sapply(common_cts, function(ct) {
    cellex_areas <- areas_common[cellex_binary[, ct] > 0]
    spatial_areas <- areas_common[spatial_cont[, ct] > 0.15]
    
    if (length(cellex_areas) == 0 & length(spatial_areas) == 0) return(TRUE)
    if (length(cellex_areas) == 0 | length(spatial_areas) == 0) return(FALSE)
    
    length(intersect(cellex_areas, spatial_areas)) > 0
  })
  
  mean_agree <- mean(agreements)
  n_pred <- sum(predictions$prediction != "negative")
  
  all_results <- rbind(all_results, data.frame(
    strategy = "relative_threshold",
    parameter = rel_thresh,
    correlation = mean_cor,
    mean_agreement = mean_agree,
    n_predictions = n_pred
  ))
  
  cat(sprintf("cor=%.3f, agree=%.2f, n=%d\n", mean_cor, mean_agree, n_pred))
}

# STRATEGY 3: Top N areas
cat("\nStrategy 3: Top N scoring areas\n")
for (top_n in top_n_values) {
  cat(sprintf("  Testing top %d... ", top_n))
  
  predictions <- data.frame(
    cell_type = colnames(cellex_spatial),
    stringsAsFactors = FALSE
  )
  
  for (i in 1:ncol(cellex_spatial)) {
    ct <- colnames(cellex_spatial)[i]
    scores <- cellex_spatial[, ct]
    scores[is.na(scores)] <- 0
    
    # Get top N areas (only those with positive scores)
    positive_scores <- scores[scores > 0]
    
    if (length(positive_scores) == 0) {
      predictions$prediction[i] <- "negative"
    } else {
      # Sort and take top N
      sorted_scores <- sort(positive_scores, decreasing = TRUE)
      n_to_take <- min(top_n, length(sorted_scores))
      top_areas <- names(sorted_scores)[1:n_to_take]
      
      predictions$prediction[i] <- paste(sort(top_areas), collapse = "+")
    }
  }
  
  # Same comparison
  common_cts <- intersect(spatial_stats_summary$cell_type, predictions$cell_type)
  
  cellex_binary <- matrix(0, nrow = length(areas_common), ncol = length(common_cts),
                         dimnames = list(areas_common, common_cts))
  
  for (ct in common_cts) {
    pred <- predictions$prediction[predictions$cell_type == ct]
    if (!is.na(pred) && pred != "negative") {
      pred_areas <- strsplit(pred, "\\+")[[1]]
      cellex_binary[pred_areas, ct] <- 1
    }
  }
  
  spatial_cont <- matrix(0, nrow = length(areas_common), ncol = length(common_cts),
                        dimnames = list(areas_common, common_cts))
  
  for (area in areas_common) {
    col_name <- paste0("prop_", area)
    for (ct in common_cts) {
      val <- spatial_stats_summary %>% filter(cell_type == ct) %>% pull(!!col_name)
      spatial_cont[area, ct] <- ifelse(length(val) > 0, val, 0)
    }
  }
  
  cors <- sapply(areas_common, function(area) {
    if (sum(cellex_binary[area, ]) > 0) {
      cor(cellex_binary[area, ], spatial_cont[area, ], method = "spearman")
    } else {
      NA
    }
  })
  
  mean_cor <- mean(cors, na.rm = TRUE)
  
  agreements <- sapply(common_cts, function(ct) {
    cellex_areas <- areas_common[cellex_binary[, ct] > 0]
    spatial_areas <- areas_common[spatial_cont[, ct] > 0.15]
    
    if (length(cellex_areas) == 0 & length(spatial_areas) == 0) return(TRUE)
    if (length(cellex_areas) == 0 | length(spatial_areas) == 0) return(FALSE)
    
    length(intersect(cellex_areas, spatial_areas)) > 0
  })
  
  mean_agree <- mean(agreements)
  n_pred <- sum(predictions$prediction != "negative")
  
  all_results <- rbind(all_results, data.frame(
    strategy = "top_n",
    parameter = top_n,
    correlation = mean_cor,
    mean_agreement = mean_agree,
    n_predictions = n_pred
  ))
  
  cat(sprintf("cor=%.3f, agree=%.2f, n=%d\n", mean_cor, mean_agree, n_pred))
}

cat("\n=== Optimization Complete ===\n\n")

# Find best overall
best_idx <- which.max(all_results$correlation)
best_strategy <- all_results[best_idx, ]

cat("=== OPTIMAL STRATEGY ===\n")
cat(sprintf("Strategy: %s\n", best_strategy$strategy))
cat(sprintf("Parameter: %.2f\n", best_strategy$parameter))
cat(sprintf("Correlation: %.3f\n", best_strategy$correlation))
cat(sprintf("Agreement: %.2f%%\n", 100 * best_strategy$mean_agreement))
cat(sprintf("Predictions: %d cell types\n\n", best_strategy$n_predictions))

cat("=== Top 5 Strategy-Parameter Combinations ===\n")
top_5 <- all_results %>%
  arrange(desc(correlation)) %>%
  head(5) %>%
  mutate(rank = 1:5)
print(top_5)

# Visualizations
pdf(paste0(atlaspath, "cellex_strategy_optimization.pdf"), width = 14, height = 10)

par(mfrow = c(2, 2))

# Plot each strategy
for (strat in unique(all_results$strategy)) {
  strat_data <- all_results %>% filter(strategy == strat)
  
  # Skip if all correlations are NA/NaN
  if (all(is.na(strat_data$correlation) | is.nan(strat_data$correlation))) {
    plot(1, type = "n", xlab = "", ylab = "", main = sprintf("Strategy: %s\n(No valid results)", strat))
    text(1, 1, "No valid correlations", cex = 1.5)
    next
  }
  
  # Get finite values for ylim
  valid_cors <- strat_data$correlation[is.finite(strat_data$correlation)]
  
  if (length(valid_cors) > 0) {
    plot(strat_data$parameter, strat_data$correlation,
         type = "b", pch = 19, col = "blue", lwd = 2,
         xlab = "Parameter Value", ylab = "Correlation",
         main = sprintf("Strategy: %s", strat),
         ylim = range(valid_cors, na.rm = TRUE))
    grid()
    
    # Mark best for this strategy
    best_for_strat <- strat_data[which.max(strat_data$correlation), ]
    if (!is.na(best_for_strat$correlation)) {
      points(best_for_strat$parameter, best_for_strat$correlation,
             pch = 8, col = "red", cex = 2, lwd = 2)
    }
  }
}

# Combined comparison - only valid strategies
par(mfrow = c(1, 1))
valid_results <- all_results %>%
  filter(!is.na(correlation) & !is.nan(correlation))

if (nrow(valid_results) > 0) {
  p <- ggplot(valid_results, aes(x = parameter, y = correlation, color = strategy)) +
    geom_line(size = 1.2) +
    geom_point(size = 3) +
    theme_minimal() +
    labs(title = "CELLEX Prediction Strategy Comparison",
         x = "Parameter Value", y = "Correlation with Spatial Stats",
         color = "Strategy") +
    theme(legend.position = "bottom")
  print(p)
}

dev.off()

# Save results
write.csv(all_results, 
          paste0(atlaspath, "cellex_strategy_optimization_results.csv"),
          row.names = FALSE)

write.csv(best_strategy,
          paste0(atlaspath, "optimal_cellex_strategy.csv"),
          row.names = FALSE)

cat("\n=== Files Saved ===\n")
cat("1. cellex_strategy_optimization.pdf - Strategy comparison plots\n")
cat("2. cellex_strategy_optimization_results.csv - All tested strategies\n")
cat("3. optimal_cellex_strategy.csv - Best strategy parameters\n")

cat("\n✓ Optimization complete!\n")
cat(sprintf("\nRecommendation: Use '%s' strategy with parameter %.2f\n",
            best_strategy$strategy, best_strategy$parameter))

# Compare CELLEX cell-type predictions vs Spatial Statistics predictions
# Both methods should predict projection targets for neurons

library(dplyr)
library(ggplot2)
library(tidyr)

# Read the CELLEX spatial enrichment data
cellex_spatial <- read.csv(paste0(spatialpath,"combined_spatial_cellex.csv"), 
                           row.names = 1, check.names = FALSE)

# Proper comparison: Aggregate spatial predictions to cell type level
# CELLEX predicts: "Does cell type X project to area Y?" (cell type level)
# Spatial stats predicts: "Does this specific cell project to area Y?" (single cell level)
# We need to aggregate spatial stats to cell type level for fair comparison

cat("=== Cell Type Level Comparison ===\n")
cat("CELLEX: Enrichment of AAV in each cell type → proxy for projection\n")
cat("Spatial Stats: Per-cell predictions → aggregate to cell type level\n\n")

# Get data for neurons at projection sites only
neuron_data <- data.frame(
  cell_type = seurat_merged$plot,
  ground_truth = seurat_merged$projection,
  spatial_pred = seurat_merged$projection_spatial_pred,
  neurotransmitter = seurat_merged$neurotransmitter,
  injection_filter = seurat_merged$injection_filter,
  run_index = seurat_merged$run_index,
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(neurotransmitter) & neurotransmitter != "Glia",
         !is.na(injection_filter) & injection_filter == "projection")

cat(sprintf("Total neurons at projection sites: %d\n", nrow(neuron_data)))
cat(sprintf("Unique cell types: %d\n", length(unique(neuron_data$cell_type))))

# Areas
areas <- c("PAG", "BST", "PBN", "LHA", "CEA", "PVT")

# 1. CELLEX METHOD: Cell type → Area enrichment scores
cellex_matrix <- cellex_spatial
cat("\n=== CELLEX Matrix (areas × cell types) ===\n")
cat(sprintf("Dimensions: %d areas × %d cell types\n", nrow(cellex_matrix), ncol(cellex_matrix)))

# 2. SPATIAL STATS METHOD: Aggregate per-cell predictions to cell type level
# For each cell type, calculate: what proportion project to each area?

spatial_celltype_summary <- neuron_data %>%
  filter(ground_truth != "negative") %>%  # only cells that actually project
  group_by(cell_type) %>%
  summarize(
    n_cells = n(),
    # For each area, calculate proportion of this cell type that projects there
    prop_PAG = sum(grepl("PAG", spatial_pred)) / n(),
    prop_BST = sum(grepl("BST", spatial_pred)) / n(),
    prop_PBN = sum(grepl("PBN", spatial_pred)) / n(),
    prop_LHA = sum(grepl("LHA", spatial_pred)) / n(),
    prop_CEA = sum(grepl("CEA", spatial_pred)) / n(),
    prop_PVT = sum(grepl("PVT", spatial_pred)) / n(),
    .groups = "drop"
  )

cat("\n=== Spatial Stats Aggregated (cell types) ===\n")
cat(sprintf("Cell types with projecting cells: %d\n", nrow(spatial_celltype_summary)))

# 3. GROUND TRUTH: What proportion of each cell type actually projects to each area?
ground_truth_summary <- neuron_data %>%
  filter(ground_truth != "negative") %>%
  group_by(cell_type) %>%
  summarize(
    n_cells = n(),
    gt_prop_PAG = sum(grepl("PAG", ground_truth)) / n(),
    gt_prop_BST = sum(grepl("BST", ground_truth)) / n(),
    gt_prop_PBN = sum(grepl("PBN", ground_truth)) / n(),
    gt_prop_LHA = sum(grepl("LHA", ground_truth)) / n(),
    gt_prop_CEA = sum(grepl("CEA", ground_truth)) / n(),
    gt_prop_PVT = sum(grepl("PVT", ground_truth)) / n(),
    .groups = "drop"
  )

# Convert to matrices for comparison
# Match cell types across all three
common_cell_types <- intersect(
  intersect(colnames(cellex_matrix), spatial_celltype_summary$cell_type),
  ground_truth_summary$cell_type
)

cat(sprintf("\nCell types present in all datasets: %d\n", length(common_cell_types)))

# Create matrices: rows = areas, columns = cell types
cellex_mat <- as.matrix(cellex_matrix[areas, common_cell_types])
cellex_mat[is.na(cellex_mat)] <- 0

spatial_mat <- matrix(NA, nrow = length(areas), ncol = length(common_cell_types),
                      dimnames = list(areas, common_cell_types))
for (area in areas) {
  col_name <- paste0("prop_", area)
  for (ct in common_cell_types) {
    val <- spatial_celltype_summary %>% 
      filter(cell_type == ct) %>% 
      pull(!!col_name)
    spatial_mat[area, ct] <- ifelse(length(val) > 0, val, 0)
  }
}

gt_mat <- matrix(NA, nrow = length(areas), ncol = length(common_cell_types),
                 dimnames = list(areas, common_cell_types))
for (area in areas) {
  col_name <- paste0("gt_prop_", area)
  for (ct in common_cell_types) {
    val <- ground_truth_summary %>% 
      filter(cell_type == ct) %>% 
      pull(!!col_name)
    gt_mat[area, ct] <- ifelse(length(val) > 0, val, 0)
  }
}

# 4. CORRELATION ANALYSIS
# For each area, correlate CELLEX enrichment with:
# a) Ground truth proportion
# b) Spatial stats proportion

cat("\n=== Correlation Analysis (per area) ===\n")
correlation_results <- data.frame(
  area = areas,
  cellex_vs_gt_cor = NA,
  cellex_vs_gt_pval = NA,
  spatial_vs_gt_cor = NA,
  spatial_vs_gt_pval = NA,
  cellex_vs_spatial_cor = NA,
  cellex_vs_spatial_pval = NA
)

for (i in seq_along(areas)) {
  area <- areas[i]
  
  # Get non-zero entries (cell types with data)
  has_data <- !is.na(cellex_mat[area, ]) & 
              !is.na(spatial_mat[area, ]) & 
              !is.na(gt_mat[area, ])
  
  if (sum(has_data) > 3) {  # Need at least 3 points for correlation
    # CELLEX vs Ground Truth
    cor_test <- cor.test(cellex_mat[area, has_data], 
                         gt_mat[area, has_data], 
                         method = "spearman")
    correlation_results$cellex_vs_gt_cor[i] <- cor_test$estimate
    correlation_results$cellex_vs_gt_pval[i] <- cor_test$p.value
    
    # Spatial vs Ground Truth
    cor_test <- cor.test(spatial_mat[area, has_data], 
                         gt_mat[area, has_data], 
                         method = "spearman")
    correlation_results$spatial_vs_gt_cor[i] <- cor_test$estimate
    correlation_results$spatial_vs_gt_pval[i] <- cor_test$p.value
    
    # CELLEX vs Spatial
    cor_test <- cor.test(cellex_mat[area, has_data], 
                         spatial_mat[area, has_data], 
                         method = "spearman")
    correlation_results$cellex_vs_spatial_cor[i] <- cor_test$estimate
    correlation_results$cellex_vs_spatial_pval[i] <- cor_test$p.value
  }
}

print(correlation_results)

# 5. VISUALIZATIONS
pdf(paste0(atlaspath, "celltype_level_comparison.pdf"), width = 16, height = 12)

# Heatmaps
par(mfrow = c(2, 2))

# Ground Truth
pheatmap(gt_mat, 
         main = "Ground Truth: Proportion of Each Cell Type Projecting to Area",
         cluster_rows = FALSE, cluster_cols = TRUE,
         color = colorRampPalette(c("white", "red"))(100),
         display_numbers = FALSE,
         fontsize = 8)

# CELLEX
pheatmap(cellex_mat, 
         main = "CELLEX: AAV Enrichment in Cell Types",
         cluster_rows = FALSE, cluster_cols = TRUE,
         color = colorRampPalette(c("white", "blue"))(100),
         display_numbers = FALSE,
         fontsize = 8)

# Spatial Stats
pheatmap(spatial_mat, 
         main = "Spatial Stats: Predicted Projection Proportions",
         cluster_rows = FALSE, cluster_cols = TRUE,
         color = colorRampPalette(c("white", "green"))(100),
         display_numbers = FALSE,
         fontsize = 8)

dev.off()

# Scatter plots per area
pdf(paste0(atlaspath, "celltype_correlations_by_area.pdf"), width = 14, height = 20)

par(mfrow = c(6, 2))

for (area in areas) {
  has_data <- !is.na(cellex_mat[area, ]) & 
              !is.na(spatial_mat[area, ]) & 
              !is.na(gt_mat[area, ])
  
  if (sum(has_data) > 3) {
    # CELLEX vs Ground Truth
    plot(cellex_mat[area, has_data], gt_mat[area, has_data],
         xlab = "CELLEX Enrichment", ylab = "Ground Truth Proportion",
         main = sprintf("%s: CELLEX vs Ground Truth\n(ρ=%.3f, p=%.3f)", 
                       area, 
                       cor(cellex_mat[area, has_data], gt_mat[area, has_data], 
                           method = "spearman", use = "complete.obs"),
                       cor.test(cellex_mat[area, has_data], gt_mat[area, has_data], 
                               method = "spearman")$p.value),
         pch = 19, col = "blue")
    abline(lm(gt_mat[area, has_data] ~ cellex_mat[area, has_data]), col = "red")
    
    # Spatial vs Ground Truth
    plot(spatial_mat[area, has_data], gt_mat[area, has_data],
         xlab = "Spatial Stats Proportion", ylab = "Ground Truth Proportion",
         main = sprintf("%s: Spatial Stats vs Ground Truth\n(ρ=%.3f, p=%.3f)", 
                       area,
                       cor(spatial_mat[area, has_data], gt_mat[area, has_data], 
                           method = "spearman", use = "complete.obs"),
                       cor.test(spatial_mat[area, has_data], gt_mat[area, has_data], 
                               method = "spearman")$p.value),
         pch = 19, col = "darkgreen")
    abline(lm(gt_mat[area, has_data] ~ spatial_mat[area, has_data]), col = "red")
  }
}

dev.off()

# Summary barplot
pdf(paste0(atlaspath, "method_comparison_summary.pdf"), width = 10, height = 6)

correlation_long <- correlation_results %>%
  select(area, cellex_vs_gt_cor, spatial_vs_gt_cor) %>%
  pivot_longer(cols = c(cellex_vs_gt_cor, spatial_vs_gt_cor),
               names_to = "method", values_to = "correlation") %>%
  mutate(method = ifelse(method == "cellex_vs_gt_cor", "CELLEX", "Spatial Stats"))

p <- ggplot(correlation_long, aes(x = area, y = correlation, fill = method)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(title = "Spearman Correlation with Ground Truth (Cell Type Level)",
       subtitle = "Higher = better prediction of which cell types project to each area",
       x = "Projection Area", y = "Spearman ρ", fill = "Method") +
  scale_fill_manual(values = c("CELLEX" = "steelblue", "Spatial Stats" = "darkgreen")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p)

dev.off()

# Save results
write.csv(correlation_results, paste0(atlaspath, "celltype_correlation_results.csv"), row.names = FALSE)
write.csv(spatial_celltype_summary, paste0(atlaspath, "spatial_celltype_aggregated.csv"), row.names = FALSE)
write.csv(ground_truth_summary, paste0(atlaspath, "ground_truth_celltype_summary.csv"), row.names = FALSE)

cat("\n=== Overall Method Performance ===\n")
cat(sprintf("Mean CELLEX correlation: %.3f\n", mean(correlation_results$cellex_vs_gt_cor, na.rm = TRUE)))
cat(sprintf("Mean Spatial Stats correlation: %.3f\n", mean(correlation_results$spatial_vs_gt_cor, na.rm = TRUE)))

cat("\n=== Files saved ===\n")
cat("1. celltype_level_comparison.pdf - Heatmaps of all three matrices\n")
cat("2. celltype_correlations_by_area.pdf - Scatter plots per area\n")
cat("3. method_comparison_summary.pdf - Bar plot of correlations\n")
cat("4. celltype_correlation_results.csv - Correlation statistics\n")
cat("5. spatial_celltype_aggregated.csv - Spatial stats aggregated to cell type level\n")
cat("6. ground_truth_celltype_aggregated.csv - Ground truth aggregated to cell type level\n")

cat("\n✓ Cell type level comparison complete!\n")

# Comprehensive CELLEX Comparison
# Compare 4 approaches against aggregated spatial stats (ground truth)
# 1. Spatial CELLEX
# 2. Single-cell CELLEX  
# 3. Spatial CELLEX
# 4. Combined (average of single-cell + spatial)

library(dplyr)
library(ggplot2)
library(tidyr)
library(pheatmap)

cat("=== Comprehensive CELLEX Method Comparison ===\n\n")

# Read data
cellex_spatial <- read.csv(paste0(spatialpath,"combined_spatial_cellex.csv"), 
                           row.names = 1, check.names = FALSE)
cellex_sc <- read.csv(paste0(spatialpath,"df_combined_new.csv"), 
                      header = FALSE, row.names = NULL)


# Fix single-cell CELLEX format (rows should be areas, first row is column names)
colnames(cellex_sc) <- as.character(cellex_sc[1, ])
cellex_sc <- cellex_sc[-1, ]
rownames(cellex_sc) <- c("BST", "CEA", "LHA", "PAG", "PBN", "PVH", "PVT")
cellex_sc <- cellex_sc[, -1]  # Remove first empty column if exists

# Convert to numeric
cellex_sc[] <- lapply(cellex_sc, function(x) as.numeric(as.character(x)))

cat("Single-cell CELLEX dimensions:", nrow(cellex_sc), "×", ncol(cellex_sc), "\n")
cat("Spatial CELLEX dimensions:", nrow(cellex_spatial), "×", ncol(cellex_spatial), "\n")

# We only compare the 6 areas that are in both datasets
areas_common <- c("PAG", "BST", "PBN", "LHA", "CEA", "PVT")

# Subset matrices to common areas
cellex_sc <- cellex_sc[areas_common, ]
cellex_spatial <- cellex_spatial[areas_common, ]

cat("\nAreas for comparison:", paste(areas_common, collapse = ", "), "\n")

# Get ground truth from spatial stats (aggregated to cell type level)
neuron_data <- data.frame(
  cell_type = seurat_merged$plot,
  ground_truth = seurat_merged$projection,
  neurotransmitter = seurat_merged$neurotransmitter,
  injection_filter = seurat_merged$injection_filter,
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(neurotransmitter) & neurotransmitter != "Glia",
         !is.na(injection_filter) & injection_filter == "projection",
         ground_truth != "negative")

# Aggregate to cell type level - proportion projecting to each area
ground_truth_summary <- neuron_data %>%
  group_by(cell_type) %>%
  summarize(
    n_cells = n(),
    gt_prop_PAG = sum(grepl("PAG", ground_truth)) / n(),
    gt_prop_BST = sum(grepl("BST", ground_truth)) / n(),
    gt_prop_PBN = sum(grepl("PBN", ground_truth)) / n(),
    gt_prop_LHA = sum(grepl("LHA", ground_truth)) / n(),
    gt_prop_CEA = sum(grepl("CEA", ground_truth)) / n(),
    gt_prop_PVT = sum(grepl("PVT", ground_truth)) / n(),
    .groups = "drop"
  )

cat(sprintf("\nGround truth cell types: %d\n", nrow(ground_truth_summary)))

# Find common cell types across all datasets
common_cell_types <- Reduce(intersect, list(
  colnames(cellex_sc),
  colnames(cellex_spatial),
  ground_truth_summary$cell_type
))

cat(sprintf("Common cell types across all datasets: %d\n", length(common_cell_types)))

# Create matrices for comparison (all: areas × cell types)
# 1. Single-cell CELLEX
sc_mat <- as.matrix(cellex_sc[areas_common, common_cell_types])
sc_mat[is.na(sc_mat)] <- 0

# 2. Spatial CELLEX
spatial_mat <- as.matrix(cellex_spatial[areas_common, common_cell_types])
spatial_mat[is.na(spatial_mat)] <- 0

# 3. Combined CELLEX (average of single-cell and spatial)
combined_mat <- (sc_mat + spatial_mat) / 2

# 4. Ground truth
gt_mat <- matrix(NA, nrow = length(areas_common), ncol = length(common_cell_types),
                 dimnames = list(areas_common, common_cell_types))
for (area in areas_common) {
  col_name <- paste0("gt_prop_", area)
  for (ct in common_cell_types) {
    val <- ground_truth_summary %>% 
      filter(cell_type == ct) %>% 
      pull(!!col_name)
    gt_mat[area, ct] <- ifelse(length(val) > 0, val, 0)
  }
}

# CORRELATION ANALYSIS
cat("\n=== Correlation Analysis (Spearman) ===\n\n")

methods <- c("Single-cell CELLEX", "Spatial CELLEX", "Combined CELLEX")
mat_list <- list(sc_mat, spatial_mat, combined_mat)

# Create results dataframe
correlation_results <- data.frame(
  area = rep(areas_common, each = 3),
  method = rep(methods, times = length(areas_common)),
  correlation = NA,
  p_value = NA,
  stringsAsFactors = FALSE
)

idx <- 1
for (area in areas_common) {
  for (i in 1:3) {
    mat <- mat_list[[i]]
    
    # Only use cell types with data in both matrices
    has_data <- !is.na(mat[area, ]) & !is.na(gt_mat[area, ])
    
    if (sum(has_data) > 3) {
      cor_test <- cor.test(mat[area, has_data], 
                          gt_mat[area, has_data], 
                          method = "spearman")
      correlation_results$correlation[idx] <- cor_test$estimate
      correlation_results$p_value[idx] <- cor_test$p.value
    }
    idx <- idx + 1
  }
}

# Print results table
cat("Correlation with Ground Truth (Aggregated Spatial Stats):\n\n")
results_wide <- correlation_results %>%
  select(area, method, correlation) %>%
  pivot_wider(names_from = method, values_from = correlation)
print(results_wide)

# Summary statistics
cat("\n=== Mean Correlation by Method ===\n")
mean_cors <- correlation_results %>%
  group_by(method) %>%
  summarize(
    mean_correlation = mean(correlation, na.rm = TRUE),
    sd_correlation = sd(correlation, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_correlation))
print(mean_cors)

# ADDITIONAL COMPARISON: Single-cell vs Spatial CELLEX
cat("\n=== Direct Comparison: Single-cell vs Spatial CELLEX ===\n")

sc_vs_spatial <- data.frame(
  area = areas_common,
  correlation = NA,
  p_value = NA
)

for (i in seq_along(areas_common)) {
  area <- areas_common[i]
  has_data <- !is.na(sc_mat[area, ]) & !is.na(spatial_mat[area, ])
  
  if (sum(has_data) > 3) {
    cor_test <- cor.test(sc_mat[area, has_data], 
                        spatial_mat[area, has_data], 
                        method = "spearman")
    sc_vs_spatial$correlation[i] <- cor_test$estimate
    sc_vs_spatial$p_value[i] <- cor_test$p.value
  }
}

print(sc_vs_spatial)

# VISUALIZATIONS
pdf(paste0(atlaspath, "comprehensive_cellex_comparison.pdf"), width = 16, height = 20)

# 1. Heatmaps of all methods
layout(matrix(c(1,2,3,4), nrow=2, byrow=TRUE))

# Use base R heatmap instead of pheatmap for better PDF compatibility
heatmap(gt_mat, 
        main = "Ground Truth (Spatial Stats Aggregated)",
        Rowv = NA, Colv = NA, scale = "none",
        col = colorRampPalette(c("white", "darkred"))(100),
        margins = c(8, 8), cexRow = 0.7, cexCol = 0.7)

heatmap(sc_mat, 
        main = "Single-cell CELLEX",
        Rowv = NA, Colv = NA, scale = "none",
        col = colorRampPalette(c("white", "blue"))(100),
        margins = c(8, 8), cexRow = 0.7, cexCol = 0.7)

heatmap(spatial_mat, 
        main = "Spatial CELLEX",
        Rowv = NA, Colv = NA, scale = "none",
        col = colorRampPalette(c("white", "darkgreen"))(100),
        margins = c(8, 8), cexRow = 0.7, cexCol = 0.7)

heatmap(combined_mat, 
        main = "Combined CELLEX (SC + Spatial Average)",
        Rowv = NA, Colv = NA, scale = "none",
        col = colorRampPalette(c("white", "purple"))(100),
        margins = c(8, 8), cexRow = 0.7, cexCol = 0.7)

# 2. Bar plot comparing methods
par(mfrow = c(1, 1))  # Reset layout
p1 <- ggplot(correlation_results, 
             aes(x = area, y = correlation, fill = method)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(title = "Correlation with Ground Truth by Area and Method",
       x = "Area", y = "Spearman Correlation", fill = "Method") +
  scale_fill_manual(values = c("Single-cell CELLEX" = "steelblue",
                               "Spatial CELLEX" = "darkgreen",
                               "Combined CELLEX" = "purple")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p1)

# 3. Mean correlation comparison
p2 <- ggplot(mean_cors, aes(x = reorder(method, mean_correlation), 
                            y = mean_correlation, fill = method)) +
  geom_col() +
  geom_errorbar(aes(ymin = mean_correlation - sd_correlation,
                    ymax = mean_correlation + sd_correlation),
                width = 0.2) +
  coord_flip() +
  theme_minimal() +
  labs(title = "Overall Method Performance",
       subtitle = "Mean ± SD Spearman Correlation with Ground Truth",
       x = "Method", y = "Mean Correlation") +
  scale_fill_manual(values = c("Single-cell CELLEX" = "steelblue",
                               "Spatial CELLEX" = "darkgreen",
                               "Combined CELLEX" = "purple")) +
  theme(legend.position = "none")
print(p2)

# 4. Scatter plots - each method vs ground truth (for best area)
best_area <- results_wide %>%
  pivot_longer(cols = -area, names_to = "method", values_to = "cor") %>%
  filter(method == "Combined CELLEX") %>%
  arrange(desc(cor)) %>%
  slice(1) %>%
  pull(area)

par(mfrow = c(2, 2))

# Ground truth vs SC
has_data <- !is.na(sc_mat[best_area, ]) & !is.na(gt_mat[best_area, ])
plot(sc_mat[best_area, has_data], gt_mat[best_area, has_data],
     xlab = "Single-cell CELLEX", ylab = "Ground Truth",
     main = sprintf("%s: SC CELLEX vs GT", best_area),
     pch = 19, col = "steelblue")
abline(lm(gt_mat[best_area, has_data] ~ sc_mat[best_area, has_data]), col = "red")

# Ground truth vs Spatial
plot(spatial_mat[best_area, has_data], gt_mat[best_area, has_data],
     xlab = "Spatial CELLEX", ylab = "Ground Truth",
     main = sprintf("%s: Spatial CELLEX vs GT", best_area),
     pch = 19, col = "darkgreen")
abline(lm(gt_mat[best_area, has_data] ~ spatial_mat[best_area, has_data]), col = "red")

# Ground truth vs Combined
plot(combined_mat[best_area, has_data], gt_mat[best_area, has_data],
     xlab = "Combined CELLEX", ylab = "Ground Truth",
     main = sprintf("%s: Combined CELLEX vs GT", best_area),
     pch = 19, col = "purple")
abline(lm(gt_mat[best_area, has_data] ~ combined_mat[best_area, has_data]), col = "red")

# SC vs Spatial
plot(sc_mat[best_area, has_data], spatial_mat[best_area, has_data],
     xlab = "Single-cell CELLEX", ylab = "Spatial CELLEX",
     main = sprintf("%s: SC vs Spatial CELLEX", best_area),
     pch = 19, col = "black")
abline(lm(spatial_mat[best_area, has_data] ~ sc_mat[best_area, has_data]), col = "red")

dev.off()

# Save results
write.csv(correlation_results, paste0(atlaspath, "comprehensive_cellex_correlations.csv"), row.names = FALSE)
write.csv(sc_vs_spatial, paste0(atlaspath, "sc_vs_spatial_cellex_correlation.csv"), row.names = FALSE)

# Save matrices
write.csv(sc_mat, paste0(atlaspath, "singlecell_cellex_matrix.csv"), row.names = TRUE)
write.csv(spatial_mat, paste0(atlaspath, "spatial_cellex_matrix.csv"), row.names = TRUE)
write.csv(combined_mat, paste0(atlaspath, "combined_cellex_matrix.csv"), row.names = TRUE)
write.csv(gt_mat, paste0(atlaspath, "ground_truth_matrix.csv"), row.names = TRUE)

cat("\n=== Files Saved ===\n")
cat("1. comprehensive_cellex_comparison.pdf - All visualizations\n")
cat("2. comprehensive_cellex_correlations.csv - Full correlation results\n")
cat("3. sc_vs_spatial_cellex_correlation.csv - Direct SC vs Spatial comparison\n")
cat("4. singlecell_cellex_matrix.csv - SC CELLEX values\n")
cat("5. spatial_cellex_matrix.csv - Spatial CELLEX values\n")
cat("6. combined_cellex_matrix.csv - Combined CELLEX values\n")
cat("7. ground_truth_matrix.csv - Ground truth proportions\n")

cat("\n=== CONCLUSION ===\n")
best_method <- mean_cors$method[1]
best_cor <- mean_cors$mean_correlation[1]
cat(sprintf("Best performing method: %s (mean ρ = %.3f)\n", best_method, best_cor))

cat("\n✓ Comprehensive comparison complete!\n")

# Optimize Combined CELLEX Metric for Cell Type Projection Prediction
# Goal: Find best combination of SC + Spatial CELLEX and optimal threshold
# Test different weighting schemes and noise thresholds

library(dplyr)
library(ggplot2)
library(tidyr)

cat("=== Optimizing Combined CELLEX Projection Predictor ===\n\n")

# Read data
cellex_spatial <- read.csv(paste0(spatialpath,"combined_spatial_cellex.csv"), 
                           row.names = 1, check.names = FALSE)
cellex_sc <- read.csv(paste0(spatialpath,"df_combined_new.csv"), 
                      header = FALSE, row.names = NULL)

# Fix single-cell CELLEX format
colnames(cellex_sc) <- as.character(cellex_sc[1, ])
cellex_sc <- cellex_sc[-1, ]
rownames(cellex_sc) <- c("BST", "CEA", "LHA", "PAG", "PBN", "PVH", "PVT")
cellex_sc[] <- lapply(cellex_sc, function(x) as.numeric(as.character(x)))

# Common areas and cell types
areas_common <- c("PAG", "BST", "PBN", "LHA", "CEA", "PVT")
cellex_sc <- cellex_sc[areas_common, ]
cellex_spatial <- cellex_spatial[areas_common, ]

common_cell_types <- intersect(colnames(cellex_sc), colnames(cellex_spatial))

# Create matrices
sc_mat <- as.matrix(cellex_sc[areas_common, common_cell_types])
sc_mat[is.na(sc_mat)] <- 0

spatial_mat <- as.matrix(cellex_spatial[areas_common, common_cell_types])
spatial_mat[is.na(spatial_mat)] <- 0

# Get spatial stats for agreement (not ground truth)
neuron_data <- data.frame(
  cell_type = seurat_merged$plot,
  projection = seurat_merged$projection,
  neurotransmitter = seurat_merged$neurotransmitter,
  injection_filter = seurat_merged$injection_filter,
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(neurotransmitter) & neurotransmitter != "Glia",
         !is.na(injection_filter) & injection_filter == "projection",
         projection != "negative")

spatial_stats_summary <- neuron_data %>%
  group_by(cell_type) %>%
  summarize(
    n_cells = n(),
    prop_PAG = sum(grepl("PAG", projection)) / n(),
    prop_BST = sum(grepl("BST", projection)) / n(),
    prop_PBN = sum(grepl("PBN", projection)) / n(),
    prop_LHA = sum(grepl("LHA", projection)) / n(),
    prop_CEA = sum(grepl("CEA", projection)) / n(),
    prop_PVT = sum(grepl("PVT", projection)) / n(),
    .groups = "drop"
  )

# Create spatial stats matrix
common_with_stats <- intersect(common_cell_types, spatial_stats_summary$cell_type)
cat(sprintf("Cell types with all data: %d\n", length(common_with_stats)))

spatial_stats_mat <- matrix(NA, nrow = length(areas_common), ncol = length(common_with_stats),
                             dimnames = list(areas_common, common_with_stats))
for (area in areas_common) {
  col_name <- paste0("prop_", area)
  for (ct in common_with_stats) {
    val <- spatial_stats_summary %>% filter(cell_type == ct) %>% pull(!!col_name)
    spatial_stats_mat[area, ct] <- ifelse(length(val) > 0, val, 0)
  }
}

# Subset to common cell types with stats
sc_mat <- sc_mat[, common_with_stats]
spatial_mat <- spatial_mat[, common_with_stats]

cat(sprintf("Testing on %d cell types across %d areas\n\n", 
            length(common_with_stats), length(areas_common)))

# OPTIMIZATION PARAMETERS
weight_sc_values <- seq(0, 1, by = 0.1)  # Weight for SC (1-weight for Spatial)
threshold_values <- seq(0, 0.5, by = 0.05)  # Noise threshold to set to 0

# Store results
optimization_results <- expand.grid(
  weight_sc = weight_sc_values,
  threshold = threshold_values,
  correlation_with_stats = NA,
  mean_max_score = NA,
  prop_clear_winner = NA,
  stringsAsFactors = FALSE
)

cat("Testing", nrow(optimization_results), "parameter combinations...\n")

for (i in 1:nrow(optimization_results)) {
  w_sc <- optimization_results$weight_sc[i]
  w_spatial <- 1 - w_sc
  thresh <- optimization_results$threshold[i]
  
  # Create combined matrix
  combined <- w_sc * sc_mat + w_spatial * spatial_mat
  
  # Apply threshold (set values below threshold to 0)
  combined[combined < thresh] <- 0
  
  # Calculate metrics
  
  # 1. Correlation with spatial stats
  correlations <- numeric(length(areas_common))
  for (j in seq_along(areas_common)) {
    area <- areas_common[j]
    has_data <- !is.na(combined[area, ]) & !is.na(spatial_stats_mat[area, ])
    if (sum(has_data) > 3) {
      correlations[j] <- cor(combined[area, has_data], 
                            spatial_stats_mat[area, has_data], 
                            method = "spearman", use = "complete.obs")
    } else {
      correlations[j] <- NA
    }
  }
  optimization_results$correlation_with_stats[i] <- mean(correlations, na.rm = TRUE)
  
  # 2. Mean maximum score (how strong are predictions?)
  max_scores <- apply(combined, 2, max, na.rm = TRUE)
  optimization_results$mean_max_score[i] <- mean(max_scores[is.finite(max_scores)])
  
  # 3. Proportion with clear winner (max > 2x second highest)
  clear_winners <- sapply(1:ncol(combined), function(j) {
    scores <- combined[, j]
    scores_sorted <- sort(scores[scores > 0], decreasing = TRUE)
    if (length(scores_sorted) >= 2) {
      scores_sorted[1] > 2 * scores_sorted[2]
    } else if (length(scores_sorted) == 1) {
      TRUE  # Only one projection is clear
    } else {
      FALSE
    }
  })
  optimization_results$prop_clear_winner[i] <- mean(clear_winners)
}

cat("✓ Optimization complete\n\n")

# Find best parameter combinations
cat("=== Top 10 Parameter Combinations ===\n")
cat("Ranked by correlation with spatial stats\n\n")

top_10 <- optimization_results %>%
  arrange(desc(correlation_with_stats)) %>%
  head(10) %>%
  mutate(
    weight_spatial = 1 - weight_sc,
    rank = 1:10
  ) %>%
  select(rank, weight_sc, weight_spatial, threshold, correlation_with_stats, 
         mean_max_score, prop_clear_winner)

print(top_10)

# Best overall parameters
best_params <- optimization_results %>%
  filter(correlation_with_stats == max(correlation_with_stats, na.rm = TRUE)) %>%
  slice(1)

cat("\n=== OPTIMAL PARAMETERS ===\n")
cat(sprintf("SC weight: %.2f\n", best_params$weight_sc))
cat(sprintf("Spatial weight: %.2f\n", 1 - best_params$weight_sc))
cat(sprintf("Noise threshold: %.3f\n", best_params$threshold))
cat(sprintf("Correlation with spatial stats: %.3f\n", best_params$correlation_with_stats))
cat(sprintf("Mean max score: %.3f\n", best_params$mean_max_score))
cat(sprintf("Proportion clear winners: %.2f%%\n", 100 * best_params$prop_clear_winner))

# Create optimal combined matrix
optimal_combined <- best_params$weight_sc * sc_mat + 
                    (1 - best_params$weight_sc) * spatial_mat
optimal_combined[optimal_combined < best_params$threshold] <- 0

# Analyze the optimal combined matrix
cat("\n=== Analysis of Optimal Combined Matrix ===\n")

# For each cell type, predict projection area(s)
predictions <- data.frame(
  cell_type = colnames(optimal_combined),
  predicted_area = NA,
  max_score = NA,
  second_score = NA,
  confidence = NA,  # ratio of max to second
  stringsAsFactors = FALSE
)

for (i in 1:ncol(optimal_combined)) {
  scores <- optimal_combined[, i]
  scores_sorted <- sort(scores, decreasing = TRUE)
  max_area <- names(scores)[which.max(scores)]
  
  predictions$predicted_area[i] <- ifelse(scores_sorted[1] > 0, max_area, "none")
  predictions$max_score[i] <- scores_sorted[1]
  predictions$second_score[i] <- ifelse(length(scores_sorted) > 1, scores_sorted[2], 0)
  predictions$confidence[i] <- ifelse(scores_sorted[2] > 0, 
                                      scores_sorted[1] / scores_sorted[2], 
                                      Inf)
}

# Add spatial stats info for comparison
predictions <- predictions %>%
  left_join(spatial_stats_summary, by = "cell_type") %>%
  rowwise() %>%
  mutate(
    spatial_stats_areas = {
      areas_vec <- c(prop_PAG, prop_BST, prop_PBN, prop_LHA, prop_CEA, prop_PVT)
      names(areas_vec) <- areas_common
      top_areas <- names(sort(areas_vec, decreasing = TRUE)[1:2])
      paste(top_areas[areas_vec[top_areas] > 0.1], collapse = "+")
    },
    agreement = grepl(predicted_area, spatial_stats_areas) | 
                grepl(spatial_stats_areas, predicted_area)
  ) %>%
  ungroup()

cat("\nPrediction summary:\n")
cat(sprintf("Cell types with predictions: %d / %d\n", 
            sum(predictions$predicted_area != "none"), 
            nrow(predictions)))
cat(sprintf("High confidence (>2x ratio): %d\n", 
            sum(predictions$confidence > 2 & is.finite(predictions$confidence))))
cat(sprintf("Agreement with spatial stats: %.1f%%\n", 
            100 * mean(predictions$agreement, na.rm = TRUE)))

# Show top predictions
cat("\n=== Top 20 Cell Type Predictions ===\n")
top_predictions <- predictions %>%
  filter(predicted_area != "none") %>%
  arrange(desc(max_score)) %>%
  select(cell_type, predicted_area, max_score, confidence, 
         spatial_stats_areas, agreement, n_cells) %>%
  head(20)
print(top_predictions)

# VISUALIZATIONS
pdf(paste0(atlaspath, "cellex_optimization_results.pdf"), width = 14, height = 18)

# 1. Heatmap of correlation vs parameters
cor_matrix <- optimization_results %>%
  select(weight_sc, threshold, correlation_with_stats) %>%
  pivot_wider(names_from = threshold, values_from = correlation_with_stats) %>%
  select(-weight_sc) %>%
  as.matrix()

rownames(cor_matrix) <- sprintf("SC=%.1f", weight_sc_values)

image(1:ncol(cor_matrix), 1:nrow(cor_matrix), t(cor_matrix),
      col = colorRampPalette(c("blue", "white", "red"))(100),
      xlab = "Noise Threshold", ylab = "SC Weight",
      main = "Correlation with Spatial Stats\n(Red = Better)",
      axes = FALSE)
axis(1, at = 1:ncol(cor_matrix), labels = colnames(cor_matrix))
axis(2, at = 1:nrow(cor_matrix), labels = rownames(cor_matrix), las = 1)
contour(1:ncol(cor_matrix), 1:nrow(cor_matrix), t(cor_matrix), add = TRUE)

# 2. Parameter effects
par(mfrow = c(2, 2))

# Effect of SC weight (at optimal threshold)
opt_thresh_data <- optimization_results %>%
  filter(threshold == best_params$threshold)

plot(opt_thresh_data$weight_sc, opt_thresh_data$correlation_with_stats,
     type = "b", pch = 19, col = "blue",
     xlab = "SC Weight", ylab = "Correlation",
     main = sprintf("Effect of SC Weight\n(threshold=%.2f)", best_params$threshold))
abline(v = best_params$weight_sc, col = "red", lty = 2)

# Effect of threshold (at optimal weight)
opt_weight_data <- optimization_results %>%
  filter(weight_sc == best_params$weight_sc)

plot(opt_weight_data$threshold, opt_weight_data$correlation_with_stats,
     type = "b", pch = 19, col = "darkgreen",
     xlab = "Threshold", ylab = "Correlation",
     main = sprintf("Effect of Threshold\n(SC weight=%.2f)", best_params$weight_sc))
abline(v = best_params$threshold, col = "red", lty = 2)

# Mean max score vs correlation
plot(optimization_results$mean_max_score, optimization_results$correlation_with_stats,
     pch = 19, col = rgb(0, 0, 1, 0.3),
     xlab = "Mean Max Score", ylab = "Correlation",
     main = "Prediction Strength vs Agreement")
points(best_params$mean_max_score, best_params$correlation_with_stats, 
       pch = 19, col = "red", cex = 2)

# Clear winners vs correlation
plot(optimization_results$prop_clear_winner, optimization_results$correlation_with_stats,
     pch = 19, col = rgb(0, 0, 1, 0.3),
     xlab = "Proportion Clear Winners", ylab = "Correlation",
     main = "Prediction Clarity vs Agreement")
points(best_params$prop_clear_winner, best_params$correlation_with_stats, 
       pch = 19, col = "red", cex = 2)

# 3. Optimal combined matrix heatmap
par(mfrow = c(1, 1))
heatmap(optimal_combined,
        main = sprintf("Optimal Combined CELLEX\n(SC=%.2f, Spatial=%.2f, threshold=%.3f)",
                      best_params$weight_sc, 1 - best_params$weight_sc, best_params$threshold),
        col = colorRampPalette(c("white", "red", "darkred"))(100),
        scale = "none", margins = c(8, 8),
        cexRow = 0.8, cexCol = 0.6)

# 4. Prediction confidence distribution
par(mfrow = c(1, 2))
hist(predictions$max_score[predictions$predicted_area != "none"],
     breaks = 30, col = "steelblue",
     xlab = "Max CELLEX Score", main = "Distribution of Prediction Scores")

conf_finite <- predictions$confidence[is.finite(predictions$confidence)]
hist(log2(conf_finite),
     breaks = 30, col = "darkgreen",
     xlab = "log2(Confidence Ratio)", 
     main = "Distribution of Prediction Confidence\n(max/second)")
abline(v = log2(2), col = "red", lty = 2, lwd = 2)
text(log2(2), par("usr")[4] * 0.9, "2x threshold", pos = 4, col = "red")

dev.off()

# Save results
write.csv(optimization_results, paste0(atlaspath, "cellex_optimization_grid.csv"), row.names = FALSE)
write.csv(predictions, paste0(atlaspath, "celltype_projection_predictions.csv"), row.names = FALSE)
write.csv(optimal_combined, paste0(atlaspath, "optimal_combined_cellex_matrix.csv"), row.names = TRUE)

# Save optimal parameters
optimal_params_summary <- data.frame(
  parameter = c("SC_weight", "Spatial_weight", "Threshold", 
                "Correlation_with_stats", "Mean_max_score", "Prop_clear_winner"),
  value = c(best_params$weight_sc, 1 - best_params$weight_sc, best_params$threshold,
            best_params$correlation_with_stats, best_params$mean_max_score, 
            best_params$prop_clear_winner)
)
write.csv(optimal_params_summary, paste0(atlaspath, "optimal_cellex_parameters.csv"), row.names = FALSE)

cat("\n=== Files Saved ===\n")
cat("1. cellex_optimization_results.pdf - Optimization visualizations\n")
cat("2. cellex_optimization_grid.csv - All tested parameter combinations\n")
cat("3. celltype_projection_predictions.csv - Predictions for each cell type\n")
cat("4. optimal_combined_cellex_matrix.csv - Final combined CELLEX matrix\n")
cat("5. optimal_cellex_parameters.csv - Best parameters summary\n")

cat("\n✓ Optimization complete!\n")

# Unbiased Analysis: Does Single-Cell CELLEX Add Independent Predictive Value?
# Problem: Spatial CELLEX vs Spatial Stats has circular dependency (same data)
# Question: What unique information does SC CELLEX provide?

library(dplyr)
library(ggplot2)
library(tidyr)

cat("=== Independent Predictive Value of Single-Cell CELLEX ===\n\n")

# Read data
cellex_spatial <- read.csv(paste0(spatialpath,"combined_spatial_cellex.csv"), 
                           row.names = 1, check.names = FALSE)
cellex_sc <- read.csv(paste0(spatialpath,"df_combined_new.csv"), 
                      header = FALSE, row.names = NULL)

# Fix SC CELLEX format
colnames(cellex_sc) <- as.character(cellex_sc[1, ])
cellex_sc <- cellex_sc[-1, ]
rownames(cellex_sc) <- c("BST", "CEA", "LHA", "PAG", "PBN", "PVH", "PVT")
cellex_sc[] <- lapply(cellex_sc, function(x) as.numeric(as.character(x)))

areas_common <- c("PAG", "BST", "PBN", "LHA", "CEA", "PVT")
cellex_sc <- cellex_sc[areas_common, ]
cellex_spatial <- cellex_spatial[areas_common, ]

common_cell_types <- intersect(colnames(cellex_sc), colnames(cellex_spatial))

# Create matrices
sc_mat <- as.matrix(cellex_sc[areas_common, common_cell_types])
sc_mat[is.na(sc_mat)] <- 0

spatial_mat <- as.matrix(cellex_spatial[areas_common, common_cell_types])
spatial_mat[is.na(spatial_mat)] <- 0

cat(sprintf("Analyzing %d cell types across %d areas\n\n", 
            length(common_cell_types), length(areas_common)))

# ANALYSIS 1: Where do the two methods DISAGREE?
# These are cases where SC CELLEX might provide independent information

cat("=== ANALYSIS 1: Where SC and Spatial CELLEX Disagree ===\n\n")

disagreement_analysis <- data.frame(
  cell_type = character(),
  sc_top_area = character(),
  sc_top_score = numeric(),
  spatial_top_area = character(),
  spatial_top_score = numeric(),
  agree = logical(),
  sc_adds_value = logical(),
  stringsAsFactors = FALSE
)

for (ct in common_cell_types) {
  sc_scores <- sc_mat[, ct]
  spatial_scores <- spatial_mat[, ct]
  
  sc_top <- names(which.max(sc_scores))
  sc_top_val <- max(sc_scores)
  
  spatial_top <- names(which.max(spatial_scores))
  spatial_top_val <- max(spatial_scores)
  
  # Do they agree on top area?
  agree <- (sc_top == spatial_top)
  
  # Does SC suggest something different with high confidence?
  # (SC top score > 0.3 AND different from spatial)
  sc_adds_value <- (!agree & sc_top_val > 0.3)
  
  disagreement_analysis <- rbind(disagreement_analysis, data.frame(
    cell_type = ct,
    sc_top_area = sc_top,
    sc_top_score = sc_top_val,
    spatial_top_area = spatial_top,
    spatial_top_score = spatial_top_val,
    agree = agree,
    sc_adds_value = sc_adds_value,
    stringsAsFactors = FALSE
  ))
}

cat("Agreement between SC and Spatial CELLEX:\n")
cat(sprintf("  Agree on top area: %d / %d (%.1f%%)\n", 
            sum(disagreement_analysis$agree), 
            nrow(disagreement_analysis),
            100 * mean(disagreement_analysis$agree)))
cat(sprintf("  SC suggests different area with confidence: %d (%.1f%%)\n",
            sum(disagreement_analysis$sc_adds_value),
            100 * mean(disagreement_analysis$sc_adds_value)))

cat("\nCases where SC CELLEX suggests different area (score > 0.3):\n")
different_predictions <- disagreement_analysis %>%
  filter(sc_adds_value) %>%
  arrange(desc(sc_top_score))
print(different_predictions)

# ANALYSIS 2: Complementary information
# Where SC is strong but Spatial is weak (and vice versa)

cat("\n\n=== ANALYSIS 2: Complementary Strengths ===\n\n")

complementary <- data.frame(
  cell_type = character(),
  area = character(),
  sc_score = numeric(),
  spatial_score = numeric(),
  pattern = character(),
  stringsAsFactors = FALSE
)

for (ct in common_cell_types) {
  for (area in areas_common) {
    sc_val <- sc_mat[area, ct]
    spatial_val <- spatial_mat[area, ct]
    
    # Categorize the pattern
    pattern <- "both_low"
    if (sc_val > 0.3 & spatial_val > 0.3) pattern <- "both_high"
    if (sc_val > 0.3 & spatial_val < 0.1) pattern <- "sc_only"
    if (spatial_val > 0.3 & sc_val < 0.1) pattern <- "spatial_only"
    
    if (pattern %in% c("sc_only", "spatial_only", "both_high")) {
      complementary <- rbind(complementary, data.frame(
        cell_type = ct,
        area = area,
        sc_score = sc_val,
        spatial_score = spatial_val,
        pattern = pattern,
        stringsAsFactors = FALSE
      ))
    }
  }
}

pattern_summary <- complementary %>%
  group_by(pattern) %>%
  summarize(
    n_cases = n(),
    mean_sc = mean(sc_score),
    mean_spatial = mean(spatial_score),
    .groups = "drop"
  )

cat("Complementary information patterns:\n")
print(pattern_summary)

cat("\nCases where SC is HIGH but Spatial is LOW:\n")
sc_only <- complementary %>%
  filter(pattern == "sc_only") %>%
  arrange(desc(sc_score))
print(head(sc_only, 10))

cat("\nCases where Spatial is HIGH but SC is LOW:\n")
spatial_only <- complementary %>%
  filter(pattern == "spatial_only") %>%
  arrange(desc(spatial_score))
print(head(spatial_only, 10))

# ANALYSIS 3: Does averaging help when both are informative?

cat("\n\n=== ANALYSIS 3: Value of Averaging When Both Methods Agree ===\n\n")

# For cell types where both methods agree, does averaging improve prediction clarity?

both_high_cases <- complementary %>%
  filter(pattern == "both_high")

cat(sprintf("Cases where both SC and Spatial are high: %d\n", nrow(both_high_cases)))

# For these cases, compare:
# 1. Spatial alone
# 2. SC alone  
# 3. Average

evaluation <- data.frame(
  cell_type = character(),
  spatial_winner = character(),
  sc_winner = character(),
  avg_winner = character(),
  spatial_clarity = numeric(),
  sc_clarity = numeric(),
  avg_clarity = numeric(),
  stringsAsFactors = FALSE
)

for (ct in unique(both_high_cases$cell_type)) {
  spatial_scores <- spatial_mat[, ct]
  sc_scores <- sc_mat[, ct]
  avg_scores <- (spatial_scores + sc_scores) / 2
  
  # Get winners
  spatial_win <- names(which.max(spatial_scores))
  sc_win <- names(which.max(sc_scores))
  avg_win <- names(which.max(avg_scores))
  
  # Calculate clarity (ratio of top to second)
  spatial_sorted <- sort(spatial_scores, decreasing = TRUE)
  sc_sorted <- sort(sc_scores, decreasing = TRUE)
  avg_sorted <- sort(avg_scores, decreasing = TRUE)
  
  spatial_clarity <- spatial_sorted[1] / (spatial_sorted[2] + 0.001)
  sc_clarity <- sc_sorted[1] / (sc_sorted[2] + 0.001)
  avg_clarity <- avg_sorted[1] / (avg_sorted[2] + 0.001)
  
  evaluation <- rbind(evaluation, data.frame(
    cell_type = ct,
    spatial_winner = spatial_win,
    sc_winner = sc_win,
    avg_winner = avg_win,
    spatial_clarity = spatial_clarity,
    sc_clarity = sc_clarity,
    avg_clarity = avg_clarity,
    stringsAsFactors = FALSE
  ))
}

cat("\nDoes averaging improve prediction clarity?\n")
cat(sprintf("Mean clarity - Spatial alone: %.2f\n", mean(evaluation$spatial_clarity)))
cat(sprintf("Mean clarity - SC alone: %.2f\n", mean(evaluation$sc_clarity)))
cat(sprintf("Mean clarity - Average: %.2f\n", mean(evaluation$avg_clarity)))

cat("\nWinner agreement:\n")
cat(sprintf("Spatial = SC: %d / %d (%.1f%%)\n",
            sum(evaluation$spatial_winner == evaluation$sc_winner),
            nrow(evaluation),
            100 * mean(evaluation$spatial_winner == evaluation$sc_winner)))
cat(sprintf("Average = Spatial: %d / %d (%.1f%%)\n",
            sum(evaluation$avg_winner == evaluation$spatial_winner),
            nrow(evaluation),
            100 * mean(evaluation$avg_winner == evaluation$spatial_winner)))
cat(sprintf("Average = SC: %d / %d (%.1f%%)\n",
            sum(evaluation$avg_winner == evaluation$sc_winner),
            nrow(evaluation),
            100 * mean(evaluation$avg_winner == evaluation$sc_winner)))

# ANALYSIS 4: Unique SC predictions validated by spatial stats

cat("\n\n=== ANALYSIS 4: SC Predictions NOT in Spatial CELLEX ===\n\n")

# Get spatial stats for validation
neuron_data <- data.frame(
  cell_type = seurat_merged$plot,
  projection = seurat_merged$projection,
  neurotransmitter = seurat_merged$neurotransmitter,
  injection_filter = seurat_merged$injection_filter,
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(neurotransmitter) & neurotransmitter != "Glia",
         !is.na(injection_filter) & injection_filter == "projection",
         projection != "negative")

spatial_stats_summary <- neuron_data %>%
  group_by(cell_type) %>%
  summarize(
    n_cells = n(),
    prop_PAG = sum(grepl("PAG", projection)) / n(),
    prop_BST = sum(grepl("BST", projection)) / n(),
    prop_PBN = sum(grepl("PBN", projection)) / n(),
    prop_LHA = sum(grepl("LHA", projection)) / n(),
    prop_CEA = sum(grepl("CEA", projection)) / n(),
    prop_PVT = sum(grepl("PVT", projection)) / n(),
    .groups = "drop"
  )

# For SC-only predictions, check if spatial stats supports them
sc_unique <- sc_only %>%
  left_join(spatial_stats_summary, by = "cell_type") %>%
  rowwise() %>%
  mutate(
    spatial_stats_score = case_when(
      area == "PAG" ~ prop_PAG,
      area == "BST" ~ prop_BST,
      area == "PBN" ~ prop_PBN,
      area == "LHA" ~ prop_LHA,
      area == "CEA" ~ prop_CEA,
      area == "PVT" ~ prop_PVT
    ),
    sc_validated = spatial_stats_score > 0.15
  ) %>%
  ungroup()

cat("SC CELLEX predictions where Spatial CELLEX is weak:\n")
cat(sprintf("Total cases: %d\n", nrow(sc_unique)))
cat(sprintf("Validated by spatial stats (>15%% projection): %d (%.1f%%)\n",
            sum(sc_unique$sc_validated, na.rm = TRUE),
            100 * mean(sc_unique$sc_validated, na.rm = TRUE)))

cat("\nValidated SC-only predictions:\n")
validated <- sc_unique %>%
  filter(sc_validated) %>%
  arrange(desc(sc_score)) %>%
  select(cell_type, area, sc_score, spatial_score, spatial_stats_score, n_cells)
print(validated)

# VISUALIZATIONS
pdf(paste0(atlaspath, "sc_cellex_independent_value.pdf"), width = 14, height = 16)

# 1. Scatter: SC vs Spatial CELLEX (all cell type-area pairs)
all_pairs <- data.frame(
  sc = as.vector(sc_mat),
  spatial = as.vector(spatial_mat),
  cell_type = rep(colnames(sc_mat), each = nrow(sc_mat)),
  area = rep(rownames(sc_mat), times = ncol(sc_mat))
)

p1 <- ggplot(all_pairs, aes(x = spatial, y = sc)) +
  geom_point(alpha = 0.3) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  geom_smooth(method = "lm", color = "blue") +
  theme_minimal() +
  labs(title = "SC vs Spatial CELLEX Scores",
       subtitle = sprintf("Correlation: %.3f", cor(all_pairs$sc, all_pairs$spatial, use = "complete.obs")),
       x = "Spatial CELLEX", y = "Single-cell CELLEX")
print(p1)

# 2. Pattern distribution
p2 <- ggplot(complementary, aes(x = pattern, fill = pattern)) +
  geom_bar() +
  theme_minimal() +
  labs(title = "Complementary Information Patterns",
       subtitle = "Where each method shows high vs low scores",
       x = "Pattern", y = "Count") +
  theme(legend.position = "none")
print(p2)

# 3. Clarity comparison (when both agree)
clarity_long <- evaluation %>%
  select(cell_type, spatial_clarity, sc_clarity, avg_clarity) %>%
  pivot_longer(cols = c(spatial_clarity, sc_clarity, avg_clarity),
               names_to = "method", values_to = "clarity")

p3 <- ggplot(clarity_long, aes(x = method, y = clarity, fill = method)) +
  geom_boxplot() +
  geom_jitter(alpha = 0.3, width = 0.2) +
  theme_minimal() +
  labs(title = "Prediction Clarity When Both Methods Agree",
       subtitle = "Ratio of top score to second score",
       x = "Method", y = "Clarity (higher = better)") +
  scale_x_discrete(labels = c("avg_clarity" = "Average", 
                              "sc_clarity" = "SC only",
                              "spatial_clarity" = "Spatial only"))
print(p3)

# 4. Disagreement heatmap
disagree_mat <- disagreement_analysis %>%
  filter(!agree) %>%
  mutate(comparison = paste(sc_top_area, "vs", spatial_top_area))

if (nrow(disagree_mat) > 0) {
  p4 <- ggplot(disagree_mat, aes(x = sc_top_area, y = spatial_top_area)) +
    geom_count(aes(size = sc_top_score), alpha = 0.6) +
    theme_minimal() +
    labs(title = "Where SC and Spatial CELLEX Disagree",
         subtitle = "Top predicted area for each method",
         x = "SC CELLEX Top Area", y = "Spatial CELLEX Top Area",
         size = "SC Score") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p4)
}

dev.off()

# SUMMARY REPORT
cat("\n\n=== SUMMARY: Independent Value of SC CELLEX ===\n\n")

cat("1. AGREEMENT:\n")
cat(sprintf("   - Methods agree on top area: %.1f%% of cell types\n", 
            100 * mean(disagreement_analysis$agree)))
cat(sprintf("   - Overall correlation: %.3f\n", 
            cor(all_pairs$sc, all_pairs$spatial, use = "complete.obs")))

cat("\n2. UNIQUE SC PREDICTIONS:\n")
cat(sprintf("   - SC suggests different area: %d cell types\n", 
            sum(disagreement_analysis$sc_adds_value)))
cat(sprintf("   - Validated by spatial stats: %d / %d\n",
            sum(sc_unique$sc_validated, na.rm = TRUE),
            nrow(sc_unique)))

cat("\n3. COMPLEMENTARY INFORMATION:\n")
cat(sprintf("   - Both methods high: %d cases\n", 
            sum(complementary$pattern == "both_high")))
cat(sprintf("   - SC high, Spatial low: %d cases\n", 
            sum(complementary$pattern == "sc_only")))
cat(sprintf("   - Spatial high, SC low: %d cases\n", 
            sum(complementary$pattern == "spatial_only")))

cat("\n4. AVERAGING BENEFIT:\n")
cat(sprintf("   - Spatial clarity: %.2f\n", mean(evaluation$spatial_clarity)))
cat(sprintf("   - SC clarity: %.2f\n", mean(evaluation$sc_clarity)))
cat(sprintf("   - Average clarity: %.2f\n", mean(evaluation$avg_clarity)))
if (mean(evaluation$avg_clarity) > mean(evaluation$spatial_clarity)) {
  cat("   ✓ Averaging IMPROVES clarity over spatial alone\n")
} else {
  cat("   ✗ Averaging does NOT improve clarity over spatial alone\n")
}

# Save results
write.csv(disagreement_analysis, paste0(atlaspath, "sc_spatial_disagreements.csv"), row.names = FALSE)
write.csv(complementary, paste0(atlaspath, "complementary_patterns.csv"), row.names = FALSE)
write.csv(sc_unique, paste0(atlaspath, "sc_unique_predictions.csv"), row.names = FALSE)

cat("\n=== Files Saved ===\n")
cat("1. sc_cellex_independent_value.pdf - Visualizations\n")
cat("2. sc_spatial_disagreements.csv - Where methods disagree\n")
cat("3. complementary_patterns.csv - Complementary information\n")
cat("4. sc_unique_predictions.csv - SC predictions not in spatial\n")

cat("\n✓ Analysis complete!\n")

# Analysis: SUN1 Positivity vs Projection Status and AAV Expression
# Question: Does SUN1+ status correlate with being called a projecting cell?
# Does SUN1+ correlate with AAV transcript levels?

library(dplyr)
library(ggplot2)
library(tidyr)

cat("=== SUN1 Analysis ===\n\n")

# Get cell-level data
cell_data <- data.frame(
  cell_id = rownames(seurat_merged@meta.data),
  sun1 = seurat_merged$sun1,
  projection = seurat_merged$projection,
  projection_spatial_pred = seurat_merged$projection_spatial_pred,
  neurotransmitter = seurat_merged$neurotransmitter,
  injection_filter = seurat_merged$injection_filter,
  run_index = seurat_merged$run_index,
  cell_type = seurat_merged$plot,
  stringsAsFactors = FALSE
)

# Get AAV expression from raw assay (Seurat v5 format)
# Need to check what the AAV genes are called in the assay
cat("Checking for AAV genes in raw assay...\n")
all_features <- rownames(seurat_merged[["raw"]])

# Look for AAV genes more specifically
aav_genes <- grep("p36_NLS_Cre|p147_mCherry_2A_iCre", all_features, value = TRUE)

if (length(aav_genes) == 0) {
  # Try with underscores replaced by hyphens
  aav_genes <- grep("p36-NLS-Cre|p147-mCherry-2A-iCre", all_features, value = TRUE)
}

if (length(aav_genes) == 0) {
  # Show what genes might be AAV-related
  cat("Could not find AAV genes with expected names.\n")
  cat("Searching for any p36 or p147 genes:\n")
  aav_genes <- grep("p36|p147", all_features, value = TRUE, ignore.case = TRUE)
}

cat("AAV-related genes found:\n")
print(aav_genes)
cat("\n")

# Get counts and normalized expression for AAV genes
if (length(aav_genes) >= 2) {
  # Assume first two are the AAV genes
  gene1 <- aav_genes[1]
  gene2 <- aav_genes[2]
  
  cat(sprintf("Using genes: %s, %s\n\n", gene1, gene2))
  
  # Get counts and normalized data from Seurat v5 assay
  # Counts are in: seurat_merged@assays$raw$counts
  # Normalized are in: seurat_merged@assays$raw$scran_normalization
  
  cat("Accessing expression data from raw assay...\n")
  
  # Get feature indices
  all_genes <- rownames(seurat_merged@assays$raw)
  gene1_idx <- which(all_genes == gene1)
  gene2_idx <- which(all_genes == gene2)
  
  # Extract counts
  counts_gene1 <- seurat_merged@assays$raw$counts[gene1_idx, ]
  counts_gene2 <- seurat_merged@assays$raw$counts[gene2_idx, ]
  
  # Extract normalized data
  norm_gene1 <- seurat_merged@assays$raw$scran_normalization[gene1_idx, ]
  norm_gene2 <- seurat_merged@assays$raw$scran_normalization[gene2_idx, ]
  
  cat("✓ Successfully extracted AAV expression data\n\n")
  
  # Add to cell data
  cell_data$aav1_counts <- as.numeric(counts_gene1)
  cell_data$aav2_counts <- as.numeric(counts_gene2)
  cell_data$aav1_norm <- as.numeric(norm_gene1)
  cell_data$aav2_norm <- as.numeric(norm_gene2)
  
  # Total AAV
  cell_data$total_aav_counts <- cell_data$aav1_counts + cell_data$aav2_counts
  cell_data$total_aav_norm <- cell_data$aav1_norm + cell_data$aav2_norm
  
} else {
  cat("Warning: Could not find both AAV genes. Proceeding without AAV expression data.\n")
  cell_data$aav1_counts <- 0
  cell_data$aav2_counts <- 0
  cell_data$aav1_norm <- 0
  cell_data$aav2_norm <- 0
  cell_data$total_aav_counts <- 0
  cell_data$total_aav_norm <- 0
}

# Summary statistics
cat("=== Overall Statistics ===\n")
cat(sprintf("Total cells: %d\n", nrow(cell_data)))
cat(sprintf("SUN1+ cells: %d (%.1f%%)\n", 
            sum(cell_data$sun1, na.rm = TRUE),
            100 * mean(cell_data$sun1, na.rm = TRUE)))
cat(sprintf("SUN1- cells: %d (%.1f%%)\n", 
            sum(!cell_data$sun1, na.rm = TRUE),
            100 * mean(!cell_data$sun1, na.rm = TRUE)))

# Filter for neurons at projection sites (exclude glia and injection sites)
neuron_data <- cell_data %>%
  filter(!is.na(neurotransmitter) & neurotransmitter != "Glia",
         !is.na(injection_filter) & injection_filter == "projection")

cat(sprintf("\nNeurons at projection sites: %d\n", nrow(neuron_data)))
cat(sprintf("SUN1+ neurons: %d (%.1f%%)\n",
            sum(neuron_data$sun1, na.rm = TRUE),
            100 * mean(neuron_data$sun1, na.rm = TRUE)))

# ANALYSIS 1: SUN1 vs Ground Truth Projection Status
cat("\n=== Analysis 1: SUN1 vs Ground Truth Projection ===\n")

projection_analysis <- neuron_data %>%
  filter(!is.na(projection)) %>%
  mutate(
    has_projection = projection != "negative"
  )

sun1_by_projection <- projection_analysis %>%
  group_by(has_projection) %>%
  summarize(
    n = n(),
    n_sun1_pos = sum(sun1, na.rm = TRUE),
    pct_sun1_pos = 100 * mean(sun1, na.rm = TRUE),
    .groups = "drop"
  )

cat("\nSUN1+ rate by projection status (ground truth):\n")
print(sun1_by_projection)

# Chi-square test
if (nrow(projection_analysis) > 0) {
  contingency <- table(projection_analysis$sun1, projection_analysis$has_projection)
  cat("\nContingency table:\n")
  print(contingency)
  
  if (all(dim(contingency) == c(2, 2))) {
    chi_test <- chisq.test(contingency)
    cat(sprintf("\nChi-square test: p = %.2e\n", chi_test$p.value))
    cat(sprintf("Interpretation: %s\n", 
                ifelse(chi_test$p.value < 0.001, 
                       "Highly significant association",
                       ifelse(chi_test$p.value < 0.05,
                              "Significant association",
                              "No significant association"))))
  }
}

# ANALYSIS 2: SUN1 vs Spatial Stats Prediction
cat("\n=== Analysis 2: SUN1 vs Spatial Stats Prediction ===\n")

spatial_analysis <- neuron_data %>%
  filter(!is.na(projection_spatial_pred)) %>%
  mutate(
    has_spatial_pred = projection_spatial_pred != "negative"
  )

sun1_by_spatial <- spatial_analysis %>%
  group_by(has_spatial_pred) %>%
  summarize(
    n = n(),
    n_sun1_pos = sum(sun1, na.rm = TRUE),
    pct_sun1_pos = 100 * mean(sun1, na.rm = TRUE),
    .groups = "drop"
  )

cat("\nSUN1+ rate by spatial prediction:\n")
print(sun1_by_spatial)

# ANALYSIS 3: SUN1 vs AAV Expression
cat("\n=== Analysis 3: SUN1 vs AAV Expression ===\n")

# Compare AAV expression in SUN1+ vs SUN1- cells
aav_by_sun1 <- neuron_data %>%
  group_by(sun1) %>%
  summarize(
    n = n(),
    mean_aav1_counts = mean(aav1_counts, na.rm = TRUE),
    median_aav1_counts = median(aav1_counts, na.rm = TRUE),
    mean_aav2_counts = mean(aav2_counts, na.rm = TRUE),
    median_aav2_counts = median(aav2_counts, na.rm = TRUE),
    mean_total_counts = mean(total_aav_counts, na.rm = TRUE),
    median_total_counts = median(total_aav_counts, na.rm = TRUE),
    pct_with_any_aav = 100 * mean(total_aav_counts > 0, na.rm = TRUE),
    .groups = "drop"
  )

cat("\nAAV expression by SUN1 status:\n")
print(aav_by_sun1)

# Statistical tests
if (length(aav_genes) >= 2) {
  cat("\nWilcoxon tests (SUN1+ vs SUN1-):\n")
  
  sun1_pos <- neuron_data$sun1 == TRUE
  sun1_neg <- neuron_data$sun1 == FALSE
  
  # Test AAV1
  if (sum(neuron_data$aav1_counts[sun1_pos] > 0) > 0 & 
      sum(neuron_data$aav1_counts[sun1_neg] > 0) > 0) {
    test_aav1 <- wilcox.test(neuron_data$aav1_counts[sun1_pos],
                             neuron_data$aav1_counts[sun1_neg])
    cat(sprintf("  AAV1 (%s): p = %.2e\n", gene1, test_aav1$p.value))
  }
  
  # Test AAV2
  if (sum(neuron_data$aav2_counts[sun1_pos] > 0) > 0 & 
      sum(neuron_data$aav2_counts[sun1_neg] > 0) > 0) {
    test_aav2 <- wilcox.test(neuron_data$aav2_counts[sun1_pos],
                             neuron_data$aav2_counts[sun1_neg])
    cat(sprintf("  AAV2 (%s): p = %.2e\n", gene2, test_aav2$p.value))
  }
  
  # Test total
  if (sum(neuron_data$total_aav_counts[sun1_pos] > 0) > 0 & 
      sum(neuron_data$total_aav_counts[sun1_neg] > 0) > 0) {
    test_total <- wilcox.test(neuron_data$total_aav_counts[sun1_pos],
                              neuron_data$total_aav_counts[sun1_neg])
    cat(sprintf("  Total AAV: p = %.2e\n", test_total$p.value))
  }
}

# ANALYSIS 4: By Run (different AAV-area mappings)
cat("\n=== Analysis 4: SUN1 by Run ===\n")

sun1_by_run <- neuron_data %>%
  group_by(run_index) %>%
  summarize(
    n = n(),
    n_sun1_pos = sum(sun1, na.rm = TRUE),
    pct_sun1_pos = 100 * mean(sun1, na.rm = TRUE),
    mean_aav_counts = mean(total_aav_counts, na.rm = TRUE),
    .groups = "drop"
  )

cat("\nSUN1+ rate by experimental run:\n")
print(sun1_by_run)

# VISUALIZATIONS
pdf(paste0(atlaspath, "sun1_analysis.pdf"), width = 16, height = 20)

# 1. SUN1 vs Projection Status - improved
par(mfrow = c(2, 3), mar = c(5, 5, 4, 2))

if (nrow(sun1_by_projection) > 0) {
  # Ground truth
  bp1 <- barplot(sun1_by_projection$pct_sun1_pos,
          names.arg = c("Non-projecting", "Projecting"),
          col = c("gray80", "steelblue"),
          ylab = "% SUN1+", ylim = c(0, 30),
          main = "SUN1+ Rate by Projection Status\n(Ground Truth)",
          cex.names = 1.2, cex.axis = 1.2, cex.lab = 1.2)
  text(bp1, sun1_by_projection$pct_sun1_pos + 1.5, 
       sprintf("%.1f%%", sun1_by_projection$pct_sun1_pos),
       cex = 1.3, font = 2)
  text(bp1, sun1_by_projection$pct_sun1_pos + 4,
       sprintf("n=%d", sun1_by_projection$n),
       cex = 1, col = "gray40")
  grid()
}

# 2. SUN1 vs Spatial Prediction
if (nrow(sun1_by_spatial) > 0) {
  bp2 <- barplot(sun1_by_spatial$pct_sun1_pos,
          names.arg = c("Non-projecting", "Projecting"),
          col = c("gray80", "darkgreen"),
          ylab = "% SUN1+", ylim = c(0, 20),
          main = "SUN1+ Rate by Spatial Prediction",
          cex.names = 1.2, cex.axis = 1.2, cex.lab = 1.2)
  text(bp2, sun1_by_spatial$pct_sun1_pos + 1, 
       sprintf("%.1f%%", sun1_by_spatial$pct_sun1_pos),
       cex = 1.3, font = 2)
  text(bp2, sun1_by_spatial$pct_sun1_pos + 3,
       sprintf("n=%d", sun1_by_spatial$n),
       cex = 1, col = "gray40")
  grid()
}

# 3. Fold enrichment comparison
if (nrow(sun1_by_projection) > 0) {
  fold_enrichment <- sun1_by_projection$pct_sun1_pos[2] / sun1_by_projection$pct_sun1_pos[1]
  barplot(c(1, fold_enrichment),
          names.arg = c("Non-projecting\n(baseline)", "Projecting"),
          col = c("gray80", "coral"),
          ylab = "Fold Enrichment of SUN1+",
          main = sprintf("SUN1 Enrichment\n(%.1f-fold)", fold_enrichment),
          cex.names = 1.2, cex.axis = 1.2, cex.lab = 1.2,
          ylim = c(0, max(fold_enrichment * 1.2, 2)))
  abline(h = 1, lty = 2, col = "red", lwd = 2)
  text(2, fold_enrichment + 0.5, sprintf("%.1f×", fold_enrichment),
       cex = 1.5, font = 2, col = "coral")
  grid()
}

# 4. AAV expression by SUN1 - log scale boxplots
if (length(aav_genes) >= 2) {
  # AAV1 - focused on positive values
  aav1_positive <- neuron_data$aav1_counts[neuron_data$aav1_counts > 0]
  sun1_positive <- neuron_data$sun1[neuron_data$aav1_counts > 0]
  
  boxplot(log10(aav1_positive + 1) ~ sun1_positive,
          names = c("SUN1-", "SUN1+"),
          col = c("gray70", "coral"),
          ylab = "log10(Counts + 1)", 
          main = sprintf("AAV1 Expression\n(%s)\n(cells with AAV1 > 0)", gene1),
          outline = FALSE, cex.names = 1.2, cex.axis = 1.2, cex.lab = 1.2)
  grid()
  
  # AAV2
  aav2_positive <- neuron_data$aav2_counts[neuron_data$aav2_counts > 0]
  sun1_positive2 <- neuron_data$sun1[neuron_data$aav2_counts > 0]
  
  if (length(aav2_positive) > 10) {
    boxplot(log10(aav2_positive + 1) ~ sun1_positive2,
            names = c("SUN1-", "SUN1+"),
            col = c("gray70", "coral"),
            ylab = "log10(Counts + 1)", 
            main = sprintf("AAV2 Expression\n(%s)\n(cells with AAV2 > 0)", gene2),
            outline = FALSE, cex.names = 1.2, cex.axis = 1.2, cex.lab = 1.2)
    grid()
  }
  
  # Detection rate comparison
  detection_rates <- neuron_data %>%
    group_by(sun1) %>%
    summarize(
      pct_aav1_detected = 100 * mean(aav1_counts > 0),
      pct_aav2_detected = 100 * mean(aav2_counts > 0),
      pct_any_detected = 100 * mean(total_aav_counts > 0),
      .groups = "drop"
    )
  
  bp_detect <- barplot(as.matrix(t(detection_rates[, c("pct_aav1_detected", "pct_any_detected")])),
          beside = TRUE, 
          col = c("gray70", "coral"),
          names.arg = c("AAV1", "Any AAV"),
          ylab = "% Cells with Detection",
          main = "AAV Detection Rate by SUN1 Status",
          ylim = c(0, max(detection_rates$pct_any_detected) * 1.2),
          cex.names = 1.2, cex.axis = 1.2, cex.lab = 1.2,
          legend.text = c("SUN1-", "SUN1+"),
          args.legend = list(x = "topright", cex = 1.2))
  
  # Add values on bars
  for (i in 1:ncol(bp_detect)) {
    text(bp_detect[, i], 
         as.matrix(t(detection_rates[, c("pct_aav1_detected", "pct_any_detected")]))[, i] + 3,
         sprintf("%.1f%%", as.matrix(t(detection_rates[, c("pct_aav1_detected", "pct_any_detected")]))[, i]),
         cex = 1.1)
  }
  grid()
}

# 5. AAV expression distributions - better visualization
par(mfrow = c(2, 2), mar = c(5, 5, 4, 2))

if (length(aav_genes) >= 2) {
  # AAV1 violin-style with better bins
  sun1_neg_aav1 <- log10(neuron_data$aav1_counts[neuron_data$sun1 == FALSE] + 1)
  sun1_pos_aav1 <- log10(neuron_data$aav1_counts[neuron_data$sun1 == TRUE] + 1)
  
  h1 <- hist(sun1_neg_aav1, breaks = 100, plot = FALSE)
  h2 <- hist(sun1_pos_aav1, breaks = 100, plot = FALSE)
  
  plot(h1, col = rgb(0.5, 0.5, 0.5, 0.5), border = NA,
       xlab = "log10(Counts + 1)", main = "AAV1 Distribution",
       xlim = c(0, max(c(sun1_neg_aav1, sun1_pos_aav1))),
       cex.lab = 1.2, cex.axis = 1.2, cex.main = 1.3)
  plot(h2, col = rgb(1, 0.5, 0, 0.6), border = NA, add = TRUE)
  legend("topright", c("SUN1-", "SUN1+"), 
         fill = c(rgb(0.5, 0.5, 0.5, 0.5), rgb(1, 0.5, 0, 0.6)),
         cex = 1.2, bty = "n")
  grid()
  
  # AAV total
  sun1_neg_total <- log10(neuron_data$total_aav_counts[neuron_data$sun1 == FALSE] + 1)
  sun1_pos_total <- log10(neuron_data$total_aav_counts[neuron_data$sun1 == TRUE] + 1)
  
  h3 <- hist(sun1_neg_total, breaks = 100, plot = FALSE)
  h4 <- hist(sun1_pos_total, breaks = 100, plot = FALSE)
  
  plot(h3, col = rgb(0.5, 0.5, 0.5, 0.5), border = NA,
       xlab = "log10(Total Counts + 1)", main = "Total AAV Distribution",
       xlim = c(0, max(c(sun1_neg_total, sun1_pos_total))),
       cex.lab = 1.2, cex.axis = 1.2, cex.main = 1.3)
  plot(h4, col = rgb(1, 0.5, 0, 0.6), border = NA, add = TRUE)
  legend("topright", c("SUN1-", "SUN1+"), 
         fill = c(rgb(0.5, 0.5, 0.5, 0.5), rgb(1, 0.5, 0, 0.6)),
         cex = 1.2, bty = "n")
  grid()
  
  # Mean expression comparison
  mean_expr <- neuron_data %>%
    group_by(sun1) %>%
    summarize(
      mean_aav1 = mean(aav1_counts),
      mean_total = mean(total_aav_counts),
      .groups = "drop"
    )
  
  bp_mean <- barplot(as.matrix(t(mean_expr[, c("mean_aav1", "mean_total")])),
          beside = TRUE,
          col = c("gray70", "coral"),
          names.arg = c("AAV1", "Total AAV"),
          ylab = "Mean Counts",
          main = "Mean AAV Expression",
          cex.names = 1.2, cex.axis = 1.2, cex.lab = 1.2,
          legend.text = c("SUN1-", "SUN1+"),
          args.legend = list(x = "topright", cex = 1.2))
  
  # Add values
  for (i in 1:ncol(bp_mean)) {
    vals <- as.matrix(t(mean_expr[, c("mean_aav1", "mean_total")]))[, i]
    text(bp_mean[, i], vals + max(vals) * 0.05,
         sprintf("%.1f", vals), cex = 1.1)
  }
  grid()
  
  # Scatter with density coloring
  sun1_colors <- ifelse(neuron_data$sun1, 
                       rgb(1, 0.3, 0, 0.4), 
                       rgb(0.3, 0.3, 0.3, 0.1))
  
  plot(jitter(log10(neuron_data$aav1_counts + 1), amount = 0.05),
       jitter(log10(neuron_data$aav2_counts + 1), amount = 0.05),
       pch = 16, cex = 0.4, col = sun1_colors,
       xlab = "log10(AAV1 + 1)", ylab = "log10(AAV2 + 1)",
       main = "AAV1 vs AAV2 Expression",
       cex.lab = 1.2, cex.axis = 1.2, cex.main = 1.3)
  legend("topright", c("SUN1-", "SUN1+"),
         pch = 16, col = c(rgb(0.3, 0.3, 0.3, 0.5), rgb(1, 0.3, 0, 0.8)),
         cex = 1.2, bty = "n")
  grid()
}

# 6. SUN1 by cell type - improved
par(mfrow = c(1, 1), mar = c(5, 12, 4, 2))

if (nrow(sun1_by_celltype) > 0) {
  # Color by percentage
  colors <- colorRampPalette(c("lightblue", "darkblue"))(nrow(sun1_by_celltype))
  
  bp_ct <- barplot(sun1_by_celltype$pct_sun1,
          names.arg = sun1_by_celltype$cell_type,
          horiz = TRUE, las = 1,
          col = colors,
          xlab = "% SUN1+",
          main = "SUN1+ Rate by Cell Type\n(Top 20, n≥50)",
          xlim = c(0, max(sun1_by_celltype$pct_sun1) * 1.2),
          cex.names = 1, cex.axis = 1.2, cex.lab = 1.2, cex.main = 1.3)
  
  # Add percentage labels
  text(sun1_by_celltype$pct_sun1 + 0.5, bp_ct,
       sprintf("%.1f%%", sun1_by_celltype$pct_sun1),
       pos = 4, cex = 1)
  
  # Add cell counts
  text(0.2, bp_ct,
       sprintf("n=%d", sun1_by_celltype$n),
       pos = 4, cex = 0.8, col = "white", font = 2)
  grid()
}

dev.off()

# Save results
sun1_summary <- data.frame(
  analysis = c("Overall", "Neurons_projection_sites", "Ground_truth_projecting", 
               "Ground_truth_negative", "Spatial_pred_projecting", "Spatial_pred_negative"),
  n_cells = c(
    nrow(cell_data),
    nrow(neuron_data),
    sum(projection_analysis$has_projection, na.rm = TRUE),
    sum(!projection_analysis$has_projection, na.rm = TRUE),
    sum(spatial_analysis$has_spatial_pred, na.rm = TRUE),
    sum(!spatial_analysis$has_spatial_pred, na.rm = TRUE)
  ),
  n_sun1_pos = c(
    sum(cell_data$sun1, na.rm = TRUE),
    sum(neuron_data$sun1, na.rm = TRUE),
    sum(projection_analysis$sun1[projection_analysis$has_projection], na.rm = TRUE),
    sum(projection_analysis$sun1[!projection_analysis$has_projection], na.rm = TRUE),
    sum(spatial_analysis$sun1[spatial_analysis$has_spatial_pred], na.rm = TRUE),
    sum(spatial_analysis$sun1[!spatial_analysis$has_spatial_pred], na.rm = TRUE)
  ),
  pct_sun1_pos = c(
    100 * mean(cell_data$sun1, na.rm = TRUE),
    100 * mean(neuron_data$sun1, na.rm = TRUE),
    100 * mean(projection_analysis$sun1[projection_analysis$has_projection], na.rm = TRUE),
    100 * mean(projection_analysis$sun1[!projection_analysis$has_projection], na.rm = TRUE),
    100 * mean(spatial_analysis$sun1[spatial_analysis$has_spatial_pred], na.rm = TRUE),
    100 * mean(spatial_analysis$sun1[!spatial_analysis$has_spatial_pred], na.rm = TRUE)
  )
)

write.csv(sun1_summary, paste0(atlaspath, "sun1_analysis_summary.csv"), row.names = FALSE)
write.csv(aav_by_sun1, paste0(atlaspath, "aav_expression_by_sun1.csv"), row.names = FALSE)
write.csv(sun1_by_celltype, paste0(atlaspath, "sun1_by_celltype.csv"), row.names = FALSE)

cat("\n=== Files Saved ===\n")
cat("1. sun1_analysis.pdf - Comprehensive visualizations\n")
cat("2. sun1_analysis_summary.csv - Summary statistics\n")
cat("3. aav_expression_by_sun1.csv - AAV levels by SUN1 status\n")
cat("4. sun1_by_celltype.csv - SUN1+ rates by cell type\n")

cat("\n✓ SUN1 analysis complete!\n")

unique(seurat_merged$region)

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

names(colors_detailed) <- cells_detailed
names(colors_detailed_light) <- cells_detailed
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

seurat_merged$detailed_celltype <- seurat_merged$detailed_celltype_annotation
seurat_merged@meta.data[["detailed_celltype"]] <- ifelse(seurat_merged@meta.data[["score_detailed_celltype_annotation"]] > 0.4, 
                                                         seurat_merged@meta.data[["detailed_celltype"]], "no_label")
glia <- c('04_Crnde_Lmx1a', '27_Lmx1a_Ttc6', '32_Gucy2c_Ttc6', '36_Ebf2_Crhbp')
seurat_merged$plot <- seurat_merged$detailed_celltype
seurat_merged@meta.data[["plot"]] <- ifelse(seurat_merged@meta.data[["plot"]] %in% glia, seurat_merged@meta.data[["plot"]], "no_label")

plot <- ImageDimPlot(seurat_merged, 
             border.color = NA, 
             boundaries = "segmentation", 
             group.by = "plot", 
             cols = colors_detailed_light,
             fov = "S3C2",
             dark.background = FALSE,
            )+ theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      legend.key.size = unit(0.5, "lines"),
                      legend.text = element_text(size = 9), 
                      plot.background = element_rect(fill = "white"),
                      panel.background = element_rect(fill = "white"), 
                      plot.margin = unit(c(0, 0, 0, 0), "cm"))
show(plot)
resolve_save(plot, paste0(atlaspath, "lmx1a_bright.pdf"))

plot <- ImageDimPlot(seurat_merged, 
             border.color = NA, 
             boundaries = "segmentation", 
             group.by = "plot", 
             cols = colors_detailed_light,
             fov = "S3C2",
             molecules = c("Lmx1b"),
             mols.cols = c("yellow"),
             mols.size = 1,
             mols.alpha = 0.25,
             nmols=1000,
             dark.background = FALSE,
            )+ theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      legend.key.size = unit(0.5, "lines"),
                      legend.text = element_text(size = 9), 
                      plot.background = element_rect(fill = "white"),
                      panel.background = element_rect(fill = "white"), 
                      plot.margin = unit(c(0, 0, 0, 0), "cm"))
show(plot)
resolve_save(plot, paste0(atlaspath, "lmx1a-b_bright.pdf"))

table(seurat_merged$old.ident, seurat_merged$bregma)

seurat_merged$detailed_celltype <- seurat_merged$detailed_celltype_annotation
seurat_merged@meta.data[["detailed_celltype"]] <- ifelse(seurat_merged@meta.data[["score_detailed_celltype_annotation"]] > 0.4, 
                                                         seurat_merged@meta.data[["detailed_celltype"]], "no_label")
glia <- c('04_Crnde_Lmx1a', '27_Lmx1a_Ttc6', '32_Gucy2c_Ttc6', '36_Ebf2_Crhbp')
seurat_merged$plot <- seurat_merged$detailed_celltype
seurat_merged@meta.data[["plot"]] <- ifelse(seurat_merged@meta.data[["plot"]] %in% glia, "Lmx1a+", "Lmx1a-")

plot <- ImageDimPlot(seurat_merged, 
             border.color = NA, 
             boundaries = "segmentation", 
             group.by = "plot", 
             cols = c("#d4d4d4", "#452C2C"),
             fov = "S1C1",
             dark.background = FALSE,
            )+ theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      legend.key.size = unit(0.5, "lines"),
                      legend.text = element_text(size = 9), 
                      plot.background = element_rect(fill = "white"),
                      panel.background = element_rect(fill = "white"), 
                      plot.margin = unit(c(0, 0, 0, 0), "cm"))
show(plot)
resolve_save(plot, paste0(atlaspath, "lmx1a_S1C1.pdf"))

plot <- ImageDimPlot(seurat_merged, 
             border.color = NA, 
             boundaries = "segmentation", 
             group.by = "plot", 
             cols = c("#d4d4d4", "#452C2C"),
             fov = "S2C2",
             dark.background = FALSE,
            )+ theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      legend.key.size = unit(0.5, "lines"),
                      legend.text = element_text(size = 9), 
                      plot.background = element_rect(fill = "white"),
                      panel.background = element_rect(fill = "white"), 
                      plot.margin = unit(c(0, 0, 0, 0), "cm"))
show(plot)
resolve_save(plot, paste0(atlaspath, "lmx1a_S2C2.pdf"))

plot <- ImageDimPlot(seurat_merged, 
             border.color = NA, 
             boundaries = "segmentation", 
             group.by = "plot", 
             cols = c("#d4d4d4", "#452C2C"),
             fov = "S3C2",
             dark.background = FALSE,
            )+ theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      legend.key.size = unit(0.5, "lines"),
                      legend.text = element_text(size = 9), 
                      plot.background = element_rect(fill = "white"),
                      panel.background = element_rect(fill = "white"), 
                      plot.margin = unit(c(0, 0, 0, 0), "cm"))
show(plot)
resolve_save(plot, paste0(atlaspath, "lmx1a_S3C2.pdf"))

plot <- ImageDimPlot(seurat_merged, 
             border.color = NA, 
             boundaries = "segmentation", 
             group.by = "plot", 
             cols = c("#452C2C", "#d4d4d4"),
             fov = "S2C1",
             dark.background = FALSE,
            )+ theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      legend.key.size = unit(0.5, "lines"),
                      legend.text = element_text(size = 9), 
                      plot.background = element_rect(fill = "white"),
                      panel.background = element_rect(fill = "white"), 
                      plot.margin = unit(c(0, 0, 0, 0), "cm"))
show(plot)
resolve_save(plot, paste0(atlaspath, "lmx1a_S2C1.pdf"))

plot <- ImageDimPlot(seurat_merged, 
             border.color = NA, 
             boundaries = "segmentation", 
             group.by = "plot", 
             cols = c("#d4d4d4", "#452C2C"),
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
resolve_save(plot, paste0(atlaspath, "lmx1a_S3C1.pdf"))

library(ggplot2)
library(dplyr)

# Define cell types of interest
cell_types <- c("04_Crnde_Lmx1a", "27_Lmx1a_Ttc6", "32_Gucy2c_Ttc6", "36_Ebf2_Crhbp")

# Build summary dataframe from the seurat object
fov_summary <- seurat_merged@meta.data %>%
  filter(plot %in% cell_types) %>%
  group_by(old.ident) %>%
  summarise(total_cells = n(), .groups = "drop") %>%
  arrange(desc(total_cells)) %>%
  mutate(old.ident = factor(old.ident, levels = old.ident))  # preserve rank order

# Plot
ggplot(fov_summary, aes(x = old.ident, y = total_cells, fill = total_cells)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = total_cells), vjust = -0.4, size = 3) +
  scale_fill_gradient(low = "#a8d8ea", high = "#d62728") +
  labs(
    title = "FOVs Ranked by Total Cells Across 4 Cell Types of Interest",
    subtitle = "Cell types: 04_Crnde_Lmx1a | 27_Lmx1a_Ttc6 | 32_Gucy2c_Ttc6 | 36_Ebf2_Crhbp",
    x = "FOV (old.ident)",
    y = "Total Cell Count",
    fill = "Total Cells"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    legend.position = "none"
  )

library(ggplot2)
library(dplyr)

# Define cell types of interest
cell_types <- c("04_Crnde_Lmx1a", "27_Lmx1a_Ttc6", "32_Gucy2c_Ttc6", "36_Ebf2_Crhbp")

# Build summary dataframe grouped by bregma, filtered to coronal sections only
fov_summary <- seurat_merged@meta.data %>%
  filter(section_pane == "coronal",
         plot %in% cell_types) %>%
  group_by(old.ident, bregma) %>%
  summarise(total_cells = n(), .groups = "drop") %>%
  arrange(desc(total_cells)) %>%
  mutate(old.ident = factor(old.ident, levels = old.ident))

# Plot
ggplot(fov_summary, aes(x = old.ident, y = total_cells, fill = total_cells)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = total_cells), vjust = -0.4, size = 3) +
  scale_fill_gradient(low = "#a8d8ea", high = "#d62728") +
  scale_x_discrete(labels = setNames(fov_summary$bregma, fov_summary$old.ident)) +
  labs(
    title = "Coronal FOVs Ranked by Total Cells Across 4 Cell Types of Interest",
    subtitle = "Cell types: 04_Crnde_Lmx1a | 27_Lmx1a_Ttc6 | 32_Gucy2c_Ttc6 | 36_Ebf2_Crhbp",
    x = "Bregma Level",
    y = "Total Cell Count"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    legend.position = "none"
  )
ggsave(paste0(atlaspath, "lmx1a_bregma.pdf"))

table(seurat_merged$plot, seurat_merged$projection)

library(ggplot2)
library(dplyr)

# Define cell types of interest
cell_types <- c("04_Crnde_Lmx1a", "27_Lmx1a_Ttc6", "32_Gucy2c_Ttc6", "36_Ebf2_Crhbp")

# Extract matched colors directly from the palette
cell_colors <- setNames(
  colors_detailed[match(cell_types, cells_detailed)],
  cell_types
)

# Build summary dataframe
projection_summary <- seurat_merged@meta.data %>%
  filter(plot %in% cell_types,
         !projection %in% c("negative", "no_label"),
         injection_filter == "projection") %>%
  group_by(projection, plot) %>%
  summarise(total_cells = n(), .groups = "drop")

# Rank projections by total cells across all 4 types
projection_order <- projection_summary %>%
  group_by(projection) %>%
  summarise(grand_total = sum(total_cells)) %>%
  arrange(desc(grand_total)) %>%
  pull(projection)

projection_summary <- projection_summary %>%
  mutate(projection = factor(projection, levels = projection_order),
         plot = factor(plot, levels = cell_types))

# Plot
ggplot(projection_summary, aes(x = projection, y = total_cells, fill = plot)) +
  geom_bar(stat = "identity") +
  geom_text(
    data = projection_summary %>%
      group_by(projection) %>%
      summarise(grand_total = sum(total_cells)),
    aes(x = projection, y = grand_total, label = grand_total),
    inherit.aes = FALSE,
    vjust = -0.4, size = 3
  ) +
  scale_fill_manual(values = cell_colors, name = "Cell Type") +
  labs(
    title = "Projections Ranked by Total Cells Across 4 Cell Types of Interest",
    subtitle = "Excluding injection regions & 'negative' projection calls | All sections",
    x = "Projection Target",
    y = "Total Cell Count"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    legend.position = "right"
  )
ggsave(paste0(atlaspath, "lmx1a_projections.pdf"))

library(ggplot2)
library(dplyr)

# Define cell types of interest
cell_types <- c("04_Crnde_Lmx1a", "27_Lmx1a_Ttc6", "32_Gucy2c_Ttc6", "36_Ebf2_Crhbp")

# Extract matched colors directly from the palette
cell_colors <- setNames(
  colors_detailed[match(cell_types, cells_detailed)],
  cell_types
)

# Build summary dataframe grouped by bregma + plot, filtered to coronal sections only
fov_summary <- seurat_merged@meta.data %>%
  filter(section_pane == "coronal",
         plot %in% cell_types) %>%
  group_by(old.ident, bregma, plot) %>%
  summarise(total_cells = n(), .groups = "drop")

# Rank FOVs by grand total across all 4 cell types
fov_order <- fov_summary %>%
  group_by(old.ident, bregma) %>%
  summarise(grand_total = sum(total_cells), .groups = "drop") %>%
  arrange(desc(grand_total))

fov_summary <- fov_summary %>%
  mutate(old.ident = factor(old.ident, levels = fov_order$old.ident),
         plot = factor(plot, levels = cell_types))

# Plot
ggplot(fov_summary, aes(x = old.ident, y = total_cells, fill = plot)) +
  geom_bar(stat = "identity") +
  geom_text(
    data = fov_order,
    aes(x = old.ident, y = grand_total, label = grand_total),
    inherit.aes = FALSE,
    vjust = -0.4, size = 3
  ) +
  scale_fill_manual(values = cell_colors, name = "Cell Type") +
  scale_x_discrete(labels = setNames(fov_order$bregma, fov_order$old.ident)) +
  labs(
    title = "Coronal FOVs Ranked by Total Cells Across 4 Cell Types of Interest",
    subtitle = "Cell types: 04_Crnde_Lmx1a | 27_Lmx1a_Ttc6 | 32_Gucy2c_Ttc6 | 36_Ebf2_Crhbp",
    x = "Bregma Level",
    y = "Total Cell Count"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    legend.position = "right"
  )

ggsave(paste0(atlaspath, "lmx1a_bregma.pdf"))

library(ggplot2)
library(dplyr)

# Define cell types of interest
cell_types <- c("04_Crnde_Lmx1a", "27_Lmx1a_Ttc6", "32_Gucy2c_Ttc6", "36_Ebf2_Crhbp")

# Extract matched colors directly from the palette
cell_colors <- setNames(
  colors_detailed[match(cell_types, cells_detailed)],
  cell_types
)

# Build summary dataframe grouped by bregma + plot, filtered to coronal sections only
fov_summary <- seurat_merged@meta.data %>%
  filter(section_pane == "coronal",
         plot %in% cell_types) %>%
  group_by(old.ident, bregma, plot) %>%
  summarise(total_cells = n(), .groups = "drop")

# Rank FOVs by bregma position: anterior (most positive) → posterior (most negative)
fov_order <- fov_summary %>%
  group_by(old.ident, bregma) %>%
  summarise(grand_total = sum(total_cells), .groups = "drop") %>%
  arrange(desc(as.numeric(bregma)), desc(grand_total))

fov_summary <- fov_summary %>%
  mutate(old.ident = factor(old.ident, levels = fov_order$old.ident),
         plot = factor(plot, levels = cell_types))

# Plot
ggplot(fov_summary, aes(x = old.ident, y = total_cells, fill = plot)) +
  geom_bar(stat = "identity") +
  geom_text(
    data = fov_order,
    aes(x = old.ident, y = grand_total, label = grand_total),
    inherit.aes = FALSE,
    vjust = -0.4, size = 3
  ) +
  scale_fill_manual(values = cell_colors, name = "Cell Type") +
  scale_x_discrete(labels = setNames(fov_order$bregma, fov_order$old.ident)) +
  labs(
    title = "Coronal FOVs Ranked by Total Cells Across 4 Cell Types of Interest",
    subtitle = "Cell types: 04_Crnde_Lmx1a | 27_Lmx1a_Ttc6 | 32_Gucy2c_Ttc6 | 36_Ebf2_Crhbp",
    x = "Bregma Level",
    y = "Total Cell Count"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    legend.position = "right"
  )

ggsave(paste0(atlaspath, "lmx1a_bregma.pdf"))

# Build summary dataframe grouped by bregma + plot, filtered to coronal sections only
fov_summary <- seurat_merged@meta.data %>%
  filter(section_pane == "coronal",
         plot %in% cell_types) %>%
  group_by(old.ident, bregma, plot) %>%
  summarise(total_cells = n(), .groups = "drop")

# Merge duplicate bregma levels by averaging across FOVs
fov_summary <- fov_summary %>%
  group_by(bregma, plot) %>%
  summarise(total_cells = mean(total_cells), .groups = "drop")

# Order bregma anterior → posterior and set factor levels
bregma_levels <- fov_summary %>%
  distinct(bregma) %>%
  arrange(desc(as.numeric(bregma))) %>%
  pull(bregma)

fov_summary <- fov_summary %>%
  mutate(bregma = factor(bregma, levels = bregma_levels),
         plot   = factor(plot, levels = cell_types))

# Compute grand totals AFTER factor is set, so x aligns with bars
fov_order <- fov_summary %>%
  group_by(bregma) %>%
  summarise(grand_total = sum(total_cells), .groups = "drop")

ggplot(fov_summary, aes(x = bregma, y = total_cells, fill = plot)) +
  geom_bar(stat = "identity") +
  geom_text(
    data = fov_order,
    aes(x = bregma, y = grand_total, label = round(grand_total)),
    inherit.aes = FALSE,
    vjust = -0.4, size = 3
  ) +
  scale_fill_manual(values = cell_colors, name = "Cell Type") +
  labs(
    title = "Number of mapped cells for each coronal FOV \n (filtered to 4 cell type of interest)", 
    x = "Bregma Level",
    y = "Total Cell Count"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    legend.position = "right"
  )

ggsave(paste0(atlaspath, "lmx1a_bregma_merged.pdf"))

library(ggplot2)
library(dplyr)
library(tidyr)

# Build summary - rename projection column to avoid clash with injection_filter values
qc_heatmap <- seurat_merged@meta.data %>%
  filter(plot != "no_label",
         projection != "negative",
         !is.na(injection_filter)) %>%
  group_by(plot, projection, injection_filter) %>%
  summarise(total_cells = n(), .groups = "drop") %>%
  rename(proj_target = projection) %>%  # avoid name clash with "projection" level in injection_filter
  pivot_wider(names_from = injection_filter,
              values_from = total_cells,
              values_fill = 0)

# Compute Pearson residuals
qc_residuals <- qc_heatmap %>%
  mutate(
    total        = injection + projection,
    prop_inj     = sum(injection) / (sum(injection) + sum(projection)),
    expected_inj = total * prop_inj,
    residual     = (injection - expected_inj) / sqrt(expected_inj)
  ) %>%
  group_by(plot) %>%
  mutate(mean_residual = mean(residual)) %>%
  ungroup() %>%
  mutate(plot       = reorder(plot, mean_residual),
         proj_target = reorder(proj_target, residual, FUN = mean))

# Symmetric color scale limits
lim <- max(abs(qc_residuals$residual), na.rm = TRUE)

# Plot
ggplot(qc_residuals, aes(x = proj_target, y = plot, fill = residual)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(
    low      = "#4e79a7",
    mid      = "white",
    high     = "#e15759",
    midpoint = 0,
    limits   = c(-lim, lim),
    name     = "Pearson\nResidual"
  ) +
  labs(
    title    = "Projection Call Residuals: Injection vs. Projection Sites",
    subtitle = "Red = over-represented in injection sites | Blue = over-represented in projection sites",
    x        = "Projection Target",
    y        = "Cell Type"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x   = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y   = element_text(size = 8),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    legend.position = "right"
  )

ggsave(paste0(atlaspath, "lmx1a_projection_residuals_heatmap.pdf"), width = 10, height = 14)


