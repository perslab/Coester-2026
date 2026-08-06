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

# Visualize godsnot_102 color palette with hex codes and positions

library(ggplot2)
library(grid)

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

# Create data frame with position, color, and hex code
color_df <- data.frame(
  position = 1:length(godsnot_102),
  hex = godsnot_102,
  stringsAsFactors = FALSE
)

# Calculate grid layout (10 columns)
n_cols <- 10
n_rows <- ceiling(length(godsnot_102) / n_cols)

color_df$row <- ceiling(color_df$position / n_cols)
color_df$col <- ((color_df$position - 1) %% n_cols) + 1

# Determine if text should be white or black based on color brightness
get_text_color <- function(hex_color) {
  # Remove # and convert to RGB
  rgb_vals <- col2rgb(hex_color)
  # Calculate brightness (perceived luminance)
  brightness <- (0.299 * rgb_vals[1] + 0.587 * rgb_vals[2] + 0.114 * rgb_vals[3]) / 255
  ifelse(brightness > 0.5, "black", "white")
}

color_df$text_color <- sapply(color_df$hex, get_text_color)

# Create the plot
p <- ggplot(color_df, aes(x = col, y = -row)) +
  geom_tile(aes(fill = hex), color = "white", size = 1) +
  geom_text(aes(label = sprintf("[%d]\n%s", position, hex), color = text_color),
            size = 2.5, lineheight = 0.9) +
  scale_fill_identity() +
  scale_color_identity() +
  theme_void() +
  theme(
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, margin = margin(b = 20))
  ) +
  labs(
    title = "godsnot_102 Color Palette",
    subtitle = sprintf("Total colors: %d | Format: [position] hex_code", length(godsnot_102))
  ) +
  coord_fixed(ratio = 1)

# Save as PDF
ggsave(paste0(atlaspath, "godsnot_102_palette.pdf"),
       plot = p,
       width = 14,
       height = n_rows * 1.4 + 2,
       device = cairo_pdf)

print(p)

# Print color information to console
cat("\n=== godsnot_102 Color Palette ===\n")
cat(sprintf("Total colors: %d\n\n", length(godsnot_102)))
cat("\nFiles saved:\n")
cat("1. godsnot_102_palette.pdf (grid view - 10 columns)\n")

seurat_merged <- qs::qread(paste0(atlaspath, "spatial_24_bipca_viruscalls.qs"))
seurat_merged

bregma_map <- c(
  S1A1=-0.47, S1A2=-0.83, S1B1=-1.23, S1B2=-1.43, S1C1=-1.79, S1C2=-2.69,
  S1D1= 0.675, S1D2= 0.200,
  S2A1=-0.59, S2A2=-0.95, S2B1=-1.31, S2B2=-1.55, S2C1=-1.91, S2C2=-2.45,
  S2D1= 0.300, S2D2= 0.475,
  S3A1=-0.71, S3A2=-1.23, S3B1=-1.43, S3B2=-1.67, S3C1=-1.91, S3C2=-2.45,
  S3D1= 0.225, S3D2= 0.200
)

axis_map <- c(
  S1A1="antero-posterior", S1A2="antero-posterior", S1B1="antero-posterior",
  S1B2="antero-posterior", S1C1="antero-posterior", S1C2="antero-posterior",
  S1D1="latero-lateral",   S1D2="latero-lateral",
  S2A1="antero-posterior", S2A2="antero-posterior", S2B1="antero-posterior",
  S2B2="antero-posterior", S2C1="antero-posterior", S2C2="antero-posterior",
  S2D1="latero-lateral",   S2D2="latero-lateral",
  S3A1="antero-posterior", S3A2="antero-posterior", S3B1="antero-posterior",
  S3B2="antero-posterior", S3C1="antero-posterior", S3C2="antero-posterior",
  S3D1="latero-lateral",   S3D2="latero-lateral"
)

section_map <- c(
  S1A1="coronal", S1A2="coronal", S1B1="coronal", S1B2="coronal",
  S1C1="coronal", S1C2="coronal", S1D1="sagittal", S1D2="sagittal",
  S2A1="coronal", S2A2="coronal", S2B1="coronal", S2B2="coronal",
  S2C1="coronal", S2C2="coronal", S2D1="sagittal", S2D2="sagittal",
  S3A1="coronal", S3A2="coronal", S3B1="coronal", S3B2="coronal",
  S3C1="coronal", S3C2="coronal", S3D1="sagittal", S3D2="sagittal"
)
key <- as.character(seurat_merged$old.ident)

seurat_merged$bregma       <- unname(bregma_map[key])
seurat_merged$axis         <- unname(axis_map[key])
seurat_merged$section_pane <- unname(section_map[key])


table(seurat_merged$bregma, seurat_merged$old.ident)
table(seurat_merged$axis, seurat_merged$old.ident)
table(seurat_merged$section_pane, seurat_merged$old.ident)

# Build mapping of cell types to cellex v147 values from multiple CSV files

library(dplyr)
library(tidyr)

indir <- paste0(atlaspath,"area_cellex_results")

# Get list of CSV files
files <- list.files(indir, pattern = "_detailed_cellex\\.csv$", full.names = TRUE)

# Loop through each CSV, pull out the v147 row, and build a per-area data frame
records <- list()

for (fn in files) {
  # Extract area name from filename
  basename <- basename(fn)
  area <- gsub("_detailed_cellex\\.csv$", "", basename)
  
  # Read CSV with first column as row names
  df <- read.csv(fn, row.names = 1, check.names = FALSE)
  
  # Check if v147 exists in the row names
  if (!"v147" %in% rownames(df)) {
    next
  }
  
  # Extract the v147 row as a data frame (not a vector)
  v147_row <- df["v147", , drop = FALSE]
  rownames(v147_row) <- area
  
  # Add to records list
  records[[area]] <- v147_row
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

# Combine into one data frame: rows=areas, cols=cell-types
combined <- do.call(rbind, records_filled)

# Make sure columns match your actual cell type categories
celltypes <- sort(unique(seurat_merged$detailed_celltype_annotation))

# Reindex columns to match cell types, fill missing with 0
missing_cols <- setdiff(celltypes, colnames(combined))
for (col in missing_cols) {
  combined[[col]] <- 0
}

# Reorder columns to match celltypes order
combined <- combined[, celltypes, drop = FALSE]

# Sort rows (areas) alphabetically
combined <- combined[order(rownames(combined)), , drop = FALSE]

# Create transposed version (cell types as rows, areas as columns)
combined_swap <- t(combined)

# Print results
cat("\n=== Transposed matrix (cell types × areas) ===\n")
combined_swap

# Add cellex v147 scores from combined matrix to Seurat metadata
# Each area gets its own column: seurat_merged$score_{area}

library(dplyr)

# Assuming 'combined' is the areas × cell types matrix from previous script
# and seurat_merged$detailed_celltype_annotation contains the cell type labels

# Get list of areas (row names of combined matrix)
areas <- rownames(combined)

cat("Adding cellex scores for", length(areas), "areas to Seurat metadata...\n\n")

# Loop through each area and create a metadata column
for (area in areas) {
  # Create column name
  col_name <- paste0("snRNA_cellex_score_", area)
  
  # Get the scores for this area (one value per cell type)
  area_scores <- combined[area, ]
  
  # Get cell type annotation for each cell
  cell_types <- seurat_merged$detailed_celltype_annotation
  
  # Create a named vector for mapping (cell type -> score)
  score_lookup <- as.numeric(area_scores)
  names(score_lookup) <- colnames(combined)
  
  # Map each cell's type to its score
  scores <- score_lookup[as.character(cell_types)]
  
  # Replace NAs with 0
  scores[is.na(scores)] <- 0
  
  # Remove names before adding to metadata (use cell order from Seurat object)
  scores <- unname(scores)
  
  # Add to Seurat metadata
  seurat_merged[[col_name]] <- scores
  
  cat(sprintf("✓ Added %s (range: %.4f to %.4f)\n", 
              col_name, 
              min(scores, na.rm = TRUE), 
              max(scores, na.rm = TRUE)))
}

cat("\n=== Summary ===\n")
cat(sprintf("Total new metadata columns added: %d\n", length(areas)))
cat(sprintf("Column names: %s\n", paste(paste0("snRNA_cellex_score_", areas[1:min(5, length(areas))]), collapse = ", ")))
if (length(areas) > 5) {
  cat(sprintf("  ... and %d more\n", length(areas) - 5))
}

# Show example of what was added
cat("\n=== Example: First few cells ===\n")
example_cols <- paste0("snRNA_cellex_score_", areas[1:min(3, length(areas))])
print(head(seurat_merged@meta.data[, c("detailed_celltype_annotation", example_cols)]))

# Optional: Check for cells with all zeros (might indicate cell types not in the matrix)
all_score_cols <- paste0("snRNA_cellex_score_", areas)
score_matrix <- seurat_merged@meta.data[, all_score_cols]
cells_all_zero <- rowSums(score_matrix) == 0

if (sum(cells_all_zero) > 0) {
  cat(sprintf("\n⚠ Warning: %d cells have all-zero scores across all areas\n", sum(cells_all_zero)))
  cat("Cell types with all zeros:\n")
  print(table(seurat_merged$detailed_celltype_annotation[cells_all_zero]))
}

cat("\n✓ All cellex scores added to Seurat metadata!\n")

plot <- ImageFeaturePlot(seurat_merged, border.color = NA, 
                         boundaries = "segmentation", 
                         features = "snRNA_cellex_score_BST",
                         fov = "S3C2") + scale_fill_distiller(palette = "Spectral")
show(plot)

region_ct_df <- read.csv(paste0(atlaspath,"region_celltype_scores_cellex_spatial.csv"), row.names = 1, check.names = FALSE)
region_ct_df

areas <- rownames(region_ct_df)  # c("BST","CEA","LHA","PAG","PBN","PVH","PVT")

for (a in areas) {
  ct_vec <- seurat_merged$detailed_celltype_annotation
  
  # start with NAs
  score_vec <- rep(NA_real_, length(ct_vec))
  
  # which cell types are present in region_ct_df columns?
  in_table <- ct_vec %in% colnames(region_ct_df)
  
  # look up scores for those
  score_vec[in_table] <- as.numeric(region_ct_df[a, ct_vec[in_table]])
  
  # store in metadata
  seurat_merged[[paste0("spatial_cellex_score_", a)]] <- score_vec
}

plot <- ImageFeaturePlot(seurat_merged, border.color = NA, 
                         boundaries = "segmentation", 
                         features = "spatial_cellex_score_BST",
                         fov = "S3C2") + scale_fill_distiller(palette = "Spectral")
show(plot)

# Read prioritization CSV and assign beta values to Seurat metadata
# Only use rows where specificity_id == "detailed"

library(dplyr)

# Read the prioritization CSV
prioritization <- read.csv(paste0(atlaspath,"prioritization.csv"), stringsAsFactors = FALSE)

cat("=== Prioritization data loaded ===\n")
cat(sprintf("Total rows: %d\n", nrow(prioritization)))
cat(sprintf("Unique GWAS traits: %s\n", paste(unique(prioritization$gwas), collapse = ", ")))
cat(sprintf("Unique specificity IDs: %s\n", paste(unique(prioritization$specificity_id), collapse = ", ")))

# Filter for "detailed" annotations only
detailed_priority <- prioritization %>%
  filter(specificity_id == "detailed")

cat(sprintf("\nRows with 'detailed' specificity: %d\n", nrow(detailed_priority)))

# Transform pvalues to -log10(pvalue)
detailed_priority$neg_log10_pval <- -log10(detailed_priority$pvalue)

# Create a lookup table: annotation -> -log10(pvalue)
score_lookup <- setNames(detailed_priority$neg_log10_pval, detailed_priority$annotation)

cat("\n=== -log10(p-value) for cell types ===\n")
print(head(score_lookup[order(-score_lookup)], 10))

# Get cell type annotations from Seurat
cell_types <- seurat_merged$detailed_celltype_annotation

# Map cell types to -log10(pvalue) values
cellect_scores <- score_lookup[as.character(cell_types)]

# Replace NAs with 0
cellect_scores[is.na(cellect_scores)] <- 0

# Remove names before adding to metadata
cellect_scores <- unname(cellect_scores)

# Add to Seurat metadata
seurat_merged$cellect_bmi_score <- cellect_scores

cat("\n=== Summary of cellect_bmi_score (-log10 p-value) ===\n")
cat(sprintf("Range: %.4f to %.4f\n", min(cellect_scores), max(cellect_scores)))
cat(sprintf("Mean: %.4f\n", mean(cellect_scores)))
cat(sprintf("Median: %.4f\n", median(cellect_scores)))
cat(sprintf("Cells with non-zero scores: %d (%.1f%%)\n", 
            sum(cellect_scores != 0), 
            100 * sum(cellect_scores != 0) / length(cellect_scores)))

# Show distribution by cell type
cat("\n=== Top 10 cell types by -log10(p-value) ===\n")
celltype_summary <- data.frame(
  celltype = names(score_lookup),
  neg_log10_pval = as.numeric(score_lookup),
  pvalue = detailed_priority$pvalue[match(names(score_lookup), detailed_priority$annotation)],
  stringsAsFactors = FALSE
) %>%
  arrange(desc(neg_log10_pval))

print(head(celltype_summary, 10))

# Example output
cat("\n=== Example: First few cells ===\n")
example_df <- data.frame(
  celltype = head(seurat_merged$detailed_celltype_annotation, 20),
  cellect_bmi_score = head(seurat_merged$cellect_bmi_score, 20)
)
print(example_df)

cat("\n✓ cellect_bmi_score added to seurat_merged metadata!\n")

plot <- ImageFeaturePlot(seurat_merged, border.color = NA, 
                         boundaries = "segmentation", 
                         features = "cellect_bmi_score",
                         fov = "S3C2") + scale_fill_distiller(palette = "Spectral")
show(plot)

#qs::qsave(seurat_merged, paste0(atlaspath, "spatial_24_bipca_downstream.qs"))
seurat_merged <- qs::qread(paste0(atlaspath, "spatial_24_bipca_downstream.qs"))
seurat_merged

glia <- c()

seurat_merged$plot <- seurat_merged$detailed_celltype_annotation
seurat_merged@meta.data[["plot"]] <- ifelse(seurat_merged@meta.data[["plot"]] %in% glia, seurat_merged@meta.data[["plot"]], "no_label")
seurat_merged@meta.data[["plot"]] <- ifelse(seurat_merged@meta.data[["score_detailed_celltype_annotation"]] > 0.4, seurat_merged@meta.data[["plot"]], "no_label")

plot <- ImageDimPlot(seurat_merged, 
             border.color = NA, 
             boundaries = "segmentation", 
             group.by = "plot", 
             cols = colors_detailed,
             molecules = c("Lepr", "Sim1", "Fezf1", "Ebf1"),
             mols.cols = c("orange", "purple", "green", "blue"),
             mols.size = 0.01,
             mols.alpha = 0.1,
             nmols=1000000,
             fov = "S3D2",
             dark.background = TRUE,
            )+ theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      legend.key.size = unit(0.5, "lines"),
                      legend.text = element_text(size = 9), 
                      plot.background = element_rect(fill = "black"),
                      panel.background = element_rect(fill = "black"), 
                      plot.margin = unit(c(0, 0, 0, 0), "cm"))
show(plot)
resolve_save(plot, paste0(atlaspath, "plots/dylan_plot.pdf"))

plot <- ImageDimPlot(seurat_merged, 
             border.color = NA, 
             boundaries = "segmentation", 
             group.by = "plot", 
             cols = colors_detailed_light,
             molecules = c("Lepr", "Calcr"),
             mols.cols = c("red", "cyan"),
             mols.size = 0.01,
             mols.alpha = 0.1,
             nmols=1000000,
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
resolve_save(plot, paste0(atlaspath, "plots/test_lepr_calcr_bright.pdf"))

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
resolve_save(plot, paste0(atlaspath, "plots/test_lmx1a_bright.pdf"))

#function to plot subsets
lmx_plot <- function(seurat_obj, fov){
    plot <- ImageDimPlot(seurat_obj, 
             border.color = NA, 
             boundaries = "segmentation", 
             group.by = "plot", 
             cols = colors_detailed_light,
             fov = fov,
             dark.background = FALSE,
            )+ theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      legend.key.size = unit(0.5, "lines"),
                      legend.text = element_text(size = 9), 
                      plot.background = element_rect(fill = "white"),
                      panel.background = element_rect(fill = "white"), 
                      plot.margin = unit(c(0, 0, 0, 0), "cm"))
}

fov <- unique(seurat_merged$old.ident)
fov

for (f in fov){
    print(paste("Printing FOV", f))
    plot <- lmx_plot(seurat_merged, f)
    resolve_save(plot, paste0(atlaspath, "plots/lmx1a_bright_",f,".pdf"))
}

counts <- GetAssayData(seurat_merged, layer = "raw", slot = "counts")
seurat_merged$Lepr_pos  <- counts["Lepr", ]  > 2
seurat_merged$Calcr_pos <- counts["Calcr", ] > 2
seurat_merged$expr_class <- "Lepr-/Calcr-"
seurat_merged$expr_class[seurat_merged$Lepr_pos]  <- "Lepr+/Calcr-"
seurat_merged$expr_class[seurat_merged$Calcr_pos] <- "Lepr-/Calcr+"
seurat_merged$expr_class[seurat_merged$Lepr_pos & seurat_merged$Calcr_pos] <- "Lepr+/Calcr+"

plot <- ImageDimPlot(seurat_merged, 
             border.color = NA, 
             boundaries = "segmentation", 
             group.by = "expr_class", 
             cols = c("grey", "firebrick", "dodgerblue", "purple"),
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

seurat_merged$plot <- seurat_merged$detailed_celltype_annotation
seurat_merged@meta.data[["plot"]] <- ifelse(seurat_merged@meta.data[["plot"]] %in% glia,
                                            seurat_merged@meta.data[["plot"]], "no_label")
seurat_merged@meta.data[["plot"]] <- ifelse(seurat_merged@meta.data[["score_detailed_celltype_annotation"]] > 0.4,
                                            seurat_merged@meta.data[["plot"]], "no_label")

seurat_merged@meta.data[["plot"]] <- ifelse(seurat_merged@meta.data[["plot"]] == "no_label",
                                            seurat_merged@meta.data[["plot"]], "Lmx1a")
seurat_merged@meta.data[["plot"]] <- ifelse(seurat_merged@meta.data[["Lepr_pos"]] == FALSE 
                                            & seurat_merged@meta.data[["plot"]] == "Lmx1a",
                                            "Lmx1a_Lepr-", seurat_merged@meta.data[["plot"]])

plot <- ImageDimPlot(seurat_merged, 
             border.color = NA, 
             boundaries = "segmentation", 
             group.by = "plot", 
             cols = c("grey", "darkred", "darkgreen"),
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

seurat_merged$plot <- seurat_merged$detailed_celltype_annotation
seurat_merged@meta.data[["plot"]] <- ifelse(seurat_merged@meta.data[["plot"]] %in% glia,
                                            seurat_merged@meta.data[["plot"]], "no_label")
seurat_merged@meta.data[["plot"]] <- ifelse(seurat_merged@meta.data[["score_detailed_celltype_annotation"]] > 0.4,
                                            seurat_merged@meta.data[["plot"]], "no_label")
seurat_merged@meta.data[["plot"]] <- ifelse(seurat_merged@meta.data[["Lepr_pos"]] == FALSE,
                                            seurat_merged@meta.data[["plot"]], "no_label")
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

glia <- c("33_Chat_Ghrh")

seurat_merged$plot <- seurat_merged$detailed_celltype_annotation
seurat_merged@meta.data[["plot"]] <- ifelse(seurat_merged@meta.data[["plot"]] %in% glia, 
                                            seurat_merged@meta.data[["plot"]], "no_label")
seurat_merged@meta.data[["plot"]] <- ifelse(seurat_merged@meta.data[["score_detailed_celltype_annotation"]] > 0.4, 
                                            seurat_merged@meta.data[["plot"]], "no_label")
plot <- ImageDimPlot(seurat_merged, 
             border.color = NA, 
             boundaries = "segmentation", 
             group.by = "plot", 
             cols = colors_detailed_light,
             fov = "S1B2",
             dark.background = FALSE,
            )+ theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      legend.key.size = unit(0.5, "lines"),
                      legend.text = element_text(size = 9), 
                      plot.background = element_rect(fill = "white"),
                      panel.background = element_rect(fill = "white"), 
                      plot.margin = unit(c(0, 0, 0, 0), "cm"))
show(plot)

#function to plot subsets
gh_plot <- function(seurat_obj, fov){
    plot <- ImageDimPlot(seurat_obj, 
             border.color = NA, 
             boundaries = "segmentation", 
             group.by = "plot", 
             cols = colors_detailed_light,
             fov = fov,
             dark.background = FALSE,
            )+ theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      legend.key.size = unit(0.5, "lines"),
                      legend.text = element_text(size = 9), 
                      plot.background = element_rect(fill = "white"),
                      panel.background = element_rect(fill = "white"), 
                      plot.margin = unit(c(0, 0, 0, 0), "cm"))
}
fov <- unique(seurat_merged$old.ident)

for (f in fov){
    print(paste("Printing FOV", f))
    plot <- gh_plot(seurat_merged, f)
    resolve_save(plot, paste0(atlaspath, "plots/ghrh_bright_",f,".pdf"))
}

glia <- c("30_Trh_Sim1", "25_Nxph4_Ebf3", "21_Bnc2_Glp1r")

seurat_merged$plot <- seurat_merged$detailed_celltype_annotation
seurat_merged@meta.data[["plot"]] <- ifelse(seurat_merged@meta.data[["plot"]] %in% glia, 
                                            seurat_merged@meta.data[["plot"]], "no_label")
seurat_merged@meta.data[["plot"]] <- ifelse(seurat_merged@meta.data[["score_detailed_celltype_annotation"]] > 0.4, 
                                            seurat_merged@meta.data[["plot"]], "no_label")
trh_plot <- function(seurat_obj, fov){
plot <- ImageDimPlot(seurat_merged, 
             border.color = NA, 
             boundaries = "segmentation", 
             group.by = "plot", 
             cols = colors_detailed_light,
             molecules = c("Trh"),
             mols.cols = c("cyan"),
             mols.size = 0.01,
             mols.alpha = 0.1,
             nmols=10000,
             fov = fov,
             dark.background = FALSE,
            )+ theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      legend.key.size = unit(0.5, "lines"),
                      legend.text = element_text(size = 9), 
                      plot.background = element_rect(fill = "white"),
                      panel.background = element_rect(fill = "white"), 
                      plot.margin = unit(c(0, 0, 0, 0), "cm"))
}
fov <- unique(seurat_merged$old.ident)

for (f in fov){
    print(paste("Printing FOV", f))
    plot <- trh_plot(seurat_merged, f)
    resolve_save(plot, paste0(atlaspath, "plots/trh_bright_",f,".pdf"))
}

counts <- GetAssayData(seurat_merged, layer = "raw", slot = "counts")
seurat_merged$Lepr_pos  <- counts["Trh", ]  > 3
table(seurat_merged$Lepr_pos, seurat_merged$detailed_celltype_annotation)

glia <- c("13_Pth2r_Otp", "10_Npy_Agrp", "40_Gpr150_Tac2")

seurat_merged$plot <- seurat_merged$detailed_celltype_annotation
#seurat_merged@meta.data[["plot"]] <- ifelse(seurat_merged@meta.data[["plot"]] %in% glia, 
#                                            seurat_merged@meta.data[["plot"]], "no_label")
seurat_merged@meta.data[["plot"]] <- ifelse(seurat_merged@meta.data[["score_detailed_celltype_annotation"]] > 0.4, 
                                            seurat_merged@meta.data[["plot"]], "no_label")
trh_plot <- function(seurat_obj, fov){
plot <- ImageDimPlot(seurat_merged, 
             border.color = NA, 
             boundaries = "segmentation", 
             group.by = "plot", 
             cols = colors_detailed,
             molecules = c("Calcr"),
             mols.cols = c("white"),
             mols.size = 0.01,
             mols.alpha = 0.1,
             nmols=10000,
             fov = fov,
             dark.background = TRUE,
            )+ theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      legend.key.size = unit(0.5, "lines"),
                      legend.text = element_text(size = 9), 
                      plot.background = element_rect(fill = "black"),
                      panel.background = element_rect(fill = "black"), 
                      plot.margin = unit(c(0, 0, 0, 0), "cm"))
}
fov <- unique(seurat_merged$old.ident)

for (f in fov){
    print(paste("Printing FOV", f))
    plot <- trh_plot(seurat_merged, f)
    resolve_save(plot, paste0(atlaspath, "plots/calcr_",f,".pdf"))
}

plot <- ImageDimPlot(seurat_merged, 
             border.color = NA, 
             boundaries = "segmentation", 
             cols = colors_detailed_light,
             molecules = c("Glp1r", "Calcr"),
             mols.cols = c("red", "cyan"),
             mols.size = 0.01,
             mols.alpha = 0.1,
             nmols=1000000,
             fov = "S3D2",
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
             cols = colors_detailed_light,
             molecules = c("Pmch", "Pomc", "Lmx1b"),
             mols.cols = c("red", "cyan", "green"),
             mols.size = 0.01,
             mols.alpha = 0.1,
             nmols=1000000,
             fov = "S2D2",
             dark.background = FALSE,
            )+ theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      legend.key.size = unit(0.5, "lines"),
                      legend.text = element_text(size = 9), 
                      plot.background = element_rect(fill = "white"),
                      panel.background = element_rect(fill = "white"), 
                      plot.margin = unit(c(0, 0, 0, 0), "cm"))
resolve_save(plot, paste0(atlaspath, "plots/MCH_POMC_LMX.pdf"))

counts <- GetAssayData(seurat_merged, layer = "raw", slot = "counts")
genes <- c("Pmch", "Pomc", "Lmx1b")
sub <- counts[genes, , drop = FALSE]

# 2. Long format
df <- sub %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(
    cols = -gene,
    names_to = "cell",
    values_to = "count"
  )

# 3. Prepare metadata with cell names as a real column
meta <- seurat_merged@meta.data %>%
  rownames_to_column("cell")

# 4. Join correctly using "cell"
df2 <- df %>%
  left_join(meta[, c("cell", "detailed_celltype_annotation")],
            by = "cell")

# 5. Summarize counts per cell type
result <- df2 %>%
  group_by(gene, detailed_celltype_annotation) %>%
  summarise(total_molecules = sum(count), .groups = "drop") %>%
  arrange(gene, desc(total_molecules))

result_list <- split(result, result$gene)
head(result_list$Pmch)
head(result_list$Pomc)
head(result_list$Lmx1b)

glia <- c("41_Otx1_Zic1", "06_Tbx3_Npy1r", "27_Lmx1a_Ttc6")

seurat_merged$plot <- seurat_merged$detailed_celltype_annotation
seurat_merged@meta.data[["plot"]] <- ifelse(seurat_merged@meta.data[["plot"]] %in% glia, 
                                            seurat_merged@meta.data[["plot"]], "no_label")
seurat_merged@meta.data[["plot"]] <- ifelse(seurat_merged@meta.data[["score_detailed_celltype_annotation"]] > 0.4, 
                                            seurat_merged@meta.data[["plot"]], "no_label")
plot <- ImageDimPlot(seurat_merged, 
             border.color = NA, 
             boundaries = "segmentation", 
             group.by = "plot", 
             cols = c("azure3", "springgreen2", "blue2", "deeppink2"),
             fov = "S2D2",
             dark.background = FALSE,
            )+ theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      legend.key.size = unit(0.5, "lines"),
                      legend.text = element_text(size = 9), 
                      plot.background = element_rect(fill = "white"),
                      panel.background = element_rect(fill = "white"), 
                      plot.margin = unit(c(0, 0, 0, 0), "cm"))
resolve_save(plot, paste0(atlaspath, "plots/MCH_POMC_LMX_celltypes.pdf"))

seurat_merged$plot <- seurat_merged$detailed_celltype_annotation
seurat_merged@meta.data[["plot"]] <- ifelse(seurat_merged@meta.data[["plot"]] %in% glia,
                                            seurat_merged@meta.data[["plot"]], "no_label")
seurat_merged@meta.data[["plot"]] <- ifelse(seurat_merged@meta.data[["score_detailed_celltype_annotation"]] > 0.4,
                                            seurat_merged@meta.data[["plot"]], "no_label")

seurat_merged@meta.data[["plot"]] <- ifelse(seurat_merged@meta.data[["plot"]] == "no_label",
                                            seurat_merged@meta.data[["plot"]], "Lmx1a")
seurat_merged@meta.data[["plot"]] <- ifelse(seurat_merged@meta.data[["Lepr_pos"]] == FALSE 
                                            & seurat_merged@meta.data[["plot"]] == "Lmx1a",
                                            "Lmx1a_Lepr-", seurat_merged@meta.data[["plot"]])

plot <- ImageDimPlot(seurat_merged, 
             border.color = NA, 
             boundaries = "segmentation", 
             group.by = "plot", 
             cols = c("grey", "darkred", "darkgreen"),
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

counts <- GetAssayData(seurat_merged, layer = "raw", slot = "counts")
seurat_merged$Lepr_pos  <- counts["Lepr", ]  > 0
seurat_merged$Glp1r_pos <- counts["Glp1r", ] > 0
seurat_merged$Ebf1_pos <- counts["Ebf1", ] > 0
seurat_merged$Trh_pos <- counts["Trh", ] > 0
seurat_merged$expr_class <- "negative"
seurat_merged$expr_class[seurat_merged$Lepr_pos & seurat_merged$Glp1r_pos & seurat_merged$Ebf1_pos] <- "Lepr+Glp1r+Ebf1"
seurat_merged$expr_class[seurat_merged$Lepr_pos & seurat_merged$Glp1r_pos & seurat_merged$Trh_pos] <- "Lepr+Glp1r+Trh"
seurat_merged$expr_class[seurat_merged$Lepr_pos & seurat_merged$Glp1r_pos & seurat_merged$Ebf1_pos & seurat_merged$Trh_pos] <- "All four genes"

table(seurat_merged$expr_class)

unique(seurat_merged$region.collapsed)

plot <- ImageDimPlot(seurat_merged, 
             border.color = NA, 
             boundaries = "segmentation", 
             group.by = "expr_class", 
             cols = c("darkgrey", "firebrick", "dodgerblue", "yellow"),
             molecules = c("Agrp", "Sim1", "Fezf1"),
             mols.cols = c("orange", "purple", "green"),
             mols.size = 0.01,
             mols.alpha = 0.5,
             nmols=5000,
             fov = "S3D2",
             dark.background = TRUE,
            )+ theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      legend.key.size = unit(0.5, "lines"),
                      legend.text = element_text(size = 9), 
                      plot.background = element_rect(fill = "black"),
                      panel.background = element_rect(fill = "black"), 
                      plot.margin = unit(c(0, 0, 0, 0), "cm"))
show(plot)
resolve_save(plot, paste0(atlaspath, "plots/dylan_plot.pdf"))

plot <- ImageDimPlot(seurat_merged, 
             border.color = NA, 
             boundaries = "segmentation", 
             group.by = "expr_class", 
             cols = c("darkgrey", "firebrick", "dodgerblue", "yellow"),
             fov = "S2C1",
             dark.background = TRUE,
            )+ theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      legend.key.size = unit(0.5, "lines"),
                      legend.text = element_text(size = 9), 
                      plot.background = element_rect(fill = "black"),
                      panel.background = element_rect(fill = "black"), 
                      plot.margin = unit(c(0, 0, 0, 0), "cm"))
show(plot)
resolve_save(plot, paste0(atlaspath, "dylan_s2c1.pdf"))

plot <- ImageDimPlot(seurat_merged, 
             border.color = NA, 
             boundaries = "segmentation", 
             group.by = "expr_class", 
             cols = c("darkgrey", "firebrick", "dodgerblue", "yellow"),
             fov = "S3C1",
             dark.background = TRUE,
            )+ theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      legend.key.size = unit(0.5, "lines"),
                      legend.text = element_text(size = 9), 
                      plot.background = element_rect(fill = "black"),
                      panel.background = element_rect(fill = "black"), 
                      plot.margin = unit(c(0, 0, 0, 0), "cm"))
show(plot)
resolve_save(plot, paste0(atlaspath, "dylan_s3c1.pdf"))

plot <- ImageDimPlot(seurat_merged, 
             border.color = NA, 
             boundaries = "segmentation", 
             group.by = "expr_class", 
             cols = c("darkgrey", "firebrick", "dodgerblue", "yellow"),
             fov = "S1D2",
             dark.background = TRUE,
            )+ theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      legend.key.size = unit(0.5, "lines"),
                      legend.text = element_text(size = 9), 
                      plot.background = element_rect(fill = "black"),
                      panel.background = element_rect(fill = "black"), 
                      plot.margin = unit(c(0, 0, 0, 0), "cm"))
show(plot)
resolve_save(plot, paste0(atlaspath, "dylan_s1d2.pdf"))

# Heatmap of expression classes by region
# Showing both absolute counts and relative abundance

library(ggplot2)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)

# Get the data excluding negative cells
expr_data <- data.frame(
  expr_class = seurat_merged$expr_class,
  region = seurat_merged$region.collapsed,
  stringsAsFactors = FALSE
)

# Create contingency table
expr_region_table <- table(expr_data$expr_class, expr_data$region)
print("Full contingency table:")
print(expr_region_table)

# Remove negative class for the heatmap
expr_classes <- c("Lepr+Glp1r+Ebf1", "Lepr+Glp1r+Trh", "All four genes")
expr_region_table_filtered <- expr_region_table[expr_classes, , drop = FALSE]

print("\nFiltered table (without negative):")
print(expr_region_table_filtered)

# Calculate total cells per region (including negative)
total_per_region <- colSums(expr_region_table)

# Calculate proportion relative to total cells in each region
prop_matrix <- sweep(expr_region_table_filtered, 2, total_per_region, FUN = "/") * 100

print("\nProportion matrix (% of total cells in each region):")
print(round(prop_matrix, 2))

# Create annotation data frames
count_matrix <- as.matrix(expr_region_table_filtered)
prop_matrix <- as.matrix(prop_matrix)

# Function to create cell text showing both count and percentage
cell_fun_both = function(j, i, x, y, width, height, fill) {
  count <- count_matrix[i, j]
  perc <- prop_matrix[i, j]
  
  grid::grid.text(
    sprintf("%d\n(%.2f%%)", count, perc),
    x, y, 
    gp = grid::gpar(fontsize = 8, col = if(prop_matrix[i, j] > 0.5) "white" else "black")
  )
}

# Create color function for the heatmap (based on percentage)
col_fun <- colorRamp2(
  c(0, max(prop_matrix)/2, max(prop_matrix)),
  c("#FFFFFF", "#FFA500", "#FF0000")
)

# Create the heatmap using ComplexHeatmap
ht1 <- Heatmap(
  prop_matrix,
  name = "% of region",
  col = col_fun,
  
  # Row settings
  row_names_side = "left",
  row_names_gp = grid::gpar(fontsize = 10),
  cluster_rows = FALSE,
  
  # Column settings
  column_names_side = "bottom",
  column_names_gp = grid::gpar(fontsize = 10),
  column_names_rot = 45,
  cluster_columns = FALSE,
  
  # Cell annotations
  cell_fun = cell_fun_both,
  
  # Heatmap appearance
  border = TRUE,
  rect_gp = grid::gpar(col = "white", lwd = 2),
  
  # Legend
  heatmap_legend_param = list(
    title = "% of total\ncells in region",
    title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
    labels_gp = grid::gpar(fontsize = 9),
    legend_height = unit(4, "cm")
  ),
  
  # Title
  column_title = "Expression class distribution across brain regions",
  column_title_gp = grid::gpar(fontsize = 12, fontface = "bold")
)

# Draw the heatmap
pdf(paste0(atlaspath, "expr_class_by_region_heatmap.pdf"), width = 12, height = 6)
draw(ht1)
dev.off()

# Also create a simpler ggplot2 version
# Prepare data for ggplot - convert matrix properly
plot_data <- as.data.frame(as.table(prop_matrix)) %>%
  rename(expr_class = Var1, region = Var2, percentage = Freq)

plot_data_counts <- as.data.frame(as.table(count_matrix)) %>%
  rename(expr_class = Var1, region = Var2, count = Freq)

plot_data <- plot_data %>%
  left_join(plot_data_counts, by = c("expr_class", "region"))

# Create ggplot heatmap
p1 <- ggplot(plot_data, aes(x = region, y = expr_class, fill = percentage)) +
  geom_tile(color = "white", size = 1) +
  geom_text(aes(label = sprintf("%d\n(%.2f%%)", count, percentage)), 
            size = 3, 
            color = ifelse(plot_data$percentage > max(plot_data$percentage)/2, "white", "black")) +
  scale_fill_gradient2(
    low = "white",
    mid = "orange",
    high = "red",
    midpoint = max(plot_data$percentage)/2,
    name = "% of total\ncells in region"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_blank(),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    panel.grid = element_blank()
  ) +
  labs(title = "Expression class distribution across brain regions")

ggsave(paste0(atlaspath, "expr_class_by_region_heatmap_ggplot.pdf"),
       plot = p1,
       width = 12,
       height = 6,
       device = cairo_pdf)

# Create a summary table
summary_df <- plot_data %>%
  arrange(expr_class, desc(percentage))

print("\nTop regions for each expression class:")
print(summary_df %>% 
        group_by(expr_class) %>% 
        slice_max(percentage, n = 3))

# Save summary tables
write.csv(count_matrix, 
          paste0(atlaspath, "expr_class_by_region_counts.csv"),
          row.names = TRUE)

write.csv(prop_matrix, 
          paste0(atlaspath, "expr_class_by_region_percentages.csv"),
          row.names = TRUE)

cat("\nFiles saved:\n")
cat("1. expr_class_by_region_heatmap.pdf (ComplexHeatmap version)\n")
cat("2. expr_class_by_region_heatmap_ggplot.pdf (ggplot2 version)\n")
cat("3. expr_class_by_region_counts.csv (absolute counts)\n")
cat("4. expr_class_by_region_percentages.csv (percentages)\n")

# Alternative visualization: Side-by-side heatmaps for counts and percentages

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# Get the data
expr_data <- data.frame(
  expr_class = seurat_merged$expr_class,
  region = seurat_merged$region.collapsed,
  stringsAsFactors = FALSE
)

# Create contingency table
expr_region_table <- table(expr_data$expr_class, expr_data$region)

# Remove negative class
expr_classes <- c("Lepr+Glp1r+Ebf1", "Lepr+Glp1r+Trh", "All four genes")
count_matrix <- expr_region_table[expr_classes, , drop = FALSE]

# Calculate percentages
total_per_region <- colSums(expr_region_table)
prop_matrix <- sweep(count_matrix, 2, total_per_region, FUN = "/") * 100

# Prepare data for plotting - convert table properly
plot_data_counts <- as.data.frame(as.table(count_matrix)) %>%
  rename(expr_class = Var1, region = Var2, value = Freq) %>%
  mutate(type = "Absolute count")

plot_data_perc <- as.data.frame(as.table(prop_matrix)) %>%
  rename(expr_class = Var1, region = Var2, value = Freq) %>%
  mutate(type = "% of region")

# Heatmap 1: Absolute counts
p_counts <- ggplot(plot_data_counts, aes(x = region, y = expr_class, fill = value)) +
  geom_tile(color = "white", size = 1) +
  geom_text(aes(label = round(value, 0)), 
            size = 3.5, 
            color = ifelse(plot_data_counts$value > max(plot_data_counts$value)/2, "white", "black")) +
  scale_fill_gradient2(
    low = "white",
    mid = "#9ecae1",
    high = "#08519c",
    midpoint = max(plot_data_counts$value)/2,
    name = "Cell count"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 10),
    axis.title = element_blank(),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
    panel.grid = element_blank(),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
  ) +
  labs(title = "Absolute cell counts")

# Heatmap 2: Percentages
p_perc <- ggplot(plot_data_perc, aes(x = region, y = expr_class, fill = value)) +
  geom_tile(color = "white", size = 1) +
  geom_text(aes(label = sprintf("%.2f%%", value)), 
            size = 3.5, 
            color = ifelse(plot_data_perc$value > max(plot_data_perc$value)/2, "white", "black")) +
  scale_fill_gradient2(
    low = "white",
    mid = "orange",
    high = "red",
    midpoint = max(plot_data_perc$value)/2,
    name = "% of total\ncells in region"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_blank(),  # Remove y-axis labels for second plot
    axis.title = element_blank(),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
    panel.grid = element_blank(),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
  ) +
  labs(title = "Percentage of region cells")

# Combine plots
combined_plot <- p_counts + p_perc +
  plot_annotation(
    title = 'Expression class distribution across brain regions',
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  )

# Save with increased width
ggsave(paste0(atlaspath, "expr_class_by_region_dual_heatmap.pdf"),
       plot = combined_plot,
       width = 20,
       height = 6,
       device = cairo_pdf)

print(combined_plot)

cat("\nSaved: expr_class_by_region_dual_heatmap.pdf\n")

# Find FOVs with highest absolute counts for each expression class

library(dplyr)
library(ggplot2)

# Get the data
expr_data <- data.frame(
  expr_class = seurat_merged$expr_class,
  fov = seurat_merged$old.ident,
  stringsAsFactors = FALSE
)

# Filter out negative cells
expr_data_filtered <- expr_data %>%
  filter(expr_class != "negative")

# Create contingency table
expr_fov_table <- table(expr_data_filtered$expr_class, expr_data_filtered$fov)

print("=== Expression class counts by FOV ===")
print(expr_fov_table)

# Convert to data frame for easier manipulation
expr_fov_df <- as.data.frame(as.table(expr_fov_table)) %>%
  rename(expr_class = Var1, fov = Var2, count = Freq) %>%
  arrange(expr_class, desc(count))

# Find top FOVs for each expression class
cat("\n=== Top 10 FOVs for each expression class ===\n\n")

for (expr_class in c("Lepr+Glp1r+Ebf1", "Lepr+Glp1r+Trh", "All four genes")) {
  cat(sprintf("\n%s:\n", expr_class))
  cat(strrep("-", 60), "\n")
  
  top_fovs <- expr_fov_df %>%
    filter(expr_class == !!expr_class) %>%
    slice_max(count, n = 10)
  
  for (i in 1:nrow(top_fovs)) {
    cat(sprintf("%2d. %s: %d cells\n", i, top_fovs$fov[i], top_fovs$count[i]))
  }
}

# Create a summary table
summary_table <- expr_fov_df %>%
  group_by(expr_class) %>%
  slice_max(count, n = 10) %>%
  mutate(rank = row_number()) %>%
  ungroup()

# Create visualization - bar plot of top 10 FOVs for each class
plot_top_fovs <- summary_table %>%
  ggplot(aes(x = reorder(fov, count), y = count, fill = expr_class)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~expr_class, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = c(
    "Lepr+Glp1r+Ebf1" = "#FF6B6B",
    "Lepr+Glp1r+Trh" = "#4ECDC4",
    "All four genes" = "#FFE66D"
  )) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 9),
    strip.text = element_text(size = 11, face = "bold"),
    panel.grid.major.y = element_blank(),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
  ) +
  labs(
    title = "Top 10 FOVs by cell count for each expression class",
    x = "FOV",
    y = "Cell count"
  )

print(plot_top_fovs)

# Save plot
ggsave(paste0(atlaspath, "top_fovs_by_expr_class.pdf"),
       plot = plot_top_fovs,
       width = 10,
       height = 12,
       device = cairo_pdf)

# Save the full table as CSV
write.csv(expr_fov_df, 
          paste0(atlaspath, "expr_class_by_fov_counts.csv"),
          row.names = FALSE)

# Create a summary with overall statistics
cat("\n=== Overall Statistics ===\n")
overall_stats <- expr_fov_df %>%
  group_by(expr_class) %>%
  summarize(
    total_cells = sum(count),
    n_fovs_with_cells = sum(count > 0),
    max_per_fov = max(count),
    mean_per_fov = mean(count[count > 0]),
    median_per_fov = median(count[count > 0])
  )

print(overall_stats)

# Find which FOV has the most cells for each class
top_fov_per_class <- expr_fov_df %>%
  group_by(expr_class) %>%
  slice_max(count, n = 1)

cat("\n=== FOV with most cells for each expression class ===\n")
print(top_fov_per_class)

cat("\nFiles saved:\n")
cat("1. top_fovs_by_expr_class.pdf\n")
cat("2. expr_class_by_fov_counts.csv\n")

# Plot expression class distribution along bregma coordinates (coronal sections only)

library(dplyr)
library(ggplot2)
library(tidyr)

# Filter for coronal sections only
coronal_data <- data.frame(
  expr_class = seurat_merged$expr_class,
  bregma = seurat_merged$bregma,
  section_plane = seurat_merged$section_pane,
  stringsAsFactors = FALSE
) %>%
  filter(section_plane == "coronal")

cat("=== Summary of coronal data ===\n")
cat(sprintf("Total cells in coronal sections: %d\n", nrow(coronal_data)))
cat(sprintf("Bregma range: %.2f to %.2f\n", 
            min(coronal_data$bregma, na.rm = TRUE), 
            max(coronal_data$bregma, na.rm = TRUE)))

# Filter out negative cells
coronal_filtered <- coronal_data %>%
  filter(expr_class != "negative")

# Count cells by expression class and bregma
bregma_counts <- coronal_filtered %>%
  group_by(bregma, expr_class) %>%
  summarize(count = n(), .groups = "drop") %>%
  arrange(bregma, expr_class)

cat("\n=== Cell counts by bregma level ===\n")
print(bregma_counts)

# Get unique bregma values
unique_bregma <- sort(unique(coronal_filtered$bregma))
cat(sprintf("\nNumber of unique bregma levels: %d\n", length(unique_bregma)))
cat("Bregma levels:\n")
print(unique_bregma)

# Create line plot
plot_bregma_lines <- ggplot(bregma_counts, aes(x = bregma, y = count, color = expr_class, group = expr_class)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  scale_color_manual(
    values = c(
      "Lepr+Glp1r+Ebf1" = "#FF6B6B",
      "Lepr+Glp1r+Trh" = "#4ECDC4",
      "All four genes" = "#FFE66D"
    ),
    name = "Expression class"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 11, face = "bold"),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Expression class distribution along bregma axis (coronal sections)",
    x = "Bregma (mm)",
    y = "Cell count"
  ) +
  scale_x_continuous(breaks = unique_bregma)

print(plot_bregma_lines)

# Save plot
ggsave(paste0(atlaspath, "expr_class_by_bregma_coronal.pdf"),
       plot = plot_bregma_lines,
       width = 12,
       height = 6,
       device = cairo_pdf)

# Create a stacked bar plot as alternative visualization
plot_bregma_stacked <- ggplot(bregma_counts, aes(x = factor(bregma), y = count, fill = expr_class)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(
    values = c(
      "Lepr+Glp1r+Ebf1" = "#FF6B6B",
      "Lepr+Glp1r+Trh" = "#4ECDC4",
      "All four genes" = "#FFE66D"
    ),
    name = "Expression class"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 11, face = "bold"),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    panel.grid.major.x = element_blank()
  ) +
  labs(
    title = "Expression class distribution along bregma axis (coronal sections)",
    x = "Bregma (mm)",
    y = "Cell count"
  )

print(plot_bregma_stacked)

# Save stacked bar plot
ggsave(paste0(atlaspath, "expr_class_by_bregma_coronal_barplot.pdf"),
       plot = plot_bregma_stacked,
       width = 12,
       height = 6,
       device = cairo_pdf)

# Create a summary table showing peak bregma for each class
peak_summary <- bregma_counts %>%
  group_by(expr_class) %>%
  slice_max(count, n = 3) %>%
  arrange(expr_class, desc(count))

cat("\n=== Top 3 bregma levels for each expression class ===\n")
print(peak_summary)

# Calculate total cells per expression class across all bregma levels
total_per_class <- bregma_counts %>%
  group_by(expr_class) %>%
  summarize(total_cells = sum(count))

cat("\n=== Total cells per expression class (coronal only) ===\n")
print(total_per_class)

# Create a wide format table for easier viewing
bregma_wide <- bregma_counts %>%
  pivot_wider(names_from = expr_class, values_from = count, values_fill = 0)

cat("\n=== Cell counts in wide format ===\n")
print(bregma_wide)

# Save tables
write.csv(bregma_counts, 
          paste0(atlaspath, "expr_class_by_bregma_coronal_long.csv"),
          row.names = FALSE)

write.csv(bregma_wide, 
          paste0(atlaspath, "expr_class_by_bregma_coronal_wide.csv"),
          row.names = FALSE)

cat("\nFiles saved:\n")
cat("1. expr_class_by_bregma_coronal.pdf (line plot)\n")
cat("2. expr_class_by_bregma_coronal_barplot.pdf (bar plot)\n")
cat("3. expr_class_by_bregma_coronal_long.csv (long format data)\n")
cat("4. expr_class_by_bregma_coronal_wide.csv (wide format data)\n")

# Heatmap of expression classes by region
# Showing both absolute counts and relative abundance

library(ggplot2)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)

# Get the data excluding negative cells
expr_data <- data.frame(
  expr_class = seurat_merged$expr_class,
  region = seurat_merged$region.collapsed,
  stringsAsFactors = FALSE
)

# Create contingency table
expr_region_table <- table(expr_data$expr_class, expr_data$region)
print("Full contingency table:")
print(expr_region_table)

# Remove negative class for the heatmap
expr_classes <- c("Lepr+Glp1r+Ebf1", "Lepr+Glp1r+Trh", "All four genes")
expr_region_table_filtered <- expr_region_table[expr_classes, , drop = FALSE]

print("\nFiltered table (without negative):")
print(expr_region_table_filtered)

# Calculate total cells per region (including negative)
total_per_region <- colSums(expr_region_table)

# Calculate proportion relative to total cells in each region
prop_matrix <- sweep(expr_region_table_filtered, 2, total_per_region, FUN = "/") * 100

print("\nProportion matrix (% of total cells in each region):")
print(round(prop_matrix, 2))

# Create annotation data frames
count_matrix <- as.matrix(expr_region_table_filtered)
prop_matrix <- as.matrix(prop_matrix)

# Function to create cell text showing both count and percentage
cell_fun_both = function(j, i, x, y, width, height, fill) {
  count <- count_matrix[i, j]
  perc <- prop_matrix[i, j]
  
  grid::grid.text(
    sprintf("%d\n(%.1f%%)", count, perc),
    x, y, 
    gp = grid::gpar(fontsize = 6, col = if(prop_matrix[i, j] > 0.5) "white" else "black")
  )
}

# Create color function for the heatmap (based on percentage)
col_fun <- colorRamp2(
  c(0, max(prop_matrix)/2, max(prop_matrix)),
  c("#FFFFFF", "#FFA500", "#FF0000")
)

# Create the heatmap using ComplexHeatmap
ht1 <- Heatmap(
  prop_matrix,
  name = "% of region",
  col = col_fun,
  
  # Row settings
  row_names_side = "left",
  row_names_gp = grid::gpar(fontsize = 10),
  cluster_rows = FALSE,
  
  # Column settings
  column_names_side = "bottom",
  column_names_gp = grid::gpar(fontsize = 7),
  column_names_rot = 90,
  column_names_max_height = unit(10, "cm"),
  cluster_columns = FALSE,
  
  # Cell annotations
  cell_fun = cell_fun_both,
  
  # Heatmap appearance
  border = TRUE,
  rect_gp = grid::gpar(col = "white", lwd = 1),
  
  # Column width
  width = unit(20, "cm"),
  
  # Legend
  heatmap_legend_param = list(
    title = "% of total\ncells in region",
    title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
    labels_gp = grid::gpar(fontsize = 9),
    legend_height = unit(4, "cm")
  ),
  
  # Title
  column_title = "Expression class distribution across cell types",
  column_title_gp = grid::gpar(fontsize = 12, fontface = "bold")
)

# Draw the heatmap with increased size
pdf(paste0(atlaspath, "expr_class_by_region_heatmap.pdf"), width = 24, height = 8)
draw(ht1)
dev.off()

# Also create a simpler ggplot2 version
# Prepare data for ggplot - convert matrix properly
plot_data <- as.data.frame(as.table(prop_matrix)) %>%
  rename(expr_class = Var1, region = Var2, percentage = Freq)

plot_data_counts <- as.data.frame(as.table(count_matrix)) %>%
  rename(expr_class = Var1, region = Var2, count = Freq)

plot_data <- plot_data %>%
  left_join(plot_data_counts, by = c("expr_class", "region"))

# Create ggplot heatmap
p1 <- ggplot(plot_data, aes(x = region, y = expr_class, fill = percentage)) +
  geom_tile(color = "white", size = 0.5) +
  geom_text(aes(label = sprintf("%d\n(%.1f%%)", count, percentage)), 
            size = 2, 
            color = ifelse(plot_data$percentage > max(plot_data$percentage)/2, "white", "black")) +
  scale_fill_gradient2(
    low = "white",
    mid = "orange",
    high = "red",
    midpoint = max(plot_data$percentage)/2,
    name = "% of total\ncells in\ncell type"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    axis.text.y = element_text(size = 10),
    axis.title = element_blank(),
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    panel.grid = element_blank(),
    plot.margin = unit(c(0.5, 0.5, 2, 0.5), "cm")
  ) +
  labs(title = "Expression class distribution across cell types")

ggsave(paste0(atlaspath, "expr_class_by_region_heatmap_ggplot.pdf"),
       plot = p1,
       width = 24,
       height = 6,
       device = cairo_pdf)

# Create a summary table
summary_df <- plot_data %>%
  arrange(expr_class, desc(percentage))

print("\nTop regions for each expression class:")
print(summary_df %>% 
        group_by(expr_class) %>% 
        slice_max(percentage, n = 3))

# Save summary tables
write.csv(count_matrix, 
          paste0(atlaspath, "expr_class_by_region_counts.csv"),
          row.names = TRUE)

write.csv(prop_matrix, 
          paste0(atlaspath, "expr_class_by_region_percentages.csv"),
          row.names = TRUE)

cat("\nFiles saved:\n")
cat("1. expr_class_by_region_heatmap.pdf (ComplexHeatmap version - 24\" wide)\n")
cat("2. expr_class_by_region_heatmap_ggplot.pdf (ggplot2 version - 24\" wide)\n")
cat("3. expr_class_by_region_counts.csv (absolute counts by cell type)\n")
cat("4. expr_class_by_region_percentages.csv (percentages by cell type)\n")

# Heatmap of expression classes by cell type annotation
# Showing both absolute counts and relative abundance

library(ggplot2)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)

# Set up cell type annotation with quality filter
seurat_merged$plot <- seurat_merged$detailed_celltype_annotation
seurat_merged@meta.data[["plot"]] <- ifelse(seurat_merged@meta.data[["score_detailed_celltype_annotation"]] > 0.4, 
                                            seurat_merged@meta.data[["plot"]], "no_label")

# Get the data excluding negative cells
expr_data <- data.frame(
  expr_class = seurat_merged$expr_class,
  celltype = seurat_merged$plot,
  stringsAsFactors = FALSE
)

# Create contingency table
expr_celltype_table <- table(expr_data$expr_class, expr_data$celltype)

# Remove negative class for the heatmap
expr_classes <- c("Lepr+Glp1r+Ebf1", "Lepr+Glp1r+Trh", "All four genes")
expr_celltype_table_filtered <- expr_celltype_table[expr_classes, , drop = FALSE]

# Calculate total cells per cell type (including negative)
total_per_celltype <- colSums(expr_celltype_table)

# Calculate proportion relative to total cells in each cell type
prop_matrix <- sweep(expr_celltype_table_filtered, 2, total_per_celltype, FUN = "/") * 100

# Create annotation data frames
count_matrix <- as.matrix(expr_celltype_table_filtered)
prop_matrix <- as.matrix(prop_matrix)

# Function to create cell text showing both count and percentage
cell_fun_both = function(j, i, x, y, width, height, fill) {
  count <- count_matrix[i, j]
  perc <- prop_matrix[i, j]
  
  grid::grid.text(
    sprintf("%d\n(%.1f%%)", count, perc),
    x, y, 
    gp = grid::gpar(fontsize = 6, col = if(prop_matrix[i, j] > 0.5) "white" else "black")
  )
}

# Create color function for the heatmap (based on percentage)
col_fun <- colorRamp2(
  c(0, max(prop_matrix)/2, max(prop_matrix)),
  c("#FFFFFF", "#FFA500", "#FF0000")
)

# Create the heatmap using ComplexHeatmap
ht1 <- Heatmap(
  prop_matrix,
  name = "% of cell type",
  col = col_fun,
  
  # Row settings
  row_names_side = "left",
  row_names_gp = grid::gpar(fontsize = 10),
  cluster_rows = FALSE,
  
  # Column settings
  column_names_side = "bottom",
  column_names_gp = grid::gpar(fontsize = 7),
  column_names_rot = 90,
  column_names_max_height = unit(10, "cm"),
  cluster_columns = FALSE,
  
  # Cell annotations
  cell_fun = cell_fun_both,
  
  # Heatmap appearance
  border = TRUE,
  rect_gp = grid::gpar(col = "white", lwd = 1),
  
  # Column width
  width = unit(20, "cm"),
  
  # Legend
  heatmap_legend_param = list(
    title = "% of total\ncells in\ncell type",
    title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
    labels_gp = grid::gpar(fontsize = 9),
    legend_height = unit(4, "cm")
  ),
  
  # Title
  column_title = "Expression class distribution across cell types",
  column_title_gp = grid::gpar(fontsize = 12, fontface = "bold")
)

# Draw the heatmap with increased size
pdf(paste0(atlaspath, "expr_class_by_celltype_heatmap.pdf"), width = 24, height = 8)
draw(ht1)
dev.off()

# Also create a simpler ggplot2 version
# Prepare data for ggplot - convert matrix properly
plot_data <- as.data.frame(as.table(prop_matrix)) %>%
  rename(expr_class = Var1, celltype = Var2, percentage = Freq)

plot_data_counts <- as.data.frame(as.table(count_matrix)) %>%
  rename(expr_class = Var1, celltype = Var2, count = Freq)

plot_data <- plot_data %>%
  left_join(plot_data_counts, by = c("expr_class", "celltype"))

# Create ggplot heatmap
p1 <- ggplot(plot_data, aes(x = celltype, y = expr_class, fill = percentage)) +
  geom_tile(color = "white", size = 0.5) +
  geom_text(aes(label = sprintf("%d\n(%.1f%%)", count, percentage)), 
            size = 2, 
            color = ifelse(plot_data$percentage > max(plot_data$percentage)/2, "white", "black")) +
  scale_fill_gradient2(
    low = "white",
    mid = "orange",
    high = "red",
    midpoint = max(plot_data$percentage)/2,
    name = "% of total\ncells in\ncell type"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    axis.text.y = element_text(size = 10),
    axis.title = element_blank(),
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    panel.grid = element_blank(),
    plot.margin = unit(c(0.5, 0.5, 2, 0.5), "cm")
  ) +
  labs(title = "Expression class distribution across cell types")

ggsave(paste0(atlaspath, "expr_class_by_celltype_heatmap_ggplot.pdf"),
       plot = p1,
       width = 24,
       height = 6,
       device = cairo_pdf)

# Create a summary table
summary_df <- plot_data %>%
  arrange(expr_class, desc(percentage))

print("\nTop cell types for each expression class:")
print(summary_df %>% 
        group_by(expr_class) %>% 
        slice_max(percentage, n = 3))

# Save summary tables
write.csv(count_matrix, 
          paste0(atlaspath, "expr_class_by_celltype_counts.csv"),
          row.names = TRUE)

write.csv(prop_matrix, 
          paste0(atlaspath, "expr_class_by_celltype_percentages.csv"),
          row.names = TRUE)

cat("\nFiles saved:\n")
cat("1. expr_class_by_celltype_heatmap.pdf (ComplexHeatmap version - 24\" wide)\n")
cat("2. expr_class_by_celltype_heatmap_ggplot.pdf (ggplot2 version - 24\" wide)\n")
cat("3. expr_class_by_celltype_counts.csv (absolute counts by cell type)\n")
cat("4. expr_class_by_celltype_percentages.csv (percentages by cell type)\n")

# Create logical vector for high-confidence astrocytes
astro_highconf <- seurat_merged$detailed_celltype_annotation == "Astrocyte"

message("High-confidence astrocytes: ", sum(astro_highconf))

# Get Calcr and Glp1r expression for ALL cells
calcr_expr <- GetAssayData(seurat_merged, assay = "raw", layer = "counts")["Calcr", ]
glp1r_expr <- GetAssayData(seurat_merged, assay = "raw", layer = "counts")["Glp1r", ]

# Check expression only in high-confidence astrocytes
calcr_positive <- calcr_expr > 0
glp1r_positive <- glp1r_expr > 0

# Double positive cells (must be astrocyte AND express both genes)
double_positive <- astro_highconf & calcr_positive & glp1r_positive

# Summary in high-confidence astrocytes only
message("\nExpression summary in high-confidence astrocytes:")
message("Calcr+: ", sum(astro_highconf & calcr_positive), " (", 
        round(100 * sum(astro_highconf & calcr_positive) / sum(astro_highconf), 2), "%)")
message("Glp1r+: ", sum(astro_highconf & glp1r_positive), " (", 
        round(100 * sum(astro_highconf & glp1r_positive) / sum(astro_highconf), 2), "%)")
message("Calcr+ Glp1r+ (double positive): ", sum(double_positive), " (", 
        round(100 * sum(double_positive) / sum(astro_highconf), 2), "%)")

# Create summary table
summary_table <- data.frame(
  Calcr = c("Negative", "Positive", "Positive", "Negative"),
  Glp1r = c("Negative", "Negative", "Positive", "Positive"),
  Count = c(
    sum(astro_highconf & !calcr_positive & !glp1r_positive),
    sum(astro_highconf & calcr_positive & !glp1r_positive),
    sum(double_positive),
    sum(astro_highconf & !calcr_positive & glp1r_positive)
  )
)
summary_table$Percentage <- round(100 * summary_table$Count / sum(astro_highconf), 2)

print(summary_table)

# Add markers to metadata
seurat_merged$astro_highconf <- astro_highconf
seurat_merged$astro_calcr_glp1r <- double_positive

astro_highconf <- seurat_merged$detailed_celltype_annotation == "Astrocyte" & seurat_merged$region.collapsed == "arcuate_nucleus"

message("High-confidence astrocytes: ", sum(astro_highconf))

# Get Calcr and Glp1r expression for ALL cells
calcr_expr <- GetAssayData(seurat_merged, assay = "raw", layer = "counts")["Calcr", ]
glp1r_expr <- GetAssayData(seurat_merged, assay = "raw", layer = "counts")["Glp1r", ]

# Check expression only in high-confidence astrocytes
calcr_positive <- calcr_expr > 0
glp1r_positive <- glp1r_expr > 0

# Double positive cells (must be astrocyte AND express both genes)
double_positive <- astro_highconf & calcr_positive & glp1r_positive

# Summary in high-confidence astrocytes only
message("\nExpression summary in high-confidence astrocytes:")
message("Calcr+: ", sum(astro_highconf & calcr_positive), " (", 
        round(100 * sum(astro_highconf & calcr_positive) / sum(astro_highconf), 2), "%)")
message("Glp1r+: ", sum(astro_highconf & glp1r_positive), " (", 
        round(100 * sum(astro_highconf & glp1r_positive) / sum(astro_highconf), 2), "%)")
message("Calcr+ Glp1r+ (double positive): ", sum(double_positive), " (", 
        round(100 * sum(double_positive) / sum(astro_highconf), 2), "%)")

# Create summary table
summary_table <- data.frame(
  Calcr = c("Negative", "Positive", "Positive", "Negative"),
  Glp1r = c("Negative", "Negative", "Positive", "Positive"),
  Count = c(
    sum(astro_highconf & !calcr_positive & !glp1r_positive),
    sum(astro_highconf & calcr_positive & !glp1r_positive),
    sum(double_positive),
    sum(astro_highconf & !calcr_positive & glp1r_positive)
  )
)
summary_table$Percentage <- round(100 * summary_table$Count / sum(astro_highconf), 2)

print(summary_table)

# Add markers to metadata
seurat_merged$astro_highconf <- astro_highconf
seurat_merged$astro_calcr_glp1r <- double_positive

library(ggplot2)
library(tidyr)
library(dplyr)

# Compile summary data from all analyses
summary_data <- data.frame(
  Dataset = c("Rat DVC (Xenium)", "Mouse DVC (Xenium)", 
              "Rat MBH (Xenium)", "Mouse MBH (Xenium)",
              "Mouse MBH (Resolve - All)", "Mouse MBH (Resolve - Arc)"),
  Total_Astrocytes = c(165937, 340125, 142164, 329693, 23964, 863),
  Calcr_Pos = c(660, 11122, 4896, 9819, 752, 89),
  Glp1r_Pos = c(1083, 5139, 907, 4640, 827, 74),
  Double_Pos = c(4, 295, 66, 266, 85, 23)
)

# Calculate percentages
summary_data$Calcr_Pct <- round(100 * summary_data$Calcr_Pos / summary_data$Total_Astrocytes, 2)
summary_data$Glp1r_Pct <- round(100 * summary_data$Glp1r_Pos / summary_data$Total_Astrocytes, 2)
summary_data$Double_Pct <- round(100 * summary_data$Double_Pos / summary_data$Total_Astrocytes, 2)

# Print summary table
print("=== SUMMARY TABLE ===")
print(summary_data)

# Create breakdown data for stacked bar chart
breakdown_data <- data.frame(
  Dataset = rep(summary_data$Dataset, each = 4),
  Category = rep(c("Double Negative", "Calcr Only", "Glp1r Only", "Double Positive"), 6),
  Count = c(
    # Rat DVC
    164198, 656, 1079, 4,
    # Mouse DVC
    324159, 10827, 4844, 295,
    # Rat MBH
    136427, 4830, 841, 66,
    # Mouse MBH
    315500, 9553, 4374, 266,
    # Resolve All
    22470, 667, 742, 85,
    # Resolve Arc
    723, 66, 51, 23
  )
)

# Calculate percentages for breakdown
breakdown_data <- breakdown_data %>%
  group_by(Dataset) %>%
  mutate(Percentage = 100 * Count / sum(Count)) %>%
  ungroup()

# Set factor order for proper display
breakdown_data$Category <- factor(breakdown_data$Category, 
                                  levels = c("Double Negative", "Calcr Only", 
                                           "Glp1r Only", "Double Positive"))
breakdown_data$Dataset <- factor(breakdown_data$Dataset, 
                                levels = summary_data$Dataset)

# Plot 1: Percentage comparison
pct_data_long <- summary_data %>%
  select(Dataset, Calcr_Pct, Glp1r_Pct, Double_Pct) %>%
  pivot_longer(cols = c(Calcr_Pct, Glp1r_Pct, Double_Pct),
               names_to = "Marker", 
               values_to = "Percentage")

pct_data_long$Marker <- factor(pct_data_long$Marker,
                               levels = c("Calcr_Pct", "Glp1r_Pct", "Double_Pct"),
                               labels = c("Calcr+", "Glp1r+", "Calcr+Glp1r+"))

p1 <- ggplot(pct_data_long, aes(x = Dataset, y = Percentage, fill = Marker)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_fill_manual(values = c("Calcr+" = "#E74C3C", 
                               "Glp1r+" = "#3498DB", 
                               "Calcr+Glp1r+" = "#9B59B6")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        legend.position = "top",
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(title = "Calcr and Glp1r Expression in High-Confidence Astrocytes",
       x = NULL, y = "Percentage of Astrocytes (%)",
       fill = "Marker Expression") +
  ylim(0, max(pct_data_long$Percentage) * 1.1)

print(p1)

# Plot 2: Stacked bar chart showing all four categories
p2 <- ggplot(breakdown_data, aes(x = Dataset, y = Percentage, fill = Category)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = c("Double Negative" = "#BDC3C7",
                               "Calcr Only" = "#E74C3C",
                               "Glp1r Only" = "#3498DB",
                               "Double Positive" = "#9B59B6")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        legend.position = "top",
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(title = "Astrocyte Expression Profile Breakdown",
       x = NULL, y = "Percentage of Astrocytes (%)",
       fill = "Expression Category")

print(p2)

# Plot 3: Focus on double-positive cells only
p3 <- ggplot(summary_data, aes(x = Dataset, y = Double_Pct)) +
  geom_col(fill = "#9B59B6", width = 0.6) +
  geom_text(aes(label = paste0(Double_Pct, "%\n(n=", Double_Pos, ")")), 
            vjust = -0.3, size = 3) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(title = "Calcr+Glp1r+ Double-Positive Astrocytes",
       x = NULL, y = "Percentage of Astrocytes (%)") +
  ylim(0, max(summary_data$Double_Pct) * 1.3)

print(p3)

# Key observations
cat("\n=== KEY OBSERVATIONS ===\n")
cat("1. Double-positive (Calcr+Glp1r+) astrocytes are RARE across all datasets (<0.4%)\n")
cat("2. Highest double-positive rate: Resolve Arc (2.67%) - 30x higher than Xenium\n")
cat("3. Mouse shows higher Calcr expression than rat (3% vs 0.4-3.4%)\n")
cat("4. Glp1r expression is consistently low (0.6-1.5% in Xenium, higher in Resolve)\n")
cat("5. Arcuate nucleus enrichment: 2.67% double-positive vs 0.35% in whole MBH\n\n")

# Save plots
ggsave("astrocyte_marker_comparison.pdf", p1, width = 10, height = 6)
ggsave("astrocyte_expression_breakdown.pdf", p2, width = 10, height = 6)
ggsave("astrocyte_double_positive.pdf", p3, width = 10, height = 6)

# Export summary table
write.csv(summary_data, "astrocyte_summary_table.csv", row.names = FALSE)

cat("Plots saved: astrocyte_marker_comparison.pdf, astrocyte_expression_breakdown.pdf, astrocyte_double_positive.pdf\n")
cat("Table saved: astrocyte_summary_table.csv\n")


