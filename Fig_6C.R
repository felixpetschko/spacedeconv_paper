# DensityAggregationSection1
# Date: 2024-02-13

# ---- Libraries ----
suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(jsonlite)
  library(purrr)
  library(pbapply)
  library(ggplot2)
  library(spacedeconv)
  library(gridExtra)
  library(sf)
})

# ---- Inputs ----
a1 <- read.table("./data/V1_Section1_detections.txt", sep = "\t", header = TRUE)

# ---- Load Visium data ----
spe <- read10xVisium("./data/breastCancer/Section1/outs/")
print(spe)

# ---- Spot coordinates (fullres px) ----
coords <- spatialCoords(spe)

# ---- Diagnostics ----
coords_df <- as.data.frame(coords)
cat(paste0("Max X Cell Centroids: ", max(a1$Centroid.X.px), "\n"))
cat(paste0("Max Col Spot Coordinates: ", max(coords_df[, "pxl_col_in_fullres"]), "\n"))
cat(paste0("Max Y: ", max(a1$Centroid.Y.px), "\n"))
cat(paste0("Max Row Spot Coordinates: ", max(coords_df[, "pxl_row_in_fullres"]), "\n"))
cat(paste0("Number of cells: ", nrow(a1), "\n"))

# ---- Read Visium scalefactors & radius ----
scalefactors <- jsonlite::fromJSON("./data/breastCancer/Section1/outs/spatial/scalefactors_json.json")
radius <- scalefactors$spot_diameter_fullres / 2

# ---- Helpers ----
calculate_distances <- function(spot_x, spot_y, cell_xs, cell_ys) {
  sqrt((spot_x - cell_xs)^2 + (spot_y - cell_ys)^2)
}

count_cells_within_radius <- function(spot_row, a1, radius) {
  spot_x <- spot_row$pxl_col_in_fullres
  spot_y <- spot_row$pxl_row_in_fullres
  distances <- calculate_distances(spot_x, spot_y, a1$Centroid.X.px, a1$Centroid.Y.px)
  sum(distances <= radius)
}

# ---- Compute per-spot cell density ----
spe$cell_density <- unlist(
  pblapply(seq_len(nrow(coords_df)),
           function(i) count_cells_within_radius(coords_df[i, ], a1, radius),
           cl = 1)
)

# ---- QC plots ----
plot_umi_count(spe)
plot_celltype(spe, cell_type = "cell_density")
cat(paste0("Mean Number of cells per spot: ", mean(spe$cell_density), "\n"))

# ---- Load deconvolution results & attach counts ----
deconv <- readRDS("./data/allresults_minor_andersson.rds")
deconv$cell_counts <- spe$cell_density

# Scale cell counts for all available rctd results
for (result in available_results(deconv, "rctd")) {
  deconv <- scale_cell_counts(deconv, result, "cell_counts")
}

# ---- Plot settings ----
smooth <- FALSE

# ---- Subset a rectangular region (by fullres pixel ranges) ----
cutout <- subsetSPE(deconv, colRange = c(11200, 22000), rowRange = c(0, 15300))

# ---- Absolute & subset cell count plots ----
abs <- plot_celltype(spe, "cell_density",
                     density = FALSE, title = "Cell counts",
                     smooth = smooth, title_size = 25)
cut <- plot_celltype(cutout, "cell_counts",
                     density = FALSE, title = "Cell counts subset",
                     smooth = smooth, title_size = 25)

# ---- Draw a dashed rectangle on the full view ----
xmin <- 270; xmax <- 510; ymin <- -380; ymax <- -90
rect_df <- data.frame(
  x = c(xmin, xmax, xmax, xmin, xmin),
  y = c(ymin, ymin, ymax, ymax, ymin)
)

abs <- abs +
  geom_polygon(data = rect_df, aes(x = x, y = y),
               fill = NA, color = "#FF5733", size = 2, linetype = "dashed") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())

print(abs)

# ---- Selected cell-type maps (absolute & relative) ----
lumb <- plot_celltype(cutout, "rctd_Cancer.LumB.SC",
                      density = FALSE, title = "Cancer (LumB)",
                      smooth = smooth, title_size = 25)

lumbabs <- plot_celltype(cutout, "rctd_Cancer.LumB.SC_absolute",
                         density = FALSE, title = "Cancer (LumB) absolute",
                         smooth = smooth, title_size = 25)

caf <- plot_celltype(cutout, "rctd_CAFs.myCAF.like",
                     density = FALSE, title = "myCAFs",
                     smooth = smooth, title_size = 25)

cafabs <- plot_celltype(cutout, "rctd_CAFs.myCAF.like_absolute",
                        density = FALSE, title = "myCAFs absolute",
                        smooth = smooth, title_size = 25)

# ---- Arrange and save figure (2x3 grid) ----
grid_layout <- rbind(c(1, 2, 3),
                     c(4, 5, 6))
final <- grid.arrange(abs, lumb, caf, cut, lumbabs, cafabs,
                      layout_matrix = grid_layout)

ggsave(filename = "./export2/fig_6C.png",
       plot = final, dpi = 300, width = 14, height = 9, units = "in")
