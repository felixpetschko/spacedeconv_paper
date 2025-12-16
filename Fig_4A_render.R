# ---- Setup ----
options(stringsAsFactors = FALSE)

# ---- Packages ----
suppressPackageStartupMessages({
  library(spacedeconv)
  library(SpatialExperiment)
  library(ggplot2)
  library(decoupleR)
  library(gridExtra)
})

# ---- Load Data ----
deconv <- readRDS("./export2/paper/PAPER_4A_C2L_3.rds")

# ---- Selected Cell-type Maps ----
# Note: palette names follow scico/viridis style where applicable inside spacedeconv.
truth <- plot_celltype(deconv, cell_type = "groundTruth", density = FALSE, smooth = FALSE,
                       title = "Ground truth (human)", show_image = FALSE, palette = "Inferno", title_size = 25, shift_positive=FALSE)

l5    <- plot_celltype(deconv, "cell2location_L5.IT",   density = FALSE, smooth = TRUE, title = "L5.IT C2L",       title_size = 25)
l6    <- plot_celltype(deconv, "cell2location_L2.3.IT", density = FALSE, smooth = TRUE, title = "L2.3.IT C2L",     title_size = 25)
oligo <- plot_celltype(deconv, "cell2location_Oligo",    density = FALSE, smooth = TRUE, title = "Oligodendrocytes C2L", title_size = 25)
astro <- plot_celltype(deconv, "cell2location_Astro",    density = FALSE, smooth = TRUE, title = "Astrocytes",     title_size = 25)

# ---- Clustering & Most Abundant ----
# Avoid name collision with the function `cluster()` by using clus_* variables.
clus_expr <- spacedeconv::cluster(deconv, spmethod = "expression",
                                  clusres = c(0.2, 0.25, 0.3, 0.4, 0.45, 0.5))
plot_celltype(clus_expr, "cluster", palette = "inferno", density = FALSE, shift_positive=FALSE)

clus_kmeans <- spacedeconv::cluster(deconv, spmethod = "cell2location",
                                    ncluster = 7, method = "kmeans", dist_method = "euclidean")
kmeans_7 <- plot_celltype(clus_kmeans, "cluster_cell2location_nclusters_7", palette = "Accent",
                          density = FALSE, title = "k-means C2L 7", title_size = 25, shift_positive=FALSE)
# kmeans_8 <- plot_celltype(clus_kmeans, "cluster_cell2location_nclusters_8", palette = "Accent",
#                           density = FALSE, title = "k-means C2L 8", title_size = 25)

# Most abundant
mAbundant <- plot_most_abundant(deconv, method = "cell2location", min_abundance = 0.1, min_spot = 200,
                                remove = "cell2location_Endo", title = "Most Abundant C2L", title_size = 25) +
  scale_fill_manual(values = colorspace::sequential_hcl(4, "viridis"),
                    labels = c("Astro", "L2.3.IT", "L5.IT", "Oligo"))

# ---- decoupleR activities on clustered object ----
ref <- get_decoupleR_reference(method = "collectri", split_complexes = FALSE)
res <- spacedeconv::normalize(clus_kmeans)
res <- compute_activities(res, assay = "cpm", method = "wmean", reference = ref)

ref <- get_decoupleR_reference(method = "dorothea")
res <- compute_activities(res, assay = "cpm", method = "wmean", reference = ref)

res <- spacedeconv::cluster(res, spmethod = "cell2location", method = "kmeans",
                            ncluster = 3:8, dist_method = "euclidean")

# ---- Assemble and Save Figure ----
grid_layout <- rbind(c(1, 2, 3),
                     c(4, 5, 6))

final <- grid.arrange(truth, kmeans_7, mAbundant, oligo, l5, l6, layout_matrix = grid_layout)

# Ensure output directory exists
out_dir <- "./export2/spacedeconvPlots"
try(dir.create(out_dir, recursive = TRUE), silent = TRUE)

ggsave(filename = file.path(out_dir, "4A.png"), plot = final, dpi = 300, width = 16, height = 10, units = "in")
cat(sprintf("Saved figure to %s/4A.png\n", out_dir))