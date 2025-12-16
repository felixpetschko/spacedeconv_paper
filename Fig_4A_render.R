options(stringsAsFactors = FALSE)

library(spacedeconv)
library(SpatialExperiment)
library(ggplot2)
library(decoupleR)
library(gridExtra)

deconv <- readRDS("./export/4A_C2L.rds")

title_size <- 22
font_size <- 18
legend_size <- 20

# Note: palette names follow scico/viridis style where applicable inside spacedeconv.
truth <- plot_spatial(
    deconv,
    result = "groundTruth",
    density = FALSE,
    smooth = FALSE,
    title = "Ground truth (human)",
    show_image = FALSE,
    palette = "Inferno",
    title_size = title_size,
    font_size = font_size,
    legend_size = legend_size,
    shift_positive = FALSE
)

l5 <- plot_spatial(
    deconv,
    result = "cell2location_L5.IT",
    density = FALSE,
    smooth = TRUE,
    title = "L5.IT C2L",
    title_size = title_size,
    font_size = font_size,
    legend_size = legend_size
)

l6 <- plot_spatial(
    deconv,
    result = "cell2location_L2.3.IT",
    density = FALSE,
    smooth = TRUE,
    title = "L2.3.IT C2L",
    title_size = title_size,
    font_size = font_size,
    legend_size = legend_size
)

oligo <- plot_spatial(
    deconv,
    result = "cell2location_Oligo",
    density = FALSE,
    smooth = TRUE,
    title = "Oligodendrocytes C2L",
    title_size = title_size,
    font_size = font_size,
    legend_size = legend_size
)

astro <- plot_spatial(
    deconv,
    result = "cell2location_Astro",
    density = FALSE,
    smooth = TRUE,
    title = "Astrocytes",
    title_size = title_size,
    font_size = font_size,
    legend_size = legend_size
)

# Avoid name collision with the function `cluster()` by using clus_* variables.
clus_expr <- spacedeconv::cluster(
    deconv,
    spmethod = "expression",
    clusres = c(0.2, 0.25, 0.3, 0.4, 0.45, 0.5)
)
plot_spatial(
    clus_expr,
    result = "cluster",
    palette = "inferno",
    density = FALSE,
    shift_positive = FALSE
)

clus_kmeans <- spacedeconv::cluster(
    deconv,
    spmethod = "cell2location",
    ncluster = 7,
    method = "kmeans",
    dist_method = "euclidean"
)

kmeans_7 <- plot_spatial(
    clus_kmeans,
    result = "cluster_cell2location_nclusters_7",
    palette = "Accent",
    density = FALSE,
    title = "k-means C2L 7",
    title_size = title_size,
    font_size = font_size,
    legend_size = legend_size,
    shift_positive = FALSE
)

# Most abundant
mAbundant <- plot_most_abundant(
    deconv,
    method = "cell2location",
    min_abundance = 0.1,
    min_spot = 200,
    remove = "cell2location_Endo",
    title = "Most Abundant C2L",
    title_size = title_size,
    font_size = font_size,
    legend_size = legend_size
) +
    scale_fill_manual(
        values = colorspace::sequential_hcl(4, "viridis"),
        labels = c("Astro", "L2.3.IT", "L5.IT", "Oligo")
    )

ref <- get_decoupleR_reference(method = "collectri", split_complexes = FALSE)
res <- spacedeconv::normalize(clus_kmeans)
res <- compute_activities(res, assay = "cpm", method = "wmean", reference = ref)

ref <- get_decoupleR_reference(method = "dorothea")
res <- compute_activities(res, assay = "cpm", method = "wmean", reference = ref)

res <- spacedeconv::cluster(
    res,
    spmethod = "cell2location",
    method = "kmeans",
    ncluster = 3:8,
    dist_method = "euclidean"
)

grid_layout <- rbind(c(1, 2, 3), c(4, 5, 6))

final <- grid.arrange(
    truth, kmeans_7, mAbundant,
    oligo, l5, l6,
    layout_matrix = grid_layout
)

ggsave(
    filename = "./export/fig_4A.png",
    plot = final,
    dpi = 600,
    width = 15,
    height = 10,
    units = "in",
    bg = "white"
)
