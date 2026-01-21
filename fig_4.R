library(SpatialExperiment)
library(jsonlite)
library(purrr)
library(pbapply)
library(ggplot2)
library(spacedeconv)
library(gridExtra)
library(sf)
figures_dir <- "./export/figures"
objects_dir <- "./export/objects"
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(objects_dir, showWarnings = FALSE, recursive = TRUE)


a1 <- read.table("./data/V1_Section1_detections.txt", sep = "\t", header = TRUE)
spe <- read10xVisium("./data/breastCancer/Section1/outs/")

coords <- spatialCoords(spe)

scalefactors <- jsonlite::fromJSON(
    "./data/breastCancer/Section1/outs/spatial/scalefactors_json.json"
)

radius <- scalefactors$spot_diameter_fullres / 2

calculate_distances <- function(spot_x, spot_y, cell_xs, cell_ys) {
    sqrt((spot_x - cell_xs)^2 + (spot_y - cell_ys)^2)
}

count_cells_within_radius <- function(spot_row, a1, radius) {
    spot_x <- spot_row$pxl_col_in_fullres
    spot_y <- spot_row$pxl_row_in_fullres
    distances <- calculate_distances(
        spot_x,
        spot_y,
        a1$Centroid.X.px,
        a1$Centroid.Y.px
    )
    sum(distances <= radius)
}

coords_df <- as.data.frame(coords)

spe$cell_density <- unlist(
    pblapply(
        seq_len(nrow(coords_df)),
        function(i) count_cells_within_radius(coords_df[i, ], a1, radius),
        cl = 1
    )
)

deconv <- readRDS("./data/allresults_minor_andersson.rds")
deconv$cell_counts <- spe$cell_density

for (result in available_results(deconv, "rctd")) {
    deconv <- scale_cell_counts(deconv, result, "cell_counts")
}

smooth <- FALSE

cutout <- subsetSPE(deconv, colRange = c(11200, 22000), rowRange = c(0, 15300))

title_size <- 22
font_size <- 18
legend_size <- 20

abs <- plot_spatial(
    spe,
    result = "cell_density",
    density = FALSE,
    title = "Cell counts",
    smooth = smooth,
    title_size = title_size,
    font_size = font_size,
    legend_size = legend_size,
)

cut <- plot_spatial(
    cutout,
    result = "cell_counts",
    density = FALSE,
    title = "Cell counts ROI",
    smooth = smooth,
    title_size = title_size,
    font_size = font_size,
    legend_size = legend_size,
)

xmin <- 270
xmax <- 510
ymin <- -380
ymax <- -90

rect_df <- data.frame(
    x = c(xmin, xmax, xmax, xmin, xmin),
    y = c(ymin, ymin, ymax, ymax, ymin)
)

abs <- abs +
    geom_polygon(
        data = rect_df,
        aes(x = x, y = y),
        fill = NA,
        color = "#FF5733",
        size = 2,
        linetype = "dashed"
    ) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank())

lumb <- plot_spatial(
    cutout,
    result = "rctd_Cancer.LumB.SC",
    density = FALSE,
    title = "Cancer (LumB)",
    smooth = smooth,
    title_size = title_size,
    font_size = font_size,
    legend_size = legend_size,
)

lumbabs <- plot_spatial(
    cutout,
    result = "rctd_Cancer.LumB.SC_absolute",
    density = FALSE,
    title = "Cancer (LumB) abs.",
    smooth = smooth,
    title_size = title_size,
    font_size = font_size,
    legend_size = legend_size,
)

caf <- plot_spatial(
    cutout,
    result = "rctd_CAFs.myCAF.like",
    density = FALSE,
    title = "myCAFs",
    smooth = smooth,
    title_size = title_size,
    font_size = font_size,
    legend_size = legend_size,
)

cafabs <- plot_spatial(
    cutout,
    result = "rctd_CAFs.myCAF.like_absolute",
    density = FALSE,
    title = "myCAFs abs.",
    smooth = smooth,
    title_size = title_size,
    font_size = font_size,
    legend_size = legend_size,
)

add_enum <- function(p, label) {
    p +
        labs(tag = paste0(label)) +
        theme(
            plot.title = element_text(hjust = 0.5),
            plot.title.position = "plot",
            plot.tag = element_text(hjust = 0, vjust = 1, face = "bold", size = 24),
            plot.tag.position = c(0, 1)
        )
}

add_enum_list <- function(plots, labels) {
    Map(add_enum, plots, labels)
}

plots <- add_enum_list(
    list(abs, cut, lumb, lumbabs, caf, cafabs),
    LETTERS[1:6]
)

grid_layout <- rbind(
    c(1, 2),
    c(3, 4),
    c(5, 6)
)

final <- grid.arrange(
    grobs = plots,
    layout_matrix = grid_layout
)

ggsave(
    filename = "./export/figures/fig_4.png",
    plot = final,
    dpi = 600,
    width = 10,
    height = 15,
    units = "in",
    bg = "white"
)
