library(spacedeconv)
library(SpatialExperiment)

deconvEstimate <- readRDS("./export/3AEstimate.rds")
deconvEPIC <- readRDS("./export/3AEpic.rds")
deconvQuanTiseq <- readRDS("./export/3AQuantiseq.rds")

title_size <- 25

# ESTIMATE panels
estimate_tumor <- plot_celltype(
    deconvEstimate,
    cell_type = "estimate_tumor.purity",
    smooth = TRUE,
    density = FALSE,
    title = "ESTIMATE Tumor Purity",
    title_size = title_size
)

estimate_immune <- plot_celltype(
    deconvEstimate,
    cell_type = "estimate_immune.score",
    smooth = TRUE,
    density = FALSE,
    title = "ESTIMATE Immune Score",
    title_size = title_size
)

estimate_stroma <- plot_celltype(
    deconvEstimate,
    cell_type = "estimate_stroma.score",
    smooth = TRUE,
    density = FALSE,
    title = "ESTIMATE Stroma Score",
    title_size = title_size
)

# EPIC panels
epic_caf <- plot_celltype(
    deconvEPIC,
    cell_type = "epic_Cancer.associated.fibroblast",
    smooth = TRUE,
    density = FALSE,
    title = "EPIC CAFs",
    title_size = title_size
)

epic_b <- plot_celltype(
    deconvEPIC,
    cell_type = "epic_B.cell",
    smooth = TRUE,
    density = FALSE,
    title = "EPIC B cells",
    title_size = title_size
)

epic_endothelial <- plot_celltype(
    deconvEPIC,
    cell_type = "epic_Endothelial.cell",
    smooth = TRUE,
    density = FALSE,
    title = "EPIC Endothelial cell",
    title_size = title_size
)

# quanTIseq panels
quantiseq_cd8 <- plot_celltype(
    deconvQuanTiseq,
    cell_type = "quantiseq_T.cell.CD8.",
    smooth = TRUE,
    density = FALSE,
    title = "quanTIseq CD8 T cells",
    title_size = title_size,
    shift_positive = FALSE
)

quantiseq_cd4 <- plot_celltype(
    deconvQuanTiseq,
    cell_type = "aggCD4",
    smooth = TRUE,
    density = FALSE,
    title = "quanTIseq CD4+",
    title_size = title_size,
    shift_positive = FALSE
)

quantiseq_t <- plot_celltype(
    deconvQuanTiseq,
    cell_type = "aggT",
    smooth = TRUE,
    density = FALSE,
    title = "quanTIseq T cells",
    title_size = title_size,
    shift_positive = FALSE
)

grid_layout <- rbind(c(1, 2, 3, 4),
                     c(5, 6, 7, 8))

final_plot <- grid.arrange(
  estimate_tumor, estimate_immune, estimate_stroma, epic_caf,
  epic_b, epic_endothelial, quantiseq_t, quantiseq_cd8,
  layout_matrix = grid_layout
)

ggsave(filename = "./export/3A.png",
    plot = final_plot,
    dpi = 600,
    width = 19,
    height = 10,
    units = "in",
    bg = "white"
)