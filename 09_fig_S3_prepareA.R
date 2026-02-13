library(spacedeconv)
library(SpatialExperiment)
library(ggplot2)
figures_dir <- "./export/figures"
objects_dir <- "./export/objects"
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(objects_dir, showWarnings = FALSE, recursive = TRUE)


deconvEstimate <- readRDS("./export/objects/fig_S3A_estimate.rds")
deconvEPIC <- readRDS("./export/objects/fig_S3A_epic.rds")
deconvQuanTiseq <- readRDS("./export/objects/fig_S3A_quantiseq.rds")

title_size <- 22
font_size <- 18
legend_size <- 20

# ESTIMATE panels
estimate_tumor <- plot_spatial(
  deconvEstimate,
  result = "estimate_tumor.purity",
  smooth = TRUE,
  density = FALSE,
  title = "Tumor Purity",
  title_size = title_size,
  font_size = font_size,
  legend_size = legend_size
)

estimate_immune <- plot_spatial(
  deconvEstimate,
  result = "estimate_immune.score",
  smooth = TRUE,
  density = FALSE,
  title = "Immune Score",
  title_size = title_size,
  font_size = font_size,
  legend_size = legend_size
)

estimate_stroma <- plot_spatial(
  deconvEstimate,
  result = "estimate_stroma.score",
  smooth = TRUE,
  density = FALSE,
  title = "Stroma Score",
  title_size = title_size,
  font_size = font_size,
  legend_size = legend_size
)

# EPIC panels
epic_caf <- plot_spatial(
  deconvEPIC,
  result = "epic_Cancer.associated.fibroblast",
  smooth = TRUE,
  density = FALSE,
  title = "CAFs",
  title_size = title_size,
  font_size = font_size,
  legend_size = legend_size
)

epic_b <- plot_spatial(
  deconvEPIC,
  result = "epic_B.cell",
  smooth = TRUE,
  density = FALSE,
  title = "B cells",
  title_size = title_size,
  font_size = font_size,
  legend_size = legend_size
)

epic_endothelial <- plot_spatial(
  deconvEPIC,
  result = "epic_Endothelial.cell",
  smooth = TRUE,
  density = FALSE,
  title = "Endothelial cell",
  title_size = title_size,
  font_size = font_size,
  legend_size = legend_size
)

# quanTIseq panels
quantiseq_cd8 <- plot_spatial(
  deconvQuanTiseq,
  result = "quantiseq_T.cell.CD8.",
  smooth = TRUE,
  density = FALSE,
  title = "CD8 T cells",
  title_size = title_size,
  font_size = font_size,
  legend_size = legend_size,
  shift_positive = FALSE
)

quantiseq_t <- plot_spatial(
  deconvQuanTiseq,
  result = "aggT",
  smooth = TRUE,
  density = FALSE,
  title = "T cells",
  title_size = title_size,
  font_size = font_size,
  legend_size = legend_size,
  shift_positive = FALSE
)

plots_3a <- list(
  estimate_tumor,
  estimate_immune,
  estimate_stroma,
  epic_caf,
  epic_b,
  epic_endothelial,
  quantiseq_t,
  quantiseq_cd8
)

saveRDS(plots_3a, file = "./export/objects/fig_S3A_panels.rds")
