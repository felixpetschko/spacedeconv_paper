library(spacedeconv)
library(gridExtra)
library(ggplot2)
library(grid)

allresults_minor_andersson <- readRDS("./data/allresults_minor_andersson.rds")

# B cells naive (C2L)
p1 <- spacedeconv::plot_celltype(
  allresults_minor_andersson,
  cell_type = "c2l_B.cells.Naive",
  density = FALSE, smooth = FALSE,
  title_size = 40, font_size = 25, legend_size = 30,
  title = "C2L Naïve B cells",
  nDigits = 2
)

# B cells naive (C2L) sqrt
p2 <- spacedeconv::plot_celltype(
  allresults_minor_andersson,
  cell_type = "c2l_B.cells.Naive",
  density = FALSE, smooth = FALSE, transform_scale = "sqrt",
  title_size = 40, font_size = 25, legend_size = 30,
  title = "C2L Naïve B cells (sqrt)",
  nDigits = 2
)

# B cells naive (C2L) log2
p3 <- spacedeconv::plot_celltype(
  allresults_minor_andersson,
  cell_type = "c2l_B.cells.Naive",
  density = FALSE, smooth = FALSE, transform_scale = "log2",
  title_size = 40, font_size = 25, legend_size = 30,
  title = "C2L Naïve B cells (log2)",
  nDigits = 2, pseudocount = 0 # added pseudo_count = 0
)

# B cells naive smoothed
p4 <- spacedeconv::plot_celltype(
  allresults_minor_andersson,
  cell_type = "c2l_B.cells.Naive",
  density = FALSE, smooth = TRUE,
  title_size = 40, font_size = 25, legend_size = 30,
  title = "C2L Naïve B cells (smoothed)",
  nDigits = 2
)

# iCAF (C2L)
p5 <- spacedeconv::plot_celltype(
  allresults_minor_andersson,
  cell_type = "c2l_CAFs.MSC.iCAF.like",
  density = FALSE, smooth = TRUE,
  title_size = 40, font_size = 25, legend_size = 30,
  title = "C2L iCAF (smoothed)",
  nDigits = 2
)

# myCAF (C2L)
p6 <- spacedeconv::plot_celltype(
  allresults_minor_andersson,
  cell_type = "c2l_CAFs.myCAF.like",
  density = FALSE, smooth = TRUE,
  title_size = 40, font_size = 25, legend_size = 30,
  title = "C2L myCAF (smoothed)",
  nDigits = 2
)

# iCAF vs. myCAF (C2L)
p7 <- spacedeconv::plot_comparison(
  allresults_minor_andersson,
  cell_type_1 = "c2l_CAFs.MSC.iCAF.like",
  cell_type_2 = "c2l_CAFs.myCAF.like",
  density = FALSE, smooth = TRUE,
  title_size = 40, font_size = 25, legend_size = 30,
  title = "C2L iCAF vs. myCAF (smoothed)",
  palette = "Purple-Green", reverse_palette = TRUE,
  shift_positive = FALSE,
  nDigits = 2
)

# RCTD most abundant cancer molecular subtype
p8 <- spacedeconv::plot_celltype(
  allresults_minor_andersson,
  cell_type = "rctd_Cancer.LumB.SC",
  density = FALSE, smooth = TRUE,
  title_size = 40, font_size = 25, legend_size = 30,
  title = "RCTD LumB (smoothed)",
  nDigits = 2
)

# Create a grid layout
grid_layout <- rbind(
  c(1, 2, 3, 4),
  c(5, 6, 7, 8)
)

# Arrange the plots in a grid
grid_2 <- grid.arrange(
  p1, p2, p3, p4,
  p5, p6, p7, p8,
  layout_matrix = grid_layout
)

# Save the UMAP as a high-resolution PNG
ggplot2::ggsave(
  filename = "./export/fig_2B.png",
  plot = grid_2,
  width = 32, height = 17, units = "in",
  dpi = 300,
  limitsize = FALSE
)