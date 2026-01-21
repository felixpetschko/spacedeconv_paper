library(spacedeconv)
library(ggplot2)
library(cowplot)
figures_dir <- "./export/figures"
objects_dir <- "./export/objects"
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(objects_dir, showWarnings = FALSE, recursive = TRUE)


allresults_minor_andersson <- readRDS("./data/allresults_minor_andersson.rds")

title_size <- 22
font_size <- 18
legend_size <- 20

# B cells naive (C2L)
p1 <- spacedeconv::plot_spatial(
  allresults_minor_andersson,
  result = "c2l_B.cells.Naive",
  density = FALSE,
  smooth = FALSE,
  title_size = title_size,
  font_size = font_size,
  legend_size = legend_size,
  title = "Naïve B cells",
  nDigits = 2
)

# B cells naive (C2L) sqrt
p2 <- spacedeconv::plot_spatial(
  allresults_minor_andersson,
  result = "c2l_B.cells.Naive",
  density = FALSE,
  smooth = FALSE,
  transform_scale = "sqrt",
  title_size = title_size,
  font_size = font_size,
  legend_size = legend_size,
  title = "Naïve B cells (sqrt)",
  nDigits = 2
)

# B cells naive (C2L) log2
p3 <- spacedeconv::plot_spatial(
  allresults_minor_andersson,
  result = "c2l_B.cells.Naive",
  density = FALSE,
  smooth = FALSE,
  transform_scale = "log2",
  title_size = title_size,
  font_size = font_size,
  legend_size = legend_size,
  title = "Naïve B cells (log2)",
  nDigits = 2,
  pseudocount = 0
)

# B cells naive smoothed
p4 <- spacedeconv::plot_spatial(
  allresults_minor_andersson,
  result = "c2l_B.cells.Naive",
  density = FALSE,
  smooth = TRUE,
  title_size = title_size,
  font_size = font_size,
  legend_size = legend_size,
  title = "Naïve B cells (smoothed)",
  nDigits = 2
)

add_enum <- function(p, label) {
  p + ggplot2::ggtitle(paste0(label, ") ", p$labels$title))
}

add_enum_list <- function(plots, labels) {
  Map(add_enum, plots, labels)
}

plots <- add_enum_list(list(p1, p2, p3, p4), LETTERS[1:4])

# Arrange the plots in a grid
final_plot <- plot_grid(
  plotlist = plots,
  ncol = 4,
  align = "hv",
  axis = "tblr"
)

# Save the UMAP as a high-resolution PNG
ggplot2::ggsave(
  filename = "./export/figures/fig_S1.png",
  plot = final_plot,
  width = 20,
  height = 5,
  units = "in",
  dpi = 600,
  bg = "white"
)
