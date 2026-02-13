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

# iCAF (C2L)
p5 <- spacedeconv::plot_spatial(
  allresults_minor_andersson,
  result = "c2l_CAFs.MSC.iCAF.like",
  density = FALSE,
  smooth = TRUE,
  title_size = title_size,
  font_size = font_size,
  legend_size = legend_size,
  title = "iCAF (sm.)",
  nDigits = 2
)

# myCAF (C2L)
p6 <- spacedeconv::plot_spatial(
  allresults_minor_andersson,
  result = "c2l_CAFs.myCAF.like",
  density = FALSE,
  smooth = TRUE,
  title_size = title_size,
  font_size = font_size,
  legend_size = legend_size,
  title = "myCAF (sm.)",
  nDigits = 2
)

# iCAF vs. myCAF (C2L)
p7 <- spacedeconv::plot_comparison(
  allresults_minor_andersson,
  cell_type_1 = "c2l_CAFs.MSC.iCAF.like",
  cell_type_2 = "c2l_CAFs.myCAF.like",
  density = FALSE,
  smooth = TRUE,
  title_size = title_size,
  font_size = font_size,
  legend_size = legend_size,
  title = "iCAF vs. myCAF (sm.)",
  palette = "Purple-Green",
  reverse_palette = TRUE,
  shift_positive = FALSE,
  nDigits = 2
)

# RCTD most abundant cancer molecular subtype
p8 <- spacedeconv::plot_spatial(
  allresults_minor_andersson,
  result = "rctd_Cancer.LumB.SC",
  density = FALSE,
  smooth = TRUE,
  title_size = title_size,
  font_size = font_size,
  legend_size = legend_size,
  title = "LumB (sm.)",
  nDigits = 2
)

add_enum <- function(p, label) {
  p +
    labs(tag = paste0(label)) +
    theme(
      plot.title = element_text(hjust = 0.5, vjust = 1),
      plot.title.position = "plot",
      plot.tag = element_text(hjust = 0, vjust = 1, face = "bold", size = 24),
      plot.tag.position = c(0, 1),
      plot.margin = margin(6, 6, 6, 6)
    )
}

add_enum_list <- function(plots, labels) {
  Map(add_enum, plots, labels)
}

plots <- add_enum_list(list(p5, p6, p7, p8), LETTERS[1:4])

# Arrange the plots in a grid
final_plot <- plot_grid(
  plotlist = plots,
  ncol = 2,
  align = "hv",
  axis = "tblr"
)

# Save the UMAP as a high-resolution PNG
ggplot2::ggsave(
  filename = "./export/figures/fig_2.png",
  plot = final_plot,
  width = 10,
  height = 10,
  units = "in",
  dpi = 600,
  bg = "white"
)
