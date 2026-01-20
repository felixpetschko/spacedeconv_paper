library(ggplot2)
library(gridExtra)

plots_3a <- readRDS("./export/fig_3A_panels.rds")
plots_3b <- readRDS("./export/fig_3B_panels.rds")

plots <- c(plots_3a, plots_3b)

add_enum <- function(p, label) {
  p + ggplot2::ggtitle(paste0(label, ") ", p$labels$title))
}

add_enum_list <- function(plots) {
  Map(function(p, i) add_enum(p, LETTERS[i]), plots, seq_along(plots))
}

plots <- add_enum_list(plots)

layout_matrix <- rbind(
  c(1, 2, 3, 4),
  c(5, 6, 7, 8),
  c(9, 10, 11, 12),
  c(13, 14, 15, 16),
  c(17, 17, 18, 18)
)

final_plot <- grid.arrange(
  grobs = plots,
  layout_matrix = layout_matrix
)

ggplot2::ggsave(
  filename = "./export/fig_3.png",
  plot = final_plot,
  width = 20,
  height = 25,
  units = "in",
  dpi = 600,
  bg = "white"
)
