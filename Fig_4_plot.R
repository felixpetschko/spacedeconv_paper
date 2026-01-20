library(ggplot2)
library(gridExtra)

plots <- readRDS("./export/fig_4_panels.rds")

add_enum <- function(p, label) {
  p + ggplot2::ggtitle(paste0(label, ") ", p$labels$title))
}

add_enum_list <- function(plots) {
  Map(function(p, i) add_enum(p, LETTERS[i]), plots, seq_along(plots))
}

plots <- add_enum_list(plots)

grid_layout <- rbind(c(1, 2, 3), c(4, 5, 6))

final <- grid.arrange(
  grobs = plots,
  layout_matrix = grid_layout
)

ggplot2::ggsave(
  filename = "./export/fig_4.png",
  plot = final,
  dpi = 600,
  width = 15,
  height = 10,
  units = "in",
  bg = "white"
)
