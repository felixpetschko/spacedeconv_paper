library(ggplot2)
library(gridExtra)

plots <- readRDS("./export/fig_S2_panels.rds")

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

add_enum_list <- function(plots) {
  Map(function(p, i) add_enum(p, LETTERS[i]), plots, seq_along(plots))
}

plots <- add_enum_list(plots)

grid_layout <- rbind(
    c(1, 1, 1, 2, 2, 2),
    c(3, 3, 4, 4, 5, 5),
    c(6, 6, 7, 7, 8, 8)
)

final <- grid.arrange(
    grobs = plots,
    layout_matrix = grid_layout
)

ggsave(
    filename = "./export/fig_S2.png",
    plot = final,
    dpi = 600,
    width = 15,
    height = 15,
    units = "in",
    bg = "white"
)
