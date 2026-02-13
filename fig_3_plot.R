library(ggplot2)
library(gridExtra)
figures_dir <- "./export/figures"
objects_dir <- "./export/objects"
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(objects_dir, showWarnings = FALSE, recursive = TRUE)


plots <- readRDS("./export/objects/fig_3_panels.rds")
plots <- plots[-5]

add_enum <- function(p, label) {
  p +
    labs(tag = paste0(label)) +
    theme(
      plot.title = element_text(hjust = 0.5, vjust = 1),
      plot.title.position = "plot",
      plot.tag = element_text(hjust = 0, vjust = 1, face = "bold", size = 24),
      plot.tag.position = c(0, 1),
      plot.margin = margin(0, 0, 0, 0)
    )
}

add_enum_list <- function(plots) {
  Map(function(p, i) add_enum(p, LETTERS[i]), plots, seq_along(plots))
}

plots <- add_enum_list(plots)
plots[[1]] <- plots[[1]] + theme(plot.margin = margin(0, 250, 0, 0))

grid_layout <- rbind(
    c(1, 1, 1, 1),
    c(2, 2, 3, 3),
    c(4, 4, 5, 5),
    c(6, 6, 7, 7)
)

final <- grid.arrange(
    grobs = plots,
    layout_matrix = grid_layout
)

ggsave(
    filename = "./export/figures/fig_3.png",
    plot = final,
    dpi = 600,
    width = 10,
    height = 20,
    units = "in",
    bg = "white"
)
