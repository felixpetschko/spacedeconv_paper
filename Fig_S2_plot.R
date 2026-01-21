library(ggplot2)
library(gridExtra)

plots <- readRDS("./export/fig_S2_panels.rds")

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
