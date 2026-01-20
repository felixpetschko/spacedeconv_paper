library(ggplot2)
library(gridExtra)

plots <- readRDS("./export/fig_5A_panels.rds")

grid_layout <- rbind(c(1, 2, 3))
final <- grid.arrange(grobs = plots, layout_matrix = grid_layout)

ggsave(
    filename = "./export/fig_5A.png",
    plot = final,
    dpi = 600,
    width = 15,
    height = 5,
    units = "in",
    bg = "white"
)
