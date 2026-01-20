library(ggplot2)
library(gridExtra)

plots <- readRDS("./export/fig_5A_panels.rds")

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

plots <- add_enum_list(plots, LETTERS[1:length(plots)])

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
