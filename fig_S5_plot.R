library(ggplot2)
library(cowplot)
figures_dir <- "./export/figures"
objects_dir <- "./export/objects"
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(objects_dir, showWarnings = FALSE, recursive = TRUE)


plots <- readRDS("./export/objects/fig_S5_panels.rds")

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

final <- plot_grid(plotlist = plots, ncol = 3, align = "hv", axis = "tblr")

ggsave(
    filename = "./export/figures/fig_S5.png",
    plot = final,
    dpi = 600,
    width = 15,
    height = 5,
    units = "in",
    bg = "white"
)
