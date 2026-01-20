library(ggplot2)
library(cowplot)

plots <- readRDS("./export/fig_4_panels.rds")

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

final <- plot_grid(
  plotlist = plots,
  ncol = 3,
  align = "hv",
  axis = "tblr"
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
