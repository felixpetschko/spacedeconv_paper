# ---- setup ----
knitr::opts_chunk$set(echo = TRUE)

library(spacedeconv)
library(SpatialExperiment)

allresults_minor_1142243F <- readRDS("./data/allresults_minor_1142243F.rds")
allresults_minor_1160920F <- readRDS("./data/allresults_minor_1160920F.rds")
allresults_minor_4465 <- readRDS("./data/allresults_minor_4465.rds")
allresults_minor_44971 <- readRDS("./data/allresults_minor_44971.rds")

# lumA
allresults_minor_4290 <- readRDS("./data/allresults_minor_4290.rds")
# lumB
allresults_minor_4535 <- readRDS("./data/allresults_minor_4535.rds")


# ---- plot_celltype calls ----
# most abundant cancer cell types, with density
x1 <- plot_celltype(allresults_minor_1142243F, "rctd_Cancer.Basal.SC", density = FALSE,
                    title_size = 5, font_size = 25, legend_size = 30,
                    smooth = TRUE, limits = c(0,1), title = "")
x2 <- plot_celltype(allresults_minor_1142243F, "rctd_Cancer.LumA.SC", density = FALSE,
                    title_size = 5, font_size = 25, legend_size = 30,
                    smooth = TRUE, limits = c(0,1), title = "")
x3 <- plot_celltype(allresults_minor_1142243F, "rctd_Cancer.LumB.SC", density = FALSE,
                    title_size = 5, font_size = 25, legend_size = 30,
                    smooth = TRUE, limits = c(0,1), title = "")
x4 <- plot_celltype(allresults_minor_1142243F, "rctd_Cancer.Her2.SC", density = FALSE,
                    title_size = 5, font_size = 25, legend_size = 30,
                    smooth = TRUE, limits = c(0,1), title = "")
# ======================================================================================
x5 <- plot_celltype(allresults_minor_1160920F, "rctd_Cancer.Basal.SC", density = FALSE,
                    title_size = 5, font_size = 25, legend_size = 30,
                    smooth = TRUE, limits = c(0,1), title = "") 
x6 <- plot_celltype(allresults_minor_1160920F, "rctd_Cancer.LumA.SC", density = FALSE,
                    title_size = 5, font_size = 25, legend_size = 30,
                    smooth = TRUE, limits = c(0,1), title = "") 
x7 <- plot_celltype(allresults_minor_1160920F, "rctd_Cancer.LumB.SC", density = FALSE,
                    title_size = 5, font_size = 25, legend_size = 30,
                    smooth = TRUE, limits = c(0,1), title = "") 
x8 <- plot_celltype(allresults_minor_1160920F, "rctd_Cancer.Her2.SC", density = FALSE,
                    title_size = 5, font_size = 25, legend_size = 30,
                    smooth = TRUE, limits = c(0,1), title = "") 
# ======================================================================================
x9 <- plot_celltype(allresults_minor_4465, "rctd_Cancer.Basal.SC", density = FALSE,
                    title_size = 5, font_size = 25, legend_size = 30,
                    smooth = TRUE, limits = c(0,1), title = "") 
x10 <- plot_celltype(allresults_minor_4465, "rctd_Cancer.LumA.SC", density = FALSE,
                     title_size = 5, font_size = 25, legend_size = 30,
                     smooth = TRUE, limits = c(0,1), title = "") 
x11 <- plot_celltype(allresults_minor_4465, "rctd_Cancer.LumB.SC", density = FALSE,
                     title_size = 5, font_size = 25, legend_size = 30,
                     smooth = TRUE, limits = c(0,1), title = "") 
x12 <- plot_celltype(allresults_minor_4465, "rctd_Cancer.Her2.SC", density = FALSE,
                     title_size = 5, font_size = 25, legend_size = 30,
                     smooth = TRUE, limits = c(0,1), title = "") 
# =====================================================================================
x13 <- plot_celltype(allresults_minor_44971, "rctd_Cancer.Basal.SC", density = FALSE,
                     title_size = 5, font_size = 25, legend_size = 30,
                     smooth = TRUE, limits = c(0,1), title = "") 
x14 <- plot_celltype(allresults_minor_44971, "rctd_Cancer.LumA.SC", density = FALSE,
                     title_size = 5, font_size = 25, legend_size = 30,
                     smooth = TRUE, limits = c(0,1), title = "") 
x15 <- plot_celltype(allresults_minor_44971, "rctd_Cancer.LumB.SC", density = FALSE,
                     title_size = 5, font_size = 25, legend_size = 30,
                     smooth = TRUE, limits = c(0,1), title = "") 
x16 <- plot_celltype(allresults_minor_44971, "rctd_Cancer.Her2.SC", density = FALSE,
                     title_size = 5, font_size = 25, legend_size = 30,
                     smooth = TRUE, limits = c(0,1), title = "") 
# ===================================================================================
x17 <- plot_celltype(allresults_minor_4290, "rctd_Cancer.Basal.SC", density = FALSE,
                     title_size = 5, font_size = 25, legend_size = 30,
                     smooth = TRUE, limits = c(0,1), title = "")
x18 <- plot_celltype(allresults_minor_4290, "rctd_Cancer.LumA.SC", density = FALSE,
                     title_size = 5, font_size = 25, legend_size = 30,
                     smooth = TRUE, limits = c(0,1), title = "")
x19 <- plot_celltype(allresults_minor_4290, "rctd_Cancer.LumB.SC", density = FALSE,
                     title_size = 5, font_size = 25, legend_size = 30,
                     smooth = TRUE, limits = c(0,1), title = "")
x20 <- plot_celltype(allresults_minor_4290, "rctd_Cancer.Her2.SC", density = FALSE,
                     title_size = 5, font_size = 25, legend_size = 30,
                     smooth = TRUE, limits = c(0,1), title = "")
# ===================================================================================
x21 <- plot_celltype(allresults_minor_4535, "rctd_Cancer.Basal.SC", density = FALSE,
                     title_size = 5, font_size = 25, legend_size = 30,
                     smooth = TRUE, limits = c(0,1), title = "")
x22 <- plot_celltype(allresults_minor_4535, "rctd_Cancer.LumA.SC", density = FALSE,
                     title_size = 5, font_size = 25, legend_size = 30,
                     smooth = TRUE, limits = c(0,1), title = "")
x23 <- plot_celltype(allresults_minor_4535, "rctd_Cancer.LumB.SC", density = FALSE,
                     title_size = 5, font_size = 25, legend_size = 30,
                     smooth = TRUE, limits = c(0,1), title = "")
x24 <- plot_celltype(allresults_minor_4535, "rctd_Cancer.Her2.SC", density = FALSE,
                     title_size = 5, font_size = 25, legend_size = 30,
                     smooth = TRUE, limits = c(0,1), title = "")


# ---- grid arrangement and saving ----
library(gridExtra)
library(ggplot2)
library(grid)
library(ggpubr)

# Create a grid layout
grid_layout <- rbind(
  c(1, 2, 3, 4),
  c(5, 6, 7, 8),
  c(9, 10, 11, 12),
  c(13, 14, 15, 16),
  c(17, 18, 19, 20),
  c(21, 22, 23, 24)
)

# Arrange the plots in a grid
grid_1 <- grid.arrange(
  x1, x2, x3, x4,
  x5, x6, x7, x8,
  x9, x10, x11, x12,
  x13, x14, x15, x16,
  x17, x18, x19, x20,
  x21, x22, x23, x24,
  layout_matrix = grid_layout
)

ggsave("./export2/suppl1_300dpi_extra.png", plot = grid_1, dpi = 300,
       width = 28, height = 42, units = "in")

ggsave("./export2/suppl1_300dpi_extra.pdf", plot = grid_1, dpi = 300,
       width = 28, height = 42, units = "in")