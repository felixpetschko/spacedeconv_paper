library(spacedeconv)
library(SpatialExperiment)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)

decouple <- readRDS("./export/3BdecoupleR.rds")
tcr <- readRDS("./export/3Btcr.rds")

title_size <- 22
font_size <- 18
legend_size <- 20

tgfb <- plot_spatial(
  decouple,
  result = "progeny_TGFb",
  title = "TGFb",
  density = FALSE,
  smooth = TRUE,
  title_size = title_size,
  font_size = font_size,
  legend_size = legend_size,
  spot_size = 1.03
)

wnt <- plot_spatial(
  decouple,
  result = "progeny_WNT",
  title = "WNT",
  density = FALSE,
  smooth = TRUE,
  title_size = title_size,
  font_size = font_size,
  legend_size = legend_size,
  spot_size = 1.03
)

hypox <- plot_spatial(
  decouple,
  result = "progeny_Hypoxia",
  title = "Hypoxia",
  density = FALSE,
  smooth = TRUE,
  title_size = title_size,
  font_size = font_size,
  legend_size = legend_size,
  spot_size = 1.03
)

myc <- plot_spatial(
  decouple,
  result = "collectri_MYC",
  title = "MYC",
  density = FALSE,
  smooth = TRUE,
  title_size = title_size,
  font_size = font_size,
  legend_size = legend_size,
  spot_size = 1.03
)

jak <- plot_spatial(
  decouple,
  result = "progeny_JAK.STAT",
  title = "JAK-STAT",
  density = FALSE,
  smooth = TRUE,
  title_size = title_size,
  font_size = font_size,
  legend_size = legend_size,
  spot_size = 1.03
)

stat <- plot_spatial(
  decouple,
  result = "collectri_STAT1",
  title = "STAT1",
  density = FALSE,
  smooth = TRUE,
  title_size = title_size,
  font_size = font_size,
  legend_size = legend_size,
  spot_size = 1.03
)

tcr$UMIB0 <- colData(tcr)$UMITCR > 0
tcr$UMIB <- colData(tcr)$UMITCR >= 5
tcr <- preprocess(tcr, min_umi = 87)

tcr5 <- plot_spatial(
  tcr,
  result = "UMIB",
  density = FALSE,
  title = "TCR",
  smooth = FALSE,
  show_image = TRUE,
  palette = "Oslo",
  show_legend = FALSE,
  zoom = TRUE,
  palette_type = "discrete",
  image_id = "hires",
  title_size = title_size,
  font_size = font_size,
  legend_size = legend_size,
  spot_size = 1.03
)

spe <- read10xVisium("./data/sudmeier/750/")
rownames(spe) <- rowData(spe)$symbol
spe <- subsetSPE(spe, colRange = c(0, 10000))
spe <- spacedeconv::normalize(spe)

deconvEPIC <- readRDS("./export/3Bepic.rds")
deconvQuanTIseq <- readRDS("./export/3Bquantiseq.rds")

assays(deconvEPIC)[["RNA"]] <- assays(deconvEPIC)[["counts"]]
assays(deconvQuanTIseq)[["RNA"]] <- assays(deconvQuanTIseq)[["counts"]]

deconvQuanTIseq <- deconvQuanTIseq[!duplicated(rownames(deconvQuanTIseq)), ]

cluster1 <- spacedeconv::cluster(
  deconvQuanTIseq,
  data = "expression",
  clusres = 0.1
)

cluster2 <- spacedeconv::cluster(
  deconvQuanTIseq,
  data = "expression",
  clusres = 0.2
)

cluster3 <- spacedeconv::cluster(
  deconvQuanTIseq,
  data = "expression",
  clusres = 0.3
)

clusterX <- spacedeconv::cluster(deconvEPIC, spmethod = "epic", nclusters = 3)

cl1 <- "cluster_expression_res_0.1"
cl2 <- "cluster_expression_res_0.2"
cl3 <- "cluster_expression_res_0.3"
clX <- "cluster_epic_nclusters_3"

SummarizedExperiment::colData(cluster1)[[cl1]] <- as.character(SummarizedExperiment::colData(cluster1)[[cl1]])
SummarizedExperiment::colData(cluster2)[[cl2]] <- as.character(SummarizedExperiment::colData(cluster2)[[cl2]])
SummarizedExperiment::colData(cluster3)[[cl3]] <- as.character(SummarizedExperiment::colData(cluster3)[[cl3]])
SummarizedExperiment::colData(clusterX)[[clX]] <- as.character(SummarizedExperiment::colData(clusterX)[[clX]])

pclus1 <- plot_spatial(
  cluster1,
  result = "cluster",
  palette = "Accent",
  title = "Clustering 0.1",
  density = FALSE,
  title_size = title_size,
  font_size = font_size,
  legend_size = legend_size,
  spot_size = 1.03,
  image_id = NULL,
  show_image = FALSE,
  palette_type = "discrete"
)

pclus2 <- plot_spatial(
  cluster2,
  result = "cluster",
  palette = "inferno",
  title = "Clustering 0.2",
  density = FALSE,
  title_size = title_size,
  font_size = font_size,
  legend_size = legend_size,
  spot_size = 1.03,
  image_id = NULL,
  show_image = FALSE
)

pclus3 <- plot_spatial(
  cluster3,
  result = "cluster",
  palette = "Accent",
  title = "Clustering 0.3",
  density = FALSE,
  title_size = title_size,
  font_size = font_size,
  legend_size = legend_size,
  spot_size = 1.03,
  image_id = NULL,
  show_image = FALSE
)

plot_spatial(
  clusterX,
  result = "cluster",
  palette = "inferno",
  title = "Clustering EPIC",
  density = FALSE,
  title_size = title_size,
  font_size = font_size,
  legend_size = legend_size,
  spot_size = 1.03,
  image_id = NULL,
  show_image = FALSE
)

mAbundantSTD <- plot_most_abundant(
  deconvEPIC,
  method = "epic",
  palette = "inferno",
  title_size = title_size,
  font_size = font_size,
  legend_size = legend_size,
  min_abundance = 0.05,
  title = "Most Abundant",
  spot_size = 1.03
)

mAbundantNoCancer <- plot_most_abundant(
  deconvEPIC,
  method = "epic",
  palette = "inferno",
  title_size = title_size,
  font_size = font_size,
  legend_size = legend_size,
  remove = "epic_uncharacterized.cell",
  min_abundance = 0.05,
  title = "Without Tumor",
  spot_size = 1.03
)

mAbundantSTD <- mAbundantSTD +
  ggplot2::scale_fill_manual(
    values = RColorBrewer::brewer.pal(7, "Accent"),
    labels = c(
      "B cells",
      "CAFs",
      "Endothelial",
      "Macrophage",
      "CD4",
      "CD8",
      "Tumor"
    )
  )

mAbundantNoCancer <- mAbundantNoCancer +
  ggplot2::scale_fill_manual(
    values = c(RColorBrewer::brewer.pal(6, "Accent"), "#D3D3D3"),
    labels = c(
      "B cells",
      "CAFs",
      "Endothelial",
      "Macrophage",
      "CD4",
      "CD8",
      "Undefined"
    )
  )

spe <- get_lr(spe, resource = "Consensus", method = "min", organism = "human")

col <- plot_spatial(
  spe,
  result = "lr_COL18A1.ITGB1",
  density = FALSE,
  smooth = TRUE,
  title = "COL18A1-ITGB1",
  title_size = title_size,
  font_size = font_size,
  legend_size = legend_size,
  spot_size = 1.03
)

grid_layout <- rbind(
  c(1, 2, 3, 4),
  c(5, 6, 7, 8),
  c(9, 9, 10, 10)
)

tmp <- grid.arrange(
  wnt, hypox, myc, tgfb,
  jak, col, pclus1, pclus3,
  mAbundantSTD, mAbundantNoCancer,
  layout_matrix = grid_layout
)

ggplot2::ggsave(
  filename = "./export/fig_3B.png",
  plot = tmp,
  dpi = 600,
  width = 20,
  height = 15,
  units = "in",
  bg = "white"
)