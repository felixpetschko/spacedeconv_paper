# Figure 6B
# Author: Constantin Zackl
# Date: Sys.Date()

# ---- Libraries ----
suppressPackageStartupMessages({
  library(spacedeconv)
  library(SpatialExperiment)
  library(gridExtra)
  library(ggplot2)
})

# ---- Load Visium data ----
spe <- read10xVisium("./data/renal/frozen_b1/")

# ---- Add pathology annotation ----
annot <- read.csv("./data/renal/frozen_b1/outs/spatial/GSM5924046_frozen_b_1_TLS_annotation.csv")
spe$tls <- annot$TLS_2_cat

# ---- Preprocess & normalize ----
spe <- preprocess(spe)
spe <- spacedeconv::normalize(spe)
rownames(spe) <- rowData(spe)$symbol

# ---- Quick pathology plot ----
plot_celltype(
  spe,
  cell_type = "tls",
  offset_rotation = TRUE,
  density = FALSE,
  legend_size = 50,
  font_size = 30,
  palette = "BluYl",
  title = "TLS"
)

# ---- Deconvolution (quanTIseq) ----
deconv <- deconvolute(
  spe,
  method = "quantiseq",
  assay_sp = "cpm",
  tumor = TRUE
)

# ---- Optional alt. deconv methods (disabled, mirror Rmd) ----
# deconv <- deconvolute(deconv, method = "mcp_counter", assay_sp = "cpm")

# ref <- get_decoupleR_reference("progeny")
# deconv <- compute_activities(deconv, ref)
# ref <- get_decoupleR_reference("collectri")
# deconv <- compute_activities(deconv, ref)

# ---- TLS-related gene set & score ----
geneset <- c(
  "IGHA1","IGHG1","IGHG2","IGHG3","IGHG4","IGHGP","IGHM",
  "IGKG","IGLC1","IGLC2","IGLC3","JCHAIN","CD79A","FCRL5",
  "MZB1","SSR4","XBP1","TRBC2","IL7R","CXCL12","LUM",
  "C1QA","C7","CD52","APOE","PTLP","PTGSD","PIM2","DERL3"
)

deconv <- gene_set_score(deconv, genes = geneset, assay = "cpm")

# ---- Optional clustering (disabled, mirror Rmd) ----
# deconv <- cluster(deconv, spmethod = "expression", clusres = 0.3)
# cluster_plot <- plot_celltype(
#   deconv, density = FALSE,
#   cell_type = "cluster_expression_res_0.3",
#   offset_rotation = TRUE, title = "Cluster 0.3"
# )

# ---- Plots to combine ----
patho <- plot_celltype(
  deconv, "tls",
  density = FALSE, palette = "BluYl",
  offset_rotation = TRUE,
  title = "Pathology annotation",
  title_size = 25, spot_size = 1.03
)

geneSet <- plot_celltype(
  deconv, "geneSet",
  density = FALSE, offset_rotation = TRUE,
  title = "TLS score",
  title_size = 25, spot_size = 1.03
)

qBcells <- plot_celltype(
  deconv, "quantiseq_B.cell",
  density = FALSE, offset_rotation = TRUE, smooth = TRUE,
  title = "quanTIseq B cells smoothed",
  title_size = 25, spot_size = 1.03
)

# ---- Arrange and save figure ----
grid_layout <- rbind(c(1, 2, 3))
final <- grid.arrange(patho, geneSet, qBcells, layout_matrix = grid_layout)

# Match the Rmd chunk fig.width/fig.height (14.8 x 5 inches)
ggsave(
  filename = "./export2/fig_6B.png",
  plot = final, dpi = 300, width = 14.8, height = 5, units = "in"
)
