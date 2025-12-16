library(spacedeconv)
library(SpatialExperiment)
library(gridExtra)
library(ggplot2)

spe <- read10xVisium("./data/renal/frozen_b1/")

annot <- read.csv(
    "./data/renal/frozen_b1/outs/spatial/GSM5924046_frozen_b_1_TLS_annotation.csv"
)

spe$tls <- annot$TLS_2_cat

spe <- preprocess(spe)
spe <- spacedeconv::normalize(spe)
rownames(spe) <- rowData(spe)$symbol

plot_spatial(
    spe,
    result = "tls",
    offset_rotation = TRUE,
    density = FALSE,
    legend_size = 50,
    font_size = 30,
    palette = "BluYl",
    title = "TLS"
)

deconv <- deconvolute(
    spe,
    method = "quantiseq",
    assay_sp = "cpm",
    tumor = TRUE
)

geneset <- c(
  "IGHA1","IGHG1","IGHG2","IGHG3","IGHG4","IGHGP","IGHM",
  "IGKG","IGLC1","IGLC2","IGLC3","JCHAIN","CD79A","FCRL5",
  "MZB1","SSR4","XBP1","TRBC2","IL7R","CXCL12","LUM",
  "C1QA","C7","CD52","APOE","PTLP","PTGSD","PIM2","DERL3"
)

deconv <- gene_set_score(deconv, genes = geneset, assay = "cpm")

patho <- plot_spatial(
    deconv,
    result = "tls",
    density = FALSE,
    palette = "BluYl",
    offset_rotation = TRUE,
    title = "Pathology annotation",
    title_size = 25,
    spot_size = 1.03
)

geneSet <- plot_spatial(
    deconv,
    result = "geneSet",
    density = FALSE,
    offset_rotation = TRUE,
    title = "TLS score",
    title_size = 25,
    spot_size = 1.03
)

qBcells <- plot_spatial(
    deconv,
    result = "quantiseq_B.cell",
    density = FALSE,
    offset_rotation = TRUE,
    smooth = TRUE,
    title = "quanTIseq B cells smoothed",
    title_size = 25,
    spot_size = 1.03
)

grid_layout <- rbind(c(1, 2, 3))
final <- grid.arrange(patho, geneSet, qBcells, layout_matrix = grid_layout)

ggsave(
    filename = "./export/fig_6B.png",
    plot = final,
    dpi = 600,
    width = 15,
    height = 5,
    units = "in",
    bg = "white"
)
