# PAPER_6A
# Author: Constantin Zackl
# Date: Sys.Date()

# ---- Libraries ----
suppressPackageStartupMessages({
  library(spacedeconv)
  library(SpatialExperiment)
  library(readr)
  library(progress)
  library(gridExtra)
  library(ggplot2)
})

# ---- Load Visium data ----
spe <- read10xVisium("./data/sudmeier/750/")
rownames(spe) <- rowData(spe)$symbol
spe <- preprocess(spe, min_umi = 87)
spe <- spacedeconv::normalize(spe)

# ---- Build comparison stats from metadata ----
metadata <- readr::read_tsv("./data/TCR/IEDB_VDJDB_table_anno.tsv")
metadata$barcode <- paste0(metadata$barcode, "-1")
metadata <- metadata[metadata$Sample_x == "16", ]

pb <- progress::progress_bar$new(total = length(colnames(spe)))

TCRstat <- data.frame(row.names = colnames(spe))

for (spot in colnames(spe)) {
  tmp <- metadata[metadata$barcode == spot, ]

  # IE Organism
  unknownIE <- tmp[tmp$Category_IEDB == "Unknown", ]
  otherIE <- tmp[tmp$Category_IEDB == "viral", ]
  sumUnknownIE <- sum(unknownIE$UMIs)
  sumOther <- sum(otherIE$UMIs)
  TCRstat[spot, "IEUnknown"] <- sumUnknownIE
  TCRstat[spot, "IEViral"] <- sumOther

  # VDJ
  unknownVDJ <- tmp[tmp$Category_VDJDB == "Unknown", ]
  otherVDJ <- tmp[tmp$Category_VDJDB == "viral", ]
  sumUnknownVDJ <- sum(unknownVDJ$UMIs)
  sumOtherVDJ <- sum(otherVDJ$UMIs)
  TCRstat[spot, "VDJUnknown"] <- sumUnknownVDJ
  TCRstat[spot, "VDJViral"] <- sumOtherVDJ

  pb$tick()
}

colData(spe) <- cbind(colData(spe), TCRstat)
pb$terminate()

# ---- Plot comparisons ----
ie <- plot_comparison(
  spe,
  cell_type_1 = "IEUnknown",
  cell_type_2 = "IEViral",
  density = FALSE,
  smooth = FALSE,
  title = "IE Unknown vs Viral",
  palette = "Vik",
  reverse_palette = TRUE,
  title_size = 25,
  spot_size = 1.03
)

vdj <- plot_comparison(
  spe,
  cell_type_1 = "VDJUnknown",
  cell_type_2 = "VDJViral",
  density = FALSE,
  smooth = FALSE,
  title = "VDJUnknown vs Viral",
  palette = "Vik",
  title_size = 25,
  spot_size = 1.03,
  reverse_palette = TRUE
)

# ---- Load tcr object and plots ----
tcr <- readRDS("./export/paper/3Btcr.rds")

tcr$UMIB0 <- colData(tcr)$UMITCR > 0
tcr$UMIB <- colData(tcr)$UMITCR >= 5

tcr <- preprocess(tcr, min_umi = 87)

tcr5 <- plot_celltype(
  tcr,
  cell_type = "UMIB",
  density = FALSE,
  title = "TCR UMI > 5",
  smooth = FALSE,
  show_image = TRUE,
  palette = "Oslo",
  show_legend = FALSE,
  zoom = TRUE,
  palette_type = "discrete",
  image_id = "hires",
  title_size = 25,
  spot_size = 1.03
)

tcr0 <- plot_celltype(
  tcr,
  cell_type = "UMIB0",
  density = FALSE,
  title = "TCR All",
  smooth = FALSE,
  show_image = TRUE,
  palette = "Oslo",
  show_legend = FALSE,
  zoom = TRUE,
  palette_type = "discrete",
  image_id = "hires",
  title_size = 25,
  spot_size = 1.03
)

tcrUMI <- plot_celltype(
  tcr,
  cell_type = "UMITCR",
  density = FALSE,
  title = "TCR UMI",
  smooth = FALSE,
  show_image = FALSE,
  show_legend = TRUE,
  zoom = TRUE,
  image_id = "hires",
  title_size = 25,
  spot_size = 1.03
)

# ---- Arrange and save plot ----
grid_layout <- rbind(c(1, 2, 3))
final <- grid.arrange(tcrUMI, tcr5, ie, layout_matrix = grid_layout)
ggsave(filename = "./export2/fig_6A.png",
       plot = final, dpi = 300, width = 15, height = 7, units = "in")