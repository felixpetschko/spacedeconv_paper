library(spacedeconv)
library(SpatialExperiment)
library(readr)
library(progress)
library(gridExtra)
library(ggplot2)

spe <- read10xVisium("./data/sudmeier/750/")
rownames(spe) <- rowData(spe)$symbol
spe <- preprocess(spe, min_umi = 87)
spe <- spacedeconv::normalize(spe)

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

title_size <- 22
font_size <- 18
legend_size <- 20

ie <- plot_comparison(
    spe,
    cell_type_1 = "IEUnknown",
    cell_type_2 = "IEViral",
    density = FALSE,
    smooth = FALSE,
    title = "IE Unknown vs Viral",
    palette = "Vik",
    reverse_palette = TRUE,
    title_size = title_size,
    font_size = font_size,
    legend_size = legend_size,
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
    title_size = title_size,
    font_size = font_size,
    legend_size = legend_size,
    spot_size = 1.03,
    reverse_palette = TRUE
)

tcr <- readRDS("./export/3Btcr.rds")

tcr$UMIB0 <- colData(tcr)$UMITCR > 0
tcr$UMIB <- colData(tcr)$UMITCR >= 5

tcr <- preprocess(tcr, min_umi = 87)

tcr5 <- plot_spatial(
    tcr,
    result = "UMIB",
    density = FALSE,
    title = "TCR UMI > 5",
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

tcr0 <- plot_spatial(
    tcr,
    result = "UMIB0",
    density = FALSE,
    title = "TCR All",
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

tcrUMI <- plot_spatial(
    tcr,
    result = "UMITCR",
    density = FALSE,
    title = "TCR UMI",
    smooth = FALSE,
    show_image = FALSE,
    show_legend = TRUE,
    zoom = TRUE,
    image_id = "hires",
    title_size = title_size,
    font_size = font_size,
    legend_size = legend_size,
    spot_size = 1.03
)

grid_layout <- rbind(c(1, 2, 3))

final <- grid.arrange(tcrUMI, tcr5, ie, layout_matrix = grid_layout)

ggsave(
    filename = "./export/fig_6A.png",
    plot = final,
    dpi = 600,
    width = 15,
    height = 5,
    units = "in",
    bg = "white"
)
