# Fig 3B Deconvolution
# Author: Constantin Zackl
# Date: Sys.Date()

# ---- Libraries ----
suppressPackageStartupMessages({
  library(spacedeconv)
  library(SpatialExperiment)
  library(decoupleR)
  library(dplyr)
})

# ---- Load Visium data ----
spe <- read10xVisium("./data/sudmeier/750/")

# ---- Preprocessing ----
rownames(spe) <- rowData(spe)$symbol
spe <- preprocess(spe, min_umi = 87) # remove low count cells
spe <- spacedeconv::normalize(spe)

# ---- TCR metadata ----
metadata <- read.csv("./data/TCR/Spatial_final_table.tsv", sep = "\t")

# filter only TCR clonotype
metadata <- dplyr::filter(metadata, grepl("TRB", metadata$Clonoytpe_ids))
metadata$barcode <- paste0(metadata$barcode, "-1")
metadata$Clonoytpe_ids <- make.names(metadata$Clonoytpe_ids)
metadata <- metadata[metadata$Sample == "16", ]

# add metadata to object
for (clonotype in unique(metadata$Clonoytpe_ids)) {
  tmp <- metadata[metadata$Clonoytpe_ids == clonotype, ]
  spe <- annotate_spots(spe, tmp$barcode, name = clonotype)
}

# ---- decoupleR activities ----
ref <- get_decoupleR_reference("progeny")
deconv <- spacedeconv::compute_activities(spe, ref)

ref <- get_decoupleR_reference("collectri")
deconv <- spacedeconv::compute_activities(deconv, ref)

saveRDS(deconv, file = "./export2/paper/3BdecoupleR.rds")

# ---- Plots for activities ----
plot_celltype(deconv, cell_type = "progeny_TGFb", smooth = TRUE, density = FALSE, title = "TGFb")
plot_celltype(deconv, cell_type = "progeny_WNT", smooth = TRUE, density = FALSE, title = "WNT")

plot_celltype(deconv, cell_type = "collectri_BATF", smooth = TRUE, density = FALSE, title = "BATF")
plot_celltype(deconv, cell_type = "collectri_EOMES", smooth = TRUE, density = FALSE, title = "EOMES")

# ---- TCR UMI plots ----
speTCR <- read10xVisium("./data/sudmeier/750/", images = "hires")
rownames(speTCR) <- rowData(speTCR)$symbol

for (spot in colnames(speTCR)) {
  tmp <- metadata[metadata$barcode == spot, ]
  umi <- sum(tmp$UMIs)
  colData(speTCR)[spot, "UMITCR"] <- umi
}

library(DelayedArray)
library(S4Vectors)

# Alles aus HDF5/DelayedArray in den RAM laden
for (nm in assayNames(speTCR)) {
  assay(speTCR, nm) <- as.matrix(assay(speTCR, nm))
}

# Option: Bilddaten entfernen (oft groß und OOM)
if ("imgData" %in% slotNames(speTCR)) {
  imgData(speTCR) <- DataFrame()
}

# Jetzt normales Speichern möglich
saveRDS(speTCR, file = "./export2/paper/3Btcr.rds")

###

plot_celltype(speTCR, cell_type = "UMITCR", density = FALSE, title = "TCR UMI", smooth = TRUE, show_image = FALSE, image_id = NULL)

speTCR$UMIB <- colData(speTCR)$UMITCR >= 5
plot_celltype(speTCR, cell_type = "UMIB", density = FALSE, title = "TCR UMI >= 5",
              smooth = FALSE, show_image = TRUE, palette = "Oslo", image_id = NULL)

# ---- Repeat preprocessing for deconvolution ----
spe <- read10xVisium("./data/sudmeier/750/")
rownames(spe) <- rowData(spe)$symbol
spe <- preprocess(spe, min_umi = 87)
spe <- spacedeconv::normalize(spe)

# ---- Deconvolution (EPIC & quanTIseq) ----
deconvEPIC <- deconvolute(spe, method = "epic", assay_sp = "cpm", tumor = TRUE)
deconvQuanTIseq <- deconvolute(spe, method = "quantiseq", assay_sp = "cpm", tumor = TRUE)

saveRDS(deconvEPIC, file = "./export2/paper/3Bepic.rds")
saveRDS(deconvQuanTIseq, file = "./export2/paper/3Bquantiseq.rds")