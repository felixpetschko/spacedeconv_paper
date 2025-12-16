library(spacedeconv)
library(SpatialExperiment)
library(decoupleR)
library(dplyr)

spe <- read10xVisium("./data/sudmeier/750/")

rownames(spe) <- rowData(spe)$symbol
spe <- preprocess(spe, min_umi = 87) # remove low count cells
spe <- spacedeconv::normalize(spe)

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

ref <- get_decoupleR_reference("progeny")
deconv <- spacedeconv::compute_activities(spe, ref)

ref <- get_decoupleR_reference("collectri")
deconv <- spacedeconv::compute_activities(deconv, ref)

saveRDS(deconv, file = "./export/3BdecoupleR.rds")


speTCR <- read10xVisium("./data/sudmeier/750/", images = "hires")
rownames(speTCR) <- rowData(speTCR)$symbol

for (spot in colnames(speTCR)) {
  tmp <- metadata[metadata$barcode == spot, ]
  umi <- sum(tmp$UMIs)
  colData(speTCR)[spot, "UMITCR"] <- umi
}

library(DelayedArray)
library(S4Vectors)

for (nm in assayNames(speTCR)) {
  assay(speTCR, nm) <- as.matrix(assay(speTCR, nm))
}

if ("imgData" %in% slotNames(speTCR)) {
  imgData(speTCR) <- DataFrame()
}

saveRDS(speTCR, file = "./export/3Btcr.rds")


speTCR$UMIB <- colData(speTCR)$UMITCR >= 5

spe <- read10xVisium("./data/sudmeier/750/")
rownames(spe) <- rowData(spe)$symbol
spe <- preprocess(spe, min_umi = 87)
spe <- spacedeconv::normalize(spe)

deconvEPIC <- deconvolute(spe, method = "epic", assay_sp = "cpm", tumor = TRUE)
deconvQuanTIseq <- deconvolute(
  spe,
  method = "quantiseq",
  assay_sp = "cpm",
  tumor = TRUE
)

saveRDS(deconvEPIC, file = "./export/3Bepic.rds")
saveRDS(deconvQuanTIseq, file = "./export/3Bquantiseq.rds")
