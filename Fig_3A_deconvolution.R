library(spacedeconv)
library(SpatialExperiment)
library(ggplot2)
library(gridExtra)
library(dplyr)

spe <- read10xVisium("./data/sudmeier/750/")

rownames(spe) <- rowData(spe)$symbol
spe <- preprocess(spe, min_umi = 87) # remove low count spots
spe <- spacedeconv::normalize(spe)

metadata <- read.csv("./data/IEDB_VDJDB_table_anno.tsv", sep = "\t")

# filter only TCR clonotype
metadata <- dplyr::filter(metadata, grepl("TRB", metadata$Clonoytpe_ids))
metadata$barcode <- paste0(metadata$barcode, "-1")

metadata$Clonoytpe_ids <- make.names(metadata$Clonoytpe_ids)
metadata <- metadata[metadata$Sample == "16", ] # patient 16

# add metadata to object
for (clonotype in unique(metadata$Clonoytpe_ids)) {
  tmp <- metadata[metadata$Clonoytpe_ids == clonotype, ]
  spe <- annotate_spots(spe, tmp$barcode, name = clonotype)
}

deconvEstimate <- deconvolute(spe, method = "estimate", assay_sp = "cpm")
saveRDS(deconvEstimate, file = "./export/3AEstimate.rds")

deconvEPIC <- deconvolute(spe, method = "epic", assay_sp = "cpm", tumor = TRUE)
saveRDS(deconvEPIC, file = "./export/3AEpic.rds")

deconvQuanTiseq <- deconvolute(
  spe,
  method = "quantiseq",
  assay_sp = "cpm",
  tumor = TRUE
)

deconvQuanTiseq <- aggregate_results(
  deconvQuanTiseq,
  cell_type_1 = "quantiseq_T.cell.CD4...non.regulatory.",
  cell_type_2 = "quantiseq_T.cell.regulatory..Tregs.",
  name = "aggCD4"
)

deconvQuanTiseq <- aggregate_results(
  deconvQuanTiseq,
  cell_type_1 = "aggCD4",
  cell_type_2 = "quantiseq_T.cell.CD8.",
  name = "aggT"
)

saveRDS(deconvQuanTiseq, file = "./export/3AQuantiseq.rds")