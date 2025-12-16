# PAPER 4B deconvolution
# Author: Constantin Zackl
# Date: Automatically generated on: 
print(Sys.Date())

# ---- Load Packages ----
library(spacedeconv)
library(SpatialExperiment)

# ---- Load Spatial Data ----
spe = read10xVisium("./data/cell2location/visium/48/")
rownames(spe) <- rowData(spe)$symbol
spe = preprocess(spe)
spe = spacedeconv::normalize(spe)
rownames(spe) <- make.names(rownames(spe), unique = TRUE)

# ---- Load Reference Data ----
sce = readRDS("./data/cell2location/sce.rds")
rownames(sce) <- make.names(rownames(sce), unique = TRUE)
sce = sce[, !grepl("Unk", sce$annotation_1)]
sce = sce[, !grepl("LowQ", sce$annotation_1)]
sce = subsetSCE(sce, ncells = 20000, cell_type_col = "annotation_1")

# ---- Deconvolution with Cell2location ----
sce = spacedeconv::normalize(sce)
signature = build_model(
  sce,
  cell_type_col = "annotation_1",
  method = "cell2location",
  epochs = 250,
  gpu = TRUE,
  batch_id_col = "sample"
)

deconv = deconvolute(
  spe,
  method = "cell2location",
  signature = signature,
  epochs = 30000,
  gpu = TRUE,
  values = "relative"
)

saveRDS(deconv, file = "./export2/paper/PAPER_4B_C2L_3.rds")

# ---- Deconvolution with RCTD ----
deconv = deconvolute(
  spe,
  method = "rctd",
  single_cell_obj = sce,
  cell_type_col = "annotation_1",
  n_cores = 8
)

saveRDS(deconv, file = "./export2/paper/PAPER_4B_RCTD_3.rds")