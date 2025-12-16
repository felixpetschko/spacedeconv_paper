# ---- Setup ----
options(stringsAsFactors = FALSE)

# ---- Packages ----
suppressPackageStartupMessages({
  library(spacedeconv)
  library(SpatialExperiment)
})

# ---- Spatial Data ----
spe <- read10xVisium("./data/maynard/")
rownames(spe) <- rowData(spe)$symbol

# 151673
spatialLIBDAnnotation <- readRDS("./data/maynard/outs/spatial/spatial.rds")

# ---- Prepare Ground Truth ----
spatialLIBDAnnotation <- filter_sample_id(spatialLIBDAnnotation, "151673")
groundTruth <- spatialLIBDAnnotation$spatialLIBD
spe$groundTruth <- groundTruth

spe <- preprocess(spe)

rownames(spe) <- make.names(rownames(spe), unique = TRUE)

# ---- Reference Data ----
sce <- readRDS("./data/brainAtlas/allenBrain.rds")
sce <- subsetSCE(sce, cell_type_col = "subclass_label", ncells = 40000)
sce$subclass_label <- make.names(sce$subclass_label)

# ---- Deconvolution (cell2location) ----
# Note: 'gpu = TRUE' assumes you have a compatible GPU environment configured.
signature <- build_model(
  sce,
  cell_type_col = "subclass_label",
  method = "cell2location",
  epochs = 250,
  batch_id_col = "external_donor_name_label",
  gpu = TRUE,
  assay_sc = "counts"
)

deconvC2L <- deconvolute(
  spe,
  method = "cell2location",
  signature = signature,
  epochs = 30000,
  gpu = TRUE,
  assay_sc = "counts",
  assay_sp = "counts"
)

# ---- Save Output ----
saveRDS(deconvC2L, file = "./export2/paper/PAPER_4A_C2L_3.rds")
cat("Saved deconvolution result to ./export2/paper/PAPER_4A_C2L_3.rds\n")