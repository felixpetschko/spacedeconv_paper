library(spacedeconv)
library(SpatialExperiment)
library(corrplot)
library(gridExtra)
library(RColorBrewer)
library(ggplot2)

c2l <- readRDS("./export/4B_C2L.rds")
rctd <- readRDS("./export/4B_RCTD.rds")

spot_to_remove <- c(
  "AATCGAGGTCTCAAGG-1", "ACCACACGGTTGATGG-1",
  "ACCGAGTCTCCTTATT-1", "AGCAGTCGAAGCATGC-1",
  "AGTACAGAAGCTTATA-1", "CCTACAGTTGAGGGAG-1",
  "CGCTTTCCGCCAAGGT-1", "CGGTATGGGCACTCTG-1",
  "CGGTTGGGCAGGGTCC-1", "CTATCACAACGCTGGA-1",
  "CTCTAACACCGGCAGC-1", "CTTACATAGATTTCTT-1",
  "CTTCATCACCAGGGCT-1", "GAAGGGTCATTAAGAC-1",
  "GACCACACTTCCCTTT-1", "GCACGCCGATTCCCGC-1",
  "GCGGCTTTAGCAAGTT-1", "GCGTCGTAACATGGTC-1",
  "GGATTAATCATGGACC-1", "GTACTACGGCCTCGTT-1",
  "TACACAGCCGTGGTGC-1", "TCAGCGCACGCCGTTT-1",
  "TGCCCGATAGTTAGAA-1", "TGTAGCCAATTCCGTT-1"
)

duplicate <- c(
    "Gm16701",
    "Gm2464",
    "Aldoa",
    "Gm35558",
    "Gcat",
    "Pick1",
    "Atp5o",
    "Ints5"
)

c2l <- c2l[, !colnames(c2l) %in% spot_to_remove]
c2l <- c2l[!rownames(c2l) %in% duplicate, ]

cluster <- spacedeconv::cluster(c2l, data = "expression", clusres = c(0.1, 0.6))

get_cluster_features(
    cluster,
    clusterid = "cluster_expression_res_0.1",
    spmethod = "cell2location",
    topn = 5
)

get_cluster_features(
    cluster,
    clusterid = "cluster_expression_res_0.6",
    spmethod = "cell2location",
    topn = 5
)

c2l2 <- c2l

# aggregate the Inh types
c2l2 <- aggregate_results(
    c2l2,
    cell_types = c(
        "cell2location_Inh_1",
        "cell2location_Inh_2",
        "cell2location_Inh_3",
        "cell2location_Inh_5",
        "cell2location_Inh_4",
        "cell2location_Inh_6",
        "cell2location_Inh_Lamp5",
        "cell2location_Inh_Pvalb",
        "cell2location_Inh_Sst",
        "cell2location_Inh_Meis2_1",
        "cell2location_Inh_Meis2_2",
        "cell2location_Inh_Meis2_3",
        "cell2location_Inh_Meis2_4",
        "cell2location_Inh_Vip"
    ),
    name = "cell2location_Inh",
    remove = TRUE
)

# aggregate the Astro Types
c2l2 <- aggregate_results(
    c2l2,
    cell_types = c(
        "cell2location_Astro_CTX",
        "cell2location_Astro_HPC",
        "cell2location_Astro_HYPO",
        "cell2location_Astro_THAL_lat",
        "cell2location_Astro_STR",
        "cell2location_Astro_WM",
        "cell2location_Astro_THAL_hab",
        "cell2location_Astro_THAL_med",
        "cell2location_Astro_AMY",
        "cell2location_Astro_AMY_CTX"
    ),
    name = "cell2location_Astro",
    remove = TRUE
)

# aggregate HPC types (Ext group)
c2l2 <- aggregate_results(
    c2l2,
    cell_types = c(
        "cell2location_Ext_Hpc_CA1",
        "cell2location_Ext_Hpc_CA2",
        "cell2location_Ext_Hpc_CA3",
        "cell2location_Ext_Hpc_DG1",
        "cell2location_Ext_Hpc_DG2",
        "cell2location_Ext_L5_1",
        "cell2location_Ext_L5_2",
        "cell2location_Ext_L56",
        "cell2location_Ext_L5_3",
        "cell2location_Ext_L6",
        "cell2location_Ext_L6B",
        "cell2location_Ext_Thal_1",
        "cell2location_Ext_Thal_2",
        "cell2location_Ext_L25",
        "cell2location_Ext_L23",
        "cell2location_Ext_Amy_1",
        "cell2location_Ext_Amy_2",
        "cell2location_Ext_ClauPyr",
        "cell2location_Ext_Med",
        "cell2location_Ext_Pir"
    ),
    name = "cell2location_Ext",
    remove = TRUE
)

# aggregate Oligos
c2l2 <- aggregate_results(
    c2l2,
    cell_types = c("cell2location_Oligo_1", "cell2location_Oligo_2"),
    name = "cell2location_Oligo",
    remove = TRUE
)

# aggregate OPC
c2l2 <- aggregate_results(
    c2l2,
    cell_types = c("cell2location_OPC_1", "cell2location_OPC_2"),
    name = "cell2location_OPC",
    remove = TRUE
)

# aggregate NB
c2l2 <- aggregate_results(
    c2l2,
    cell_types = c("cell2location_Nb_1", "cell2location_Nb_2"),
    name = "cell2location_Nb",
    remove = TRUE
)

title_size <- 22
font_size <- 18
legend_size <- 20

mAbundantAgg <- plot_most_abundant(
    c2l2,
    method = "cell2location",
    title = "Most Abundant",
    min_abundance = 0.05,
    palette = "Dark2",
    title_size = title_size,
    font_size = font_size,
    legend_size = legend_size,
    spot_size = 1.03
) +
    scale_fill_manual(
        values = RColorBrewer::brewer.pal(6, "Dark2"),
        labels = c("Astro", "Endo", "Ext", "Inh", "Nb", "Oligo")
    )

clusterC2L <- spacedeconv::cluster(
    spe = c2l,
    data = "deconvolution",
    spmethod = "cell2location",
    nclusters = c(6, 8, 10, 12, 14, 16, 18, 20, 22),
    method = "hclust"
)

clusterC2L <- spacedeconv::cluster(
    spe = c2l,
    spmethod = "cell2location",
    nclusters = c(6, 8, 10, 12, 14, 16, 18, 20, 22),
    method = "kmeans"
)

hypo <- plot_spatial(
    c2l,
    result = "cell2location_Astro_HYPO",
    density = FALSE,
    smooth = TRUE,
    title = "Astro HYPO C2L",
    title_size = title_size,
    font_size = font_size,
    legend_size = legend_size,
    spot_size = 1.03,
    shift_positive = FALSE
)

thal <- plot_spatial(
    c2l,
    result = "cell2location_Astro_THAL_lat",
    density = FALSE,
    smooth = TRUE,
    title = "Astro THAL C2L",
    title_size = title_size,
    font_size = font_size,
    legend_size = legend_size,
    spot_size = 1.03,
    shift_positive = FALSE
)

wm <- plot_spatial(
    c2l,
    result = "cell2location_Astro_WM",
    density = FALSE,
    smooth = TRUE,
    title = "Astro WM C2L (sqrt)",
    transform_scale = "sqrt",
    title_size = title_size,
    font_size = font_size,
    legend_size = legend_size,
    spot_size = 1.03,
    shift_positive = FALSE
)

clus01 <- plot_spatial(
    cluster,
    result = "cluster_expression_res_0.1",
    density = FALSE,
    smooth = FALSE,
    title = "Cluster 0.1",
    palette = "Accent",
    title_size = title_size,
    font_size = font_size,
    legend_size = legend_size,
    spot_size = 1.03,
    shift_positive = FALSE
)

clus06 <- plot_spatial(
    cluster,
    result = "cluster_expression_res_0.6",
    density = FALSE,
    smooth = FALSE,
    title = "Cluster 0.6",
    palette = "Accent",
    title_size = title_size,
    font_size = font_size,
    legend_size = legend_size,
    spot_size = 1.03,
    shift_positive = FALSE
)

c2l3 <- aggregate_results(
    c2l,
    cell_types = c("cell2location_Ext_L5_1", "cell2location_Ext_L5_2"),
    remove = TRUE,
    name = "cell2location_L5"
)

l5l6 <- plot_comparison(
    c2l3,
    cell_type_1 = "cell2location_L5",
    cell_type_2 = "cell2location_Ext_L6",
    density = FALSE,
    smooth = TRUE,
    title = "L5 vs L6",
    palette = "Purple-Green",
    title_size = title_size,
    font_size = font_size,
    legend_size = legend_size,
    spot_size = 1.03,
    shift_positive = FALSE
)

an1 <- read.csv("./data/20200904_regionAnnotation_per_location.csv")
rownames(an1) <- an1$spot_id
an1 <- an1[grepl("ST8059048", an1$spot_id), ]
rownames(an1) <- gsub(".*_", "", rownames(an1))
an1 <- an1[colnames(c2l), ]
c2l$brain_region <- an1$location

truth <- plot_spatial(
    c2l,
    result = "brain_region",
    density = FALSE,
    palette = "Accent",
    title = "Ground Truth (mouse)",
    title_size = title_size,
    font_size = font_size,
    legend_size = legend_size,
    spot_size = 1.03,
    shift_positive = FALSE
)

grid_layout <- rbind(
    c(1, 1, 1, 2, 2, 2),
    c(3, 3, 4, 4, 5, 5),
    c(6, 6, 7, 7, 8, 8)
)

final <- grid.arrange(
    truth, mAbundantAgg,
    hypo, thal, wm,
    l5l6, clus01, clus06,
    layout_matrix = grid_layout
)

ggsave(
    filename = "./export/fig_4B.png",
    plot = final,
    dpi = 600,
    width = 15,
    height = 15,
    units = "in",
    bg = "white"
)
