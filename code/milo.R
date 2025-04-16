library(miloR)
library(SingleCellExperiment)
library(dplyr)
library(patchwork)
library(Seurat)
library(ggplot2)


load("./seurat_object1.Rda")
embryo_data <- as.SingleCellExperiment(seurat_object1)
embryo_milo <- Milo(embryo_data)

embryo_milo <- buildGraph(embryo_milo, k = 30, d = 30, reduced.dim = "INTEGRATED.RPCA")
embryo_milo <- makeNhoods(embryo_milo, prop = 0.1, k = 30, d=30, refined = TRUE, reduced_dims = "INTEGRATED.RPCA")
embryo_milo <- countCells(embryo_milo, meta.data = data.frame(colData(embryo_milo)), sample="patients")
embryo_design <- data.frame(colData(embryo_milo))[,c("patients", "endotype")]
desired_order <- c("SRS1","SRS2")
embryo_design$endotype <- as.factor(embryo_design$endotype)
levels(embryo_design$endotype) <- desired_order
embryo_design <- distinct(embryo_design)
rownames(embryo_design) <- embryo_design$patients
embryo_milo <- calcNhoodDistance(embryo_milo, d=30, reduced.dim = "INTEGRATED.RPCA")
da_results <- testNhoods(embryo_milo, design = ~ 0 + endotype, design.df = embryo_design)
embryo_milo <- buildNhoodGraph(embryo_milo)
plotNhoodGraphDA(embryo_milo, da_results, layout="UMAP",alpha=0.05)
da_results1 <- annotateNhoods(embryo_milo, da_results, coldata_col = "celltype")
plotDAbeeswarm(da_results1, group.by = "celltype")

