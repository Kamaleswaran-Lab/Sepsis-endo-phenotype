library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)


load("./seurat_object1.Rda")
load("./seurat_object2.Rda")
load("./seurat_object3.Rda")
load("./seurat_object4.Rda")
load("./seurat_object5.Rda")
load("./seurat_object6.Rda")


merged_obj <- merge(x = seurat_object1, y = c(seurat_object2, seurat_object3, seurat_object4, seurat_object5,seurat_object6))
merged_obj <- NormalizeData(merged_obj)
merged_obj <- FindVariableFeatures(merged_obj)
merged_obj <- ScaleData(merged_obj)
merged_obj <- RunPCA(merged_obj)
merged_obj <- IntegrateLayers(object = merged_obj, method = RPCAIntegration, orig.reduction = "pca", new.reduction = "integrated.rpca",
                              verbose = FALSE)
merged_obj[["RNA"]] <- JoinLayers(merged_obj)

merged_obj <- FindNeighbors(merged_obj, reduction = "integrated.rpca", dims = 1:30)
merged_obj <- FindClusters(merged_obj, resolution = 1)
merged_obj <- RunUMAP(merged_obj, dims = 1:30, reduction = "integrated.rpca")

DimPlot(merged_obj, reduction = "integrated.rpca")

marker = c("CD79A","CD79B","XBP1",#b plasma
            "FCGR3A","CD14","S100A9",#CD14 CD16 Mono
            "CD3D","CD3G", "CCR7","LEF1","TCF7","IL7R","LTB",#cd4
            "CD8A","CCL5","GZMH","TRGC2",#cd8
            "FCER1A","CST3","JCHAIN","MZB1",#DC
            "HBB",#Erythro
            "HLA-DQA1",#M1 Macro
            "S100A8","MMP8","MMP9","NAMPT","NEAT1","LCN2","PADI4","ITGAM",#Neutrophil
            "MPO","ELANE","AZU1",#Neutrophil_progenitors
            "GNLY","NKG7","PRF1","FGFBP2","GZMA",#nk
            "ITGA2B","ITGB3","GP1BA","SELP"#Platelets
)


#cell annotation
DotPlot(object = merged_obj, features = marker) + RotatedAxis()

save(merged_obj, file ="/integrated_seurat_object.Rda")


