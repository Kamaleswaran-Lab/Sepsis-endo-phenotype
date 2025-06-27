library(Seurat)
library(DESeq2)
library(dplyr)
library(ggplot2)

options(future.globals.maxSize = 100e+09)

new_counts <- Read10X(data.dir = "./whole_blood/",gene.column=1)
seurat_object = CreateSeuratObject(counts = new_counts, min.cells = 0, min.features = 0, project = 'kwok')
seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")

new_counts_adt <- Read10X(data.dir = "./whole_blood/adt/",gene.column=1)
seurat_object1 = CreateAssayObject(counts = new_counts_adt, min.cells = 0, min.features = 0, project = 'kwok_adt')
seurat_object[["ADT"]] <- seurat_object1

METADATA <- read.csv("./whole_blood/Kwok-whole-blood_cell-metadata.tsv", row.names = 1, sep = "\t")

seurat_object <- AddMetaData(seurat_object, metadata = METADATA["sample_id"])


alldata <- load("./all_integrated_harmony.Rda")

cite_cells <- colnames(alldata)[alldata$dataset == "Kwok_whole_blood"]
adt_counts <- GetAssayData(seurat_object, assay = "ADT", slot = "counts")[, cite_cells]

adt_assay <- CreateAssayObject(counts = adt_counts)
alldata[["ADT"]] <- adt_assay

alldata <- NormalizeData(alldata, assay = "ADT", normalization.method = "CLR")

save(alldata, file ="./kwok_ADT_harmony.Rda")

DefaultAssay(alldata) <- "ADT"
adt_features <- rownames(alldata[["ADT"]])

data_df <- FetchData(alldata, vars = c(adt_features, "celltype_sub"), layer = "data", assay = "ADT")
colnames(data_df)

plot_list <- list()
for (feature in adt_features) {
  p <- ggplot(data_df, aes(x = celltype_sub, y = !!sym(feature))) +
    geom_boxplot(outlier.shape = NA) +
    theme_minimal(base_size = 10) +
    labs(title = feature, x = "", y = "CLR expression") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  plot_list[[feature]] <- p
}


pdf("ADT_boxplots_celltypesub_harmony.pdf", width = 12, height = 10)
n <- length(plot_list)
for (i in seq(1, n, by = 6)) {
  print(wrap_plots(plot_list[i:min(i+5, n)], ncol = 2))
}
dev.off()





