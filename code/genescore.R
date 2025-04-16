library(tidyverse)
library(data.table)
library(pals)
library(rawr)
library(pROC)
library(ggpubr)
library(scales)
library(ggbeeswarm)
library(ggrastr)
library(Seurat)
library(ComplexHeatmap)


geomMean <- function(x, na.rm = FALSE){
  if (!is.numeric(x) && !is.complex(x) && !is.logical(x)) {
    warning("argument is not numeric or logical: returning NA")
    return(as.numeric(NA))
  }
  if (na.rm){
    x <- x[!is.na(x)]
  }
  if (any(x < 0)){
    stop("'x' contains negative value(s)")
  }
  if(any(x == 0)){ #if I don't this, having any 0s will result in the geometric mean being 0
    x = x[-c(which(x == 0))]
    if(length(x) == 0){ #this means all values were 0
      return(0)
    }
  }
  return(exp(sum(log(x))/length(x)))
}

load("./integrated_seurat_object.Rda")
meta = alldata@meta.data 
umap_data <- read.csv(file = './UMAP.csv',row.names = 1)
umap_data$cell = rownames(umap_data)

genes.df<-alldata[["RNA"]]$data
minvalue<-min(genes.df)
genes.df<-genes.df-minvalue

endotype_genes <- read.csv("./genelist.csv")
genes_of_interest <- c(endotype_genes)
genes <- genes_of_interest[c(genes_of_interest) %in% rownames(genes.df)] 
genes.df = genes.df[genes,,drop = F] + 1
endotype_expr <- 2^genes.df[(row.names(genes.df) %in% endotype_genes),]


endotype_missing <- function(geneMtx, genes, out.missing=TRUE){
  if(out.missing){
    missinggenes = genes[!(genes %in% rownames(geneMtx))]
    missing = c(missinggenes)
    if(length(missing)>0){
      cat("Missing these genes:",missing,"\n")
    }
  }
}

endotype_missing(endotype_expr,genes = endotype_genes)

endotype_expr <- as.matrix(endotype_expr)
endotype_int1 <- data.frame(t(endotype_expr)) 
endotype_int2 <- endotype_int1%>%
  mutate(geomean = apply(endotype_int1, 1, FUN = geomMean, na.rm = TRUE))

endotype_table <- (endotype_int2 - endotype_int2[,"geomean"])^2 %>%
  mutate(endotype_score = rowSums(.)/1000000,
         cell = rownames(endotype_int2))%>%
  select("cell","endotype_score")

meta$cell = rownames(meta)
srt_meta <- meta

gene_counts <- t(as.matrix(genes.df))
gene_counts = data.frame(gene_counts)
gene_counts$cell = rownames(gene_counts)

srt_final <- left_join(srt_meta,endotype_table, by = "cell") %>%
  left_join(.,gene_counts, by = "cell") %>%
  left_join(.,umap_data, by = "cell")

srt_final <- srt_final[sample(1:nrow(srt_final)),]

srt_cell_prop <- srt_final%>%
  filter(!is.na(celltype))%>%
  summarize(`Neutrophil Immature` = log(mean(celltype == "Neutrophil Immature", na.rm = T)*1000),
            `Neutrophil Mature` = log(mean(celltype == "Neutrophil Mature", na.rm = T)*1000),
            `Neutrophil Progenitor` = log(mean(celltype == "Neutrophil Progenitor", na.rm = T)*1000),
            `Monocyte CD14` = log(mean(celltype == "Monocyte CD14", na.rm = T)*1000),
            `Monocyte CD16` = log(mean(celltype == "Monocyte CD16", na.rm = T)*1000),
            `B-Lymphocytes` = log(mean(celltype == "B-Lymphocytes", na.rm = T)*1000),
            `Dendritic cell` = log(mean(celltype == "Dendritic cell", na.rm = T)*1000),
            `T Lymphocyte CD4` = log(mean(celltype == "T Lymphocyte CD4", na.rm = T)*1000),
            `T-Lymphocyte CD8` = log(mean(celltype == "T-Lymphocyte CD8", na.rm = T)*1000),
            `Natural Killer Cell CD56` = log(mean(celltype == "Natural Killer Cell CD56", na.rm = T)*1000),
            `Macrophage M1` = log(mean(celltype == "Macrophage M1", na.rm = T)*1000),
            `Plasma Cells` = log(mean(celltype == "Plasma Cells", na.rm = T)*1000),
            `Mast Cell/Eosinophils` = log(mean(celltype == "Mast Cell/Eosinophils", na.rm = T)*1000))

endotype_genes_intersect <- intersect(colnames(srt_final),endotype_genes)

srt_endotype <- srt_final%>%
  dplyr::select(celltype,endotype_genes_intersect)%>%
  group_by(celltype) %>%
  summarize(across(everything(), mean, na.rm = T))

mat.expr = as.matrix(srt_endotype[,-1])
rownames(mat.expr) = srt_endotype$celltype
Endotype_up_score <- apply(mat.expr, 1, mean, na.rm = TRUE)
Endotype_up_score <- data.frame(Endotype_up_score)
Endotype_up_score = scale(Endotype_up_score)
Endotype_up_score = t(Endotype_up_score)
Endotype_up_score <- as.matrix(Endotype_up_score)

#same way for Endotype_down_score


