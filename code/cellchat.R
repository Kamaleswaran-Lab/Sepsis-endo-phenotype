library(CellChat)
library(patchwork)
library(NMF)
library(ggalluvial)


load("./seurat_object1.Rda")
SRS1 <- subset(seurat_object1, subset = endotype =="SRS1")
data.SRS1 <- SRS1[["RNA"]]$data # normalized data matrix
cellchat_SRS1 <- createCellChat(object = data.SRS1, meta = SRS1[[]], group.by = "celltype")
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling
cellchat_SRS1@DB <- CellChatDB.use
cellchat_SRS1 <- subsetData(cellchat_SRS1)# This step is necessary even if using the whole database
cellchat_SRS1 <- identifyOverExpressedGenes(cellchat_SRS1)
cellchat_SRS1 <- identifyOverExpressedInteractions(cellchat_SRS1)
cellchat_SRS1 <- computeCommunProb(cellchat_SRS1, type = "triMean")
cellchat_SRS1 <- filterCommunication(cellchat_SRS1, min.cells = 10)
cellchat_SRS1 <- computeCommunProbPathway(cellchat_SRS1)
cellchat_SRS1 <- aggregateNet(cellchat_SRS1)
groupSize <- as.numeric(table(cellchat_SRS1@idents))
cellchat_SRS1 <- netAnalysis_computeCentrality(cellchat_SRS1, slot.name = "netP")
selectK(cellchat_SRS1, pattern = "outgoing")
nPatterns = 2
cellchat_SRS1 <- identifyCommunicationPatterns(cellchat_SRS1, pattern = "outgoing", k = nPatterns)
selectK(cellchat_SRS1, pattern = "incoming")
nPatterns = 2
cellchat_SRS1 <- identifyCommunicationPatterns(cellchat_SRS1, pattern = "incoming", k = nPatterns)
cellchat_SRS1 <- computeNetSimilarity(cellchat_SRS1, type = "functional")
cellchat_SRS1 <- netEmbedding(cellchat_SRS1, type = "functional")
cellchat_SRS1 <- netClustering(cellchat_SRS1, type = "functional")
cellchat_SRS1 <- computeNetSimilarity(cellchat_SRS1, type = "structural")
cellchat_SRS1 <- netEmbedding(cellchat_SRS1, type = "structural")
cellchat_SRS1 <- netClustering(cellchat_SRS1, type = "structural")

#same step for cellchat_SRS2


#compare SRS1 and SRS2

object.list <- list(SRS2 = cellchat_SRS2, SRS1 = cellchat_SRS1)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
compareInteractions(cellchat, show.legend = F, group = c(1,2))
compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_heatmap(cellchat)

weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
patchwork::wrap_plots(plots = gg)


cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)

rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)

i = 1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))


ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

gg1 <- netVisual_bubble(cellchat, sources.use = c(1,2), targets.use = c(1,2,9,10),  comparison = c(1, 2), max.dataset = 2, title.name = paste0("Increased signaling in ", names(object.list)[2]), angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, sources.use = c(1,2), targets.use = c(1,2,9,10),  comparison = c(1, 2), max.dataset = 1, title.name = paste0("Decreased signaling in ", names(object.list)[2]), angle.x = 45, remove.isolate = T)
gg1 + gg2

















