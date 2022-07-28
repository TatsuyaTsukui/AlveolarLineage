library(monocle3)
library(dplyr)
library(Seurat)
library(patchwork)
library(SeuratWrappers)
library(ggplot2)

##load seurat object
MG <- readRDS(file = "MGcleaned.rds")

##subset alveolar, inflammatory, stress-activated, fibrotic, and proliferating to run Monocle3
Trc <- subset(MG, idents = c("Alveolar1", "Alveolar2", "Inflammatory", "Stress-activated", "Fibrotic", "Proliferating"))
Trc <- NormalizeData(Trc)
Trc <- FindVariableFeatures(Trc)
Trc <- RunFastMNN(object.list = SplitObject(Trc, split.by = "sample"))
Trc <- RunUMAP(Trc, reduction = "mnn", dims = 1:38)
Trc <- FindNeighbors(Trc, reduction = "mnn", dims = 1:38)
Trc <- FindClusters(Trc, resolution = 0.3)
DimPlot(Trc, label = T)
annotation <- c("Alveolar1", "Fibrotic", "Alveolar2", "Inflammatory", "Stress-ativated", "Proliferating")
names(annotation) <- levels(Trc)
Trc <- RenameIdents(Trc, annotation)

#export seurat object for Monocle3
Trc.cds <- as.cell_data_set(Trc)
Trc.cds <- cluster_cells(cds = Trc.cds, reduction_method = "UMAP")
Trc.cds <- learn_graph(Trc.cds, use_partition = TRUE)

get_earliest_principal_node <- function(Trc.cds, time_bin="Day0"){
  cell_ids <- which(colData(Trc.cds)[, "Day"] == time_bin)
  
  closest_vertex <-
    Trc.cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(Trc.cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(Trc.cds)[["UMAP"]])$name[as.numeric(names
                                                                  (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
Trc.cds <- order_cells(Trc.cds, root_pr_nodes=get_earliest_principal_node(Trc.cds))
plot_cells(Trc.cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

saveRDS(Trc, file="Trc.rds")
saveRDS(Trc.cds, file="Trc_cds.rds")


##subset alveolar, inflammatory, and fibrotic clusters to run Monocle3
ltr <- subset(MG, idents = c("Alveolar1", "Alveolar2", "Inflammatory", "Fibrotic"))
ltr <- FindVariableFeatures(ltr)
ltr <- RunFastMNN(object.list = SplitObject(ltr, split.by = "sample"))
ltr <- RunUMAP(ltr, reduction = "mnn", dims = 1:38)
ltr <- FindNeighbors(ltr, reduction = "mnn", dims = 1:38)
ltr <- FindClusters(ltr, resolution = 0.25)
DimPlot(ltr, label = T)

FeaturePlot(ltr, features = c("Ces1d", "Saa3", "Cthrc1", "Spp1"))

annotation <- c("Alveolar", "Fibrotic", "Alveolar", "Inflammatory")
names(annotation) <- levels(ltr)
ltr <- RenameIdents(ltr, annotation)

ltr.cds <- as.cell_data_set(ltr)
ltr.cds <- estimate_size_factors(ltr.cds)
ltr.cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(ltr[["RNA"]])
ltr.cds <- cluster_cells(cds = ltr.cds, reduction_method = "UMAP")
ltr.cds <- learn_graph(ltr.cds, use_partition = TRUE)

get_earliest_principal_node <- function(ltr.cds, time_bin="Day0"){
  cell_ids <- which(colData(ltr.cds)[, "Day"] == time_bin)
  
  closest_vertex <-
    ltr.cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(ltr.cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(ltr.cds)[["UMAP"]])$name[as.numeric(names
                                                                  (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
ltr.cds <- order_cells(ltr.cds, root_pr_nodes=get_earliest_principal_node(ltr.cds))
plot_cells(ltr.cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

saveRDS(ltr, file = "ltr.rds")
saveRDS(ltr.cds, file = "ltr.cds.rds")

##generate a graph of gene expression change along pseudotime
pt <- pseudotime(ltr.cds)
ltr@meta.data$pseudotime <- pt
ltr@meta.data$CellType <- Idents(ltr)
ltr <- ScaleData(ltr)
data.use <- data.frame(FetchData(ltr, vars = c("pseudotime", "CellType")), 
                       t(ltr@assays$RNA@scale.data[c("Pdgfra", "Tcf21", "Saa3", "Lcn2", "Spp1", "Cthrc1"), ]))
data.use$Cell <- rownames(data.use)
data.usemod <- data.use %>% tidyr::gather(key="genes", value="value", -pseudotime, -CellType, -Cell)
color1 <- "green4" 
color2 <- "red"
color3 <- "blue" 
color4 <- "lightgrey" 

p1 <- ggplot(data.use) +
  theme_bw() +
  geom_smooth(aes(pseudotime, Pdgfra), span = 0.75, method = "loess", color = color2) +
  geom_smooth(aes(pseudotime, Tcf21), span = 0.75, method = "loess", color = color2) +
  geom_smooth(aes(pseudotime, Saa3), span = 0.75, method = "loess", color = color3) +
  geom_smooth(aes(pseudotime, Lcn2), span = 0.75, method = "loess", color = color3) +
  geom_smooth(aes(pseudotime, Spp1), span = 0.75, method = "loess", color = color1) +
  geom_smooth(aes(pseudotime, Cthrc1), span = 0.75, method = "loess", color = color1) +
  scale_y_continuous(limits = c(-2, NA)) + xlab("Pseudotime") + ylab("expression") +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


##generate a heatmap of gene expression change along pseudotime by SlingShot
library(monocle3)
library(dplyr)
library(Seurat)
library(patchwork)
library(SeuratWrappers)
library(slingshot)
library(SingleCellExperiment)
library(scater)
library(RColorBrewer)
library(grDevices)
library(tradeSeq)

ltr <- readRDS(file = "ltr.rds")

ltr <- FindVariableFeatures(ltr, nfeatures = 3000)
HVG <- ltr@assays$RNA@var.features
ltrs <- subset(x = ltr, features = HVG)
ltr.sce <- as.SingleCellExperiment(ltrs)
ltr.sce <- slingshot(ltr.sce, clusterLabels = 'ident', reducedDim = 'UMAP', start.clus = "Alveolar", end.clus = "Fibrotic")
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(ltr.sce$slingPseudotime_1, breaks=100)]
plot(reducedDims(ltr.sce)$UMAP, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(ltr.sce), lwd=2, col='black')

ltr.sce <- fitGAM(ltr.sce)
saveRDS(ltr.sce, file = "ltr.sce.rds")
ATres <- associationTest(ltr.sce)

library(ComplexHeatmap)
genes <- c("Pdgfra","Tcf21","Npnt","Inmt","Ces1d","Saa3","Lcn2","Sod2","Hp","Cxcl12","Sfrp1","Spp1", "Fst","S100a4", "Cd9","Postn","Cthrc1","Grem1")
pst.ord <- order(ltr.sce$slingPseudotime_1, na.last = NA)
heatdata <- GetAssayData(ltr, slot = "scale.data")[genes, pst.ord]
heatclus <- ltr.sce$ident[pst.ord]
ha = HeatmapAnnotation( Cluster = heatclus,  col= list(Cluster = c("Alveolar" = "firebrick1", "Inflammatory" = "dodgerblue3", "Fibrotic" = "green3")))
Heatmap(heatdata, name = "Expression",  column_dend_reorder = FALSE, cluster_rows = FALSE, row_dend_reorder = FALSE, cluster_columns = FALSE, show_column_names = FALSE, top_annotation = ha)
