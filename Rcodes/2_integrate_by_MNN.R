library(dplyr)
library(Seurat)
library(patchwork)
library(SeuratWrappers)

##load seurat objects that were processed with DoubletFinder and add metadata
d0_1 <- readRDS(file = "d0_1.rds")
d0_2 <- readRDS(file = "d0_2.rds")
d0_3 <- readRDS(file = "d0_3.rds")
d0_1@meta.data$sample <- "d0_1"
d0_2@meta.data$sample <- "d0_2"
d0_3@meta.data$sample <- "d0_3"
d0 <- merge(d0_1, y = c(d0_2, d0_3), project = "Day0")
d0@meta.data$Day <- "Day0"

d7_1 <- readRDS(file = "d7_1.rds")
d7_2 <- readRDS(file = "d7_2.rds")
d7_3 <- readRDS(file = "d7_3.rds")
d7_1@meta.data$sample <- "d7_1"
d7_2@meta.data$sample <- "d7_2"
d7_3@meta.data$sample <- "d7_3"
d7 <- merge(d7_1, y = c(d7_2, d7_3), project = "Day7")
d7@meta.data$Day <- "Day7"

d14_1 <- readRDS(file = "d14_1.rds")
d14_2 <- readRDS(file = "d14_2.rds")
d14_3 <- readRDS(file = "d14_3.rds")
d14_1@meta.data$sample <- "d14_1"
d14_2@meta.data$sample <- "d14_2"
d14_3@meta.data$sample <- "d14_3"
d14 <- merge(d14_1, y = c(d14_2, d14_3), project = "Day14")
d14@meta.data$Day <- "Day14"

d21_1 <- readRDS(file = "d21_1.rds")
d21_2 <- readRDS(file = "d21_2.rds")
d21_3 <- readRDS(file = "d21_3.rds")
d21_1@meta.data$sample <- "d21_1"
d21_2@meta.data$sample <- "d21_2"
d21_3@meta.data$sample <- "d21_3"
d21 <- merge(d21_1, y = c(d21_2, d21_3), project = "Day21")
d21@meta.data$Day <- "Day21"

d0@meta.data$treatment <- "Untreated"
d7@meta.data$treatment <- "Bleomycin"
d14@meta.data$treatment <- "Bleomycin"
d21@meta.data$treatment <- "Bleomycin"

MG <- merge(d0, y = c(d7, d14, d21), project = "combine")


##integrate samples by FastMNN
MG <- NormalizeData(MG)
MG <- FindVariableFeatures(MG)
MG <- RunFastMNN(object.list = SplitObject(MG, split.by = "sample"))
MG <- RunUMAP(MG, reduction = "mnn", dims = 1:35)
MG <- FindNeighbors(MG, reduction = "mnn", dims = 1:35)
MG <- FindClusters(MG, resolution = 0.8)
DimPlot(MG, label = T)
DimPlot(MG, split.by = "Day", label = T)
saveRDS(MG, file = "MG.rds")

MG.markers <- FindAllMarkers(MG, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
MG.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5
all.genes <- rownames(MG)
MG <- ScaleData(MG, features = all.genes)
DoHeatmap(MG, features = top5$gene) + NoLegend()

table(MG@meta.data$seurat_clusters)
####
#0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18 
#6538 5279 4807 4729 4450 3968 3138 2536 2289 1929 1918 1604 1297 1256  705  694  339  168  165 
####
##cluster 17 are lineage+ cells. cluster 18 seems doublets

#exclude cluster 17 and 18
MG <- subset(MG, idents = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16))
MG <- DietSeurat(MG, counts = TRUE, data = TRUE, scale.data = FALSE)

#re-clustering
MG <- NormalizeData(MG)
MG <- FindVariableFeatures(MG)
MG <- RunFastMNN(object.list = SplitObject(MG, split.by = "sample"))
MG <- RunUMAP(MG, reduction = "mnn", dims = 1:38)
MG <- FindNeighbors(MG, reduction = "mnn", dims = 1:38)
MG <- FindClusters(MG, resolution = 0.3)
DimPlot(MG, label = T)
DimPlot(MG, split.by = "Day", label = T)

annotation <- c("Alveolar", "Fibrotic", "Alveolar", "Inflammatory", "Adventitial", "Smooth muscle",
                     "Peribronchial", "Stress-activated", "Proliferating", "Mesothelial", "Pericyte")
names(annotation) <- levels(MG)
MG <- RenameIdents(MG, annotation)
saveRDS(MG, file = "MGcleaned.rds")

##define tdTomato+ cells as tdTomato expression level>3.5
poscells <- WhichCells(MG, expression = rna_tdTomato > 3.5)
MG@meta.data$lineage<- ifelse(colnames(MG) %in% poscells, "tdTomato+", "Negative")
saveRDS(MG, file = "MGcleaned.rds")

##export metadata for tdTomato+ cell quantification
meta <- as.data.frame(Idents(MG))
meta$sample <- MG@meta.data$sample
meta$lineage <- MG@meta.data$lineage
write.table(meta, file="meta.txt")

##examine frequency of clusters at each time point
count_ID <- meta %>% count(sample, Idents.MG.)

##identify cluster markers
MG.markers <- FindAllMarkers(MG, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
MG <- ScaleData(MG, features = top10$gene)
DoHeatmap(MG, features = top10$gene) + NoLegend()

##run gsfisher
library(kableExtra)
library(knitr)
library(gsfisher)

MG.markers <- FindAllMarkers(MG, only.pos = TRUE, 
                                  min.pct = 0.25, logfc.threshold = 0.25)
saveRDS(MG.markers, "MGmarkers.rds")

getExpressedGenesFromSeuratObject <- function(seurat_object, 
                                              clusters, 
                                              min.pct=0.1)
{
  expressed <- c()
  for(cluster in clusters)
  {
    
    cluster_cells <- names(MG$seurat_clusters[MG$seurat_clusters==cluster])
    clust_pcts <- apply(MG@assays$RNA@counts[,cluster_cells], 
                        1, function(x) sum(x>0)/length(x))
    
    detected_in_clust <- names(clust_pcts[clust_pcts>min.pct])
    
    
    other_cells <- names(MG$seurat_clusters[MG$seurat_clusters!=cluster])
    other_pcts <- apply(MG@assays$RNA@counts[,other_cells], 
                        1, function(x) sum(x>0)/length(x))
    
    detected_in_other_cells <- names(other_pcts[other_pcts>min.pct])
    
    expressed <- c(expressed, detected_in_clust, detected_in_other_cells)
  }
  expressed <- unique(expressed)
  expressed
}

expressed_genes <- getExpressedGenesFromSeuratObject(
  MG,unique(MG$seurat_clusters), min.pct=0.25)

annotation <- fetchAnnotation(species="mm",
                              ensembl_version=NULL,
                              ensembl_host = NULL)

MG.markers$entrez_id <- as.character(annotation$entrez_id[
  match(MG.markers$gene, annotation$gene_name)])
MG.markers <- MG.markers[!is.na(MG.markers$entrez_id),]
background_entrez <- as.character(annotation$entrez_id[
  match(expressed_genes, annotation$gene_name)])
background_entrez <- background_entrez[!is.na(background_entrez)]

MG.markers.filtered <- MG.markers[MG.markers$p_val_adj < 0.05,]

go.results <- runGO.all(results=MG.markers.filtered, 
                        species = "mm",
                        background_ids = background_entrez,
                        gene_id_col="entrez_id",
                        gene_id_type="entrez",
                        sample_col="cluster",
                        p_col="p_val_adj",
                        p_threshold=0.05)
saveRDS(go.results, "go.results.MG.rds")

go.results.filtered <- go.results[go.results$ontology=="BP",]
go.results.filtered <- filterGenesets(go.results.filtered, 
                                      min_foreground_genes = 3,
                                      max_genes_geneset = 500, 
                                      min_odds_ratio = 2, 
                                      p_col = "p.val",
                                      padjust_method = "BH",
                                      use_adjusted_pvalues = TRUE,
                                      pvalue_threshold = 0.05)

go.results.top <- go.results.filtered %>% 
  group_by(cluster) %>% 
  top_n(n=5, -p.val)

sampleEnrichmentDotplot(go.results.top, 
                        selection_col = "description",
                        selected_genesets = unique(go.results.top$description),
                        sample_id_col = "cluster",
                        fill_var = "odds.ratio",
                        maxl=50,
                        title="GO biological pathway")
