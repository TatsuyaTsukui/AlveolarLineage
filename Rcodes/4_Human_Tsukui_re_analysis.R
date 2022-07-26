library(dplyr)
library(Seurat)
library(patchwork)
library(SeuratWrappers)
library(monocle3)
library(kableExtra)
library(knitr)
library(gsfisher)
library(tidyr)
library(scCustomize)
library(biomaRt)


#load seurat object of human mesenchymal cells from Tsukui et al. 2020
#Email me tatsuya.tsukui@ucsf.edu if you need this seurat object
hlin <- readRDS(file="hLung0827.rds")

##subset alveolar and pathologic fibroblasts
alpfb <- subset(hlin, idents = c(1, 3))

##re-clustering 
alpfb <- NormalizeData(alpfb)
alpfb <- FindVariableFeatures(alpfb)
alpfb <- RunFastMNN(object.list = SplitObject(alpfb, split.by = "ID"))
alpfb <- RunUMAP(alpfb, reduction = "mnn", dims = 1:20)
alpfb <- FindNeighbors(alpfb, reduction = "mnn", dims = 1:20)
alpfb <- FindClusters(alpfb, resolution = 0.3)
DimPlot(alpfb, label = T)
DimPlot(alpfb, group.by = "ID", label = T)

alpfb.markers <- FindAllMarkers(alpfb, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

annotation <- c("Alveolar", "Alveolar", "Inflammatory1", "Fibrotic", "Inflammatory2")
names(annotation) <- levels(alpfb)
alpfb <- RenameIdents(alpfb, annotation)
saveRDS(file = "alpfb.rds")

##run Monocle 3
alpfb.cds <- as.cell_data_set(alpfb)
alpfb.cds <- estimate_size_factors(alpfb.cds)
alpfb.cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(alpfb[["RNA"]])
alpfb.cds <- cluster_cells(cds = alpfb.cds, reduction_method = "UMAP")
alpfb.cds <- learn_graph(alpfb.cds, use_partition = TRUE)

get_earliest_principal_node <- function(alpfb.cds, time_bin="Normal"){
  cell_ids <- which(colData(alpfb.cds)[, "disease"] == time_bin)
  
  closest_vertex <-
    alpfb.cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(alpfb.cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(alpfb.cds)[["UMAP"]])$name[as.numeric(names
                                                                    (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
alpfb.cds <- order_cells(alpfb.cds, root_pr_nodes=get_earliest_principal_node(alpfb.cds))

plot_cells(alpfb.cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

saveRDS(alpfb, file = "alpfb.rds")
saveRDS(alpfb.cds, file = "alpfb.cds.rds")

##run gsfisher
getExpressedGenesFromSeuratObject <- function(seurat_object, 
                                              clusters, 
                                              min.pct=0.1)
{
  expressed <- c()
  for(cluster in clusters)
  {
    
    cluster_cells <- names(alpfb$seurat_clusters[alpfb$seurat_clusters==cluster])
    clust_pcts <- apply(alpfb@assays$RNA@counts[,cluster_cells], 
                        1, function(x) sum(x>0)/length(x))
    
    detected_in_clust <- names(clust_pcts[clust_pcts>min.pct])
    
    
    other_cells <- names(alpfb$seurat_clusters[alpfb$seurat_clusters!=cluster])
    other_pcts <- apply(alpfb@assays$RNA@counts[,other_cells], 
                        1, function(x) sum(x>0)/length(x))
    
    detected_in_other_cells <- names(other_pcts[other_pcts>min.pct])
    
    expressed <- c(expressed, detected_in_clust, detected_in_other_cells)
  }
  expressed <- unique(expressed)
  expressed
}

expressed_genes <- getExpressedGenesFromSeuratObject(
  alpfb,unique(alpfb$seurat_clusters), min.pct=0.25)

annotation <- fetchAnnotation(species="hs",
                              ensembl_version=NULL,
                              ensembl_host = NULL)


alpfb.markers$entrez_id <- as.character(annotation$entrez_id[
  match(alpfb.markers$gene, annotation$gene_name)])
alpfb.markers <- alpfb.markers[!is.na(alpfb.markers$entrez_id),]
background_entrez <- as.character(annotation$entrez_id[
  match(expressed_genes, annotation$gene_name)])
background_entrez <- background_entrez[!is.na(background_entrez)]

alpfb.markers.filtered <- alpfb.markers[alpfb.markers$p_val_adj < 0.05,]

go.results <- runGO.all(results=alpfb.markers.filtered, 
                          species = "hs",
                          background_ids = background_entrez,
                          gene_id_col="entrez_id",
                          gene_id_type="entrez",
                          sample_col="cluster",
                          p_col="p_val_adj",
                          p_threshold=0.05)

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
                        
                        
##Correlation alaysis for human and mouse emergent clusters
library(dplyr)
library(Seurat)
library(patchwork)
library(SeuratWrappers)
library(ggplot2)
library(pheatmap)
                        
MG <- readRDS(file = "MGcleaned.rds") #mouse seurat object
alpfb <- readRDS(file = "alpfb.rds")
MGpfb <- subset(MG, idents = c("Alveolar", "Inflammatory", "Stress-activated", "Fibrotic"))
                        
alpfb <- ScaleData(alpfb)
MGpfb <- ScaleData(MGpfb)
AEm <- AverageExpression(MGpfb, assays = "RNA", slot = "scale.data")
AEh <- AverageExpression(alpfb, assays = "RNA", slot = "scale.data")
AEM <- as.data.frame(AEm$RNA)
AEH <- as.data.frame(AEh$RNA)
hgenes <- rownames(AEH)

human = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host="https://dec2021.archive.ensembl.org/")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl" , host="https://dec2021.archive.ensembl.org/")

test <- getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", 
               values = hgenes , mart = human, attributesL = c("mgi_symbol"), 
               martL = mouse, uniqueRows=T)

id <- match(hgenes , test$HGNC.symbol)
AEH$mgene <- test$MGI.symbol[id]
non <- which(duplicated(AEH$mgene))
AEH <- AEH[-non,]
which(is.na(AEH$mgene))
AEH <- AEH[-4,]
rownames(AEH) <- AEH$mgene
AEH[,5] <- NULL
                        
hAllgene <- rownames(AEH)
mousegene <- rownames(AEM)
genes.use <- intersect(hAllgene, mousegene)

AEH <- AEH[genes.use, , drop=FALSE]
AEM <- AEM[genes.use,]
AEH <- AEH[, c(1,2,4,3)]
AEM <- AEM[, c(1,3,4,2)]

comball <- cbind(AEH, AEM)
hcorrelation <- cor(comball, method = "spearman")
trim <- hcorrelation[1:4, 5:8]
pheatmap(trim, cluster_cols = F, cluster_rows = F)
write.table(hcorrelation, file="hcorrelation.csv")
