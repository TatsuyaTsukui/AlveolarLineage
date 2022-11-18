library(dplyr)
library(kableExtra)
library(knitr)
library(Seurat)
library(patchwork)
library(gsfisher)
library(tidyr)
library(SeuratWrappers)
library(scCustomize)

##load alveolar inflammatory fbirotic cluster seurat object re-analyzed from Tsukui et al. 2020. add metadata
alpfb <- readRDS(file="alpfb.rds")
alpfb@meta.data$DataSet <- "Tsukui"
alpfb@meta.data$Tsukui_clusters <- Idents(alpfb)
alpfb@meta.data$disease <- gsub("Normal", "Control", alpfb@meta.data$disease)

##load Adams et al. mesenchyme data set

RC <- read.table(file = "GSE147066_RawCounts_triples.txt")
names(RC) <- RC[1,]
RC <- RC[-1,]
adams <- acast(RC, row~col, value.var='UMI_Counts', fill=0)

MD <- read.table(file = "GSE147066_Cell_MetaData.txt")
names(MD) <- MD[1,]
MD <- MD[-1,]
rownames(MD)<- MD[,1]
MD <- MD[,-1]

##create seurat object and quality check
Adams <- CreateSeuratObject(adams, project = "SeuratProject", assay = "RNA",
                            min.cells = 3, min.features = 200, names.field = 1,
                            names.delim = "_", meta.data = MD)
Adams[["percent.mt"]] <- PercentageFeatureSet(Adams, pattern = "^MT-")
VlnPlot(Adams, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
saveRDS(Adams, file = "adams.mesenchyme.rds")

##clustering by FastMNN
Adams <- NormalizeData(Adams)
Adams <- FindVariableFeatures(Adams)
Adams <- RunFastMNN(object.list = SplitObject(Adams, split.by = "subject.ident"))
Adams <- RunUMAP(Adams, reduction = "mnn", dims = 1:25)
Adams <- FindNeighbors(Adams, reduction = "mnn", dims = 1:25)
Adams <- FindClusters(Adams, resolution = 0.8)
DimPlot(Adams, label = T)
FeaturePlot(Adams, features = c("CTHRC1", "COL1A1", "PI16", "SFRP2"))
DimPlot(Adams, split.by = "disease.ident", label = T)
saveRDS(Adams, file = "adams.mesenchyme.rds")

##subset alveolar and pathologic fibroblast. add meta data
Adams.alvp <- subset(x = Adams, idents = c(0,1,3,8))
Adams.alvp@meta.data$disease <- Adams.alvp@meta.data$disease.ident
Adams.alvp@meta.data$DataSet <- "Adams"
Adams.alvp@meta.data$ID <- Adams.alvp@meta.data$subject.ident
Adams.alvp@meta.data$Tsukui_clusters <- "Adams"
saveRDS(Adams.alvp, file = "Adams.alvp.rds")

##load mesenchymal cell seurat object from Habermann et al.
NVFB <- readRDS(file="NvFB.rds") ##this object was made by subsetting mesenchymal cells from the full size seurat object from GSE135893
Idents(NVFB) <- "celltype"

##subset "Myofibroblasts" cluster that contains alveolar and pathologic fibroblasts. add metadata
myo <- subset(NVFB, idents = "Myofibroblasts")
myo@meta.data$disease <- myo@meta.data$Diagnosis
myo@meta.data$ID <- myo@meta.data$Sample_Name
myo@meta.data$DataSet <- "Habermann"
myo@meta.data$Tsukui_clusters <- "Habermann"

saveRDS(myo, file = "NVmyosubset.rds")

##integrate the three data sets by FastMNN with split.by = ID (individual patients or donors)
HFBmg <- merge(alpfb, y = c(myo, Adams.alvp), project = "HumanFBmg")
HFBmg <- NormalizeData(HFBmg)
HFBmg <- FindVariableFeatures(HFBmg)
HFBmg <- RunFastMNN(object.list = SplitObject(HFBmg, split.by = "ID"))
HFBmg <- RunUMAP(HFBmg, reduction = "mnn", dims = 1:21)
HFBmg <- FindNeighbors(HFBmg, reduction = "mnn", dims = 1:21)
HFBmg <- FindClusters(HFBmg, resolution = 0.8)
DimPlot(HFBmg, label = T)
DimPlot(HFBmg, split.by = "DataSet", label = T)
FeaturePlot(HFBmg, features = c("CTHRC1", "NPNT", "SFRP2", "SFRP4"))
FeaturePlot(HFBmg, features = "nFeature_RNA")
HFBmg[["percent.mt"]] <- PercentageFeatureSet(HFBmg, pattern = "^MT-")
saveRDS(HFBmg, file = "HFBmg.rds")


##exclude cluster 9 and 10, which are unique to Adams et al. and seem driven by technical reasons
HFBmgss <- subset(x = HFBmg, idents = c(9, 10), invert = TRUE)

##re-clustering
HFBmgss <- NormalizeData(HFBmgss)
HFBmgss <- FindVariableFeatures(HFBmgss)
HFBmgss <- RunFastMNN(object.list = SplitObject(HFBmgss, split.by = "ID"))
HFBmgss <- RunUMAP(HFBmgss, reduction = "mnn", dims = 1:18)
HFBmgss <- FindNeighbors(HFBmgss, reduction = "mnn", dims = 1:18)
HFBmgss <- FindClusters(HFBmgss, resolution = 0.4)
DimPlot(HFBmgss, label = T)
FeaturePlot(HFBmgss, features = c("CTHRC1", "NPNT", "SFRP2", "SFRP4"))

##highlight transferred cluster meta data from Tsukui et al.
Meta_Highlight_Plot(seurat_object = HFBmgss, meta_data_column = "Tsukui_clusters", meta_data_highlight = "Fibrotic", highlight_color = "#E76BF3", background_color = "lightgray")
Meta_Highlight_Plot(seurat_object = HFBmgss, meta_data_column = "Tsukui_clusters", meta_data_highlight = "Inflammatory1", highlight_color = "#00BF7D", background_color = "lightgray")
Meta_Highlight_Plot(seurat_object = HFBmgss, meta_data_column = "Tsukui_clusters", meta_data_highlight = "Inflammatory2", highlight_color = "#00B0F6", background_color = "lightgray")
Meta_Highlight_Plot(seurat_object = HFBmgss, meta_data_column = "Tsukui_clusters", meta_data_highlight = "Alveolar", highlight_color = "#F8766D", background_color = "lightgray")

##annotation
new.cluster.ids <- c("Alveolar", "Inflammatory1", "Fibrotic", "Inflammatory1", "Inflammatory2", "Alveolar")
names(new.cluster.ids) <- levels(HFBmgss)
HFBmgss <- RenameIdents(HFBmgss, new.cluster.ids)
saveRDS(HFBmgss, file = "HFBmgss.rds")
