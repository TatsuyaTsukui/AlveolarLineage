library(dplyr)
library(Seurat)
library(patchwork)
library(DoubletFinder)

d0_1.data <- Read10X(data.dir = "E:\\Tatsuya\\TT3635_cellranger\\per_sample_outs\\Day0-1\\count\\matrix")
d0_1 <- CreateSeuratObject(counts = d0_1.data$`Gene Expression`, project = "d0_1", min.cells = 3, min.features = 200)
d0_1[["percent.mt"]] <- PercentageFeatureSet(d0_1, pattern = "^mt-")
d0_1 <- subset(d0_1, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 15)
d0_1 <- NormalizeData(d0_1)
d0_1 <- FindVariableFeatures(d0_1, selection.method = "vst", nfeatures = 2000)
d0_1 <- ScaleData(d0_1)
d0_1 <- RunPCA(d0_1)
d0_1 <- RunUMAP(d0_1, dims = 1:15)
d0_1 <- FindNeighbors(d0_1, dims = 1:15)
d0_1 <- FindClusters(d0_1, resolution = 0.4)

sweep.res.list_d0_1 <- paramSweep_v3(d0_1, PCs = 1:15, sct = FALSE)
sweep.stats_d0_1 <- summarizeSweep(sweep.res.list_d0_1, GT = FALSE)
bcmvn_d0_1 <- find.pK(sweep.stats_d0_1)

annotations <- d0_1@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations) 
nExp_poi <- round(0.02*nrow(d0_1@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

d0_1 <- doubletFinder_v3(d0_1, PCs = 1:15, pN = 0.25, pK = 0.28, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
d0_1 <- doubletFinder_v3(d0_1, PCs = 1:15, pN = 0.25, pK = 0.28, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.28_81", sct = FALSE)

d0_1 <- subset(d0_1, subset = DF.classifications_0.25_0.28_62 == "Singlet")
saveRDS(d0_1, file = "d0_1.rds")


d0_2 <- Read10X(data.dir = "E:\\Tatsuya\\TT3635_cellranger\\per_sample_outs\\Day0-2\\count\\matrix")
d0_2 <- CreateSeuratObject(counts = d0_2$`Gene Expression`, project = "d0_2", min.cells = 3, min.features = 200)
d0_2[["percent.mt"]] <- PercentageFeatureSet(d0_2, pattern = "^mt-")
d0_2 <- subset(d0_2, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 15)
d0_2 <- NormalizeData(d0_2)
d0_2 <- FindVariableFeatures(d0_2, selection.method = "vst", nfeatures = 2000)
d0_2 <- ScaleData(d0_2)
d0_2 <- RunPCA(d0_2)
ElbowPlot(d0_2)

d0_2 <- RunUMAP(d0_2, dims = 1:17)
d0_2 <- FindNeighbors(d0_2, dims = 1:17)
d0_2 <- FindClusters(d0_2, resolution = 0.4)

sweep.res.list_d0_2 <- paramSweep_v3(d0_2, PCs = 1:17, sct = FALSE)
sweep.stats_d0_2 <- summarizeSweep(sweep.res.list_d0_2, GT = FALSE)
bcmvn_d0_2 <- find.pK(sweep.stats_d0_2)

annotations <- d0_2@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations) 
nExp_poi <- round(0.02*nrow(d0_2@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

d0_2 <- doubletFinder_v3(d0_2, PCs = 1:17, pN = 0.25, pK = 0.01, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
d0_2 <- doubletFinder_v3(d0_2, PCs = 1:17, pN = 0.25, pK = 0.01, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.01_84", sct = FALSE)

DimPlot(d0_2, group.by = "DF.classifications_0.25_0.01_64")

d0_2 <- subset(d0_2, subset = DF.classifications_0.25_0.01_64 == "Singlet")
saveRDS(d0_2, file = "d0_2.rds")


d0_3 <- Read10X(data.dir = "E:\\Tatsuya\\TT3635_cellranger\\per_sample_outs\\Day0-3\\count\\matrix")
d0_3 <- CreateSeuratObject(counts = d0_3$`Gene Expression`, project = "d0_3", min.cells = 3, min.features = 200)
d0_3[["percent.mt"]] <- PercentageFeatureSet(d0_3, pattern = "^mt-")
d0_3 <- subset(d0_3, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 15)
d0_3 <- NormalizeData(d0_3)
d0_3 <- FindVariableFeatures(d0_3, selection.method = "vst", nfeatures = 2000)
d0_3 <- ScaleData(d0_3)
d0_3 <- RunPCA(d0_3)
ElbowPlot(d0_3)

d0_3 <- RunUMAP(d0_3, dims = 1:17)
d0_3 <- FindNeighbors(d0_3, dims = 1:17)
d0_3 <- FindClusters(d0_3, resolution = 0.4)

sweep.res.list_d0_3 <- paramSweep_v3(d0_3, PCs = 1:17, sct = FALSE)
sweep.stats_d0_3 <- summarizeSweep(sweep.res.list_d0_3, GT = FALSE)
bcmvn_d0_3 <- find.pK(sweep.stats_d0_3)

annotations <- d0_3@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations) 
nExp_poi <- round(0.02*nrow(d0_3@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

d0_3 <- doubletFinder_v3(d0_3, PCs = 1:17, pN = 0.25, pK = 0.01, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
d0_3 <- doubletFinder_v3(d0_3, PCs = 1:17, pN = 0.25, pK = 0.01, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.01_69", sct = FALSE)

DimPlot(d0_3, group.by = "DF.classifications_0.25_0.01_50")

d0_3 <- subset(d0_3, subset = DF.classifications_0.25_0.01_50 == "Singlet")
saveRDS(d0_3, file = "d0_3.rds")


d7_1 <- Read10X(data.dir = "E:\\Tatsuya\\TT3635_cellranger\\per_sample_outs\\Day7-1\\count\\matrix")
d7_1 <- CreateSeuratObject(counts = d7_1$`Gene Expression`, project = "d7_1", min.cells = 3, min.features = 200)
d7_1[["percent.mt"]] <- PercentageFeatureSet(d7_1, pattern = "^mt-")
d7_1 <- subset(d7_1, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 15)
d7_1 <- NormalizeData(d7_1)
d7_1 <- FindVariableFeatures(d7_1, selection.method = "vst", nfeatures = 2000)
d7_1 <- ScaleData(d7_1)
d7_1 <- RunPCA(d7_1)
ElbowPlot(d7_1)

d7_1 <- RunUMAP(d7_1, dims = 1:20)
d7_1 <- FindNeighbors(d7_1, dims = 1:20)
d7_1 <- FindClusters(d7_1, resolution = 0.4)
DimPlot(d7_1)

sweep.res.list_d7_1 <- paramSweep_v3(d7_1, PCs = 1:20, sct = FALSE)
sweep.stats_d7_1 <- summarizeSweep(sweep.res.list_d7_1, GT = FALSE)
bcmvn_d7_1 <- find.pK(sweep.stats_d7_1)

annotations <- d7_1@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations) 
nExp_poi <- round(0.02*nrow(d7_1@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

d7_1 <- doubletFinder_v3(d7_1, PCs = 1:20, pN = 0.25, pK = 0.2, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
d7_1 <- doubletFinder_v3(d7_1, PCs = 1:20, pN = 0.25, pK = 0.2, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.2_39", sct = FALSE)

DimPlot(d7_1, group.by = "DF.classifications_0.25_0.2_33")

d7_1 <- subset(d7_1, subset = DF.classifications_0.25_0.2_33 == "Singlet")
saveRDS(d7_1, file = "d7_1.rds")


d7_2 <- Read10X(data.dir = "E:\\Tatsuya\\TT3635_cellranger\\per_sample_outs\\Day7-2\\count\\matrix")
d7_2 <- CreateSeuratObject(counts = d7_2$`Gene Expression`, project = "d7_2", min.cells = 3, min.features = 200)
d7_2[["percent.mt"]] <- PercentageFeatureSet(d7_2, pattern = "^mt-")
d7_2 <- subset(d7_2, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 15)
d7_2 <- NormalizeData(d7_2)
d7_2 <- FindVariableFeatures(d7_2, selection.method = "vst", nfeatures = 2000)
d7_2 <- ScaleData(d7_2)
d7_2 <- RunPCA(d7_2)
ElbowPlot(d7_2)

d7_2 <- RunUMAP(d7_2, dims = 1:20)
d7_2 <- FindNeighbors(d7_2, dims = 1:20)
d7_2 <- FindClusters(d7_2, resolution = 0.4)
DimPlot(d7_2)

sweep.res.list_d7_2 <- paramSweep_v3(d7_2, PCs = 1:20, sct = FALSE)
sweep.stats_d7_2 <- summarizeSweep(sweep.res.list_d7_2, GT = FALSE)
bcmvn_d7_2 <- find.pK(sweep.stats_d7_2)
bcmvn_d7_2

annotations <- d7_2@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations) 
nExp_poi <- round(0.02*nrow(d7_2@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

d7_2 <- doubletFinder_v3(d7_2, PCs = 1:20, pN = 0.25, pK = 0.005, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
d7_2 <- doubletFinder_v3(d7_2, PCs = 1:20, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_60", sct = FALSE)

DimPlot(d7_2, group.by = "DF.classifications_0.25_0.005_51")

d7_2 <- subset(d7_2, subset = DF.classifications_0.25_0.005_51 == "Singlet")
saveRDS(d7_2, file = "d7_2.rds")


d7_3 <- Read10X(data.dir = "E:\\Tatsuya\\TT3635_cellranger\\per_sample_outs\\Day7-3\\count\\matrix")
d7_3 <- CreateSeuratObject(counts = d7_3$`Gene Expression`, project = "d7_3", min.cells = 3, min.features = 200)
d7_3[["percent.mt"]] <- PercentageFeatureSet(d7_3, pattern = "^mt-")
d7_3 <- subset(d7_3, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 15)
d7_3 <- NormalizeData(d7_3)
d7_3 <- FindVariableFeatures(d7_3, selection.method = "vst", nfeatures = 2000)
d7_3 <- ScaleData(d7_3)
d7_3 <- RunPCA(d7_3)
ElbowPlot(d7_3)

d7_3 <- RunUMAP(d7_3, dims = 1:20)
d7_3 <- FindNeighbors(d7_3, dims = 1:20)
d7_3 <- FindClusters(d7_3, resolution = 0.4)
DimPlot(d7_3)

sweep.res.list_d7_3 <- paramSweep_v3(d7_3, PCs = 1:20, sct = FALSE)
sweep.stats_d7_3 <- summarizeSweep(sweep.res.list_d7_3, GT = FALSE)
bcmvn_d7_3 <- find.pK(sweep.stats_d7_3)
bcmvn_d7_3

annotations <- d7_3@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations) 
nExp_poi <- round(0.02*nrow(d7_3@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

d7_3 <- doubletFinder_v3(d7_3, PCs = 1:20, pN = 0.25, pK = 0.06, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
d7_3 <- doubletFinder_v3(d7_3, PCs = 1:20, pN = 0.25, pK = 0.06, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.06_86", sct = FALSE)

DimPlot(d7_3, group.by = "DF.classifications_0.25_0.06_73")

d7_3 <- subset(d7_3, subset = DF.classifications_0.25_0.06_73 == "Singlet")
saveRDS(d7_3, file = "d7_3.rds")


d14_1 <- Read10X(data.dir = "E:\\Tatsuya\\TT3635_cellranger\\per_sample_outs\\Day14-1\\count\\matrix")
d14_1 <- CreateSeuratObject(counts = d14_1$`Gene Expression`, project = "d14_1", min.cells = 3, min.features = 200)
d14_1[["percent.mt"]] <- PercentageFeatureSet(d14_1, pattern = "^mt-")
d14_1 <- subset(d14_1, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 15)
d14_1 <- NormalizeData(d14_1)
d14_1 <- FindVariableFeatures(d14_1, selection.method = "vst", nfeatures = 2000)
d14_1 <- ScaleData(d14_1)
d14_1 <- RunPCA(d14_1)
ElbowPlot(d14_1, ndims = 30)

d14_1 <- RunUMAP(d14_1, dims = 1:21)
d14_1 <- FindNeighbors(d14_1, dims = 1:21)
d14_1 <- FindClusters(d14_1, resolution = 0.4)
DimPlot(d14_1)

sweep.res.list_d14_1 <- paramSweep_v3(d14_1, PCs = 1:21, sct = FALSE)
sweep.stats_d14_1 <- summarizeSweep(sweep.res.list_d14_1, GT = FALSE)
bcmvn_d14_1 <- find.pK(sweep.stats_d14_1)
bcmvn_d14_1

annotations <- d14_1@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations) 
nExp_poi <- round(0.02*nrow(d14_1@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

d14_1 <- doubletFinder_v3(d14_1, PCs = 1:21, pN = 0.25, pK = 0.05, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
d14_1 <- doubletFinder_v3(d14_1, PCs = 1:21, pN = 0.25, pK = 0.05, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.05_86", sct = FALSE)

DimPlot(d14_1, group.by = "DF.classifications_0.25_0.05_62")

d14_1 <- subset(d14_1, subset = DF.classifications_0.25_0.05_62 == "Singlet")
saveRDS(d14_1, file = "d14_1.rds")


d14_2 <- Read10X(data.dir = "E:\\Tatsuya\\TT3635_cellranger\\per_sample_outs\\Day14-2\\count\\matrix")
d14_2 <- CreateSeuratObject(counts = d14_2$`Gene Expression`, project = "d14_2", min.cells = 3, min.features = 200)
d14_2[["percent.mt"]] <- PercentageFeatureSet(d14_2, pattern = "^mt-")
d14_2 <- subset(d14_2, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 15)
d14_2 <- NormalizeData(d14_2)
d14_2 <- FindVariableFeatures(d14_2, selection.method = "vst", nfeatures = 2000)
d14_2 <- ScaleData(d14_2)
d14_2 <- RunPCA(d14_2)
ElbowPlot(d14_2, ndims = 30)

d14_2 <- RunUMAP(d14_2, dims = 1:21)
d14_2 <- FindNeighbors(d14_2, dims = 1:21)
d14_2 <- FindClusters(d14_2, resolution = 0.4)
DimPlot(d14_2)

sweep.res.list_d14_2 <- paramSweep_v3(d14_2, PCs = 1:21, sct = FALSE)
sweep.stats_d14_2 <- summarizeSweep(sweep.res.list_d14_2, GT = FALSE)
bcmvn_d14_2 <- find.pK(sweep.stats_d14_2)
bcmvn_d14_2

annotations <- d14_2@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations) 
nExp_poi <- round(0.02*nrow(d14_2@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

d14_2 <- doubletFinder_v3(d14_2, PCs = 1:21, pN = 0.25, pK = 0.07, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
d14_2 <- doubletFinder_v3(d14_2, PCs = 1:21, pN = 0.25, pK = 0.07, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.07_83", sct = FALSE)

DimPlot(d14_2, group.by = "DF.classifications_0.25_0.07_64")

d14_2 <- subset(d14_2, subset = DF.classifications_0.25_0.07_64 == "Singlet")
saveRDS(d14_2, file = "d14_2.rds")


d14_3 <- Read10X(data.dir = "E:\\Tatsuya\\TT3635_cellranger\\per_sample_outs\\Day14-3\\count\\matrix")
d14_3 <- CreateSeuratObject(counts = d14_3$`Gene Expression`, project = "d14_3", min.cells = 3, min.features = 200)
d14_3[["percent.mt"]] <- PercentageFeatureSet(d14_3, pattern = "^mt-")
d14_3 <- subset(d14_3, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 15)
d14_3 <- NormalizeData(d14_3)
d14_3 <- FindVariableFeatures(d14_3, selection.method = "vst", nfeatures = 2000)
d14_3 <- ScaleData(d14_3)
d14_3 <- RunPCA(d14_3)
ElbowPlot(d14_3, ndims = 30)

d14_3 <- RunUMAP(d14_3, dims = 1:23)
d14_3 <- FindNeighbors(d14_3, dims = 1:23)
d14_3 <- FindClusters(d14_3, resolution = 0.4)
DimPlot(d14_3)

sweep.res.list_d14_3 <- paramSweep_v3(d14_3, PCs = 1:23, sct = FALSE)
sweep.stats_d14_3 <- summarizeSweep(sweep.res.list_d14_3, GT = FALSE)
bcmvn_d14_3 <- find.pK(sweep.stats_d14_3)
bcmvn_d14_3

annotations <- d14_3@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations) 
nExp_poi <- round(0.02*nrow(d14_3@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

d14_3 <- doubletFinder_v3(d14_3, PCs = 1:23, pN = 0.25, pK = 0.07, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
d14_3 <- doubletFinder_v3(d14_3, PCs = 1:23, pN = 0.25, pK = 0.07, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.07_98", sct = FALSE)

DimPlot(d14_3, group.by = "DF.classifications_0.25_0.07_80")

d14_3 <- subset(d14_3, subset = DF.classifications_0.25_0.07_80 == "Singlet")
saveRDS(d14_3, file = "d14_3.rds")


d21_1 <- Read10X(data.dir = "E:\\Tatsuya\\TT3635_cellranger\\per_sample_outs\\Day21-1\\count\\matrix")
d21_1 <- CreateSeuratObject(counts = d21_1$`Gene Expression`, project = "d21_1", min.cells = 3, min.features = 200)
d21_1[["percent.mt"]] <- PercentageFeatureSet(d21_1, pattern = "^mt-")
d21_1 <- subset(d21_1, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 15)
d21_1 <- NormalizeData(d21_1)
d21_1 <- FindVariableFeatures(d21_1, selection.method = "vst", nfeatures = 2000)
d21_1 <- ScaleData(d21_1)
d21_1 <- RunPCA(d21_1)
ElbowPlot(d21_1, ndims = 30)

d21_1 <- RunUMAP(d21_1, dims = 1:22)
d21_1 <- FindNeighbors(d21_1, dims = 1:22)
d21_1 <- FindClusters(d21_1, resolution = 0.4)
DimPlot(d21_1)

sweep.res.list_d21_1 <- paramSweep_v3(d21_1, PCs = 1:22, sct = FALSE)
sweep.stats_d21_1 <- summarizeSweep(sweep.res.list_d21_1, GT = FALSE)
bcmvn_d21_1 <- find.pK(sweep.stats_d21_1)
bcmvn_d21_1

annotations <- d21_1@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations) 
nExp_poi <- round(0.02*nrow(d21_1@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

d21_1 <- doubletFinder_v3(d21_1, PCs = 1:22, pN = 0.25, pK = 0.23, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
d21_1 <- doubletFinder_v3(d21_1, PCs = 1:22, pN = 0.25, pK = 0.23, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.23_98", sct = FALSE)

DimPlot(d21_1, group.by = "DF.classifications_0.25_0.23_79")

d21_1 <- subset(d21_1, subset = DF.classifications_0.25_0.23_79 == "Singlet")
saveRDS(d21_1, file = "d21_1.rds")


d21_2 <- Read10X(data.dir = "E:\\Tatsuya\\TT3635_cellranger\\per_sample_outs\\Day21-2\\count\\matrix")
d21_2 <- CreateSeuratObject(counts = d21_2$`Gene Expression`, project = "d21_2", min.cells = 3, min.features = 200)
d21_2[["percent.mt"]] <- PercentageFeatureSet(d21_2, pattern = "^mt-")
d21_2 <- subset(d21_2, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 15)
d21_2 <- NormalizeData(d21_2)
d21_2 <- FindVariableFeatures(d21_2, selection.method = "vst", nfeatures = 2000)
d21_2 <- ScaleData(d21_2)
d21_2 <- RunPCA(d21_2)
ElbowPlot(d21_2, ndims = 30)

d21_2 <- RunUMAP(d21_2, dims = 1:22)
d21_2 <- FindNeighbors(d21_2, dims = 1:22)
d21_2 <- FindClusters(d21_2, resolution = 0.4)
DimPlot(d21_2)

sweep.res.list_d21_2 <- paramSweep_v3(d21_2, PCs = 1:22, sct = FALSE)
sweep.stats_d21_2 <- summarizeSweep(sweep.res.list_d21_2, GT = FALSE)
bcmvn_d21_2 <- find.pK(sweep.stats_d21_2)
bcmvn_d21_2

annotations <- d21_2@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations) 
nExp_poi <- round(0.02*nrow(d21_2@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

d21_2 <- doubletFinder_v3(d21_2, PCs = 1:22, pN = 0.25, pK = 0.3, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
d21_2 <- doubletFinder_v3(d21_2, PCs = 1:22, pN = 0.25, pK = 0.3, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.3_97", sct = FALSE)

DimPlot(d21_2, group.by = "DF.classifications_0.25_0.3_77")

d21_2 <- subset(d21_2, subset = DF.classifications_0.25_0.3_77 == "Singlet")
saveRDS(d21_2, file = "d21_2.rds")


d21_3 <- Read10X(data.dir = "E:\\Tatsuya\\TT3635_cellranger\\per_sample_outs\\Day21-3\\count\\matrix")
d21_3 <- CreateSeuratObject(counts = d21_3$`Gene Expression`, project = "d21_3", min.cells = 3, min.features = 200)
d21_3[["percent.mt"]] <- PercentageFeatureSet(d21_3, pattern = "^mt-")
d21_3 <- subset(d21_3, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 15)
d21_3 <- NormalizeData(d21_3)
d21_3 <- FindVariableFeatures(d21_3, selection.method = "vst", nfeatures = 2000)
d21_3 <- ScaleData(d21_3)
d21_3 <- RunPCA(d21_3)
ElbowPlot(d21_3, ndims = 30)

d21_3 <- RunUMAP(d21_3, dims = 1:23)
d21_3 <- FindNeighbors(d21_3, dims = 1:23)
d21_3 <- FindClusters(d21_3, resolution = 0.4)
DimPlot(d21_3)

sweep.res.list_d21_3 <- paramSweep_v3(d21_3, PCs = 1:23, sct = FALSE)
sweep.stats_d21_3 <- summarizeSweep(sweep.res.list_d21_3, GT = FALSE)
bcmvn_d21_3 <- find.pK(sweep.stats_d21_3)
bcmvn_d21_3

annotations <- d21_3@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations) 
nExp_poi <- round(0.02*nrow(d21_3@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

d21_3 <- doubletFinder_v3(d21_3, PCs = 1:23, pN = 0.25, pK = 0.05, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
d21_3 <- doubletFinder_v3(d21_3, PCs = 1:23, pN = 0.25, pK = 0.05, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.05_92", sct = FALSE)

DimPlot(d21_3, group.by = "DF.classifications_0.25_0.05_75")

d21_3 <- subset(d21_3, subset = DF.classifications_0.25_0.05_75 == "Singlet")
saveRDS(d21_3, file = "d21_3.rds")


