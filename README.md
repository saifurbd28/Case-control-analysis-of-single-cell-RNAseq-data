# Dissecting cell-type specific pathophysiological changes in human type-2 diabetic pancreatic islets using single-cell transcriptomic data
#R-Code of scRNA-seq data analysis between No Diabets (ND), pre-diabetic (PD), and type 2 diabetic (T2D) states (n=3/group) 

# Dataset: Publicly available data, extracted from NCBI, GEO Accession GSE221156

# Environment
library(dplyr)

library(Seurat)

library(patchwork)

library(cowplot)

library(magrittr)

library(ggplot2)

# Setting up the working directory
setwd("~/GSE221156_human type-2 diabetic pancreatic islets")

#  Load the dataset: Group-1
ND_1 <- Read10X("~/GSE221156_human type-2 diabetic pancreatic islets/GSE221156/ND/Islet29")

ND_2 <- Read10X("~/GSE221156_human type-2 diabetic pancreatic islets/GSE221156/ND/Islet34")

ND_3 <- Read10X("~/GSE221156_human type-2 diabetic pancreatic islets/GSE221156/ND/Islet37")

# Create Seurat objects for Group-1
ND_1 <- CreateSeuratObject(counts = ND_1, project = "ND", min.cells = 3, min.features = 200)

ND_1 <- PercentageFeatureSet(ND_1, pattern = "^MT", col.name = "percent.mt")

ND_2 <- CreateSeuratObject(counts = ND_2, project = "ND", min.cells = 3, min.features = 200)

ND_2 <- PercentageFeatureSet(ND_2, pattern = "^MT", col.name = "percent.mt")

ND_3 <- CreateSeuratObject(counts = ND_3, project = "ND", min.cells = 3, min.features = 200)

ND_3 <- PercentageFeatureSet(ND_3, pattern = "^MT", col.name = "percent.mt")

# Merge the all Seurat objects of Group-1

ND <- merge(x = ND_1, y = c(ND_2, ND_3), add.cell.ids = ls()[1:3], project = 'ND')

# Load the dataset: Group-2
PD_1 <- Read10X("~/GSE221156_human type-2 diabetic pancreatic islets/GSE221156/PD/Islet40")

PD_2 <- Read10X("~/GSE221156_human type-2 diabetic pancreatic islets/GSE221156/PD/Islet50")

PD_3 <- Read10X("~/GSE221156_human type-2 diabetic pancreatic islets/GSE221156/PD/Islet53")

# Create Seurat objects for Group-2
PD_1 <- CreateSeuratObject(counts = PD_1, project = "PD", min.cells = 3, min.features = 200)

PD_1 <- PercentageFeatureSet(PD_1, pattern = "^MT", col.name = "percent.mt")

PD_2 <- CreateSeuratObject(counts = PD_2, project = "PD", min.cells = 3, min.features = 200)

PD_2 <- PercentageFeatureSet(PD_2, pattern = "^MT", col.name = "percent.mt")

PD_3 <- CreateSeuratObject(counts = PD_3, project = "PD", min.cells = 3, min.features = 200)

PD_3 <- PercentageFeatureSet(PD_3, pattern = "^MT", col.name = "percent.mt")

# Merge the all Seurat objects of Group-2
PD <- merge(x = PD_1, y = c(PD_2, PD_3), add.cell.ids = ls()[1:3], project = 'PD')

# Load the dataset: Group-3
T2D_1 <- Read10X("~/GSE221156_human type-2 diabetic pancreatic islets/GSE221156/T2D/Islet28")

T2D_2 <- Read10X("~/GSE221156_human type-2 diabetic pancreatic islets/GSE221156/T2D/Islet30")

T2D_3 <- Read10X("~/GSE221156_human type-2 diabetic pancreatic islets/GSE221156/T2D/Islet31")

# Create Seurat objects for Group-3
T2D_1 <- CreateSeuratObject(counts = T2D_1, project = "T2D", min.cells = 3, min.features = 200)

T2D_1 <- PercentageFeatureSet(T2D_1, pattern = "^MT", col.name = "percent.mt")

T2D_2 <- CreateSeuratObject(counts = T2D_2, project = "T2D", min.cells = 3, min.features = 200)

T2D_2 <- PercentageFeatureSet(T2D_2, pattern = "^MT", col.name = "percent.mt")

T2D_3 <- CreateSeuratObject(counts = T2D_3, project = "T2D", min.cells = 3, min.features = 200)

T2D_3 <- PercentageFeatureSet(T2D_3, pattern = "^MT", col.name = "percent.mt")

# Merge the all Seurat objects of Group-3
T2D <- merge(x = T2D_1, y = c(T2D_2, T2D_3), add.cell.ids = ls()[1:3], project = 'T2D')

# Merge all Seurat objects 
merged <- merge(x = ND, y = c(PD, T2D), add.cell.ids = c("ND", "PD", "T2D"))

pbmc = merged #Rename

# This is a great place to start QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

pbmc <- subset(pbmc, subset = nFeature_RNA > 50 & nFeature_RNA < 10000 & percent.mt < 10)

VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
![image](https://github.com/saifurbd28/Case-control-analysis-of-single-cell-RNAseq-data/assets/100442163/b8509304-d0da-4e47-9d41-060c84ac9cca)

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")

plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

plot1 + plot2

#Normalization
pbmc <- NormalizeData(pbmc)

# Feature selection
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
![image](https://github.com/saifurbd28/Case-control-analysis-of-single-cell-RNAseq-data/assets/100442163/11d13ec1-4785-424a-99f7-c376bc67bae4)

#Scaling
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

#Perform linear dimensional reduction
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
![image](https://github.com/saifurbd28/Case-control-analysis-of-single-cell-RNAseq-data/assets/100442163/30648f9b-eb75-4430-9bdc-e9c0d71f0e39)

DimPlot(pbmc, reduction = "pca") + NoLegend()
![image](https://github.com/saifurbd28/Case-control-analysis-of-single-cell-RNAseq-data/assets/100442163/2dc2ef18-deeb-46b8-9908-1964d706a932)

DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 2, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:20, cells = 500, balanced = TRUE)
![image](https://github.com/saifurbd28/Case-control-analysis-of-single-cell-RNAseq-data/assets/100442163/3f7f785e-38c3-4f64-8cbe-027f9f47f1b0)


ElbowPlot(pbmc)
![image](https://github.com/saifurbd28/Case-control-analysis-of-single-cell-RNAseq-data/assets/100442163/47187378-c6a7-46ef-a373-dba6734c6eda)

pbmc <- FindNeighbors(pbmc, dims = 1:20)
pbmc <- FindClusters(pbmc, resolution = 0.1)
pbmc <- RunUMAP(pbmc, dims = 1:20)
DimPlot(pbmc, reduction = "umap", label = TRUE)

![image](https://github.com/saifurbd28/Case-control-analysis-of-single-cell-RNAseq-data/assets/100442163/2ff7068d-60db-407e-9aef-56eed1460010)


#To visualize the two conditions side-by-side
DimPlot(pbmc, reduction = "umap", split.by = "orig.ident", label = TRUE)

![image](https://github.com/saifurbd28/Case-control-analysis-of-single-cell-RNAseq-data/assets/100442163/ffca10e0-af56-4d37-b3e8-197517ca2ee7)

# Compute t-SNE reduction
pbmc <- RunTSNE(pbmc)
# Visualize t-SNE plot

DimPlot(pbmc, reduction = "tsne")

DimPlot(pbmc, reduction = "tsne", label = TRUE)

DimPlot(pbmc, reduction = "tsne", split.by = "orig.ident")

# Acinar cell + beta cell + alpha cell + Delta cell + Ductal cell + Apical cell + Stellate cell +
# pancreatic stem cell + Endothelial cell + Macrophage
features <- c("CPA1", "PRSS1",
              "INS", "CDKN1C", "C-peptide",
              "GCG", "GBA",
              "SST",
              "KRT19",
              "ZO1",
              "COL1A1",
              "NANOG", 
              "VWF",
              "SDS", "CD68")
DotPlot(pbmc, features = features, cols = c("green", "brown"),col.min = -5, col.max = 5 ) 
![image](https://github.com/saifurbd28/Case-control-analysis-of-single-cell-RNAseq-data/assets/100442163/73c725a8-ee74-4074-baf5-c49c0d3c6f42)


# Remane the clusters based on cell markers using RenameIdents function
![image](https://github.com/saifurbd28/Case-control-analysis-of-single-cell-RNAseq-data/assets/100442163/b0370451-779e-46a0-bb83-28f9a03baddc)


####################################################

cluster_counts <- table(Idents(pbmc))#Count the number of cells in each cluster

print(cluster_counts)

![image](https://github.com/saifurbd28/Case-control-analysis-of-single-cell-RNAseq-data/assets/100442163/7ad478a0-8581-4fdd-ac24-78e362f95bd0)


#To visualize one condition from the all conditions 

pbmc_meta.data <- pbmc@meta.data

write.csv(pbmc_meta.data, file = "pbmc_cluster cell count.csv", row.names = TRUE) #from this excel file you will get the cluster cell number for each cluster and each condition

################## The cluster-specific markers ################################
###############################################################################

pbmc1 <- JoinLayers(pbmc, features = c("assay1", "assay2", "assay3"))

# find top 10 markers for every cluster compared, report only the positive ones

pbmc.markers <- FindAllMarkers(pbmc1, only.pos = TRUE)

data.markers<-pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

write.csv(data.markers, file = "cluster.markers.csv", row.names = TRUE)

#Heatmap for all clusters

pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top10

DoHeatmap(pbmc, features = top10$gene) + NoLegend()

#colour option-1
DoHeatmap(object = pbmc, features = top10$gene) + 
  scale_fill_gradientn(colors = RColorBrewer::brewer.pal(n = 10, name = "RdBu"))

![image](https://github.com/saifurbd28/Case-control-analysis-of-single-cell-RNAseq-data/assets/100442163/70f475e2-0a40-49ec-9447-aaf8d53373ef)


#colour option-2
DoHeatmap(object = pbmc, features = top10$gene) + scale_fill_gradientn(colors = colorRampPalette(c("#0200ad", "#fbfcbd", "#ff0000"))(256))

![image](https://github.com/saifurbd28/Case-control-analysis-of-single-cell-RNAseq-data/assets/100442163/0c2953e8-b0df-46f5-836e-566fb4772620)

#colour option-2
DoHeatmap(object = pbmc, features = top10$gene) + scale_fill_gradientn(colors = colorRampPalette(c("dodgerblue3", "snow1", "brown3"))(256))

![image](https://github.com/saifurbd28/Case-control-analysis-of-single-cell-RNAseq-data/assets/100442163/a78061bd-9ed8-4b30-a671-45b4a5a9ad12)


########### Find differential expression between two conditions (e.g., T2D and ND) in a specific cluster-1###################
# Example with Cluster_0

Cluster_0 <- subset(x = pbmc, idents = c("0"))

DimPlot(Cluster_0, reduction = "umap", label = TRUE)

Idents(Cluster_0) <- "orig.ident" 

# T2D vs ND

merged <- JoinLayers(object = Cluster_0, features = c("assay3", "assay1"))

DE.Regression <- FindMarkers(merged, ident.1 = "T2D", ident.2 = "ND")

write.csv(DE.Regression, file = "DE.T2D vs ND _clustr0_alpha cell.csv", row.names = TRUE)


# T2D vs PD

merged <- JoinLayers(object = Cluster_0, features = c("assay3", "assay2"))

DE.Regression <- FindMarkers(merged, ident.1 = "T2D", ident.2 = "PD")

write.csv(DE.Regression, file = "DE.T2D vs PD _clustr0_alpha cell.csv", row.names = TRUE)


# PD vs ND

merged <- JoinLayers(object = Cluster_0, features = c("assay2", "assay1"))

DE.Regression <- FindMarkers(merged, ident.1 = "PD", ident.2 = "ND")

write.csv(DE.Regression, file = "DE.PD vs ND _clustr0_alpha cell.csv", row.names = TRUE)








