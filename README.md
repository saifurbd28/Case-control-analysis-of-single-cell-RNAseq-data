# Case-control-analysis-of-single-cell-RNAseq-data
R-Code of scRNA-seq data analysis between WT and KO mice (n=3/group) 

# Dataset: Publicly available data, extracted from NCBI, GEO Accession GSE212546
# Citation: Liu XH, Zhou JT, Yan CX, Cheng C, Fan JN, Xu J, Zheng Q, Bai Q, Li Z, Li S, Li X. Single-cell RNA sequencing reveals a novel inhibitory effect of ApoA4 on NAFL mediated by liver-specific subsets of myeloid cells. Front Immunol. 2022 Nov 8;13:1038401. doi: 10.3389/fimmu.2022.1038401. PMID: 36426356; PMCID: PMC9678944.

# Summary of the dataset: scRNA-seq on liver immune cells from WT and ApoA4-deficient mice administered a high-fat diet. 
# Overall design: 	Liver immune cells isolated from WT and KO mice (n=3/group) were subjected to scRNA-seq

# Environment
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(cowplot)
library(tidyr)

# Setting up the working directory
setwd("~/GSE212546_Liver_WT_KO_NAFLD")

# Load the datasets
control <- Read10X("~/GSE212546_Liver_WT_KO_NAFLD/WT")
case <- Read10X("~/GSE212546_Liver_WT_KO_NAFLD/KO")

# Create Seurat objects for control and case datasets
control <- CreateSeuratObject(counts = control, project = "Control", min.cells = 3, min.features = 200)
case <- CreateSeuratObject(counts = case, project = "Case", min.cells = 3, min.features = 200)

# Merge the two Seurat objects
merged <- merge(x = control, y = case, add.cell.ids = c("Control", "Case"))

# Perform QC and store the 'percent.mt' in the metadata
merged[["percent.mt"]] <- PercentageFeatureSet(merged, pattern = "^Mt-") # In case of human data, please use "^MT-" instead of "^Mt-"

# VlnPlot
VlnPlot(merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

![image](https://github.com/saifurbd28/Case-control-analysis-of-single-cell-RNAseq-data/assets/100442163/0e40d81f-22bc-4b62-a6ed-4fc5986e48d5)

# FeatureScatter plots
plot1 <- FeatureScatter(merged, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(merged, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# Combine the two FeatureScatter plots using 'patchwork'
plot1 + plot2

![image](https://github.com/saifurbd28/Case-control-analysis-of-single-cell-RNAseq-data/assets/100442163/ebd87539-f289-4d98-9a1f-b04cb54a51fe)

# Filtering -----------------
merged <- subset(merged, subset = nFeature_RNA > 50 & nFeature_RNA < 2500 & 
                 percent.mt < 0.05)

# Normalizing the data-------------
merged <- NormalizeData(merged)

#Identification of highly variable features------------
merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(merged), 10)
plot1 <- VariableFeaturePlot(merged)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

![image](https://github.com/saifurbd28/Case-control-analysis-of-single-cell-RNAseq-data/assets/100442163/d1d14b0b-841a-42f7-bf72-4a9133441371)

#Scaling the data-------------
merged <- ScaleData(merged)

#Perform linear dimensional reduction---------------
merged <- RunPCA(merged, features = VariableFeatures(merged))
print(merged[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(merged, dims = 1:2, reduction = "pca")

![image](https://github.com/saifurbd28/Case-control-analysis-of-single-cell-RNAseq-data/assets/100442163/325d4da8-8330-452c-a88e-d44ce033467c)

DimPlot(merged, reduction = "pca")

![image](https://github.com/saifurbd28/Case-control-analysis-of-single-cell-RNAseq-data/assets/100442163/0d3bd939-c11e-4b50-b8f6-898bdfcfb3d5)
DimHeatmap(merged, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(merged, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(merged)

![image](https://github.com/saifurbd28/Case-control-analysis-of-single-cell-RNAseq-data/assets/100442163/1133c586-fa5a-4831-afa9-97f71f34ed8a)

#Cluster the cells----------------------
merged <- FindNeighbors(merged, dims = 1:10)
merged <- FindClusters(merged, resolution = 0.1)
head(Idents(merged), 5)

merged <- RunUMAP(merged, dims = 1:10)
DimPlot(merged, reduction = "umap")

![image](https://github.com/saifurbd28/Case-control-analysis-of-single-cell-RNAseq-data/assets/100442163/167ec997-fedf-44b3-826a-c3efb55f37f9)
DimPlot(merged, reduction = "umap", label = TRUE)

![image](https://github.com/saifurbd28/Case-control-analysis-of-single-cell-RNAseq-data/assets/100442163/57828fb8-4a59-465d-aa5c-c7785286fbc3)

# find all markers of cluster 5
cluster5.markers <- FindMarkers(merged, ident.1 = 5, min.pct = 0.25)
head(cluster5.markers, n = 1000)
write.csv(cluster5.markers, file = "cluster5_markers.csv", row.names = TRUE)
# visualize top 2 markers in cluster5------------------
VlnPlot(merged, features = c(row.names(cluster5.markers)[1], row.names(cluster5.markers)[2]))

![image](https://github.com/saifurbd28/Case-control-analysis-of-single-cell-RNAseq-data/assets/100442163/9976638e-762f-4c4a-b3ee-43ab69b08292)

# find all markers distinguishing cluster 5 from all other clusters (i.e., 0, 1,2,3,4)
cluster5vsall.markers <- FindMarkers(merged, ident.1 = 5, ident.2 = c(0,1,2,3,4), min.pct = 0.25)
head(cluster5vsall.markers, n = 10)
write.csv(cluster5vsall.markers, file = "cluster5vsall.markers.csv", row.names = TRUE)

# vinplot for any markers
VlnPlot(merged, features = c("Cd68", "Csf1r"))

![image](https://github.com/saifurbd28/Case-control-analysis-of-single-cell-RNAseq-data/assets/100442163/92d01331-65ef-4828-80d9-2104ef8a343f)

FeaturePlot(merged, features = c("Cd68", "Csf1r"))

![image](https://github.com/saifurbd28/Case-control-analysis-of-single-cell-RNAseq-data/assets/100442163/7c6421d2-d84e-4c61-95e8-0796eafdc67e)

# Find the markers for cluster 5 and heatmap
cluster5.markers <- FindMarkers(object = merged, ident.1 = "5", min.pct = 0.25)
print(colnames(cluster5.markers))# Print the column names of the cluster6.markers data frame
top_genes_cluster5 <- cluster5.markers %>%
  top_n(n = 100, wt = avg_log2FC) %>%
  rownames()# Extract the gene names using row names
seurat_obj_cluster5 <- subset(merged, idents = "5", features = top_genes_cluster5)# Subset your Seurat object to cluster 5 and selected genes
heatmap_plot_cluster5 <- DoHeatmap(seurat_obj_cluster5, features = top_genes_cluster5) + NoLegend()# Create a heatmap for cluster 5 using top markers
print(heatmap_plot_cluster5)# Print the heatmap for cluster 5

![image](https://github.com/saifurbd28/Case-control-analysis-of-single-cell-RNAseq-data/assets/100442163/ed310255-59c7-45bf-9a26-6862e575b795)






