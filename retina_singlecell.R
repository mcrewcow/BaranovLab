#Install and load the packages for the analysis

library(Seurat)
library(patchwork)
library(dplyr)
library(ggplot2)
library(cowplot)
library(reticulate)

library(EnhancedVolcano)
library(DoubletFinder)
library(CytoTRACE)
library(biomaRt)
library(SeuratDisk)
library(SeuratWrappers)
library(Matrix)

data_dir <- 'C:/Users/Emil/BaranovLab/Retina/D59_fetal_filtered_gene_bc_matrices/GRCh38'
list.files(data_dir)
D59fetal <- Read10X(data.dir = data_dir)
D59fetalS = CreateSeuratObject(counts = D59fetal)

s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

#Create the function for creating and processing the Seurat objects

ProcessSeu <- function(Seurat){
Seurat <- NormalizeData(Seurat)
Seurat <- FindVariableFeatures(Seurat, selection.method = "vst", nfeatures = 2000)
Seurat <- ScaleData(Seurat)
Seurat <- RunPCA(Seurat)
Seurat <- FindNeighbors(Seurat, dims = 1:10)
Seurat <- FindClusters(Seurat, resolution = 0.5)
Seurat <- RunUMAP(Seurat, dims = 1:10)
Seurat <- RunTSNE(Seurat,  dims.use = 1:10 )
DimPlot(object = Seurat, reduction = "umap")
return (Seurat)
}

D59fetalS[["percent.rb"]] <- PercentageFeatureSet(D59fetalS, pattern = "^RPS|^RPL|^MRPS|^MRPL", assay = 'RNA')
D59fetalS[["percent.mt"]] <- PercentageFeatureSet(D59fetalS, pattern = "^MT-")
D59fetalS <- CellCycleScoring(D59fetalS, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
D59fetalS <- subset(D59fetalS, subset = nCount_RNA > 1000 & nCount_RNA < 20000 & nFeature_RNA > 600 & nFeature_RNA < 6000 & percent.mt < 15 & percent.rb < 30)
D59fetalS <- ScaleData(D59fetalS, verbose = T, vars.to.regress = c('nCount_RNA', 'percent.mt', "percent.rb","S.Score","G2M.Score"))
D59fetalS <- ProcessSeu(D59fetalS)

FeaturePlot(D59fetalS, features = c('RBPMS'))
DimPlot(D59fetalS, reduction = 'umap', label = TRUE, repel = TRUE)


###Metadata input
#D59metadata <- read.table(file = "clipboard", sep = "\t", header=TRUE)
#rownames(D59metadata) <- D59metadata$barcode
#D59fetalS1 <- AddMetaData(object = D59fetalS1, metadata = D59metadata)
#D59fetalS2 <- subset(D59fetalS2, subset = type != 'NA')
#D59fetalS2 <- RenameIdents(object = D59fetalS2, '0' = 'RGC', '1' = 'T1', '2' = 'RGC','3' = 'Progenitors', '4' = 'T1', '5' = 'Amacrine', '6' = 'RGC', '7' = 'Progenitors', '8' = 'Progenitors')

D59fetalS.markers <- FindAllMarkers(D59fetalS, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

D59fetalS.markers %>%
     group_by(cluster) %>%
     slice_max(n = 2, order_by = avg_log2FC)

D59fetalS.markers %>%
     group_by(cluster) %>%
     top_n(n = 10, wt = avg_log2FC) -> top10
 DoHeatmap(D59fetalS, features = top10$gene) + NoLegend()

#For merging P and C datasets
D82PCfetalS <- merge(D82cfetalS , y = D82pfetalS)

#FD59 clusters annotation
#D59fetalS1 <- RenameIdents(object = D59fetalS1, '0' = 'RGC', '1' = 'T1A', '2' = 'RGC','3' = 'Progenitors', '4' = 'T1B', '5' = 'AC/HC', '6' = 'RGC', '7' = 'Progenitors', '8' = 'Progenitors')

ProcessInt(my_data_frame) #This line is to choose the dataset for Integrated Processing in Seurat

ProcessInt <- function(data.integrated){
data.integrated <- ScaleData(data.integrated, verbose = FALSE)
data.integrated <- RunPCA(data.integrated, npcs = 30, verbose = FALSE)
data.integrated <- FindNeighbors(data.integrated, dims = 1:20)
data.integrated <- FindClusters(data.integrated, resolution = 0.5)
data.integrated <- RunUMAP(data.integrated, reduction = "pca", dims = 1:20)
data.integrated <- RunTSNE(data.integrated,  dims.use = 1:10 )
}

#Choose the objects for integration

integration_list <- list(poc, pocGH, soc, socGH, pocV, socV)

features <- SelectIntegrationFeatures(object.list = integration_list)
data.anchors <- FindIntegrationAnchors(object.list = integration_list, anchor.features = features)

data.combined <- IntegrateData(anchorset = data.anchors)

pocpocGHsocsocGHpocVsocV <- ProcessInt(data.combined)
