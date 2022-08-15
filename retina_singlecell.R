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

library(dittoSeq)

data_dir <- 'C:/Users/Emil/BaranovLab/Retina/D59_fetal_filtered_gene_bc_matrices/GRCh38'
list.files(data_dir)
D59fetal <- Read10X(data.dir = data_dir)
D59fetalS = CreateSeuratObject(counts = D59fetal)

memory.limit(size=56000)

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

#For saving the Seurat objects
SaveH5Seurat(D125PCfetalS1, 'C://Users/Emil/10X/scretina/FD125.h5Seurat', overwrite = TRUE)
seuratobject <- LoadH5Seurat('C://Users/Emil/10X/scretina/FD125.h5Seurat')

#For subsetting the clusters / identification
D82PCfetalS11 <- subset(D82PCfetalS, subset = seurat_clusters == 11)
names <- colnames(D82PCfetalS11)
D82PCfetalS <- SetIdent(D82PCfetalS, cells = names, value = 'Amacrine')
DimPlot(D82PCfetalS, reduction = 'umap', label = TRUE, label.box = TRUE)

#For subcluster/expression analysis
HLAexpr <- GetAssayData(adultretinaSsub, assay = 'RNA', slot = 'data')['HLA-DRA',]
poshla <- names(which(HLAexpr>0))
neghla <- names(which(HLAexpr=0))
poscells <- subset(adultretinaSsub, cells = poshla)
negcells <- subset(adultretinaSsub, cells = neghla)
adultretinaS1 <- SetIdent(adultretinaS1, cells = poshla, value = 'Microglia')

library(htmlwidgets)
library(plotly)
library(ggplotly)

#For Violin Plots
#Subset RGCs and merge the data adding the $stage metadata
RGCCXCR <- merge(D45orgSrgc, y = c(D60orgSrgc, D59fetalSRGC, D82PCfetalSRGC, D125PCfetalSRGC, adultretinaSRGC))
RGCCXCR$stage <- factor(RGCCXCR$stage, levels = c('FD59','FD82','FD125','Adult','OD45','OD60'))
library(ggplot2)
vln_df = data.frame(CXCR4 = RGCCXCR[["RNA"]]@data["CXCR4",], cluster = RGCCXCR$stage)
ggplot(vln_df, aes(x = cluster, y = CXCR4)) + geom_violin(aes(fill = cluster), trim=TRUE, scale = "width") + geom_jitter(width = 0.2) + theme_minimal()

#For bars
ggplot(bars, aes(x = Stage, y = CXCR4, fill = Stage)) + geom_bar(stat = 'identity') + theme_minimal()
ggplot(bars, aes(x = Stage, y = CXCR4, color = Stage)) + geom_bar(stat = 'identity', fill = 'white', size = 1.2) +
 geom_text(aes(label = round(CXCR4, digits = 1)), vjust = -1, color = 'black', size = 3.5) + theme_minimal() 
ggplot(bars, aes(x = Stage, y = CXCR4, color = Stage)) + geom_bar(stat = 'identity', fill = 'white', size = 1.2) +
 geom_text(aes(label = Amount), vjust = -1, color = 'black', size = 3.5) + theme_minimal()

#For the dotplots merge the data with the $analysis metadata
#The genes include: CXCR4, DCC, PTCH1, ROBO1 etc
analysis <- RenameIdents(analysis, 'Rod' = 'PR','Cone' = 'PR', 'ON-bipolar' = 'Bipolar', 'OFF-cone bipolar' = 'Bipolar', 'GABA-Amacrine' = 'Amacrine')
CXCR.FD59 <- analysis[['RNA']]@data["CXCR4",] * (analysis$an == "FD59")
CXCR.FD82 <- analysis[['RNA']]@data["CXCR4",] * (analysis$an == "FD82")
CXCR.FD125 <- analysis[['RNA']]@data["CXCR4",] * (analysis$an == "FD125")
CXCR.Adult <- analysis[['RNA']]@data["CXCR4",] * (analysis$an == "Adult")
analysis[['NEW']] <- CreateAssayObject(data = rbind(CXCR.FD59, CXCR.FD82, CXCR.FD125, CXCR.Adult))
DefaultAssay(analysis) <- "NEW"
DotPlot(analysis, features = c("CXCR.FD59", "CXCR.FD82",'CXCR.FD125','CXCR.Adult'), dot.scale = 10) 

analysisrgcs <- subset(analysis, idents = c('RGC'))
#analysisrgcs$an <- factor(RGCCXCR$stage, levels = c('FD59','FD82','FD125','Adult'))
DotPlot(analysisrgcs, features = c('ACKR3','ADGRG1','ADGRL3','CCR4','CELSR1','CELSR2','CELSR3',
'CXCR4','DCC','DRD1','DRD2','ERBB4','ESR2','FGFR1','FZD3','GFRA3','GPR173','IL1R1','NR2F1','NR2F2',
'NR4A2','NSMF','NTRK2','PTCH1','PTPRZ1','ROBO1','ROBO2','ROBO3','UNC5C','UNC5D'), split.by = 'an', cols = c('blue','blue','blue','blue')) 
#Or cols = 'RdBu' for red/blue scale colour

#For FeatureScatter
RGCCXCR <- merge(D59fetalSRGC, y = c(D82PCfetalSRGC, D125PCfetalSRGC, adultretinaSRGC))
RGCCXCR$stage <- factor(RGCCXCR$stage, levels = c('FD59','FD82','FD125','Adult'))
FeatureScatter(RGCCXCR, feature1 = 'CXCR4', feature2 = 'RBPMS', group.by = 'stage', pt.size = 5)

#For the multiple gene analysis
mp_genes <- read.table(file = "clipboard", sep = "\t", header=TRUE)
mp_genes <- as.list(mp_genes)
mp_genes <- lapply(mp_genes, toupper)
mp_genes <- unlist(mp_genes)
mp <- unique(mp_genes)
DotPlot(analysisrgcs, features = mp, split.by = 'an', cols = c('blue','blue','blue','blue'))





