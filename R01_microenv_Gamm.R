expression_matrix <- read.table('G://Cones_GAMM//GSM4271906_WT_d170-3.dge.txt', , header = TRUE, row.names = 1)
expr1 <- CreateSeuratObject(counts = expression_matrix, project = '1')

s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

ProcessSeu <- function(Seurat){
  Seurat <- NormalizeData(Seurat)
  Seurat <- FindVariableFeatures(Seurat, selection.method = "vst", nfeatures = 2000)
  Seurat <- ScaleData(Seurat)
  Seurat <- RunPCA(Seurat)
  Seurat <- FindNeighbors(Seurat, dims = 1:20)
  Seurat <- FindClusters(Seurat, resolution = 0.5)
  Seurat <- RunUMAP(Seurat, dims = 1:20)
  Seurat <- RunTSNE(Seurat,  dims.use = 1:20 )
  DimPlot(object = Seurat, reduction = "umap")
  return (Seurat)
}

expr1 [["percent.rb"]] <- PercentageFeatureSet(expr1 , pattern = "^RPS|^RPL|^MRPS|^MRPL", assay = 'RNA')
expr1 <- CellCycleScoring(expr1, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE, nbin = 3)
expr1[["percent.mt"]] <- PercentageFeatureSet(expr1, pattern = "^MT-")
VlnPlot(expr1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", 'percent.rb'), ncol = 4)

expr1 <- subset(expr1, subset = nCount_RNA > 100 & nCount_RNA < 10000 & nFeature_RNA > 100 & nFeature_RNA < 1500 & percent.mt < 15 & percent.rb < 35)

expr1 <- ScaleData(expr1, verbose = T, vars.to.regress = c('nCount_RNA', 'nFeature_RNA', 'percent.mt', "percent.rb","S.Score","G2M.Score"))

expr1 <- ProcessSeu(expr1)

RDoublet <- function(tmp){
  sweep.res.list <- paramSweep_v3(tmp, PCs = 1:30, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pKopt <- as.numeric(as.character(bcmvn$pK[bcmvn$BCmetric == max(bcmvn$BCmetric)]))
  pKopt <- pKopt[order(pKopt, decreasing = TRUE) ]
  pKopt <- pKopt[1]
  homotypic.prop <- modelHomotypic(tmp$seurat_clusters) 
  nExp_poi <- round(0.1*length(colnames(tmp)))  ## Assuming 10% doublet formation rate 
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  tmp <- doubletFinder_v3(tmp, PCs = 1:30, pN = 0.25, pK = pKopt, nExp = nExp_poi, reuse.pANN = FALSE)
  tmp <- doubletFinder_v3(tmp, PCs = 1:30, pN = 0.25, pK = pKopt, nExp = nExp_poi.adj, reuse.pANN = paste("pANN_0.25",pKopt,nExp_poi, sep="_"))
  return (tmp) 
}
library(DoubletFinder)

expr1 <- RDoublet(expr1)
expr1 <- subset(expr1, cells = colnames(expr1)[which(expr1[[]][12] == 'Singlet')])
expr1 <- subset(expr1, cells = colnames(expr1)[which(expr1[[]][13] == 'Singlet')])

#Wt_d100_1 <- expr1
#Wt_d100_2 <- expr1
#Wt_d100_3 <- expr1
#Wt_d100_4 <- expr1
#Wt_d170_1 <- expr1
#Wt_d170_2 <- expr1
#Wt_d170_3 <- expr1

Wt_d100_1 <- ProcessSeu(Wt_d100_1)
Wt_d100_1$stage <- 'Day100'
Wt_d100_2 <- ProcessSeu(Wt_d100_2)
Wt_d100_2$stage <- 'Day100'
Wt_d100_3 <- ProcessSeu(Wt_d100_3)
Wt_d100_3$stage <- 'Day100'
Wt_d100_4 <- ProcessSeu(Wt_d100_4)
Wt_d100_4$stage <- 'Day100'
Wt_d170_1 <- ProcessSeu(Wt_d170_1)
Wt_d170_1$stage <- 'Day170'
Wt_d170_2 <- ProcessSeu(Wt_d170_2)
Wt_d170_2$stage <- 'Day170'
Wt_d170_3 <- ProcessSeu(Wt_d170_3)
Wt_d170_3$stage <- 'Day170'

ProcessInt <- function(data.integrated){
  data.integrated <- ScaleData(data.integrated, verbose = T,  vars.to.regress = c('nCount_RNA', 'nFeature_RNA', 'percent.mt', "percent.rb","S.Score","G2M.Score"))
  data.integrated <- RunPCA(data.integrated, npcs = 30, verbose = FALSE)
  data.integrated <- FindNeighbors(data.integrated, dims = 1:20)
  data.integrated <- FindClusters(data.integrated, resolution = 0.5)
  data.integrated <- RunUMAP(data.integrated, reduction = "pca", dims = 1:20)
  data.integrated <- RunTSNE(data.integrated,  dims.use = 1:20 )
}

#Choose the objects for integration

integration_list <- list(Wt_d100_1, Wt_d100_2, Wt_d100_3, Wt_d100_4, Wt_d170_1, Wt_d170_2, Wt_d170_3)

features <- SelectIntegrationFeatures(object.list = integration_list)
data.anchors <- FindIntegrationAnchors(object.list = integration_list, anchor.features = features)

data.combined <- IntegrateData(anchorset = data.anchors)

data.combined <- ProcessInt(data.combined)

data.combined.markers <- FindAllMarkers(data.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

data.combined.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

data.combined.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(data.combined, features = top10$gene) + NoLegend()

data.combined <- RenameIdents(data.combined, '0' = 'Cones', '1' = 'Glia', '2' = 'Rods', '3' = 'Rods',
                              '4' = 'Cones','5'='Cones','6'='PR precursors', '7' = 'Glia', '8' = 'Progenitors',
                              '9' = 'PR precursors','10'='Cones','11'='Cones','12'='Glia','13'='Progenitors','14'='Cones',
                              '15'='Progenitors','16'='Glia','17'='Glia','18'='Glia','19'='Glia','20'='Glia')

FD125@active.ident <- FD125$seurat_clusters
FD125 <- RenameIdents(FD125, '0' = 'Rods', '1' = 'Progenitors', '2' = 'AC/HC', '3' = 'Rods', '4'='Progenitors',
                      '5'='Glia','6'='RGC','7'='Cones','8'='Bipolar','9'='Horizontal','10'='Muller glia','11'='Amacrine',
                      '12'='Cones','13'='Glia','14'='Glia')
