ONC_1d <-  read.csv('C://Bioinf/Zhigang_ONC/GSE137398_ONCRGCs_1w_afterCrush_count_mat.csv')
rownames(ONC_1d) <- ONC_1d$X
ONC_1d <- ONC_1d[,-1]
ONC_1d <- CreateSeuratObject(counts = ONC_1d, min.cells = 3, min.features = 100)
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

convertHumanGeneList <- function(x){
  require("biomaRt")
  human = useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", mirror = "www", host = "dec2021.archive.ensembl.org")
  mouse = useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl", mirror = "www", host = "dec2021.archive.ensembl.org")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}
m.s.genes <- convertHumanGeneList(cc.genes.updated.2019$s.genes)
m.g2m.genes <- convertHumanGeneList(cc.genes.updated.2019$g2m.genes)

ProcessSeu <- function(Seurat){
  Seurat <- NormalizeData(Seurat)
  Seurat <- FindVariableFeatures(Seurat, selection.method = "vst", nfeatures = 2000)
  Seurat <- ScaleData(Seurat, verbose = T, vars.to.regress = c('nCount_RNA', 'nFeature_RNA', 'percent.mt',"percent.rb","S.Score","G2M.Score"))
  Seurat <- RunPCA(Seurat)
  Seurat <- FindNeighbors(Seurat, dims = 1:30)
  Seurat <- FindClusters(Seurat, resolution = 1)
  Seurat <- RunUMAP(Seurat, dims = 1:30)
  Seurat <- RunTSNE(Seurat,  dims.use = 1:30 )
  DimPlot(object = Seurat, reduction = "umap")
  return (Seurat)
}
ONC_1d [["percent.rb"]] <- PercentageFeatureSet(ONC_1d , pattern = "^Rps|^Rpl|^Mrps|^Mrpl", assay = 'RNA')
ONC_1d <- CellCycleScoring(ONC_1d, s.features = m.s.genes, g2m.features = m.g2m.genes, set.ident = FALSE)
ONC_1d[["percent.mt"]] <- PercentageFeatureSet(ONC_1d, pattern = "^mt-")
VlnPlot(ONC_1d, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", 'percent.rb'), ncol = 4)
ONC_1d <- subset(ONC_1d, subset = nCount_RNA > 100 & nCount_RNA < 70000 & nFeature_RNA > 100 & nFeature_RNA < 9000 & percent.mt < 20 & percent.rb < 35)
ONC_1d <- ScaleData(ONC_1d, verbose = T, vars.to.regress = c('nCount_RNA', 'nFeature_RNA','percent.mt', "percent.rb","S.Score","G2M.Score"))
gc()
ONC_1d <- ProcessSeu(ONC_1d)
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

ONC_1d <- RDoublet(ONC_1d)
ONC_1d <- subset(ONC_1d, cells = colnames(ONC_1d)[which(ONC_1d[[]][12] == 'Singlet')])
ONC_1d <- subset(ONC_1d, cells = colnames(ONC_1d)[which(ONC_1d[[]][13] == 'Singlet')])
gc()
ONC_1d <- ProcessSeu(ONC_1d)
ONC_1d$stage <- 'ONC 1 week ZH'
SaveH5Seurat(ONC_1d, 'C://Bioinf/Zhigang_ONC/ONC_1w_EK.h5Seurat')
