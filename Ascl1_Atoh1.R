#dogs Levi late dataset
levi.data <- Read10X_h5("G://Levi_Jon/late/Ascl1_Atoh1_late_timepoint/outs/filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
levi <- CreateSeuratObject(counts = levi.data, project = "late", min.cells = 3, min.features = 200)
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

levi [["percent.rb"]] <- PercentageFeatureSet(levi , pattern = "^Rps|^Rpl|^Mrps|^Mrpl", assay = 'RNA')

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

levi <- CellCycleScoring(levi, s.features = m.s.genes, g2m.features = m.g2m.genes, set.ident = FALSE)
levi[["percent.mt"]] <- PercentageFeatureSet(levi, pattern = "^mt-")

levi <- subset(levi, subset = nCount_RNA > 300 & nCount_RNA < 10000 & nFeature_RNA > 600 & nFeature_RNA < 4000 & percent.mt < 15 & percent.rb < 15)

levi <- ScaleData(levi, verbose = T, vars.to.regress = c('nCount_RNA', 'nFeature_RNA', 'percent.mt', "percent.rb","S.Score","G2M.Score"))

levi <- ProcessSeu(levi)

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

remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
library(DoubletFinder)

levi <- RDoublet(levi)
levi <- subset(levi, subset = DF.classifications_0.25_0.04_241 == 'Singlet')
levi <- subset(levi, subset = DF.classifications_0.25_0.04_191 == 'Singlet')

levi <- ProcessSeu(levi)

levi_late <- levi

#repeat for early timepoint

levi_early <- levi


