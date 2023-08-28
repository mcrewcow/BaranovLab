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
FD125$stage <- 'FD125'
adult <- RenameIdents(adult, 'Cone' = 'Cones', 'Rod' = 'Rods')
adult <- subset(adult, subset = annotation == c('Cones','Rods'))

integration_list <- list(Gamm170, FD125, adult)

features <- SelectIntegrationFeatures(object.list = integration_list)
data.anchors <- FindIntegrationAnchors(object.list = integration_list, anchor.features = features)

data.combined <- IntegrateData(anchorset = data.anchors)

data.combined <- ProcessInt(data.combined)
DimPlot(data.combined, group.by = 'annotation', label = T, repel = T, label.box = T)

data.combined.cones <- subset(data.combined, subset = annotation == 'Cones')
data.combined.rods <- subset(data.combined, subset = annotation == 'Rods')

AverageExpression(data.combined.cones, features = c('GFRA1','RET','NTRK2','DCC','SDC1','GFRA4','CD44','ACVR2A','ACVR2B','TGFBR1','TGFBR2','BMPR1A','FGFR8','FGFR1','FGFR2','NGFR','IGF1R','IGF2R'), group.by = c('stage'), assay = 'RNA')

mp_genes <- read.table(file = "clipboard", sep = "\t", header=TRUE) #Read the file with list of the genes
mp_genes <- as.list(mp_genes)
mp_genes <- lapply(mp_genes, toupper) 
mp_genes <- unlist(mp_genes)
mp <- unique(mp_genes)
[1] "ACKR1"   "ACKR2"   "ACKR3"   "ADGRL3"  "ASTN1"   "ASTN2"   "CCR10"   "CCR6"    "CCR9"   
[10] "CELSR2"  "CHL1"    "CNTN4"   "CTNNA2"  "CTNNB1"  "CX3CR1"  "CXCR4"   "CXCR6"   "FAT3"   
[19] "FGFR1"   "FZD3"    "FZD5"    "LRP12"   "NCAM2"   "NLGN1"   "NLGN2"   "NLGN3"   "NLGN4X" 
[28] "NRCAM"   "NRXN1"   "NRXN2"   "NRXN3"   "NTRK3"   "PLXNB2"  "RBM15"   "ROBO1"   "SEMA6A" 
[37] "SLC12A2" "STK39"   "UNC5C"   "UNC5D"   "WBP1L"   "WNK1"  
AverageExpression(data.combined.rods, features = mp, group.by = c('stage'), assay = 'RNA')

df <- read.table(file = "clipboard", sep = "\t", header=TRUE)
ggplot(df, aes(x = reorder(GENE, -Value), y = Stage, fill= Value)) +
geom_tile(color = "black",
lwd = 1,linetype = 1) +
scale_fill_gradient2(low = "#075AFF",
mid = "white",
high = "#FF0000", midpoint =0.5) + theme_bw() + coord_fixed() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab('Gene expression in the integated fetal, adult, and organoids cones')

df$Stage <- factor(df$Stage, levels = c('Adult','OD170','FD125'))
avg_expression <- df %>%
  group_by(GENE) %>%
  summarize(avg_expression = mean(Value)) %>%
  arrange(desc(avg_expression))

df$GENE <- factor(df$GENE, levels = avg_expression$GENE)

ggplot(df, aes(x = GENE, y = Stage, fill= Value)) +
  geom_tile(color = "black",
            lwd = 1,linetype = 1) +
  scale_fill_gradient2(low = "#075AFF",
                       mid = "white",
                       high = "#FF0000", midpoint =0.5) + theme_bw() + coord_fixed() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab('Gene expression in the integated fetal, adult, and organoids rods')

AMD1 <- readRDS('G://human_AMD/Menon_2019.rds')
AMD1$nCount_RNA = colSums(x = AMD1, slot = "data")  # nCount_RNA
AMD1$nFeature = colSums(x = GetAssayData(object = AMD1, slot = "data") > 0)
AMD1 <- subset(AMD1, subset = cell_type == c('retinal bipolar neuron','retinal cone cell','retina horizontal cell',
                                             'macroglial cell','microglial cell','retinal rod cell'))
AMD1 <- subset(AMD1, subset = cell_type == c('retinal ganglion cell'), invert = T) #'amacrine cell','retinal blood vessel endothelial cell','retinal ganglion cell'
AMD1$cell_type = droplevels(AMD1$cell_type, exclude = setdiff(levels(AMD1$cell_type),unique(AMD1$cell_type)))
AMD1$cell_type1 <- paste(AMD1$cell_type,'Menon')

AMD2 <- readRDS('G://human_AMD/Orozco_2020.rds')
AMD2$nCount_RNA = colSums(x = AMD2, slot = "data")  # nCount_RNA
AMD2$nFeature = colSums(x = GetAssayData(object = AMD2, slot = "data") > 0)
AMD2 <- subset(AMD2, subset = cell_type == c('retinal bipolar neuron','retinal cone cell','retina horizontal cell',
                                             'Mueller cell','retinal rod cell'))
AMD2 <- subset(AMD2, subset = cell_type == c('retinal pigment epithelial cell'), invert = T) #,'native cell','astrocyte','amacrine cell',
#'retinal ganglion cell','retinal pigment epithelial cell','myeloid cell'
AMD2$cell_type = droplevels(AMD2$cell_type, exclude = setdiff(levels(AMD2$cell_type),unique(AMD2$cell_type)))
AMD2$cell_type1 <- paste(AMD2$cell_type,'Orozco')
integration_list <- list(AMD1, AMD2)

features <- SelectIntegrationFeatures(object.list = integration_list)
data.anchors <- FindIntegrationAnchors(object.list = integration_list, anchor.features = features)

data.combined <- IntegrateData(anchorset = data.anchors)
ProcessInt <- function(data.integrated){
  data.integrated <- ScaleData(data.integrated, verbose = T)
  data.integrated <- RunPCA(data.integrated, npcs = 30, verbose = FALSE)
  data.integrated <- FindNeighbors(data.integrated, dims = 1:20)
  data.integrated <- FindClusters(data.integrated, resolution = 0.5)
  data.integrated <- RunUMAP(data.integrated, reduction = "pca", dims = 1:20)
  data.integrated <- RunTSNE(data.integrated,  dims.use = 1:20 )
}

data.combined <- ProcessInt(data.combined)
DimPlot(data.combined, group.by = 'cell_type1', label = T, repel = T, label.box = T)
data.combined <- SetIdent(data.combined, value = 'cell_type1')
DimPlot(data.combined, label = T, repel = T, label.box = T)
data.combined <- SetIdent(data.combined, value = 'cell_type1')
data.combined <- RenameIdents(data.combined, 'retinal rod cell Menon' = 'Rods', 'retinal rod cell Orozco' = 'Rods',
                              'microglial cell Menon' = 'Microglia', 'Mueller cell Orozco' = 'Mueller glia',
                              'retina horizontal cell Menon' = 'Horizontal', 'retina horizontal cell Orozco' = 'Horizontal',
                              'retinal cone cell Orozco' = 'Cones', 'retinal cone cell Menon' = 'Cones',
                              'retinal bipolar neuron Orozco' = 'Bipolar', 'retinal bipolar neuron Menon' = 'Bipolar',
                              'macroglial cell Menon' = 'Glia')

DimPlot(data.combined, label = T, repel = T, label.box = T)
data.combined$EK_anno <- data.combined@active.ident
data.combined$stage <- 'AMD'
data.combined$annotation <- data.combined@active.ident

DefaultAssay(FD125) <- 'RNA'
DefaultAssay(AMD1) <- 'RNA'
AMD2$stage <- 'AMD'
AMD2$annotation <- AMD2$cell_type
DefaultAssay(data.combined) <- 'integrated'
DefaultAssay(adult) <- 'RNA'
integration_list <- list(AMD1, FD125, adult)
features <- SelectIntegrationFeatures(object.list = integration_list)
data.anchors <- FindIntegrationAnchors(object.list = integration_list, anchor.features = features)
data.combined <- IntegrateData(anchorset = data.anchors)
data.combined <- ProcessInt(data.combined)

AMD1 <- readRDS('G://human_AMD/Menon_2019.rds')
write.csv(AMD1@meta.data,"G://human_AMD/Menon/seurat_metadata.csv")
expression_matrix <- read.delim('G://human_AMD/Orozco/GSM3988006_SAM24362284.txt.gz', header = T, stringsAsFactors = F)
expression_matrix <- ReadMtx(mtx = "G://human_AMD/Menon/GSE137537_counts.mtx", features = "G://human_AMD/Menon/test.tsv",cells = "G://human_AMD/Menon/seurat_metadata.tsv")

library(Seurat)
library(patchwork)
library(dplyr)
library(biomaRt)
library(org.Hs.eg.db)
library(ggplot2)
library(Matrix)

ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
bm <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), values=rownames(expression_matrix), mart=ensembl)
print("Mitochondrial genes before conversion:")
print(grep("^MT-", rownames(expression_matrix), value = TRUE))
hgnc.symbols <- bm$hgnc_symbol[match(rownames(expression_matrix), bm$ensembl_gene_id)]
expression_matrix <- as.matrix(expression_matrix)
rownames(expression_matrix) <- hgnc.symbols
print("Mitochondrial genes after conversion:")
print(grep("^MT-", rownames(expression_matrix), value = TRUE))

AMD1 <- CreateSeuratObject(counts = expression_matrix, min.cells = 3, min.features = 100)
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

ProcessSeu <- function(Seurat){
  Seurat <- NormalizeData(Seurat)
  Seurat <- FindVariableFeatures(Seurat, selection.method = "vst", nfeatures = 2000)
  Seurat <- ScaleData(Seurat, verbose = T, vars.to.regress = c('nCount_RNA', 'nFeature_RNA', "percent.rb","S.Score","G2M.Score"))
  Seurat <- RunPCA(Seurat)
  Seurat <- FindNeighbors(Seurat, dims = 1:20)
  Seurat <- FindClusters(Seurat, resolution = 0.5)
  Seurat <- RunUMAP(Seurat, dims = 1:20)
  Seurat <- RunTSNE(Seurat,  dims.use = 1:20 )
  DimPlot(object = Seurat, reduction = "umap")
  return (Seurat)
}

AMD1 [["percent.rb"]] <- PercentageFeatureSet(AMD1 , pattern = "^RPS|^RPL|^MRPS|^MRPL", assay = 'RNA')
AMD1 <- CellCycleScoring(AMD1, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE, nbin = 3)
AMD1[["percent.mt"]] <- PercentageFeatureSet(AMD1, pattern = "^MT-")
VlnPlot(AMD1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", 'percent.rb'), ncol = 4)

AMD1 <- subset(AMD1, subset = nCount_RNA > 100 & nCount_RNA < 10000 & nFeature_RNA > 100 & nFeature_RNA < 4000 & percent.mt < 15 & percent.rb < 35)

AMD1 <- ScaleData(AMD1, verbose = T, vars.to.regress = c('nCount_RNA', 'nFeature_RNA', "percent.rb","S.Score","G2M.Score"))

AMD1 <- ProcessSeu(AMD1)

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

AMD1 <- RDoublet(AMD1)
AMD1 <- subset(AMD1, cells = colnames(AMD1)[which(AMD1[[]][12] == 'Singlet')])
AMD1 <- subset(AMD1, cells = colnames(AMD1)[which(AMD1[[]][13] == 'Singlet')])

AMD1 <- ProcessSeu(AMD1)
DimPlot(AMD1)
AMD1.markers <- FindAllMarkers(AMD1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

AMD1.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

AMD1.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(AMD1, features = top10$gene) + NoLegend()


# Install and load the biomaRt package
if (!requireNamespace("biomaRt", quietly = TRUE)) {
  install.packages("biomaRt")
}
library(biomaRt)
# Connect to the Ensembl database
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# List of gene names
gene_names <- c("NLGN1", "NLGN2", "NCAM1", "FGF1", "FGF2", "FGF3", "FGF4", "FGF5", "FGF6", "FGF7", "FGF8", "FGF9",
                "FGF10", "FGF15", "FGF16", "FGF17", "FGF18", "FGF20", "FGF21", "FGF22", "FGF23", "CCL2", "CCL3",
                "CCL3L1", "CCL4", "CCL5", "CCL7", "CCL8", "CCL11", "CCL13", "CCL22", "CCL24", "CCL26", "CCL14",
                "CCL17", "CCL27", "CCL28", "CXCL11", "CXCL12", "MIF", "CX3CL1")
# Get Ensembl IDs for gene names
ensembl_ids <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                     filters = "external_gene_name",
                     values = gene_names,
                     mart = ensembl)
# Print the Ensembl IDs
print(ensembl_ids)

write.table(as.matrix(GetAssayData(object = AMDfull, slot = "data")), 
            'G://Human_AMD/merged_counts.csv', 
            sep = ',', row.names = T, col.names = T, quote = F)
AverageExpression(AMDfull, features = ensembl_ids$ensembl_gene_id, group.by = c('stage'), assay = 'RNA')

df <- read.table(file = "clipboard", sep = "\t", header=TRUE)
ggplot(df, aes(x = reorder(GENE, -Value), y = Stage, fill= Value)) +
  geom_tile(color = "black",
            lwd = 1,linetype = 1) +
  scale_fill_gradient2(low = "#075AFF",
                       mid = "white",
                       high = "#FF0000", midpoint =0.5) + theme_bw() + coord_fixed() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab('Gene expression in the integated fetal, adult, and organoids cones')

df$Stage <- factor(df$Stage, levels = c('Adult','OD170','FD125'))
avg_expression <- df %>%
  group_by(GENE) %>%
  summarize(avg_expression = mean(Value)) %>%
  arrange(desc(avg_expression))

df$GENE <- factor(df$GENE, levels = avg_expression$GENE)

ggplot(df, aes(x = GENE, y = Stage, fill= Value)) +
  geom_tile(color = "black",
            lwd = 1,linetype = 1) +
  scale_fill_gradient2(low = "#075AFF",
                       mid = "white",
                       high = "#FF0000", midpoint =0.5) + theme_bw() + coord_fixed() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab('Gene expression in the integated fetal, adult, and organoids rods')

expression_matrix <- read.csv('G://human_AMD//merged_counts.csv', header = T, stringsAsFactors = F)

library(Seurat)
library(patchwork)
library(dplyr)
library(biomaRt)
library(org.Hs.eg.db)
library(ggplot2)
library(Matrix)

ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
bm <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), values=rownames(expression_matrix), mart=ensembl)
print("Mitochondrial genes before conversion:")
print(grep("^MT-", rownames(expression_matrix), value = TRUE))
hgnc.symbols <- bm$hgnc_symbol[match(rownames(expression_matrix), bm$ensembl_gene_id)]
expression_matrix <- as.matrix(expression_matrix)
rownames(expression_matrix) <- hgnc.symbols
print("Mitochondrial genes after conversion:")
print(grep("^MT-", rownames(expression_matrix), value = TRUE))
expression_matrix <- expression_matrix[!(is.na(rownames(expression_matrix)) | rownames(expression_matrix) == ' '), ] #and ""
unique_rows <- !duplicated(rownames(expression_matrix))
expr1 <- expression_matrix[unique_rows, ]

AMD1 <- CreateSeuratObject(counts = expr1)




orgreh$stage <- 'OD4560'
orgreh60$stage <- 'OD4560'
FD59$stage <- 'FD59'
adult$stage <- 'Adult'
FD59$annotation <- FD59@active.ident
adult$annotation <- adult@active.ident
orgreh$annotation <- orgreh@active.ident
orgreh60$annotation <- orgreh60@active.ident
FD59 <- subset(FD59, subset = annotation == 'RGC')
adult <- subset(adult, subset = annotation == 'RGC')
orgreh <- subset(orgreh, subset = annotation == c('T1'), invert = T)
orgreh$annotation = droplevels(orgreh$annotation, exclude = setdiff(levels(orgreh$annotation),unique(orgreh$annotation)))

orgreh60 <- subset(orgreh60, subset = annotation == 'RGC')   

integration_list <- list(FD59, adult, orgreh, orgreh60)
features <- SelectIntegrationFeatures(object.list = integration_list)
data.anchors <- FindIntegrationAnchors(object.list = integration_list, anchor.features = features)
RGCs_merged_orgreh4560_fd59_adult <- IntegrateData(anchorset = data.anchors)
RGCs_merged_orgreh4560_fd59_adult <- ProcessInt(RGCs_merged_orgreh4560_fd59_adult)

DimPlot(RGCs_merged_orgreh4560_fd59_adult, group.by = 'stage')
mp_genes <- read.table(file = "clipboard", sep = "\t", header=TRUE) #Read the file with list of the genes
mp_genes <- as.list(mp_genes)
mp_genes <- lapply(mp_genes, toupper) 
mp_genes <- unlist(mp_genes)
mp <- unique(mp_genes)
df <- AverageExpression(RGCs_merged_orgreh4560_fd59_adult, features = mp, group.by = c('stage'), assay = 'RNA')

df <- read.table(file = "clipboard", sep = "\t", header=TRUE)
ggplot(df, aes(x = reorder(GENE, -Value), y = Stage, fill= Value)) +
  geom_tile(color = "black",
            lwd = 1,linetype = 1) +
  scale_fill_gradient2(low = "#075AFF",
                       mid = "white",
                       high = "#FF0000", midpoint =0.5) + theme_bw() + coord_fixed() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab('Gene expression in the integated fetal, adult, and organoids cones')

df$Stage <- factor(df$Stage, levels = c('Adult','OD45_60','FD59'))
avg_expression <- df %>%
  group_by(GENE) %>%
  summarize(avg_expression = mean(Value)) %>%
  arrange(desc(avg_expression))

df$GENE <- factor(df$GENE, levels = avg_expression$GENE)

ggplot(df, aes(x = GENE, y = Stage, fill= Value)) +
  geom_tile(color = "black",
            lwd = 1,linetype = 1) +
  scale_fill_gradient2(low = "#075AFF",
                       mid = "white",
                       high = "#FF0000", midpoint =0.5) + theme_bw() + coord_fixed() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab('Gene expression in the integated fetal, adult, and organoids RGCs')

expression_matrix <- read.csv('G://Mandeep_dataset/counts.csv')
View(expression_matrix)
expr1 <- CreateSeuratObject(counts = expression_matrix)
metadata <- read.csv('G://Mandeep_dataset/GSE197847_metadata.csv', row.names = 1)
expr1 <- AddMetaData(expr1, metadata = metadata)
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
expr1 [["percent.rb"]] <- PercentageFeatureSet(expr1 , pattern = "^RPS|^RPL|^MRPS|^MRPL", assay = 'RNA')
expr1 <- CellCycleScoring(expr1, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE, nbin = 3)
expr1[["percent.mt"]] <- PercentageFeatureSet(expr1, pattern = "^MT-")
VlnPlot(expr1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", 'percent.rb'), ncol = 4)
expr1 <- subset(expr1, subset = nCount_RNA > 100 & nCount_RNA < 35000 & nFeature_RNA > 100 & nFeature_RNA < 7000 & percent.mt < 20 & percent.rb < 35)
expr1 <- ScaleData(expr1, verbose = T, vars.to.regress = c('nCount_RNA', 'nFeature_RNA', 'percent.mt', "percent.rb","S.Score","G2M.Score"))
expr1 <- ProcessSeu(expr1)
library(DoubletFinder)
expr1 <- RDoublet(expr1)

FD125@active.ident <- FD125$seurat_clusters
FD125 <- RenameIdents(FD125, '0' = 'Rods', '1' = 'Progenitors', '2' = 'AC/HC', '3' = 'Rods', '4'='Progenitors',
                      '5'='Glia','6'='RGC','7'='Cones','8'='Bipolar','9'='Horizontal','10'='Muller glia','11'='Amacrine',
                      '12'='Cones','13'='Glia','14'='Glia')
